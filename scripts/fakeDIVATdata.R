#### Recreating data like the Kidney data in Foucher et al. (2010)
setwd("H:/Multistate models/SemiMarkovMultistate")
rm(list = ls())
library(igraph)
library(parallel)
library(cubature)
library(MASS)

source("sm_msm_preprocessing_func.R")
source("sm_msm_likelihood_func.R")

nn <- 839

## 2 absorbing states!
gg = graph_from_literal("1"--+"2"--+"3", "1"--+"4", "2"--+"4","1"--+"3")

## Generate covariates
XX <- matrix(NA,nn,9)
colnames(XX) <- c("donor_gender","patient_gender","cold_isch_time","donor_age","patient_age",
                  "num_HLAincomp","ind_treat","panel_react_antib","del_graft_func")
Sig <- matrix(NA,9,9)
Sig[lower.tri(Sig,diag=TRUE)] <- c(1,0,0,0,0,0.1,0.1,0.1,0.1,1,-0.1,0,0,-0.3,0.2,0.1,0.1,
                                   1,0,0.4,0.2,0.1,-0.1,0.1,1,0,0,0,0,0,
                                   1,0.6,0.4,0.3,0.5,1,0.2,0.3,0.4,1,0.2,0.2,1,0.3,1)
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
Sig <- makeSymm(Sig)
ZZ <- mvrnorm(nn,rep(0,9),Sig)
XX[,1] <- as.numeric(ZZ[,1]>0)
XX[,2] <- as.numeric(ZZ[,2]>-0.1)
XX[,3] <- as.numeric(ZZ[,3]>1)
XX[,4] <- as.numeric(ZZ[,4]>0.1)
XX[,5] <- as.numeric(ZZ[,5]>-0.3)
XX[,6] <- as.numeric(ZZ[,6]>0)
XX[,7] <- as.numeric(ZZ[,7]>0.1)
XX[,8] <- as.numeric(ZZ[,8]>0.1)
XX[,9] <- as.numeric(ZZ[,9]>0.1)

## Generate durations
qGenWei <- function(uval,nu,sig,the){
  m <- length(uval)
  dval <- rep(NA,m)
  for (i in 1:m){
    dval[i] <- sig*((1-log(1-uval[i]))^the-1)^(1/nu)
  }
  return(dval)
}
rGenWei <- function(n,nu,sig,the){
  uu <- runif(n)
  qGenWei(uu,nu,sig,the)
}
dd1 <- rGenWei(nn,0.77,68.8,0.25)
hist(dd1)

qGenWeiX <- function(uval,nu,sig,the,beta,x){
  m <- length(uval)
  dval <- rep(NA,m)
  for (i in 1:m){
    dval[i] <- sig*((1-exp(-beta%*%x)*log(1-uval[i]))^the-1)^(1/nu)
  }
  return(dval)
}
rGenWeiX <- function(n,nu,sig,the,beta,x){
  uu <- runif(n)
  qGenWeiX(uu,nu,sig,the,beta,x)
}

dd1 <- rGenWeiX(nn,0.77,68.8,0.25,c(0.36,-0.26,0.96),XX[1,c("patient_gender","donor_age","ind_treat")])
hist(dd1)

## Make transition probabilities
pMat <- matrix(NA,nn,5)
colnames(pMat) <- c("p12","p13","p14","p23","p24")
pMat[,"p12"] <- exp(1.87+0.73*XX[,"del_graft_func"]-2.03*XX[,"donor_age"])
pMat[,"p13"] <- exp(-2.81)
pMat[,"p14"] <- exp(0)
rsums <- rowSums(pMat[,c("p12","p13","p14")])
pMat[,"p12"] <- pMat[,"p12"]/rsums
pMat[,"p13"] <- pMat[,"p13"]/rsums
pMat[,"p14"] <- pMat[,"p14"]/rsums

pMat[,"p23"] <- exp(1.13+0.93*XX[,"num_HLAincomp"])
pMat[,"p24"] <- exp(0)
rsums <- rowSums(pMat[,c("p23","p24")])
pMat[,"p23"] <- pMat[,"p23"]/rsums
pMat[,"p24"] <- pMat[,"p24"]/rsums


## Generate transition history for each patient
tMat <- matrix(NA,nn,5)
colnames(tMat) <- c("t12","t13","t14","t23","t24")

for (i in 1:nn){
  draw1 <- rmultinom(1,1,as.numeric(pMat[i,c("p12","p13","p14")]))
  if (draw1[1,1]==1){
    dd12 <- rGenWeiX(1,0.77,68.8,0.25,c(-0.26,0.96,0.36),XX[i,c("patient_gender","donor_age","ind_treat")])
    tMat[i,"t12"] <- dd12
    
    draw2 <- rbinom(1,1,pMat[i,"p23"])
    if (draw2==1){
      dd23 <- rGenWeiX(1,1,12.66,1,c(0.88,1.09),XX[i,c("num_HLAincomp","panel_react_antib")])
      tMat[i,"t23"] <- dd23
    }else if (draw2==0){
      dd24 <- rGenWeiX(1,1,6.54,1,c(2.03,1.54),XX[i,c("del_graft_func","patient_gender")])
      tMat[i,"t24"] <- dd24
    }
  }else if (draw1[2,1]==1){
    dd13 <- rGenWeiX(1,1,42.84,1,5.02,XX[i,"cold_isch_time",drop=F])
    tMat[i,"t13"] <- dd13
  }else if (draw1[3,1]==1){
    dd14 <- rGenWei(1,1,109.83,1)
    tMat[i,"t14"] <- dd14
  }
}

## Generate observation times and states
tMat_full <- matrix(NA,nn*30,3)
colnames(tMat_full) <- c("patient","time","state")
ii <- 1
for (i in 1:nn){
  mm <- sample(5:55,1)
  #tobs <- c(0,cumsum(rgamma(mm-1,0.5,1)))
  #tobs <- c(0,sort(rgamma(mm-1,0.5,0.3)))
  tobs <- c(0,sort(rlnorm(mm-1,-1,1)))
  iistart <- ii
  for (j in 1:mm){
    if (tobs[j] <= max(tMat[i,c("t12","t13","t14")],na.rm=T)){
      tMat_full[ii,c("time","state")] <- c(tobs[j],1)
      ii <- ii+1
    }else if (!is.na(tMat[i,c("t13")])){
      tMat_full[ii,c("time","state")] <- c(tobs[j],3)
      ii <- ii+1
      break
    }else if (!is.na(tMat[i,c("t14")])){
      tMat_full[ii,c("time","state")] <- c(tobs[j],4)
      ii <- ii+1
      break
    }else if (tobs[j] > tMat[i,c("t12")] & tobs[j] <= tMat[i,c("t12")]+max(tMat[i,c("t23","t24")],na.rm=T)){
      tMat_full[ii,c("time","state")] <- c(tobs[j],2)
      ii <- ii+1
    }else if(!is.na(tMat[i,c("t23")])){
      tMat_full[ii,c("time","state")] <- c(tobs[j],3)
      ii <- ii+1
      break
    }else if(!is.na(tMat[i,c("t24")])){
      tMat_full[ii,c("time","state")] <- c(tobs[j],4)
      ii <- ii+1
      break
    }
  }
  tMat_full[iistart:(ii-1),"patient"] <- i
}
tMat_full <- tMat_full[-which(is.na(tMat_full[,"patient"])),]

dim(tMat_full)[1]/nn
tMat_agg1 <- aggregate(tMat_full[,"state"],by=list(tMat_full[,"patient"]),FUN=unique)
tMat_agg2 <- aggregate(tMat_full[,"time"],by=list(tMat_full[,"patient"]),FUN=max)
tMat_agg <- merge(tMat_agg1,tMat_agg2,by="Group.1")
colnames(tMat_agg) <- c("patient","states","mean_time")
tMat_agg$states <- as.character(tMat_agg$states)
table(tMat_agg$states)

aggregate(tMat_agg$mean_time,by=list(tMat_agg$states),FUN=mean)

### Fit true model
ddt <- arrange_data(data.frame(tMat_full),gg, abs_int_cens = NULL)
table(ddt[,"obs_type"])

