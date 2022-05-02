### Fitting four-state illness death to the CAV data
setwd("H:/Multistate models")
source("weibull_4state_illness_death.R")

## CAV data
library(msm)
dd <- cav
dd <- dd[!is.na(dd$pdiag),]
dim(dd[dd$firstobs==1,]) # 614 patients
#first obs (at years=0) is always in state 1
id_wrong <- unique(dd$PTNUM[which(dd$state!=dd$statemax)])  # observations where the patient appears to go back to a previous state
dd <- dd[-which(dd$PTNUM %in% id_wrong),]

## Keep only relevant time-points
inds <- unique(dd$PTNUM)
nn <- length(inds)
ddr <- dd
idd <- NULL
pt_types <- list()

for (i in 1:nn){
  ddi <- dd[which(dd$PTNUM==inds[i]),]
  states <- ddi$state
  rlei <- matrix(0,2,4)
  rlei[1,] <- 1:4
  rlei[2,rle(states)$values] <- rle(states)$lengths
  id_unrelevant <- NULL
  if (rlei[2,1]>=2){
    id_unrelevant <- 1:(rlei[2,1]-1)
  }
  if (rlei[2,2]>=3){
    id_unrelevant <- c(id_unrelevant,(rlei[2,1]+2):(rlei[2,1]+rlei[2,2]-1))
  }
  if (rlei[2,3]>=3){
    id_unrelevant <- c(id_unrelevant,(rlei[2,1]+rlei[2,2]+2):(sum(rlei[2,1:3])-1))
  }
  if (rlei[2,4]>=2){
    id_unrelevant <- c(id_unrelevant,(sum(rlei[2,1:3])+2):sum(rlei[2,1:4]))
  }
  idd <- c(idd,which(dd$PTNUM==inds[i])[id_unrelevant])
  #pt_types[[i]] <- idtypes(rle(states)$values)
}

ddr <- ddr[-idd,]

# Make matrix with patients and time-points
idtypes <- function(tpoints){ # Identifies one of eight observation types
  # timepointss: a vector of length 6 with NAs indicating unobserved states
  kk <- sum(!is.na(tpoints))
  observed <- which(!is.na(tpoints))
  otype <- NA
  if (kk == 1){ #Type "0"
    otype <- 1 
  }else if (kk==2){#Types "0(12)3", "0(1)3" or "03"
    otype <- 4
  }else if (kk==3){
    if (sum(observed==c(1,2,3))==kk){ #Type "01"
      otype <- 2 
    }else if (sum(observed==c(1,4,5))==kk){ #Type "0(1)2"
      otype <- 3
    }
  }else if (kk==4){
    if (sum(observed==c(1,2,3,6))==kk){ #Types "01(2)3" or "013"
      otype <- 6 
    }else if (sum(observed==c(1,4,5,6))==kk){ #Type "0(1)23"
      otype <- 7 
    }
  }else if(kk==5){ #Type "012"
    otype <- 5 
  }else if (kk == 6){ #Type "0123"
    otype <- 8 
  }
  return(otype)
}

inds <- unique(ddr$PTNUM)
nn <- length(inds)
timepoints <- matrix(NA,nn,6)
colnames(timepoints) <- c("t0M","t1m","t1M","t2m","t2M","t3m")
otypes <- rep(NA,nn)

for (i in 1:nn){
  ddi <- ddr[which(ddr$PTNUM==inds[i]),]
  tti <- ddi$years[pmatch(c(1,2,2,3,3,4),ddi$state)]
  if (is.na(tti[3]) & !is.na(tti[2])) tti[3] <- tti[2]
  if (is.na(tti[5]) & !is.na(tti[4])) tti[5] <- tti[4]
  timepoints[i,] <- tti
  otypes[i] <- idtypes(timepoints[i,])
}

# Make log-likelihood
mll <- function(param,timeMat){
  # timeMat: matrix of relevant time-points for each patient
  nn <- dim(timeMat)[1]
  lli <- rep(NA,nn)
  for (i in 1:nn){
    otype <- idtypes(timepoints[i,])
    if (otype==1){
      lli[i] <- otype_1(param,timepoints[i,])
    }else if(otype==2){
      lli[i] <- otype_2(param,timepoints[i,])
    }else if(otype==3){
      lli[i] <- otype_3(param,timepoints[i,])
    }else if(otype==4){
      lli[i] <- otype_4(param,timepoints[i,])
    }else if(otype==5){
      lli[i] <- otype_5(param,timepoints[i,])
    }else if(otype==6){
      lli[i] <- otype_6(param,timepoints[i,])
    }else if(otype==7){
      lli[i] <- otype_7(param,timepoints[i,])
    }else if(otype==8){
      lli[i] <- otype_8(param,timepoints[i,])
    }
  }
  return(sum(lli))
}

aa <- c(1,12,1,3,1,3,1,25,1,16)
mll(aa,timepoints)

params <- c(1.43,8.72,1.21,3.05,1.12,3.75,0.44,418.32,0.71,6.78)
mll(params,timepoints)

system.time({
  mll(params,timepoints)
})


oo <- nlminb(aa,mll,timeMat=timepoints,lower=rep(0.01,10))

oo <- optim(aa,mll,timeMat=timepoints,method = "L-BFGS-B",lower=rep(0.0001,10),hessian=T)

aa2 <- c(1.4343907,8.7152152,1.2118153,3.0472622,1.1245528,3.7461502,0.4360518,418.3155272,0.7088931,6.7814314)

hss <- pracma::hessian(mll,aa2,timeMat=timepoints)
round((oo$par-1)/sqrt(diag(solve(oo$hessian))),3) #5.063  17.074   1.396   6.495   0.859   5.736 -11.897   1.544  -1.762   3.017. 
# seems that T12 could have been expo, and also T23.

2*mll(aa2,timepoints)+2*10 # 2829.4, better than all the models Marthe tried?

## Occupancy probabilities
aa <- c(1.4,8.6,1.3,2.7,0.9,3,0.39,1166,0.25,284.16)
aa <- aa2
tval <- seq(0.01,20,length=100)
pi0f <- pi0(aa,tval)
pi1f <- pi1(aa,tval)
pi2f <- pi2(aa,tval)
pi3f <- pi3(aa,tval)
plot(tval,pi0f,type="l",ylim=c(0,1))
lines(tval,pi1f,col="red")
lines(tval,pi2f,col="green")
lines(tval,pi3f,col="blue")

plot(tval,pi0f+pi1f+pi2f+pi3f,type="l")

## Embedded transition probabilities


# Homogenous markov
twoway4.q <- rbind(c(0, 0.25, 0, 0.25), c(0.166, 0, 0.166, 0.166),c(0, 0.25, 0, 0.25), c(0, 0, 0, 0))
rownames(twoway4.q) <- colnames(twoway4.q) <- c("Well", "Mild","Severe", "Death")
cav.msm <- msm(state ~ years, subject = PTNUM, data = dd,qmatrix = twoway4.q, death = 4)

prev <- prevalence.msm(cav.msm,times=tval)


pdf("cav_weibull_prevalence.pdf", width=11, height=7)
par(mfrow=c(2,2))
par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1.2)
plot(tval,prev$`Observed percentages`[,1],type="l",col="black",ylim=c(0,100),lwd=2,xlab=" ",
     ylab="prevalence (%)",main="State 0")
lines(tval,prev$`Expected percentages`[,1],col="#ca0020",lwd=3)
lines(tval,pi0f*100,col="#0571b0",lwd=3)
legend("topright",legend=c("observed","Markov model","MCN model"),col=c("black","#ca0020","#0571b0"),lwd=3,bty="n")

plot(tval,prev$`Observed percentages`[,2],type="l",col="black",ylim=c(0,100),lwd=2,xlab=" ",
     ylab="prevalence (%)",main="State 1")
lines(tval,prev$`Expected percentages`[,2],col="#ca0020",lwd=3)
lines(tval,pi1f*100,col="#0571b0",lwd=3)
legend("topright",legend=c("observed","Markov model","MCN model"),col=c("black","#ca0020","#0571b0"),lwd=3,bty="n")

plot(tval,prev$`Observed percentages`[,3],type="l",col="black",ylim=c(0,100),lwd=2,xlab="years after transplantation",
     ylab="prevalence (%)",main="State 2")
lines(tval,prev$`Expected percentages`[,3],col="#ca0020",lwd=3)
lines(tval,pi2f*100,col="#0571b0",lwd=3)
legend("topright",legend=c("observed","Markov model","MCN model"),col=c("black","#ca0020","#0571b0"),lwd=3,bty="n")

plot(tval,prev$`Observed percentages`[,4],type="l",col="black",ylim=c(0,100),lwd=2,xlab="years after transplantation",
     ylab="prevalence (%)",main="State 3")
lines(tval,prev$`Expected percentages`[,4],col="#ca0020",lwd=3)
lines(tval,pi3f*100,col="#0571b0",lwd=3)
legend("bottomright",legend=c("observed","Markov model","MCN model"),col=c("black","#ca0020","#0571b0"),lwd=3,bty="n")
dev.off()

## Transition probabilities From State A at time t1 to State B at time t2
tval <- seq(0,20,length=100)
uu <- 10
p01f <- p01(aa,t1=uu,t2=uu+tval)
p12f <- p12(aa,t1=uu,t2=uu+tval)
p23f <- p23(aa,t1=uu,t2=uu+tval)
p13f <- p13(aa,t1=uu,t2=uu+tval)
p03f <- p03(aa,t1=uu,t2=uu+tval)
plot(tval,p01f,type="l",ylim=c(0,1))
lines(tval,p12f,col="red")
lines(tval,p23f,col="green")
lines(tval,p13f,col="blue")
lines(tval,p03f,col="magenta")

## Plot fitted hazard functions
tval <- seq(0,20,length=500)
aa2 <- c(1,1/(0.08118),1,1/(0.33),1,1/(0.2889),1,1/(0.0445),1,1/0.06349)

pdf("cav_weibull_hazard.pdf", width=11, height=7)
par(mfrow=c(2,2))
par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1.2)
plot(tval,h_01(tval,aa),type="l",ylim=c(0,0.8),col="#0571b0",lwd=3,xlab=" ",ylab="hazard",main="0 -> 1")
lines(tval,h_01(tval,aa2),col="#ca0020",lwd=3)
legend("topright",legend=c("Markov model","MCN model"),col=c("#ca0020","#0571b0"),lwd=3,bty="n")

plot(tval,h_12(tval,aa),type="l",ylim=c(0,0.8),col="#0571b0",lwd=3,xlab=" ",ylab="hazard",main="1 -> 2")
lines(tval,h_12(tval,aa2),col="#ca0020",lwd=3)
legend("topright",legend=c("Markov model","MCN model"),col=c("#ca0020","#0571b0"),lwd=3,bty="n")

plot(tval,h_23(tval,aa),type="l",ylim=c(0,0.8),col="#0571b0",lwd=3,xlab="years after transplantation",ylab="hazard",main="2 -> 3")
lines(tval,h_23(tval,aa2),col="#ca0020",lwd=3)
legend("topright",legend=c("Markov model","MCN model"),col=c("#ca0020","#0571b0"),lwd=3,bty="n")


plot(tval,h_03(tval,aa),type="l",ylim=c(0,0.8),col="#0571b0",lwd=3,xlab="years after transplantation",ylab="hazard",main="0 -> 3 and 1 -> 3")
lines(tval,h_03(tval,aa2),col="#ca0020",lwd=3)
lines(tval,h_13(tval,aa),lty=2,col="#0571b0",lwd=3)
lines(tval,h_13(tval,aa2),col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("Markov model, 0 -> 3","MCN model, 0 -> 3","Markov model, 1 -> 3","MCN model, 1 -> 3"),
       col=c("#ca0020","#0571b0","#ca0020","#0571b0"),lwd=3,bty="n",lty=c(1,1,2,2))

plot(tval,h_13(tval,aa),type="l",ylim=c(0,0.8),col="#0571b0",lwd=3,xlab="years after transplantation",ylab="hazard",main="1 -> 3")
lines(tval,h_13(tval,aa2),col="#ca0020",lwd=3)
legend("topright",legend=c("Markov model","MCN model"),col=c("#ca0020","#0571b0"),lwd=3,bty="n")
dev.off()

## Plot fitted overall survival function against KM
tval <- seq(0,20,length=100)
pdf("cav_weibull_survival.pdf", width=11, height=7)
par(mfrow=c(1,1))
par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1.5)
plot.survfit.msm(cav.msm, col.surv="black",lwd.surv=2,xlab="years after transplantation",
                 ylab="survival",main=" ",legend.pos=None,col="#ca0020",lwd=3)
lines(tval,S_total(tval,aa),col="#0571b0",lwd=3)
legend("topright",legend=c("Kaplan-Meier","Markov model","MCN model"),col=c("black","#ca0020","#0571b0"),lwd=3,bty="n")
dev.off()


aa2 <- c(1,1/(0.08118),1,1/(0.33),1,1/(0.2889),1,1/(0.0445),1,1/0.06349)

plot.survfit.msm(cav.msm, col.surv="black",lwd.surv=2,xlab="years after transplantation",
                 ylab="survival",main=" ",legend.pos=None,col="#ca0020",lwd=3)
lines(tval,S_total(tval,aa2),col="green",lty=2)
lines(tval,S_total(tval,aa),col="#0571b0",lwd=3)

