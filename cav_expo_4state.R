### Fitting four-state illness death to the CAV data

source("expo_4state_illness_death.R")
#source("expo_4state_illness_death_marthes.R")

## CAV data
library(msm)
dd <- cav
dd <- dd[!is.na(dd$pdiag),]
dim(dd[dd$firstobs==1,]) # 614 patients
#first obs (at years=0) is always in state 1
id_wrong <- unique(dd$PTNUM[which(dd$state!=dd$statemax)])  # observations where the patient appears to go back to a previous state
dd <- dd[-which(dd$PTNUM %in% id_wrong),]

aa <- c(0.5,1,1.1,0.4,1.7)

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
mll2 <- function(param,timeMat){
  # timeMat: matrix of relevant time-points for each patient
  nn <- dim(timeMat)[1]
  lli <- rep(NA,nn)
  for (i in 1:nn){
    otype <- idtypes(timepoints[i,])
    if (otype==1){
      lli[i] <- mtype_1(param,timepoints[i,])
    }else if(otype==2){
      lli[i] <- mtype_2(param,timepoints[i,])
    }else if(otype==3){
      lli[i] <- mtype_3(param,timepoints[i,])
    }else if(otype==4){
      lli[i] <- mtype_4(param,timepoints[i,])
    }else if(otype==5){
      lli[i] <- mtype_5(param,timepoints[i,])
    }else if(otype==6){
      lli[i] <- mtype_6(param,timepoints[i,])
    }else if(otype==7){
      lli[i] <- mtype_7(param,timepoints[i,])
    }else if(otype==8){
      lli[i] <- mtype_8(param,timepoints[i,])
    }
  }
  return(sum(lli))
}

mll2_2 <- function(param,timeMat,otypes){
  # timeMat: matrix of relevant time-points for each patient
  # otypes: vector with the observation type of each patient
  tt1 <- timeMat[otypes==1,]
  tt2 <- timeMat[otypes==2,]
  tt3 <- timeMat[otypes==3,]
  tt4 <- timeMat[otypes==4,]
  tt5 <- timeMat[otypes==5,]
  tt6 <- timeMat[otypes==6,]
  tt7 <- timeMat[otypes==7,]
  tt8 <- timeMat[otypes==8,]
  lli <- rep(NA,8)
  lli[1] <- mtype_1(param,tt1)
  lli[2] <- mtype_2(param,tt2)
  lli[3] <- mtype_3(param,tt3)
  lli[4] <- mtype_4(param,tt4)
  lli[5] <- mtype_5(param,tt5)
  lli[6] <- mtype_6(param,tt6)
  lli[7] <- mtype_7(param,tt7)
  lli[8] <- mtype_8(param,tt8)
  return(sum(lli))
}
#timepoints[which(timepoints[,1]==0),1] <- 0.0000001

mll(aa,timepoints)
mll(c(0.08,0.33,0.29,0.04,0.06),timepoints)
mll2(c(0.08,0.33,0.29,0.04,0.06),timepoints)
mll2_2(c(0.08,0.33,0.29,0.04,0.06),timepoints,otypes)
mll(rep(10^(-3),5),timepoints)
mll(rep(10^(-2),5),timepoints)

system.time(mll(c(0.08,0.33,0.29,0.04,0.06),timepoints)) 
system.time(mll2(c(0.08,0.33,0.29,0.04,0.06),timepoints)) 
system.time(mll2_2(c(0.08,0.33,0.29,0.04,0.06),timepoints,otypes))

system.time(oo <- nlminb(aa,mll,timeMat=timepoints,lower=rep(0,5)) ) #Time:  61.681 
system.time(oo <- nlminb(aa,mll2,timeMat=timepoints,lower=rep(0,5)) ) #Time: 65.453
system.time(oo <- nlminb(aa,mll2_2,timeMat=timepoints,otypes=otypes,lower=rep(0,5)) ) #Time: 65.117
# Marthes actually seems a bit slower

oo <- nlminb(c(0.08,0.33,0.29,0.04,0.06),mll,timeMat=timepoints,lower=rep(10^(-3),5))

oo <- optim(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1), mll,timeMat=timepoints, method = "L-BFGS-B", hessian = FALSE
            ,lower=rep(10^(-5),5))

# Make plots of transition probabilities and prevalence (state occupancy)
aa <- c(0.0812,0.330,0.289,0.0445,0.0635)
tval <- seq(0,20,length=100)
plot(tval,pi0(aa,tval),type="l",ylim=c(0,1))
lines(tval,pi1(aa,tval),col="red")
lines(tval,pi2(aa,tval),col="green")
lines(tval,pi3(aa,tval),col="blue")

plot(tval,pi0(aa,tval)+pi1(aa,tval)+pi2(aa,tval)+pi3(aa,tval),type="l")

# Homogenous markov
twoway4.q <- rbind(c(0, 0.25, 0, 0.25), c(0.166, 0, 0.166, 0.166),c(0, 0.25, 0, 0.25), c(0, 0, 0, 0))
rownames(twoway4.q) <- colnames(twoway4.q) <- c("Well", "Mild","Severe", "Death")
cav.msm <- msm(state ~ years, subject = PTNUM, data = dd,qmatrix = twoway4.q, death = 4)

prev <- prevalence.msm(cav.msm,times=tval)

plot(tval,prev$`Observed percentages`[,1],type="l",col="blue",ylim=c(0,100))
lines(tval,prev$`Expected percentages`[,1],col="red")
lines(tval,pi0(aa,tval)*100,col="green",lty=2)

plot(tval,prev$`Observed percentages`[,2],type="l",col="blue",ylim=c(0,100))
lines(tval,prev$`Expected percentages`[,2],col="red")
lines(tval,pi1(aa,tval)*100,col="green",lty=2)

plot(tval,prev$`Observed percentages`[,3],type="l",col="blue",ylim=c(0,100))
lines(tval,prev$`Expected percentages`[,3],col="red")
lines(tval,pi2(aa,tval)*100,col="green",lty=2)

plot(tval,prev$`Observed percentages`[,4],type="l",col="blue",ylim=c(0,100))
lines(tval,prev$`Expected percentages`[,4],col="red")
lines(tval,pi3(aa,tval)*100,col="green",lty=2)

# Transition probabilities
aa <- c(0.0812,0.330,0.289,0.0445,0.0635)
tval <- seq(0,20,length=100)
uu <- 20
plot(tval,p01(aa,t1=uu,t2=uu+tval),type="l",ylim=c(0,1))
lines(tval,p12(aa,t1=uu,t2=uu+tval),col="red")
lines(tval,p23(aa,t1=uu,t2=uu+tval),col="green")
lines(tval,p13(aa,t1=uu,t2=uu+tval),col="blue")
lines(tval,p03(aa,t1=uu,t2=uu+tval),col="magenta")

plot(tval,p01(aa,t1=uu,t2=uu+tval)+p03(aa,t1=uu,t2=uu+tval),type="l",ylim=c(0,1))

plot(tval,p12(aa,t1=uu,t2=uu+tval)+p13(aa,t1=uu,t2=uu+tval),type="l",ylim=c(0,1))



############# old code

inds <- unique(ddr$PTNUM)
nn <- length(inds)
otypes <- rep(NA,nn)

for (i in 1:nn){
  ddi <- ddr[which(ddr$PTNUM==inds[i]),]
  otypes[i] <- idtypes(timepoints[i,])
}

####

inds <- unique(ddr$PTNUM)
nn <- length(inds)
lli <- rep(NA,nn)

for (i in 1:nn){
  ddi <- ddr[which(ddr$PTNUM==inds[i]),]
  otype <- idtypes(timepoints[i,])
  if (otype==1){
    lli[i] <- otype_1(aa,timepoints[i,])
  }else if(otype==2){
    lli[i] <- otype_2(aa,timepoints[i,])
  }else if(otype==3){
    lli[i] <- otype_3(aa,timepoints[i,])
  }else if(otype==4){
    lli[i] <- otype_4(aa,timepoints[i,])
  }else if(otype==5){
    lli[i] <- otype_5(aa,timepoints[i,])
  }else if(otype==6){
    lli[i] <- otype_6(aa,timepoints[i,])
  }else if(otype==7){
    lli[i] <- otype_7(aa,timepoints[i,])
  }else if(otype==8){
    lli[i] <- otype_8(aa,timepoints[i,])
  }
}

## Investigating some likelihood contributions
pp1 <- dd[dd$PTNUM=="100021",] #stays for a while in state 1 and then jumps to 4 (so can be type 03 or 03(1) or 03(12) - 16, 21 or 5)
pp2 <- dd[dd$PTNUM=="100067",] #straight to state 4 (so can be type 3(0) or 3(1) or 3(1,2) - 17, 20 or 14)

t_i16 <- pp1$years[3]
t_i21 <- t_i5 <- matrix(c(pp1$years[2:3]),1,2)
type_16(aa)
type_21(aa)
type_5(aa)

t_i17 <- t_i20 <- t_i14 <- matrix(c(pp2$years[2]),1,1)
type_17(aa)
type_20(aa)
type_14(aa)

t_i16 <- pp2$years[2]
t_i21 <- t_i5 <- matrix(c(pp2$years[1:2]),1,2)
type_16(aa)
type_21(aa)
type_5(aa)


idtypes <- function(obs_states,times){ # Identifies one of eight observation types
  # obs_states: a vector indicating which states have been observed, f.eks. c(1,3,3,4)
  # times: a vector of same lenght with the relevant times (after filtrating out the unrelevant ones). At most 6 values
  states <- rle(obs_states)$values
  kk <- length(states) 
  mm <- length(obs_states)
  otype <- NA
  if (kk == 1){
    if (sum(states==1)==kk){ #Type "0"
      lims <- c(times["t0M"],rep(NA,9))
      otype <- 1 #list("0" = lims)
    }
  }else if (kk == 2){
    if (sum(states==c(1,2))==kk){ #Type "01"
      lims <- c(times["t0M"],times["t1m"],NA,times["t1M"],rep(NA,3))
      otype <- 2 #list("01" = lims)
    }else if (sum(states==c(1,3))==kk){ #Type "0(1)2"
      lims <- c(times["t0M"],times["t2m"],NA,0,times["t2m"],NA,times["t2M"])
      otype <- 3 #list("012" = lims)
    }else if (sum(states==c(1,4))==kk){ #Types "0(12)3", "0(1)3" or "03"
      lims1 <- c(times["t0M"],times["t3m"],NA,0,times["t3m"],NA,times["t3m"])
      lims2 <- c(times["t0M"],times["t3m"],NA,NA,NA,times["t3m"],NA)
      lims3 <- c(NA,NA,times["t3m"],rep(NA,4))
      otype <- 4 #list("0123" = lims1,"013" = lims2,"03" = lims3)
    }
  }else if (kk == 3){
    if (sum(states==c(1,2,3))==kk){  #Type "012"
      lims <- c(times["t0M"],times["t1m"],NA,times["t1M"],times["t2m"],NA,times["t2M"])
      otype <- 5 #list("012" = lims)
    }else if (sum(states==c(1,2,4))==kk){ #Types "01(2)3" or "013"
      lims1 <- c(times["t0M"],times["t1m"],NA,times["t1M"],times["t3m"],NA,times["t3m"])
      lims2 <- c(times["t0M"],times["t1m"],NA,NA,NA,times["t3m"],NA)
      otype <- 6 #list("0123"=lims1,"013"=lims2)
    }else if (sum(states==c(1,3,4))==kk){ #Type "0(1)23"
      lims <- c(times["t0M"],times["t2m"],NA,0,times["t2m"],NA,times["t3m"])
      otype <- 7 #list("0123" = lims)
    }
  }else if (kk == 4){ #Type "0123"
    if (sum(states==c(1,2,3,4))==kk){
      lims <- c(times["t0M"],times["t1m"],NA,times["t1M"],times["t2m"],NA,times["t3m"])
      otype <- 8 #list("0123" = lims)
    }
  }
  return(otype)
}
