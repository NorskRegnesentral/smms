## Four-state illness-death, exact time of transition to absorbing state known
## Expo transition times

## Survival functions
S_01 = function(tt,param){1-pexp(tt,param[1])}
S_12 = function(tt,param){1-pexp(tt,param[2])}
S_23 = function(tt,param){1-pexp(tt,param[3])}
S_03 = function(tt,param){1-pexp(tt,param[4])}
S_13 = function(tt,param){1-pexp(tt,param[5])}

## Density functions 
f_01 = function(tt,param){dexp(tt,param[1])}
f_12 = function(tt,param){dexp(tt,param[2])}
f_23 = function(tt,param){dexp(tt,param[3])}
f_03 = function(tt,param){dexp(tt,param[4])}
f_13 = function(tt,param){dexp(tt,param[5])}

## Hazard functions
h_01 = function(tt,param){f_01(tt,param)/S_01(tt,param)}
h_12 = function(tt,param){f_12(tt,param)/S_12(tt,param)}
h_23 = function(tt,param){f_23(tt,param)/S_23(tt,param)}
h_03 = function(tt,param){f_03(tt,param)/S_03(tt,param)}
h_13 = function(tt,param){f_13(tt,param)/S_13(tt,param)}

# Inner functions
f01_S12_S03_S13 = function(ss,teval,param,teval2=NA){
  f_01(ss,param)*S_03(ss,param)*S_12(teval-ss,param)*S_13(teval-ss,param)
}
f01_f12_S23_S03_S13 = function(uu,ss,teval,param,teval2=NA){
    f_01(ss,param)*f_12(uu,param)*S_23(teval-uu-ss,param)*S_03(ss,param)*S_13(uu,param)
}
f03_S01 <- function(ss,teval,param,teval2=NA){
  f_03(ss,param)*S_01(ss,param)
}
f01_f13_S12_S03 = function(uu,ss,teval,param,teval2=NA){
  f_01(ss,param)*S_03(ss,param)*f_13(uu,param)*S_12(uu,param)
}
f01_f12_S23_S23_S03_S13 = function(uu,ss,teval,teval2=NA,param){
  if (is.na(teval2)){
    up <- 1
  }else{
    up <- S_23(teval2-uu-ss,param)
  }
  f_01(ss,param)*f_12(uu,param)*(up-S_23(teval-uu-ss,param))*S_03(ss,param)*S_13(uu,param)
}
f01_f12_f23_S03_S13 = function(uu,ss,teval,param,teval2=NA){ # if T23 is observed
  f_01(ss,param)*f_12(uu,param)*f_23(teval-uu-ss,param)*S_03(ss,param)*S_13(uu,param)
}
f01_f13_S12_S03_exact = function(ss,teval,param,teval2=NA){ # if T13 is observed
  f_01(ss,param)*S_03(ss,param)*f_13(teval-ss,param)*S_12(teval-ss,param)
}

# Integrals over functions of 3 variables
int3 <- function(uu,ss,innerfunc,teval,param,tlim1_l3,tlim2_l3){ #integrate over rr
  mm <- length(uu)
  out <- rep(NA,mm)
  for (i in 1:mm){
    out[i] <- integrate(innerfunc,lower=max(tlim1_l3-uu[i]-ss,0),upper=tlim2_l3-uu[i]-ss,teval=teval,uu=uu[i], 
                        ss=ss,param=param)$value
  }
  return(out)
}
# Integrals over functions of 2 variables
int2 <- function(ss,innerfunc,teval,param,tlim1_l2,tlim2_l2,teval2=NA){ #integrate over uu
  mm <- length(ss)
  out <- rep(NA,mm)
  for (i in 1:mm){
    out[i] <- integrate(innerfunc,lower=max(tlim1_l2-ss[i],0),upper=tlim2_l2-ss[i],teval=teval,teval2=teval2,
                          ss=ss[i],param=param)$value
  }
  return(out)
} 
# Integral over functions of 1 variable
int1 <- function(innerfunc,teval,param,tlim1,tlim2,tlim1_l2,tlim2_l2,teval2=NA){ #integrate over ss
  if (is.na(tlim1_l2)){
    out <- integrate(innerfunc,lower=tlim1,upper=tlim2,teval=teval,param=param,teval2=teval2)$value
  }else{
    out <- integrate(int2,innerfunc=innerfunc,lower=max(tlim1,0),upper=tlim2,teval=teval,param=param,tlim1_l2=tlim1_l2,tlim2_l2=tlim2_l2,
                     teval2=teval2)$value
  }
  return(out)
} 


## Obs type 1 ("0")
otype_1 = function(param,times){-log(S_01(times["t0M"],param)) - log(S_03(times["t0M"],param))}

## Obs type 2 ("01")
otype_2 = function(param,times){
  -log(int1(f01_S12_S03_S13,teval=times["t1M"],param=param,tlim1=times["t0M"],tlim2=times["t1m"],tlim1_l2=NA,tlim2_l2=NA))
}

## Obs type 3 ("0(1)2")
otype_3 = function(param,times){
  -log(int1(f01_f12_S23_S03_S13,teval=times["t2M"],param=param,tlim1=times["t0M"],tlim2=times["t2m"],tlim1_l2=times["t0M"],tlim2_l2=times["t2m"]))
}

## Obs type 4 (0(12)3", "0(1)3" or "03") 
otype_4 = function(param,times){
  -log(int1(f01_f12_f23_S03_S13,teval=times["t3m"],param=param,tlim1=times["t0M"],tlim2=times["t3m"],tlim1_l2=times["t0M"],tlim2_l2=times["t3m"])
    + int1(f01_f13_S12_S03_exact,teval=times["t3m"],param=param,tlim1=times["t0M"],tlim2=times["t3m"],tlim1_l2=NA,tlim2_l2=NA)
    + f_03(times["t3m"],param)*S_01(times["t3m"],param))
}

## Obs type 5 ("012")
otype_5 = function(param,times){
  -log(int1(f01_f12_S23_S03_S13,teval=times["t2M"],param=param,tlim1=times["t0M"],tlim2=times["t1m"],tlim1_l2=times["t1M"],tlim2_l2=times["t2m"]))
}

## Obs type 6 (01(2)3"or "013") 
otype_6 = function(param,times){
  -log(int1(f01_f12_f23_S03_S13,teval=times["t3m"],param=param,tlim1=times["t0M"],tlim2=times["t1m"],tlim1_l2=times["t1M"],tlim2_l2=times["t3m"])
       + int1(f01_f13_S12_S03_exact,teval=times["t3m"],param=param,tlim1=times["t0M"],tlim2=times["t1m"],tlim1_l2=NA,tlim2_l2=NA))
}

## Obs type 7 ("0(1)23") 
otype_7 = function(param,times){
  -log(int1(f01_f12_f23_S03_S13,teval=times["t3m"],param=param,tlim1=times["t0M"],tlim2=times["t2m"],tlim1_l2=times["t0M"],tlim2_l2=times["t2m"]))
}

## Obs type 8 ("01234")
otype_8 = function(param,times){
  -log(int1(f01_f12_f23_S03_S13,teval=times["t3m"],param=param,tlim1=times["t0M"],tlim2=times["t1m"],tlim1_l2=times["t1M"],tlim2_l2=times["t2m"]))
}

## Occupancy probabilities
pi0 = function(param,tt){S_01(tt,param)*S_03(tt,param)}

pi1 = function(param,tt){
  nn <- length(tt)
  pis <- rep(NA,nn)
  for (i in 1:nn){
    pis[i] <- integrate(f01_S12_S03_S13, lower = 0, upper = tt[i], teval = tt[i], param = param)$value
  }
  return(pis)
}

pi2 = function(param,tt){
  nn <- length(tt)
  pis <- rep(NA,nn)
  for (i in 1:nn){
    pis[i] <- int1(f01_f12_S23_S03_S13,teval=tt[i],param=param,tlim1=0,tlim2=tt[i],tlim1_l2=0,tlim2_l2=tt[i])
  }
  return(pis)
}

pi3 = function(param,tt){
  nn <- length(tt)
  pis <- rep(NA,nn)
  for (i in 1:nn){
    pis[i] <- (int1(f01_f12_S23_S23_S03_S13,teval=tt[i],param=param,tlim1=0,tlim2=tt[i],tlim1_l2=0,tlim2_l2=tt[i])
               + int1(f01_f13_S12_S03,teval=NA,param=param,tlim1=0,tlim2=tt[i],tlim1_l2=0,tlim2_l2=tt[i])
               + int1(f03_S01,teval=NA,param=param,tlim1=0,tlim2=tt[i],tlim1_l2=NA,tlim2_l2=NA))
  }
  return(pis)
}


## Transition probabilities From State A at time t1 to State B at time t2
p01 = function(param,t1,t2){
  nn <- length(t2)
  pp <- rep(NA,nn)
  for (i in 1:nn){
    pp[i] <- int1(f01_S12_S03_S13,teval=t2[i],param=param,tlim1=t1,tlim2=t2[i],tlim1_l2=NA,tlim2_l2=NA)
  }
  dn <- pi0(param=param,t1)
  return(pp/dn)
}
p12 = function(param,t1,t2){
  nn <- length(t2)
  pp <- rep(NA,nn)
  for (i in 1:nn){
    pp[i] <- int1(f01_f12_S23_S03_S13,teval=t2[i],param=param,tlim1=0,tlim2=t1,tlim1_l2=t1,tlim2_l2=t2[i])
  }
  dn <- pi1(param=param,t1)
  return(pp/dn)
}
p23 = function(param,t1,t2){
  nn <- length(t2)
  pp <- rep(NA,nn)
  for (i in 1:nn){
    pp[i] <- int1(f01_f12_S23_S23_S03_S13,teval=t2[i],teval2=t1,param=param,tlim1=0,tlim2=t1,tlim1_l2=0,tlim2_l2=t1)
  }
  dn <- pi2(param=param,t1)
  return(pp/dn)
}
p13 = function(param,t1,t2){
  nn <- length(t2)
  pp <- rep(NA,nn)
  for (i in 1:nn){
    pp[i] <- int1(f01_f13_S12_S03,teval=NA,param=param,tlim1=0,tlim2=t1,tlim1_l2=t1,tlim2_l2=t2[i])
  }
  dn <- pi1(param=param,t1)
  return(pp/dn)
}
p03 = function(param,t1,t2){
  nn <- length(t2)
  pp <- rep(NA,nn)
  for (i in 1:nn){
    pp[i] <- int1(f03_S01,teval=NA,param=param,tlim1=t1,tlim2=t2[i],tlim1_l2=NA,tlim2_l2=NA)
  }
  dn <- pi0(param=param,t1)
  return(pp/dn)
}

## Sub-survival functions
f01_S03 <- function(ss,param,teval2=NA){
  f_01(ss,param)*S_03(ss,param)
}

S01_sub <- function(tt,param){
  nn <- length(tt)
  pp <- rep(NA,nn)
  uu <- integrate(f01_S03,0,Inf,param=param)$value
  for (i in 1:nn){
    pp[i] <- integrate(f01_S03,tt[i],Inf,param=param)$value
  }
  return(pp)
  #exp(-(param[1]+param[4])*tt)
}

S03_sub <- function(tt,param){
  nn <- length(tt)
  pp <- rep(NA,nn)
  uu <- integrate(f03_S01,0,Inf,param=param)$value
  for (i in 1:nn){
    pp[i] <- integrate(f03_S01,tt[i],Inf,param=param)$value
  }
  return(pp)
  #exp(-(param[1]+param[4])*tt)
}

# tval <- seq(0,30,length=100)
# plot(tval,S01_sub(tval,aa),type="l",ylim=c(0,1))
# lines(tval,S03_sub(tval,aa),col="red")
# 
# plot(tval,S01_sub(tval,aa)+S03_sub(tval,aa),type="l",ylim=c(0,1))

## Overall survival probability (stay in state 0, move from 0 to 1 and stay in state 1, move from 0 to 2 and stay in state 2)
S_total <- function(tt,param){
  nn <- length(tt)
  pp <- rep(NA,nn)
  for (i in 1:nn){
    pp[i] <- (S_01(tt[i],param)*S_03(tt[i],param) + 
                int1(f01_S12_S03_S13,teval=tt[i],param=param,tlim1=0,tlim2=tt[i],tlim1_l2=NA,tlim2_l2=NA) +
                int1(f01_f12_S23_S03_S13,teval=tt[i],param=param,tlim1=0,tlim2=tt[i],tlim1_l2=0,tlim2_l2=tt[i])) 
  }
  return(pp)
}
