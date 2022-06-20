#### CAV example for testing
### Expo model with covariates (to see that we get the same as Jackson)
setwd("H:/Multistate models/SemiMarkovMultistate")
rm(list = ls())
devtools::load_all() 
library(igraph)
library(parallel)
library(cubature)
## For test-data
library(msm)


## Preprocessing
gg = graph_from_literal("1"--+"2"--+"3"--+"4", "1"--+"4", "2"--+"4")

# Testing for CAV-data
dd = cav
dd = dd[!is.na(dd$pdiag),]
dim(dd[dd$firstobs==1,]) # 614 patients
#first obs (at years=0) is always in state 1
id_wrong = unique(dd$PTNUM[which(dd$state!=dd$statemax)])  # observations where the patient appears to go back to a previous state
dd = dd[-which(dd$PTNUM %in% id_wrong),]
## Only relevant parts 
dd = dd[ ,-c(2, 5, 7, 9, 10)]
dd$dage_st = (dd$dage-mean(dd$dage))/sd(dd$dage)
dd$ihd = (dd$pdiag=="IHD")
colnames(dd)[1:2] <- c("patient","time")

# Model:
S_01 = function(param, x, t){1-pweibull(t,exp(param[1]),exp(param[2]+param[3]*x[1]+param[4]*x[2]))}
S_12 = function(param, x, t){1-pweibull(t,exp(param[5]),exp(param[6]+param[7]*x[1]+param[8]*x[2]))}
S_23 = function(param, x, t){1-pweibull(t,exp(param[9]),exp(param[10]+param[11]*x[1]+param[12]*x[2]))}
S_03 = function(param, x, t){1-pweibull(t,exp(param[13]),exp(param[14]+param[15]*x[1]+param[16]*x[2]))}
S_13 = function(param, x, t){1-pweibull(t,exp(param[17]),exp(param[18]+param[19]*x[1]+param[20]*x[2]))}

f_01 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[1]),exp(param[2]+param[3]*x[1]+param[4]*x[2]))}
f_12 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[5]),exp(param[6]+param[7]*x[1]+param[8]*x[2]))}
f_23 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[9]),exp(param[10]+param[11]*x[1]+param[12]*x[2]))}
f_03 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[13]),exp(param[14]+param[15]*x[1]+param[16]*x[2]))}
f_13 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[17]),exp(param[18]+param[19]*x[1]+param[20]*x[2]))}


## Part 3: From the time points of a given patient to an integral
all_data_set = arrange_data(dd, gg)
X_data_set = aggregate(dd[,c("dage_st","ihd")],by=list(dd$patient),FUN=median)
X_data_set = as.matrix(X_data_set[,2:3])

## Finding all types, integrand, time points ++
formula_obs_types = all_types(gg)
edge_mats <- edge_matrices(gg)
state_ord = state_ordering(gg)
absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
names_surv_dens = names_of_survival_density(gg)

timepointMat <- all_data_set[,1:(dim(all_data_set)[2]-1)]

observation_type = rep(NA, nrow(all_data_set))
all_integral_limits = list()
integrand = list()

for(i in 1:nrow(all_data_set)){
  observation_type[i] = all_data_set[i,"obs_type"]
  f_types = names(which(formula_obs_types[, observation_type[i]] == 1))
  integrand_mellomregn = list()
  integral_mellomregn= list()
  for(j in 1:length(f_types)){
    integrand_mellomregn[[j]] = eval(parse(text=type_to_integrand(f_types[j], edge_mats, names_surv_dens)))
    integral_mellomregn[[j]] = finding_limits(timepointMat[i,],f_types[j],edge_mats, absorbing_states)
  }
  all_integral_limits[[i]] = integral_mellomregn
  integrand[[i]] = integrand_mellomregn
}

params <- c(0.41,2.3,-0.14,-0.327,0.25,1.11,0.17,-0.22,-0.01,0.92,0.09,0.43,-0.92,7.6,-1.59,-0.36,-1.08,4.88,1.89,-1.16)


mloglikelihood(params,integrand,all_integral_limits,method1 = "hcubature",X=X_data_set,mc_cores=10)

system.time({
  oo <- nlminb(params,mloglikelihood,integrand = integrand,limits = all_integral_limits,
               mc_cores=10,X=X_data_set)
})
2*oo$objective


# With covariates on the scale parameter
2*oo$objective #2736.213
## time: 3728 seconds (1 hour)
#oo$par
#[1]  0.395176238  2.306533314 -0.142396794 -0.326719502  0.252982108
#[6]  1.106276893  0.171426567 -0.218790051 -0.007948915  0.916087018
#[11]  0.092580630  0.426935854 -0.900624794  7.580816624 -1.589690902
#[16] -0.355544744 -1.066852225  4.882702792  1.891081007 -1.159004675

aa <- c(0.4,2.31,-0.14,-0.33,0.25,1.11,0.17,-0.22,-0.01,0.92,0.09,0.43,
        -0.9,7.58,-1.59,-0.36,-1.07,4.88,1.89,-1.16)


######################## Plots #################################

## Hazard functions
h_01 = function(tt,x,param){f_01(param,x,tt)/S_01(param,x,tt)}
h_12 = function(tt,x,param){f_12(param,x,tt)/S_12(param,x,tt)}
h_23 = function(tt,x,param){f_23(param,x,tt)/S_23(param,x,tt)}
h_03 = function(tt,x,param){f_03(param,x,tt)/S_03(param,x,tt)}
h_13 = function(tt,x,param){f_13(param,x,tt)/S_13(param,x,tt)}


## Plot fitted hazard functions
tval <- seq(0,20,length=500)

pdf("cav_weibull_cov_hazard.pdf", width=11, height=7)
par(mfrow=c(2,3))
par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1)
plot(tval,h_01(tval,c(-1,0),aa),type="l",ylim=c(0,1),col="#0571b0",lwd=3,xlab=" ",ylab="hazard",main="0 -> 1")
lines(tval,h_01(tval,c(-1,1),aa),col="#ca0020",lwd=3)
lines(tval,h_01(tval,c(1,0),aa),col="#0571b0",lwd=3,lty=2)
lines(tval,h_01(tval,c(1,1),aa),col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)

plot(tval,h_12(tval,c(-1,0),aa),type="l",ylim=c(0,1),col="#0571b0",lwd=3,xlab=" ",ylab="hazard",main="1 -> 2")
lines(tval,h_12(tval,c(-1,1),aa),col="#ca0020",lwd=3)
lines(tval,h_12(tval,c(1,0),aa),col="#0571b0",lwd=3,lty=2)
lines(tval,h_12(tval,c(1,1),aa),col="#ca0020",lwd=3,lty=2)
legend("bottomright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)

plot(tval,h_23(tval,c(-1,0),aa),type="l",ylim=c(0,1),col="#0571b0",lwd=3,xlab="years after entering state",ylab="hazard",main="2 -> 3")
lines(tval,h_23(tval,c(-1,1),aa),col="#ca0020",lwd=3)
lines(tval,h_23(tval,c(1,0),aa),col="#0571b0",lwd=3,lty=2)
lines(tval,h_23(tval,c(1,1),aa),col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)


plot(tval,h_03(tval,c(-1,0),aa),type="l",ylim=c(0,0.8),col="#0571b0",lwd=3,xlab="years after entering state",ylab="hazard",main="0 -> 3")
lines(tval,h_03(tval,c(-1,1),aa),col="#ca0020",lwd=3)
lines(tval,h_03(tval,c(1,0),aa),col="#0571b0",lwd=3,lty=2)
lines(tval,h_03(tval,c(1,1),aa),col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)


plot(tval,h_13(tval,c(-1,0),aa),type="l",ylim=c(0,1),col="#0571b0",lwd=3,xlab="years after entering state",ylab="hazard",main="1 -> 3")
lines(tval,h_13(tval,c(-1,1),aa),col="#ca0020",lwd=3)
lines(tval,h_13(tval,c(1,0),aa),col="#0571b0",lwd=3,lty=2)
lines(tval,h_13(tval,c(1,1),aa),col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)

dev.off()



## Occupancy probabilities

tval <- seq(0.01,30,length=100)
p0_y0 <- occupancy_prob("0",tval,params,gg,xval=c(-1,0))
p0_y1 <- occupancy_prob("0",tval,params,gg,xval=c(-1,1))
p0_o0 <- occupancy_prob("0",tval,params,gg,xval=c(1,0))
p0_o1 <- occupancy_prob("0",tval,params,gg,xval=c(1,1))
p1_y0 <- occupancy_prob("1",tval,params,gg,xval=c(-1,0))
p1_y1 <- occupancy_prob("1",tval,params,gg,xval=c(-1,1))
p1_o0 <- occupancy_prob("1",tval,params,gg,xval=c(1,0))
p1_o1 <- occupancy_prob("1",tval,params,gg,xval=c(1,1))
p2_y0 <- occupancy_prob("2",tval,params,gg,xval=c(-1,0))
p2_y1 <- occupancy_prob("2",tval,params,gg,xval=c(-1,1))
p2_o0 <- occupancy_prob("2",tval,params,gg,xval=c(1,0))
p2_o1 <- occupancy_prob("2",tval,params,gg,xval=c(1,1))
p3_y0 <- occupancy_prob("3",tval,params,gg,xval=c(-1,0))
p3_y1 <- occupancy_prob("3",tval,params,gg,xval=c(-1,1))
p3_o0 <- occupancy_prob("3",tval,params,gg,xval=c(1,0))
p3_o1 <- occupancy_prob("3",tval,params,gg,xval=c(1,1))

plot(tval,p0_y0+p1_y0+p2_y0+p3_y0,type="l")
plot(tval,p0_y1+p1_y1+p2_y1+p3_y1,type="l")
plot(tval,p0_o0+p1_o0+p2_o0+p3_o0,type="l")
plot(tval,p0_o1+p1_o1+p2_o1+p3_o1,type="l")

# msm package with covariates
twoway4.q <- rbind(c(0, 0.25, 0, 0.25), c(0.166, 0, 0.166, 0.166),c(0, 0.25, 0, 0.25), c(0, 0, 0, 0))
rownames(twoway4.q) <- colnames(twoway4.q) <- c("Well", "Mild","Severe", "Death")
cav.msm <- msm(state ~ time, subject = patient, data = dd,covariates=~dage_st+ihd,
               qmatrix = twoway4.q, death = 4,method = "BFGS", control = list(fnscale = 4000, maxit = 10000))

prev_y0 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=-1,ihdTRUE=0))
prev_y1 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=-1,ihdTRUE=1))
prev_o0 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=1,ihdTRUE=0))
prev_o1 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=1,ihdTRUE=1))


pdf("cav_weibull_cov_prevalence.pdf", width=11, height=7)
par(mfrow=c(2,2))
par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1)
plot(tval,prev_y0$`Observed percentages`[,1],type="l",ylim=c(0,100),lwd=3,xlab=" ",col="dark grey",
     ylab="prevalence (%)",main="State 0")
lines(tval,p0_y0*100,col="#0571b0",lwd=3)
lines(tval,p0_y1*100,col="#ca0020",lwd=3)
lines(tval,p0_o0*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p0_o1*100,col="#ca0020",lwd=3,lty=2)
lines(tval,prev_y0$`Expected percentages`[,1],col="#0571b0",lwd=1)
lines(tval,prev_y1$`Expected percentages`[,1],col="#ca0020",lwd=1)
lines(tval,prev_o0$`Expected percentages`[,1],col="#0571b0",lwd=1,lty=2)
lines(tval,prev_o1$`Expected percentages`[,1],col="#ca0020",lwd=1,lty=2)

plot(tval,prev_y0$`Observed percentages`[,2],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab=" ",
     ylab="prevalence (%)",main="State 1")
lines(tval,p1_y0*100,col="#0571b0",lwd=3)
lines(tval,p1_y1*100,col="#ca0020",lwd=3)
lines(tval,p1_o0*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p1_o1*100,col="#ca0020",lwd=3,lty=2)
lines(tval,prev_y0$`Expected percentages`[,2],col="#0571b0",lwd=1)
lines(tval,prev_y1$`Expected percentages`[,2],col="#ca0020",lwd=1)
lines(tval,prev_o0$`Expected percentages`[,2],col="#0571b0",lwd=1,lty=2)
lines(tval,prev_o1$`Expected percentages`[,2],col="#ca0020",lwd=1,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)

#legend("topright",legend=c("no IHD","IHD","older donor, no IHD","older donor, IHD"),
#       col=c("#0571b0","#ca0020","black","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)

plot(tval,prev_y0$`Observed percentages`[,3],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="State 2")
lines(tval,p2_y0*100,col="#0571b0",lwd=3)
lines(tval,p2_y1*100,col="#ca0020",lwd=3)
lines(tval,p2_o0*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p2_o1*100,col="#ca0020",lwd=3,lty=2)
lines(tval,prev_y0$`Expected percentages`[,3],col="#0571b0",lwd=1)
lines(tval,prev_y1$`Expected percentages`[,3],col="#ca0020",lwd=1)
lines(tval,prev_o0$`Expected percentages`[,3],col="#0571b0",lwd=1,lty=2)
lines(tval,prev_o1$`Expected percentages`[,3],col="#ca0020",lwd=1,lty=2)

plot(tval,prev_y0$`Observed percentages`[,4],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="State 3")
lines(tval,p3_y0*100,col="#0571b0",lwd=3)
lines(tval,p3_y1*100,col="#ca0020",lwd=3)
lines(tval,p3_o0*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p3_o1*100,col="#ca0020",lwd=3,lty=2)
lines(tval,prev_y0$`Expected percentages`[,4],col="#0571b0",lwd=1)
lines(tval,prev_y1$`Expected percentages`[,4],col="#ca0020",lwd=1)
lines(tval,prev_o0$`Expected percentages`[,4],col="#0571b0",lwd=1,lty=2)
lines(tval,prev_o1$`Expected percentages`[,4],col="#ca0020",lwd=1,lty=2)
dev.off()

#### Overall survival
S_total <- function(tt,param,x){
  nn <- length(tt)
  pp <- rep(NA,nn)
  for (i in 1:nn){
    pp[i] <- (S_01(param,x,tt[i])*S_03(param,x,tt[i]) + 
                int1(f01_S12_S03_S13,teval=tt[i],param=param,x=x,tlim1=0,tlim2=tt[i],tlim1_l2=NA,tlim2_l2=NA) +
                int1(f01_f12_S23_S03_S13,teval=tt[i],param=param,x=x,tlim1=0,tlim2=tt[i],tlim1_l2=0,tlim2_l2=tt[i])) 
  }
  return(pp)
}

tval <- seq(0,30,length=100)
Sy0 <- S_total(tval,aa,c(-1,0))
Sy1 <- S_total(tval,aa,c(-1,1))
So0 <- S_total(tval,aa,c(1,0))
So1 <- S_total(tval,aa,c(1,1))

plot.survfit.msm(cav.msm, col.surv="black",lwd.surv=2,xlab="years after transplantation",
                 ylab="survival",main=" ",legend.pos=None,col="#ca0020",lwd=3,covariates=list(dage_st=-1,ihdTRUE=0))

pdf("cav_weibull_cov_survival.pdf", width=11, height=7)
par(mfrow=c(1,1))
par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1.2)
plot(tval,Sy0,type="l",col="#0571b0",ylim=c(0,1),lwd=3,xlab="years after transplantation",
     ylab="survival")
lines(tval,Sy1,col="#ca0020",lwd=3)
lines(tval,So0,col="#0571b0",lwd=3,lty=2)
lines(tval,So1,col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.9)
dev.off()

