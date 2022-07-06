#### CAV example for testing
### Generalised Weibull model with covariates (to see that we get the same as Jackson)
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

# Density and survival function
pXweibull <- function(tt,a,b,th){
  pp <- 1-exp(1-(1+(tt/b)^a)^(1/th))
  pp[tt<0] <- 0
  return(pp)
}
dXweibull <- function(tt,a,b,th){
  (1/th)*(a/b)*(tt/b)^(a-1)*(1+(tt/b)^a)^(1/th-1)*exp(1-(1+(tt/b)^a)^(1/th))
}

# Model:
S_01 = function(param, x, t){(1-pXweibull(t,exp(param[1]),exp(param[2]+param[3]*x[1]+param[4]*x[2]),exp(param[19])))}
S_12 = function(param, x, t){(1-pXweibull(t,exp(param[5]),exp(param[6]+param[7]*x[1]+param[8]*x[2]),exp(param[20])))}
S_23 = function(param, x, t){(1-pXweibull(t,exp(param[9]),exp(param[10]+param[11]*x[1]+param[12]*x[2]),exp(param[22])))}
S_03 = function(param, x, t){(1-pXweibull(t,exp(param[13]),exp(param[14]+param[15]*x[1]),exp(param[22])))}
S_13 = function(param, x, t){(1-pXweibull(t,exp(param[16]),exp(param[17]+param[18]*x[1]),exp(param[23])))}

f_01 = function(param, x, t){as.numeric(t>=0)*dXweibull(t,exp(param[1]),exp(param[2]+param[3]*x[1]+param[4]*x[2]),exp(param[19]))}
f_12 = function(param, x, t){as.numeric(t>=0)*dXweibull(t,exp(param[5]),exp(param[6]+param[7]*x[1]+param[8]*x[2]),exp(param[20]))}
f_23 = function(param, x, t){as.numeric(t>=0)*dXweibull(t,exp(param[9]),exp(param[10]+param[11]*x[1]+param[12]*x[2]),exp(param[21]))}
f_03 = function(param, x, t){as.numeric(t>=0)*dXweibull(t,exp(param[13]),exp(param[14]+param[15]*x[1]),exp(param[22]))}
f_13 = function(param, x, t){as.numeric(t>=0)*dXweibull(t,exp(param[16]),exp(param[17]+param[18]*x[1]),exp(param[23]))}


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

params <- log(c(1.4,3,0.21,1.65,1.2,3,0.8,1.2,
                1.1,3.7,0.9,0.6,
                0.4,4,1.6,1.3,
                0.71,3,0.34,3,rep(1,5)))

params <- c(0.34784701,2.52621502,-0.14892845,-0.29490633,-0.05209313,4.84311946,
            0.19968936,-0.28632690,-0.49958056,6.64122178,0.12169109,0.40170468,
            -0.54523069,0.66891295,-1.73124666,-2.04755639,4.59824638,
            4.75956540,-0.25620822,-3.18974878,-3.18017775,2.51101361,
            0.90323822)
mloglikelihood(params,integrand,all_integral_limits,method1 = "hcubature",X=X_data_set,mc_cores=1)

system.time({
  oo <- nlminb(params,mloglikelihood,integrand = integrand,limits = all_integral_limits,
               mc_cores=10,X=X_data_set,lower=rep(-50,length(params)),upper=rep(50,length(params)))
})
2*oo$objective
# 2725.521

#> round(exp(oo$par),3)
#[1]   1.416  12.506   0.862   0.745   0.949 126.864   1.221   0.751   0.607
#[10] 766.030   1.129   1.494   0.580   1.952   0.177   0.443   0.129  99.310
#[19] 116.695   0.194   0.774   0.041   0.042  12.317   2.468

hessian = pracma::hessian(mloglikelihood, params, integrand = integrand,limits = all_integral_limits,
                                      mc_cores=10,X=X_data_set)

aa <- c(0.34784701,2.52621502,-0.14892845,-0.29490633,-0.05209313,4.84311946,
        0.19968936,-0.28632690,-0.49958056,6.64122178,0.12169109,0.40170468,
        -0.54523069,0.66891295,-1.73124666,-0.81331241,-2.04755639,4.59824638,
        4.75956540,-0.25620822,-3.18974878,-3.18017775,2.51101361,
        0.90323822)
        


######################## Plots #################################

## Hazard functions
h_01 = function(tt,x,param){f_01(param,x,tt)/S_01(param,x,tt)}
h_12 = function(tt,x,param){f_12(param,x,tt)/S_12(param,x,tt)}
h_23 = function(tt,x,param){f_23(param,x,tt)/S_23(param,x,tt)}
h_03 = function(tt,x,param){f_03(param,x,tt)/S_03(param,x,tt)}
h_13 = function(tt,x,param){f_13(param,x,tt)/S_13(param,x,tt)}


## Plot fitted hazard functions
tval <- seq(0,20,length=500)

pdf("cav_Xweibull_cov_hazard.pdf", width=11, height=7)
par(mfrow=c(2,3))
par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1)
plot(tval,h_01(tval,c(-1,0),aa),type="l",ylim=c(0,1),col="#0571b0",lwd=3,xlab=" ",ylab="hazard",main="0 -> 1")
lines(tval,h_01(tval,c(-1,1),aa),col="#ca0020",lwd=3)
lines(tval,h_01(tval,c(1,0),aa),col="#0571b0",lwd=3,lty=2)
lines(tval,h_01(tval,c(1,1),aa),col="#ca0020",lwd=3,lty=2)
#legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
#       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)

plot(tval,h_12(tval,c(-1,0),aa),type="l",ylim=c(0,3.5),col="#0571b0",lwd=3,xlab=" ",ylab="hazard",main="1 -> 2")
lines(tval,h_12(tval,c(-1,1),aa),col="#ca0020",lwd=3)
lines(tval,h_12(tval,c(1,0),aa),col="#0571b0",lwd=3,lty=2)
lines(tval,h_12(tval,c(1,1),aa),col="#ca0020",lwd=3,lty=2)
#legend("bottomright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
#       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)

plot(tval,h_23(tval,c(-1,0),aa),type="l",ylim=c(0,2),col="#0571b0",lwd=3,xlab="years after entering state",ylab="hazard",main="2 -> 3")
lines(tval,h_23(tval,c(-1,1),aa),col="#ca0020",lwd=3)
lines(tval,h_23(tval,c(1,0),aa),col="#0571b0",lwd=3,lty=2)
lines(tval,h_23(tval,c(1,1),aa),col="#ca0020",lwd=3,lty=2)
#legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
#       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)


plot(tval,h_03(tval,c(-1,0),aa),type="l",ylim=c(0,0.8),col="#0571b0",lwd=3,xlab="years after entering state",ylab="hazard",main="0 -> 3")
lines(tval,h_03(tval,c(-1,1),aa),col="#ca0020",lwd=3)
lines(tval,h_03(tval,c(1,0),aa),col="#0571b0",lwd=3,lty=2)
lines(tval,h_03(tval,c(1,1),aa),col="#ca0020",lwd=3,lty=2)
#legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
#       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)


plot(tval,h_13(tval,c(-1,0),aa),type="l",ylim=c(0,1),col="#0571b0",lwd=3,xlab="years after entering state",ylab="hazard",main="1 -> 3")
lines(tval,h_13(tval,c(-1,1),aa),col="#ca0020",lwd=3)
lines(tval,h_13(tval,c(1,0),aa),col="#0571b0",lwd=3,lty=2)
lines(tval,h_13(tval,c(1,1),aa),col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)

dev.off()



## Occupancy probabilities
tval <- seq(0.01,30,length=100)
p0_y0 <- occupancy_prob("0",tval,aa,gg,xval=c(-1,0))
p0_y1 <- occupancy_prob("0",tval,aa,gg,xval=c(-1,1))
p0_o0 <- occupancy_prob("0",tval,aa,gg,xval=c(1,0))
p0_o1 <- occupancy_prob("0",tval,aa,gg,xval=c(1,1))
p1_y0 <- occupancy_prob("1",tval,aa,gg,xval=c(-1,0))
p1_y1 <- occupancy_prob("1",tval,aa,gg,xval=c(-1,1))
p1_o0 <- occupancy_prob("1",tval,aa,gg,xval=c(1,0))
p1_o1 <- occupancy_prob("1",tval,aa,gg,xval=c(1,1))
p2_y0 <- occupancy_prob("2",tval,aa,gg,xval=c(-1,0))
p2_y1 <- occupancy_prob("2",tval,aa,gg,xval=c(-1,1))
p2_o0 <- occupancy_prob("2",tval,aa,gg,xval=c(1,0))
p2_o1 <- occupancy_prob("2",tval,aa,gg,xval=c(1,1))
p3_y0 <- occupancy_prob("3",tval,aa,gg,xval=c(-1,0))
p3_y1 <- occupancy_prob("3",tval,aa,gg,xval=c(-1,1))
p3_o0 <- occupancy_prob("3",tval,aa,gg,xval=c(1,0))
p3_o1 <- occupancy_prob("3",tval,aa,gg,xval=c(1,1))

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


pdf("cav_Xweibull_cov_prevalence.pdf", width=13, height=7)
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
#lines(tval,prev_y0$`Expected percentages`[,1],col="#0571b0",lwd=1)
# lines(tval,prev_y1$`Expected percentages`[,1],col="#ca0020",lwd=1)
# lines(tval,prev_o0$`Expected percentages`[,1],col="#0571b0",lwd=1,lty=2)
# lines(tval,prev_o1$`Expected percentages`[,1],col="#ca0020",lwd=1,lty=2)

plot(tval,prev_y0$`Observed percentages`[,2],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab=" ",
     ylab="prevalence (%)",main="State 1")
lines(tval,p1_y0*100,col="#0571b0",lwd=3)
lines(tval,p1_y1*100,col="#ca0020",lwd=3)
lines(tval,p1_o0*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p1_o1*100,col="#ca0020",lwd=3,lty=2)
# lines(tval,prev_y0$`Expected percentages`[,2],col="#0571b0",lwd=1)
# lines(tval,prev_y1$`Expected percentages`[,2],col="#ca0020",lwd=1)
# lines(tval,prev_o0$`Expected percentages`[,2],col="#0571b0",lwd=1,lty=2)
# lines(tval,prev_o1$`Expected percentages`[,2],col="#ca0020",lwd=1,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=1)

#legend("topright",legend=c("no IHD","IHD","older donor, no IHD","older donor, IHD"),
#       col=c("#0571b0","#ca0020","black","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.7)

plot(tval,prev_y0$`Observed percentages`[,3],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="State 2")
lines(tval,p2_y0*100,col="#0571b0",lwd=3)
lines(tval,p2_y1*100,col="#ca0020",lwd=3)
lines(tval,p2_o0*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p2_o1*100,col="#ca0020",lwd=3,lty=2)
# lines(tval,prev_y0$`Expected percentages`[,3],col="#0571b0",lwd=1)
# lines(tval,prev_y1$`Expected percentages`[,3],col="#ca0020",lwd=1)
# lines(tval,prev_o0$`Expected percentages`[,3],col="#0571b0",lwd=1,lty=2)
# lines(tval,prev_o1$`Expected percentages`[,3],col="#ca0020",lwd=1,lty=2)

plot(tval,prev_y0$`Observed percentages`[,4],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="State 3")
lines(tval,p3_y0*100,col="#0571b0",lwd=3)
lines(tval,p3_y1*100,col="#ca0020",lwd=3)
lines(tval,p3_o0*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p3_o1*100,col="#ca0020",lwd=3,lty=2)
# lines(tval,prev_y0$`Expected percentages`[,4],col="#0571b0",lwd=1)
# lines(tval,prev_y1$`Expected percentages`[,4],col="#ca0020",lwd=1)
# lines(tval,prev_o0$`Expected percentages`[,4],col="#0571b0",lwd=1,lty=2)
# lines(tval,prev_o1$`Expected percentages`[,4],col="#ca0020",lwd=1,lty=2)
dev.off()


#### Overall survival: 
tval <- seq(0.01,30,length=100)
Sy0 <- overall_survival(tval,aa,gg,c(-1,0))
Sy1 <- overall_survival(tval,aa,gg,c(-1,1))
So0 <- overall_survival(tval,aa,gg,c(1,0))
So1 <- overall_survival(tval,aa,gg,c(1,1))

plot.survfit.msm(cav.msm, col.surv="black",lwd.surv=2,xlab="years after transplantation",
                 ylab="survival",main=" ",legend.pos=None,col="grey",lwd=3,covariates=list(dage_st=-1,ihdTRUE=0))

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