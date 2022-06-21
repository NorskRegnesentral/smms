#### CAV example for testing
### Expo model with covariates (to see that we get the same as Jackson)
setwd("~/Multistate models/smms")
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
S_01 = function(param, x, t){(1-pexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2])))}
S_12 = function(param, x, t){(1-pexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2])))}
S_23 = function(param, x, t){(1-pexp(t,exp(param[7]+param[8]*x[1]+param[9]*x[2])))}
S_03 = function(param, x, t){(1-pexp(t,exp(param[10]+param[11]*x[1]+param[12]*x[2])))}
S_13 = function(param, x, t){(1-pexp(t,exp(param[13]+param[14]*x[1]+param[15]*x[2])))}

f_01 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2]))}
f_12 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2]))}
f_23 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[7]+param[8]*x[1]+param[9]*x[2]))}
f_03 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[10]+param[11]*x[1]+param[12]*x[2]))}
f_13 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[13]+param[14]*x[1]+param[15]*x[2]))}


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

params <- c(log(0.08),0.21,log(1.65),log(0.33),log(0.8),log(1.2),log(0.29),log(0.9),log(0.6),
            log(0.044),log(1.6),log(1.3),log(0.06),log(0.34),log(3))

mloglikelihood(params,integrand,all_integral_limits,method1 = "hcubature",X=X_data_set,mc_cores=1)

system.time({
  oo <- nlminb(params,mloglikelihood,integrand = integrand,limits = all_integral_limits,
               mc_cores=10,X=X_data_set)
})
2*oo$objective
# 2821.21 we find exactly the same result as in Marthes thesis (with the MSM package)

aa <- oo$par

hessian = numDeriv::hessian(mloglikelihood, aa, integrand = integrand,limits = all_integral_limits,
                           mc_cores=10,X=X_data_set)
hessian2 = pracma::hessian(mloglikelihood, aa, integrand = integrand,limits = all_integral_limits,
                           mc_cores=10,X=X_data_set)
# numDeriv and pracma give same result here. But pracma is faster.
load("cav_expo_cov_optims")

est_ci(aa,hessian)

exp(est_ci(aa,hessian))

est_ci(aa,hessian,log=F)

## Occupancy probabilities
tval <- seq(0.01,30,length=50)
p0_y0_ci <- occupancy_prob_ci_band("0",tval,aa,gg,xval=c(-1,0),hessian)
p0_y1_ci <- occupancy_prob_ci_band("0",tval,aa,gg,xval=c(-1,1),hessian)
p0_o0_ci <- occupancy_prob_ci_band("0",tval,aa,gg,xval=c(1,0),hessian)
p0_o1_ci <- occupancy_prob_ci_band("0",tval,aa,gg,xval=c(1,1),hessian)
p1_y0_ci <- occupancy_prob_ci_band("1",tval,aa,gg,xval=c(-1,0),hessian)
p1_y1_ci <- occupancy_prob_ci_band("1",tval,aa,gg,xval=c(-1,1),hessian)
p1_o0_ci <- occupancy_prob_ci_band("1",tval,aa,gg,xval=c(1,0),hessian)
p1_o1_ci <- occupancy_prob_ci_band("1",tval,aa,gg,xval=c(1,1),hessian)
p2_y0_ci <- occupancy_prob_ci_band("2",tval,aa,gg,xval=c(-1,0),hessian)
p2_y1_ci <- occupancy_prob_ci_band("2",tval,aa,gg,xval=c(-1,1),hessian)
p2_o0_ci <- occupancy_prob_ci_band("2",tval,aa,gg,xval=c(1,0),hessian)
p2_o1_ci <- occupancy_prob_ci_band("2",tval,aa,gg,xval=c(1,1),hessian)
p3_y0_ci <- occupancy_prob_ci_band("3",tval,aa,gg,xval=c(-1,0),hessian)
p3_y1_ci <- occupancy_prob_ci_band("3",tval,aa,gg,xval=c(-1,1),hessian)
p3_o0_ci <- occupancy_prob_ci_band("3",tval,aa,gg,xval=c(1,0),hessian)
p3_o1_ci <- occupancy_prob_ci_band("3",tval,aa,gg,xval=c(1,1),hessian)


# msm package with covariates
oneway4.q <- rbind(c(0, 0.25, 0, 0.25), c(0, 0, 0.25, 0.25),c(0, 0, 0, 0.5), c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well", "Mild","Severe", "Death")
cav.msm <- msm(state ~ time, subject = patient, data = dd,covariates=~dage_st+ihd,
               qmatrix = oneway4.q, death = 4,method = "BFGS", control = list(fnscale = 4000, maxit = 10000))

prev_y0 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=-1,ihdTRUE=0))
prev_y1 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=-1,ihdTRUE=1))
prev_o0 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=1,ihdTRUE=0))
prev_o1 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=1,ihdTRUE=1))

png("cav_expo_cov_prevalence.png", width=700, height=400)
par(mfrow=c(2,2))
par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1)
plot(tval,prev_y0$`Observed percentages`[,1],type="l",ylim=c(0,100),lwd=3,xlab=" ",col="dark grey",
     ylab="prevalence (%)",main="State 0")
polygon(c(tval, rev(tval)), c(p0_y0_ci$upper*100, rev(p0_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p0_y1_ci$upper*100, rev(p0_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p0_o0_ci$upper*100, rev(p0_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p0_o1_ci$upper*100, rev(p0_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p0_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p0_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p0_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p0_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)


plot(tval,prev_y0$`Observed percentages`[,2],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab=" ",
     ylab="prevalence (%)",main="State 1")
polygon(c(tval, rev(tval)), c(p1_y0_ci$upper*100, rev(p1_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p1_y1_ci$upper*100, rev(p1_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p1_o0_ci$upper*100, rev(p1_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p1_o1_ci$upper*100, rev(p1_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p1_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p1_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p1_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p1_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=1)


plot(tval,prev_y0$`Observed percentages`[,3],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="State 2")
polygon(c(tval, rev(tval)), c(p2_y0_ci$upper*100, rev(p2_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p2_y1_ci$upper*100, rev(p2_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p2_o0_ci$upper*100, rev(p2_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p2_o1_ci$upper*100, rev(p2_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p2_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p2_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p2_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p2_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)

plot(tval,prev_y0$`Observed percentages`[,4],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="State 3")
polygon(c(tval, rev(tval)), c(p3_y0_ci$upper*100, rev(p3_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p3_y1_ci$upper*100, rev(p3_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p3_o0_ci$upper*100, rev(p3_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p3_o1_ci$upper*100, rev(p3_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p3_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p3_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p3_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p3_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)
dev.off()



#### Overall survival: 
tval <- seq(0.01,30,length=50)
Sy0 <- overall_survival_ci_band(tval,aa,gg,c(-1,0),hessian)
Sy1 <- overall_survival_ci_band(tval,aa,gg,c(-1,1),hessian)
So0 <- overall_survival_ci_band(tval,aa,gg,c(1,0),hessian)
So1 <- overall_survival_ci_band(tval,aa,gg,c(1,1),hessian)

png("cav_expo_cov_survival.png", width=700, height=400)
par(mfrow=c(1,1))
par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1.2)
plot.survfit.msm(cav.msm, col.surv="black",lwd.surv=2,xlab="years after transplantation",
                 ylab="survival",main=" ",legend.pos=None,col=NULL,lwd=3,covariates=list(dage_st=-1,ihdTRUE=0))
polygon(c(tval, rev(tval)), c(Sy0$upper, rev(Sy0$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(Sy1$upper, rev(Sy1$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(So0$upper, rev(So0$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(So1$upper, rev(So1$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,Sy0$est,col="#0571b0",lwd=3)
lines(tval,Sy1$est,col="#ca0020",lwd=3)
lines(tval,So0$est,col="#0571b0",lwd=3,lty=2)
lines(tval,So1$est,col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=0.9)
dev.off()



