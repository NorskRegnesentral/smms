#### CAV dataset
### Expo model (i.e. time-homogeneous Markov) without covariates 
setwd("~/Multistate models/smms")
rm(list = ls())
devtools::load_all() 


## Preprocessing
gg = igraph::graph_from_literal("well"--+"mild"--+"severe"--+"death", "well"--+"death", "mild"--+"death")

# Testing for CAV-data
library(msm)
dd = cav
dd = dd[!is.na(dd$pdiag),]
dim(dd[dd$firstobs==1,]) # 614 patients
#first obs (at years=0) is always in state 1
id_wrong = unique(dd$PTNUM[which(dd$state!=dd$statemax)])  # observations where the patient appears to go back to a previous state
dd = dd[-which(dd$PTNUM %in% id_wrong),]
## Only relevant parts 
dd = dd[ ,-c(2, 5, 7, 9, 10)]
colnames(dd)[1:2] <- c("patient","time") # rename relevant columns (necessary in current version)

# Model:
names_of_survival_density(gg)

S_01 = function(param, x, t){1-pexp(t,exp(param[1]))}
S_12 = function(param, x, t){1-pexp(t,exp(param[2]))}
S_23 = function(param, x, t){1-pexp(t,exp(param[3]))}
S_03 = function(param, x, t){1-pexp(t,exp(param[4]))}
S_13 = function(param, x, t){1-pexp(t,exp(param[5]))}

f_01 = function(param, x, t){dexp(t,exp(param[1]))}
f_12 = function(param, x, t){dexp(t,exp(param[2]))}
f_23 = function(param, x, t){dexp(t,exp(param[3]))}
f_03 = function(param, x, t){dexp(t,exp(param[4]))}
f_13 = function(param, x, t){dexp(t,exp(param[5]))}


# Fitting the model 
startval <- c(-2.5,-1.1,-1.2,-3.1,-2.8)

mlo <- smms(startval,dd,gg, mc_cores = 1)
ho <- hessian_matrix(mlo$par,dd,gg,mc_cores=1) #calculates hessian (can be done inside smms() too)

# Compute AIC (higher values are better with this definition)
aic <- (-2*mlo$objective)-2*length(startval) #-2887.1

# Look at estimates and 95% confidence intervals
round(est_ci(mlo$par,ho),2) #parameters living on the real line

round(exp(est_ci(mlo$par,ho)),2) #parameters living on the positive halv-line (more interpretable for some parameters)

#### Occupancy probabilities (with confidence bands), for some specific covariate values
tval <- seq(0.01,30,length=50)
p0_ci <- occupancy_prob_ci_band("well",tval,mlo$par,gg,hessian=ho)
p1_ci <- occupancy_prob_ci_band("mild",tval,mlo$par,gg,hessian=ho)
p2_ci <- occupancy_prob_ci_band("severe",tval,mlo$par,gg,hessian=ho)
p3_ci <- occupancy_prob_ci_band("death",tval,mlo$par,gg,hessian=ho)


# msm package with covariates (to get non-parameteric prevalence curves)
oneway4.q <- rbind(c(0, 0.25, 0, 0.25), c(0, 0, 0.25, 0.25),c(0, 0, 0, 0.5), c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well", "Mild","Severe", "Death")
cav.msm <- msm(state ~ time, subject = patient, data = dd,qmatrix = oneway4.q, death = 4,method = "BFGS")

prev_np <- prevalence.msm(cav.msm,times=tval)

par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1)
plot(tval,prev_np$`Observed percentages`[,1],type="l",ylim=c(0,100),lwd=3,xlab=" ",col="dark grey",
     ylab="prevalence (%)",main="well")
polygon(c(tval, rev(tval)), c(p0_ci$upper*100, rev(p0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
lines(tval,p0_ci$est*100,col="#0571b0",lwd=3)

plot(tval,prev_np$`Observed percentages`[,2],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab=" ",
     ylab="prevalence (%)",main="mild")
polygon(c(tval, rev(tval)), c(p1_ci$upper*100, rev(p1_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
lines(tval,p1_ci$est*100,col="#0571b0",lwd=3)


plot(tval,prev_np$`Observed percentages`[,3],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="State 2")
polygon(c(tval, rev(tval)), c(p2_ci$upper*100, rev(p2_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
lines(tval,p2_ci$est*100,col="#0571b0",lwd=3)

plot(tval,prev_np$`Observed percentages`[,4],type="l",col="dark grey",ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="State 3")
polygon(c(tval, rev(tval)), c(p3_ci$upper*100, rev(p3_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
lines(tval,p3_ci$est*100,col="#0571b0",lwd=3)


#### Overall survival curves (with confidence bands), for some specific covariate values
tval <- seq(0.01,30,length=50)
So <- overall_survival_ci_band(tval,mlo$par,gg,hessian=ho)


par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1.2)
plot.survfit.msm(cav.msm, col.surv="black",lwd.surv=2,xlab="years after transplantation",
                 ylab="survival",main=" ",legend.pos=None,col=NULL,lwd=3)
# using msm package for the nonparametric estimate (for now)
polygon(c(tval, rev(tval)), c(So$upper, rev(So$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
lines(tval,So$est,col="#0571b0",lwd=3)



#### Transition probabilities (with confidence bands)
tval <- seq(10,40,length=50)
vt <- 10
t01_ci <- transition_prob_ci_band("well-mild",tval,vt,mlo$par,gg,hessian=ho)
t12_ci <- transition_prob_ci_band("mild-severe",tval,vt,mlo$par,gg,hessian=ho)
t23_ci <- transition_prob_ci_band("severe-death",tval,vt,mlo$par,gg,hessian=ho)
t13_ci <- transition_prob_ci_band("mild-death",tval,vt,mlo$par,gg,hessian=ho)
t03_ci <- transition_prob_ci_band("well-death",tval,vt,mlo$par,gg,hessian=ho)


par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1.2)
plot(tval,t01_ci$est,type="l",ylim=c(0,1),col="#5ab4ac",lwd=3,ylab="transition probability",xlab="time")
polygon(c(tval, rev(tval)), c(t01_ci$upper, rev(t01_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t12_ci$upper, rev(t12_ci$lower)),
        col = adjustcolor("#5ab4ac",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t23_ci$upper, rev(t23_ci$lower)),
        col = adjustcolor("#f6e8c3",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t13_ci$upper, rev(t13_ci$lower)),
        col = adjustcolor("#d8b365",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t03_ci$upper, rev(t03_ci$lower)),
        col = adjustcolor("#8c510a",alpha.f=0.2), border=NA)
lines(tval,t12_ci$est,lwd=3,col="#5ab4ac")
lines(tval,t23_ci$est,lwd=3,col="#f6e8c3")
lines(tval,t13_ci$est,lwd=3,col="#d8b365")
lines(tval,t03_ci$est,lwd=3,col="#8c510a")
legend("right",legend=c("well-mild","mild-severe","severe-death","mild-death","well-death"),
       col=c("#5ab4ac","#01665e","#f6e8c3","#d8b365","#8c510a"),lwd=2,bty="n",cex=1)

