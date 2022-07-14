#### CAV dataset
### Expo model (i.e. time-homogeneous Markov) with covariates 
setwd("~/Multistate models/smms")
rm(list = ls())
devtools::load_all() 


## Preprocessing
gg = igraph::graph_from_literal("1"--+"2"--+"3"--+"4", "1"--+"4", "2"--+"4")

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
dd$dage_st = (dd$dage-mean(dd$dage))/sd(dd$dage)
dd$ihd = (dd$pdiag=="IHD")
colnames(dd)[1:2] <- c("patient","time")

X_data = aggregate(dd[,c("dage_st","ihd")],by=list(dd$patient),FUN=median) #FUN does not matter
X_data = as.matrix(X_data[,2:3])

# Model:
names_of_survival_density(gg)

S_01 = function(param, x, t){(1-pexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2])))}
S_12 = function(param, x, t){(1-pexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2])))}
S_23 = function(param, x, t){(1-pexp(t,exp(param[7]+param[8]*x[1]+param[9]*x[2])))}
S_03 = function(param, x, t){(1-pexp(t,exp(param[10]+param[11]*x[1]+param[12]*x[2])))}
S_13 = function(param, x, t){(1-pexp(t,exp(param[13]+param[14]*x[1]+param[15]*x[2])))}

f_01 = function(param, x, t){dexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2]))}
f_12 = function(param, x, t){dexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2]))}
f_23 = function(param, x, t){dexp(t,exp(param[7]+param[8]*x[1]+param[9]*x[2]))}
f_03 = function(param, x, t){dexp(t,exp(param[10]+param[11]*x[1]+param[12]*x[2]))}
f_13 = function(param, x, t){dexp(t,exp(param[13]+param[14]*x[1]+param[15]*x[2]))}


# Fitting the model (this may take some time, when using only one core)
startval <- c(-2.5,0.21,0.5,-1.1,-0.2,0.2,-1.2,-0.1,-0.5,
              -3.1,0.5,0.3,-2.8,-1.1,1.1)

mlo <- smms(startval,dd,gg,X_data, mc_cores = 1)
ho <- hessian_matrix(mlo$par,dd,gg,X_data,mc_cores=1) #calculates hessian, also a bit slow with 1 core (can be done inside smms() too)

# Compute AIC (higher values are better with this definition)
aic <- (-2*mlo$objective)-2*length(startval) #-2851.2

# Look at estimates and 95% confidence intervals
round(est_ci(mlo$par,ho),2) #parameters living on the real line

round(exp(est_ci(mlo$par,ho)),2) #parameters living on the positive halv-line (more interpretable for some parameters)

#### Occupancy probabilities (with confidence bands), for some specific covariate values
tval <- seq(0.01,30,length=50)
p0_y0_ci <- occupancy_prob_ci_band("0",tval,mlo$par,gg,xval=c(-1,0),ho)
p0_y1_ci <- occupancy_prob_ci_band("0",tval,mlo$par,gg,xval=c(-1,1),ho)
p0_o0_ci <- occupancy_prob_ci_band("0",tval,mlo$par,gg,xval=c(1,0),ho)
p0_o1_ci <- occupancy_prob_ci_band("0",tval,mlo$par,gg,xval=c(1,1),ho)
p1_y0_ci <- occupancy_prob_ci_band("1",tval,mlo$par,gg,xval=c(-1,0),ho)
p1_y1_ci <- occupancy_prob_ci_band("1",tval,mlo$par,gg,xval=c(-1,1),ho)
p1_o0_ci <- occupancy_prob_ci_band("1",tval,mlo$par,gg,xval=c(1,0),ho)
p1_o1_ci <- occupancy_prob_ci_band("1",tval,mlo$par,gg,xval=c(1,1),ho)
p2_y0_ci <- occupancy_prob_ci_band("2",tval,mlo$par,gg,xval=c(-1,0),ho)
p2_y1_ci <- occupancy_prob_ci_band("2",tval,mlo$par,gg,xval=c(-1,1),ho)
p2_o0_ci <- occupancy_prob_ci_band("2",tval,mlo$par,gg,xval=c(1,0),ho)
p2_o1_ci <- occupancy_prob_ci_band("2",tval,mlo$par,gg,xval=c(1,1),ho)
p3_y0_ci <- occupancy_prob_ci_band("3",tval,mlo$par,gg,xval=c(-1,0),ho)
p3_y1_ci <- occupancy_prob_ci_band("3",tval,mlo$par,gg,xval=c(-1,1),ho)
p3_o0_ci <- occupancy_prob_ci_band("3",tval,mlo$par,gg,xval=c(1,0),ho)
p3_o1_ci <- occupancy_prob_ci_band("3",tval,mlo$par,gg,xval=c(1,1),ho)


# msm package with covariates (to get non-parameteric prevalence curves)
oneway4.q <- rbind(c(0, 0.25, 0, 0.25), c(0, 0, 0.25, 0.25),c(0, 0, 0, 0.5), c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well", "Mild","Severe", "Death")
cav.msm <- msm(state ~ time, subject = patient, data = dd,covariates=~dage_st+ihd,
               qmatrix = oneway4.q, death = 4,method = "BFGS", control = list(fnscale = 4000, maxit = 10000))

prev_y0 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=-1,ihdTRUE=0))
prev_y1 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=-1,ihdTRUE=1))
prev_o0 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=1,ihdTRUE=0))
prev_o1 <- prevalence.msm(cav.msm,times=tval,covariates=list(dage_st=1,ihdTRUE=1))

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
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=1)

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
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=1)

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
legend("bottomright",legend=c("youger donor, no IHD","youger donor, IHD","older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",lty=c(1,1,2,2),cex=1)



#### Overall survival curves (with confidence bands), for some specific covariate values
tval <- seq(0.01,30,length=50)
Sy0 <- overall_survival_ci_band(tval,mlo$par,gg,c(-1,0),ho)
Sy1 <- overall_survival_ci_band(tval,mlo$par,gg,c(-1,1),ho)
So0 <- overall_survival_ci_band(tval,mlo$par,gg,c(1,0),ho)
So1 <- overall_survival_ci_band(tval,mlo$par,gg,c(1,1),ho)


par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1.2)
plot.survfit.msm(cav.msm, col.surv="black",lwd.surv=2,xlab="years after transplantation",
                 ylab="survival",main=" ",legend.pos=None,col=NULL,lwd=3,covariates=list(dage_st=-1,ihdTRUE=0))
# using msm package for the nonparametric estimate (for now)
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


#### Transition probabilities (with confidence bands), for some specific covariate values
tval <- seq(2,30,length=50)
vt <- 1
t_y0_ci <- transition_prob("1","0",tval,vt,mlo$par,gg,xval=c(-1,0))
p0_y1_ci <- occupancy_prob_ci_band("0",tval,mlo$par,gg,xval=c(-1,1),ho)
p0_o0_ci <- occupancy_prob_ci_band("0",tval,mlo$par,gg,xval=c(1,0),ho)
p0_o1_ci <- occupancy_prob_ci_band("0",tval,mlo$par,gg,xval=c(1,1),ho)
p1_y0_ci <- occupancy_prob_ci_band("1",tval,mlo$par,gg,xval=c(-1,0),ho)
p1_y1_ci <- occupancy_prob_ci_band("1",tval,mlo$par,gg,xval=c(-1,1),ho)
p1_o0_ci <- occupancy_prob_ci_band("1",tval,mlo$par,gg,xval=c(1,0),ho)
p1_o1_ci <- occupancy_prob_ci_band("1",tval,mlo$par,gg,xval=c(1,1),ho)
p2_y0_ci <- occupancy_prob_ci_band("2",tval,mlo$par,gg,xval=c(-1,0),ho)
p2_y1_ci <- occupancy_prob_ci_band("2",tval,mlo$par,gg,xval=c(-1,1),ho)
p2_o0_ci <- occupancy_prob_ci_band("2",tval,mlo$par,gg,xval=c(1,0),ho)
p2_o1_ci <- occupancy_prob_ci_band("2",tval,mlo$par,gg,xval=c(1,1),ho)
p3_y0_ci <- occupancy_prob_ci_band("3",tval,mlo$par,gg,xval=c(-1,0),ho)
p3_y1_ci <- occupancy_prob_ci_band("3",tval,mlo$par,gg,xval=c(-1,1),ho)
p3_o0_ci <- occupancy_prob_ci_band("3",tval,mlo$par,gg,xval=c(1,0),ho)
p3_o1_ci <- occupancy_prob_ci_band("3",tval,mlo$par,gg,xval=c(1,1),ho)


