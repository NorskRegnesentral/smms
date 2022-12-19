#### Competing risk model with 2 causes of death: Boag (1949) dataset
### Expo model without covariates 

library(smms)
library(igraph) # For specifying the multi-state graph 

# Graph:
gg = graph_from_literal("0"--+"1", "0"--+"2")
par(mar=c(1,1,1,1))
plot(gg,layout=layout_with_sugiyama(gg,layers=c(1,2,2))$layout,vertex.size=20)

#write_loglikelihood(gg,abs_exact=T)

# Data:
dd <- read.table(file="boag1949_breastCancer.csv",header=T,sep=";",dec=",")
colnames(dd) <- c("time","state")

dd <- as.data.frame(dd)
dd$time <- dd$time/max(dd$time)

n <- dim(dd)[1]

dd$patient <- 1:n

# Specify model:
S_01 = function(param, x, t){(1-pexp(t,exp(param[1])))}
S_02 = function(param, x, t){(1-pexp(t,exp(param[2])))}

f_01 = function(param, x, t){dexp(t,exp(param[1]))}
f_02 = function(param, x, t){dexp(t,exp(param[2]))}

# Optimize:
startval <- c(0.8,-0.5)
mlo <- smms(startval,dd,gg, mc_cores = 1, hessian_matrix = T)

# Prevalence plot
tval <- seq(0.01,1,length=50) 
# a sequence of time-points over which to compute the state occupancies
p0_ci <- occupancy_prob_ci_band("0",tval,mlo$opt$par,gg,hessian=mlo$hess) 
#for computing the confidence bands, the hessian needs to be provided (but there
# exists a function computing only the occupancy probabilities too)
p1_ci <- occupancy_prob_ci_band("1",tval,mlo$opt$par,gg,hessian=mlo$hess)
p2_ci <- occupancy_prob_ci_band("2",tval,mlo$opt$par,gg,hessian=mlo$hess)

library("cmprsk")
ff <- cuminc(dd[,"time"],dd[,"state"])
prvnp1 <- rep(NA,50)
for (i in 1:50){
  prvnp1[i] <- 100*ff$`1 1`$est[max(which(ff$`1 1`$time<=tval[i]))]
}
prvnp2 <- rep(NA,50)
for (i in 1:50){
  prvnp2[i] <- 100*ff$`1 2`$est[max(which(ff$`1 2`$time<=tval[i]))]
}
prvnp0 <- rep(100,50)-prvnp1-prvnp2

par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1)
plot(tval,prvnp0,type="l",ylim=c(0,100),lwd=3,
     xlab="time",col="dark grey", ylab="prevalence (%)")
polygon(c(tval, rev(tval)), c(p0_ci$upper*100, rev(p0_ci$lower*100)),
        col = adjustcolor("#b2e2e2",alpha.f=0.2), border=NA)
lines(tval,p0_ci$est*100,col="#b2e2e2",lwd=3)

lines(tval,prvnp1,col="dark grey",lwd=3)
polygon(c(tval, rev(tval)), c(p1_ci$upper*100, rev(p1_ci$lower*100)),
        col = adjustcolor("#fdae61",alpha.f=0.2), border=NA)
lines(tval,p1_ci$est*100,col="#fdae61",lwd=3)

lines(tval,prvnp2,col="dark grey",lwd=3)
polygon(c(tval, rev(tval)), c(p2_ci$upper*100, rev(p2_ci$lower*100)),
        col = adjustcolor("#a50026",alpha.f=0.2), border=NA)
lines(tval,p2_ci$est*100,col="#a50026",lwd=3)
legend("topright",legend=c("non-parametric","alive","death by cause 1","death by cause 2"),
       col=c("dark grey","#b2e2e2","#fdae61","#a50026"),
       lwd=2,bty="n",cex=1)
# the simple expo model doesn't seem to fit so well based on these plots

# Overall survival
tval <- seq(0.01,1,length=50)
So <- overall_survival_ci_band(tval,mlo$opt$par,gg,hessian=mlo$hess)

library(survival)
km <- survfit(Surv(dd[,"time"],dd["state"]>0)~1)

par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1)
plot(km, col="dark grey",lwd=2,xlab="time",ylab="survival")
# using msm package for the nonparametric estimate (for now)
polygon(c(tval, rev(tval)), c(So$upper, rev(So$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
lines(tval,So$est,col="#0571b0",lwd=3)
# the simple expo model doesn't seem to fit very well based on these plots
