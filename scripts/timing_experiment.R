### Timing experiment
#### up to 7 state progressive with absorbing not observed exactly
setwd("~/Multistate models/smms")
rm(list = ls())
devtools::load_all() 
library(igraph)
library(parallel)
library(cubature)
library(microbenchmark)


gg1 = graph_from_literal("1"--+"2")
gg2 = graph_from_literal("1"--+"2"--+"3")
gg3 = graph_from_literal("1"--+"2"--+"3"--+"4")
gg4 = graph_from_literal("1"--+"2"--+"3"--+"4"--+"5")
gg5 = graph_from_literal("1"--+"2"--+"3"--+"4"--+"5"--+"6")
gg6 = graph_from_literal("1"--+"2"--+"3"--+"4"--+"5"--+"6"--+"7")

ctimes <- data.frame(int_dim=c(rep(0:5,each=3)),num_par=rep(c("1","+1","2"),6),otime=rep(NA,18))
nn <- 10 #size of datasets
mm <- 50 # number of times to repeat the timing process
mc_cores <- 1 # number of cores to use - not sure what is best here?
graphs <- list(gg1,gg2,gg3,gg4,gg5,gg6)

## Expo transition times
S_01 = function(param, x, t){ (1-pexp(t,exp(param[1])))}
S_12 = function(param, x, t){ (1-pexp(t,exp(param[2])))}
S_23 = function(param, x, t){ (1-pexp(t,exp(param[3])))}
S_34 = function(param, x, t){ (1-pexp(t,exp(param[4])))}
S_45 = function(param, x, t){ (1-pexp(t,exp(param[5])))}
S_56 = function(param, x, t){ (1-pexp(t,exp(param[6])))}

f_01 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[1]))}
f_12 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[2]))}
f_23 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[3]))}
f_34 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[4]))}
f_45 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[5]))}
f_56 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[6]))}

for (k in 1:6){
  tt <- rep(NA,mm)
  
  for (m in 1:mm){
    param <- rep(1,k)
    int_dim <- k-1
    dd <- data.frame(patient=rep(1:nn,each=(k+1)),time=rep(NA,nn*(k+1)),state=rep(NA,nn*(k+1)))
    for (h in 1:nn){
      tts <- c(0,cumsum(rexp(k,exp(param[1]))))
      dd[dd$patient==h,"time"] <- c(runif(k,tts[1:k],tts[2:(k+1)]),max(tts)+runif(1,0.001,0.2))
      dd[dd$patient==h,"state"] <- 1:(k+1)
    }
    gg <- graphs[[k]]
    dda <- arrange_data(dd,gg)
    
    formula_obs_types = all_types(gg)
    edge_mats <- edge_matrices(gg)
    state_ord = state_ordering(gg)
    absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
    names_surv_dens = names_of_survival_density(gg)
    
    timepointMat <- dda[,1:(dim(dda)[2]-1)]
    
    observation_type = rep(NA, nrow(dda))
    all_integral_limits = list()
    integrand = list()
    
    for(i in 1:nrow(dda)){
      observation_type[i] = dda[i,"obs_type"]
      f_types = names(which(formula_obs_types[, observation_type[i]] == 1))
      integrand_mellomregn = list()
      integral_mellomregn= list()
      for(j in 1:length(f_types)){
        integrand_mellomregn[[j]] = eval(parse(text=type_to_integrand(f_types[j], edge_mats, names_surv_dens,abs_exact = F)))
        integral_mellomregn[[j]] = finding_limits(timepointMat[i,],f_types[j],edge_mats, absorbing_states,abs_exact = F)
      }
      all_integral_limits[[i]] = integral_mellomregn
      integrand[[i]] = integrand_mellomregn
    }
    mt <- microbenchmark(oo<-nlminb(param,mloglikelihood,integrand = integrand,limits = all_integral_limits,
                                    mc_cores=mc_cores,X=NULL),times=1)
    tt[m] <- mt$time/10^9
  }
  ctimes[which(ctimes$int_dim==(k-1) & ctimes$num_par=="1"),3] <- mean(tt)
}


## One extra param: Weibull first transition
S_01 = function(param, x, t){as.numeric(t>=0)*(1-pweibull(t,exp(param[1]),exp(param[2])))}
S_12 = function(param, x, t){as.numeric(t>=0)*(1-pexp(t,exp(param[3])))}
S_23 = function(param, x, t){as.numeric(t>=0)*(1-pexp(t,exp(param[4])))}
S_34 = function(param, x, t){as.numeric(t>=0)*(1-pexp(t,exp(param[5])))}
S_45 = function(param, x, t){as.numeric(t>=0)*(1-pexp(t,exp(param[6])))}
S_56 = function(param, x, t){as.numeric(t>=0)*(1-pexp(t,exp(param[7])))}

f_01 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[1]),exp(param[2]))}
f_12 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[3]))}
f_23 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[4]))}
f_34 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[5]))}
f_45 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[6]))}
f_56 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[7]))}


for (k in 1:6){
  tt <- rep(NA,mm)
  
  for (m in 1:mm){
    param <- c(1,0.5,rep(1,k-1))
    int_dim <- k-1
    dd <- data.frame(patient=rep(1:nn,each=(k+1)),time=rep(NA,nn*(k+1)),state=rep(NA,nn*(k+1)))
    for (h in 1:nn){
      if(k == 1){
        tts <- c(0,cumsum(c(rweibull(1,exp(param[1]),exp(param[2])))))
      }
      else if(k > 1){
      tts <- c(0,cumsum(c(rweibull(1,exp(param[1]),exp(param[2])),rexp(k-1,exp(param[3])))))
      }
      dd[dd$patient==h,"time"] <- c(runif(k,tts[1:k],tts[2:(k+1)]),max(tts)+runif(1,0.001,0.2))
      dd[dd$patient==h,"state"] <- 1:(k+1)
    }
    gg <- graphs[[k]]
    dda <- arrange_data(dd,gg)
    
    formula_obs_types = all_types(gg)
    edge_mats <- edge_matrices(gg)
    state_ord = state_ordering(gg)
    absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
    names_surv_dens = names_of_survival_density(gg)
    
    timepointMat <- dda[,1:(dim(dda)[2]-1)]
    
    observation_type = rep(NA, nrow(dda))
    all_integral_limits = list()
    integrand = list()
    
    for(i in 1:nrow(dda)){
      observation_type[i] = dda[i,"obs_type"]
      f_types = names(which(formula_obs_types[, observation_type[i]] == 1))
      integrand_mellomregn = list()
      integral_mellomregn= list()
      for(j in 1:length(f_types)){
        integrand_mellomregn[[j]] = eval(parse(text=type_to_integrand(f_types[j], edge_mats, names_surv_dens,abs_exact = F)))
        integral_mellomregn[[j]] = finding_limits(timepointMat[i,],f_types[j],edge_mats, absorbing_states,abs_exact = F)
      }
      all_integral_limits[[i]] = integral_mellomregn
      integrand[[i]] = integrand_mellomregn
    }
    mt <- microbenchmark(oo<-nlminb(param,mloglikelihood,integrand = integrand,limits = all_integral_limits,
                                    mc_cores=mc_cores,X=NULL),times=1)
    tt[m] <- mt$time/10^9
  }
  ctimes[which(ctimes$int_dim==(k-1) & ctimes$num_par=="+1"),3] <- mean(tt)
}



## Weibull transition times
S_01 = function(param, x, t){as.numeric(t>=0)*(1-pweibull(t,exp(param[1]),exp(param[2])))}
S_12 = function(param, x, t){ as.numeric(t>=0)*(1-pweibull(t,exp(param[3]),exp(param[4])))}
S_23 = function(param, x, t){as.numeric(t>=0)*(1-pweibull(t,exp(param[5]),exp(param[6])))}
S_34 = function(param, x, t){as.numeric(t>=0)*(1-pweibull(t,exp(param[7]),exp(param[8])))}
S_45 = function(param, x, t){as.numeric(t>=0)*(1-pweibull(t,exp(param[9]),exp(param[10])))}
S_56 = function(param, x, t){as.numeric(t>=0)*(1-pweibull(t,exp(param[11]),exp(param[12])))}

f_01 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[1]),exp(param[2]))}
f_12 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[3]),exp(param[4]))}
f_23 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[5]),exp(param[6]))}
f_34 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[7]),exp(param[8]))}
f_45 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[9]),exp(param[10]))}
f_56 = function(param, x, t){as.numeric(t>=0)*dweibull(t,exp(param[11]),exp(param[12]))}


for (k in 1:6){
  tt <- rep(NA,mm)
  
  for (m in 1:mm){
    param <- c(rep(c(1,0.5),k))
    int_dim <- k-1
    dd <- data.frame(patient=rep(1:nn,each=(k+1)),time=rep(NA,nn*(k+1)),state=rep(NA,nn*(k+1)))
    for (h in 1:nn){
      tts <- c(0,cumsum(rweibull(k,exp(param[1]),exp(param[2]))))
      dd[dd$patient==h,"time"] <- c(runif(k,tts[1:k],tts[2:(k+1)]),max(tts)+runif(1,0.001,0.2))
      dd[dd$patient==h,"state"] <- 1:(k+1)
    }
    gg <- graphs[[k]]
    dda <- arrange_data(dd,gg)
    
    formula_obs_types = all_types(gg)
    edge_mats <- edge_matrices(gg)
    state_ord = state_ordering(gg)
    absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
    names_surv_dens = names_of_survival_density(gg)
    
    timepointMat <- dda[,1:(dim(dda)[2]-1)]
    
    observation_type = rep(NA, nrow(dda))
    all_integral_limits = list()
    integrand = list()
    
    for(i in 1:nrow(dda)){
      observation_type[i] = dda[i,"obs_type"]
      f_types = names(which(formula_obs_types[, observation_type[i]] == 1))
      integrand_mellomregn = list()
      integral_mellomregn= list()
      for(j in 1:length(f_types)){
        integrand_mellomregn[[j]] = eval(parse(text=type_to_integrand(f_types[j], edge_mats, names_surv_dens,abs_exact = F)))
        integral_mellomregn[[j]] = finding_limits(timepointMat[i,],f_types[j],edge_mats, absorbing_states,abs_exact = F)
      }
      all_integral_limits[[i]] = integral_mellomregn
      integrand[[i]] = integrand_mellomregn
    }
    mt <- microbenchmark(oo<-nlminb(param,mloglikelihood,integrand = integrand,limits = all_integral_limits,
                                    mc_cores=mc_cores,X=NULL),times=1)
    tt[m] <- mt$time/10^9
  }
  ctimes[which(ctimes$int_dim==(k-1) & ctimes$num_par=="2"),3] <- mean(tt)
}

save(ctimes,file="~/H/Paperskriving/Masteroppgave_gamma_prosesser/Simulations/comp_times_exp1.RData")
