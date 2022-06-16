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
S_01 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2])))}
S_12 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2])))}
S_23 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,exp(param[7]+param[8]*x[1]+param[9]*x[2])))}
S_03 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,exp(param[10]+param[11]*x[1]+param[12]*x[2])))}
S_13 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,exp(param[13]+param[14]*x[1]+param[15]*x[2])))}

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
    integral_mellomregn[[j]] = finding_limits(timepointMat[i,],f_types[j],absorbing_states)
  }
  all_integral_limits[[i]] = integral_mellomregn
  integrand[[i]] = integrand_mellomregn
}

params <- c(log(0.08),0.21,log(1.65),log(0.33),log(0.8),log(1.2),log(0.29),log(0.9),log(0.6),
            log(0.044),log(1.6),log(1.3),log(0.06),log(0.34),log(3))

mloglikelihood(params,integrand,all_integral_limits,method1 = "hcubature",X=X_data_set,mc_cores=1)

system.time({
  oo <- nlminb(params,mloglikelihood,integrand = integrand,limits = all_integral_limits,
               mc_cores=1,X=X_data_set)
})
2*oo$objective
# 2821.21 we find exactly the same result as in Marthes thesis (with the MSM package)
