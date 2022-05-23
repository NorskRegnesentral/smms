#### CAV example for testing
### Expo model with covariates (to see that we get the same as Jackson)
setwd("H:/Multistate models/SemiMarkovMultistate")
rm(list = ls())
library(igraph)
library(parallel)
library(cubature)
## For test-data
library(msm)

source("sm_msm_preprocessing_func.R")
source("sm_msm_likelihood_func.R")

## Preprocessing
gg = graph_from_literal("1"--+"2"--+"3", "1"--+"2", "1"--+"3")
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
S_01 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,param[1]*exp(param[2]*x[1]+param[3]*x[2])))}
S_12 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,param[4]*exp(param[5]*x[1]+param[6]*x[2])))}
S_23 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,param[7]*exp(param[8]*x[1]+param[9]*x[2])))}
S_03 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,param[10]*exp(param[11]*x[1]+param[12]*x[2])))}
S_13 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,param[13]*exp(param[14]*x[1]+param[15]*x[2])))}

f_01 = function(param, x, t){as.numeric(t>=0)*dexp(t,param[1]*exp(param[2]*x[1]+param[3]*x[2]))}
f_12 = function(param, x, t){as.numeric(t>=0)*dexp(t,param[4]*exp(param[5]*x[1]+param[6]*x[2]))}
f_23 = function(param, x, t){as.numeric(t>=0)*dexp(t,param[7]*exp(param[8]*x[1]+param[9]*x[2]))}
f_03 = function(param, x, t){as.numeric(t>=0)*dexp(t,param[10]*exp(param[11]*x[1]+param[12]*x[2]))}
f_13 = function(param, x, t){as.numeric(t>=0)*dexp(t,param[13]*exp(param[14]*x[1]+param[15]*x[2]))}

densities <- c(f_01,f_03,f_12,f_13,f_23)
names(densities) <- c("01","03","12","13","23")

survival_functions <- c(S_01,S_03,S_12,S_13,S_23)
names(survival_functions) <- c("01","03","12","13","23") 

## The edges for the transition to the absorbing state written in the "new" way
edge_abs <- c("03","13","23")
edge_mats <- edge_matrices(gg)


## Part 3: From the time points of a given patient to an integral

all_edges = get.edgelist(gg)
all_absorbing = which(!(all_edges[,2] %in% all_edges[,1]))
absorbing_state = unique(sapply(all_absorbing, function(p) all_edges[p,2]))
all_states_old = state_ordering(gg)[,1]
row_absorbing = sapply(1:length(absorbing_state), function(p) which(all_states_old == absorbing_state[p], arr.ind = TRUE))
absorbing_state_new = sapply(row_absorbing, function(p) state_ordering(gg)[p,2])
all_data_set = arrange_data(data_set = dd, gg)
X_data_set = aggregate(dd[,c("dage_st","ihd")],by=list(dd$patient),FUN=median)
X_data_set = as.matrix(X_data_set[,2:3])

## Finding all types, integrand, time points ++
formula_obs_types = all_types(gg)
all_data_set = arrange_data(data_set = dd, gg)
observation_type = rep(NA, nrow(all_data_set))
all_integral_limits = list()
integrand = list()
integrand2 = list()
all_types = list()
for(i in 1:nrow(all_data_set)){
  observation_type[i] = all_data_set[i,"obs_type"]
  type = names(which(formula_obs_types[, observation_type[i]] == 1))
  integrand_mellomregn = list()
  integrand_mellomregn2 = list()
  integral_mellomregn= list()
  type_1 = list()
  for(j in 1:length(type)){
    integrand_mellomregn[[j]] = eval(parse(text=type_to_integrand_absExact(type[j], edge_mats, edge_abs)))
    integrand_mellomregn2[[j]] = eval(parse(text=type_to_integrand_absExact_v2(type[j], edge_mats, edge_abs)))
    integral_mellomregn[[j]] = finding_limits_integral(i, type[j], gg, all_edges, absorbing_state_new, all_data_set)
    type_1[[j]] = type[j]
  }
  all_integral_limits[[i]] = integral_mellomregn
  integrand[[i]] = integrand_mellomregn
  integrand2[[i]] = integrand_mellomregn2
  all_types[[i]] = type_1
}

params <- c(0.08,0.21,log(1.65),0.33,log(0.8),log(1.2),0.29,log(0.9),log(0.6),0.044,log(1.6),log(1.3),0.06,log(0.34),log(3))

from_time_point_to_integral(params,method1 = "hcubature", integrand = integrand,integrand2 = integrand2, 
                            all_integral_limits = all_integral_limits,X=X_data_set,mc_cores=1)

system.time({
  from_time_point_to_integral(params,method1 = "hcubature", integrand = integrand,integrand2 = integrand2,
                              all_integral_limits = all_integral_limits,X=X_data_set,mc_cores=1)
})

system.time({
  oo <- nlminb(params,from_time_point_to_integral,integrand = integrand,integrand2 = integrand2, 
               all_integral_limits = all_integral_limits,X=X_data_set,mc_cores=1,lower=rep(c(0.0001,-Inf,-Inf),5))
})
2*oo$objective
# we find exactly the same result as in Marthes thesis (with the MSM package)
