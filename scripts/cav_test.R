#### CAV example for testing
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
dd = dd[ ,-c(2, 4, 5, 6, 7, 9, 10)]
colnames(dd)[1:2] <- c("patient","time")

# Model:
S_01 = function(param, x, t){(as.numeric(t>=0))* (1-pweibull(t,param[1],param[2]))}
S_12 = function(param, x, t){(as.numeric(t>=0))* (1-pweibull(t,param[3],param[4]))}
S_23 = function(param, x, t){(as.numeric(t>=0))* (1-pweibull(t,param[5],param[6]))}
S_03 = function(param, x, t){(as.numeric(t>=0))* (1-pweibull(t,param[7],param[8]))}
S_13 = function(param, x, t){(as.numeric(t>=0))* (1-pweibull(t,param[9],param[10]))}

f_01 = function(param, x, t){as.numeric(t>=0)*dweibull(t,param[1],param[2])}
f_12 = function(param, x, t){as.numeric(t>=0)*dweibull(t,param[3],param[4])}
f_23 = function(param, x, t){as.numeric(t>=0)*dweibull(t,param[5],param[6])}
f_03 = function(param, x, t){as.numeric(t>=0)*dweibull(t,param[7],param[8])}
f_13 = function(param, x, t){as.numeric(t>=0)*dweibull(t,param[9],param[10])}

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

params <- c(1.43,8.72,1.21,3.05,1.12,3.75,0.44,418.32,0.71,6.78)

from_time_point_to_integral(params,method1 = "hcubature", integrand = integrand,integrand2 = integrand2, 
                            all_integral_limits = all_integral_limits,mc_cores=1)

system.time({
  from_time_point_to_integral(params,method1 = "hcubature", integrand = integrand,integrand2 = integrand2,
                              all_integral_limits = all_integral_limits,mc_cores=1)
})

system.time({
  oo <- nlminb(params,from_time_point_to_integral,integrand = integrand,integrand2 = integrand2, 
              all_integral_limits = all_integral_limits,mc_cores=1,lower=rep(0.0001,10))
})

