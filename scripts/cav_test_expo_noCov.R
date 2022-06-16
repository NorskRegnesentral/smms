#### CAV example for testing
### Expo model without covariates (to see that we get the same as Jackson)
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
dd = dd[ ,-c(2, 4, 5, 6, 7, 9, 10)]
colnames(dd)[1:2] <- c("patient","time")

# Model:
S_01 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,exp(param[1])))}
S_12 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,exp(param[2])))}
S_23 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,exp(param[3])))}
S_03 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,exp(param[4])))}
S_13 = function(param, x, t){(as.numeric(t>=0))* (1-pexp(t,exp(param[5])))}

f_01 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[1]))}
f_12 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[2]))}
f_23 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[3]))}
f_03 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[4]))}
f_13 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[5]))}


all_data_set = arrange_data(dd, gg)

## Part 3: From the time points of a given patient to an integral
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

params <- c(0.1,0.3,0.3,0.1,0.1)

mloglikelihood(log(params),integrand,all_integral_limits,method1 = "hcubature",X=NULL,mc_cores=1)

system.time({
  oo <- nlminb(log(params),mloglikelihood,integrand = integrand,limits = all_integral_limits,mc_cores=1,X=NULL)
})
2*oo$objective
# 2877.069. we find exactly the same result as in Marthes thesis (with the MSM package)

