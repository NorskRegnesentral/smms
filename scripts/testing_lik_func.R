## Testing likelihood functions

rm(list = ls())
devtools::load_all() 

## Testing finding_limits
gg1 = graph_from_literal("0"--+"1"--+"2"--+"3"--+"5","0"--+"4"--+"5")
dd1 <- data.frame(patient=c(rep(1,10)),time=c(0,0.5,1,1.5,2,2.5,3,3.5,4,5),state=c(0,0,5,5,5,5,5,5,5,5))

dd1a = arrange_data(dd1,gg1)

edge_mats=edge_matrices(gg1)
state_ord = state_ordering(gg1)
absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])

ll = finding_limits(dd1a[1,1:10],"0",edge_mats,absorbing_states)


## Testing integrand functions
gg = graph_from_literal("1"--+"2"--+"3"--+"4", "1"--+"4", "2"--+"4")
edge_mats <- edge_matrices(gg)
names_surv_dens <- names_of_survival_density(gg)

type_to_integrand("0", edge_mats,names_surv_dens)
# [1] "function(times,tt,param,x){ ss<-times[1]\n   S_01 (param,x, tt ) *S_03 (param,x, tt )  }"
type_to_integrand("01", edge_mats, names_surv_dens)
# [1] "function(times,tt,param,x){ ss<-times[1]\n f_01 (param,x, ss ) * S_03 (param,x, ss ) * S_12 (param,x, tt- ss ) *S_13 (param,x, tt- ss )  }"
type_to_integrand("012", edge_mats, names_surv_dens)
# [1] "function(times,tt,param,x){ ss<-times[1]\nuu<-times[2]\n 
# f_01 (param,x, ss ) *f_12 (param,x, uu-ss ) * S_03 (param,x, ss ) *S_13 (param,x, uu-ss ) * S_23 (param,x, tt- uu )  }"
type_to_integrand("0123", edge_mats, names_surv_dens)
# [1] "function(times,tt,param,x){ ss<-times[1]\nuu<-times[2]\n 
# f_01 (param,x, ss ) *f_12 (param,x, uu-ss ) *f_23 (param,x, tt-uu ) * S_03 (param,x, ss ) *S_13 (param,x, uu-ss )   }"
type_to_integrand("013", edge_mats, names_surv_dens)
# [1] "function(times,tt,param,x){ ss<-times[1]\n f_01 (param,x, ss ) *f_13 (param,x, tt-ss ) * S_03 (param,x, ss ) *S_12 (param,x, tt-ss )   }"
type_to_integrand("03", edge_mats, names_surv_dens)
# [1] "function(times,tt,param,x){ ss<-times[1]\n f_03 (param,x, tt ) * S_01 (param,x, tt )   }"

type_to_integrand("0123", edge_mats, names_surv_dens,abs_exact = F)
# [1] "function( ss ,uu ,tt2=NULL , tt,param,x){ f_01 (param,x, ss ) *f_12 (param,x, uu ) *
# ( S_23 (param,x, tt2- ss- uu ) - S_23 (param,x, tt- ss- uu ))* S_03 (param,x, ss ) *S_13 (param,x, uu )   }"
type_to_integrand("013", edge_mats, names_surv_dens,abs_exact = F)
# [1] "function( ss ,uu , tt,param,x){ f_01 (param,x, ss ) *f_13 (param,x, uu ) * S_03 (param,x, ss ) *S_12 (param,x, uu )   }"
type_to_integrand("03", edge_mats, names_surv_dens,abs_exact = F)
# [1] "function( ss , tt,param,x){ f_03 (param,x, ss ) * S_01 (param,x, ss )   }"

## Testing occupancy probabilities
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

#### 5 state progressive with absorbing not observed exactly
gg1 = graph_from_literal("0"--+"1"--+"2"--+"3"--+"4")
dd1 <- data.frame(patient=c(rep(1,10)),time=c(0,0.5,1,1.5,2,2.5,3,3.5,4,5),state=c(0,0,1,1,2,2,3,3,4,4))
dd2 <- data.frame(patient=rep(2,8),time=c(0,0.5,1,1.5,2,2.5,4,5),state=c(0,0,1,1,2,2,4,4))
dd3 <- data.frame(patient=rep(3,2),time=c(0,5),state=c(0,4))
dd4 <- data.frame(patient=rep(4,4),time=c(0,0.5,2,3),state=c(0,0,2,3))
dd <- rbind(dd1,dd2,dd3,dd4)

dda = arrange_data(dd,gg1)

formula_obs_types = all_types(gg1)
edge_mats <- edge_matrices(gg1)
state_ord = state_ordering(gg1)
absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
names_surv_dens = names_of_survival_density(gg1)

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

integrand
all_integral_limits

S_01 = function(param, x, t){ (1-pexp(t,exp(param[1])))}
S_12 = function(param, x, t){ (1-pexp(t,exp(param[2])))}
S_23 = function(param, x, t){ (1-pexp(t,exp(param[3])))}
S_34 = function(param, x, t){ (1-pexp(t,exp(param[4])))}

f_01 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[1]))}
f_12 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[2]))}
f_23 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[3]))}
f_34 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[4]))}

aa <- c(1.2,0.4,1.5,2,1)

mloglikelihood(aa,integrand,all_integral_limits,method1 = "hcubature",mc_cores=1)
# seems correct! (only compared for "01234")

#### 6 state progressive with absorbing not observed exactly
gg1 = graph_from_literal("0"--+"1"--+"2"--+"3"--+"4"--+"5")
dd1 <- data.frame(patient=c(rep(1,11)),time=c(0,0.5,1,1.5,2,2.5,3,3.5,4,5,5.1),state=c(0,0,1,1,2,2,3,3,4,4,5))
dd2 <- data.frame(patient=rep(2,8),time=c(0,0.5,1,1.5,2,2.5,4,5),state=c(0,0,1,1,2,2,5,5))
dd3 <- data.frame(patient=rep(3,2),time=c(0,5),state=c(0,4))
dd4 <- data.frame(patient=rep(4,4),time=c(0,0.5,2,3),state=c(0,0,2,3))
dd <- rbind(dd1,dd2,dd3,dd4)

dda = arrange_data(dd,gg1)

formula_obs_types = all_types(gg1)
edge_mats <- edge_matrices(gg1)
state_ord = state_ordering(gg1)
absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
names_surv_dens = names_of_survival_density(gg1)

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

integrand
all_integral_limits

S_01 = function(param, x, t){ (1-pexp(t,exp(param[1])))}
S_12 = function(param, x, t){ (1-pexp(t,exp(param[2])))}
S_23 = function(param, x, t){ (1-pexp(t,exp(param[3])))}
S_34 = function(param, x, t){ (1-pexp(t,exp(param[4])))}
S_45 = function(param, x, t){ (1-pexp(t,exp(param[5])))}

f_01 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[1]))}
f_12 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[2]))}
f_23 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[3]))}
f_34 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[4]))}
f_45 = function(param, x, t){as.numeric(t>=0)*dexp(t,exp(param[5]))}

aa <- c(1.2,0.4,1.5,2,1,2)

mloglikelihood(aa,integrand,all_integral_limits,method1 = "hcubature",mc_cores=1)

#### 7 state progressive with absorbing not observed exactly
gg1 = graph_from_literal("0"--+"1"--+"2"--+"3"--+"4"--+"5"--+"6")
dd1 <- data.frame(patient=c(rep(1,11)),time=c(0,0.5,1,1.5,2,2.5,3,3.5,4,5,5.1),state=c(0,0,1,1,2,2,3,3,4,4,6))
dd2 <- data.frame(patient=rep(2,8),time=c(0,0.5,1,1.5,2,2.5,4,5),state=c(0,0,1,1,2,2,6,6))
dd3 <- data.frame(patient=rep(3,2),time=c(0,5),state=c(0,4))
dd4 <- data.frame(patient=rep(4,4),time=c(0,0.5,2,3),state=c(0,0,2,3))
dd <- rbind(dd1,dd2,dd3,dd4)

dda = arrange_data(dd,gg1)

formula_obs_types = all_types(gg1)
edge_mats <- edge_matrices(gg1)
state_ord = state_ordering(gg1)
absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
names_surv_dens = names_of_survival_density(gg1)

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

integrand
all_integral_limits

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

aa <- c(1.2,0.4,1.5,2,1,2,0.3)

mloglikelihood(aa,integrand,all_integral_limits,method1 = "hcubature",mc_cores=1)

