library(numDeriv)
## NOT DONE
# start_params = function(data_set, graph){
#   names_surv_dens = names_of_survial_density(graph)
#   which_abs = which(names_surv_dens[,"type"] == "abs")
#   edge_abs = names_surv_dens[which_abs, "all_transitions"]
#   all_data_set = arrange_data(data_set, graph)
#   midpoint = matrix(NA, nrow = nrow(all_data_set), ncol = 6)
#   for(i in 1:nrow(all_data_set)){
#     observation_type[i] = all_data_set[i,"obs_type"]
#     type = names(which(formula_obs_types[, observation_type[i]] == 1))
#     time_points = all_data_set[i,1:(ncol(all_data_set)-1)]
#     if( observational_type[i] == type){
#       integral_limits = finding_limits_integral(i, type, graph, all_edges, absorbing_state_new, all_data_set)
#       na_omit_intergal_limits = integral_limits[!is.na(integral_limits)]
#       if("tmin_int1" %in% names(na_omit_intergal_limits)){
#         midpoint[i, 1] = (na_omit_intergal_limits["tmax_int1"]-na_omit_intergal_limits["tmin_int1"])/2
#       }
#       if("tmin_int2" %in% names(na_omit_intergal_limits)){
#         midpoint[i, 2] = midpoint[i, 1] + (na_omit_intergal_limits["tmax_int2"]-na_omit_intergal_limits["tmin_int2"])/2
#       }
#       if("tmin_int3" %in% names(na_omit_intergal_limits)){
#         midpoint[i, 3] = midpoint[i, 2] + (na_omit_intergal_limits["tmax_int3"]-na_omit_intergal_limits["tmin_int3"])/2
#       }
#       if("tmin_int4" %in% names(na_omit_intergal_limits)){
#         midpoint[i, 4] = midpoint[i, 3] + (na_omit_intergal_limits["tmax_int4"]-na_omit_intergal_limits["tmin_int4"])/2
#       }
#       midpoint[i, 5] = na_omit_intergal_limits["tmax"]
#     } else{
#       time_points = all_data_set[i,1:(ncol(all_data_set)-1)]
#     }
#   }
# }



optim_func = function(params, data_set, graph, X = NULL, mc_cores = 3){
  formula_obs_types = all_types(graph)
  edge_mats <- edge_matrices(graph)
  all_transitions = state_ordering(gg)
  names_surv_dens = names_of_survial_density(graph)
  which_abs = which(names_surv_dens[,"type"] == "abs")
  edge_abs = names_surv_dens[which_abs, "all_transitions"]
  absorbing_state_new = all_transitions[which(all_transitions[,"type"] == "abs"), "order"]
  all_data_set = arrange_data(data_set, graph)
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
    for(j in 1:length(type)){
      integrand_mellomregn[[j]] = eval(parse(text=type_to_integrand_absExact(type[j], edge_mats, edge_abs)))
      integrand_mellomregn2[[j]] = eval(parse(text=type_to_integrand_absExact_v2(type[j], edge_mats, edge_abs)))
      integral_mellomregn[[j]] = finding_limits_integral(i, type[j], graph, all_edges, absorbing_state_new, all_data_set)
    }
    all_integral_limits[[i]] = integral_mellomregn
    integrand[[i]] = integrand_mellomregn
    integrand2[[i]] = integrand_mellomregn2
  }

    
  optimizer <- optim(params,from_time_point_to_integral, method = "L-BFGS", integrand = integrand,integrand2 = integrand2, 
               all_integral_limits = all_integral_limits,mc_cores=mc_cores,X=X,lower=rep(0.0001,5), hessian = FALSE)
  hessian_optimizer = hessian(from_time_point_to_integral, optimizer$par, integrand = integrand,integrand2 = integrand2, 
                              all_integral_limits = all_integral_limits,mc_cores=mc_cores,X=X)
  return(list(optimizer, hessian_optimizer))
}