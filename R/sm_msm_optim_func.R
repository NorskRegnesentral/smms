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

#' Fit a semi-Markovian multistate model
#'
#' The function the user should interact with. Assumes that the user has provided densities and 
#' survival functions of the right format (ideally with all parameters living on -Inf to Inf - 
#' to avoid having to set lower limits for some parameters).
#'
#' @param startval The parameter values where the optimisation routine will start. The dimension
#' and ordering is given by the user-specified densities.
#' @param data A data frame where the number of rows correspond to the total number of observations for all patients,
#' and  with 3 columns: patient (numbering or names for each patient), time (the time when a patient was observed),
#' state (the state which the patient occupies at the observation time). FOR NOW!!
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param X A matrix with the covariates. Number of rows equal to the number of patients and one column
#' for each covariate. The covariate specification is given by the user-specified densities.
#' @param mc_cores The number of cores to use (for parallelisation). The function uses the mclapply()
#' function from the parallel package.
#' @return The result from optimising the log-likelihood: parameter estimates with corresponding variance-covariance
#' matrix, and the maximum log-likelihood value.
#' 

# This function needs to be updated so that it takes column names as input: to know which column corresponds
# to "patient" (numbering or names for each patient), "time" (the time when a patient was observed),
# "state" (the state which the patient occupies at the observation time).
smms = function(startval, data, graph, X = NULL, mc_cores = 3){
  formula_obs_types = all_types(graph)
  edge_mats = edge_matrices(graph)
  state_ord = state_ordering(graph)
  names_surv_dens = names_of_survival_density(graph)
  absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
  
  all_data_set = arrange_data(data, graph)
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
      integral_mellomregn[[j]] = finding_limits(timepointMat[i,],f_types[j],edge_mats,absorbing_states)
    }
    all_integral_limits[[i]] = integral_mellomregn
    integrand[[i]] = integrand_mellomregn
  }

    
  optimizer <- optim(startval,mloglikelihood,integrand = integrand,limits = all_integral_limits,X=X, method = "L-BFGS",
               mc_cores=mc_cores,hessian = FALSE)
  hessian_optimizer = numDeriv::hessian(mloglikelihood, optimizer$par, integrand = integrand,limits = all_integral_limits,
                              mc_cores=mc_cores,X=X)
  return(list(optimizer, hessian_optimizer))
}

#optim_func(params = params, dd, gg)



#' Estimates and confidence intervals
#'
#' Provide estimates and approximate confidence intervals for all parameters 
#' (on the "original" scale - meaning that all parameters live on -Inf to +Inf).
#'
#' @param param  The parameter estimates.
#' @param hessian The hessian matrix.
#' @param level The level of confidence. Default value is 0.95.
#' @param log Boolean denoting whether the intervals should be on the log scale (default), or on the exp/hazard scale.
#' @return A vector of the same length as time.
#' 
est_ci = function(param, hessian,level=0.95,log=TRUE){
  zz <- qnorm(0.5*(1-level),lower.tail = F)
  varCov <- solve(hessian)
  if (log){
    tt <- data.frame(estimate=param,lower.ci=param-zz*sqrt(diag(varCov)),upper.ci=param+zz*sqrt(diag(varCov)))
  }else{
    ses <- sqrt(diag(varCov))*exp(param)
    tt <- data.frame(estimate=exp(param),lower.ci=exp(param)-zz*ses,upper.ci=exp(param)+zz*ses)
  }
  
  return(tt)
}

###################### Functions for diagnostic plots ##############


#' Compute the occupancy probability over time
#'
#' The last state is often heavy to compute, but recall that is has to be equal to 1 minus the others
#'
#' @param state A string with a number indicating the state for which to compute the probability (in the ordered naming system). 
#' @param time A vector of timepoints for which to compute the occupancy probability.
#' @param param  The parameter values in which the probabilities should be evaluated. The dimension
#' and ordering is given by the user-specified densities.
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param xval A vector of covariate values.
#' @return A vector of the same length as time.
#' 
occupancy_prob = function(state, time, param, graph, xval = NULL){
  edge_mats = edge_matrices(graph)
  names_surv_dens = names_of_survival_density(graph)
  
  state_ord = state_ordering(graph)
  absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
  
  f_types = construct_formula_types(graph)
  f_types_match = f_types[which(substr(f_types,nchar(f_types),nchar(f_types)) == state)]
  mm <- length(f_types_match)
  
  kk <- length(time)
  op <- rep(0,kk)
  for (i in 1:mm){
    f_type = f_types_match[i]
    integrand = eval(parse(text=type_to_integrand(f_type,edge_mats, names_surv_dens,abs_exact=FALSE))) # always with abs_exact=FALSE
    
    choice <- (length(which(edge_mats$passedBy[f_type,]==max(edge_mats$traveled[f_type,])))>0)  #TRUE: if absorbing state is reached from a state where there was an option to go a different way
    
    if (substr(f_type,nchar(f_type),nchar(f_type)) %in% absorbing_states & choice==FALSE){ #if patient is observed in absorbing state
      dim_int <- nchar(f_type)-2  
    }else{ #if patient has not reached absorbing state, or has reached it from a node with a choice
      dim_int <- nchar(f_type)-1 
    }
    
    opi <- rep(NA,kk)
    for (j in 1:kk){
      lower <- rep(0,dim_int)
      upper <- rep(time[j],dim_int)
      tmax <- c(time[j],-1)
      
      if(length(lower) == 0){
        opi[j] = integrand(times=1,tt = tmax[1], tt2=tmax[2],param = param, x = xval) # the value of times does not matter
      }else if(length(lower)>0){
        if(length(lower)<=2 ){
          #opi[j] = repintegrate(integrand,tt = tmax[1], tt2=tmax[2],lower=lower,upper = upper, param = param, x = xval)
          
          opi[j] = tryCatch({
            repintegrate(integrand,tt=tmax[1],tt2=tmax[2],lower=lower,upper = upper, param = param,
                         x = xval)
          },error=function(cond){
            integrand2 <- change_integrand(integrand)
            llij = cubintegrate(integrand2, lower = lower,upper = upper, method = "divonne", maxEval = 500,
                                  tt = tmax[1], tt2=tmax[2],param = param, x = xval)$integral
            return(llij)
          })
          
        }else if (length(lower)>2){
          opi[j] = cubintegrate(integrand, lower = lower,upper = upper, method = "divonne", maxEval = 500,
                                tt = tmax[1], tt2=tmax[2],param = param, x = xval)$integral
        }
      }
    }
    op = op + opi
  }
  return(op)
}

#' Compute uncertainty bands
#' 
occupancy_prob_ci_band <- function(state,time,param,graph,xval,hessian,level=0.95){
  est <- occupancy_prob(state=state,time=time,param=param,graph=graph,xval=xval)
  delta <- occupancy_prob_delta(state=state,time=time,param=param,graph=graph,xval=xval)
  varCov <- solve(hessian)
  kk <- length(time)
  lci <- rep(NA,kk)
  uci <- rep(NA,kk)
  zz <- qnorm(0.5*(1-level),lower.tail = F)
  for (i in 1:kk){
    sei <- sqrt(delta[i,,drop=F]%*%varCov%*%t(delta[i,,drop=F]))
    lci[i] <- est[i]-zz*sei
    uci[i] <- est[i]+zz*sei
  }
  return(list(est=est,lower=lci,upper=uci))
}

occupancy_prob_delta = function(state, time, param, graph, xval = NULL){
  edge_mats = edge_matrices(graph)
  names_surv_dens = names_of_survival_density(graph)
  
  state_ord = state_ordering(graph)
  absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
  
  f_types = construct_formula_types(graph)
  f_types_match = f_types[which(substr(f_types,nchar(f_types),nchar(f_types)) == state)]
  mm <- length(f_types_match)
  
  kk <- length(time)
  op <-  matrix(0,kk,length(param))
  for (i in 1:mm){
    f_type = f_types_match[i]
    integrand = eval(parse(text=type_to_integrand(f_type,edge_mats, names_surv_dens,abs_exact=FALSE))) # always with abs_exact=FALSE
    
    choice <- (length(which(edge_mats$passedBy[f_type,]==max(edge_mats$traveled[f_type,])))>0)  #TRUE: if absorbing state is reached from a state where there was an option to go a different way
    
    if (substr(f_type,nchar(f_type),nchar(f_type)) %in% absorbing_states & choice==FALSE){ #if patient is observed in absorbing state
      dim_int <- nchar(f_type)-2  
    }else{ #if patient has not reached absorbing state, or has reached it from a node with a choice
      dim_int <- nchar(f_type)-1 
    }
    
    opi <-  matrix(0,kk,length(param))
    for (j in 1:kk){
      lower <- rep(0,dim_int)
      upper <- rep(time[j],dim_int)
      tmax <- c(time[j],-1)
      
      if(length(lower) == 0){
        opi[j,] = pracma::grad(integrand,x0=param,times=1,tt = tmax[1], tt2=tmax[2],x = xval) # the value of times does not matter
      }else if(length(lower)>0){
        if(length(lower)<=2 ){
          #opi[j,] = pracma::grad(repintegrate,x0=param,innerfunc=integrand,tt = tmax[1], tt2=tmax[2],lower=lower,upper = upper,x = xval)
          
          opi[j,] = tryCatch({
            pracma::grad(repintegrate,x0=param,innerfunc=integrand,tt=tmax[1],tt2=tmax[2],lower=lower,upper = upper,x = xval)
          },error=function(cond){
            integrand2 <- change_integrand(integrand)
            llij = pracma::grad(cubint,x0=param,integrand=integrand2,lower = lower,upper = upper, tmax=tmax,xval=xval)
            return(llij)
          })
          
        }else if (length(lower)>2){
          opi[j,] = pracma::grad(cubint,x0=param,integrand=integrand,lower = lower,upper = upper, tmax=tmax,xval=xval)
          
        }
      }
    }
    op = op + opi
  }
  return(op)
}


#' Compute the overall survival
#'
#' ...
#'
#' @param time A vector of timepoints for which to compute the occupancy probability.
#' @param param  The parameter values in which the probabilities should be evaluated. The dimension
#' and ordering is given by the user-specified densities.
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param xval A vector of covariate values.
#' @return A vector of the same length as time.
#' 
overall_survival = function(time, param, graph, xval = NULL){
  edge_mats = edge_matrices(graph)
  names_surv_dens = names_of_survival_density(graph)
  
  state_ord = state_ordering(graph)
  absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
  
  f_types = construct_formula_types(graph)
  f_types_match = f_types[which(!(substr(f_types,nchar(f_types),nchar(f_types)) %in% absorbing_states))]
  mm <- length(f_types_match)
  
  kk <- length(time)
  op <- rep(0,kk)
  for (i in 1:mm){
    f_type = f_types_match[i]
    integrand = eval(parse(text=type_to_integrand(f_type,edge_mats, names_surv_dens,abs_exact=FALSE))) # always with abs_exact=FALSE

    dim_int <- nchar(f_type)-1 
    
    opi <- rep(NA,kk)
    for (j in 1:kk){
      lower <- rep(0,dim_int)
      upper <- rep(time[j],dim_int)
      tmax <- c(time[j],-1)
      
      if(length(lower) == 0){
        opi[j] = integrand(times=1,tt = tmax[1], tt2=tmax[2],param = param, x = xval) # the value of times does not matter
      }else if(length(lower)>0){
        if(length(lower)<=2 ){
          opi[j] = repintegrate(integrand,tt = tmax[1], tt2=tmax[2],lower=lower,upper = upper, param = param, x = xval)
          
          # opi[j] = tryCatch({
          #   repintegrate(integrand,tt=tmax[1],tt2=tmax[2],lower=lower,upper = upper, param = param, 
          #                x = xval)
          # },error=function(cond){
          #   integrand2 <- change_integrand(integrand)
          #   llij = cubintegrate(integrand2, lower = lower,upper = upper, method = "divonne", maxEval = 500,
          #                         tt = tmax[1], tt2=tmax[2],param = param, x = xval)$integral
          #   return(llij)
          # })
          
        }else if (length(lower)>2){
          opi[j] = cubintegrate(integrand, lower = lower,upper = upper, method = "divonne", maxEval = 500,
                                tt = tmax[1], tt2=tmax[2],param = param, x = xval)$integral
        }
      }
    }
    op = op + opi
  }
  return(op)
}

#' Compute uncertainty bands
#' 
overall_survival_ci_band <- function(time,param,graph,xval,hessian,level=0.95){
  est <- overall_survival(time=time,param=param,graph=graph,xval=xval)
  delta <- overall_survival_delta(time=time,param=param,graph=graph,xval=xval)
  varCov <- solve(hessian)
  kk <- length(time)
  lci <- rep(NA,kk)
  uci <- rep(NA,kk)
  zz <- qnorm(0.5*(1-level),lower.tail = F)
  for (i in 1:kk){
    sei <- sqrt(delta[i,,drop=F]%*%varCov%*%t(delta[i,,drop=F]))
    lci[i] <- est[i]-zz*sei
    uci[i] <- est[i]+zz*sei
  }
  return(list(est=est,lower=lci,upper=uci))
}

overall_survival_delta = function(time, param, graph, xval = NULL){
  edge_mats = edge_matrices(graph)
  names_surv_dens = names_of_survival_density(graph)
  
  state_ord = state_ordering(graph)
  absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
  
  f_types = construct_formula_types(graph)
  f_types_match = f_types[which(!(substr(f_types,nchar(f_types),nchar(f_types)) %in% absorbing_states))]
  mm <- length(f_types_match)
  
  kk <- length(time)
  op <- matrix(0,kk,length(param))
  for (i in 1:mm){
    f_type = f_types_match[i]
    integrand = eval(parse(text=type_to_integrand(f_type,edge_mats, names_surv_dens,abs_exact=FALSE))) # always with abs_exact=FALSE
    
    dim_int <- nchar(f_type)-1 
    
    opi <- matrix(0,kk,length(param))
    for (j in 1:kk){
      lower <- rep(0,dim_int)
      upper <- rep(time[j],dim_int)
      tmax <- c(time[j],-1)
      
      if(length(lower) == 0){
        opi[j,] = pracma::grad(integrand,x0=param,times=1,tt = tmax[1], tt2=tmax[2], x=xval)
      }else if(length(lower)>0){
        if(length(lower)<=2 ){
          opi[j,] = pracma::grad(repintegrate,x0=param,innerfunc=integrand,tt = tmax[1], tt2=tmax[2],lower=lower,upper = upper,x = xval)
          
        }else if (length(lower)>2){
          #opi[j,] = pracma::grad(cubintegrate,x0=param,f=integrand, lower = lower,upper = upper, method = "divonne", maxEval = 500,
          #                      tt = tmax[1], tt2=tmax[2], x = xval)
          #needs to put cubintegrate in a separate function
          print("not ready yet")
        }
      }
    }
    op = op + opi
  }
  return(op)
}


