### Multi-state functions part 3 - general smms function, and functions for quantities of interest


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
#' @param abs_exact A boolean indicating whether the time of entrance into absorbing states is observed
#' exactly (TRUE) or not (FALSE). Default value is TRUE.
#' @param mc_cores The number of cores to use (for parallelisation). The function uses the mclapply()
#' function from the parallel package. Defaults to 1.
#' @param hessian_matrix Whether the hessian matrix (observed Fisher information matrix) for the parameter estimates 
#' should be calculated or not. 
#' @param cmethod The integration method of choice for the cubintegrate() function. Only for integrals 
#' of higher dimension than 2. Defaults to "hcubature".
#' @return The result from optimising the log-likelihood: parameter estimates with corresponding variance-covariance
#' matrix if variance_matrix = TRUE, and the maximum log-likelihood value.
#' 
smms = function(startval, data, graph, X = NULL, abs_exact=TRUE, mc_cores = 1, hessian_matrix = FALSE,cmethod = "hcubature"){
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
      integrand_mellomregn[[j]] = eval(parse(text=type_to_integrand(f_types[j], edge_mats, names_surv_dens,abs_exact=abs_exact)))
      integral_mellomregn[[j]] = finding_limits(timepointMat[i,],f_types[j],edge_mats,absorbing_states,abs_exact=abs_exact)
    }
    all_integral_limits[[i]] = integral_mellomregn
    integrand[[i]] = integrand_mellomregn
  }
  
  
  # optimizer <- stats::optim(startval,mloglikelihood,integrand = integrand,limits = all_integral_limits,X=X, method = "L-BFGS",
  #                    mc_cores=mc_cores,hessian = FALSE)
  
  optimizer <- stats::nlminb(startval,mloglikelihood,integrand = integrand,limits = all_integral_limits,X=X, cmethod=cmethod,
                            mc_cores=mc_cores,hessian = FALSE, lower=rep(-50,length(startval)),upper=rep(50,length(startval)))
  
  if(hessian_matrix == TRUE){
    hessian_optimizer = numDeriv::hessian(mloglikelihood, optimizer$par, integrand = integrand,limits = all_integral_limits, cmethod=cmethod,
                                          mc_cores=mc_cores,X=X)
    return(list(optimizer, hessian_optimizer))
  } else{
    return(optimizer)
  }
}

#' Calculates the hessian matrix (observed Fisher information matrix).
#'
#'
#' @param param  The parameter estimates from optimizer$par in smms.
#' @param data A data frame where the number of rows correspond to the total number of observations for all patients,
#' and  with 3 columns: patient (numbering or names for each patient), time (the time when a patient was observed),
#' state (the state which the patient occupies at the observation time). FOR NOW!!
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param X A matrix with the covariates. Number of rows equal to the number of patients and one column
#' for each covariate. The covariate specification is given by the user-specified densities.
#' @param mc_cores The number of cores to use (for parallelisation). The function uses the mclapply()
#' function from the parallel package. Defaults to 1.
#' @param cmethod The integration method of choice for the cubintegrate() function. Only for integrals 
#' of higher dimension than 2. Defaults to "hcubature".
#' @return A matrix with the hessian matrix.
#' 
hessian_matrix = function(param, data, graph, X = NULL, mc_cores = 1,cmethod = "hcubature"){
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
  hessian_optimizer = numDeriv::hessian(mloglikelihood, param, integrand = integrand,limits = all_integral_limits,
                                        mc_cores=mc_cores,X=X,cmethod=cmethod)
  return(hessian_optimizer)
}

#' Estimates and confidence intervals
#'
#' Provide estimates and approximate confidence intervals for all parameters 
#' (on the "original" scale - meaning that all parameters live on -Inf to +Inf).
#'
#' @param param  The parameter estimates.
#' @param hessian The hessian matrix.
#' @param level The level of confidence. Default value is 0.95.
#' @param pos Boolean denoting whether the estimates and intervals should live on the real line (default), or on the postive half-line (pos=TRUE).
#' @return A table with one row per estimate, and 3 columns (estimate, lower.ci, upper.ci).
#' 
est_ci = function(param, hessian,level=0.95,pos=FALSE){
  zz <- stats::qnorm(0.5*(1-level),lower.tail = F)
  varCov <- solve(hessian)
  if (pos==FALSE){
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
#' @param state A string with the name of a state for which to compute the probability (user-defined names). 
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
  state <- state_ord$ord[state_ord$state==state] #find state numbering
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
            llij = cubature::cubintegrate(integrand2, lower = lower,upper = upper, method = "divonne", maxEval = 500,
                                  tt = tmax[1], tt2=tmax[2],param = param, x = xval)$integral
            return(llij)
          })
          
        }else if (length(lower)>2){
          opi[j] = cubature::cubintegrate(integrand, lower = lower,upper = upper, method = "divonne", maxEval = 500,
                                tt = tmax[1], tt2=tmax[2],param = param, x = xval)$integral
        }
      }
    }
    op = op + opi
  }
  return(op)
}


#' Compute the occupancy probability over time with uncertainty bands
#'
#' The last state is often heavy to compute, but recall that is has to be equal to 1 minus the others
#'
#' @param state A string with the name of a state for which to compute the probability (user-defined names). 
#' @param time A vector of timepoints for which to compute the occupancy probability.
#' @param param  The parameter values in which the probabilities should be evaluated. The dimension
#' and ordering is given by the user-specified densities.
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param xval A vector of covariate values.
#' @param hessian The hessian matrix from the optimisation of the negative log-likelihood.
#' @param level The confidence level for the uncertainty bands, defaults to 0.95.
#' @return A vector of the same length as time.
#' 
occupancy_prob_ci_band <- function(state,time,param,graph,xval,hessian,level=0.95){
  est <- occupancy_prob(state=state,time=time,param=param,graph=graph,xval=xval)
  delta <- occupancy_prob_delta(state=state,time=time,param=param,graph=graph,xval=xval)
  varCov <- solve(hessian)
  kk <- length(time)
  lci <- rep(NA,kk)
  uci <- rep(NA,kk)
  zz <- stats::qnorm(0.5*(1-level),lower.tail = F)
  for (i in 1:kk){
    sei <- sqrt(delta[i,,drop=F]%*%varCov%*%t(delta[i,,drop=F]))
    lci[i] <- est[i]-zz*sei
    uci[i] <- est[i]+zz*sei
  }
  return(list(est=est,lower=lci,upper=uci))
}

#' Calculating partial derivatives of occupancy probability functions wrt parameters
#'
#' Helper function for computing uncertainty bands.
#'
#' @param state A string with the name of a state for which to compute the probability (user-defined names). 
#' @param time A vector of timepoints for which to compute the occupancy probability.
#' @param param  The parameter values in which the probabilities should be evaluated. The dimension
#' and ordering is given by the user-specified densities.
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param xval A vector of covariate values.
#' @return A vector of the same length as time.
#' 
occupancy_prob_delta = function(state, time, param, graph, xval = NULL){
  edge_mats = edge_matrices(graph)
  names_surv_dens = names_of_survival_density(graph)
  
  state_ord = state_ordering(graph)
  state <- state_ord$ord[state_ord$state==state] #find state numbering
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
          opi[j] = cubature::cubintegrate(integrand, lower = lower,upper = upper, method = "divonne", maxEval = 500,
                                tt = tmax[1], tt2=tmax[2],param = param, x = xval)$integral
        }
      }
    }
    op = op + opi
  }
  return(op)
}

#' Compute the overall survival with uncertainty bands
#'
#' ...
#'
#' @param time A vector of timepoints for which to compute the occupancy probability.
#' @param param  The parameter values in which the probabilities should be evaluated. The dimension
#' and ordering is given by the user-specified densities.
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param xval A vector of covariate values.
#' @param hessian The hessian matrix from the optimisation of the negative log-likelihood.
#' @param level The confidence level for the uncertainty bands, defaults to 0.95.
#' @return A vector of the same length as time.
#' 
overall_survival_ci_band <- function(time,param,graph,xval,hessian,level=0.95){
  est <- overall_survival(time=time,param=param,graph=graph,xval=xval)
  delta <- overall_survival_delta(time=time,param=param,graph=graph,xval=xval)
  varCov <- solve(hessian)
  kk <- length(time)
  lci <- rep(NA,kk)
  uci <- rep(NA,kk)
  zz <- stats::qnorm(0.5*(1-level),lower.tail = F)
  for (i in 1:kk){
    sei <- sqrt(delta[i,,drop=F]%*%varCov%*%t(delta[i,,drop=F]))
    lci[i] <- est[i]-zz*sei
    uci[i] <- est[i]+zz*sei
  }
  return(list(est=est,lower=lci,upper=uci))
}

#' Calculating partial derivatives of overall survival function wrt parameters
#'
#' Helper function for computing uncertainty bands.
#'
#' @param time A vector of timepoints for which to compute the occupancy probability.
#' @param param  The parameter values in which the probabilities should be evaluated. The dimension
#' and ordering is given by the user-specified densities.
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param xval A vector of covariate values.
#' @return A vector of the same length as time.
#' 
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
          opi[j,] = pracma::grad(cubint,x0=param,integrand=integrand,lower = lower,upper = upper, tmax=tmax,xval=xval)
        }
      }
    }
    op = op + opi
  }
  return(op)
}

#' Compute the transition probability over time
#'
#' The probability that a patient is found in state i at a time-point t given that she was in state j at 
#' a prior time-point v. Returns a vector of zeros if the transition from j to i is not possible (given the graph).
#' 
#' Only computes probabilities for direct transitions (a specific edge in the graph)
#'
#' @param trans_ji A string with two state names separated by "-" indicating the prior state j and the current state i (using user-defined state names). 
#' @param time_t A vector of timepoints denoting the "current" time.
#' @param time_v A number indicating the single prior time-point (must be smaller than all time_t).
#' @param param  The parameter values in which the probabilities should be evaluated. The dimension
#' and ordering is given by the user-specified densities.
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param xval A vector of covariate values.
#' @return A vector of the same length as time_t.
#' 
transition_prob = function(trans_ji, time_t,time_v, param, graph, xval = NULL){
  edge_mats = edge_matrices(graph)
  names_surv_dens = names_of_survival_density(graph)
  
  state_ord = state_ordering(graph)
  absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
  
  f_types = construct_formula_types(graph)
  state_j_name <- unlist(strsplit(trans_ji,"-"))[1]
  state_i_name <- unlist(strsplit(trans_ji,"-"))[2]
  state_j <- state_ord$order[state_ord$state==state_j_name]
  state_i <- state_ord$order[state_ord$state==state_i_name]
  
  f_types_match = f_types[which(substr(f_types,nchar(f_types)-1,nchar(f_types)) == paste(state_j,state_i,sep="") )]
  
  # Make a time-point format:
  dd <- data.frame(patient=c(rep(1,2)),time=c("v","t"),state=c(state_j_name,state_i_name))
  timepoints <- arrange_data(dd,graph)
  timepoints <- timepoints[1:(dim(timepoints)[2]-1)]
  #
  mm <- length(f_types_match) #0 or 1
  
  kk <- length(time_t)
  op <- rep(0,kk)
  if (mm==0){
    return(op)
  }else{
    # Denominator:
    dn <- occupancy_prob(state_j_name,time_v,param,graph,xval)
    
    # Numerator:
    f_type = f_types_match
    integrand = eval(parse(text=type_to_integrand(f_type,edge_mats, names_surv_dens,abs_exact=FALSE))) # always with abs_exact=FALSE
    
    choice <- (length(which(edge_mats$passedBy[f_type,]==max(edge_mats$traveled[f_type,])))>0)  #TRUE: if absorbing state is reached from a state where there was an option to go a different way
    
    if (substr(f_type,nchar(f_type),nchar(f_type)) %in% absorbing_states & choice==FALSE){ #if patient is observed in absorbing state
      dim_int <- nchar(f_type)-2  
    }else{ #if patient has not reached absorbing state, or has reached it from a node with a choice
      dim_int <- nchar(f_type)-1 
    }
    
    for (j in 1:kk){
      timepointsj <- timepoints
      timepointsj[which(timepointsj=="v")] <- time_v
      timepointsj[which(timepointsj=="t")] <- time_t[j]
      tt <- as.numeric(timepointsj)
      names(tt) <- names(timepointsj)
      
      lims <- finding_limits(tt,f_type,edge_mats,absorbing_states,abs_exact=F)
      lower <- lims$lower
      upper <- lims$upper
      tmax <- lims$tmax
      
      if(length(lower) == 0){
        op[j] = integrand(times=1,tt = tmax[1], tt2=tmax[2],param = param, x = xval) # the value of times does not matter
      }else if(length(lower)>0){
        if(length(lower)<=2 ){
          #opi[j] = repintegrate(integrand,tt = tmax[1], tt2=tmax[2],lower=lower,upper = upper, param = param, x = xval)
          
          op[j] = tryCatch({
            repintegrate(integrand,tt=tmax[1],tt2=tmax[2],lower=lower,upper = upper, param = param,
                         x = xval)
          },error=function(cond){
            integrand2 <- change_integrand(integrand)
            llij = cubature::cubintegrate(integrand2, lower = lower,upper = upper, method = "divonne", maxEval = 500,
                                          tt = tmax[1], tt2=tmax[2],param = param, x = xval)$integral
            return(llij)
          })
          
        }else if (length(lower)>2){
          op[j] = cubature::cubintegrate(integrand, lower = lower,upper = upper, method = "divonne", maxEval = 500,
                                          tt = tmax[1], tt2=tmax[2],param = param, x = xval)$integral
        }
      }
    }
    return(op/dn)
  }
}
#' Compute the transition probability over time with uncertainty bands
#'
#'
#' @param trans_ji A string with two state names separated by "-" indicating the prior state j and the current state i (using user-defined state names). 
#' @param time_t A vector of timepoints denoting the "current" time.
#' @param time_v A number indicating the single prior time-point (must be smaller than all time_t).
#' @param param  The parameter values in which the probabilities should be evaluated. The dimension
#' and ordering is given by the user-specified densities.
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param xval A vector of covariate values.
#' @param hessian The hessian matrix from the optimisation of the negative log-likelihood.
#' @param level The confidence level for the uncertainty bands, defaults to 0.95.
#' @return A vector of the same length as time.
#' 
transition_prob_ci_band <- function(trans_ji, time_t,time_v, param,graph,xval,hessian,level=0.95){
  est <- transition_prob(trans_ji=trans_ji,time_t=time_t,time_v=time_v,param=param,graph=graph,xval=xval)
  delta <- transition_prob_delta(trans_ji=trans_ji,time_t=time_t,time_v=time_v,param=param,graph=graph,xval=xval)
  varCov <- solve(hessian)
  kk <- length(time_t)
  lci <- rep(NA,kk)
  uci <- rep(NA,kk)
  zz <- stats::qnorm(0.5*(1-level),lower.tail = F)
  for (i in 1:kk){
    sei <- sqrt(delta[i,,drop=F]%*%varCov%*%t(delta[i,,drop=F]))
    lci[i] <- est[i]-zz*sei
    uci[i] <- est[i]+zz*sei
  }
  return(list(est=est,lower=lci,upper=uci))
}

#' Calculating partial derivatives of transition probability functions wrt parameters
#'
#' Helper function for computing uncertainty bands.
#'
#' @param trans_ji A string with two state names separated by "-" indicating the prior state j and the current state i (using user-defined state names). 
#' @param time_t A vector of timepoints denoting the "current" time.
#' @param time_v A number indicating the single prior time-point (must be smaller than all time_t).
#' @param param  The parameter values in which the probabilities should be evaluated. The dimension
#' and ordering is given by the user-specified densities.
#' @param graph A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).
#' @param xval A vector of covariate values.
#' @return A vector of the same length as time.
#' 
transition_prob_delta = function(trans_ji, time_t,time_v, param, graph, xval = NULL){
  edge_mats = edge_matrices(graph)
  names_surv_dens = names_of_survival_density(graph)
  
  state_ord = state_ordering(graph)
  absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
  
  f_types = construct_formula_types(graph)
  state_j_name <- unlist(strsplit(trans_ji,"-"))[1]
  state_i_name <- unlist(strsplit(trans_ji,"-"))[2]
  state_j <- state_ord$order[state_ord$state==state_j_name]
  state_i <- state_ord$order[state_ord$state==state_i_name]
  
  f_types_match = f_types[which(substr(f_types,nchar(f_types)-1,nchar(f_types)) == paste(state_j,state_i,sep="") )]
  
  # Make a time-point format:
  dd <- data.frame(patient=c(rep(1,2)),time=c("v","t"),state=c(state_j_name,state_i_name))
  timepoints <- arrange_data(dd,graph)
  timepoints <- timepoints[1:(dim(timepoints)[2]-1)]
  #
  mm <- length(f_types_match) #0 or 1
  
  kk <- length(time_t)
  n_deriv <- matrix(0,kk,length(param))
  if (mm==0){
    return(n_deriv)
  }else{
    # Denominator:
    dn <- occupancy_prob(state_j_name,time_v,param,graph,xval)
    dn_deriv <- occupancy_prob_delta(state_j_name,time_v,param,graph,xval)
    
    # Numerator:
    f_type = f_types_match
    integrand = eval(parse(text=type_to_integrand(f_type,edge_mats, names_surv_dens,abs_exact=FALSE))) # always with abs_exact=FALSE
    
    choice <- (length(which(edge_mats$passedBy[f_type,]==max(edge_mats$traveled[f_type,])))>0)  #TRUE: if absorbing state is reached from a state where there was an option to go a different way
    
    if (substr(f_type,nchar(f_type),nchar(f_type)) %in% absorbing_states & choice==FALSE){ #if patient is observed in absorbing state
      dim_int <- nchar(f_type)-2  
    }else{ #if patient has not reached absorbing state, or has reached it from a node with a choice
      dim_int <- nchar(f_type)-1 
    }
    
    for (j in 1:kk){
      timepointsj <- timepoints
      timepointsj[which(timepointsj=="v")] <- time_v
      timepointsj[which(timepointsj=="t")] <- time_t[j]
      tt <- as.numeric(timepointsj)
      names(tt) <- names(timepointsj)
      
      lims <- finding_limits(tt,f_type,edge_mats,absorbing_states,abs_exact=F)
      lower <- lims$lower
      upper <- lims$upper
      tmax <- lims$tmax
      
      if(length(lower) == 0){
        n_deriv[j,] = pracma::grad(integrand,x0=param,times=1,tt = tmax[1], tt2=tmax[2], x=xval)
      }else if(length(lower)>0){
        if(length(lower)<=2 ){
          #opi[j] = repintegrate(integrand,tt = tmax[1], tt2=tmax[2],lower=lower,upper = upper, param = param, x = xval)
          
          n_deriv[j,] = tryCatch({
            pracma::grad(repintegrate,x0=param,innerfunc=integrand,tt=tmax[1],tt2=tmax[2],lower=lower,upper = upper,x = xval)
          },error=function(cond){
            integrand2 <- change_integrand(integrand)
            llij = pracma::grad(cubint,x0=param,integrand=integrand2,lower = lower,upper = upper, tmax=tmax,xval=xval)
            return(llij)
          })
          
        }else if (length(lower)>2){
          n_deriv[j,] = pracma::grad(cubint,x0=param,integrand=integrand2,lower = lower,upper = upper, tmax=tmax,xval=xval)
        }
      }
    }
    nn <- transition_prob(trans_ji=trans_ji,time_t=time_t,time_v=time_v,param=param,graph=graph,xval=xval)*dn
    deriv <- (n_deriv*dn-nn%*%dn_deriv)/dn^2
    return(deriv)
  }
}


#' Helper function for differentiating cubintegrate functions
#'
#' @param integrand The function to integrate.
#' @param lower  Lower limits.
#' @param upper Upper limits.
#' @param tmax Time-points in which to evaluate the integrand.
#' @param param The parameter value in which to evaluate the integrand.
#' @param xval A vector of covariate values.
#' @return The value after integration.
#' 
cubint <- function(integrand,lower,upper,tmax,param,xval){
  cubature::cubintegrate(integrand, lower = lower,upper = upper, method = "divonne", maxEval = 500,
               tt = tmax[1], tt2=tmax[2],param = param, x = xval)$integral
}
