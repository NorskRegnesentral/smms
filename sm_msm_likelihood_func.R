### Multi-state functions part 1 - likelihood


## Chapter 3 - from a given formula type to an integrand
## 3.1 Connecting the edges with the survival and density functions
## Input: gg (the graph)
names_of_survial_density = function(gg){
  data_frame_possible_travels = get.data.frame(gg, what= "edges" )
  all_transitions = state_ordering(gg)
  data_frame_possible_travels
  all_transitions
  matrix_names = matrix(ncol = 5, nrow = nrow(data_frame_possible_travels))
  colnames(matrix_names) = c("survival_name", "density_name", "from_prev", "to_prev", "all_transitions")
  for(i in 1:nrow(data_frame_possible_travels)){
    kk_from = which(all_transitions[, "prev"] %in% data_frame_possible_travels[i, "from"]) 
    kk_to = which(all_transitions[, "prev"] %in% data_frame_possible_travels[i, "to"])
    matrix_names[i, "survival_name"] = paste(c("S_", all_transitions[kk_from, "new"], all_transitions[kk_to, "new"]), collapse = "")
    matrix_names[i, "density_name"] = paste(c("f_", all_transitions[kk_from, "new"], all_transitions[kk_to, "new"]), collapse = "")
    matrix_names[i, "from_prev"] = all_transitions[kk_from, "prev"]
    matrix_names[i, "to_prev"] = all_transitions[kk_to, "prev"]
    matrix_names[i, "all_transitions"] = paste(c(all_transitions[kk_from, "new"], all_transitions[kk_to, "new"]), collapse = "")
  }
  return(matrix_names)
}


## 3.2 Survival functions
## Input: gg (the graph), type (formula type), possible_passedBy (survival functions from possible/passedBy)
survival_function = function(gg, type, possible_passedBy){
  nrow = which(rownames(possible_passedBy) == type)
  all_possible_passedBy = which(possible_passedBy[nrow, ] != 0)
  all_colnames = names(all_possible_passedBy)
  if(length(all_colnames) > 0){
    which_transition = rep(NA, length(all_colnames))
    survival_function = rep(NA, length(all_colnames))
    for(i in 1:length(all_colnames)){
      which_transition[i] = which(names_of_survial_density(gg)[, "all_transitions"] == all_colnames[i])
      survival_function[i] = names_of_survial_density(gg)[which_transition[i], "survival_name"]
    }
  }
  else
    survival_function = "Error: no survival function for the chosen possible/passedBy for this type"
  return(survival_function)
}


## Density functions
## Input: gg (the graph), type (formula type), travelled (density functions from travelled)
density_function = function(gg, type, travelled){
  nrow = which(rownames(travelled) == type)
  all_travelled = which(travelled[nrow, ] != 0)
  all_colnames = names(all_travelled)
  if(length(all_colnames) > 0){
    which_transition = rep(NA, length(all_colnames))
    density_function = rep(NA, length(all_colnames))
    for(i in 1:length(all_colnames)){
      which_transition[i] = which(names_of_survial_density(gg)[, "all_transitions"] == all_colnames[i])
      density_function[i] = names_of_survial_density(gg)[which_transition[i], "density_name"]
    }
  }
  else
    density_function = "Error: no density function for travelled for this type"
  return(density_function)
}


# this version of the function is for datasets (or observations) where the time into absorbing state 
# is observed exactly
# for now it is limited to integrals of dimension 3 (but can be extended)
type_to_integrand_absExact = function(type,edge_mats,edge_abs){
  # inputs: type, edge_mats, (names of density functions, survival functions,) vector of edges into absorbing state
  # output: an integrand function written as a string per type
  
  travi <- edge_mats$traveled[type,]
  passi <- edge_mats$passedBy[type,]
  posi <- edge_mats$possible[type,]
  
  density_names <- c("f_01","f_03","f_12","f_13","f_23") #make as input, should be sorted in the same way as columns in edge mats
  surv_names <- c("S_01","S_03","S_12","S_13","S_23") #should be sorted in the same way as columns in edge mats
  
  # Write variable calls
  dim_integral <- max(travi)
  #if (names(travi)[which.max(travi)] %in% edge_abs) dim_integral <- dim_integral-1
  variable_names <- c("ss","uu","rr","vv")
  var_def <- rep(NA,dim_integral)
  for (i in 1:dim_integral){
    var_def[i] <- paste(variable_names[i],"<-","xx[",i,"]","\n",sep="")
  }
  
  # If the absorbing state is reached (and we know exactly when): change the last one
  if (names(travi)[which.max(travi)] %in% edge_abs){
    var_def[dim_integral+1] <- paste(variable_names[dim_integral],"<-","tt","\n",sep="")
  }
  
  # Write "travelled" functions
  variable_names_eval <- c("ss","uu-ss","rr-uu","vv-rr")
  n_trav <- sum(travi>0)
  f_trav <- rep(NA,n_trav)
  if (n_trav>0){
    for (i in 1:n_trav){
      f_trav[i] <- paste(density_names[which(travi==i)],"(param,",variable_names_eval[i],")","*")
    }
  }
  
  
  # Write "passedBy" functions
  n_passed <- sum(passi>0)
  S_passed <- rep(NA,n_passed)
  ii <- 1
  if (n_passed>0){
    for (i in 1:max(passi)){
      n_pass_i <- sum(passi==i)
      if (n_pass_i==0) next
      for (j in 1:n_pass_i){
        S_passed[ii] <- paste(surv_names[which(passi==i)][j],"(param,",variable_names_eval[i],")","*")
        ii <- ii+1
      }
    }
  }
  
  # Write "possible" functions
  n_posi <- sum(posi>0)
  S_posi <- rep(NA,n_posi)
  if (n_posi>0){
    ii <- 1
    if (dim_integral==0){
      eval <- "tt"
    }else{
      eval <- paste("tt-",variable_names[dim_integral])
    }
    
    for (i in 1:max(posi)){
      n_posi_i <- sum(posi==i)
      if (n_posi_i==0) next
      for (j in 1:n_posi_i){
        S_posi[ii] <- paste(surv_names[which(posi==i)][j],"(param,",eval,")","*")
        ii <- ii+1
      }
    }
    
  }
  
  # remove the last "*" in the last function
  if (length(S_passed)==0 & length(S_posi)==0){
    f_trav[n_trav] <- gsub('.{1}$', '', f_trav[n_trav]) 
  }else if (length(S_posi)==0){
    S_passed[n_passed] <- gsub('.{1}$', '', S_passed[n_passed])
  }else if (length(S_posi)!=0){
    S_posi[n_posi] <- gsub('.{1}$', '', S_posi[n_posi]) 
  }
  
  # paste everything together
  f_intro <- "function(xx,tt,param){"
  f_end <- "}"
  f_prod <- paste(paste(f_trav,collapse=""),paste(S_passed,collapse=""),paste(S_posi,collapse=""))
  ff <- paste(f_intro,paste(var_def,collapse=""),f_prod,f_end)
  
  return(ff)
}


finding_limits_integral = function(i, type, gg, all_edges, absorbing_state_new, all_data_set){
  ## Finding the absorbing state(s) 
  ## Possibly make a separate function to find the initial, transient and absorbing states 
  
  ## Arranging the data set and extracting the row. 
  time_points = all_data_set[i,1:(ncol(all_data_set)-1)]
  ## All the types
  all_time_points = substr(names(time_points), 2,2)
  splitted_type = unlist(strsplit(type, ""))
  ## Only include the time points which is found in the formula type
  all_included_time_points = c()
  for(j in 1:length(all_time_points)){
    if(all_time_points[j] %in% splitted_type){
      all_included_time_points = c(all_included_time_points, time_points[j])
    }
  }
  
  if(any(is.na(all_included_time_points))){
    included_time_points_na = which(is.na(all_included_time_points))
    unique_na = unique(substr(names(included_time_points_na), 2,2))
    upper_bound = rep(NA, length(unique_na))
    lower_bound = rep(NA, length(unique_na))
    for(j in 1:length(unique_na)){
      ## In case we are in a situation: 0(1)2(3)4
      col_numbers = which(substr(names(included_time_points_na), 2,2) == unique_na[j])
      upper_bound[j] = min(as.numeric(na.omit(all_included_time_points[(max(
        included_time_points_na[max(col_numbers)])+1):length(all_included_time_points)])))
      lower_bound[j] = min(as.numeric(na.omit(all_included_time_points[1:(min(included_time_points_na[min(col_numbers)])-1)])))
      which_to_replace = which(all_time_points == unique_na[j])
      all_included_time_points[which_to_replace[1]] = upper_bound[j]
      all_included_time_points[which_to_replace[2]] = lower_bound[j]
    }
  }
  all_included_time_points = as.numeric(unname(all_included_time_points))
  
  tmin_int1 = NA; tmax_int1 = NA; tmin_int2 = NA; tmax_int2 = NA 
  tmin_int3 = NA; tmax_int3 = NA; tmin_int4 = NA; tmax_int4 = NA
  tmax = max(all_included_time_points)
  
  for(j in 1:length(splitted_type)){
    if(!(splitted_type[length(splitted_type)] %in% absorbing_state_new)){
      if(j == 2){
        tmin_int1 = all_included_time_points[1]
        tmax_int1 = all_included_time_points[2]
      } else if(j == 3){
        tmin_int2 = all_included_time_points[3]
        tmax_int2 = all_included_time_points[4]
      } else if(j == 4){
        tmin_int3 = all_included_time_points[5]
        tmax_int3 = all_included_time_points[6]
      } else if(j == 5){
        tmin_int4 = all_included_time_points[7]
        tmax_int4 = all_included_time_points[8]
      }
    } else if((splitted_type[length(splitted_type)] %in% absorbing_state_new)){
      if(j == 3){
        tmin_int1 = all_included_time_points[1]
        tmax_int1 = all_included_time_points[2]
      } else if(j == 4){
        tmin_int2 = all_included_time_points[3]
        tmax_int2 = all_included_time_points[4]
      } else if(j == 5){
        tmin_int3 = all_included_time_points[5]
        tmax_int3 = all_included_time_points[6]
      } else if(j == 6){
        tmin_int4 = all_included_time_points[7]
        tmax_int4 = all_included_time_points[8]
      }
    }
  }
  all_limits = c(tmin_int1, tmax_int1, tmin_int2, tmax_int2, tmin_int3, tmax_int3, 
                 tmin_int4, tmax_int4, tmax)
  names(all_limits) = c("tmin_int1", "tmax_int1", "tmin_int2", "tmax_int2", "tmin_int3", "tmax_int3", 
                        "tmin_int4", "tmax_int4", "tmax")
  return(all_limits)
}


## Make additional input for cubintegrate
from_time_point_to_integral = function(param, method1 = "hcubature", integrand = integrand, all_data_set = all_data_set, 
                                       all_integral_limits = all_integral_limits, mc_cores = 2){
  ## Including all types and the data set
  final_integral = sum(unlist(mclapply(1:nrow(all_data_set), function(i){
    ## Using the observation type to find all possible formula types
    ## If we only have one type for this observation type
    if(length(all_integral_limits[[i]]) == 1){
      lower_integral_mellomregn = c(as.numeric(all_integral_limits[[i]][[1]][1]), as.numeric(all_integral_limits[[i]][[1]][3]),
                                    as.numeric(all_integral_limits[[i]][[1]][5]), as.numeric(all_integral_limits[[i]][[1]][7]))
      lower_integral = lower_integral_mellomregn[!is.na(lower_integral_mellomregn)]
      upper_integral_mellomregn = c(as.numeric(all_integral_limits[[i]][[1]][2]), as.numeric(all_integral_limits[[i]][[1]][4]),
                                    as.numeric(all_integral_limits[[i]][[1]][6]), as.numeric(all_integral_limits[[i]][[1]][8]))
      upper_integral = upper_integral_mellomregn[!is.na(upper_integral_mellomregn)]
      tmax = as.numeric(all_integral_limits[[i]][[1]][9])
      
      length_lower_integral = length(lower_integral)
      length_unique_lower_integral = length(unique(lower_integral))
      
      if(length_lower_integral == 0){
        log(integrand[[i]][[1]](xx=NULL,tt = tmax, param = param))
      } else if(length_unique_lower_integral != length_lower_integral){
        log(cubintegrate(integrand[[i]][[1]], lower = lower_integral, upper = upper_integral, method = "vegas",tt = tmax, param = param, maxEval = 10^3)$integral)
      } else if(length_unique_lower_integral == length_lower_integral){
        log(cubintegrate(integrand[[i]][[1]], lower = lower_integral, upper = upper_integral, method = method1, tt = tmax, param = param)$integral)
      }
    } else{ #if(length(all_integral_limits[[i]]) > 1){ ## If we have more than one type for this observation type
      calculate_integral_multiple_types = rep(NA, length(all_integral_limits[[i]]))
      for(j in 1:length(all_integral_limits[[i]])){
        lower_integral_mellomregn = c(as.numeric(all_integral_limits[[i]][[j]][1]), as.numeric(all_integral_limits[[i]][[j]][3]), 
                                      as.numeric(all_integral_limits[[i]][[j]][5]), as.numeric(all_integral_limits[[i]][[j]][7]))
        lower_integral = lower_integral_mellomregn[!is.na(lower_integral_mellomregn)]
        upper_integral_mellomregn = c(as.numeric(all_integral_limits[[i]][[j]][2]), as.numeric(all_integral_limits[[i]][[j]][4]),
                                      as.numeric(all_integral_limits[[i]][[j]][6]), as.numeric(all_integral_limits[[i]][[j]][8]))
        upper_integral = upper_integral_mellomregn[!is.na(upper_integral_mellomregn)]
        tmax = as.numeric(all_integral_limits[[i]][[j]][9])
        
        length_lower_integral = length(lower_integral)
        length_unique_lower_integral = length(unique(lower_integral))
        
        if(length_lower_integral == 0){
          calculate_integral_multiple_types[j] = integrand[[i]][[j]](xx=NULL,tt = tmax, param = param) # the value of xx does not matter
        } else if(length_unique_lower_integral != length_lower_integral){
          calculate_integral_multiple_types[j] = cubintegrate(integrand[[i]][[j]], lower = lower_integral, upper = upper_integral, method = "vegas", maxEval = 1000,
                                                              tt = tmax, param = param)$integral
          # needs to check that it does not matter that the dimension of xx is higher then lower and upper
        } else if(length_unique_lower_integral == length_lower_integral){
          calculate_integral_multiple_types[j] = cubintegrate(integrand[[i]][[j]], lower = lower_integral, upper = upper_integral, 
                                                              method = method1, tt = tmax, param = param)$integral
        } 
      }
      ## Taking the log of the sum of the likelihood contribution for each possible type
      log(sum(calculate_integral_multiple_types))
    }
  }, mc.cores = mc_cores)))
  return(final_integral)
}
