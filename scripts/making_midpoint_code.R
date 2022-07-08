devtools::load_all()
library(data.table)
## This is a work-in-progress script (not been tested) which automatically
## finds the midpoint.
## Should only be used to find the startparameters to smms. 

arrange_data_midtpoint = function(data, graph, abs_exact = TRUE){
  set.seed(6789)
  all_data_set = arrange_data(data, graph)
  state_ord = state_ordering(graph)
  init <- sort(state_ord$order[which(state_ord$type=="init")])
  trans <- sort(state_ord$order[which(state_ord$type=="trans")])
  abs <- sort(state_ord$order[which(state_ord$type=="abs")])
  inds = unique(data$patient)
  nn = length(inds)
  max_depth = max(nchar(construct_formula_types(gg)))
  
  timepoints_midpoint = as.data.frame(matrix(NA,nn,max_depth-1))
  timepoints_midpoint = cbind(timepoints_midpoint, all_data_set$obs_type)
  names_midpoint = paste(rep("mid", max_depth-1),
                        1:(max_depth-1), sep="")
  colnames(timepoints_midpoint) = c(names_midpoint, "obs_type")
  
  ## Filling out the timepoints matrix
  observation_type = rep(NA, nrow(all_data_set))
  for (i in 1:nn){
    observation_type[i] = all_data_set[i,"obs_type"]
    formula_obs_types = all_types(graph)
    type = names(which(formula_obs_types[, observation_type[i]] == 1))
    time_points = all_data_set[i,1:(ncol(all_data_set)-1)]
    if(!(abs %in% unlist(strsplit(observation_type[i], "")))){
      ## if a patient stays in a transient state
      excluded_timepoint = c()
      included_names = c(init)
      for(j in 1:length(trans)){
        id_dupl <- which(as.character(c(init, rep(trans, each = 2),abs)) %in% trans[j])
        if(!(max(id_dupl) == max(which(!is.na(time_points))))){
          excluded_timepoint = c(excluded_timepoint, max(id_dupl))
          included_names = c(included_names, trans[j])
        } else{
          included_names = c(included_names, rep(trans[j],2))
        }
      }
      timepoints_mid_length = cbind(time_points[-excluded_timepoint])
      colnames(timepoints_mid_length) = c(included_names, abs)
    } else{  
      excluded_timepoint = c()
      for(j in 1:length(trans)){
        id_dupl <- which(as.character(c(init, rep(trans, each = 2),abs)) %in% trans[j])
        excluded_timepoint = c(excluded_timepoint, max(id_dupl))
      }
      timepoints_mid_length = time_points[-excluded_timepoint]
      colnames(timepoints_mid_length) = c(init, trans, abs)
    }
    if(length(type != 1)){
      ## Draw if we have more than one possuble type
      type_sample = sample(1:length(type), 1)
      type_new = type[type_sample]
      type = type_new
    }
    if(type %in% init){
      ## Only in transient state
        timepoints_midpoint[i, 1] = time_points[1]
      } else if(observation_type[i] == type){
        ## Exclude all NA and the double timepoints
        na_omit_intergal_limits = time_points[!is.na(time_points)]
        timepoints_mid_length = na_omit_intergal_limits[-excluded_timepoint]
        for(j in 1:(length(timepoints_mid_length)-1)){
          ## If stays in a state or is observed in absorbing exact
          if(j == (length(timepoints_mid_length)-1) & j != 1 & (abs_exact == TRUE | !(abs %in% unlist(strsplit(observation_type[i], ""))))){
            timepoints_midpoint[i, j] = timepoints_mid_length[j+1] - timepoints_midpoint[i, j-1]
          } else if (j == (length(timepoints_mid_length)-1) & j == 1 & abs_exact == TRUE){ ## If observed exact in absorbing (direct)
            timepoints_midpoint[i, j] = timepoints_mid_length[j+1]
          } else if((j == 1 & j != (length(timepoints_mid_length)-1)) | (j == (length(timepoints_mid_length)-1) & j == 1 & abs_exact == FALSE)){ ## First midpoint 
            timepoints_midpoint[i, j] = (timepoints_mid_length[j] + timepoints_mid_length[j+1])/2   
          }else if((j != 1 & j != (length(timepoints_mid_length)-1)) | (j == (length(timepoints_mid_length)-1) & j != 1 & abs_exact == FALSE)){ ## Second midpoint
            timepoints_midpoint[i, j] = (-timepoints_midpoint[i, j-1] + timepoints_mid_length[j+1])/2   
          }
        }
      } else if(observation_type[i] != type){
        timepoints_midpoint$obs_type[i] = type
        if(length(which(!(colnames(timepoints_mid_length) %in% unlist(strsplit(type, ""))))) != 0){
          rr = which(!(colnames(timepoints_mid_length) %in% unlist(strsplit(type, ""))))
          timepoints_mid_length = timepoints_mid_length[-rr]
        }
        which_na = which(is.na(timepoints_mid_length))
        which_not_na = which(!is.na(timepoints_mid_length))
        between_1 = matrix(NA, nrow = length(which_na), ncol = 2)
        for(j in 1:length(which_na)){
          between_1[j, 1] = which_not_na[max(which(which_na[j] > which_not_na))]
          between_1[j, 2] = which_not_na[min(which(which_na[j] < which_not_na))]
        }
        between_1 = as.data.table(between_1)
        names(between_1) = c("lower", "upper")
        between_1[, count := .N, by = c("lower", "upper")]
        between_1[, mm := .GRP, by = c("lower", "upper")]
        unique_between = unique(between_1)
        for(j in 1:nrow(unique_between)){
          sample_1 = sort(sample(1:100,unique_between[j, count]))
          sample_from = seq(from = unlist(timepoints_mid_length[unique_between[j, lower]]), 
                            to = unlist(timepoints_mid_length[unique_between[j, upper]]), length.out = 100)
          for(k in 1:length(sample_1)){
            timepoints_mid_length[which_na[k]] = sample_from[sample_1[k]]
          }
        }
        for(j in 1:(length(timepoints_mid_length)-1)){
          ## If stays in a state or is observed in absorbing exact
          if(j == (length(timepoints_mid_length)-1) & j != 1 & (abs_exact == TRUE | !(abs %in% unlist(strsplit(observation_type[i], ""))))){
            timepoints_midpoint[i, j] = timepoints_mid_length[j+1] - timepoints_midpoint[i, j-1]
          } else if (j == (length(timepoints_mid_length)-1) & j == 1 & abs_exact == TRUE){ ## If observed exact in absorbing (direct)
            timepoints_midpoint[i, j] = timepoints_mid_length[j+1]
          } else if((j == 1 & j != (length(timepoints_mid_length)-1)) | (j == (length(timepoints_mid_length)-1) & j == 1 & abs_exact == FALSE)){ ## First midpoint 
            timepoints_midpoint[i, j] = (timepoints_mid_length[j] + timepoints_mid_length[j+1])/2   
          }else if((j != 1 & j != (length(timepoints_mid_length)-1)) | (j == (length(timepoints_mid_length)-1) & j != 1 & abs_exact == FALSE)){ ## Second midpoint
            timepoints_midpoint[i, j] = (-timepoints_midpoint[i, j-1] + timepoints_mid_length[j+1])/2   
          }
        }
      }
    }  
  return(timepoints_midpoint)
}

type_to_integrand_midpoint = function(form_type,edge_mats,names_surv_dens,abs_exact=TRUE){
  travi <- edge_mats$traveled[form_type,,drop=F]
  passi <- edge_mats$passedBy[form_type,,drop=F]
  posi <- edge_mats$possible[form_type,,drop=F]
  
  density_names <- names_surv_dens$density_name
  surv_names <- names_surv_dens$survival_name
  edge_abs <- names_surv_dens$edge_name[names_surv_dens$type=="abs"]
  
  # Determine the integral dimension
  dim_integral <- max(c(travi, posi, passi))
  intoAbs <- colnames(travi)[which.max(travi)] %in% edge_abs
  choice <- (length(which(passi==max(travi)))>0)   #TRUE: if absorbing state is reached from a state where there was an option to go a different way
  
  variable_names <- c("ss","uu","rr","vv","ww","tt")
  var_def <- rep(NA,dim_integral)
    # Write variable calls
    for (i in 1:dim_integral){
      var_def[i] <- paste(variable_names[i],"<-","times[",i,"]","\n",sep="")
    }
    
    # Write "travelled" functions
    variable_names_eval <- c("ss","uu","rr","vv","ww","tt")
    n_trav <- sum(travi>0)
    
    f_trav <- rep(NA,n_trav)
    if (n_trav>0){
      for (i in 1:n_trav){
        f_trav[i] <- paste(density_names[which(travi==i)],"(param,x,",variable_names_eval[i],")","*")
      }
    }
    if (intoAbs==TRUE & abs_exact==FALSE & choice==FALSE){
      f_trav[n_trav] <- paste("(",surv_names[which(travi==n_trav)],"(param,x,",variable_names_eval2[n_trav],")","-",
                              surv_names[which(travi==n_trav)],"(param,x,",variable_names_eval[n_trav],"))*")
      # this solution requires that if tt2=-1 (which is a value it should be given when not relevant), 
      # the survival function S_23 will return 1
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
          S_passed[ii] <- paste(surv_names[which(passi==i)][j],"(param,x,",variable_names_eval[i],")","*")
          ii <- ii+1
        }
      }
    }
    
    # Write "possible" functions
    n_posi <- sum(posi>0)
    S_posi <- rep(NA,n_posi)
    if (n_posi>0){
      ii <- 1
      for (i in 1:max(posi)){
        n_posi_i <- sum(posi==i)
        if (n_posi_i==0) next
        for (j in 1:n_posi_i){
          S_posi[ii] <- paste(surv_names[which(posi==i)][j],"(param,x,",variable_names_eval[i],")","*")
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
    f_intro <- "function(times,param,x,...){"
    if (intoAbs==TRUE & abs_exact==FALSE & choice==FALSE){
      f_intro <- "function(times,param,x,tt2=-1){" #perhaps tt2=NULL?
    } 
    f_end <- "}"
    f_prod <- paste(paste(f_trav,collapse=""),paste(S_passed,collapse=""),paste(S_posi,collapse=""))
    ff <- paste(f_intro,paste(var_def,collapse=""),f_prod,f_end)
  return(ff)
}

loglikelihood_midpoint = function(param, integrand, limits, X = NULL){
  vv = rep(NA, length(integrand))
  for(i in 1:length(integrand)){
    vv[i] = log(integrand[[i]](limits[[i]], param = param, x = X))
  }
  return(-sum(vv))
}

optim_midpoint = function(startval, data, graph, X = NULL){
  data_mid = arrange_data_midtpoint(data, graph)
  observation_type = rep(NA, nrow(data_mid))
  limits = data_mid[-ncol(data_mid)]
  edge_mats <- edge_matrices(gg)
  names_surv_dens = names_of_survival_density(gg)
  integrand = list()
  limit2 = list()
  for(i in 1:nrow(data_mid)){
    observation_type[i] = as.character(data_mid[i,]$obs_type)
  f_types = observation_type[i]
  limit1 = limits[i,]
  limit2[[i]] = as.vector(unlist(limit1[which(!is.na(limit1))]))
  integrand[[i]] = eval(parse(text=type_to_integrand_midpoint(f_types, edge_mats, names_surv_dens)))
  }
  t = nlminb(startval, loglikelihood_midpoint, integrand = integrand, limits = limit2, X = X)
  startval = t$par
  return(startval)
}




