### Multi-state functions part 2 - likelihood

#' Produce a table with each edge of the graph and the corresponding names of survival functions 
#' and densities
#'
#' @param graph A directed, acyclic graph in the igraph format (igraph package).
#' @return A data frame with one row for each edge in the graph. First column gives the standard edge
#' name, second column indicates the names of the survival functions, third column gives the names of the 
#' density functions. The two next columns give the names of the source node and the target node *in the
#' original user-defined names*. The last column gives the type of edge: either "abs" for absorbing if the target node
#' is an absorbing state, or "trans" for transient if not.
#' @export
names_of_survival_density = function(graph){
  all_edges = igraph::get.edgelist(graph)
  # Update with state ordering as node names:
  state_ord = state_ordering(graph)
  all_edges_new=all_edges
  all_edges_new[,1] <- state_ord$order[match(all_edges[,1],state_ord$state)]
  all_edges_new[,2] <- state_ord$order[match(all_edges[,2],state_ord$state)]
  edge_names = apply(all_edges_new,1,paste,collapse="")
  
  matrix_names = as.data.frame(matrix(ncol = 6, nrow = length(edge_names)))
  colnames(matrix_names) = c("edge_name","survival_name", "density_name", "from_prev", "to_prev", "type")
  for(i in 1:length(edge_names)){
    kk_to = which(state_ord$state %in% all_edges[i,2])
    matrix_names[i, "survival_name"] = paste(c("S_", edge_names[i]), collapse = "")
    matrix_names[i, "density_name"] = paste(c("f_", edge_names[i]), collapse = "")
    matrix_names[i, "from_prev"] = all_edges[i,1]
    matrix_names[i, "to_prev"] = all_edges[i,2]
    matrix_names[i, "edge_name"] = edge_names[i]
    matrix_names[i, "type"] = state_ord[kk_to, "type"]
  }
  return(matrix_names)
}

#' Write out the integrand as a string
#'
#' For a given formula type the appropriate integrand is generated as a text string. The integrand
#' takes the form of an R type function. For now this function is limited to integrals up to dimension 5.
#'
#' @param form_type A formula type (as a string).
#' @param edge_mats A list with 3 matrices produced by the edge_matrices() function.
#' @param names_surv_dens A data frame indicating the names of survival functions and densities. 
#' The output of the names_of_survival_density() function.
#' @param abs_exact A boolean indicating whether the time of entrance into absorbing states is observed
#' exactly (TRUE) or not (FALSE). Default value is TRUE.
#' @return A text string giving the R syntax of the integrand.
#' @export
type_to_integrand = function(form_type,edge_mats,names_surv_dens,abs_exact=TRUE){
  
  travi <- edge_mats$traveled[form_type,,drop=F]
  passi <- edge_mats$passedBy[form_type,,drop=F]
  posi <- edge_mats$possible[form_type,,drop=F]

  density_names <- names_surv_dens$density_name
  surv_names <- names_surv_dens$survival_name
  edge_abs <- names_surv_dens$edge_name[names_surv_dens$type=="abs"]
  
  # Determine the integral dimension
  dim_integral <- max(travi)
  if (max(travi)>0){
    intoAbs <- colnames(travi)[which.max(travi)] %in% edge_abs
  }else{
    intoAbs <- FALSE
  }
  choice <- (length(which(passi==max(travi)))>0)   #TRUE: if absorbing state is reached from a state where there was an option to go a different way
  if (intoAbs){
    if(abs_exact==TRUE) dim_integral <- dim_integral-1
    if (abs_exact==FALSE & choice==FALSE) dim_integral <- dim_integral-1
  } 
  
  variable_names <- c("ss","uu","rr","vv","ww")
  var_def <- rep(NA,dim_integral)
  
  if (dim_integral<=2 & dim_integral>0){  # here we write the integrand in the form which suits the integrate() function
    
    # Write variable calls
    for (i in 1:dim_integral){
      var_def[i] <- paste(variable_names[i],",")
    }
    
    if (intoAbs==TRUE & abs_exact==FALSE & choice==FALSE) var_def <- c(var_def,"tt2=-1 ,") #perhaps tt2=NULL?
    
    
    # Write "travelled" functions
    variable_names_eval <- c("ss","uu") #c("ss","uu","rr","vv","ww")
    n_trav <- sum(travi>0)
    if (intoAbs & abs_exact==TRUE){
      variable_names_eval[n_trav] <-  paste("tt-", paste0(variable_names_eval[1:(n_trav-1)],collapse="",sep="- "))
      variable_names_eval[n_trav] <- gsub('.{2}$', '', variable_names_eval[n_trav]) 
    }
    if (intoAbs==TRUE & abs_exact==FALSE & choice==FALSE){
      variable_names_eval[n_trav] <-  paste("tt-", paste0(variable_names_eval[1:(n_trav-1)],collapse="",sep="- "))
      variable_names_eval[n_trav] <- gsub('.{2}$', '', variable_names_eval[n_trav]) 
      
      variable_names_eval2 <- c("ss","uu")
      variable_names_eval2[n_trav] <-  paste("tt2-", paste0(variable_names_eval[1:(n_trav-1)],collapse="",sep="- "))
      variable_names_eval2[n_trav] <- gsub('.{2}$', '', variable_names_eval2[n_trav]) 
    }
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
    variable_names_eval_pos <- c("tt","tt-ss","tt-ss-uu") #,"tt-rr-uu-ss","tt-rr-uu-ss-vv")
    n_posi <- sum(posi>0)
    S_posi <- rep(NA,n_posi)
    if (n_posi>0){
      ii <- 1
      eval <- variable_names_eval_pos[dim_integral+1]
      
      for (i in 1:max(posi)){
        n_posi_i <- sum(posi==i)
        if (n_posi_i==0) next
        for (j in 1:n_posi_i){
          S_posi[ii] <- paste(surv_names[which(posi==i)][j],"(param,x,",eval,")","*")
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
    f_intro <- paste("function(",paste(var_def,collapse=""),"tt,param,x,...){")
    f_end <- "}"
    f_prod <- paste(paste(f_trav,collapse=""),paste(S_passed,collapse=""),paste(S_posi,collapse=""))
    ff <- paste(f_intro,f_prod,f_end)
    
  }else{ # here we write the integrand in the form suitable for cubintegrate()
    
    # Write variable calls
    for (i in 1:dim_integral){
      var_def[i] <- paste(variable_names[i],"<-","times[",i,"]","\n",sep="")
    }
    
    # Write "travelled" functions
    variable_names_eval <- c("ss","uu-ss","rr-uu","vv-rr","ww-vv","tt-ww")
    n_trav <- sum(travi>0)
    if (intoAbs & abs_exact==TRUE){
      variable_names_eval[n_trav] <- sub("^.{2}", "tt", variable_names_eval[n_trav])
    }
    if (intoAbs==TRUE & abs_exact==FALSE & choice==FALSE){
      variable_names_eval[n_trav] <- sub("^.{2}", "tt", variable_names_eval[n_trav]) 
      
      variable_names_eval2 <- c("ss","uu-ss","rr-uu","vv-rr","ww-vv","tt-ww")
      variable_names_eval2[n_trav] <- sub("^.{2}", "tt2", variable_names_eval[n_trav]) 
    }

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
      if (dim_integral==0){
        eval <- "tt"
      }else{
        eval <- paste("tt-",variable_names[dim_integral],sep="")
      }
      
      for (i in 1:max(posi)){
        n_posi_i <- sum(posi==i)
        if (n_posi_i==0) next
        for (j in 1:n_posi_i){
          S_posi[ii] <- paste(surv_names[which(posi==i)][j],"(param,x,",eval,")","*")
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
    f_intro <- "function(times,tt,param,x,...){"
    if (intoAbs==TRUE & abs_exact==FALSE & choice==FALSE){
      f_intro <- "function(times,tt,param,x,tt2=-1){" #perhaps tt2=NULL?
    } 
    f_end <- "}"
    f_prod <- paste(paste(f_trav,collapse=""),paste(S_passed,collapse=""),paste(S_posi,collapse=""))
    ff <- paste(f_intro,paste(var_def,collapse=""),f_prod,f_end)
  }
  return(ff)
}



#' Integrate functions of dimension 2
#'
#' The following two functions are helper functions for integrating the inner function
#' in an integral of dimension 2.
#'
#' @param ss The variable to integrate over in the outer integral.
#' @param innerfunc The integrand (an R function).
#' @param tt An input to innerfunc (a time-point in which to evaluate the integrand).
#' @param param The parameter value in which to evaluate the integrand.
#' @param x A vector of covariate values.
#' @param lower2  Lower limits.
#' @param upper2 Upper limits.
#' @param tt2 An input to innerfunc (a time-point in which to evaluate the integrand).
#' @param ... further arguments to innerfunc.
#' @return The value of the integral.
#'
repint2 <- function(ss,innerfunc,tt,param,lower2,upper2,x,tt2=-1,...){ #integrate over uu
  mm <- length(ss)
  out <- rep(NA,mm)
  for (i in 1:mm){
    out[i] <- stats::integrate(innerfunc,lower=max(lower2-ss[i],0),upper=upper2-ss[i],tt=tt,x=x,
                        ss=ss[i],param=param,tt2=tt2)$value
  }
  return(out)
} 

#' Integrate functions of dimension 1 
#'
#' The following is a helper function for integrating in dimension 1.
#'
#' @param innerfunc The integrand (an R function).
#' @param tt An input to innerfunc (a time-point in which to evaluate the integrand).
#' @param param The parameter value in which to evaluate the integrand.
#' @param x A vector of covariate values.
#' @param lower  Lower limits.
#' @param upper Upper limits.
#' @param tt2 An input to innerfunc (a time-point in which to evaluate the integrand).
#' @param ... further arguments to innerfunc.
#' @return The value of the integral.
#'
repintegrate <- function(innerfunc,tt,param,x,lower,upper,tt2=-1,...){ #integrate over ss
  if (length(lower)==1){
    out <- stats::integrate(innerfunc,lower=lower,upper=upper,tt=tt,param=param,x=x,tt2=tt2)$value
  }else{
    out <- stats::integrate(repint2,innerfunc=innerfunc,lower=max(lower[1],0),upper=upper[1],tt=tt,param=param,x=x,
                     lower2=lower[2],upper2=upper[2],tt2=tt2)$value
  }
  return(out)
} 

#' Change integrand 
#'
#' Change integrand syntax from the form appropriate for repintegrate() to the form appropriate 
#' for cubintegrate().
#'
#' @param integr An R function from eval-parse and type_to_integrand().
#' @return An R function
#' 
change_integrand <- function(integr){
  int <- paste0(deparse(integr), collapse = " ")
  m <- regexpr("ss,.+tt,",int)
  call <- regmatches(int,m)
  
  if (grepl("tt2",call)){
    int <- sub("ss,.+tt,","times,tt2=-1,tt,",int)
    if (grepl("uu",call)){
      int <- sub("\\{.+?f_","\\{\nss<-times[1]\nuu<-times[2]\nf_",int) #sub("\\{.+f_01","\\{\nss<-times[1]\nuu<-times[2]\nf_01",int)
      int <- gsub(", uu",", uu-ss",int)
      int <- gsub("tt - uu - ss","tt- uu",int)
      int <- gsub("tt - ss - uu","tt- uu",int)
      int <- gsub("tt2 - uu - ss","tt2- uu",int)
      int <- gsub("tt2 - ss - uu","tt2- uu",int)
    }else{
      int <- sub("\\{.+?f_","\\{\nss<-times[1]\nf_",int) #sub("\\{.+f_01","\\{\nss<-times[1]\nf_01",int)
    }
  }else{
    int <- sub("ss,.+tt,","times,tt,",int)
    if (grepl("uu",call)){
      int <- sub("\\{.+?f_","\\{\nss<-times[1]\nuu<-times[2]\nf_",int) #sub("\\{.+f_01","\\{\nss<-times[1]\nuu<-times[2]\nf_01",int)
      int <- gsub(", uu",", uu-ss",int)
      int <- gsub("tt - uu - ss","tt- uu",int)
      int <- gsub("tt - ss - uu","tt- uu",int)
    }else{
      int <- sub("\\{.+?f_","\\{\nss<-times[1]\nf_",int) #sub("\\{.+f_","\\{\nss<-times[1]\nf_01",int)
    }
  }
  return(eval(parse(text=int)))
}


#' Find the limits of integration 
#'
#' For a given formula type, a vector of (relevant) timepoints is sorted into the appropriate 
#' integral limits.
#'
#' @param timepoints A named vector of relevant timepoints (belonging to a patient), i.e. a single
#' row from the output of the arrange_data() function.
#' @param form_type A formula type (as a string).
#' @param edge_mats A list with 3 matrices produced by the edge_matrices() function.
#' @param absorbing_states A vector denoting the absorbing states (in the ordered naming system).
#' @param abs_exact A boolean indicating whether the time of entrance into absorbing states is observed
#' exactly (TRUE) or not (FALSE). Default value is TRUE.
#' @return A list with "lower" a vector of the lower limits of integration, "upper" a vector of the 
#' upper limits of integration, and "tmax" the last timepoint (in which the integrand often needs 
#' to be evaluated).
#' @export
finding_limits <- function(timepoints,form_type,edge_mats,absorbing_states,abs_exact=TRUE){
  splitted_f_type = unlist(strsplit(form_type, ""))
  
  #unobs_states <- splitted_f_type[which(!(splitted_f_type %in% unlist(strsplit(obs_type,""))))]
  
  choice <- (length(which(edge_mats$passedBy[form_type,]==max(edge_mats$traveled[form_type,])))>0)  #TRUE: if absorbing state is reached from a state where there was an option to go a different way
  intoAbs <- splitted_f_type[length(splitted_f_type)] %in% absorbing_states
  seclast <- timepoints[which(substr(names(timepoints),2,2)%in%splitted_f_type[length(splitted_f_type)-1])]
  
  if (intoAbs==TRUE){ #if patient is observed in absorbing state
    if (abs_exact==TRUE){
      dim_int <- length(splitted_f_type)-2 
      tmax = c(max(timepoints,na.rm=T),-1)
    }else{
      if (choice==TRUE){
        dim_int <- length(splitted_f_type)-1
        tmax = c(max(timepoints,na.rm=T),-1)
      }else{
        dim_int <- length(splitted_f_type)-2
        if(is.na(seclast[length(seclast)])){
          tmax = c(max(timepoints,na.rm=T),-1)
        }else{
          tmax = c(max(timepoints,na.rm=T),unlist(seclast[length(seclast)]))
        }
      }
    }
  }else{ #if patient has not reached absorbing state
    dim_int <- length(splitted_f_type)-1 
    tmax = c(max(timepoints,na.rm=T),-1)
  }
  
  if (dim_int==0){
    limit <- list(lower=NULL,upper=NULL,tmax=tmax)
  }else{
    # Lower limits are the d first "M" timepoints matching the formula type.
    # NAs are filled first with the previous values.
    M_times <- timepoints[grep("M",names(timepoints))]
    M_times <- M_times[which(substr(names(M_times),2,2)%in%splitted_f_type)]
    #id_na <- which(is.na(M_times) & substr(names(M_times),2,2)%in%unobs_states) #in unobs unecessary?
    id_na <- which(is.na(M_times))
    if (length(id_na)>0){
      for (j in 1:length(id_na)){  ## CHECK
        M_times[id_na[j]] <- M_times[id_na[j]-1]
      }
    }
    lower <- M_times[1:dim_int]
    #id_na <- which(is.na(lower))
    #if (length(id_na)>0) lower[id_na] <- lower[id_na-1]
    
    # Upper limits are the d first "m" timepoints matching the formula type.
    # NAs are filled first with the next values.
    m_times <- timepoints[grep("m",names(timepoints))]
    m_times <- m_times[which(substr(names(m_times),2,2)%in%splitted_f_type)]
    #id_na <- which(is.na(m_times) & substr(names(m_times),2,2)%in%unobs_states)
    id_na <- which(is.na(m_times))
    if (length(id_na)>0){
      for (j in length(id_na):1){  ## CHECK
        m_times[id_na[j]] <- m_times[id_na[j]+1]
      }
    }
    upper <- m_times[1:dim_int]
    #if (length(upper)!=dim_int){
    ## CHECK - not completly correct I think
    # upper <- c(upper[length(upper)],upper)
    #upper <- rep(unlist(m_times[which(substr(names(m_times),2,2)%in%splitted_f_type)[dim_int+1]]),dim_int)
    #}
    #if (length(id_na)>0 & length(id_na)<dim_int) upper[id_na] <- upper[id_na+1]
    
    # Test output:
    if(any(is.na(lower)) | any(is.na(upper))) stop("output looks wrong")
    
    limit <- list(lower=unname(lower),upper=unname(upper),tmax=tmax)
  }
  
  
  return(limit)
}

#' Calculate the negative log-likelihood
#'
#' Computes the negative log-likelihood for a set of patients, given by two nested lists "integrand"
#' and "limits".
#'
#' @param param The parameter values in which the log-likelihood should be evaluated. The dimension
#' and ordering is given by the user-specified densities.
#' @param integrands A nested list with the integrand(s) for a set of patients. The length of the 
#' outer list corresponds to the number of patients, the lengths of the inner lists correspond to the
#' the number of different paths the patient may have travelled.
#' @param limits A nested list with the limits of integration for a set of patients. The length of the 
#' outer list corresponds to the number of patients, the lengths of the inner lists correspond to the
#' the number of different paths the patient may have travelled.
#' @param X A matrix with the covariates. Number of rows equal to the number of patients and one column
#' for each covariate. The covariate specification is given by the user-specified densities.
#' @param cmethod The integration method of choice for the cubintegrate() function. Only for integrals 
#' of higher dimension than 2. Defaults to "hcubature".
#' @param mc_cores The number of cores to use (for parallelisation). The function uses the mclapply()
#' function from the parallel package.
#' @return The value of the negative log-likelihood corresponding to the dataset and parameters provided.
#' @export
mloglikelihood <-  function(param,integrands,limits, X = NULL,cmethod = "hcubature",mc_cores = 2){
  # Test that limits and integrand have same length
  
  final_integral = sum(unlist(parallel::mclapply(1:length(integrands), function(i){
    mm <- length(limits[[i]])
    lli <- rep(NA,mm)
    for (j in 1:mm){
      lower <- unlist(limits[[i]][[j]]$lower)
      upper <- unlist(limits[[i]][[j]]$upper)
      tmax <- unlist(limits[[i]][[j]]$tmax)
      
      
      if(length(lower) == 0){
        lli[j] = integrands[[i]][[j]](times=1,tt = tmax[1], tt2=tmax[2],param = param, x = X[i,]) # the value of times does not matter
      }else if(length(lower)>0){
        if(length(lower)<=2 ){
          # If integrate fails, use cubintegrate:
            lli[j] = tryCatch({
              repintegrate(integrands[[i]][[j]],tt=tmax[1],tt2=tmax[2],lower=lower,upper = upper, param = param, 
                                    x = X[i,])
            },error=function(cond){
              integrand2 <- change_integrand(integrands[[i]][[j]])
              if (length(unique(lower)) != length(lower)){
                llij = cubature::cubintegrate(integrand2, lower = lower,upper = upper, method = "divonne", maxEval = 500,
                                      tt = tmax[1], tt2=tmax[2],param = param, x = X[i,])$integral
              }else if (length(unique(lower)) == length(lower)){
                llij = cubature::cubintegrate(integrand2, lower = lower,upper = upper,maxEval = 500,
                                      method = cmethod, tt = tmax[1], tt2=tmax[2],param = param, x = X[i,])$integral
              }
              return(llij)
            })
            
        }else if (length(lower)>2){
          if (length(unique(lower)) != length(lower)){
            lli[j] = cubature::cubintegrate(integrands[[i]][[j]], lower = lower,upper = upper, method = "divonne", maxEval = 500,
                                  tt = tmax[1], tt2=tmax[2],param = param, x = X[i,])$integral
          }else if (length(unique(lower)) == length(lower)){
            lli[j] = cubature::cubintegrate(integrands[[i]][[j]], lower = lower,upper = upper,maxEval = 500,
                                  method = cmethod, tt = tmax[1], tt2=tmax[2],param = param, x = X[i,])$integral
          }
        }
      }
      
    }
    -log(sum(lli))
  }, mc.cores = mc_cores)))
  return(final_integral)
}
