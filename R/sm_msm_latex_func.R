### Multi-state functions part 4 - Writing out likelihood as latex code

# 
write_type <- function(obs_type,form_type,edge_mats,names_surv_dens,absorbing_states,abs_exact=TRUE){
  integrand = type_to_integrand(form_type,edge_mats, names_surv_dens,abs_exact=abs_exact)
  
  # Find integral dimension and tmax2 
  splitted_f_type = unlist(strsplit(form_type, ""))
  choice <- (length(which(edge_mats$passedBy[form_type,]==max(edge_mats$traveled[form_type,])))>0)  #TRUE: if absorbing state is reached from a state where there was an option to go a different way
  intoAbs <- splitted_f_type[length(splitted_f_type)] %in% absorbing_states
  
  if (intoAbs==TRUE){ #if patient is observed in absorbing state
    if (abs_exact==TRUE){
      dim_int <- length(splitted_f_type)-2 
      tmax2 = -1
    }else{
      if (choice==TRUE){
        dim_int <- length(splitted_f_type)-1
        tmax2 = -1
      }else{
        dim_int <- length(splitted_f_type)-2
        if(is.na(seclast[length(seclast)])){
          tmax2 = -1
        }else{
          tmax2 = paste("t_{",dim_int,"M}",sep="")
        }
      }
    }
  }else{ #if patient has not reached absorbing state
    dim_int <- length(splitted_f_type)-1 
    tmax2 = -1
  }
  
  # Change integrand to (sort of) repintegrate format if necessary
  if (grepl("times",integrand)){
    m <- regexpr("\\{ ss<-times\\[1\\]\n.+?[a-zA-Z]_",integrand)
    call <- regmatches(integrand,m)
    callspl <- unlist(strsplit(call,""))
    int <- sub("\\{ ss<-times\\[1\\]\n.+?[a-zA-Z]_",paste("{ ",paste(callspl[(length(callspl)-1):length(callspl)],collapse=""),sep=""),integrand)
    
    int <- gsub(", uu-ss",", uu",int)
    int <- gsub(", rr-uu",", rr",int)
    int <- gsub(", vv-rr",", vv",int)
    int <- gsub(", ww-vv",", ww",int)
    int <- gsub("tt-uu","tt - ss - uu",int)
    int <- gsub("tt-rr","tt - ss - uu - rr",int)
    int <- gsub("tt-vv","tt - ss - uu - rr - vv",int)
    int <- gsub("tt2-uu","tt2 - ss - uu",int)
    int <- gsub("tt2-rr","tt2 - ss - uu - rr",int)
    int <- gsub("tt2-vv","tt2 - ss - uu - rr - vv",int)
  }
  
  # Make integrals
  state_ord <- state_ordering(graph)
  otype_splt <- unlist(strsplit(obs_type, ""))
  timepoints = c()
  for (i in 1:length(otype_splt)){
    st <- as.numeric(otype_splt[i])
    if (state_ord$type[which(state_ord$order==st)]=="trans"){
      ti <- c(paste("t_{",st,"m}",sep=""),paste("t_{",st,"M}",sep=""))
    }else if(state_ord$type[which(state_ord$order==st)]=="init"){
      ti <- paste("t_{",st,"M}",sep="")
    }else if(state_ord$type[which(state_ord$order==st)]=="abs"){
      ti <- paste("t_{",st,"m}",sep="")
    }
    timepoints = c(timepoints,ti)
  }
  
  ints
  for (j in dim_int){
    
  }
  
}