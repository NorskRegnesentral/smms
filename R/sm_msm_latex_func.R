### Multi-state functions part 4 - Writing out likelihood as latex code

#' Write likelihood contribution for one observation type to latex type format
#'
#' @param obs_type The observation type (as a string).
#' @param graph A directed, acyclic graph in the igraph format (igraph package).
#' @param abs_exact  A boolean indicating whether the time of entrance into absorbing states is observed
#' exactly (TRUE) or not (FALSE). Default value is TRUE.
#' @return A text string with the likelihood contribution formula written in a latex type format.
#' 
write_type <- function(obs_type,graph,abs_exact=TRUE){
  
  formula_obs_types = all_types(graph)
  edge_mats <- edge_matrices(graph)
  state_ord <- state_ordering(graph)
  names_surv_dens = names_of_survival_density(graph)
  
  absorbing_states <- sort(state_ord$order[which(state_ord$type=="abs")])
  
  init <- sort(state_ord$order[which(state_ord$type=="init")])
  trans <- sort(state_ord$order[which(state_ord$type=="trans")])
  abs <- sort(state_ord$order[which(state_ord$type=="abs")])
  
  f_types = names(which(formula_obs_types[, obs_type] == 1))
  lik_parts <- rep(NA,length(f_types))
  
  for (k in 1:length(f_types)){
    form_type=f_types[k]
    
    integrand = type_to_integrand(form_type,edge_mats, names_surv_dens,abs_exact=abs_exact)
    
    # Find integral dimension and tmax2 
    splitted_f_type = unlist(strsplit(form_type, ""))
    choice <- (length(which(edge_mats$passedBy[form_type,]==max(edge_mats$traveled[form_type,])))>0)  #TRUE: if absorbing state is reached from a state where there was an option to go a different way
    intoAbs <- splitted_f_type[length(splitted_f_type)] %in% absorbing_states
    
    otype_splt <- unlist(strsplit(obs_type, ""))
    jump_states <- splitted_f_type[which(!(splitted_f_type %in% otype_splt))]
    
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
          if(splitted_f_type[length(splitted_f_type)-1] %in% jump_states){
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
      int <- gsub("tt-ww","tt - ss - uu - rr - vv - ww",int)
      int <- gsub("tt2-uu","tt2 - ss - uu",int)
      int <- gsub("tt2-rr","tt2 - ss - uu - rr",int)
      int <- gsub("tt2-vv","tt2 - ss - uu - rr - vv",int)
      int <- gsub("tt2-ww","tt2 - ss - uu - rr - vv - ww",int)
    }else{
      int <- integrand
    }
    
    # Remove stuff and change:
    int <- sub("function\\(.+,x,.+\\)\\{","",int)
    int <- sub("}$","",int)
    int <- gsub("*","",int,fixed=T)
    int <- gsub("tt2","t2",int)
    int <- gsub("uu","u",int)
    int <- gsub("ss","s",int)
    int <- gsub("rr","r",int)
    int <- gsub("vv","v",int)
    int <- gsub("ww","w",int)
    int <- gsub("param,x,","",int)
    
    m <- gregexpr("[fS]_[[:digit:]]{2}",int)
    call <- regmatches(int,m)
    for (i in 1:length(call[[1]])){
      callspl <- unlist(strsplit(call[[1]][i],""))
      int <- gsub(call[[1]][i],paste(callspl[1],"_{",paste(callspl[3:4],collapse=""),"}",sep=""),int)
    }
    
    # Make timepoints (belonging to the obs_type)
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
    
    # Make integral limits
    if (dim_int>0){
      lower <- rep(NA,dim_int)
      upper <- rep(NA,dim_int)
      lower[1] <- timepoints[1]
      upper[1] <- timepoints[2]
      
      if (dim_int>1){
        add_ons <- c("-s","-s-u","-s-u-r","-s-u-r-v","-s-u-r-v-w")
        for (j in 2:dim_int){
          current_state <- splitted_f_type[j]
          next_state <- splitted_f_type[j+1]
          if (next_state %in% jump_states & (current_state %in% jump_states)){
            lower[j] <- 0
            id_up <- which(1:length(timepoints)%in%seq(2,20,2) & as.numeric(substr(timepoints,4,4))>as.numeric(next_state))
            upper[j] <- paste(timepoints[min(id_up)],add_ons[j-1],sep="")
          }else if (next_state %in% jump_states & !(current_state %in% jump_states)){
            id_low <- which(1:length(timepoints)%in%seq(1,19,2) & as.numeric(substr(timepoints,4,4))==as.numeric(current_state))
            id_up <- which(1:length(timepoints)%in%seq(2,20,2) & as.numeric(substr(timepoints,4,4))>as.numeric(next_state))
            lower[j] <- paste(timepoints[id_low],add_ons[j-1],sep="")
            upper[j] <- paste(timepoints[min(id_up)],add_ons[j-1],sep="")
          }else if (!(next_state %in% jump_states) & (current_state %in% jump_states)){
            lower[j] <- 0
            id_up <- which(1:length(timepoints)%in%seq(2,20,2) & as.numeric(substr(timepoints,4,4))==as.numeric(next_state))
            upper[j] <- paste(timepoints[min(id_up)],add_ons[j-1],sep="")
          }else{
            id_low <- which(1:length(timepoints)%in%seq(1,19,2) & as.numeric(substr(timepoints,4,4))==as.numeric(current_state))
            id_up <- which(1:length(timepoints)%in%seq(2,20,2) & as.numeric(substr(timepoints,4,4))==as.numeric(next_state))
            lower[j] <- paste(timepoints[id_low],add_ons[j-1],sep="")
            upper[j] <- paste(timepoints[id_up],add_ons[j-1],sep="")
          }
        }
      }
      
      ints <- paste("\\int_{",lower,"}^{",upper,"}",sep="", collapse=" ")
      
      endings <- c("\\,ds","\\,du \\,ds","\\,dr \\,du \\,ds","\\,dv \\,dr \\,du \\,ds","\\,dw \\,dv \\,dr \\,du \\,ds")
      ends <- endings[dim_int]
    }else{
      ints = " "
      ends = " "
    }
    
    # Insert timepoints in which to evaluate functions:
    int <- gsub("tt",paste(timepoints[length(timepoints)],sep=""),int)
    if (tmax2==-1){
      int <- sub("S_\\{[[:digit:]]{2}\\} \\( t2.+?\\)","1",int)
    }else{
      int <- sub("t2",timepoints[length(timepoints)-1],int)
    }
    
    # Paste everything together:
    lik_parts[k] <- paste(ints,int,ends,collapse=" ")
  }
  
  lik <- paste("\\sum_{k=0}^{n_{",obs_type,"}} \\log \\{",paste(lik_parts,collapse=" + "),"\\}",sep="")
  return(lik)
}

#' Write log-likelihood for a graph in latex format.
#'
#' @param graph A directed, acyclic graph in the igraph format (igraph package).
#' @param abs_exact  A boolean indicating whether the time of entrance into absorbing states is observed
#' exactly (TRUE) or not (FALSE). Default value is TRUE.
#' @return Writes a txt file into the working directory within which one will find the entire log-likelihood
#' formula belonging to the graph, in latex format.
#' @export
write_loglikelihood <- function(graph,abs_exact=TRUE){
  o_types <- construct_obs_types(graph)
  all_parts <- rep(NA,length(o_types))
  
  for (i in 1:length(all_parts)){
    all_parts[i] <- write_type(obs_type=o_types[i],graph=graph,abs_exact=abs_exact)
  }
  eq <- paste("\\begin{align*} \n \\ell_n (\\theta)= ",paste(all_parts,collapse=" + "),". \n \\end{align*}")
  writeLines(eq, "likelihood_latex.txt")
}
