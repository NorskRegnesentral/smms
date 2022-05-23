### Multi-state functions part 1 - preprocessing

# Refer to functions with `igraph::fun()`

#' Produce an ordering of the states from a multi-state graph
#'
#' @param graph A directed, acyclic graph in the igraph format (igraph package).
#' @return A data frame with 3 columns giving a mapping from the state names (first column) to a partial ordering of the states
#' (second column, from 0 to k, and where a higher number indicates that a state is further removed from the initial state).
#' The last column indicates whether each state is initial, absorbing or transient.
state_ordering = function(graph){
  ## Calculating the initial state
  all_edges = get.edgelist(graph)
  id_initial = which(!(all_edges[,1] %in% all_edges[,2]))
  initial_states = unique(all_edges[id_initial,1])

  id_absorbing = which(!(all_edges[,2] %in% all_edges[,1]))
  absorbing_states = unique(all_edges[id_absorbing,2])

  k = length(as_ids(V(graph))) #number of states

  state_num = data.frame(state=as_ids(V(graph)),order=NA,type=NA)
  state_num$order[which(state_num$state %in% initial_states)] <- 0:(length(initial_states)-1)
  state_num$type[which(state_num$state %in% initial_states)] <- "init"

  # For the transient states (not initial or absorbing):
  # the ones on the same path should be ordered according to distance from the initial state
  dists <- distances(graph,mode="in")
  dists[which(is.infinite(dists))] <- NA
  transient_states <-  state_num$state[-which(state_num$state %in% initial_states | state_num$state %in% absorbing_states)]
  m_dists <- rep(NA,length(transient_states))
  for (i in 1:length(transient_states)){
    m_dists[i] <-  max(dists[transient_states[i],],na.rm=T) #min(dists[transient_states[i],initial_states],na.rm=T)
  }
  state_num$order[which(state_num$state %in% transient_states)] <- order(m_dists)+max(state_num$order,na.rm=T)
  state_num$type[which(state_num$state %in% transient_states)] <- "trans"

  state_num$order[which(state_num$state %in% absorbing_states)] <- (k-1-length(absorbing_states)+1):(k-1)
  state_num$type[which(state_num$state %in% absorbing_states)] <- "abs"
  return(state_num)
}

#' Filter away redundant observations in a multi-state dataset
#'
#' Only the last timepoint in initial states are relevant. Only the first timepoint in absorbing states. For the
#' transient states the first and last timepoints are relevant.
#'
#' @param data A data frame with 3 columns: patient (numbering or names for each patient), time (the time when a patient was observed),
#' state (the state which the patient occupies at the observation time).
#' @param graph A directed, acyclic graph in the igraph format (igraph package).
#' @return A data frame with the same 3 columns as 'data', but with fewer rows (potentially).
#'
relevant_timepoints = function(formula,subject,data, graph){
  inds = unique(data$patient)
  nn = length(inds)
  ddr = data
  # Add state ordering instead of state names
  state_ord = state_ordering(graph)
  ddr$state = with(state_ord, order[match(ddr$state, state)])

  # Identify initial and absorbing states
  init <- state_ord$order[which(state_ord$type=="init")]
  trans <- state_ord$order[which(state_ord$type=="trans")]
  abs <- state_ord$order[which(state_ord$type=="abs")]
  dists <- distances(graph,mode="in")

  # If the patient does not start in an initial state - add the appropriate initial state at time=0 (if the initial state is not
  # uniquely defined one just has to add one of the initial states)
  for (i in 1:nn){
    ddi = ddr[which(ddr$patient==inds[i]),]
    if (sum(ddi$state %in% init)==0){
      min_state <- state_ord$state[which(state_ord$order==min(ddi$state))]
      id_init <- which.min(dists[min_state,-which(colnames(dists)==min_state)])[1]
      ddi_1 <- c(ddi$patient[i],0,state_ord$order[which(state_ord$state==names(id_init))])
      ddi <- rbind(ddi_1,ddi)
      ddr <- rbind(ddi,ddr[-which(ddr$patient==inds[i]),])
    }else{
      next
    }
  }

  # Identify redundant timepoints
  idd = NULL
  for (i in 1:nn){
    ddi = ddr[which(ddr$patient==inds[i]),]

    rle_i <- rle(ddi$state)
    id_unrelevant = NULL
    for(j in 1:length(rle_i$values)){
      val <- rle_i$values[j]
      num <- rle_i$lengths[j]
      if(val%in%init & num >= 2){ #initial states
        id_unrelevant = which(ddi$state==val)[1:(num-1)]  # assumes that the timepoints are sorted?
      } else if(val%in%trans & num >= 3){ #transient states
        id_unrelevant = c(id_unrelevant,which(ddi$state==val)[-c(1,num)])
      } else if(val%in%abs & num >= 2){ #absorbing states
        id_unrelevant = c(id_unrelevant,which(ddi$state==val)[2:num])
      }
    }
    idd = c(idd,which(ddr$patient==inds[i])[id_unrelevant])
  }
  ddr = ddr[-idd,]
  return(ddr)
}

#' Find all formula types
#'
#' The formula types are the combination of states a patient can be last observed in and possible paths from
#' the initial states to these “last observed” states. The list of states is the path actually traveled by the patient,
#' but some of the transient state might actually be unobserved.
#'
#' @param graph A directed, acyclic graph in the igraph format (igraph package).
#' @return A vector with string elements indicating the states that were visited/traveled.
#'
construct_formula_types = function(graph){
  state_ord = state_ordering(graph)
  init <- state_ord$order[which(state_ord$type=="init")]
  trans <- state_ord$order[which(state_ord$type=="trans")]
  abs <- state_ord$order[which(state_ord$type=="abs")]

  k <- dim(state_ord)[1]

  ## Find all combination of initial states and other states
  subsets = matrix(nrow = length(init)*k, ncol = 2)
  subsets[,1] <- rep(init,each=k)
  subsets[,2] <- rep(c(init,trans,abs),times=length(init))

  ## Determine which subsets containt none, one or multiple paths
  form_types = c()
  for(p in 1:nrow(subsets)){
    paths = all_simple_paths(graph, state_ord$state[which(state_ord$order==subsets[p,1])],
                             state_ord$state[which(state_ord$order==subsets[p,2])])
    ## If only observed in initial state
    if(length(paths) == 0){
      form_types = c(form_types, as.character(subsets[p, 1]))
    } else{ ## For all the other possible roads to travel
      for(i in 1:length(paths)){
        st = sort(state_ord$order[which(state_ord$state%in%as_ids(paths[[i]]))])
        data_frame_subset_type_as_char = paste(st, collapse = "")
        form_types = c(form_types, data_frame_subset_type_as_char)
      }
    }
  }
  form_types = unique(form_types)
  return(form_types)
}


#' Find all observation types
#'
#' The observation types are the combination of states in which a patient can
#' actually observed. It does not necessarily indicate which path was traveled:
#' some observation types may correspond to multiple paths.
#'
#' @param graph A directed, acyclic graph in the igraph format (igraph package).
#' @return A vector with string elements indicating the states in which the patient is observed.
#'
construct_obs_types = function(graph){
  form_types = construct_formula_types(graph)
  obs_types = c()
  for (i in 1:length(form_types)){
    st = strsplit(form_types[i],"")[[1]]
    obs_types = c(obs_types,form_types[i])
    if (length(st)>2){
      ot = sapply(2:(length(st)-1), function(r) combn(st[1:length(st)],r),simplify=F)
      ot = lapply(ot,function(m) m[,which(m[1,]==st[1])])
      obs_types = c(obs_types,unlist(lapply(ot,function(m) apply(m,2,paste,collapse=""))))
    }
  }
  obs_types = unique(obs_types)
  return(obs_types)
}

#' Find links between formula and observation types
#'
#' Some observation types may belong to several different formula types, and some formula types
#' may give rise to several observation types.
#'
#' @param graph A directed, acyclic graph in the igraph format (igraph package).
#' @return A matrix with formula types as rows and observation types as columns: "1" indicates
#' a link, "0" indicates no link.
#'
all_types = function(graph){
  formula_types = construct_formula_types(graph)
  observation_types = construct_obs_types(graph)
  matrix_all_types = matrix(data = 0, nrow = length(formula_types), ncol = length(observation_types))
  rownames(matrix_all_types) = formula_types
  colnames(matrix_all_types) = observation_types
  for(i in 1:length(formula_types)){
    formula_types_split = unlist(strsplit(formula_types[i], ""))
    for(j in 1:length(observation_types)){
      observation_types_split = unlist(strsplit(observation_types[j], ""))
      if(formula_types_split[1] == observation_types_split[1] &
         formula_types_split[length(formula_types_split)] == observation_types_split[length(observation_types_split)] &
         all(observation_types_split %in% formula_types_split)){
        matrix_all_types[i, j] = 1
      }
    }
  }
  return(matrix_all_types)
}


## 2.6 Arrange data
## Input: data (the dataset to consider), graph (the graph)
arrange_data = function(data, graph, abs_int_cens = NULL){
  ## Finding the initial, transient and absorbing states - in the "new" format with 0, 1, ...
  all_edges = get.edgelist(graph)
  all_initial = which(!(all_edges[,1] %in% all_edges[,2]))
  initial_state = unique(sapply(all_initial, function(p) all_edges[p,1]))
  all_absorbing = which(!(all_edges[,2] %in% all_edges[,1]))
  absorbing_state = unique(sapply(all_absorbing, function(p) all_edges[p,2]))

  all_states_old = state_ordering(graph)[,1]
  all_states_new = state_ordering(graph)[,2]

  row_initial = sapply(1:length(initial_state), function(p) which(all_states_old == initial_state[p], arr.ind = TRUE))
  row_absorbing = sapply(1:length(absorbing_state), function(p) which(all_states_old == absorbing_state[p], arr.ind = TRUE))


  transient_state_new = c()
  for(p in 1:length(all_states_old)){
    if(!(all_states_old[p] %in% absorbing_state) & !(all_states_old[p] %in% initial_state)){
      transient_state_new  = c(transient_state_new, all_states_new[p])
    }
  }
  initial_state_new = sapply(row_initial, function(p) state_ordering(graph)[p,2])
  absorbing_state_new = sapply(row_absorbing, function(p) state_ordering(graph)[p,2])


  if(!is.null(abs_int_cens)){
    row_absorbing = sapply(1:length(absorbing_state), function(p) which(all_states_old == absorbing_state[p], arr.ind = TRUE))
    absorbing_state_int_cens = sapply(row_absorbing, function(p) state_ordering(graph)[p,2])
    absorbing_state_r_cens = which(!(absorbing_state_new %in% absorbing_state_int_cens))
  }


  ## Back to the timepoints - making an empty matrix with correct time points
  ddr = relevant_timepoints(data, graph)
  inds = unique(ddr$patient)
  nn = length(inds)
  if(is.null(abs_int_cens)){
    timepoints = matrix(NA,nn,length(initial_state) + length(absorbing_state) + 2*length(transient_state_new) +1)
    names_middle = paste(rep("t", length(initial_state) + 2*length(transient_state_new) + length(absorbing_state)),
                         c(initial_state_new, rep(transient_state_new,each = 2), absorbing_state_new),
                         c("M",rep(c("m","M"),length(transient_state_new)),rep("m",length(absorbing_state))), sep="")
    ### C has fixed above

  } else{
    timepoints = matrix(NA,nn,length(initial_state) + length(absorbing_state)- length(abs_int_cens) + 2*length(transient_state_new)+2*length(abs_int_cens) +1)
    names_middle = paste(rep("t", length(initial_state) + 2*length(transient_state_new) + 2*length(abs_int_cens) + length(absorbing_state_new)-length(abs_int_cens)),
                         c(initial_state_new, rep(transient_state_new,each = 2), rep(abs_int_cens, each = 2), absorbing_state_r_cens),
                         rep(c("M", "m"),
                             (length(initial_state) + 2*length(transient_state_new) + 2*length(absorbing_state_int_cens) + length(absorbing_state_r_cens))/2), sep="")
  }
  colnames(timepoints) = c(names_middle, "obs_type")

  if(is.null(abs_int_cens)){
    ## Filling out the timepoints matrix
    otypes = rep(NA,nn) # a vector indicating the observation type of each patient
    for (i in 1:nn){
      ddi = ddr[which(ddr$patient==inds[i]),]
      tti = ddi$time[pmatch(c(initial_state_new, rep(transient_state_new, each = 2),
                               absorbing_state_new),ddi$state)]
      names(tti) = as.character(c(initial_state_new, rep(transient_state_new, each = 2),
                                  absorbing_state_new),ddi$state)

      for(j in 2:length(tti)){
        ## If the initial state(s) column number is odd/even, then all the "M"'s are also odd/even.
        ## It is the "M"'s which potentially are NA. Do not want to fill them if both are NA
        if((length(initial_state) %% 2) == (j %% 2) & j > length(initial_state) &
           j <= (length(initial_state) + 2*length(transient_state_new))){
          if (is.na(tti[j]) & !is.na(tti[j-1])){
            tti[j] = tti[j-1]
          }
        }
      }
      timepoints[i,1:(ncol(timepoints)-1)] = tti
      ## Which observed type the individual is
      timepoints[i, ncol(timepoints)] = paste(unique(names(tti[!(is.na(tti))])), collapse ="")
    }
  } else{
    ## Filling out the timepoints matrix
    otypes = rep(NA,nn) # a vector indicating the observation type of each patient
    for (i in 1:nn){
      ddi = ddr[which(ddr$patient==inds[i]),]
      tti = ddi$time[pmatch(c(initial_state_new, rep(transient_state_new, each = 2),
                               rep(absorbing_state_int_cens, each = 2), absorbing_state_r_cens),ddi$state)]
      names(tti) = as.character(c(initial_state_new, rep(transient_state_new, each = 2),
                                  rep(absorbing_state_int_cens, each = 2), absorbing_state_r_cens),ddi$state)

      for(j in 2:length(tti)){
        ## If the initial state(s) column number is odd/even, then all the "M"'s are also odd/even.
        ## It is the "M"'s which potentially are NA. Do not want to fill them if both are NA
        if((length(initial_state) %% 2) == (j %% 2) & j > length(initial_state) &
           j <= (length(initial_state) + 2*length(transient_state_new) + 2*length(absorbing_state_int_cens))){
          if (is.na(tti[j]) & !is.na(tti[j-1])){
            tti[j] = tti[j-1]
          }
        }
      }
      timepoints[i,1:(ncol(timepoints)-1)] = tti
      ## Which observed type the individual is
      timepoints[i, ncol(timepoints)] = paste(unique(names(tti[!(is.na(tti))])), collapse ="")
    }
  }

  return(timepoints)
}


## 2.7 Make edge matrices
## Input: graph (the graph)
edge_matrices = function(graph){
  ## Make a duplicate to not change the original graph
  graph2 = graph
  ## The ordering of V(graph) can be different from state_ordering
  new_graph_name = c()
  for(i in 1:nrow(state_ordering(graph2))){
    r = match(V(graph2)$name[i], state_ordering(graph2)[,1])
    new_graph_name = c(new_graph_name, state_ordering(graph2)[r,2])
  }
  V(graph2)$name = new_graph_name
  data_frame_possible_travels = get.data.frame(graph2, what= "edges" )

  ## Make all of the possible travels into characters
  data_frame_travels = c()
  for(i in 1:nrow(data_frame_possible_travels)){
    data_frame_travels = c(data_frame_travels, paste(data_frame_possible_travels[i,], collapse = ""))
  }
  formula_types = construct_formula_types(graph)

  ## Make a matrix for all the edges traveled
  matrix_travelled = matrix(data = 0, nrow = length(formula_types), ncol = length(data_frame_travels))
  rownames(matrix_travelled) = formula_types
  colnames(matrix_travelled) = data_frame_travels
  for(i in 1:length(formula_types)){
    formula_types_split = unlist(strsplit(formula_types[i], ""))
    for(j in 1:nrow(data_frame_possible_travels)){
      for(p in 1:(length(formula_types_split)-1)){
        if(length(formula_types_split) > 1){
          if(formula_types_split[p] == data_frame_possible_travels[j, 1] &
             formula_types_split[p+1] == data_frame_possible_travels[j, 2]){
            matrix_travelled[i, j] = length(0:formula_types_split[p])
          }
        }
      }
    }
  }

  ## Make a matrix for all the edges which potentially can be traveled
  matrix_possible_to_travel = matrix(data = 0, nrow = length(formula_types), ncol = length(data_frame_travels))
  rownames(matrix_possible_to_travel) = formula_types
  colnames(matrix_possible_to_travel) = data_frame_travels
  for(i in 1:length(formula_types)){
    formula_types_split = unlist(strsplit(formula_types[i], ""))
    for(j in 1:nrow(data_frame_possible_travels)){
      if(formula_types_split[length(formula_types_split)] == data_frame_possible_travels[j, 1]){
        matrix_possible_to_travel[i, j] = length(formula_types_split)
      }
    }
  }

  ## Make a matrix for all the edges which were not traveled
  matrix_did_not_travel = matrix(data = 0, nrow = length(formula_types), ncol = length(data_frame_travels))
  rownames(matrix_did_not_travel) = formula_types
  colnames(matrix_did_not_travel) = data_frame_travels
  for(i in 1:length(formula_types)){
    formula_types_split = unlist(strsplit(formula_types[i], ""))
    for(j in 1:nrow(data_frame_possible_travels)){
      for(p in 1:(length(formula_types_split)-1)){
        if(length(formula_types_split) > 1)
          if(formula_types_split[p] == data_frame_possible_travels[j, 1] & formula_types_split[p+1] != data_frame_possible_travels[j, 2]){
            matrix_did_not_travel[i, j] = length(0:formula_types_split[p])
          }
      }
    }
  }
  list_all_edges = list("traveled" = matrix_travelled , "passedBy" = matrix_did_not_travel, "possible" = matrix_possible_to_travel)
  return(list_all_edges)
}
