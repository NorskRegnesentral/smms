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
  all_edges = igraph::get.edgelist(graph)
  id_initial = which(!(all_edges[,1] %in% all_edges[,2]))
  initial_states = unique(all_edges[id_initial,1])

  id_absorbing = which(!(all_edges[,2] %in% all_edges[,1]))
  absorbing_states = unique(all_edges[id_absorbing,2])

  k = length(igraph::as_ids(igraph::V(graph))) #number of states

  state_num = data.frame(state=igraph::as_ids(igraph::V(graph)),order=NA,type=NA)
  state_num$order[which(state_num$state %in% initial_states)] <- 0:(length(initial_states)-1)
  state_num$type[which(state_num$state %in% initial_states)] <- "init"

  # For the transient states (not initial or absorbing):
  # the ones on the same path should be ordered according to distance from the initial state
  # intuition: the ordering of the transient states must be such that all edge pairs go from a lower number to a higher number
  transient_states <-  state_num$state[-which(state_num$state %in% initial_states | state_num$state %in% absorbing_states)]
  transient_edges <- all_edges[which(all_edges[,1] %in% transient_states & all_edges[,2] %in% transient_states),,drop=F]
  m_dists <- rep(NA,length(transient_states))
  names(m_dists) <- transient_states
  low <- unique(transient_edges[which(!(transient_edges[,1] %in% transient_edges[,2])),1])
  high <- unique(transient_edges[which(!(transient_edges[,2] %in% transient_edges[,1])),2])
  m_dists[low] <- 1
  m_dists[high] <- length(transient_states)
  
  mid <-  unique(transient_edges[-which(transient_edges %in% low | transient_edges %in% high)])
  if (length(mid)==1) m_dists[mid] <- 2
  ii <- 2
  while(length(mid)>1){
    mid_edges <- transient_edges[which(transient_edges[,1] %in% mid & transient_edges[,2] %in% mid),,drop=F]
    low <- unique(mid_edges[which(!(mid_edges[,1] %in% mid_edges[,2])),1])
    high <- unique(mid_edges[which(!(mid_edges[,2] %in% mid_edges[,1])),2])
    m_dists[low] <- ii
    m_dists[high] <- length(transient_states)-ii+1
    ii <- ii+1
    
    mid <-  unique(mid_edges[-which(mid_edges %in% low | mid_edges %in% high)])
    if (length(mid)==1) m_dists[mid] <- ii
  }
  state_num$order[which(state_num$state %in% transient_states)] <- rank(m_dists,ties.method="first")+max(state_num$order,na.rm=T)
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
relevant_timepoints = function(data, graph){
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
  dists <- igraph::distances(graph,mode="in")

  # If the patient does not start in an initial state - add the appropriate initial state at time=0 (if the initial state is not
  # uniquely defined one just has to add one of the initial states)
  for (i in 1:nn){
    ddi = ddr[which(ddr$patient==inds[i]),]
    if (sum(ddi$state %in% init)==0){
      min_state <- state_ord$state[which(state_ord$order==min(ddi$state))]
      dists_init <- dists[,colnames(dists)%in%state_ord$state[which(state_ord$type=="init")],drop=F]
      id_init_min <- colnames(dists_init)[which.min(dists_init[min_state,])]
      ddi_1 <- c(inds[i],0,state_ord$order[which(state_ord$state==id_init_min)])
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
  if (length(idd)>0) ddr = ddr[-idd,]
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
    paths = igraph::all_simple_paths(graph, state_ord$state[which(state_ord$order==subsets[p,1])],
                             state_ord$state[which(state_ord$order==subsets[p,2])])
    ## If only observed in initial state
    if(length(paths) == 0){
      form_types = c(form_types, as.character(subsets[p, 1]))
    } else{ ## For all the other possible roads to travel
      for(i in 1:length(paths)){
        st = sort(state_ord$order[which(state_ord$state%in%igraph::as_ids(paths[[i]]))])
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
      ot = sapply(2:(length(st)-1), function(r) utils::combn(st[1:length(st)],r),simplify=F)
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

#' Arrange data set
#'
#' Discard unrelevant time-points and arrange data into one row per patient format. The column names for the relevant time-points
#' are on the format t_{im}/t_{iM}, where i is the state numbering according to state_ordering() and m/M indicate the minimum and maximum 
#' time-point in each state. For initial states only the maximum time-point is relevant, and for absorbing states only the minimum.
#'
#' @param data A data frame with 3 columns: patient (numbering or names for each patient), time (the time when a patient was observed),
#' state (the state which the patient occupies at the observation time).
#' @param graph A directed, acyclic graph in the igraph format (igraph package).
#' @return A data frame with as many rows as there are patients and with one column per relevant time-point pluss one column indicating 
#' the observation type of each patient.
#'
#'
arrange_data = function(data, graph){
  state_ord = state_ordering(graph)
  init <- sort(state_ord$order[which(state_ord$type=="init")])
  trans <- sort(state_ord$order[which(state_ord$type=="trans")])
  abs <- sort(state_ord$order[which(state_ord$type=="abs")])

  ## Back to the timepoints - making an empty matrix with correct time points
  ddr = relevant_timepoints(data, graph)
  inds = unique(ddr$patient)
  nn = length(inds)
  timepoints = as.data.frame(matrix(NA,nn,length(init) + length(abs) + 2*length(trans) +1))
  names_middle = paste(rep("t", length(init) + length(abs) + 2*length(trans)),
                       c(init, rep(trans,each = 2), abs),
                       c(rep("M",length(init)),rep(c("m","M"),length(trans)),rep("m",length(abs))), sep="")
  colnames(timepoints) = c(names_middle, "obs_type")

  ## Filling out the timepoints matrix
  for (i in 1:nn){
    ddi = ddr[which(ddr$patient==inds[i]),]
    tti = ddi$time[pmatch(c(init, rep(trans, each = 2),abs),ddi$state)]
    names(tti) = as.character(c(init, rep(trans, each = 2),abs),ddi$state)
    
    # For transient states where only a single time-point has been observed, the second 
    # time-point is set equal to the first t_{im}=t_{iM}
    id_dupl <- which(names(tti) %in% trans & duplicated(names(tti)) & is.na(tti))
    tti[id_dupl] <- tti[id_dupl-1]
    
    timepoints[i,1:(ncol(timepoints)-1)] = tti
    ## Which observed type the individual is
    timepoints[i, ncol(timepoints)] = paste(unique(names(tti[!(is.na(tti))])), collapse ="")
  }
  return(timepoints)
}

#' Make edge matrices
#'
#'Make three matrices with the formula types as rows, and edges in the graph as columns. For each formula type, 
#'non-zero numbers indicate that edges have been "travelled", or "passed by", or if they are "possible next", 
#'meaning that they may be travelled in the next step. The non-zero numbers also indicate in what order the 
#'edges have been encountered. Zeros indicate edges that were not travelled/passed by/ possible next.
#'
#' @param graph A directed, acyclic graph in the igraph format (igraph package).
#' @return A list of three matrices: travelled, passedBy, and possible next.
#'
#'
edge_matrices = function(graph){
  all_edges = igraph::get.edgelist(graph)
  # Update with state ordering as node names:
  state_ord = state_ordering(graph)
  all_edges[,1] <- state_ord$order[match(all_edges[,1],state_ord$state)]
  all_edges[,2] <- state_ord$order[match(all_edges[,2],state_ord$state)]
  edge_names = apply(all_edges,1,paste,collapse="")
  
  formula_types = construct_formula_types(graph)

  ## Make a matrix for all the edges traveled
  matrix_travelled = matrix(data = 0, nrow = length(formula_types), ncol = length(edge_names))
  rownames(matrix_travelled) = formula_types
  colnames(matrix_travelled) = edge_names
  for(i in 1:length(formula_types)){
    formula_str_pairs = substring(formula_types[i], first = 1:(nchar(formula_types[i]) - 1), last = 2:nchar(formula_types[i]))
    matches = which(edge_names %in% formula_str_pairs)
    matches_order = match(edge_names[matches],formula_str_pairs)
    if (length(matches)==0) next
    matrix_travelled[i,matches] <- matches_order
  }

  ## Make a matrix for all the edges which potentially can be traveled (in the next step)
  matrix_possible_next = matrix(data = 0, nrow = length(formula_types), ncol = length(edge_names))
  rownames(matrix_possible_next) = formula_types
  colnames(matrix_possible_next) = edge_names
  for(i in 1:length(formula_types)){
    formula_types_split = unlist(strsplit(formula_types[i], ""))
    id_next <- which(all_edges[,1]==formula_types_split[length(formula_types_split)])
    if (length(id_next)==0) next
    matrix_possible_next[i,id_next] <- length(formula_types_split)
  }

  ## Make a matrix for all the edges which were not traveled (i.e. where passed by)
  matrix_passed = matrix(data = 0, nrow = length(formula_types), ncol = length(edge_names))
  rownames(matrix_passed) = formula_types
  colnames(matrix_passed) = edge_names
  for(i in 1:length(formula_types)){
    formula_str_pairs = substring(formula_types[i], first = 1:(nchar(formula_types[i]) - 1), last = 2:nchar(formula_types[i]))
    for (j in 1:length(formula_str_pairs)){
      pair_split = unlist(strsplit(formula_str_pairs[j], ""))
      match_passed = which(pair_split[1]==all_edges[,1] & pair_split[2]!=all_edges[,2])
      if (length(match_passed)==0) next
      matrix_passed[i,match_passed] = j
    }
  }
  list_all_edges = list("traveled" = matrix_travelled , "passedBy" = matrix_passed, "possible" = matrix_possible_next)
  return(list_all_edges)
}
