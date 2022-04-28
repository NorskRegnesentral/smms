### Multi-state functions part 1 - preprocessing


## 2.1
## Input: gg (the graph)
state_ordering = function(gg){
  ## Calculating the initial state
  all_edges = get.edgelist(gg)
  all_initial = which(!(all_edges[,1] %in% all_edges[,2]))
  initial_state = unique(sapply(all_initial, function(r) all_edges[r,1]))
  
  state_num = matrix(nrow = 0, ncol = 2)
  ## Considering all the possible initial states
  for(i in 1:length(initial_state)){
    all_states = all_simple_paths(gg, initial_state[i])
    form_types = c()
    for(j in 1:length(all_states)){
      data_frame_subset_type_as_char = paste(as.data.frame(t(sapply(all_states[j], as_ids))), collapse = ".")
      form_types = c(form_types, data_frame_subset_type_as_char)
    }
    state_num = rbind(c(initial_state[i], as.character(i - 1)))
  }
  ## Split these into types
  str_split = matrix(ncol = max(sapply(1:length(form_types), function(j) length(unlist(strsplit(form_types[j], "[.]"))))), 
                     nrow = length(form_types))
  for(j in 1:length(form_types)){
    for(p in 1:length(unlist(strsplit(form_types[j], "[.]")))){
      str_split[j,p] = unlist(strsplit(form_types[j], "[.]"))[p]
    }
  }
  ## Arranging into a "predetermined" 
  for(j in 2:ncol(str_split)){
    for(p in 1:nrow(str_split)){
      if(j < ncol(str_split)){
        if(!(str_split[p, j] %in% state_num[,1]) & !(str_split[p, j] %in% str_split[, (j+1):ncol(str_split)]) & !is.na(str_split[p,j])){
          state_num = rbind(state_num, c(str_split[p,j], as.character(as.numeric(state_num[nrow(state_num),2]) + 1)))
        }
      } else if(j == ncol(str_split)){
        if(!(str_split[p, j] %in% state_num[,1]) & !is.na(str_split[p,j])){
          state_num = rbind(state_num, c(str_split[p,j], as.character(as.numeric(state_num[nrow(state_num),2]) + 1)))
        }
      }
    }
  }
  colnames(state_num) = c("prev", "new")
  return(state_num)
}


## 2.2 Keep only relevant timepoints
## Input: data_set (the dataset to consider), gg (the graph)
relevant_timepoints = function(data_set, gg){
  inds = unique(data_set$patient)
  nn = length(inds)
  ddr = data_set
  states_in_dataset = as.numeric(state_ordering(gg)[,1]) 
  all_states_ordered = as.numeric(state_ordering(gg)[,2])
  for(j in 1:length(states_in_dataset)){
    ddr$state = replace(ddr$state, which(ddr$state == states_in_dataset[j]), all_states_ordered[j])
  }
  idd = NULL
  for (i in 1:nn){
    ddi = data_set[which(data_set$patient==inds[i]),]
    states = ddi$state
    ddi$state = states
    rlei = matrix(0,2,length(all_states_ordered))
    rlei[1,] = all_states_ordered
    rlei[2, (rle(states)$values)] = rle(states)$lengths
    id_unrelevant = NULL
    for(j in 1:ncol(rlei)){
      if(j == 1 & rlei[2,j] >= 2){
        id_unrelevant = 1:(rlei[2,j]-1)
      } else if(j > 1 & j < ncol(rlei) & rlei[2,j] >= 3){
        id_unrelevant = c(id_unrelevant,(sum(rlei[2,1:(j-1)])+2):(sum(rlei[2,1:j])-1))
      } else if(j == ncol(rlei) & rlei[2,j] >= 2){
        id_unrelevant = c(id_unrelevant,(sum(rlei[2,1:(j-1)])+2):sum(rlei[2,1:j]))
      }
    }
    idd = c(idd,which(ddr$patient==inds[i])[id_unrelevant])
  }
  ddr = ddr[-idd,]
  return(ddr)
}


## 2.3 Construct formula types
## Input: gg (the graph)
construct_formula_types = function(gg){
  ## Tranforming the graph to the "new" states
  state_start = min(as.numeric(state_ordering(gg)[,2]))
  state_end = max(as.numeric(state_ordering(gg)[,2]))
  subsets = matrix(nrow = 0, ncol = 2)
  for(i in 1:length(state_start)){
    for(j in 1:length(state_end)){
      #subsets_expand_grid = expand.grid(state_start[i], c(state_start[i]:state_end[j]))
      subsets = rbind(subsets, expand.grid(state_start[i], c(state_start[i]:state_end[j])))
    }
  }
  
  ## The ordering of V(gg) can be different from state_ordering
  new_gg_name = c()
  for(i in 1:nrow(state_ordering(gg))){
    r = match(V(gg)$name[i], state_ordering(gg)[,1])
    new_gg_name = c(new_gg_name, state_ordering(gg)[r,2])
  }
  V(gg)$name = new_gg_name
  
  characters_subset_vector = NULL
  form_types = c()
  for(p in 1:nrow(subsets)){
    ## Which list of possible combination we consider
    characters_subset_matrix = matrix(NA, ncol = ncol(subsets), nrow = nrow(subsets))
    ## Make the states into characters in order to use all_simple_paths
    for(i in 1:ncol(subsets)){
      characters_subset_vector[i] = as.character(subsets[p, i])
    }
    subset_size = all_simple_paths(gg, characters_subset_vector[1], 
                                   characters_subset_vector[ncol(subsets)])
    ## If only observed in initial state
    if(length(subset_size) == 0){
      form_types = c(form_types, as.character(subsets[p, 1])) #characters_subset_matrix[p,1])
    } else{ ## For all the other possible roads to travel
      for(i in 1:length(subset_size)){
        data_frame_subset_type_as_char = paste(as.data.frame(t(sapply(subset_size[i], as_ids))), collapse = "")
        form_types = c(form_types, data_frame_subset_type_as_char)
      }
    }
  }
  return(form_types)
}


## 2.4 Find observation times
## Input: gg (the graph)
construct_obs_types = function(gg){
  ## Finding the initial state
  all_edges = get.edgelist(gg)
  all_initial = which(!(all_edges[,1] %in% all_edges[,2]))
  initial_state = unique(sapply(1:length(all_initial), function(r) all_edges[r,1]))
  
  num_states = length(as.numeric(state_ordering(gg)[,2]))
  num_init = length(initial_state)
  states = 0:(num_states-1) # assumes that the states are enumerated from 0 to k, with the initial state in the beginning
  num_obs_types = num_init + num_init*sum(choose((num_states-num_init),1:(num_states-num_init)))
  obs_types = rep(NA,num_obs_types)
  obs_types[1:num_init] = states[1:num_init]
  states_rest = states[-(1:num_init)]
  ii = which(is.na(obs_types))[1]
  for (i in 1:length(states_rest)){
    mat = matrix(0,(i+num_init),choose((num_states-num_init),i))
    mat[-1,] = combn(states_rest,i)    # here I've assumed a single initial state (FIX later)
    otypes = apply(mat,2,paste,collapse="")
    obs_types[ii:(ii+length(otypes)-1)] = otypes
    ii = ii+length(otypes)
  }
  return(obs_types)
}


## 2.5 Links between formula and observation types
## Input: gg (the graph) 
all_types = function(gg){
  formula_types = construct_formula_types(gg)
  observation_types = construct_obs_types(gg)
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
## Input: data_set (the dataset to consider), gg (the graph)
arrange_data = function(data_set, gg, abs_int_cens = NULL){
  ## Finding the initial, transient and absorbing states - in the "new" format with 0, 1, ...
  all_edges = get.edgelist(gg)
  all_initial = which(!(all_edges[,1] %in% all_edges[,2]))
  initial_state = unique(sapply(all_initial, function(p) all_edges[p,1]))
  all_absorbing = which(!(all_edges[,2] %in% all_edges[,1]))
  absorbing_state = unique(sapply(all_absorbing, function(p) all_edges[p,2]))
  
  all_states_old = state_ordering(gg)[,1]
  all_states_new = state_ordering(gg)[,2]
  
  row_initial = sapply(1:length(initial_state), function(p) which(all_states_old == initial_state[p], arr.ind = TRUE))
  row_absorbing = sapply(1:length(absorbing_state), function(p) which(all_states_old == absorbing_state[p], arr.ind = TRUE))
  
  
  transient_state_new = c()  
  for(p in 1:length(all_states_old)){
    if(!(all_states_old[p] %in% absorbing_state) & !(all_states_old[p] %in% initial_state)){
      transient_state_new  = c(transient_state_new, all_states_new[p])
    }
  }
  initial_state_new = sapply(row_initial, function(p) state_ordering(gg)[p,2])
  absorbing_state_new = sapply(row_absorbing, function(p) state_ordering(gg)[p,2])
  
  
  if(!is.null(abs_int_cens)){
    row_absorbing = sapply(1:length(absorbing_state), function(p) which(all_states_old == absorbing_state[p], arr.ind = TRUE))
    absorbing_state_int_cens = sapply(row_absorbing, function(p) state_ordering(gg)[p,2])
    absorbing_state_r_cens = which(!(absorbing_state_new %in% absorbing_state_int_cens))
  }
  
  
  ## Back to the timepoints - making an empty matrix with correct time points
  ddr = relevant_timepoints(data_set, gg)
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
## Input: gg (the graph)
edge_matrices = function(gg){
  ## Make a duplicate to not change the original graph
  gg2 = gg
  ## The ordering of V(gg) can be different from state_ordering
  new_gg_name = c()
  for(i in 1:nrow(state_ordering(gg2))){
    r = match(V(gg2)$name[i], state_ordering(gg2)[,1])
    new_gg_name = c(new_gg_name, state_ordering(gg2)[r,2])
  }
  V(gg2)$name = new_gg_name
  data_frame_possible_travels = get.data.frame(gg2, what= "edges" )
  
  ## Make all of the possible travels into characters
  data_frame_travels = c()
  for(i in 1:nrow(data_frame_possible_travels)){
    data_frame_travels = c(data_frame_travels, paste(data_frame_possible_travels[i,], collapse = ""))
  }
  formula_types = construct_formula_types(gg)
  
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
