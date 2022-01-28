## Code for multi-state models for interval-censored data
## Celine Cunen and Marthe Aastveit

rm(list = ls())
library(igraph)
library(parallel)
library(cubature)
## For test-data
library(msm)
## Initial graph for CAV-dataset
gg = graph_from_literal("1"--+"2"--+"3"--+"4", "1"--+"4", "2"--+"4")

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
state_ordering(gg)

## Testing for CAV-data
dd = cav
dd = dd[!is.na(dd$pdiag),]
dim(dd[dd$firstobs==1,]) # 614 patients
#first obs (at years=0) is always in state 1
id_wrong = unique(dd$PTNUM[which(dd$state!=dd$statemax)])  # observations where the patient appears to go back to a previous state
dd = dd[-which(dd$PTNUM %in% id_wrong),]
## Only relevant parts 
dd = dd[ ,-c(2, 4, 5, 6, 7, 9, 10)]


## 2.2 Keep only relevant timepoints
## Input: data_set (the dataset to consider), gg (the graph)
relevant_timepoints = function(data_set, gg){
  inds = unique(data_set$PTNUM)
  nn = length(inds)
  ddr = data_set
  states_in_dataset = as.numeric(state_ordering(gg)[,1]) 
  all_states_ordered = as.numeric(state_ordering(gg)[,2])
  for(j in 1:length(states_in_dataset)){
    ddr$state = replace(ddr$state, which(ddr$state == states_in_dataset[j]), all_states_ordered[j])
  }
  idd = NULL
  for (i in 1:nn){
    ddi = data_set[which(data_set$PTNUM==inds[i]),]
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
    idd = c(idd,which(dd$PTNUM==inds[i])[id_unrelevant])
  }
  ddr = ddr[-idd,]
  return(ddr)
}
relevant_timepoints(dd, gg)

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
construct_formula_types(gg)

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
construct_obs_types(gg)

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

all_types(gg)

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
  inds = unique(ddr$PTNUM)
  nn = length(inds)
  if(is.null(abs_int_cens)){
  timepoints = matrix(NA,nn,length(initial_state) + length(absorbing_state) + 2*length(transient_state_new) +1)
  names_middle = paste(rep("t", length(initial_state) + 2*length(transient_state_new) + length(absorbing_state)),
                       c(initial_state_new, rep(transient_state_new,each = 2), absorbing_state_new), 
                       rep(c("M", "m"), 
                           (length(initial_state) + 2*length(transient_state_new) + length(absorbing_state))/2), sep="")
  
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
    ddi = ddr[which(ddr$PTNUM==inds[i]),]
    tti = ddi$years[pmatch(c(initial_state_new, rep(transient_state_new, each = 2),
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
      ddi = ddr[which(ddr$PTNUM==inds[i]),]
      tti = ddi$years[pmatch(c(initial_state_new, rep(transient_state_new, each = 2),
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
arrange_data(dd, gg)

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
edge_matrices(gg)

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
names_of_survial_density(gg)

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
passedBy = edge_matrices(gg)$passedBy
possible = edge_matrices(gg)$possible
survival_function(gg, "012", passedBy)

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
travelled = edge_matrices(gg)$traveled
density_function(gg, "012", travelled)


# S_01 = function(param, t){1-pexp(t,param[1])}
# S_12 = function(param, t){1-pexp(t,param[2])}
# S_23 = function(param, t){1-pexp(t,param[3])}
# S_03 = function(param, t){1-pexp(t,param[4])}
# S_13 = function(param, t){1-pexp(t,param[5])}
# 
# f_01 = function(param, t){dexp(t,param[1])}
# f_12 = function(param, t){dexp(t,param[2])}
# f_23 = function(param, t){dexp(t,param[3])}
# f_03 = function(param, t){dexp(t,param[4])}
# f_13 = function(param, t){dexp(t,param[5])}

S_01 = function(param, t){(as.numeric(t>=0))* (1-pweibull(t,param[1],param[2]))}
S_12 = function(param, t){(as.numeric(t>=0))* (1-pweibull(t,param[3],param[4]))}
S_23 = function(param, t){(as.numeric(t>=0))* (1-pweibull(t,param[5],param[6]))}
S_03 = function(param, t){(as.numeric(t>=0))* (1-pweibull(t,param[7],param[8]))}
S_13 = function(param, t){(as.numeric(t>=0))* (1-pweibull(t,param[9],param[10]))}

f_01 = function(param, t){as.numeric(t>=0)*dweibull(t,param[1],param[2])}
f_12 = function(param, t){as.numeric(t>=0)*dweibull(t,param[3],param[4])}
f_23 = function(param, t){as.numeric(t>=0)*dweibull(t,param[5],param[6])}
f_03 = function(param, t){as.numeric(t>=0)*dweibull(t,param[7],param[8])}
f_13 = function(param, t){as.numeric(t>=0)*dweibull(t,param[9],param[10])}



densities <- c(f_01,f_03,f_12,f_13,f_23)
names(densities) <- c("01","03","12","13","23")

survival_functions <- c(S_01,S_03,S_12,S_13,S_23)
names(survival_functions) <- c("01","03","12","13","23") 

## The edges for the transition to the absorbing state written in the "new" way
edge_abs <- c("03","13","23")
edge_mats <- edge_matrices(gg)


# this version of the function is for datasets (or observations) where the time into absorbing state 
# is observed exactly
# for now it is limited to integrals of dimension 3 (but can be extended)
type_to_integrand_absExact = function(type,edge_mats,densfunc,survfunc,edge_abs){
  # inputs: type, edge_mats, density functions, survival functions, vector of edges into absorbing state
  # output: a integrand per type
  one <- function(param,t){1}
  
  placehold0 <- function(tt,param,ft1=one){
    # Spass,Spos: lists of functions
    out <- ft1(param,tt)
    npa <- length(Spass)
    out <- out*prod(sapply(1:npa,function(i) match.fun(Spass[[i]])(param,tt)))
    
    npo <- length(Spos)
    out <- out*prod(sapply(1:npo,function(i) match.fun(Spos[[i]])(param,tt)))
    return(out)
  }
  
  placehold1 <- function(xx,tt,param,ft1,ft2=one){
    # Spass1: passedBy in node 0, Spass2: list of passedBy in node 1
    ss <- xx[1]
    out <- ft1(param,ss)*ft2(param,tt-ss)
    npa1 <- length(Spass1)
    out <- out*prod(sapply(1:npa1,function(i) match.fun(Spass1[[i]])(param,ss)))
    
    npa2 <- length(Spass2)
    out <- out*prod(sapply(1:npa2,function(i) match.fun(Spass2[[i]])(param,tt-ss)))
    
    npo <- length(Spos)
    out <- out*prod(sapply(1:npo,function(i) match.fun(Spos[[i]])(param,tt-ss)))
    return(out)
  }
  
  
  placehold2 <- function(xx,tt,param,ft1, ft2,ft3=one){
    # Spass1: passedBy in node 0, Spass2: list of passedBy in node 1, Spass3: list of passedBy in node 2
    ss <- xx[1]
    uu <- xx[2]
    out <- ft1(param,ss)*ft2(param,uu-ss)*ft3(param,tt-uu)
    npa1 <- length(Spass1)
    out <- out*prod(sapply(1:npa1,function(i) match.fun(Spass1[[i]])(param,ss)))
    
    npa2 <- length(Spass2)
    out <- out*prod(sapply(1:npa2,function(i) match.fun(Spass2[[i]])(param,uu-ss)))
    
    npa3 <- length(Spass3)
    out <- out*prod(sapply(1:npa3,function(i) match.fun(Spass3[[i]])(param,tt-uu)))
    
    npo <- length(Spos)
    out <- out*prod(sapply(1:npo,function(i) match.fun(Spos[[i]])(param,tt-uu)))
    return(out)
  }
  
  placehold3 <- function(xx,tt,param,ft1, ft2, ft3, ft4=one){
    # Spass1: passedBy in node 0, Spass2: list of passedBy in node 1, Spass3: list of passedBy in node 2, Spass4: list of passedBy in node 3
    ss <- xx[1]
    uu <- xx[2]
    rr <- xx[3]
    out <- ft1(param,ss)*ft2(param,uu-ss)*ft3(param,rr-uu)*ft4(param,tt-rr)
    npa1 <- length(Spass1)
    out <- out*prod(sapply(1:npa1,function(i) match.fun(Spass1[[i]])(param,ss)))
    
    npa2 <- length(Spass2)
    out <- out*prod(sapply(1:npa2,function(i) match.fun(Spass2[[i]])(param,uu-ss)))
    
    npa3 <- length(Spass3)
    out <- out*prod(sapply(1:npa3,function(i) match.fun(Spass3[[i]])(param,rr-uu)))
    
    npa4 <- length(Spass4)
    out <- out*prod(sapply(1:npa4,function(i) match.fun(Spass4[[i]])(param,tt-rr)))
    
    npo <- length(Spos)
    out <- out*prod(sapply(1:npo,function(i) match.fun(Spos[[i]])(param,tt-rr)))
    return(out)
  }
  
  
  travi <- edge_mats$traveled[type,]
  passi <- edge_mats$passedBy[type,]
  posi <- edge_mats$possible[type,]
  
  
  if (max(travi)==0){
    Spass=list(one)
    Spos=survival_functions[names(which(posi==1))]
    func_out <- function(tt,param){
      placehold0(tt,param)
    }
  }else if (max(travi)==1 & names(travi)[which.max(travi)] %in% edge_abs){
    # has travelled 1 edge, may have passed by X edges (X may be 0), may have 0 possible edges in the next step (since has come to absorbing)
    Spass=survival_functions[names(which(passi==1))]
    Spos=list(one)
    ft1=densities[names(which(travi==1))][[1]]
    func_out <- function(tt,param){
      placehold0(tt,param,ft1=ft1)
    }
  }else if (max(travi)==1 & !names(travi)[which.max(travi)] %in% edge_abs){
    # has travelled 1 edge, may have passed by X edges (X>=0), may have Y possible edges in the next step (Y>0)
    if (sum(passi>0)==0){
      Spos=survival_functions[names(which(posi==2))]
      Spass1=list(one)
      Spass2=list(one)
    }else{
      Spos=survival_functions[names(which(posi==2))]
      Spass1 =survival_functions[names(which(passi==1))]
      Spass2=list(one)
    }
    ft1 = densities[names(which(travi==1))][[1]]
    func_out <- function(xx,tt,param){
      placehold1(xx,tt,param,ft1=ft1)
    }
    
  }else if (max(travi)==2 & names(travi)[which.max(travi)] %in% edge_abs){
    # has traveled 2 edges, may have passed by X edges (X>=0) in each node, may have 0 possible edges in the next step (because reaches absorbing)
    if(sum(passi==1)==0){
      Spass1=list(one)
    }else{
      Spass1 =survival_functions[names(which(passi==1))]
    }
    if(sum(passi==2)==0){
      Spass2=list(one)
    }else{
      Spass2 =survival_functions[names(which(passi==2))]
    }
    Spos=list(one)
    ft1=densities[names(which(travi==1))][[1]]
    ft2=densities[names(which(travi==2))][[1]]
    func_out <- function(xx,tt,param){
      placehold1(xx,tt,param,ft1=ft1,ft2=ft2)
    }
    
    
  }else if (max(travi)==2 & !names(travi)[which.max(travi)] %in% edge_abs){
    # has traveled 2 edges, may have passed by X edges (X>=0) in each node, may have Y possible edges in the next step (Y>0)
    if(sum(passi==1)==0){
      Spass1=list(one)
    }else{
      Spass1 =survival_functions[names(which(passi==1))]
    }
    if(sum(passi==2)==0){
      Spass2=list(one)
    }else{
      Spass2 =survival_functions[names(which(passi==2))]
    }
    if(sum(passi==3)==0){
      Spass3=list(one)
    }else{
      Spass3 =survival_functions[names(which(passi==3))]
    }
    Spos=survival_functions[names(which(posi==3))]
    ft1=densities[names(which(travi==1))][[1]]
    ft2=densities[names(which(travi==2))][[1]]
    func_out <- function(xx,tt,param){
      placehold2(xx,tt,param,ft1=ft1,ft2=ft2)
    }
    
  }else if (max(travi)==3 & names(travi)[which.max(travi)] %in% edge_abs){
    # has traveled 3 edges, may have passed by X edges (X>=0) in each node, may have 0 possible edges in the next step (because reaches absorbing)
    if(sum(passi==1)==0){
      Spass1=list(one)
    }else{
      Spass1 =survival_functions[names(which(passi==1))]
    }
    if(sum(passi==2)==0){
      Spass2=list(one)
    }else{
      Spass2 =survival_functions[names(which(passi==2))]
    }
    if(sum(passi==3)==0){
      Spass3=list(one)
    }else{
      Spass3 =survival_functions[names(which(passi==3))]
    }
    Spos=list(one)
    ft1=densities[names(which(travi==1))][[1]]
    ft2=densities[names(which(travi==2))][[1]]
    ft3=densities[names(which(travi==3))][[1]]
    func_out <- function(xx,tt,param){
      placehold2(xx,tt,param,ft1=ft1,ft2=ft2,ft3=ft3)
    }
    
    
  }else if (max(travi)==3 & !names(travi)[which.max(travi)] %in% edge_abs){
    # has traveled 3 edges, may have passed by X edges (X>=0) in each node, may have Y possible edges in the next step (Y>0)
    if(sum(passi==1)==0){
      Spass1=list(one)
    }else{
      Spass1 =survival_functions[names(which(passi==1))]
    }
    if(sum(passi==2)==0){
      Spass2=list(one)
    }else{
      Spass2 =survival_functions[names(which(passi==2))]
    }
    if(sum(passi==3)==0){
      Spass3=list(one)
    }else{
      Spass3 =survival_functions[names(which(passi==3))]
    }
    if(sum(passi==4)==0){
      Spass4=list(one)
    }else{
      Spass4 =survival_functions[names(which(passi==4))]
    }
    Spos=survival_functions[names(which(posi==3))]
    ft1=densities[names(which(travi==1))][[1]]
    ft2=densities[names(which(travi==2))][[1]]
    ft3=densities[names(which(travi==3))][[1]]
    func_out <- function(xx,tt,param){
      placehold3(xx,tt,param,ft1=ft1,ft2=ft2,ft3=ft3)
    }
    
  }else if (max(travi)==4  & names(travi)[which.max(travi)] %in% edge_abs){
    # has traveled 4 edges, may have passed by X edges (X>=0) in each node,  may have 0 possible edges in the next step (because reaches absorbing)
    if(sum(passi==1)==0){
      Spass1=list(one)
    }else{
      Spass1 =survival_functions[names(which(passi==1))]
    }
    if(sum(passi==2)==0){
      Spass2=list(one)
    }else{
      Spass2 =survival_functions[names(which(passi==2))]
    }
    if(sum(passi==3)==0){
      Spass3=list(one)
    }else{
      Spass3 =survival_functions[names(which(passi==3))]
    }
    if(sum(passi==4)==0){
      Spass4=list(one)
    }else{
      Spass4 =survival_functions[names(which(passi==4))]
    }
    Spos=list(one)
    ft1=densities[names(which(travi==1))][[1]]
    ft2=densities[names(which(travi==2))][[1]]
    ft3=densities[names(which(travi==3))][[1]]
    ft4=densities[names(which(travi==4))][[1]]
    func_out <- function(xx,tt,param){
      placehold3(xx,tt,param,ft1=ft1,ft2=ft2,ft3=ft3,ft4=ft4)
    }
  }#and so on 
  
  return(func_out)
}

## Part 3: From the time points of a given patient to an integral

all_edges = get.edgelist(gg)
all_absorbing = which(!(all_edges[,2] %in% all_edges[,1]))
absorbing_state = unique(sapply(all_absorbing, function(p) all_edges[p,2]))
all_states_old = state_ordering(gg)[,1]
row_absorbing = sapply(1:length(absorbing_state), function(p) which(all_states_old == absorbing_state[p], arr.ind = TRUE))
absorbing_state_new = sapply(row_absorbing, function(p) state_ordering(gg)[p,2])
all_data_set = arrange_data(data_set = dd, gg)

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

## Finding all types, integrand, time points ++
formula_obs_types = all_types(gg)
all_data_set = arrange_data(data_set = dd, gg)
observation_type = rep(NA, nrow(all_data_set))
all_integral_limits = list()
integrand = list()
all_types = list()
for(i in 1:nrow(all_data_set)){
  observation_type[i] = all_data_set[i,"obs_type"]
  type = names(which(formula_obs_types[, observation_type[i]] == 1))
  integrand_mellomregn = list()
  integral_mellomregn= list()
  type_1 = list()
  for(j in 1:length(type)){
    integrand_mellomregn[[j]] = type_to_integrand_absExact(type[j], edge_mats, densities, survival_functions, edge_abs)
    integral_mellomregn[[j]] = finding_limits_integral(i, type[j], gg, all_edges, absorbing_state_new, all_data_set)
    type_1[[j]] = type[j]
  }
  all_integral_limits[[i]] = integral_mellomregn
  integrand[[i]] = integrand_mellomregn
  all_types[[i]] = type_1
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
        log(integrand[[i]][[1]](tt = tmax, param = param))
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
          calculate_integral_multiple_types[j] = integrand[[i]][[j]](tt = tmax, param = param)
        } else if(length_unique_lower_integral != length_lower_integral){
          calculate_integral_multiple_types[j] = cubintegrate(integrand[[i]][[j]], lower = lower_integral, upper = upper_integral, method = "vegas", maxEval = 1000,
                                                              tt = tmax, param = param)$integral
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


## Optimization
system.time({
optim(rep(1,10), from_time_point_to_integral, method1 = "hcubature", integrand = integrand, all_data_set = all_data_set, 
      all_integral_limits = all_integral_limits, mc_cores = 5, method = "L-BFGS-B", lower = rep(0.00001,10),
            control = list(fnscale = -1))
})

