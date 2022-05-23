# Testing preprocessing functions

source("sm_msm_preprocessing_func.R")
library(igraph)

gg1 = graph_from_literal("1"--+"2"--+"3"--+"4", "1"--+"4", "2"--+"4")
gg2 = graph_from_literal("1"--+"2"--+"3", "1"--+"4", "2"--+"4","1"--+"3")
gg3 = graph_from_literal("A"--+"B","C"--+"D","B"--+"E","D"--+"B")
gg4 = graph_from_literal("A"--+"B"--+"C"--+"E","A"--+"E","A"--+"D"--+"E")

state_ordering(gg1)

state_ordering(gg2)

state_ordering(gg3) 

state_ordering(gg4) 

construct_formula_types(gg1)
construct_formula_types(gg2)
construct_formula_types(gg3) 
construct_formula_types(gg4)

construct_obs_types(gg1)
construct_obs_types(gg2)
construct_obs_types(gg3) 
construct_obs_types(gg4)

all_types(gg1)
all_types(gg2)
all_types(gg3) # we assume that we now which initial state patients come from!
all_types(gg4)

##
dd1 <- data.frame(patient=c(rep(1,4),rep(2,5)),time=c(1,1.5,2,2.5,0,1,2,3,4),state=c(2,2,2,3,1,1,3,4,4))

relevant_timepoints(dd1,gg1)

dd2 <- data.frame(patient=c(rep(1,4),rep(2,5)),time=c(1,1.5,2,2.5,0.5,1,2,3,4),state=c(1,1,2,2,2,2,2,4,4))

relevant_timepoints(dd2,gg2) 

dd3 <- data.frame(patient=c(rep(1,4),rep(2,5)),time=c(1,1.5,2,2.5,0,1,2,3,4),state=c("D","D","D","B","A","A","B","E","E"))

relevant_timepoints(dd3,gg3)

##
