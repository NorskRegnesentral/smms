# Testing preprocessing functions

rm(list = ls())
devtools::load_all() 
library(igraph)

gg1 = graph_from_literal("1"--+"2"--+"3"--+"4", "1"--+"4", "2"--+"4")
gg2 = graph_from_literal("1"--+"2"--+"3", "1"--+"4", "2"--+"4","1"--+"3")
gg3 = graph_from_literal("A"--+"B","C"--+"D","B"--+"E","D"--+"B","D"--+"E","C"--+"B")
gg4 = graph_from_literal("A"--+"B"--+"C"--+"E","A"--+"E","A"--+"D"--+"E")
gg5 = graph_from_literal("A"--+"B"--+"C"--+"D"--+"E","A"--+"F"--+"G"--+"E","F"--+"B","F"--+"D")

state_ordering(gg1)
#state order  type
#1     1     0  init
#2     2     1 trans
#3     3     2 trans
#4     4     3   abs

state_ordering(gg2)
#state order  type
#1     1     0  init
#2     2     1 trans
#3     3     2   abs
#4     4     3   abs

state_ordering(gg3) 

state_ordering(gg4) 
#state order  type
#1     A     0  init
#2     B     1 trans
#3     C     3 trans
#4     E     4   abs
#5     D     2 trans

state_ordering(gg5) 

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

dd3a <- data.frame(patient=c(rep(1,4),rep(2,5)),time=c(1,1.5,2,2.5,0,1,2,3,4),state=c("B","B","B","B","A","A","B","E","E"))

relevant_timepoints(dd3a,gg3)

dd4 <- data.frame(patient=c(rep(1,4),rep(2,5)),time=c(1,1.5,2,2.5,0,1,2,3,4),state=c("B","C","C","E","A","A","D","E","E"))

relevant_timepoints(dd4,gg4)

##
arrange_data(dd1,gg1)
arrange_data(dd2,gg2)
arrange_data(dd3,gg3) 
arrange_data(dd4,gg4)

## 
edge_matrices(gg1)
edge_matrices(gg2)
edge_matrices(gg3)
edge_matrices(gg4) 
edge_matrices(gg5) 
