# types_to_integrands for all formula types at once

# Define placeholder functions (numbering indicates integral dimension)
# here I assume a maximum of 3 edges going out of each node (i.e. max two possible passedBy from each node traveled,
# and 3 possible possible nodes from each node)
one <- function(param,t){1}

placehold0 <- function(tt,param,ft1=one, Spa1=one, Spa2=one, Spo1=one, Spo2=one, Spo3=one){
  ft1(param,tt)*Spa1(param,tt)*Spa2(param,tt)*Spo1(param,tt)*Spo2(param,tt)*Spo3(param,tt)
}

placehold1 <- function(xx,tt,param,ft1,ft2=one,Spa1_1=one, Spa1_2=one, Spa2_1=one, Spa2_2=one, 
                       Spo1=one, Spo2=one, Spo3=one){
  ss <- xx[1]
  (ft1(param,ss)*ft2(param,tt-ss)*Spa1_1(param,ss)*Spa1_2(param,ss)*Spa2_1(param,tt-ss)*Spa2_2(param,tt-ss)*
      Spo1(param,tt-ss)*Spo2(param,tt-ss)*Spo3(param,tt-ss))
}

placehold2 <- function(xx,tt,param,ft1, ft2,ft3=one, Spa1_1=one, Spa1_2=one, Spa2_1=one, Spa2_2=one, 
                       Spa3_1=one, Spa3_2=one, Spo1=one, Spo2=one, Spo3=one){
  ss <- xx[1]
  uu <- xx[2]
  (ft1(param,ss)*ft2(param,uu-ss)*ft3(param,tt-uu)*Spa1_1(param,ss)*Spa1_2(param,ss)*
      Spa2_1(param,uu-ss)*Spa2_2(param,uu-ss)*Spa3_1(param,tt-uu)*Spa3_2(param,tt-uu)*
      Spo1(param,tt-uu)*Spo2(param,tt-uu)*Spo3(param,tt-uu))
}

placehold3 <- function(xx,tt,param,ft1, ft2, ft3, ft4=one, Spa1_1=one, Spa1_2=one, Spa2_1=one, Spa2_2=one, 
                       Spa3_1=one, Spa3_2=one, Spa4_1=one, Spa4_2=one, Spo1=one, Spo2=one, Spo3=one){
  ss <- xx[1]
  uu <- xx[2]
  rr <- xx[3]
  (ft1(param,ss)*ft2(param,uu-ss)*ft3(param,rr-uu)*ft4(param,tt-rr)*Spa1_1(param,ss)*
      Spa1_2(param,ss)*Spa2_1(param,uu-ss)*Spa2_2(param,uu-ss)*
      Spa3_1(param,rr-uu)*Spa3_2(param,rr-uu)*Spa4_1(param,tt-rr)*Spa4_2(param,tt-rr)*
      Spo1(param,tt-rr)*Spo2(param,tt-rr)*Spo3(param,tt-rr))
}

types_to_integrands = function(edge_mats){
  # inputs: edge_mats, density functions, survival functions, vector of edges into absorbing state
  # output: a vector nf functions (where nf is the number of formula types)
  
  trav <- edge_mats$traveled
  pass <- edge_mats$passedBy
  pos <- edge_mats$possible
  
  nf <- dim(trav)[1]
  ne <- dim(trav)[2]
  
  edge_abs <- c("03","13","23")
  
  funcs_out <- list()
  
  for (i in 1:nf){
    travi <- trav[i,] # maybe sorted
    passi <- pass[i,]
    posi <- pos[i,]
    
    if (max(trav[i,])==0){
      if (sum(posi>0)==1){
        func_out <- function(tt,param){
          placehold0(tt,param,Spo1=survival_functions[names(which(posi==1))][[1]])
        }
      }else if (sum(posi>0)==2){
        func_out <- function(tt,param){
          placehold0(tt,param,Spo1=survival_functions[names(which(posi==1))][[1]],
                               Spo2=survival_functions[names(which(posi==1))][[2]])
        }
      }else if (sum(posi>0)==3){
        func_out <- function(tt,param){
          placehold0(tt,param,Spo1=survival_functions[names(which(posi==1))][[1]],
                               Spo2=survival_functions[names(which(posi==1))][[2]],
                               Spo3=survival_functions[names(which(posi==1))][[3]])
        }
      }
      
    }else if (max(trav[i,])==1 & colnames(trav)[which.max(trav[i,])] %in% edge_abs){
      # has travelled 1 edge, may have passed by X edges (X may be 0), may have 0 possible edges in the next step (since has come to absorbing)
      if (sum(passi>0)==0){
        func_out <- function(tt,param){
          placehold0(tt,param,ft1=densities[names(which(travi==1))][[1]])
        }
      }else if (sum(passi>0)==1){
        func_out <- function(tt,param){
          placehold0(tt,param,ft1=densities[names(which(travi==1))][[1]],
                     Spa1=survival_functions[names(which(passi==1))][[1]])
        }
      }else if (sum(passi>0)==2){
        func_out <- function(tt,param){
          placehold0(tt,param,ft1=densities[names(which(travi==1))][[1]],
                               Spa1=survival_functions[names(which(passi==1))][[1]],
                               Spa2=survival_functions[names(which(passi==1))][[2]])
        }
      }
      
    }else if (max(trav[i,])==1 & !colnames(trav)[which.max(trav[i,])] %in% edge_abs){
      # has travelled 1 edge, may have passed by X edges (X>=0), may have Y possible edges in the next step (Y>0)
      if (sum(passi>0)==0 & sum(posi>0)==1){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     Spo1=survival_functions[names(which(posi==2))][[1]])
        }
      }else if (sum(passi>0)==0 & sum(posi>0)==2){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     Spo1=survival_functions[names(which(posi==2))][[1]],
                     Spo2=survival_functions[names(which(posi==2))][[2]])
        }
      }else if (sum(passi>0)==0 & sum(posi>0)==3){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     Spo1=survival_functions[names(which(posi==2))][[1]],
                     Spo2=survival_functions[names(which(posi==2))][[2]],
                     Spo3=survival_functions[names(which(posi==2))][[3]])
        }
      }else if (sum(passi>0)==1 & sum(posi>0)==1){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     Spa1=survival_functions[names(which(passi==1))][[1]],
                     Spo1=survival_functions[names(which(posi==2))][[1]])
        }
      }else if (sum(passi>0)==1 & sum(posi>0)==2){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     Spa1=survival_functions[names(which(passi==1))][[1]],
                     Spo1=survival_functions[names(which(posi==2))][[1]],
                     Spo2=survival_functions[names(which(posi==2))][[2]])
        }
      }else if (sum(passi>0)==1 & sum(posi>0)==3){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     Spa1=survival_functions[names(which(passi==1))][[1]],
                     Spo1=survival_functions[names(which(posi==2))][[1]],
                     Spo2=survival_functions[names(which(posi==2))][[2]],
                     Spo3=survival_functions[names(which(posi==2))][[3]])
        }
      }else if (sum(passi>0)==2 & sum(posi>0)==1){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     Spa1=survival_functions[names(which(passi==1))][[1]],
                     Spa2=survival_functions[names(which(passi==1))][[2]],
                     Spo1=survival_functions[names(which(posi==2))][[1]])
        }
      }else if (sum(passi>0)==2 & sum(posi>0)==2){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     Spa1=survival_functions[names(which(passi==1))][[1]],
                     Spa2=survival_functions[names(which(passi==1))][[2]],
                     Spo1=survival_functions[names(which(posi==2))][[1]],
                     Spo2=survival_functions[names(which(posi==2))][[2]])
        }
      }else if (sum(passi>0)==2 & sum(posi>0)==3){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     Spa1=survival_functions[names(which(passi==1))][[1]],
                     Spa2=survival_functions[names(which(passi==1))][[2]],
                     Spo1=survival_functions[names(which(posi==2))][[1]],
                     Spo2=survival_functions[names(which(posi==2))][[2]],
                     Spo3=survival_functions[names(which(posi==2))][[3]])
        }
      }
      
    }else if (max(trav[i,])==2 & colnames(trav)[which.max(trav[i,])] %in% edge_abs){
    # has traveled 2 edges, may have passed by X edges (X>=0) in each node, may have 0 possible edges in the next step (because reaches absorbing)
      if (sum(passi==1)==0 & sum(passi==2)==0){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     ft2=densities[names(which(travi==2))][[1]])
        }
      }else if (sum(passi==1)==0 & sum(passi==2)==1){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     ft2=densities[names(which(travi==2))][[1]],
                     Spa2_1=survival_functions[names(which(passi==2))][[1]])
        }
      }else if (sum(passi==1)==0 & sum(passi==2)==2){
        func_out <- function(tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     ft2=densities[names(which(travi==2))][[1]],
                     Spa2_1=survival_functions[names(which(passi==2))][[1]],
                     Spa2_2=survival_functions[names(which(passi==2))][[2]])
        }
      }else if (sum(passi==1)==1 & sum(passi==2)==0){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     ft2=densities[names(which(travi==2))][[1]],
                     Spa1_1=survival_functions[names(which(passi==1))][[1]])
        }
      }else if (sum(passi==1)==1 & sum(passi==2)==1){
        func_out <- function(xx,tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     ft2=densities[names(which(travi==2))][[1]],
                     Spa1_1=survival_functions[names(which(passi==1))][[1]],
                     Spa2_1=survival_functions[names(which(passi==2))][[1]])
        }
      }else if (sum(passi==1)==1 & sum(passi==2)==2){
        func_out <- function(tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     ft2=densities[names(which(travi==2))][[1]],
                     Spa1_1=survival_functions[names(which(passi==1))][[1]],
                     Spa2_1=survival_functions[names(which(passi==2))][[1]],
                     Spa2_2=survival_functions[names(which(passi==2))][[2]])
        }
      }else if (sum(passi==1)==2 & sum(passi==2)==0){
        func_out <- function(tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     ft2=densities[names(which(travi==2))][[1]],
                     Spa1_1=survival_functions[names(which(passi==1))][[1]],
                     Spa1_2=survival_functions[names(which(passi==1))][[2]])
        }
      }else if (sum(passi==1)==2 & sum(passi==2)==1){
        func_out <- function(tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     ft2=densities[names(which(travi==2))][[1]],
                     Spa1_1=survival_functions[names(which(passi==1))][[1]],
                     Spa1_2=survival_functions[names(which(passi==1))][[2]],
                     Spa2_1=survival_functions[names(which(passi==2))][[1]])
        }
      }else if (sum(passi==1)==2 & sum(passi==2)==2){
        func_out <- function(tt,param){
          placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                     ft2=densities[names(which(travi==2))][[1]],
                     Spa1_1=survival_functions[names(which(passi==1))][[1]],
                     Spa1_2=survival_functions[names(which(passi==1))][[2]],
                     Spa2_1=survival_functions[names(which(passi==2))][[1]],
                     Spa2_2=survival_functions[names(which(passi==2))][[2]])
        }
      }
      
    }else if (max(trav[i,])==2 & !colnames(trav)[which.max(trav[i,])] %in% edge_abs){
      # has traveled 2 edges, may have passed by X edges (X>=0) in each node, may have Y possible edges in the next step (Y>0)
      if (sum(posi>0)==1){
        if (sum(passi==1)==0 & sum(passi==2)==0){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]])
          }
        }else if (sum(passi==1)==0 & sum(passi==2)==1){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]])
          }
        }else if (sum(passi==1)==0 & sum(passi==2)==2){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spa2_2=survival_functions[names(which(passi==2))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]])
          }
        }else if (sum(passi==1)==1 & sum(passi==2)==0){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]])
          }
        }else if (sum(passi==1)==1 & sum(passi==2)==1){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]])
          }
        }else if (sum(passi==1)==1 & sum(passi==2)==2){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spa2_2=survival_functions[names(which(passi==2))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]])
          }
        }else if (sum(passi==1)==2 & sum(passi==2)==0){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa1_2=survival_functions[names(which(passi==1))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]])
          }
        }else if (sum(passi==1)==2 & sum(passi==2)==1){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa1_2=survival_functions[names(which(passi==1))][[2]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]])
          }
        }else if (sum(passi==1)==2 & sum(passi==2)==2){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa1_2=survival_functions[names(which(passi==1))][[2]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spa2_2=survival_functions[names(which(passi==2))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]])
          }
        }
      }else if(sum(posi>0)==2){
        if (sum(passi==1)==0 & sum(passi==2)==0){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]])
          }
        }else if (sum(passi==1)==0 & sum(passi==2)==1){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]])
          }
        }else if (sum(passi==1)==0 & sum(passi==2)==2){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spa2_2=survival_functions[names(which(passi==2))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]])
          }
        }else if (sum(passi==1)==1 & sum(passi==2)==0){
          func_out <- function(xx,tt,param){
            placehold1(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]])
          }
        }else if (sum(passi==1)==1 & sum(passi==2)==1){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]])
          }
        }else if (sum(passi==1)==1 & sum(passi==2)==2){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spa2_2=survival_functions[names(which(passi==2))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]])
          }
        }else if (sum(passi==1)==2 & sum(passi==2)==0){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa1_2=survival_functions[names(which(passi==1))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]])
          }
        }else if (sum(passi==1)==2 & sum(passi==2)==1){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa1_2=survival_functions[names(which(passi==1))][[2]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]])
          }
        }else if (sum(passi==1)==2 & sum(passi==2)==2){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa1_2=survival_functions[names(which(passi==1))][[2]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spa2_2=survival_functions[names(which(passi==2))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]])
          }
        }
      }else if(sum(posi>0)==3){
        if (sum(passi==1)==0 & sum(passi==2)==0){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]],
                       Spo3=survival_functions[names(which(posi==3))][[3]])
          }
        }else if (sum(passi==1)==0 & sum(passi==2)==1){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]],
                       Spo3=survival_functions[names(which(posi==3))][[3]])
          }
        }else if (sum(passi==1)==0 & sum(passi==2)==2){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spa2_2=survival_functions[names(which(passi==2))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]],
                       Spo3=survival_functions[names(which(posi==3))][[3]])
          }
        }else if (sum(passi==1)==1 & sum(passi==2)==0){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]],
                       Spo3=survival_functions[names(which(posi==3))][[3]])
          }
        }else if (sum(passi==1)==1 & sum(passi==2)==1){
          func_out <- function(xx,tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]],
                       Spo3=survival_functions[names(which(posi==3))][[3]])
          }
        }else if (sum(passi==1)==1 & sum(passi==2)==2){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spa2_2=survival_functions[names(which(passi==2))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]],
                       Spo3=survival_functions[names(which(posi==3))][[3]])
          }
        }else if (sum(passi==1)==2 & sum(passi==2)==0){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa1_2=survival_functions[names(which(passi==1))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]],
                       Spo3=survival_functions[names(which(posi==3))][[3]])
          }
        }else if (sum(passi==1)==2 & sum(passi==2)==1){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa1_2=survival_functions[names(which(passi==1))][[2]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]],
                       Spo3=survival_functions[names(which(posi==3))][[3]])
          }
        }else if (sum(passi==1)==2 & sum(passi==2)==2){
          func_out <- function(tt,param){
            placehold2(xx,tt,param,ft1=densities[names(which(travi==1))][[1]],
                       ft2=densities[names(which(travi==2))][[1]],
                       Spa1_1=survival_functions[names(which(passi==1))][[1]],
                       Spa1_2=survival_functions[names(which(passi==1))][[2]],
                       Spa2_1=survival_functions[names(which(passi==2))][[1]],
                       Spa2_2=survival_functions[names(which(passi==2))][[2]],
                       Spo1=survival_functions[names(which(posi==3))][[1]],
                       Spo2=survival_functions[names(which(posi==3))][[2]],
                       Spo3=survival_functions[names(which(posi==3))][[3]])
          }
        }
      }
     
      
    }else if (max(trav[i,])==3 & colnames(trav)[which.max(trav[i,])] %in% edge_abs){
      # has traveled 3 edges, may have passed by X edges (X>=0) in each node, may have 0 possible edges in the next step (because reaches absorbing)
      # FIX. will need 3^3 = 27 categories
    }else if (max(trav[i,])==3 & !colnames(trav)[which.max(trav[i,])] %in% edge_abs){
      # has traveled 3 edges, may have passed by X edges (X>=0) in each node, may have Y possible edges in the next step (Y>0)
      # FIX. will need 3^4 = 81 categories
    }else if (max(trav[i,])==4  & colnames(trav)[which.max(trav[i,])] %in% edge_abs){
      # has traveled 4 edges, may have passed by X edges (X>=0) in each node,  may have 0 possible edges in the next step (because reaches absorbing)
      # FIX. will need 3^4 = 81 categories
    }else if (max(trav[i,])==4  & !colnames(trav)[which.max(trav[i,])] %in% edge_abs){
      # has traveled 4 edges, may have passed by X edges (X>=0) in each node, may have Y possible edges in the next step (Y>0)
      # FIX. will need 3^5 = 243 categories
    } #and so on
    funcs_out[[i]] <- func_out
  }
  names(funcs_out) <- rownames(trav)
  return(funcs_out)
}


#### Testing stuff

# types_to_integrands(edge_mats)
# 
# ff <- type_to_integrand("01",gg)
# 
# fo <- function(ss,tt,param){
#   f_01(param,ss)*S_03(param,ss)*S_12(param,tt-ss)*S_13(param,tt-ss)
# }
# 
# one <- function(param,tt){1}
# 
# fp <- function(tt,f1=one,f2=one,f3=one,param){
#   f1(param,tt)*f2(param,tt)*f3(param,tt)
# }
# 
# 
# s1 <- rep(NA,10^3)
# s2 <- rep(NA,10^3)
# s3 <- rep(NA,10^3)
# 
# for (i in 1:10^3){
#   ss <- runif(1,1,3)
#   ts <- ss+runif(1,0.1,2)
#   ps <- runif(5,0.5,2)
#   s1[i] <- system.time(fo(ss,ts,ps))[1]
#   s2[i] <- system.time(ff(ss,ts,ps))[1]
#   s3[i] <- system.time(func_out(ss,ts,ps))[1]
# }
# 
# sum(s1,na.rm=T)
# sum(s2,na.rm=T)
# sum(s3,na.rm=T)
