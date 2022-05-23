
gg = graph_from_literal("1"--+"2"--+"3"--+"4", "1"--+"4", "2"--+"4")
edge_mats <- edge_matrices(gg)

densities <- c(f_01,f_03,f_12,f_13,f_23)
names(densities) <- c("01","03","12","13","23")

survival_functions <- c(S_01,S_03,S_12,S_13,S_23)
names(survival_functions) <- c("01","03","12","13","23") 

edge_abs <- c("03","13","23")

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
