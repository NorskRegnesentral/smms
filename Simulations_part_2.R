rm(list = ls())
library(Rlab)
library(numDeriv)
library(parallel)

nn = 100

theta_hat_par = matrix(NA, nrow = nn, ncol = 6)
theta_hat_par2 = matrix(NA, nrow = nn, ncol = 6)
theta_hat_solve_hessian = matrix(NA, nrow = nn, ncol = 6)

theta_hat_par_no_cens = matrix(NA, nrow = nn, ncol = 6)
theta_hat_solve_hessian_no_cens = matrix(NA, nrow = nn, ncol = 6)

## Initial survival function
S_01 = function(a, t){1-pweibull(t,a[1], a[2])}
S_12 = function(a, t){1-pweibull(t,a[3], a[4])}
S_23 = function(a, t){1-pweibull(t,a[5], a[6])}

## Initial density functions - derivative of -survival function
f_01 = function(a, t){dweibull(t,a[1],a[2])}
f_12 = function(a, t){dweibull(t,a[3],a[4])}
f_23 = function(a, t){dweibull(t,a[5],a[6])}

simulation_t = function(param, n){
  T_01 = rep(NA, n)
  T_12 = rep(NA, n)
  T_23 = rep(NA, n)
  
  t_type0 = c()
  t_type1 = matrix(, nrow = 0, ncol = 2)
  t_type2 = matrix(, nrow = 0, ncol = 3)
  t_type3 = matrix(, nrow = 0, ncol = 3)
  
  for(i in 1:n){
    unif1 = runif(1, min = 0, max = 1)
    S01 = function(T)
    {
      1-pweibull(T, param[1], param[2]) - unif1
    }
    T_01[i] = uniroot(S01, interval = c(1.e-14, 1e04), tol = 1e-9)$root
    unif2 = runif(1, min = 0, max = 1)
    S12 = function(T)
    {
      1-pweibull(T, param[3], param[4]) - unif2
    }
    T_12[i] = uniroot(S12, interval = c(1.e-14, 1e04), tol = 1e-9)$root
    unif3 = runif(1, min = 0, max = 1)
    S23 = function(T)
    {
      1-pweibull(T, param[5], param[6]) - unif3
    }
    T_23[i] = uniroot(S23, interval = c(1.e-14, 1e04), tol = 1e-9)$root
    t1 = 0.000001
    t2 = runif(1, t1, 2)
    t3 = runif(1, t1, 5)
    t4 = runif(1, t1, 8)
    t5 = runif(1, t1, 10)
    t6 = runif(1, t1, 12)
    t7 = runif(1, t1, 15)
    t8 = runif(1, t1, 18)
    t9 = runif(1, t1, 20)
    t10 = runif(1, t1, 22)
    t = c(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10)
    if(all(t[1:9] - t[2:10]<0))
    {
      t = t[1:10]
    } else if (all(t[1:8] - t[2:9]<0))
    {
      t = t[1:9]
    } else if (all(t[1:7] - t[2:8]<0))
    {
      t = t[1:8]
    } else if (all(t[1:6] - t[2:7]<0))
    {
      t = t[1:7]
    } else if (all(t[1:5] - t[2:6]<0))
    {
      t = t[1:6]
    } else if (all(t[1:4] - t[2:5]<0))
    {
      t = t[1:5]
    } else if (all(t[1:3] - t[2:4]<0))
    {
      t = t[1:4]
    } else if (all(t[1:2] - t[2:3]<0))
    {
      t = t[1:3]
    } else if (all(t[1:1] - t[2:2]<0))
    {
      t = t[1:2]
    }
    t_n = length(t)
    if(t[t_n] < T_01[i]){
      t_type0 = c(t_type0, t[t_n])
    } else if(any(t > T_01[i]) & T_01[i] + T_12[i] > t[t_n]){
      upper_int = min(which(t > T_01[i]))
      t_type1 = rbind(t_type1, c(t[upper_int], t[t_n]))
    } else if(any(t > T_01[i] + T_12[i]) & T_01[i] + T_12[i] + T_23[i] > t[t_n]){
      upper_int1 = min(which(t > T_01[i]))
      upper_int2 = min(which(t > T_01[i] + T_12[i]))
      t_type2 = rbind(t_type2, c(t[upper_int1], t[upper_int2], t[t_n]))
    } else if(any(t > T_01[i] + T_12[i]) & any(t > T_01[i] + T_12[i] + T_23[i])){
      upper_int1 = min(which(t > T_01[i]))
      upper_int2 = min(which(t > T_01[i] + T_12[i]))
      t_type3 = rbind(t_type3, c(min(t[upper_int1], T_01[i] + T_12[i] + T_23[i]), min(t[upper_int2],T_01[i] + T_12[i] + T_23[i]), 
                                 T_01[i] + T_12[i] + T_23[i]))
    }
  }
  return(list(t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, t_type3 = t_type3))
}


sum_all_no_censoring = function(a, t_type0, t_type1, t_type2, t_type3, mc_cores){
  type_1_no_censoring = -sum(sapply(t_type0, function(t) (log(S_01(a, t)))))
  type_2_no_censoring = -sum(unlist(mclapply(1:nrow(t_type1), function(i) (log(f_01(a, t_type1[i,1])*S_12(a,t_type1[i,2]-t_type1[i,1]))), mc.cores = mc_cores)))
  type_3_no_censoring = -sum(unlist(mclapply(1:nrow(t_type2), function(i) {
    return(log(f_01(a, t_type2[i,1])*f_12(a,t_type2[i,2]-t_type2[i,1])*S_23(a,t_type2[i,3]-t_type2[i,2])))
  },
  mc.cores = mc_cores)))
  type_4_no_censoring = -sum(unlist(mclapply(1:nrow(t_type3), function(i) {
    return(log(f_01(a, t_type3[i,1])*f_12(a,t_type3[i,2]-t_type3[i,1])*f_23(a,t_type3[i,3]-t_type3[i,2])))}, 
    mc.cores = mc_cores)))
  return(type_1_no_censoring + type_2_no_censoring + type_3_no_censoring + type_4_no_censoring) 
}

for(mn in 1:nn){
  set.seed(mn)
  print(mn)
  eps = 0.001
  param = c(5, 1, 3, 2, 4, 8)
  n = 1000
  t = simulation_t(param, n)
  
  t_type0 = t$t_type0
  t_type1 = t$t_type1
  t_type2 = t$t_type2
  t_type3 = t$t_type3
  
  print(length(t_type0) + nrow(t_type1) + nrow(t_type2) + nrow(t_type3))
  for(i in 1:nrow(t_type2)){
    if(t_type2[i,1] == t_type2[i,2]){
      t_type2[i,1] = (t_type2[i,2]/2)
    }
  }
  for(i in 1:nrow(t_type3)){
    if(t_type3[i,1] == t_type3[i,2] & t_type3[i,2] == t_type3[i,3]){
      t_type3[i,1] = t_type3[i,3]/3
      t_type3[i,2] = 2*t_type3[i,3]/3
    }else if(t_type3[i,2] == t_type3[i,1]){
      t_type3[i,1] = t_type3[i,2]/2
    } else if(t_type3[i,2] == t_type3[i,3]){
      t_type3[i,2] = t_type3[i,3] - (t_type3[i,3]-t_type3[i,1])/2
    }
  }
  
  theta_hat_no_cens = optim(c(3, 1, 2, 1, 2, 6), sum_all_no_censoring, t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, t_type3 = t_type3, mc_cores = 1, method = "L-BFGS-B",
                            lower = c(eps, eps, eps, eps, eps, eps))
  
  theta_hat_middle_no_cens = as.matrix(theta_hat_no_cens$par)
  theta_hat_par_no_cens[mn, 1] = theta_hat_middle_no_cens[1, 1]
  theta_hat_par_no_cens[mn, 2] = theta_hat_middle_no_cens[2, 1]
  theta_hat_par_no_cens[mn, 3] = theta_hat_middle_no_cens[3, 1]
  theta_hat_par_no_cens[mn, 4] = theta_hat_middle_no_cens[4, 1]
  theta_hat_par_no_cens[mn, 5] = theta_hat_middle_no_cens[5, 1]
  theta_hat_par_no_cens[mn, 6] = theta_hat_middle_no_cens[6, 1]
  
  theta_hat_hessian_no_cens = hessian(sum_all_no_censoring, theta_hat_no_cens$par, t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, t_type3 = t_type3, mc_cores = 1)
  theta_hat_hessian_solved_no_cens = solve(theta_hat_hessian_no_cens)
  theta_hat_solve_hessian_no_cens[mn, 1] = theta_hat_hessian_solved_no_cens[1, 1]
  theta_hat_solve_hessian_no_cens[mn, 2] = theta_hat_hessian_solved_no_cens[2, 2]
  theta_hat_solve_hessian_no_cens[mn, 3] = theta_hat_hessian_solved_no_cens[3, 3]
  theta_hat_solve_hessian_no_cens[mn, 4] = theta_hat_hessian_solved_no_cens[4, 4]
  theta_hat_solve_hessian_no_cens[mn, 5] = theta_hat_hessian_solved_no_cens[5, 5]
  theta_hat_solve_hessian_no_cens[mn, 6] = theta_hat_hessian_solved_no_cens[6, 6]
  
  
}