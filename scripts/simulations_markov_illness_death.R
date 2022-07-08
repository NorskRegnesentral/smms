rm(list = ls())
library(Rlab)
library(numDeriv)
library(parallel)

nn = 100

theta_hat_par = matrix(NA, nrow = nn, ncol = 6)
theta_hat_solve_hessian = matrix(NA, nrow = nn, ncol = 6)

theta_hat_par_markov = matrix(NA, nrow = nn, ncol = 3)
theta_hat_solve_hessian_markov = matrix(NA, nrow = nn, ncol = 3)

## Initial survival function
S_01 = function(a, t){1-pweibull(t,a[1], a[2])}
S_12 = function(a, t){1-pweibull(t,a[3], a[4])}
S_02 = function(a, t){1-pweibull(t,a[5], a[6])}

## Initial density functions - derivative of -survival function
f_01 = function(a,  t){dweibull(t,a[1], a[2])}
f_12 = function(a, t){dweibull(t,a[3], a[4])}
f_02 = function(a, t){dweibull(t,a[5], a[6])}

## Initial survival function
S_01_exp = function(a, t){1-pexp(t,a[1])}
S_12_exp = function(a, t){1-pexp(t,a[2])}
S_02_exp = function(a, t){1-pexp(t,a[3])}

## Initial density functions - derivative of -survival function
f_01_exp = function(a, t){dexp(t,a[1])}
f_12_exp = function(a, t){dexp(t,a[2])}
f_02_exp = function(a, t){dexp(t,a[3])}

## Other functions
## Not exponential
f01_S12_S02 = function(a,  w, t){f_01(a, t)*S_12(a, w - t)*S_02(a,t)}
f01_f12_S02 = function(a, w, t){f_01(a, t)*f_12(a, w - t)*S_02(a,t)}

## Exponential
f01_S12_S02_exp = function(a, w, t){f_01_exp(a, t)*S_12_exp(a, w - t)*S_02_exp(a, t)}
f01_f12_S02_exp = function(a, w, t){f_01_exp(a, t)*f_12_exp(a, w - t)*S_02_exp(a, t)}

## Simulating the time points
simulation_t = function(param, n){
  T_01 = rep(NA, n)
  T_12 = rep(NA, n)
  T_02 = rep(NA, n)
  
  t_type0 = matrix(, nrow = 0, ncol = 1)
  t_type1 = matrix(, nrow = 0, ncol = 3)
  t_type2 = matrix(, nrow = 0, ncol = 3)
  t_type3 = matrix(, nrow = 0, ncol = 1)
  
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
    S02 = function(T)
    {
      1-pweibull(T, param[5], param[6]) - unif3
    }
    T_02[i] = uniroot(S02, interval = c(1.e-14, 1e04), tol = 1e-9)$root
    
    #mm = sample(2:10,1)
    #t = c(0,sort(rlnorm(mm-1,-1,1)))
    mm = sample(2:20,1)
    t = c(0,sort(rlnorm(mm-2,-1,2)))
    t_n = length(t)
    if(t[t_n] < T_01[i] & t[t_n] < T_02[i]){
      t_type0 = rbind(t_type0, t[t_n])
    } else if(T_01[i] < T_02[i]){
      if(any(t > T_01[i]) & T_01[i] + T_12[i] > t[t_n]){
        lower_int = max(which(t < T_01[i]))
        upper_int = min(which(t > T_01[i]))
        t_type1 = rbind(t_type1, c(t[lower_int], t[upper_int], t[t_n]))
      } else if(any(t > T_01[i] + T_12[i])){
        lower_int1 = max(which(t < T_01[i]))
        upper_int1 = min(which(t > T_01[i]))
        t_type2 = rbind(t_type2, c(t[lower_int1], min(T_01[i] + T_12[i],t[upper_int1]), 
                                   T_01[i] + T_12[i]))
      }
    } else if(any(t > T_02[i]) & T_02[i] < T_01[i]){
      t_type3 = rbind(t_type3, T_02[i])
    }
  }
  return(list(t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, t_type3 = t_type3))
}

sum_all = function(a, t_type0, t_type1, t_type2, t_type3, mc_cores){
  type_1 = -sum(unlist(mclapply(1:nrow(t_type0), function(i) (log(S_01(a = a,  t = t_type0[i, 1])*S_02(a = a,  t = t_type0[i, 1]))))))
  type_2 = -sum(unlist(mclapply(1:nrow(t_type1), function(i)(log(as.numeric(
    integrate(f01_S12_S02, lower = t_type1[i,1], upper = t_type1[i,2], w = t_type1[i,3], 
              a = a)[1]))),mc.cores = mc_cores), use.names = FALSE))
  type_3 = -sum(unlist(mclapply(1:nrow(t_type2), function(i)(log(as.numeric(
    integrate(f01_f12_S02, lower = t_type2[i,1], upper = t_type2[i,2], 
              w = t_type2[i,3], a = a)[1]))),mc.cores = mc_cores), use.names = FALSE))
  type_4 = -sum(unlist(mclapply(1:nrow(t_type3), function(i) (log(f_02(a = a,  t = t_type3[i, 1])*
                                                                    S_01(a = a,  t = t_type3[i, 1]))))))
  
  return(type_1 + type_2 + type_3 + type_4) 
}

sum_all_markov = function(a, t_type0, t_type1, t_type2, t_type3, mc_cores){
  type_1_exp = -sum(unlist(mclapply(1:nrow(t_type0), function(i) (log(S_01_exp(a = a,  t = t_type0[i, 1])*S_02_exp(a = a,  t = t_type0[i, 1]))))))
  type_2_exp = -sum(unlist(mclapply(1:nrow(t_type1), function(i)(log(as.numeric(
    integrate(f01_S12_S02_exp, lower = t_type1[i,1], upper = t_type1[i,2], w = t_type1[i,3], 
              a = a)[1]))),mc.cores = mc_cores), use.names = FALSE))
  type_3_exp = -sum(unlist(mclapply(1:nrow(t_type2), function(i)(log(as.numeric(
    integrate(f01_f12_S02_exp, lower = t_type2[i,1], upper = t_type2[i,2], 
              w = t_type2[i,3], a = a)[1]))),mc.cores = mc_cores), use.names = FALSE))
  type_4_exp = -sum(unlist(mclapply(1:nrow(t_type3), function(i) (log(f_02_exp(a = a,  t = t_type3[i, 1])*
                                                                    S_01_exp(a = a,  t = t_type3[i, 1]))))))
  
  return(type_1_exp + type_2_exp + type_3_exp + type_4_exp) 
}

log_lik_semi_markov = rep(NA, nn)
log_lik_markov = rep(NA, nn)

for(mn in 1:nn){
  set.seed(mn)
  print(mn)
  eps = 0.001
  param = c(5, 3, 3, 4, 2, 4)
  ## For scenario 2
  #param = c(1.1, .9, 1, 1.1, 0.95, 1.1)
  param_markov = c(1,1,1)
  n = 1000
  t = simulation_t(param, n)
  
  t_type0 = t$t_type0
  t_type1 = t$t_type1
  t_type2 = t$t_type2
  t_type3 = t$t_type3
  
  print(nrow(t_type0) + nrow(t_type1) + nrow(t_type2) + nrow(t_type3))
  
  
  theta_hat_markov = nlminb(param_markov, sum_all_markov,t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, 
                            t_type3 = t_type3, mc_cores = 3, lower = c(eps, eps, eps))
  
  theta_hat = nlminb(param, sum_all,t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, t_type3 = t_type3,
                     mc_cores = 3, lower = c(eps, eps, eps, eps, eps, eps))
  
  theta_hat_middle_markov = as.matrix(theta_hat_markov$par)
  theta_hat_par_markov[mn, 1] = theta_hat_middle_markov[1, 1]
  theta_hat_par_markov[mn, 2] = theta_hat_middle_markov[2, 1]
  theta_hat_par_markov[mn, 3] = theta_hat_middle_markov[3, 1]
  
  theta_hat_hessian_markov = hessian(sum_all_markov, theta_hat_markov$par, t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, t_type3 = t_type3,
                                     mc_cores = 3)
  theta_hat_hessian_solved_markov = solve(theta_hat_hessian_markov)
  theta_hat_solve_hessian_markov[mn, 1] = theta_hat_hessian_solved_markov[1, 1]
  theta_hat_solve_hessian_markov[mn, 2] = theta_hat_hessian_solved_markov[2, 2]
  theta_hat_solve_hessian_markov[mn, 3] = theta_hat_hessian_solved_markov[3, 3]
  
  log_lik_markov[mn] = theta_hat_markov$objective
  
  theta_hat_middle = as.matrix(theta_hat$par)
  theta_hat_par[mn, 1] = theta_hat_middle[1, 1]
  theta_hat_par[mn, 2] = theta_hat_middle[2, 1]
  theta_hat_par[mn, 3] = theta_hat_middle[3, 1]
  theta_hat_par[mn, 4] = theta_hat_middle[4, 1]
  theta_hat_par[mn, 5] = theta_hat_middle[5, 1]
  theta_hat_par[mn, 6] = theta_hat_middle[6, 1]
  
  theta_hat_hessian = hessian(sum_all, theta_hat$par, t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, t_type3 = t_type3,mc_cores = 3)
  theta_hat_hessian_solved = solve(theta_hat_hessian)
  theta_hat_solve_hessian[mn, 1] = theta_hat_hessian_solved[1, 1]
  theta_hat_solve_hessian[mn, 2] = theta_hat_hessian_solved[2, 2]
  theta_hat_solve_hessian[mn, 3] = theta_hat_hessian_solved[3, 3]
  theta_hat_solve_hessian[mn, 4] = theta_hat_hessian_solved[4, 4]
  theta_hat_solve_hessian[mn, 5] = theta_hat_hessian_solved[5, 5]
  theta_hat_solve_hessian[mn, 6] = theta_hat_hessian_solved[6, 6]
  
  log_lik_semi_markov[mn] = theta_hat$objective
}
