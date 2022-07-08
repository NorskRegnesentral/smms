rm(list = ls())
library(Rlab)
library(numDeriv)
library(parallel)

## Simulations when the competing risks approach is true. 
nn = 100

theta_hat_par_embed = matrix(NA, nrow = nn, ncol = 10)
theta_hat_solve_hessian_embed = matrix(NA, nrow = nn, ncol = 10)

theta_hat_par_semi2_comp = matrix(NA, nrow = nn, ncol = 9)
theta_hat_solve_hessian_semi2_comp = matrix(NA, nrow = nn, ncol = 9)

## Initial survival function
S_01_embed = function(a, t, x){exp(a[3])*(1-pweibull(t,exp(a[1] + a[8]*x), exp(a[2])))}
S_12_embed = function(a, t, x){(1-pweibull(t,exp(a[4] + a[9]*x), exp(a[5])))}
S_02_embed = function(a, t, x){(1-exp(a[3]))*(1-pweibull(t, exp(a[6]+a[10]*x), exp(a[7])))}

## Initial density functions - derivative of -survival function
f_01_embed = function(a,  t,x){exp(a[3])*dweibull(t,exp(a[1] + a[8]*x), exp(a[2]))}
f_12_embed = function(a, t, x){dweibull(t,exp(a[4] + a[9]*x), exp(a[5]))}
f_02_embed = function(a, t, x){(1-exp(a[3]))*dweibull(t,exp(a[6] + a[10]*x), exp(a[7]))}

## Initial survival function
S_01_comp = function(a, t, x){(1-pweibull(t,exp(a[1] + a[7]*x), exp(a[2])))}
S_12_comp = function(a, t, x){(1-pweibull(t,exp(a[3] + a[8]*x), exp(a[4])))}
S_02_comp = function(a, t, x){(1-pweibull(t, exp(a[5] + a[9]*x), exp(a[6])))}

## Initial density functions - derivative of -survival function
f_01_comp = function(a, t, x){dweibull(t,exp(a[1] + a[7]*x), exp(a[2]))}
f_12_comp = function(a, t, x){dweibull(t,exp(a[3] + a[8]*x), exp(a[4]))}
f_02_comp = function(a, t,x){dweibull(t,exp(a[5] + a[9]*x), exp(a[6]))}

## Other functions
## Semi-Markov method 1:
f01_S12_S02_embed = function(a,  w, t, x){f_01_embed(a, t, x)*S_12_embed(a, w - t, x)}
f01_f12_S02_embed = function(a, w, t, x){f_01_embed(a, t, x)*f_12_embed(a, w - t, x)}
## Semi-Markov method 2:
f01_S12_S02_comp = function(a,  w, t, x){S_02_comp(a, t,x)*f_01_comp(a, t,x)*S_12_comp(a, w - t,x)}
f01_f12_S02_comp = function(a, w, t, x){S_02_comp(a, t,x)*f_01_comp(a, t,x)*f_12_comp(a, w - t,x)}

## Simulating time points
simulation_t = function(param, n){
  T_01 = rep(NA, n)
  T_12 = rep(NA, n)
  T_02 = rep(NA, n)
  
  t_type0 = matrix(, nrow = 0, ncol = 2)
  t_type1 = matrix(, nrow = 0, ncol = 4)
  t_type2 = matrix(, nrow = 0, ncol = 4)
  t_type3 = matrix(, nrow = 0, ncol = 2)
  
  for(i in 1:n){
    x = sample(0:1, 1)
    unif1 = runif(1, min = 0, max = 1)
    S01 = function(T)
    {
      1-pweibull(T, exp(param[1] + param[7]*x), exp(param[2])) - unif1
    }
    T_01[i] = uniroot(S01, interval = c(1.e-14, 1e04), tol = 1e-9)$root
    unif2 = runif(1, min = 0, max = 1)
    S12 = function(T)
    {
      1-pweibull(T, exp(param[3]+ param[8]*x), exp(param[4])) - unif2
    }
    T_12[i] = uniroot(S12, interval = c(1.e-14, 1e04), tol = 1e-9)$root
    unif3 = runif(1, min = 0, max = 1)
    S02 = function(T)
    {
      1-pweibull(T, exp(param[5] + param[9]*x), exp(param[6])) - unif3
    }
    T_02[i] = uniroot(S02, interval = c(1.e-14, 1e04), tol = 1e-9)$root
    
    mm <- sample(2:10,1)
    t <- c(0,sort(rlnorm(mm-1,-1,1)))
    t_n = length(t)
    if(t[t_n] < T_01[i] & t[t_n] < T_02[i]){
      t_type0 = rbind(t_type0, c(t[t_n], x))
    } else if(T_01[i] < T_02[i]){
      if(any(t > T_01[i]) & T_01[i] + T_12[i] > t[t_n]){
        lower_int = max(which(t < T_01[i]))
        upper_int = min(which(t > T_01[i]))
        t_type1 = rbind(t_type1, c(t[lower_int], t[upper_int], t[t_n], x))
      } else if(any(t > T_01[i] + T_12[i])){
        lower_int1 = max(which(t < T_01[i]))
        upper_int1 = min(which(t > T_01[i]))
        t_type2 = rbind(t_type2, c(t[lower_int1], min(T_01[i] + T_12[i],t[upper_int1]), 
                                   T_01[i] + T_12[i], x))
      }
    } else if(any(t > T_02[i]) & T_02[i] < T_01[i]){
      t_type3 = rbind(t_type3, c(T_02[i], x))
    }
  }
  return(list(t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, t_type3 = t_type3))
}

## Competing risks approach
sum_all_comp = function(a, t_type0, t_type1, t_type2, t_type3, mc_cores){
  type_1_comp = -sum(unlist(mclapply(1:nrow(t_type0), function(i) (log(S_01_comp(a = a,  t = t_type0[i, 1], x = t_type0[i, 2])*
                                                                          S_02_comp(a = a,  t = t_type0[i, 1], x = t_type0[i, 2]))))))
  type_2_comp = -sum(unlist(mclapply(1:nrow(t_type1), function(i)(log(as.numeric(
    integrate(f01_S12_S02_comp, lower = t_type1[i,1], upper = t_type1[i,2], w = t_type1[i,3], 
              a = a, x = t_type1[i, 4])[1]))),mc.cores = mc_cores), use.names = FALSE))
  type_3_comp = -sum(unlist(mclapply(1:nrow(t_type2), function(i)(log(as.numeric(
    integrate(f01_f12_S02_comp, lower = t_type2[i,1], upper = t_type2[i,2], 
              w = t_type2[i,3], a = a, x = t_type2[i, 4])[1]))),mc.cores = mc_cores), use.names = FALSE))
  type_4_comp = -sum(unlist(mclapply(1:nrow(t_type3), function(i) (log(f_02_comp(a = a,  t = t_type3[i, 1], x = t_type3[i, 2])*
                                                                          S_01_comp(a = a,  t = t_type3[i, 1], x = t_type3[i, 2]))))))
  
  return(type_1_comp + type_2_comp + type_3_comp + type_4_comp) 
}

## Embedded Markov chain approach
sum_all_embed = function(a, t_type0, t_type1, t_type2, t_type3, mc_cores){
  type_1_embed = -sum(unlist(mclapply(1:nrow(t_type0), function(i) (log(S_01_embed(a = a,  t = t_type0[i, 1], x = t_type0[i, 2])+
                                                                          S_02_embed(a = a,  t = t_type0[i, 1], x = t_type0[i, 2]))))))
  type_2_embed = -sum(unlist(mclapply(1:nrow(t_type1), function(i)(log(as.numeric(
    integrate(f01_S12_S02_embed, lower = t_type1[i,1], upper = t_type1[i,2], w = t_type1[i,3], 
              a = a, x = t_type1[i, 4])[1]))),mc.cores = mc_cores), use.names = FALSE))
  type_3_embed = -sum(unlist(mclapply(1:nrow(t_type2), function(i)(log(as.numeric(
    integrate(f01_f12_S02_embed, lower = t_type2[i,1], upper = t_type2[i,2], 
              w = t_type2[i,3], a = a, x = t_type2[i, 4])[1]))),mc.cores = mc_cores), use.names = FALSE))
  type_4_embed = -sum(unlist(mclapply(1:nrow(t_type3), function(i) (log(f_02_embed(a = a,  t = t_type3[i, 1], x = t_type3[i, 2]))))))
  
  return(type_1_embed + type_2_embed + type_3_embed + type_4_embed) 
}

log_lik_semi_markov = rep(NA, nn)
log_lik_markov = rep(NA, nn)

for(mn in 1:nn){
  set.seed(mn)
  print(mn)
  eps = 0.001
  param_comp = c(0.4,  0.6,  1.6, -0.5, 0.2,  1.1, -0.9, -1.2, 0)
  param_embed = c(0.4,  0.6,-1.2, 1.6, -0.5, 0.2,  1.1, -0.9, -1.2, 0)
  n = 1000
  t = simulation_t(param_comp, n)
  
  t_type0 = t$t_type0
  t_type1 = t$t_type1
  t_type2 = t$t_type2
  t_type3 = t$t_type3
  
  print(nrow(t_type0) + nrow(t_type1) + nrow(t_type2) + nrow(t_type3))
  
  theta_hat_comp = nlminb(param_comp, sum_all_comp,t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, t_type3 = t_type3, mc_cores = 3)
  
  
  theta_hat_embed = nlminb(param_embed, sum_all_embed,t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2,
                     t_type3 = t_type3, mc_cores = 3)
  
  theta_hat_middle_comp = as.matrix(theta_hat_comp$par)
  theta_hat_par_comp[mn, 1] = theta_hat_middle_comp[1, 1]
  theta_hat_par_comp[mn, 2] = theta_hat_middle_comp[2, 1]
  theta_hat_par_comp[mn, 3] = theta_hat_middle_comp[3, 1]
  theta_hat_par_comp[mn, 4] = theta_hat_middle_comp[4, 1]
  theta_hat_par_comp[mn, 5] = theta_hat_middle_comp[5, 1]
  theta_hat_par_comp[mn, 6] = theta_hat_middle_comp[6, 1]
  theta_hat_par_comp[mn, 7] = theta_hat_middle_comp[7, 1]
  theta_hat_par_comp[mn, 8] = theta_hat_middle_comp[8, 1]
  theta_hat_par_comp[mn, 9] = theta_hat_middle_comp[9, 1]
  
  theta_hat_hessian_comp = hessian(sum_all_comp, theta_hat_comp$par, t_type0 = t_type0, t_type1 = t_type1, 
                                    t_type2 = t_type2, t_type3 = t_type3,
                                    mc_cores = 3)
  theta_hat_hessian_solved_comp = solve(theta_hat_hessian_comp)
  theta_hat_solve_hessian_comp[mn, 1] = theta_hat_hessian_solved_comp[1, 1]
  theta_hat_solve_hessian_comp[mn, 2] = theta_hat_hessian_solved_comp[2, 2]
  theta_hat_solve_hessian_comp[mn, 3] = theta_hat_hessian_solved_comp[3, 3]
  theta_hat_solve_hessian_comp[mn, 4] = theta_hat_hessian_solved_comp[4, 4]
  theta_hat_solve_hessian_comp[mn, 5] = theta_hat_hessian_solved_comp[5, 5]
  theta_hat_solve_hessian_comp[mn, 6] = theta_hat_hessian_solved_comp[6, 6]
  theta_hat_solve_hessian_comp[mn, 7] = theta_hat_hessian_solved_comp[7, 7]
  theta_hat_solve_hessian_comp[mn, 8] = theta_hat_hessian_solved_comp[8, 8]
  theta_hat_solve_hessian_comp[mn, 9] = theta_hat_hessian_solved_comp[9, 9]
  
  theta_hat_middle = as.matrix(theta_hat_embed$par)
  theta_hat_par_embed[mn, 1] = theta_hat_middle[1, 1]
  theta_hat_par_embed[mn, 2] = theta_hat_middle[2, 1]
  theta_hat_par_embed[mn, 3] = theta_hat_middle[3, 1]
  theta_hat_par_embed[mn, 4] = theta_hat_middle[4, 1]
  theta_hat_par_embed[mn, 5] = theta_hat_middle[5, 1]
  theta_hat_par_embed[mn, 6] = theta_hat_middle[6, 1]
  theta_hat_par_embed[mn, 7] = theta_hat_middle[7, 1]
  theta_hat_par_embed[mn, 8] = theta_hat_middle[8, 1]
  theta_hat_par_embed[mn, 9] = theta_hat_middle[9, 1]
  theta_hat_par_embed[mn, 10] = theta_hat_middle[10, 1]
  
  
  theta_hat_hessian_embed = hessian(sum_all_embed, theta_hat_embed$par, t_type0 = t_type0, t_type1 = t_type1, t_type2 = t_type2, t_type3 = t_type3, mc_cores = 3)
  theta_hat_hessian_solved_embed = solve(theta_hat_hessian_embed)
  theta_hat_solve_hessian_embed[mn, 1] = theta_hat_hessian_solved_embed[1, 1]
  theta_hat_solve_hessian_embed[mn, 2] = theta_hat_hessian_solved_embed[2, 2]
  theta_hat_solve_hessian_embed[mn, 3] = theta_hat_hessian_solved_embed[3, 3]
  theta_hat_solve_hessian_embed[mn, 4] = theta_hat_hessian_solved_embed[4, 4]
  theta_hat_solve_hessian_embed[mn, 5] = theta_hat_hessian_solved_embed[5, 5]
  theta_hat_solve_hessian_embed[mn, 6] = theta_hat_hessian_solved_embed[6, 6]
  theta_hat_solve_hessian_embed[mn, 7] = theta_hat_hessian_solved_embed[7, 7]
  theta_hat_solve_hessian_embed[mn, 8] = theta_hat_hessian_solved_embed[8, 8]
  theta_hat_solve_hessian_embed[mn, 9] = theta_hat_hessian_solved_embed[9, 9]
  theta_hat_solve_hessian_embed[mn, 10] = theta_hat_hessian_solved_embed[10, 10]
  
}
