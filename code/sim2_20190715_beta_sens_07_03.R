## Given some level of noise, change the beta coefficient, see why this is affecting things


# 0.7/0.3 -----------------------------------------------------------------

p_n_ind = 70
p_n_pheno = 128
p_tying_grp = c(1,2,3,65)
p_param_pi_00 = 0
p_param_pi_11 = 1
p_grped_eps_11 = c(rep(0.7,(length(p_tying_grp)-1)),0)
p_grped_eps_10 = c(rep(0.3,(length(p_tying_grp)-1)),0)
p_coeff_mu = 0
# p_coeff_beta = 0.1
p_param_gi_prob = 0.4
p_param_sigma_beta = 5 # 0.5
p_num_sims = 50
p_seed = 20

sim_results <- list()

betas_to_try <- seq(0.1,2,0.1)

n = 1
for(beta_coeff in betas_to_try){
  sim_results[[n]] <- run_sim2(
    n_ind = p_n_ind
    ,n_pheno = p_n_pheno
    ,tying_grp = p_tying_grp
    ,param_pi_00 = p_param_pi_00
    ,param_pi_11 = p_param_pi_11
    ,grped_eps_11 = p_grped_eps_11
    ,grped_eps_10 = p_grped_eps_10
    ,coeff_mu = p_coeff_mu
    ,coeff_beta = beta_coeff
    ,param_gi_prob = p_param_gi_prob
    ,param_sigma_beta = p_param_sigma_beta
    ,num_sims = p_num_sims
    ,seed = p_seed)
  n = n + 1
}

pi_mean_list <- list()
pi_sd_list <- list()
eps_11_mean_list <- list()
eps_11_sd_list <- list()
eps_10_mean_list <- list()
eps_10_sd_list <- list()
for(i in 1:length(sim_results)){
  pi_mean_list[[i]] <- mean(unlist(sim_results[[i]]$results_pi))
  pi_sd_list[[i]] <- sd(unlist(sim_results[[i]]$results_pi))

  results_eps_11_mtx <- matrix(unlist(sim_results[[i]]$results_eps_11)
                               ,nrow = p_num_sims, ncol = (p_n_pheno - 2), byrow = T)
  results_eps_10_mtx <- matrix(unlist(sim_results[[i]]$results_eps_10)
                               ,nrow = p_num_sims, ncol = (p_n_pheno - 2), byrow = T)

  eps_11_mean_list[[i]] <- apply(results_eps_11_mtx,MARGIN = 2,mean)
  eps_11_sd_list[[i]] <- apply(results_eps_11_mtx,MARGIN = 2,sd)
  eps_10_mean_list[[i]] <- apply(results_eps_10_mtx,MARGIN = 2,mean)
  eps_10_sd_list[[i]] <- apply(results_eps_10_mtx,MARGIN = 2,sd)
}

# Pi
pi_mean_results <- unlist(pi_mean_list)
pi_sd_results <- unlist(pi_sd_list)
y_axis_bounds <- c(min(pi_mean_results - 2*pi_sd_results),max(pi_mean_results + 2*pi_sd_results))
plot(x = betas_to_try
     , y = rep(p_param_pi_11, length(betas_to_try))
     , col = "black"
     , type = "l"
     , xlab = "Beta coeffs"
     , ylab = "Pi estimate"
     , main = "Sensitity of Pi estimate by beta coeff (beta noise = 5, eps_11 = 0.7, eps_10 = 0.3)"
     , xaxt = "n"
     # , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = pi_mean_results,col = "red", type = "l")
points(x = betas_to_try, y = pi_mean_results + pi_sd_results*2, col = "blue", type = "l")
points(x = betas_to_try, y = pi_mean_results - pi_sd_results*2, col = "blue", type = "l")

# Eps_11
eps_11_means <- sapply(eps_11_mean_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_11_sd <- sapply(eps_11_sd_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_11_lb <- eps_11_means - 2*eps_11_sd
eps_11_ub <- eps_11_means + 2*eps_11_sd

y_axis_bounds <- c(min(min(eps_11_lb[3,]),min(p_grped_eps_11[3]))
                   ,max(max(eps_11_ub[3,]),max(p_grped_eps_11[3])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_11[3],eps_11_lb[3,],eps_11_ub[3,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "Eps_11 estimate", main = "Eps_11 estimation sensitivity to beta - grp 2 (beta noise = 5, eps_11 = 0.7, eps_10 = 0.3)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_11[3],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_11_means[3,],col = "red", type = "l")
points(x = betas_to_try, y = eps_11_lb[3,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_11_ub[3,],col = "blue", type = "l", lty = 2)

y_axis_bounds <- c(min(min(eps_11_lb[4,]),min(p_grped_eps_11[4]))
                   ,max(max(eps_11_ub[4,]),max(p_grped_eps_11[4])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_11[4],eps_11_lb[4,],eps_11_ub[4,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "Eps_11 estimate", main = "Eps_11 estimation sensitivity to beta - grp 3 (beta noise = 5, eps_11 = 0.7, eps_10 = 0.3)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_11[4],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_11_means[4,],col = "red", type = "l")
points(x = betas_to_try, y = eps_11_lb[4,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_11_ub[4,],col = "blue", type = "l", lty = 2)

# Eps_10
eps_10_means <- sapply(eps_10_mean_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_10_sd <- sapply(eps_10_sd_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_10_lb <- eps_10_means - 2*eps_10_sd
eps_10_ub <- eps_10_means + 2*eps_10_sd

y_axis_bounds <- c(min(min(eps_10_lb[3,]),min(p_grped_eps_10[3]))
                   ,max(max(eps_10_ub[3,]),max(p_grped_eps_10[3])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_10[3],eps_10_lb[3,],eps_10_ub[3,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "eps_10 estimate", main = "eps_10 estimation sensitivity to beta - grp 2 (beta noise = 5, eps_11 = 0.7, eps_10 = 0.3)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_10[3],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_10_means[3,],col = "red", type = "l")
points(x = betas_to_try, y = eps_10_lb[3,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_10_ub[3,],col = "blue", type = "l", lty = 2)

y_axis_bounds <- c(min(min(eps_10_lb[4,]),min(p_grped_eps_10[4]))
                   ,max(max(eps_10_ub[4,]),max(p_grped_eps_10[4])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_10[4],eps_10_lb[4,],eps_10_ub[4,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "eps_10 estimate", main = "eps_10 estimation sensitivity to beta - grp 3 (beta noise = 5, eps_11 = 0.7, eps_10 = 0.3)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_10[4],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_10_means[4,],col = "red", type = "l")
points(x = betas_to_try, y = eps_10_lb[4,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_10_ub[4,],col = "blue", type = "l", lty = 2)


write.table(sim_results[[10]]$results_y_mtx[[25]], file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/beta_sens_07_03_b1.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
cat(c("sim2","A","A",sim_results[[10]]$results_g_seq[[25]]), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/beta_sens_07_03_b1.cis.geno")

write.table(sim_results[[20]]$results_y_mtx[[25]], file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/beta_sens_07_03_b2.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
cat(c("sim2","A","A",sim_results[[20]]$results_g_seq[[25]]), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/beta_sens_07_03_b2.cis.geno")

# 0.5/0.5 -----------------------------------------------------------------

p_n_ind = 70
p_n_pheno = 128
p_tying_grp = c(1,2,3,65)
p_param_pi_00 = 0
p_param_pi_11 = 1
p_grped_eps_11 = c(rep(0.5,(length(p_tying_grp)-1)),0)
p_grped_eps_10 = c(rep(0.5,(length(p_tying_grp)-1)),0)
p_coeff_mu = 0
# p_coeff_beta = 0.1
p_param_gi_prob = 0.4
p_param_sigma_beta = 5 # 0.5
p_num_sims = 50
p_seed = 20

sim_results <- list()

betas_to_try <- seq(0.1,2,0.1)

n = 1
for(beta_coeff in betas_to_try){
  sim_results[[n]] <- run_sim2(
    n_ind = p_n_ind
    ,n_pheno = p_n_pheno
    ,tying_grp = p_tying_grp
    ,param_pi_00 = p_param_pi_00
    ,param_pi_11 = p_param_pi_11
    ,grped_eps_11 = p_grped_eps_11
    ,grped_eps_10 = p_grped_eps_10
    ,coeff_mu = p_coeff_mu
    ,coeff_beta = beta_coeff
    ,param_gi_prob = p_param_gi_prob
    ,param_sigma_beta = p_param_sigma_beta
    ,num_sims = p_num_sims
    ,seed = p_seed)
  n = n + 1
}

pi_mean_list <- list()
pi_sd_list <- list()
eps_11_mean_list <- list()
eps_11_sd_list <- list()
eps_10_mean_list <- list()
eps_10_sd_list <- list()
for(i in 1:length(sim_results)){
  pi_mean_list[[i]] <- mean(unlist(sim_results[[i]]$results_pi))
  pi_sd_list[[i]] <- sd(unlist(sim_results[[i]]$results_pi))

  results_eps_11_mtx <- matrix(unlist(sim_results[[i]]$results_eps_11)
                               ,nrow = p_num_sims, ncol = (p_n_pheno - 2), byrow = T)
  results_eps_10_mtx <- matrix(unlist(sim_results[[i]]$results_eps_10)
                               ,nrow = p_num_sims, ncol = (p_n_pheno - 2), byrow = T)

  eps_11_mean_list[[i]] <- apply(results_eps_11_mtx,MARGIN = 2,mean)
  eps_11_sd_list[[i]] <- apply(results_eps_11_mtx,MARGIN = 2,sd)
  eps_10_mean_list[[i]] <- apply(results_eps_10_mtx,MARGIN = 2,mean)
  eps_10_sd_list[[i]] <- apply(results_eps_10_mtx,MARGIN = 2,sd)
}

# Pi
pi_mean_results <- unlist(pi_mean_list)
pi_sd_results <- unlist(pi_sd_list)
y_axis_bounds <- c(min(pi_mean_results - 2*pi_sd_results),max(pi_mean_results + 2*pi_sd_results))
plot(x = betas_to_try
     , y = rep(p_param_pi_11, length(betas_to_try))
     , col = "black"
     , type = "l"
     , xlab = "Beta coeffs"
     , ylab = "Pi estimate"
     , main = "Sensitity of Pi estimate by beta coeff (beta noise = 5, eps_11 = 0.5, eps_10 = 0.5)"
     , xaxt = "n"
     # , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = pi_mean_results,col = "red", type = "l")
points(x = betas_to_try, y = pi_mean_results + pi_sd_results*2, col = "blue", type = "l")
points(x = betas_to_try, y = pi_mean_results - pi_sd_results*2, col = "blue", type = "l")

# Eps_11
eps_11_means <- sapply(eps_11_mean_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_11_sd <- sapply(eps_11_sd_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_11_lb <- eps_11_means - 2*eps_11_sd
eps_11_ub <- eps_11_means + 2*eps_11_sd

y_axis_bounds <- c(min(min(eps_11_lb[3,]),min(p_grped_eps_11[3]))
                   ,max(max(eps_11_ub[3,]),max(p_grped_eps_11[3])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_11[3],eps_11_lb[3,],eps_11_ub[3,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "Eps_11 estimate", main = "Eps_11 estimation sensitivity to beta - grp 2 (beta noise = 5, eps_11 = 0.5, eps_10 = 0.5)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_11[3],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_11_means[3,],col = "red", type = "l")
points(x = betas_to_try, y = eps_11_lb[3,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_11_ub[3,],col = "blue", type = "l", lty = 2)

y_axis_bounds <- c(min(min(eps_11_lb[4,]),min(p_grped_eps_11[4]))
                   ,max(max(eps_11_ub[4,]),max(p_grped_eps_11[4])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_11[4],eps_11_lb[4,],eps_11_ub[4,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "Eps_11 estimate", main = "Eps_11 estimation sensitivity to beta - grp 3 (beta noise = 5, eps_11 = 0.5 eps_10 = 0.5)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_11[4],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_11_means[4,],col = "red", type = "l")
points(x = betas_to_try, y = eps_11_lb[4,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_11_ub[4,],col = "blue", type = "l", lty = 2)

# Eps_10
eps_10_means <- sapply(eps_10_mean_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_10_sd <- sapply(eps_10_sd_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_10_lb <- eps_10_means - 2*eps_10_sd
eps_10_ub <- eps_10_means + 2*eps_10_sd

y_axis_bounds <- c(min(min(eps_10_lb[3,]),min(p_grped_eps_10[3]))
                   ,max(max(eps_10_ub[3,]),max(p_grped_eps_10[3])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_10[3],eps_10_lb[3,],eps_10_ub[3,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "eps_10 estimate", main = "eps_10 estimation sensitivity to beta - grp 2 (beta noise = 5, eps_11 = 0.5 eps_10 = 0.5)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_10[3],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_10_means[3,],col = "red", type = "l")
points(x = betas_to_try, y = eps_10_lb[3,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_10_ub[3,],col = "blue", type = "l", lty = 2)

y_axis_bounds <- c(min(min(eps_10_lb[4,]),min(p_grped_eps_10[4]))
                   ,max(max(eps_10_ub[4,]),max(p_grped_eps_10[4])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_10[4],eps_10_lb[4,],eps_10_ub[4,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "eps_10 estimate", main = "eps_10 estimation sensitivity to beta - grp 3 (beta noise = 5, eps_11 = 0.5 eps_10 = 0.5)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_10[4],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_10_means[4,],col = "red", type = "l")
points(x = betas_to_try, y = eps_10_lb[4,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_10_ub[4,],col = "blue", type = "l", lty = 2)

# 0.75/0.75 -----------------------------------------------------------------

p_n_ind = 70
p_n_pheno = 128
p_tying_grp = c(1,2,3,65)
p_param_pi_00 = 0
p_param_pi_11 = 0.5
p_grped_eps_11 = c(rep(0.75,(length(p_tying_grp)-1)),0)
p_grped_eps_10 = c(rep(0.75,(length(p_tying_grp)-1)),0)
p_coeff_mu = 0
# p_coeff_beta = 0.1
p_param_gi_prob = 0.4
p_param_sigma_beta = 5 # 0.5
p_num_sims = 50
p_seed = 20

sim_results <- list()

betas_to_try <- seq(0.1,2,0.1)

n = 1
for(beta_coeff in betas_to_try){
  sim_results[[n]] <- run_sim2(
    n_ind = p_n_ind
    ,n_pheno = p_n_pheno
    ,tying_grp = p_tying_grp
    ,param_pi_00 = p_param_pi_00
    ,param_pi_11 = p_param_pi_11
    ,grped_eps_11 = p_grped_eps_11
    ,grped_eps_10 = p_grped_eps_10
    ,coeff_mu = p_coeff_mu
    ,coeff_beta = beta_coeff
    ,param_gi_prob = p_param_gi_prob
    ,param_sigma_beta = p_param_sigma_beta
    ,num_sims = p_num_sims
    ,seed = p_seed)
  n = n + 1
}

pi_mean_list <- list()
pi_sd_list <- list()
eps_11_mean_list <- list()
eps_11_sd_list <- list()
eps_10_mean_list <- list()
eps_10_sd_list <- list()
for(i in 1:length(sim_results)){
  pi_mean_list[[i]] <- mean(unlist(sim_results[[i]]$results_pi))
  pi_sd_list[[i]] <- sd(unlist(sim_results[[i]]$results_pi))

  results_eps_11_mtx <- matrix(unlist(sim_results[[i]]$results_eps_11)
                               ,nrow = p_num_sims, ncol = (p_n_pheno - 2), byrow = T)
  results_eps_10_mtx <- matrix(unlist(sim_results[[i]]$results_eps_10)
                               ,nrow = p_num_sims, ncol = (p_n_pheno - 2), byrow = T)

  eps_11_mean_list[[i]] <- apply(results_eps_11_mtx,MARGIN = 2,mean)
  eps_11_sd_list[[i]] <- apply(results_eps_11_mtx,MARGIN = 2,sd)
  eps_10_mean_list[[i]] <- apply(results_eps_10_mtx,MARGIN = 2,mean)
  eps_10_sd_list[[i]] <- apply(results_eps_10_mtx,MARGIN = 2,sd)
}

# Pi
pi_mean_results <- unlist(pi_mean_list)
pi_sd_results <- unlist(pi_sd_list)
y_axis_bounds <- c(min(pi_mean_results - 2*pi_sd_results),max(pi_mean_results + 2*pi_sd_results))
plot(x = betas_to_try
     , y = rep(p_param_pi_11, length(betas_to_try))
     , col = "black"
     , type = "l"
     , xlab = "Beta coeffs"
     , ylab = "Pi estimate"
     , main = "Sensitity of Pi estimate by beta coeff (beta noise = 5, eps_11 = 0.5, eps_10 = 0.5)"
     , xaxt = "n"
     # , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = pi_mean_results,col = "red", type = "l")
points(x = betas_to_try, y = pi_mean_results + pi_sd_results*2, col = "blue", type = "l")
points(x = betas_to_try, y = pi_mean_results - pi_sd_results*2, col = "blue", type = "l")

# Eps_11
eps_11_means <- sapply(eps_11_mean_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_11_sd <- sapply(eps_11_sd_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_11_lb <- eps_11_means - 2*eps_11_sd
eps_11_ub <- eps_11_means + 2*eps_11_sd

y_axis_bounds <- c(min(min(eps_11_lb[3,]),min(p_grped_eps_11[3]))
                   ,max(max(eps_11_ub[3,]),max(p_grped_eps_11[3])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_11[3],eps_11_lb[3,],eps_11_ub[3,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "Eps_11 estimate", main = "Eps_11 estimation sensitivity to beta - grp 2 (beta noise = 5, eps_11 = 0.5, eps_10 = 0.5)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_11[3],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_11_means[3,],col = "red", type = "l")
points(x = betas_to_try, y = eps_11_lb[3,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_11_ub[3,],col = "blue", type = "l", lty = 2)

y_axis_bounds <- c(min(min(eps_11_lb[4,]),min(p_grped_eps_11[4]))
                   ,max(max(eps_11_ub[4,]),max(p_grped_eps_11[4])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_11[4],eps_11_lb[4,],eps_11_ub[4,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "Eps_11 estimate", main = "Eps_11 estimation sensitivity to beta - grp 3 (beta noise = 5, eps_11 = 0.5 eps_10 = 0.5)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_11[4],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_11_means[4,],col = "red", type = "l")
points(x = betas_to_try, y = eps_11_lb[4,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_11_ub[4,],col = "blue", type = "l", lty = 2)

# Eps_10
eps_10_means <- sapply(eps_10_mean_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_10_sd <- sapply(eps_10_sd_list, function(x){
  y <- c(0,0,x)[p_tying_grp]
  return(y)
})
eps_10_lb <- eps_10_means - 2*eps_10_sd
eps_10_ub <- eps_10_means + 2*eps_10_sd

y_axis_bounds <- c(min(min(eps_10_lb[3,]),min(p_grped_eps_10[3]))
                   ,max(max(eps_10_ub[3,]),max(p_grped_eps_10[3])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_10[3],eps_10_lb[3,],eps_10_ub[3,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "eps_10 estimate", main = "eps_10 estimation sensitivity to beta - grp 2 (beta noise = 5, eps_11 = 0.5 eps_10 = 0.5)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_10[3],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_10_means[3,],col = "red", type = "l")
points(x = betas_to_try, y = eps_10_lb[3,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_10_ub[3,],col = "blue", type = "l", lty = 2)

y_axis_bounds <- c(min(min(eps_10_lb[4,]),min(p_grped_eps_10[4]))
                   ,max(max(eps_10_ub[4,]),max(p_grped_eps_10[4])))*1.1
xval = betas_to_try
which_in_bound = xval[between(p_grped_eps_10[4],eps_10_lb[4,],eps_10_ub[4,])]
plot(1,1, type = "n"
     , xlab = "Beta coeff", ylab = "eps_10 estimate", main = "eps_10 estimation sensitivity to beta - grp 3 (beta noise = 5, eps_11 = 0.5 eps_10 = 0.5)"
     , xaxt = "n"
     , xlim = c(min(betas_to_try),max(betas_to_try))
     , ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j] - 0.05, which_in_bound[j] - 0.05, which_in_bound[j] + 0.05, which_in_bound[j] + 0.05), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}
axis(1,at = betas_to_try,labels = betas_to_try,las = 2)
points(x = betas_to_try, y = rep(p_grped_eps_10[4],length(betas_to_try)),type = "l")
points(x = betas_to_try, y = eps_10_means[4,],col = "red", type = "l")
points(x = betas_to_try, y = eps_10_lb[4,],col = "blue", type = "l", lty = 2)
points(x = betas_to_try, y = eps_10_ub[4,],col = "blue", type = "l", lty = 2)
