# Generate data for EM algorithm deep dive.

# Find some good and bad simulations
  # Good: Epsilon vaguely reflects simulated proportions of transitions
  # Bad: Epsilon does not reflect

# Some pretty heavy tying in place too, for now.

## 0.7/0.3

source("code/sim2_script.R")
library(data.table)

# Params
p_n_ind = 70
p_n_pheno = 128
p_tying_grp = c(1,2,3,65)
p_param_pi_00 = 0
p_param_pi_11 = 1
p_grped_eps_11 = c(rep(0.7,(length(p_tying_grp)-1)),0)
p_grped_eps_10 = c(rep(0.3,(length(p_tying_grp)-1)),0)
p_coeff_mu = 0
p_coeff_beta = 2
p_param_gi_prob = 0.4
p_param_sigma_beta = 0.5
p_num_sims = 100
p_seed = 20

results_list <- run_sim2(
  n_ind = p_n_ind
  ,n_pheno = p_n_pheno
  ,tying_grp = p_tying_grp
  ,param_pi_00 = p_param_pi_00
  ,param_pi_11 = p_param_pi_11
  ,grped_eps_11 = p_grped_eps_11
  ,grped_eps_10 = p_grped_eps_10
  ,coeff_mu = p_coeff_mu
  ,coeff_beta = p_coeff_beta
  ,param_gi_prob = p_param_gi_prob
  ,param_sigma_beta = p_param_sigma_beta
  ,num_sims = p_num_sims
  ,seed = p_seed)

# Verify the 'badness' of this situation

results_eps_11_mtx <- matrix(unlist(results_list$results_eps_11)
                             ,nrow = p_num_sims, ncol = (p_n_pheno - 2), byrow = T)
results_eps_10_mtx <- matrix(unlist(results_list$results_eps_10)
                             ,nrow = p_num_sims, ncol = (p_n_pheno - 2), byrow = T)

eps_11_mean <- apply(results_eps_11_mtx,MARGIN = 2,mean)
eps_11_sd <- apply(results_eps_11_mtx,MARGIN = 2,sd)
eps_11_lb <- eps_11_mean - 3*eps_11_sd
eps_11_ub <- eps_11_mean + 3*eps_11_sd
eps_11_within_range <- between(results_list$param_eps_11,eps_11_lb,eps_11_ub)
table(eps_11_within_range)

eps_10_mean <- apply(results_eps_10_mtx,MARGIN = 2,mean)
eps_10_sd <- apply(results_eps_10_mtx,MARGIN = 2,sd)
eps_10_lb <- eps_10_mean - 3*eps_10_sd
eps_10_ub <- eps_10_mean + 3*eps_10_sd
eps_10_within_range <- between(results_list$param_eps_10,eps_10_lb,eps_10_ub)
table(eps_10_within_range)

# Epsilon 11
y_axis_bounds <- c(min(min(eps_11_lb),min(results_list$param_eps_11)),max(max(eps_11_ub),max(results_list$param_eps_11)))*1.1

xval = 1:p_n_pheno
which_in_bound = xval[between(results_list$param_eps_11,eps_11_lb,eps_11_ub)]

plot(1,1, type = "n", xlab = "Wavelet scale-loc", ylab = "Probability", main = "Epsilon 11", xaxt = "n"
     ,xlim = c(0,p_n_pheno)
     ,ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j]-0.5, which_in_bound[j]-0.5, which_in_bound[j]+0.5, which_in_bound[j]+0.5), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}

axis(1,at = 2^(3:10),labels = 2^(3:10),las = 2)
points(results_list$param_eps_11,type = "l")
points(eps_11_mean,col = "red", type = "l")
points(eps_11_ub,col = "blue", type = "l", lty = 2)
points(eps_11_lb,col = "blue", type = "l", lty = 2)

# Epsilon 10
y_axis_bounds <- c(min(min(eps_10_lb),min(results_list$param_eps_10)),max(max(eps_10_ub),max(results_list$param_eps_10)))*1.1

xval = 1:p_n_pheno
which_in_bound = xval[between(results_list$param_eps_10,eps_10_lb,eps_10_ub)]

plot(1,1, type = "n", xlab = "Wavelet scale-loc", ylab = "Probability", main = "Epsilon 10", xaxt = "n"
     ,xlim = c(0,p_n_pheno)
     ,ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j]-0.5, which_in_bound[j]-0.5, which_in_bound[j]+0.5, which_in_bound[j]+0.5), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}

axis(1,at = 2^(3:10),labels = 2^(3:10),las = 2)
points(results_list$param_eps_10,type = "l")
points(eps_10_mean,col = "red", type = "l")
points(eps_10_ub,col = "blue", type = "l", lty = 2)
points(eps_10_lb,col = "blue", type = "l", lty = 2)

# Try and reconcile simulated vs retrieved transition proportions
get_gamma_transition_props <- function(eps_no_scale_data, tying_grp_vect, n_pheno){
  gamma_and_parent <- matrix(c(eps_no_scale_data[get_parent_indices(1:(n_pheno-1))],eps_no_scale_data)
                             , nrow = 2, ncol = (n_pheno-1)
                             , byrow = T)
  gamma_and_parent_dt <- as.data.table(t(gamma_and_parent))
  setnames(gamma_and_parent_dt,names(gamma_and_parent_dt),c("Parent","Child"))
  gamma_and_parent_dt[,"Transition" := paste0(Child,Parent)]
  gamma_and_parent_dt[,"TreeLvl" := findInterval(.I + 1,tying_grp_vect[-1])]

  # Exclude the root
  gamma_and_parent_dt <- copy(gamma_and_parent_dt[-1])
  gamma_and_parent_stats <- dcast.data.table(
    gamma_and_parent_dt[,.N,by = .(Transition,TreeLvl)]
    ,formula = TreeLvl ~ Transition
    ,value.var = "N"
    ,fill = 0
  )

  transitions_not_there <- setdiff(c("11","01","10","00"),names(gamma_and_parent_stats))
  if(length(transitions_not_there) > 0){
    for(col in transitions_not_there){
      gamma_and_parent_stats[,(col) := 0]
    }
  }

  gamma_and_parent_stats[,"Total" := apply(.SD,1,sum),.SDcols = 2:5]
  return(gamma_and_parent_stats[])
}

# Now cycle through each of the simulations and find some 'good' and 'bad' situations
eps_11_lists <- list()
eps_10_lists <- list()
gamma_seq_lists <- list()
gamma_transition_lists <- list()
for(i in 1:p_num_sims){
  eps_11_lists[[i]] <- c(0,0,results_list$results_eps_11[[i]])[p_tying_grp]
  eps_10_lists[[i]] <- c(0,0,results_list$results_eps_10[[i]])[p_tying_grp]
  gamma_seq_lists[[i]] <- results_list$results_gamma_seq[[i]][-1]
  gamma_transition_lists[[i]] <- get_gamma_transition_props(gamma_seq_lists[[i]]
                                                            ,tying_grp_vect = p_tying_grp
                                                            ,n_pheno = p_n_pheno)
  # Add simulated data
  gamma_transition_lists[[i]][,c("sim_eps_11","sim_eps_10") :=
                          .(`11`/(`11`+`01`)
                            ,`10`/(`10`+`00`))]
  gamma_transition_lists[[i]][,c("param_eps_11","param_eps_10") :=
                          .(p_grped_eps_11[-(1:2)]
                            ,p_grped_eps_10[-(1:2)])]
  gamma_transition_lists[[i]][,c("algo_eps_11","algo_eps_10") :=
                          .(round(eps_11_lists[[i]][-(1:2)],2)
                            ,round(eps_10_lists[[i]][-(1:2)],2))]
}

distanceList <- lapply(gamma_transition_lists,function(x){
  x[,sqrt(sum((sim_eps_11-algo_eps_11)^2+(sim_eps_10-algo_eps_10)^2))]
})
best_3 <- order(unlist(distanceList))[1:3]
worst_3 <- order(-unlist(distanceList))[1:3]

gamma_transition_lists[best_3]
gamma_transition_lists[worst_3]

# But these just seem to be 'least worst' because the simulations are closer to 0.7 than others
# Can we find any where the algo has predicted somewhere around the 0.7, and 0.3 that we want?

distanceList <- lapply(gamma_transition_lists,function(x){
  x[,sqrt(sum((param_eps_11-algo_eps_11)^2+(param_eps_10-algo_eps_10)^2))]
})
best_3 <- order(unlist(distanceList))[1:3]
worst_3 <- order(-unlist(distanceList))[1:3]

gamma_transition_lists[best_3]
gamma_transition_lists[worst_3]

# No, it seems.

# Retain one good, and one bad case here:

sim_07_03_good <- list(results_pi = results_list$results_pi[[best_3[3]]]
                       ,results_eps_11 = results_list$results_eps_11[[best_3[3]]]
                       ,results_eps_10 = results_list$results_eps_10[[best_3[3]]]
                       ,results_gamma_seq = results_list$results_gamma_seq[[best_3[3]]]
                       ,results_g_seq = results_list$results_g_seq[[best_3[3]]]
                       ,results_y_mtx = results_list$results_y_mtx[[best_3[3]]])

sim_07_03_bad <- list(results_pi = results_list$results_pi[[worst_3[2]]]
                      ,results_eps_11 = results_list$results_eps_11[[worst_3[2]]]
                      ,results_eps_10 = results_list$results_eps_10[[worst_3[2]]]
                      ,results_gamma_seq = results_list$results_gamma_seq[[worst_3[2]]]
                      ,results_g_seq = results_list$results_g_seq[[worst_3[2]]]
                      ,results_y_mtx = results_list$results_y_mtx[[worst_3[2]]])

results_list_07_03 <- copy(results_list)

### Try again with different parameters: ----
# 0.9/0.1
# Do we get more variation?

# Params
p_n_ind = 70
p_n_pheno = 128
p_tying_grp = c(1,2,3,65)
p_param_pi_00 = 0
p_param_pi_11 = 1
p_grped_eps_11 = c(rep(0.9,(length(p_tying_grp)-1)),0)
p_grped_eps_10 = c(rep(0.1,(length(p_tying_grp)-1)),0)
p_coeff_mu = 0
p_coeff_beta = 2
p_param_gi_prob = 0.4
p_param_sigma_beta = 0.5
p_num_sims = 100
p_seed = 20

results_list <- run_sim2(
  n_ind = p_n_ind
  ,n_pheno = p_n_pheno
  ,tying_grp = p_tying_grp
  ,param_pi_00 = p_param_pi_00
  ,param_pi_11 = p_param_pi_11
  ,grped_eps_11 = p_grped_eps_11
  ,grped_eps_10 = p_grped_eps_10
  ,coeff_mu = p_coeff_mu
  ,coeff_beta = p_coeff_beta
  ,param_gi_prob = p_param_gi_prob
  ,param_sigma_beta = p_param_sigma_beta
  ,num_sims = p_num_sims
  ,seed = p_seed)

# Verify the 'badness' of this situation

results_eps_11_mtx <- matrix(unlist(results_list$results_eps_11)
                             ,nrow = p_num_sims, ncol = (p_n_pheno - 2), byrow = T)
results_eps_10_mtx <- matrix(unlist(results_list$results_eps_10)
                             ,nrow = p_num_sims, ncol = (p_n_pheno - 2), byrow = T)

eps_11_mean <- apply(results_eps_11_mtx,MARGIN = 2,mean)
eps_11_sd <- apply(results_eps_11_mtx,MARGIN = 2,sd)
eps_11_lb <- eps_11_mean - 3*eps_11_sd
eps_11_ub <- eps_11_mean + 3*eps_11_sd
eps_11_within_range <- between(results_list$param_eps_11,eps_11_lb,eps_11_ub)
table(eps_11_within_range)

eps_10_mean <- apply(results_eps_10_mtx,MARGIN = 2,mean)
eps_10_sd <- apply(results_eps_10_mtx,MARGIN = 2,sd)
eps_10_lb <- eps_10_mean - 3*eps_10_sd
eps_10_ub <- eps_10_mean + 3*eps_10_sd
eps_10_within_range <- between(results_list$param_eps_10,eps_10_lb,eps_10_ub)
table(eps_10_within_range)

# Epsilon 11
y_axis_bounds <- c(min(min(eps_11_lb),min(results_list$param_eps_11)),max(max(eps_11_ub),max(results_list$param_eps_11)))*1.1

xval = 1:p_n_pheno
which_in_bound = xval[between(results_list$param_eps_11,eps_11_lb,eps_11_ub)]

plot(1,1, type = "n", xlab = "Wavelet scale-loc", ylab = "Probability", main = "Epsilon 11", xaxt = "n"
     ,xlim = c(0,p_n_pheno)
     ,ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j]-0.5, which_in_bound[j]-0.5, which_in_bound[j]+0.5, which_in_bound[j]+0.5), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}

axis(1,at = 2^(3:10),labels = 2^(3:10),las = 2)
points(results_list$param_eps_11,type = "l")
points(eps_11_mean,col = "red", type = "l")
points(eps_11_ub,col = "blue", type = "l", lty = 2)
points(eps_11_lb,col = "blue", type = "l", lty = 2)

# Epsilon 10
y_axis_bounds <- c(min(min(eps_10_lb),min(results_list$param_eps_10)),max(max(eps_10_ub),max(results_list$param_eps_10)))*1.1

xval = 1:p_n_pheno
which_in_bound = xval[between(results_list$param_eps_10,eps_10_lb,eps_10_ub)]

plot(1,1, type = "n", xlab = "Wavelet scale-loc", ylab = "Probability", main = "Epsilon 10", xaxt = "n"
     ,xlim = c(0,p_n_pheno)
     ,ylim = y_axis_bounds)
if(length(which_in_bound) > 0){
  for(j in 1:length(which_in_bound)){
    polygon(c(which_in_bound[j]-0.5, which_in_bound[j]-0.5, which_in_bound[j]+0.5, which_in_bound[j]+0.5), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
            , col ="pink", border = NA)
  }
}

axis(1,at = 2^(3:10),labels = 2^(3:10),las = 2)
points(results_list$param_eps_10,type = "l")
points(eps_10_mean,col = "red", type = "l")
points(eps_10_ub,col = "blue", type = "l", lty = 2)
points(eps_10_lb,col = "blue", type = "l", lty = 2)

# Try and reconcile simulated vs retrieved transition proportions
get_gamma_transition_props <- function(eps_no_scale_data, tying_grp_vect, n_pheno){
  gamma_and_parent <- matrix(c(eps_no_scale_data[get_parent_indices(1:(n_pheno-1))],eps_no_scale_data)
                             , nrow = 2, ncol = (n_pheno-1)
                             , byrow = T)
  gamma_and_parent_dt <- as.data.table(t(gamma_and_parent))
  setnames(gamma_and_parent_dt,names(gamma_and_parent_dt),c("Parent","Child"))
  gamma_and_parent_dt[,"Transition" := paste0(Child,Parent)]
  gamma_and_parent_dt[,"TreeLvl" := findInterval(.I + 1,tying_grp_vect[-1])]

  # Exclude the root
  gamma_and_parent_dt <- copy(gamma_and_parent_dt[-1])
  gamma_and_parent_stats <- dcast.data.table(
    gamma_and_parent_dt[,.N,by = .(Transition,TreeLvl)]
    ,formula = TreeLvl ~ Transition
    ,value.var = "N"
    ,fill = 0
  )

  transitions_not_there <- setdiff(c("11","01","10","00"),names(gamma_and_parent_stats))
  if(length(transitions_not_there) > 0){
    for(col in transitions_not_there){
      gamma_and_parent_stats[,(col) := 0]
    }
  }

  gamma_and_parent_stats[,"Total" := apply(.SD,1,sum),.SDcols = 2:5]
  return(gamma_and_parent_stats[])
}

# Now cycle through each of the simulations and find some 'good' and 'bad' situations
eps_11_lists <- list()
eps_10_lists <- list()
gamma_seq_lists <- list()
gamma_transition_lists <- list()
for(i in 1:p_num_sims){
  eps_11_lists[[i]] <- c(0,0,results_list$results_eps_11[[i]])[p_tying_grp]
  eps_10_lists[[i]] <- c(0,0,results_list$results_eps_10[[i]])[p_tying_grp]
  gamma_seq_lists[[i]] <- results_list$results_gamma_seq[[i]][-1]
  gamma_transition_lists[[i]] <- get_gamma_transition_props(gamma_seq_lists[[i]]
                                                            ,tying_grp_vect = p_tying_grp
                                                            ,n_pheno = p_n_pheno)
  # Add simulated data
  gamma_transition_lists[[i]][,c("sim_eps_11","sim_eps_10") :=
                                .(`11`/(`11`+`01`)
                                  ,`10`/(`10`+`00`))]
  gamma_transition_lists[[i]][,c("param_eps_11","param_eps_10") :=
                                .(p_grped_eps_11[-(1:2)]
                                  ,p_grped_eps_10[-(1:2)])]
  gamma_transition_lists[[i]][,c("algo_eps_11","algo_eps_10") :=
                                .(round(eps_11_lists[[i]][-(1:2)],2)
                                  ,round(eps_10_lists[[i]][-(1:2)],2))]
}

distanceList <- lapply(gamma_transition_lists,function(x){
  x[,sqrt(sum((sim_eps_11-algo_eps_11)^2+(sim_eps_10-algo_eps_10)^2))]
})
best_3 <- order(unlist(distanceList))[1:3]
worst_3 <- order(-unlist(distanceList))[1:3]

gamma_transition_lists[best_3]
gamma_transition_lists[worst_3]

# But these just seem to be 'least worst' because the simulations are closer to 0.7 than others
# Can we find any where the algo has predicted somewhere around the 0.7, and 0.3 that we want?

distanceList <- lapply(gamma_transition_lists,function(x){
  x[,sqrt(sum((param_eps_11-algo_eps_11)^2+(param_eps_10-algo_eps_10)^2))]
})
best_3 <- order(unlist(distanceList))[1:3]
worst_3 <- order(-unlist(distanceList))[1:3]

gamma_transition_lists[best_3]
gamma_transition_lists[worst_3]

sim_09_01_good_1 <- list(results_pi = results_list$results_pi[[best_3[2]]]
                         ,results_eps_11 = results_list$results_eps_11[[best_3[2]]]
                         ,results_eps_10 = results_list$results_eps_10[[best_3[2]]]
                         ,results_gamma_seq = results_list$results_gamma_seq[[best_3[2]]]
                         ,results_g_seq = results_list$results_g_seq[[best_3[2]]]
                         ,results_y_mtx = results_list$results_y_mtx[[best_3[2]]])
sim_09_01_good_2 <- list(results_pi = results_list$results_pi[[best_3[3]]]
                         ,results_eps_11 = results_list$results_eps_11[[best_3[3]]]
                         ,results_eps_10 = results_list$results_eps_10[[best_3[3]]]
                         ,results_gamma_seq = results_list$results_gamma_seq[[best_3[3]]]
                         ,results_g_seq = results_list$results_g_seq[[best_3[3]]]
                         ,results_y_mtx = results_list$results_y_mtx[[best_3[3]]])

sim_09_01_bad_1 <- list(results_pi = results_list$results_pi[[worst_3[2]]]
                        ,results_eps_11 = results_list$results_eps_11[[worst_3[2]]]
                        ,results_eps_10 = results_list$results_eps_10[[worst_3[2]]]
                        ,results_gamma_seq = results_list$results_gamma_seq[[worst_3[2]]]
                        ,results_g_seq = results_list$results_g_seq[[worst_3[2]]]
                        ,results_y_mtx = results_list$results_y_mtx[[worst_3[2]]])
sim_09_01_bad_2 <- list(results_pi = results_list$results_pi[[worst_3[3]]]
                        ,results_eps_11 = results_list$results_eps_11[[worst_3[3]]]
                        ,results_eps_10 = results_list$results_eps_10[[worst_3[3]]]
                        ,results_gamma_seq = results_list$results_gamma_seq[[worst_3[3]]]
                        ,results_g_seq = results_list$results_g_seq[[worst_3[3]]]
                        ,results_y_mtx = results_list$results_y_mtx[[worst_3[3]]])
sim_09_01_bad_nan <- list(results_pi = results_list$results_pi[[worst_3[1]]]
                          ,results_eps_11 = results_list$results_eps_11[[worst_3[1]]]
                          ,results_eps_10 = results_list$results_eps_10[[worst_3[1]]]
                          ,results_gamma_seq = results_list$results_gamma_seq[[worst_3[1]]]
                          ,results_g_seq = results_list$results_g_seq[[worst_3[1]]]
                          ,results_y_mtx = results_list$results_y_mtx[[worst_3[1]]])

## Save these examples for usage.
saveRDS(sim_07_03_good,"~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/sim_07_03_good.RDS",compress = T)
saveRDS(sim_07_03_bad,"~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/sim_07_03_bad.RDS",compress = T)

saveRDS(sim_09_01_good_1,"~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/sim_09_01_good_1.RDS",compress = T)
saveRDS(sim_09_01_good_2,"~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/sim_09_01_good_2.RDS",compress = T)
saveRDS(sim_09_01_bad_1,"~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/sim_09_01_bad_1.RDS",compress = T)
saveRDS(sim_09_01_bad_2,"~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/sim_09_01_bad_2.RDS",compress = T)
saveRDS(sim_09_01_bad_nan,"~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/sim_09_01_bad_nan.RDS",compress = T)

write.table(sim_07_03_good$results_y_mtx, file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/WCs_07_03_good.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
write.table(sim_07_03_bad$results_y_mtx, file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/WCs_07_03_bad.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
write.table(sim_09_01_good_1$results_y_mtx, file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/WCs_09_01_good_1.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
write.table(sim_09_01_good_2$results_y_mtx, file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/WCs_09_01_good_2.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
write.table(sim_09_01_bad_1$results_y_mtx, file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/WCs_09_01_bad_1.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
write.table(sim_09_01_bad_2$results_y_mtx, file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/WCs_09_01_bad_2.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
write.table(sim_09_01_bad_nan$results_y_mtx, file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/WCs_09_01_bad_nan.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)

cat(c("sim2","A","A",sim_07_03_good$results_g_seq), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/g_07_03_good.cis.geno")
cat(c("sim2","A","A",sim_07_03_bad$results_g_seq), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/g_07_03_bad.cis.geno")
cat(c("sim2","A","A",sim_09_01_good_1$results_g_seq), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/g_09_01_good_1.cis.geno")
cat(c("sim2","A","A",sim_09_01_good_2$results_g_seq), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/g_09_01_good_2.cis.geno")
cat(c("sim2","A","A",sim_09_01_bad_1$results_g_seq), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/g_09_01_bad_1.cis.geno")
cat(c("sim2","A","A",sim_09_01_bad_2$results_g_seq), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/g_09_01_bad_2.cis.geno")
cat(c("sim2","A","A",sim_09_01_bad_nan$results_g_seq), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/g_09_01_bad_nan.cis.geno")


# A toy example -----------------------------------------------------------

# Normal beta

source("code/sim2_script.R")
attempt_no=15

p_n_ind = 70
p_n_pheno = 128
p_tying_grp = c(1,2,3,33,65)
p_param_pi_00 = 0
p_param_pi_11 = 1
p_grped_eps_11 = c(rep(0.5,(length(p_tying_grp)-1)),0)
p_grped_eps_10 = c(rep(0.5,(length(p_tying_grp)-1)),0)
p_coeff_mu = 0
p_coeff_beta = 2
p_param_gi_prob = 0.4
p_param_sigma_beta = 0.5
p_num_sims = 1
p_seed = 20

results_list <- run_sim2(
  n_ind = p_n_ind
  ,n_pheno = p_n_pheno
  ,tying_grp = p_tying_grp
  ,param_pi_00 = p_param_pi_00
  ,param_pi_11 = p_param_pi_11
  ,grped_eps_11 = p_grped_eps_11
  ,grped_eps_10 = p_grped_eps_10
  ,coeff_mu = p_coeff_mu
  ,coeff_beta = p_coeff_beta
  ,param_gi_prob = p_param_gi_prob
  ,param_sigma_beta = p_param_sigma_beta
  ,num_sims = p_num_sims
  ,seed = p_seed)

tree_plot(results_list$results_gamma_seq[[1]],plot_title = "gamma_seq")

write.table(results_list$results_y_mtx[[1]][,1:16], file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_16wc.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
cat(c("sim2","A","A",results_list$results_g_seq[[1]]), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_16wc.cis.geno")

write.table(results_list$results_y_mtx[[1]][,1:32], file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_32wc.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
saveRDS(results_list$results_gamma_seq[[1]],"~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_128gamma.RDS")

# The full 1024 size of the large beta
source("code/sim2_script.R")
attempt_no=15

p_n_ind = 70
p_n_pheno = 1024
p_tying_grp = c(1,2,3,33,65,129,257,513)
p_param_pi_00 = 0
p_param_pi_11 = 1
p_grped_eps_11 = c(rep(0.5,(length(p_tying_grp)-1)),0)
p_grped_eps_10 = c(rep(0.5,(length(p_tying_grp)-1)),0)
p_coeff_mu = 0
p_coeff_beta = 2
p_param_gi_prob = 0.4
p_param_sigma_beta = 0.5
p_num_sims = 1
p_seed = 20

results_list <- run_sim2(
  n_ind = p_n_ind
  ,n_pheno = p_n_pheno
  ,tying_grp = p_tying_grp
  ,param_pi_00 = p_param_pi_00
  ,param_pi_11 = p_param_pi_11
  ,grped_eps_11 = p_grped_eps_11
  ,grped_eps_10 = p_grped_eps_10
  ,coeff_mu = p_coeff_mu
  ,coeff_beta = p_coeff_beta
  ,param_gi_prob = p_param_gi_prob
  ,param_sigma_beta = p_param_sigma_beta
  ,num_sims = p_num_sims
  ,seed = p_seed)

write.table(results_list$results_y_mtx[[1]], file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_1024wc.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
cat(c("sim2","A","A",results_list$results_g_seq[[1]]), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_1024wc.cis.geno")
saveRDS(results_list$results_gamma_seq[[1]],"~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_1024gamma.RDS")

# Small beta
source("code/sim2_script.R")
attempt_no=15

p_n_ind = 70
p_n_pheno = 128
p_tying_grp = c(1,2,3,33,65)
p_param_pi_00 = 0
p_param_pi_11 = 1
p_grped_eps_11 = c(rep(0.5,(length(p_tying_grp)-1)),0)
p_grped_eps_10 = c(rep(0.5,(length(p_tying_grp)-1)),0)
p_coeff_mu = 0
p_coeff_beta = 0.1
p_param_gi_prob = 0.4
p_param_sigma_beta = 0.5
p_num_sims = 1
p_seed = 20

results_list <- run_sim2(
  n_ind = p_n_ind
  ,n_pheno = p_n_pheno
  ,tying_grp = p_tying_grp
  ,param_pi_00 = p_param_pi_00
  ,param_pi_11 = p_param_pi_11
  ,grped_eps_11 = p_grped_eps_11
  ,grped_eps_10 = p_grped_eps_10
  ,coeff_mu = p_coeff_mu
  ,coeff_beta = p_coeff_beta
  ,param_gi_prob = p_param_gi_prob
  ,param_sigma_beta = p_param_sigma_beta
  ,num_sims = p_num_sims
  ,seed = p_seed)

tree_plot(results_list$results_gamma_seq[[1]],plot_title = "gamma_seq")

write.table(results_list$results_y_mtx[[1]][,1:16], file= "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_16wc.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)
cat(c("sim2","A","A",results_list$results_g_seq[[1]]), file = "~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_16wc.cis.geno")
