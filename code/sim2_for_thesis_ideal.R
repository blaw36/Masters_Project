### Reproducing idealistic algorithm correctness simulations for thesis


# Very idealistic ---------------------------------------------------------

# logBF: -2, 2
# attempt to run on 256 nodes, if possible
# attempt to extract no tying params, and then tied params, to approximate proportions, if possible


# Preamble ----------------------------------------------------------------

# Clear environment and do some fresh simulations
rm(list = ls());gc();cat("\014");
source("code/WaveQTL/waveqtl_hmt_test_calc_sumlog_func.R")
source("code/sim2_script.R")
library(data.table)
sim_mode = T
num_sims = 100

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

if(!sim_mode){
  # Read in gamma data ------------------------------------------------------

  ### Could we show an idealistic situation, where it can retrieve all parameters (no noise, no tying)?

  # This is a simulated gamma sequence
  full_gamma <- readRDS("~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_1024gamma.RDS")[-1]

  # 256 nodes, in R
  # logBF: -1, 1

  ## Investigate empirical gamma proportions
  # num_nodes <- 256
  # tying_groups <- c(1,2,4,8,16,32,64,128)
  num_nodes <- 32
  tying_groups <- c(1,2,4,8,16)
  # tying_groups <- c(1,16)
  # tying_groups <- 1:31

  gamma_seq <- copy(full_gamma[1:(num_nodes - 1)])
  p_logbf <- tree_plot(c(0,gamma_seq),yaxis_lims = c(min(gamma_seq),max(gamma_seq)),plot_title = "logBFs")

  gamma_trans_tree_tie <- get_gamma_transition_props(gamma_seq,tying_grp_vect = c(1,tying_groups+1),n_pheno = num_nodes)
  gamma_trans_tree_tie[,c("sim_eps_11","sim_eps_10") :=
                         .(`11`/(`11`+`01`)
                           ,`10`/(`10`+`00`))]

  gamma_256_1_n1 <- copy(gamma_seq)
  gamma_256_1_n1[gamma_256_1_n1 == 0] <- -1

  # Initial conditions
  eps_tree <- generate_eps_from_empirical(gamma_trans_tree_tie,n_pheno = num_nodes,nan_default = 0.5)

  tree_tie_256_res <- list()
  bf_mult <- 0.1

  logBF_input <- gamma_256_1_n1*bf_mult
  tree_tie_256_res <-
    waveQTL_HMT_R(tying_groups = tying_groups
                  , logBFs = logBF_input
                  # ,eps_table = log(eps_tree)
                  # , init_pi = 0.99999)
    )

  plot(cbind(round(exp(tree_tie_256_res$pp_i$`1`),5),gamma_256_1_n1),xlab = "marginal posterior prob of sl", ylab = "logBF")
  print(cbind(gamma_trans_tree_tie[,.(TreeLvl,sim_eps_11,sim_eps_10)]
              ,round(exp(tree_tie_256_res$eps[c(tying_groups[-1])][,.(`11`,`10`)]),4)))

  plot(x = cbind(gamma_trans_tree_tie[,.(TreeLvl,sim_eps_11,sim_eps_10)]
                 ,round(exp(tree_tie_256_res$eps[c(tying_groups[-1])][,.(`11`,`10`)]),4))[["sim_eps_11"]]
       ,y = cbind(gamma_trans_tree_tie[,.(TreeLvl,sim_eps_11,sim_eps_10)]
                  ,round(exp(tree_tie_256_res$eps[c(tying_groups[-1])][,.(`11`,`10`)]),4))[["11"]]
       ,xlab = "", ylabs = "Actual")
  plot(x = cbind(gamma_trans_tree_tie[,.(TreeLvl,sim_eps_11,sim_eps_10)]
                 ,round(exp(tree_tie_256_res$eps[c(tying_groups[-1])][,.(`11`,`10`)]),4))[["sim_eps_10"]]
       ,y = cbind(gamma_trans_tree_tie[,.(TreeLvl,sim_eps_11,sim_eps_10)]
                  ,round(exp(tree_tie_256_res$eps[c(tying_groups[-1])][,.(`11`,`10`)]),4))[["10"]]
       ,xlab = "Simulated", ylabs = "Actual")

}

# Loop of gamma -----------------------------------------------------------

# Simulate a gamma sequence
# Set logBFs to 1, -1 depending on gamma:
# Gamma = 1 => logBF = 1.1
# Gamma = 0 => logBF = -1.1

# Run on 32 nodes, with tying groups as per above
n_pheno <- 32
param_pi_00 <- 0 # Scaling coeff
param_pi_11 <- 1 # Head of tree
sim_eps_11 <- 0.5
sim_eps_10 <- 0.5
tying_grp <- c(1,2,3,5,9,17)

simulate_gammas_ideal <- function(
  n_pheno = 32
  ,param_pi_00 = 0
  ,param_pi_11 = 1
  ,sim_eps_11 = 0.5
  ,sim_eps_10 = 0.5
  ,tying_grp = c(1,2,3,5,9,17)
  ,num_sims = 100
){

  ### Step 1: Generate Parameters ----
  tying_grp_alt <- (tying_grp-1)[-1]
  num_tying_grps <- length(tying_grp)
  grped_eps_11 <- rep(sim_eps_11,num_tying_grps)
  grped_eps_10 <- rep(sim_eps_10,num_tying_grps)

  results_list <- list()

  for(j in 1:num_sims){
    print(paste0("Iteration: ",j))

    # Eps_11 (elements 2 -> 1023 of tree)
    # grped_eps_11 <- round(runif(n=num_tying_grps)*100)/100
    param_eps_11 <- numeric()
    for(i in 1:(num_tying_grps - 1)){
      indx_start <- tying_grp[i]
      indx_end <- tying_grp[i+1]
      param_eps_11[indx_start:indx_end] <- grped_eps_11[i]
    }
    param_eps_11[tying_grp[num_tying_grps]:n_pheno] <- grped_eps_11[num_tying_grps]
    param_eps_11 <- param_eps_11[-(1:2)]

    # Eps_10 (elements 2 -> 1023 of tree)
    # grped_eps_10 <- round(runif(n=num_tying_grps)*100)/100
    param_eps_10 <- numeric()
    for(i in 1:(num_tying_grps - 1)){
      indx_start <- tying_grp[i]
      indx_end <- tying_grp[i+1]
      param_eps_10[indx_start:indx_end] <- grped_eps_10[i]
    }
    param_eps_10[tying_grp[num_tying_grps]:n_pheno] <- grped_eps_10[num_tying_grps]
    param_eps_10 <- param_eps_10[-(1:2)]


    gamma_seq <- numeric()
    rand_seq <- runif(n_pheno)

    # Scaling coefficient
    gamma_seq[1] <- ifelse(rand_seq[1] < param_pi_00, 1, 0)

    # Head of tree
    gamma_seq[2] <- ifelse(rand_seq[2] < param_pi_11, 1, 0)

    # i is the index of tree, where i = 1 is the head of the tree
    # Using this notation because that's how 'get_parent_indices' has been written
    for(i in 2:(n_pheno - 1)){
      indx <- i
      parent_indx <- get_parent_indices(indx)
      parent_gamma <- gamma_seq[parent_indx + 1]

      if(parent_gamma == 1){
        sl_prob <- param_eps_11[indx - 1]
      }else if(parent_gamma == 0){
        sl_prob <- param_eps_10[indx - 1]
      }

      gamma_seq[i+1] <- ifelse(rand_seq[i+1] < sl_prob, 1, 0)
    }

    ### Step 2: Properties of gamma ----
    gamma_seq <- copy(gamma_seq[2:(n_pheno)])
    # p_logbf <- tree_plot(c(0,gamma_seq),yaxis_lims = c(min(gamma_seq),max(gamma_seq)),plot_title = "logBFs")

    gamma_trans_tree_tie <- get_gamma_transition_props(gamma_seq,tying_grp_vect = c(1,tying_grp_alt+1),n_pheno = n_pheno)
    gamma_trans_tree_tie[,c("sim_eps_11","sim_eps_10") :=
                           .(`11`/(`11`+`01`)
                             ,`10`/(`10`+`00`))]

    gamma_256_1_n1 <- copy(gamma_seq)
    gamma_256_1_n1[gamma_256_1_n1 == 0] <- -1

    # # Initial conditions
    # eps_tree <- generate_eps_from_empirical(gamma_trans_tree_tie,n_pheno = n_pheno,nan_default = 0.5)

    tree_tie_256_res <- list()
    bf_mult <- 1

    ### Step 3: Run R script ----
    logBF_input <- gamma_256_1_n1*bf_mult
    tree_tie_256_res <-
      waveQTL_HMT_R(tying_groups = tying_grp_alt
                    , logBFs = logBF_input
                    # ,eps_table = log(eps_tree)
                    # , init_pi = 0.99999)
      )

    result_comparison <- cbind(gamma_trans_tree_tie[,.(TreeLvl,sim_eps_11,sim_eps_10)]
                               ,round(exp(tree_tie_256_res$eps[c(tying_grp_alt[-1])][,.(`11`,`10`)]),4))

    results_list[[j]] <- result_comparison

  }

  eps_11_lvl1_act_mean <- mean(unlist(lapply(results_list,function(x){x[["sim_eps_11"]][1]})),na.rm = T)
  eps_11_lvl1_act_sd <- sd(unlist(lapply(results_list,function(x){x[["sim_eps_11"]][1]})),na.rm = T)
  eps_11_lvl2_act_mean <- mean(unlist(lapply(results_list,function(x){x[["sim_eps_11"]][2]})),na.rm = T)
  eps_11_lvl2_act_sd <- sd(unlist(lapply(results_list,function(x){x[["sim_eps_11"]][2]})),na.rm = T)
  eps_11_lvl3_act_mean <- mean(unlist(lapply(results_list,function(x){x[["sim_eps_11"]][3]})),na.rm = T)
  eps_11_lvl3_act_sd <- sd(unlist(lapply(results_list,function(x){x[["sim_eps_11"]][3]})),na.rm = T)
  eps_11_lvl4_act_mean <- mean(unlist(lapply(results_list,function(x){x[["sim_eps_11"]][4]})),na.rm = T)
  eps_11_lvl4_act_sd <- sd(unlist(lapply(results_list,function(x){x[["sim_eps_11"]][4]})),na.rm = T)
  eps_11_lvl1_algo_mean <- mean(unlist(lapply(results_list,function(x){x[["11"]][1]})),na.rm = T)
  eps_11_lvl1_algo_sd <- sd(unlist(lapply(results_list,function(x){x[["11"]][1]})),na.rm = T)
  eps_11_lvl2_algo_mean <- mean(unlist(lapply(results_list,function(x){x[["11"]][2]})),na.rm = T)
  eps_11_lvl2_algo_sd <- sd(unlist(lapply(results_list,function(x){x[["11"]][2]})),na.rm = T)
  eps_11_lvl3_algo_mean <- mean(unlist(lapply(results_list,function(x){x[["11"]][3]})),na.rm = T)
  eps_11_lvl3_algo_sd <- sd(unlist(lapply(results_list,function(x){x[["11"]][3]})),na.rm = T)
  eps_11_lvl4_algo_mean <- mean(unlist(lapply(results_list,function(x){x[["11"]][4]})),na.rm = T)
  eps_11_lvl4_algo_sd <- sd(unlist(lapply(results_list,function(x){x[["11"]][4]})),na.rm = T)

  eps_10_lvl1_act_mean <- mean(unlist(lapply(results_list,function(x){x[["sim_eps_10"]][1]})),na.rm = T)
  eps_10_lvl1_act_sd <- sd(unlist(lapply(results_list,function(x){x[["sim_eps_10"]][1]})),na.rm = T)
  eps_10_lvl1_algo_mean <- mean(unlist(lapply(results_list,function(x){x[["10"]][1]})),na.rm = T)
  eps_10_lvl1_algo_sd <- sd(unlist(lapply(results_list,function(x){x[["10"]][1]})),na.rm = T)
  eps_10_lvl2_act_mean <- mean(unlist(lapply(results_list,function(x){x[["sim_eps_10"]][2]})),na.rm = T)
  eps_10_lvl2_act_sd <- sd(unlist(lapply(results_list,function(x){x[["sim_eps_10"]][2]})),na.rm = T)
  eps_10_lvl2_algo_mean <- mean(unlist(lapply(results_list,function(x){x[["10"]][2]})),na.rm = T)
  eps_10_lvl2_algo_sd <- sd(unlist(lapply(results_list,function(x){x[["10"]][2]})),na.rm = T)
  eps_10_lvl3_act_mean <- mean(unlist(lapply(results_list,function(x){x[["sim_eps_10"]][3]})),na.rm = T)
  eps_10_lvl3_act_sd <- sd(unlist(lapply(results_list,function(x){x[["sim_eps_10"]][3]})),na.rm = T)
  eps_10_lvl3_algo_mean <- mean(unlist(lapply(results_list,function(x){x[["10"]][3]})),na.rm = T)
  eps_10_lvl3_algo_sd <- sd(unlist(lapply(results_list,function(x){x[["10"]][3]})),na.rm = T)
  eps_10_lvl4_act_mean <- mean(unlist(lapply(results_list,function(x){x[["sim_eps_10"]][4]})),na.rm = T)
  eps_10_lvl4_act_sd <- sd(unlist(lapply(results_list,function(x){x[["sim_eps_10"]][4]})),na.rm = T)
  eps_10_lvl4_algo_mean <- mean(unlist(lapply(results_list,function(x){x[["10"]][4]})),na.rm = T)
  eps_10_lvl4_algo_sd <- sd(unlist(lapply(results_list,function(x){x[["10"]][4]})),na.rm = T)

  # No first transition for 10 as this is coming from a tree root of 1
  eps_10_ub <- c(
    rep(eps_10_lvl2_algo_mean+2*eps_10_lvl2_algo_sd,4)
    ,rep(eps_10_lvl3_algo_mean+2*eps_10_lvl3_algo_sd,8)
    ,rep(eps_10_lvl4_algo_mean+2*eps_10_lvl4_algo_sd,16))
  eps_10_lb <- c(
    rep(eps_10_lvl2_algo_mean-2*eps_10_lvl2_algo_sd,4)
    ,rep(eps_10_lvl3_algo_mean-2*eps_10_lvl3_algo_sd,8)
    ,rep(eps_10_lvl4_algo_mean-2*eps_10_lvl4_algo_sd,16))
  eps_10_mean <- c(
    rep(eps_10_lvl2_algo_mean,4)
    ,rep(eps_10_lvl3_algo_mean,8)
    ,rep(eps_10_lvl4_algo_mean,16))

  plot(1,1,"n"
       ,xlim = c(2,31)
       ,ylim = c(min(eps_10_lb) * 0.99, max(eps_10_ub) * 1.01)
       ,main = "Epsilon 10 parameter"
       ,xaxt = "n")
  axis(side = 1, at = tying_grp_alt[-1], labels = tying_grp_alt[-1])
  lines(eps_10_mean, col = "blue")
  lines(eps_10_ub, col = "cyan")
  lines(eps_10_lb, col = "cyan")
  lines(rep(param_eps_10,30), col = "red")
  p_e_10 <- recordPlot()

  eps_11_ub <- c(
    rep(eps_11_lvl1_algo_mean+2*eps_11_lvl1_algo_sd,2)
    ,rep(eps_11_lvl2_algo_mean+2*eps_11_lvl2_algo_sd,4)
    ,rep(eps_11_lvl3_algo_mean+2*eps_11_lvl3_algo_sd,8)
    ,rep(eps_11_lvl4_algo_mean+2*eps_11_lvl4_algo_sd,16))
  eps_11_lb <- c(
    rep(eps_11_lvl1_algo_mean-2*eps_11_lvl1_algo_sd,2)
    ,rep(eps_11_lvl2_algo_mean-2*eps_11_lvl2_algo_sd,4)
    ,rep(eps_11_lvl3_algo_mean-2*eps_11_lvl3_algo_sd,8)
    ,rep(eps_11_lvl4_algo_mean-2*eps_11_lvl4_algo_sd,16))
  eps_11_mean <- c(
    rep(eps_11_lvl1_algo_mean,2)
    ,rep(eps_11_lvl2_algo_mean,4)
    ,rep(eps_11_lvl3_algo_mean,8)
    ,rep(eps_11_lvl4_algo_mean,16))

  plot(1,1,"n"
       ,xlim = c(1,31)
       ,ylim = c(min(eps_11_lb) * 0.99, max(eps_11_ub) * 1.01)
       ,main = "Epsilon 11 parameter"
       ,xaxt = "n")
  axis(side = 1, at = tying_grp_alt, labels = tying_grp_alt)
  lines(eps_11_mean, col = "blue")
  lines(eps_11_ub, col = "cyan")
  lines(eps_11_lb, col = "cyan")
  lines(rep(param_eps_11,31), col = "red")
  p_e_11 <- recordPlot()

  return(list(
    p_e_11 = p_e_11
    ,p_e_10 = p_e_10
    ,results_list = results_list
  ))
}

eps_0505 <- simulate_gammas_ideal(
  n_pheno = 32
  ,param_pi_00 = 0 # Scaling coeff
  ,param_pi_11 = 1 # Head of tree
  ,sim_eps_11 = 0.5
  ,sim_eps_10 = 0.5
  ,tying_grp = c(1,2,3,5,9,17)
  ,num_sims = 1000
)
eps_0703 <- simulate_gammas_ideal(
  n_pheno = 32
  ,param_pi_00 = 0 # Scaling coeff
  ,param_pi_11 = 1 # Head of tree
  ,sim_eps_11 = 0.7
  ,sim_eps_10 = 0.3
  ,tying_grp = c(1,2,3,5,9,17)
  ,num_sims = 1000
)
eps_0802 <- simulate_gammas_ideal(
  n_pheno = 32
  ,param_pi_00 = 0 # Scaling coeff
  ,param_pi_11 = 1 # Head of tree
  ,sim_eps_11 = 0.8
  ,sim_eps_10 = 0.2
  ,tying_grp = c(1,2,3,5,9,17)
  ,num_sims = 1000
)

saveRDS(eps_0505,"data/20190823_sim2_ideal_bf1_1000_0505.RDS", compress = T)
saveRDS(eps_0703,"data/20190823_sim2_ideal_bf1_1000_0703.RDS", compress = T)
saveRDS(eps_0802,"data/20190823_sim2_ideal_bf1_1000_0802.RDS", compress = T)
