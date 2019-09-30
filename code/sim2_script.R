# rm(list = ls());gc();cat("\014");

get_parent_indices <- function(indx, tree_root_indx = 1){
  return_indices <- indx
  return_indices <- (return_indices - tree_root_indx + 1) %/% 2

  # root of tree doesn't return any index
  return_indices[which(indx == tree_root_indx)] <- NA_integer_
  return(return_indices)
}

tree_plot <- function(data, yaxis_lims = c(0,1),plot_title){

  # Interleaves vectors with two bookend 0s, and a 0 between each current element
  # (for spacing of elements like a tree in a plot)
  vector_centeriser <- function(vect){
    in_between_zeros <- length(vect) - 1

    res_vect <- c(0,vect[1],0)
    if(in_between_zeros > 0){
      for(i in 1:in_between_zeros){
        res_vect <- c(res_vect,vect[i + 1],0)
      }
    }
    return(res_vect)
  }

  num_lvls <- floor(log2(length(data))) + 1
  par(mfrow = c(num_lvls,1),mar = c(1,1,1,1))
  plot(vector_centeriser(data[1]),type = "h",ylab = "",axes = F,ylim=yaxis_lims,main=plot_title) # scaling coeff
  plot(vector_centeriser(data[2]),type = "h",ylab = "",axes = F,ylim=yaxis_lims) # head of tree
  for(i in 1:(num_lvls-2)){
    plot(vector_centeriser(data[((2^i)+1):(2^(i+1))]),type = "h",ylab = "",axes = F,ylim=yaxis_lims)
  }
  p <- recordPlot()
  return(p)
}

run_sim2 <- function(
  n_ind = 70
  ,n_pheno = 1024
  ,tying_grp = c(1,2,3,5,9,17,33,65,129,257,513)
  ,param_pi_00 = 0
  ,param_pi_11 = 1
  ,grped_eps_11 = c(rep(0.9,10),0)
  ,grped_eps_10 = c(rep(0.1,10),0)
  ,coeff_mu = 0
  ,coeff_beta = 2
  ,param_gi_prob = 0.4
  ,param_sigma_beta = 0.5
  ,num_sims = 100
  ,seed = 20
  ,showOutput = F){

  set.seed(seed)
  setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")

  # Helper function to get parent indices of a given index
  # From "../code/WaveQTL/waveqtl_hmt_test_calc_sumlog.R"
  get_parent_indices <- function(indx, tree_root_indx = 1){
    return_indices <- indx
    return_indices <- (return_indices - tree_root_indx + 1) %/% 2

    # root of tree doesn't return any index
    return_indices[which(indx == tree_root_indx)] <- NA_integer_
    return(return_indices)
  }

  num_tying_grps <- length(tying_grp)

  results_pi <- list()
  results_eps_11 <- list()
  results_eps_10 <- list()
  results_gamma_seq <- list()
  results_g_seq <- list()
  results_y_mtx <- list()

  # Epsilon elements are for elements 3 -> 1024 (= 2 -> 1023 of tree)
  # Generate 1024 elements as per the above tying_grp description,
  # then just cut out the first two levels of elements (as neither the
  # scaling coefficient or head of tree need epsilons)

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

  for(its in 1:num_sims){
    cat(paste0("running iteration ",its,"...\n"))

    ### Step 1: Generate Parameters ----
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

    ### Step 2: Generate beta,mu,eps,g,y ----
    beta_seq <- rep(0,n_pheno)
    beta_seq[which(gamma_seq == 1)] <- coeff_beta

    mu_seq <- rep(0,n_pheno)
    mu_seq[which(gamma_seq == 1)] <- coeff_mu
    mu_mtx <- matrix(rep(mu_seq,70),nrow = n_ind,ncol = n_pheno,byrow = T)

    eps_seq <- rnorm(n_pheno*n_ind,mean = 0,sd = param_sigma_beta)
    eps_mtx <- matrix(eps_seq,nrow = n_ind,ncol = n_pheno,byrow = T)

    g_seq <- rbinom(n_ind,size = 2,prob = param_gi_prob)
    cat(c("sim2","A","A",g_seq), file = paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/sim2.cis.geno"))
    beta_mtx <- g_seq %*% t(beta_seq)

    y_mtx <- mu_mtx + beta_mtx + eps_mtx
    write.table(y_mtx, file= paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sim2_WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
    cat(rep(1,n_pheno), file = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/use_all.txt"))
    cat(tying_grp, file = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sim_grouping.txt"))

    ### Step 3: Run through HMT ----
    command <- paste0("../../WaveQTL -gmode 1"
                      ," -group sim_grouping.txt "
                      ,"-g ../../data/dsQTL/sim2.cis.geno -p sim2_WCs.txt -u use_all.txt -o sim2_noQT -f ",n_pheno," -hmt 1")
    system(command,show.output.on.console = showOutput)

    # I DON'T ACTUALLY HAVE A WAY OF RECOVERING PI_0 ATM UNLESS I RUN IT THROUGH WAVEQTL (NO HMT) ALSO!
    ### Step 4: Get HMT results ----
    file_prefix <- "sim2_noQT"
    pi_file <- exp(as.numeric(read.table(paste0("output/",file_prefix,".fph.pi.txt"))[2]))
    eps_11_file <- exp(as.numeric(read.table(paste0("output/",file_prefix,".fph.eps_11.txt"))[4:(n_pheno + 1)]))
    eps_10_file <- exp(as.numeric(read.table(paste0("output/",file_prefix,".fph.eps_10.txt"))[4:(n_pheno + 1)]))

    results_pi[[its]] <- pi_file
    results_eps_11[[its]] <- eps_11_file
    results_eps_10[[its]] <- eps_10_file
    results_gamma_seq[[its]] <- gamma_seq
    results_g_seq[[its]] <- g_seq
    results_y_mtx[[its]] <- y_mtx

  }

  setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
  return(list(results_pi = results_pi
              ,results_eps_11 = results_eps_11
              ,results_eps_10 = results_eps_10
              ,results_gamma_seq = results_gamma_seq
              ,results_g_seq = results_g_seq
              ,results_y_mtx = results_y_mtx
              ,param_eps_11 = param_eps_11
              ,param_eps_10 = param_eps_10
              ,param_pi_11 = param_pi_11
              ,param_pi_00 = param_pi_00))
}

# saveRDS(results_pi,file = "sim2_pi.RDS",compress = T)
# saveRDS(results_eps_11,file = "sim2_eps_11.RDS",compress = T)
# saveRDS(results_eps_10,file = "sim2_eps_10.RDS",compress = T)

sim_analysis <- function(result_data_list, cred_ints = c(0.025,0.975),num_sims,num_pheno){
  pi_mean <- mean(unlist(result_data_list$results_pi))
  pi_sd <- sd(unlist(result_data_list$results_pi))
  # pi_range <- c(pi_mean - 3*pi_sd, pi_mean + 3*pi_sd)
  pi_range <- quantile(unlist(result_data_list$results_pi),probs = cred_ints)
  pi_range; result_data_list$param_pi_11

  results_eps_11_mtx <- matrix(unlist(result_data_list$results_eps_11)
                               ,nrow = num_sims, ncol = (num_pheno - 2), byrow = T)
  results_eps_10_mtx <- matrix(unlist(result_data_list$results_eps_10)
                               ,nrow = num_sims, ncol = (num_pheno - 2), byrow = T)

  eps_11_mean <- apply(results_eps_11_mtx,MARGIN = 2,mean)
  eps_11_sd <- apply(results_eps_11_mtx,MARGIN = 2,sd)
  # eps_11_lb <- eps_11_mean - 3*eps_11_sd
  # eps_11_ub <- eps_11_mean + 3*eps_11_sd
  eps_11_lb <- apply(results_eps_11_mtx,MARGIN = 2,quantile,probs = cred_ints[1])
  eps_11_ub <- apply(results_eps_11_mtx,MARGIN = 2,quantile,probs = cred_ints[2])
  eps_11_within_range <- between(result_data_list$param_eps_11,eps_11_lb,eps_11_ub)
  table(eps_11_within_range)

  eps_10_mean <- apply(results_eps_10_mtx,MARGIN = 2,mean)
  eps_10_sd <- apply(results_eps_10_mtx,MARGIN = 2,sd)
  eps_10_lb <- apply(results_eps_10_mtx,MARGIN = 2,quantile,probs = cred_ints[1])
  eps_10_ub <- apply(results_eps_10_mtx,MARGIN = 2,quantile,probs = cred_ints[2])
  eps_10_within_range <- between(result_data_list$param_eps_10,eps_10_lb,eps_10_ub)
  table(eps_10_within_range)

  # Epsilon 11
  y_axis_bounds <- c(max(min(min(eps_11_lb),min(result_data_list$param_eps_11))*0.9,0),min(max(max(eps_11_ub),max(result_data_list$param_eps_11))*1.1,1))

  xval = 1:num_pheno
  which_in_bound = xval[between(result_data_list$param_eps_11,eps_11_lb,eps_11_ub)]

  plot(1,1, type = "n", xlab = "Wavelet scale-loc", ylab = "Probability", main = "Epsilon 11", xaxt = "n"
       ,xlim = c(0,num_pheno)
       ,ylim = y_axis_bounds)
  if(length(which_in_bound) > 0){
    for(j in 1:length(which_in_bound)){
      polygon(c(which_in_bound[j]-0.5, which_in_bound[j]-0.5, which_in_bound[j]+0.5, which_in_bound[j]+0.5), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
              , col ="pink", border = NA)
    }
  }

  axis(1,at = 2^(3:10),labels = 2^(3:10),las = 2)
  points(result_data_list$param_eps_11,type = "l")
  points(eps_11_mean,col = "red", type = "l")
  points(eps_11_ub,col = "blue", type = "l", lty = 2)
  points(eps_11_lb,col = "blue", type = "l", lty = 2)
  eps11_plot <- recordPlot()

  # Epsilon 10
  y_axis_bounds <- c(max(min(min(eps_10_lb),min(result_data_list$param_eps_10))*0.9,0),min(max(max(eps_10_ub),max(result_data_list$param_eps_10))*1.1,1))

  xval = 1:num_pheno
  which_in_bound = xval[between(result_data_list$param_eps_10,eps_10_lb,eps_10_ub)]

  plot(1,1, type = "n", xlab = "Wavelet scale-loc", ylab = "Probability", main = "Epsilon 10", xaxt = "n"
       ,xlim = c(0,num_pheno)
       ,ylim = y_axis_bounds)
  if(length(which_in_bound) > 0){
    for(j in 1:length(which_in_bound)){
      polygon(c(which_in_bound[j]-0.5, which_in_bound[j]-0.5, which_in_bound[j]+0.5, which_in_bound[j]+0.5), c(y_axis_bounds[1], y_axis_bounds[2],  y_axis_bounds[2],  y_axis_bounds[1])
              , col ="pink", border = NA)
    }
  }

  axis(1,at = 2^(3:10),labels = 2^(3:10),las = 2)
  points(result_data_list$param_eps_10,type = "l")
  points(eps_10_mean,col = "red", type = "l")
  points(eps_10_ub,col = "blue", type = "l", lty = 2)
  points(eps_10_lb,col = "blue", type = "l", lty = 2)
  eps10_plot <- recordPlot()

  return(
    list(
      pi_range = pi_range
      ,pi_mean = pi_mean
      ,pi_sd = pi_sd
      ,results_eps_11_mtx = results_eps_11_mtx
      ,results_eps_10_mtx = results_eps_10_mtx
      ,tgt_pi_11 = result_data_list$param_pi_11
      ,tgt_11 = result_data_list$param_eps_11
      ,tgt_10 = result_data_list$param_eps_10
      ,num_sims = num_sims
      ,num_pheno = num_pheno
    )
  )
}


# Sim2, with a few customisations
run_sim2_custom <- function(
  n_ind = 70
  ,n_pheno = 1024
  ,tying_grp = c(1,2,3,5,9,17,33,65,129,257,513)
  ,param_pi_00 = 0
  ,param_pi_11 = 1
  ,grped_eps_11 = c(rep(0.9,10),0)
  ,grped_eps_10 = c(rep(0.1,10),0)
  ,coeff_mu = 0
  ,coeff_beta = 2
  ,param_gi_prob = 0.4
  ,param_sigma_beta = 0.5
  ,num_sims = 100
  ,seed = 20
  ,showOutput = F
  ,sigma_prior = 14){

  set.seed(seed)
  setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")

  # Helper function to get parent indices of a given index
  # From "../code/WaveQTL/waveqtl_hmt_test_calc_sumlog.R"
  get_parent_indices <- function(indx, tree_root_indx = 1){
    return_indices <- indx
    return_indices <- (return_indices - tree_root_indx + 1) %/% 2

    # root of tree doesn't return any index
    return_indices[which(indx == tree_root_indx)] <- NA_integer_
    return(return_indices)
  }

  num_tying_grps <- length(tying_grp)

  results_pi <- list()
  results_eps_11 <- list()
  results_eps_10 <- list()
  results_gamma_seq <- list()
  results_g_seq <- list()
  results_y_mtx <- list()

  # Epsilon elements are for elements 3 -> 1024 (= 2 -> 1023 of tree)
  # Generate 1024 elements as per the above tying_grp description,
  # then just cut out the first two levels of elements (as neither the
  # scaling coefficient or head of tree need epsilons)

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

  for(its in 1:num_sims){
    cat(paste0("running iteration ",its,"...\n"))

    ### Step 1: Generate Parameters ----
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

    ### Step 2: Generate beta,mu,eps,g,y ----
    beta_seq <- rep(0,n_pheno)
    beta_seq[which(gamma_seq == 1)] <- coeff_beta

    mu_seq <- rep(0,n_pheno)
    mu_seq[which(gamma_seq == 1)] <- coeff_mu
    mu_mtx <- matrix(rep(mu_seq,70),nrow = n_ind,ncol = n_pheno,byrow = T)

    eps_seq <- rnorm(n_pheno*n_ind,mean = 0,sd = param_sigma_beta)
    eps_mtx <- matrix(eps_seq,nrow = n_ind,ncol = n_pheno,byrow = T)

    g_seq <- rbinom(n_ind,size = 2,prob = param_gi_prob)
    cat(c("sim2","A","A",g_seq), file = paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/sim2.cis.geno"))
    beta_mtx <- g_seq %*% t(beta_seq)

    y_mtx <- mu_mtx + beta_mtx + eps_mtx
    write.table(y_mtx, file= paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sim2_WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
    cat(rep(1,n_pheno), file = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/use_all.txt"))
    cat(tying_grp, file = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sim_grouping.txt"))

    ### Step 3: Run through HMT ----
    command <- paste0("../../WaveQTL -gmode 1"
                      ," -group sim_grouping.txt -a ",sigma_prior
                      ," -d ",sigma_prior/4
                      ," -g ../../data/dsQTL/sim2.cis.geno -p sim2_WCs.txt -u use_all.txt -o sim2_noQT -f ",n_pheno," -hmt 1")
    system(command,show.output.on.console = showOutput)

    # I DON'T ACTUALLY HAVE A WAY OF RECOVERING PI_0 ATM UNLESS I RUN IT THROUGH WAVEQTL (NO HMT) ALSO!
    ### Step 4: Get HMT results ----
    file_prefix <- "sim2_noQT"
    pi_file <- exp(as.numeric(read.table(paste0("output/",file_prefix,".fph.pi.txt"))[2]))
    eps_11_file <- exp(as.numeric(read.table(paste0("output/",file_prefix,".fph.eps_11.txt"))[4:(n_pheno + 1)]))
    eps_10_file <- exp(as.numeric(read.table(paste0("output/",file_prefix,".fph.eps_10.txt"))[4:(n_pheno + 1)]))

    results_pi[[its]] <- pi_file
    results_eps_11[[its]] <- eps_11_file
    results_eps_10[[its]] <- eps_10_file
    results_gamma_seq[[its]] <- gamma_seq
    results_g_seq[[its]] <- g_seq
    results_y_mtx[[its]] <- y_mtx

  }

  setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
  return(list(results_pi = results_pi
              ,results_eps_11 = results_eps_11
              ,results_eps_10 = results_eps_10
              ,results_gamma_seq = results_gamma_seq
              ,results_g_seq = results_g_seq
              ,results_y_mtx = results_y_mtx
              ,param_eps_11 = param_eps_11
              ,param_eps_10 = param_eps_10
              ,param_pi_11 = param_pi_11
              ,param_pi_00 = param_pi_00))
}
