# rm(list = ls());gc();cat("\014");

run_sim2 <- function(
  n_ind = 70
  ,n_pheno = 1024
  ,tying_grp = c(1,2,3,5,9,17,33,65,129,257,513)
  ,param_pi_00 = 0
  ,param_pi_11 = 1
  ,grped_eps_11
  ,grped_eps_10
  ,coeff_mu = 0
  ,coeff_beta = 2
  ,param_gi_prob = 0.4
  ,param_sigma_beta = 0.5
  ,num_sims = 100
  ,seed = 20){

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

  # Round to nearest 100th, for simplicity
  num_tying_grps <- length(tying_grp)


  results_pi <- list()
  results_eps_11 <- list()
  results_eps_10 <- list()

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

    eps_seq <- rnorm(n_pheno*n_ind,mean = 0,sd = sqrt(param_sigma_beta))
    eps_mtx <- matrix(eps_seq,nrow = n_ind,ncol = n_pheno,byrow = T)

    g_seq <- rbinom(n_ind,size = 2,prob = param_gi_prob)
    cat(c("sim2","A","A",g_seq), file = paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/sim2.cis.geno"))
    beta_mtx <- g_seq %*% t(beta_seq)

    y_mtx <- mu_mtx + beta_mtx + eps_mtx
    write.table(y_mtx, file= paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sim2_WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
    cat(rep(1,n_pheno), file = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/use_all.txt"))


    ### Step 3: Run through HMT ----
    command <- paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim2.cis.geno -p sim2_WCs.txt -u use_all.txt -o sim2_noQT -f ",n_pheno," -hmt 1")
    system(command,show.output.on.console = F)
    # system("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim2.cis.geno -p sim2_WCs.txt -u use_all.txt -o sim2_noQT -f 1024 -hmt 1",show.output.on.console = F)

    # I DON'T ACTUALLY HAVE A WAY OF RECOVERING PI_0 ATM UNLESS I RUN IT THROUGH WAVEQTL (NO HMT) ALSO!
    ### Step 4: Get HMT results ----
    file_prefix <- "sim2_noQT"
    pi_file <- exp(as.numeric(read.table(paste0("output/",file_prefix,".fph.pi.txt"))[2]))
    eps_11_file <- exp(as.numeric(read.table(paste0("output/",file_prefix,".fph.eps_11.txt"))[4:(n_pheno + 1)]))
    eps_10_file <- exp(as.numeric(read.table(paste0("output/",file_prefix,".fph.eps_10.txt"))[4:(n_pheno + 1)]))

    results_pi[[its]] <- pi_file
    results_eps_11[[its]] <- eps_11_file
    results_eps_10[[its]] <- eps_10_file

  }

  setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
  return(list(results_pi = results_pi
              ,results_eps_11 = results_eps_11
              ,results_eps_10 = results_eps_10
              ,param_eps_11 = param_eps_11
              ,param_eps_10 = param_eps_10))

}

# saveRDS(results_pi,file = "sim2_pi.RDS",compress = T)
# saveRDS(results_eps_11,file = "sim2_eps_11.RDS",compress = T)
# saveRDS(results_eps_10,file = "sim2_eps_10.RDS",compress = T)
