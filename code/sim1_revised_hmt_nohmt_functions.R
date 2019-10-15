
# Helper functions --------------------------------------------------------

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


# Simulation functions ----------------------------------------------------

draw_samples_no_hmt <- function(
  data_path = "~/Cpp/WaveQTL/test/dsQTL/output/"
  ,dataset = "test.no.QT"
  ,geno_select = 11
  ,wmat = Wmat_1024
  ,seed_no = 10
  ,num_samples = 1000
  ,num_pheno = 1024
  ,num_indiv = 70
){

  a_1 <- as.numeric(as.matrix(read.table(paste0(data_path,dataset,".fph.phi.txt")))[geno_select,2:1025])
  mean1 <- as.numeric(as.matrix(read.table(paste0(data_path,dataset,".fph.mean1.txt")))[geno_select,2:1025])
  var1 <- as.numeric(as.matrix(read.table(paste0(data_path,dataset,".fph.var1.txt")))[geno_select,2:1025])

  set.seed(seed_no)
  beta_data_samples <- matrix(nrow = num_samples,ncol = num_pheno)

  for(j in 1:num_samples){
    ## Sampling a bernoulli (1/0) from each phi:
    gamma_seq <- rbinom(n=num_pheno,size = 1,prob = a_1)

    ## back out a, b (from t-dist) parameters:
    t_nu <- num_indiv
    t_a <- mean1
    t_b <- var1*(t_nu-2)/t_nu

    t_sample <- stats::rt(n = num_pheno, df = t_nu)
    t_sample_3p <- t_a+(sqrt(t_b)*t_sample)

    ## Simulate beta, being either 3-param t-dist, or 0
    beta_seq <- rep(0,num_pheno)
    gamma_1_indx <- which(gamma_seq == 1)
    beta_seq[gamma_1_indx] <- t_sample_3p[gamma_1_indx]

    # Transform into data space
    beta_data_samples[j,] = as.vector(-matrix(data=beta_seq, nr = 1, nc = num_pheno)%*%as.matrix(Wmat_1024))
  }

  return(beta_data_samples)
}

draw_samples_hmt <- function(
  data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
  ,dataset = "tree_tie_noQT"
  ,waveqtl_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/WaveQTL/"
  ,waveqtl_dataset = "test.no.QT"
  ,geno_select = 11
  ,wmat = Wmat_1024
  ,seed_no = 10
  ,num_samples = 1000
  ,num_pheno = 1024
  ,num_indiv = 70
){

  a_1 <- as.matrix(read.table(paste0(data_path,dataset,".fph.pp.txt")))[geno_select,]
  b_11 <- as.matrix(read.table(paste0(data_path,dataset,".fph.pp_joint_11.txt")))[geno_select,]
  b_10 <- as.matrix(read.table(paste0(data_path,dataset,".fph.pp_joint_10.txt")))[geno_select,]

  # Just take the 1023 numeric values (excl first one as it's the scaling coefficient), from cols 3:1025. Also, we exp() our values as our software returned everything in logs.
  # Keep in mind that the first value of b_11, b_10 is just a placeholder - as the top element of the tree has no parent.
  a_1 <- exp(as.numeric(a_1[3:1025]))
  b_11 <- exp(as.numeric(b_11[3:1025]))
  b_10 <- exp(as.numeric(b_10[3:1025]))

  # Grab relevant quantities from WaveQTL (scaling coefficient)
  waveqtl_phi <- as.numeric(as.matrix(read.table(paste0(waveqtl_data_path,waveqtl_dataset,".fph.phi.txt")))[geno_select,2])

  # Load mean, var outputs from HMT
  mean1 <- as.numeric(as.matrix(read.table(paste0(data_path,dataset,".fph.mean1.txt")))[geno_select,2:1025])
  var1 <- as.numeric(as.matrix(read.table(paste0(data_path,dataset,".fph.var1.txt")))[geno_select,2:1025])

  set.seed(seed_no)
  beta_data_samples <- matrix(nrow = num_samples,ncol = num_pheno)

  for(j in 1:num_samples){
    # Generate gamma
    gamma_seq <- numeric()
    rand_seq <- runif(num_pheno)

    # Scaling coefficient
    gamma_seq[1] <- ifelse(rand_seq[1] < waveqtl_phi, 1, 0)

    # Head of tree
    gamma_seq[2] <- ifelse(rand_seq[2] < a_1[1], 1, 0)

    # i is the index of tree, where i = 1 is the head of the tree
    # Using this notation because that's how 'get_parent_indices' has been written
    for(i in 2:1023){
      indx <- i
      parent_indx <- get_parent_indices(indx)
      parent_gamma <- gamma_seq[parent_indx + 1]

      if(parent_gamma == 1){
        numerator <- b_11[indx]
        denominator <- a_1[parent_indx]
      }else if(parent_gamma == 0){
        numerator <- b_10[indx]
        denominator <- 1 - a_1[parent_indx]
      }

      post_prob <- numerator/denominator

      gamma_seq[i+1] <- ifelse(rand_seq[i+1] < post_prob, 1, 0)
    }

    # Generate beta
    t_nu <- num_indiv
    t_a <- mean1
    t_b <- var1*(t_nu-2)/t_nu

    t_sample <- stats::rt(n = num_pheno, df = t_nu)
    t_sample_3p <- t_a+(sqrt(t_b)*t_sample)

    ## Simulate beta, being either 3-param t-dist, or 0
    beta_seq <- rep(0,num_pheno)
    gamma_1_indx <- which(gamma_seq == 1)
    beta_seq[gamma_1_indx] <- t_sample_3p[gamma_1_indx]

    # Transform into data space
    beta_data_samples[j,] = as.vector(-matrix(data=beta_seq, nr = 1, nc = num_pheno)%*%as.matrix(Wmat_1024))
  }
  return(beta_data_samples)
}


# Closed form functions ---------------------------------------------------

closed_form_no_hmt <- function(
  data_path = "~/Cpp/WaveQTL/test/dsQTL/output/"
  ,dataset = "test.no.QT"
  ,geno_select = 11
  ,wmat = Wmat_1024
  ,wmat2 = W2mat_1024
  ,num_pheno = 1024
){
  ## Read posterior mean and varaince of effect sizes in Wavelet space
  ## Set a path to files
  beta_mean_path= "../test/dsQTL/output/test.no.QT.fph.mean.txt"
  beta_var_path = "../test/dsQTL/output/test.no.QT.fph.var.txt"

  ## Read posterior mean in Wavelet space and transform them back to data space
  beta_mean = as.numeric(read.table(paste0(data_path,dataset,".fph.mean.txt"))[geno_select,2:1025])
  beta_mean_dataS = as.vector(-matrix(data=beta_mean, nr = 1, nc = num_pheno)%*%as.matrix(wmat))

  ## Read posterior variance in Wavelet space, transform them back to data space, and get standard deviation
  beta_var = as.numeric(read.table(paste0(data_path,dataset,".fph.var.txt"))[geno_select,2:1025])
  beta_var_dataS = as.vector(matrix(data=beta_var, nr = 1, nc = num_pheno)%*%as.matrix(wmat2))
  beta_sd_dataS = sqrt(beta_var_dataS)

  return(list(
    beta_mean_dataS = beta_mean_dataS
    ,beta_sd_dataS = beta_sd_dataS
  ))
}

closed_form_hmt <- function(
  data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
  ,dataset = "tree_tie_noQT"
  ,waveqtl_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/WaveQTL/"
  ,waveqtl_dataset = "test.no.QT"
  ,geno_select = 11
  ,wmat = Wmat_1024
  ,seed_no = 10
  ,num_samples = 1000
  ,num_pheno = 1024
  ,num_indiv = 70
){

  # For scaling coeff, use phi from waveqtl
  scaling_coeff_phi <- read.table(paste0(waveqtl_data_path,waveqtl_dataset,".fph.phi.txt"))[geno_select,2]
  remaining_phi <- exp(as.numeric(as.matrix(read.table(paste0(data_path,dataset,".fph.pp.txt")))[geno_select,3:1025]))
  all_phi <- c(scaling_coeff_phi,remaining_phi)
  t_dist_means <- as.numeric(read.table(paste0(waveqtl_data_path,waveqtl_dataset,".fph.mean1.txt"))[geno_select,2:1025])
  t_dist_vars <- as.numeric(read.table(paste0(waveqtl_data_path,waveqtl_dataset,".fph.var1.txt"))[geno_select,2:1025])

  ## Means
  beta_mean <- all_phi*t_dist_means
  beta_mean_dataS <- as.vector(-matrix(data=beta_mean, nr = 1, nc = num_pheno)%*%as.matrix(wmat))

  # Just take the 1023 numeric values (excl first one as it's the scaling coefficient), from cols 3:1025. Also, we exp() our values as our software returned everything in logs.
  # Keep in mind that the first value of b_11, b_10 is just a placeholder - as the top element of the tree has no parent.
  a_1 <- exp(as.numeric(a_1[3:1025]))
  b_11 <- exp(as.numeric(b_11[3:1025]))
  b_10 <- exp(as.numeric(b_10[3:1025]))
  ## Read posterior mean and varaince of effect sizes in Wavelet space
  ## Set a path to files
  beta_mean_path= "../test/dsQTL/output/test.no.QT.fph.mean.txt"
  beta_var_path = "../test/dsQTL/output/test.no.QT.fph.var.txt"

  ## Read posterior mean in Wavelet space and transform them back to data space
  beta_mean = as.numeric(read.table(paste0(data_path,dataset,".fph.mean.txt"))[geno_select,2:1025])
  beta_mean_dataS = as.vector(-matrix(data=beta_mean, nr = 1, nc = num_pheno)%*%as.matrix(wmat))

  ## Read posterior variance in Wavelet space, transform them back to data space, and get standard deviation
  beta_var = as.numeric(read.table(paste0(data_path,dataset,".fph.var.txt"))[geno_select,2:1025])
  beta_var_dataS = as.vector(matrix(data=beta_var, nr = 1, nc = num_pheno)%*%as.matrix(wmat2))
  beta_sd_dataS = sqrt(beta_var_dataS)

  return(list(
    beta_mean_dataS = beta_mean_dataS
    ,beta_sd_dataS = beta_sd_dataS
  ))
}

