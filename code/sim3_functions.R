## Supporting functions for sim3 family of scripts
## Sim3 is the family of scripts which does the 3rd simulation.
## The 3rd simulation was proving HMT's efficacy over non-HMT


# get_parent_indices ------------------------------------------------------

# from 'code/sim2_script.R'
get_parent_indices <- function(indx, tree_root_indx = 1){
  return_indices <- indx
  return_indices <- (return_indices - tree_root_indx + 1) %/% 2

  # root of tree doesn't return any index
  return_indices[which(indx == tree_root_indx)] <- NA_integer_
  return(return_indices)
}


# tree_plot ---------------------------------------------------------------

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


# WC cleaning -------------------------------------------------------------
# The WC cleaning functionality, functionalised up. See 'code/WaveQTL/R/WaveQTL_preprocess_example.R' for details of where all this came from.

# pheno.dat - the simulated 70 x 1024 matrix from above
# output.path - where you want all the data to be saved
# meanR.thresh - low WC filtering threshold (default is 2)
# library.read.depth - should be read in (see WaveQTL data)
# Covaraites - should be read in (see WaveQTL data)
wavelet_cleaning_wrapper_function <- function(pheno.dat, output.path, meanR.thresh = 2, library.read.depth, Covariates){
  ## read functions for WaveQTL preprocess
  source("../code/WaveQTL/R/WaveQTL_preprocess_funcs.R")

  # data.path = "../code/WaveQTL/data/dsQTL/"

  ## set seed
  set.seed(1)

  # ## read library read depth
  # library.read.depth = scan(paste0(data.path, "library.read.depth.dat"))
  #
  # ## read Covariates
  # Covariates = as.matrix(read.table(paste0(data.path, "PC4.dat")))

  ## preprocess functional phenotype
  res = WaveQTL_preprocess(Data = pheno.dat, library.read.depth=library.read.depth , Covariates = Covariates, meanR.thresh = meanR.thresh)

  ## save output as files
  cat(res$filtered.WCs, file = paste0(output.path, "use.txt"))
  write.table(res$WCs, file= paste0(output.path, "WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)

  ## produce group information and save it as a file
  group.info = generate_Group(dim(res$WCs)[2])
  cat(group.info, file = paste0(output.path, "group.txt"))

  ## for effect size estimation, we need WCs without QT.
  set.seed(1)
  res.noQT = WaveQTL_preprocess(Data = pheno.dat, library.read.depth=library.read.depth , Covariates = Covariates, meanR.thresh = meanR.thresh, no.QT = TRUE)

  write.table(res.noQT$WCs, file= paste0(output.path, "WCs.no.QT.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
}


# Literally the same thing as above, just the source is adapted for a non-RMD script
wavelet_cleaning_wrapper_function_nonRMD <- function(pheno.dat, output.path, meanR.thresh = 2, library.read.depth, Covariates){
  ## read functions for WaveQTL preprocess
  source("code/WaveQTL/R/WaveQTL_preprocess_funcs.R")

  # data.path = "../code/WaveQTL/data/dsQTL/"

  ## set seed
  set.seed(1)

  # ## read library read depth
  # library.read.depth = scan(paste0(data.path, "library.read.depth.dat"))
  #
  # ## read Covariates
  # Covariates = as.matrix(read.table(paste0(data.path, "PC4.dat")))

  ## preprocess functional phenotype
  res = WaveQTL_preprocess(Data = pheno.dat, library.read.depth=library.read.depth , Covariates = Covariates, meanR.thresh = meanR.thresh)

  ## save output as files
  cat(res$filtered.WCs, file = paste0(output.path, "use.txt"))
  write.table(res$WCs, file= paste0(output.path, "WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)

  ## produce group information and save it as a file
  group.info = generate_Group(dim(res$WCs)[2])
  cat(group.info, file = paste0(output.path, "group.txt"))

  ## for effect size estimation, we need WCs without QT.
  set.seed(1)
  res.noQT = WaveQTL_preprocess(Data = pheno.dat, library.read.depth=library.read.depth , Covariates = Covariates, meanR.thresh = meanR.thresh, no.QT = TRUE)

  write.table(res.noQT$WCs, file= paste0(output.path, "WCs.no.QT.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
}


# Select effect interval --------------------------------------------------

# effect_vect is the 'col_posi' output from the waveqtl effect size analysis
summarise_effect_intervals <- function(effect_vect){
  require(dplyr)
  convert_to_buckets <- cumsum(diff(effect_vect) != 1)

  # Add first one back in (is ommitted by the differencing operation)
  convert_to_buckets <- c(convert_to_buckets[1],convert_to_buckets)
  # Start from 1
  convert_to_buckets <- convert_to_buckets + 1

  buckets_and_posis <- data.frame(buckets = convert_to_buckets, pos = effect_vect)
  summary_table <- buckets_and_posis %>%
    group_by(buckets) %>%
    summarise(min_pos = min(pos), max_pos = max(pos), length = max(pos) - min(pos) + 1)

  return(summary_table)
}

# based on the summary of interval lengths (above), picks a compatibly sized effect
effect_length_picker <- function(interval_summary_table, length_reqd){
  # Pick a bucket which is of appropriate length
  candidate_buckets <- interval_summary_table %>% filter(length >= length_reqd) %>% pull(buckets)

  if(length(candidate_buckets) == 0){
    # If none satisfy this, just pick any bucket and trim accordingly (and display a warning)
    warning("None of the intervals in this effect size data are long enough for desired length size. Picking an interval from data space, but will be shorter than required. Picking longest remaining bucket option.")
    bucket <- which.max(interval_summary_table$length)
    # Just pick the whole interval for bucket picked
    picked_bucket <- interval_summary_table %>% filter(buckets == bucket)
    # Output effect interval
    effect_interval <- seq(picked_bucket$min_pos,picked_bucket$max_pos, by = 1)
    return(effect_interval)
  }else{
    bucket <- sample(candidate_buckets,1)
  }

  # From that bucket, pick a random region of length required.
  # Equivalent to, pick a starting point that is at least <length required> distance from the end
  picked_bucket <- interval_summary_table %>% filter(buckets == bucket)

  # Pick starting point
  pts_to_sample <- seq(picked_bucket$min_pos,picked_bucket$max_pos - length_reqd + 1)
  if(length(pts_to_sample) > 1){
    start_pt <- sample(x = pts_to_sample, size = 1)
  }else{
    start_pt <- pts_to_sample
  }

  # Output effect interval
  effect_interval <- seq(start_pt,start_pt + length_reqd - 1, by = 1)
  return(effect_interval)
}

# Extract effect size in data space ---------------------------------------

effect_size_plot <- function(y_min, y_max, beta_mean, beta_sd, x_range = c(1,1024)
                             , x_ticks = c(1,seq(128,1024,by = 128))
                             , plot_title = "Posterior mean +/-3 posterior standard deviation"){
  # ## Visualize estimated effect size in the data space
  # ymin_beta = min(beta_dataS - 3*beta_sd_dataS) - abs(min(beta_dataS - 3*beta_sd_dataS))*0.0000000001
  # ymax_beta = max(beta_dataS + 3*beta_sd_dataS) + abs(max(beta_dataS + 3*beta_sd_dataS))*0.0000000001

  beta_mean <- beta_mean[x_range[1]:x_range[2]]
  beta_sd <- beta_sd[x_range[1]:x_range[2]]

  beta_l = beta_mean - 3*beta_sd
  beta_r = beta_mean + 3*beta_sd

  wh_l = which(beta_l > 0)
  wh_r = which(beta_r < 0)
  high_wh = sort(unique(union(wh_l, wh_r)))

  xval = x_range[1]:x_range[2]
  col_posi = xval[high_wh]

  # pdf("../test/dsQTL/effectSize.pdf", width = 8, height=3)
  par(mar = c(2,4,4,2))
  plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(y_min, y_max),xlim=x_range
       ,main = plot_title, axes=FALSE)
  axis(2)
  axis(1, at = x_ticks)
  if(length(col_posi) > 0){
    for(j in 1:length(col_posi)){
      polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(y_min-2, y_max+2, y_max+2, y_min-2), col ="pink", border = NA)
    }
  }

  abline(h = 0, col = "red")
  points(xval, beta_mean, col = "blue", type="l")
  points(xval, beta_l, col = "skyblue", type="l")
  points(xval, beta_r, col = "skyblue", type="l")
  box()

  # dev.off()
  p <- recordPlot()
  return(list(p = p
              ,col_posi = col_posi)
  )
}

# data_path - where data is saved
# data_prefix - prefix used to save the output as per WaveQTL algo
# Wmat_1024 - the IDWT matrix
# W2mat_1024 - the IDWT matrix squared, ie. Wmat_1024 * Wmat_1024
# sel_geno_IX - effect size of which element in the log? (usually corresponds to which SNP in the various SNPs)
no_hmt_effect_size <- function(data_path, data_prefix, Wmat_1024, W2mat_1024, sel_geno_IX = 1
                               ,plot_title = "Posterior mean +/-3 posterior standard deviation"){

  ## Read posterior mean in Wavelet space and transform them back to data space
  beta_mean = as.numeric(read.table(paste0(data_path,data_prefix,".fph.mean.txt"))[sel_geno_IX,2:1025])
  beta_dataS = as.vector(-matrix(data=beta_mean, nr = 1, nc = 1024)%*%as.matrix(Wmat_1024))

  ## Read posterior variance in Wavelet space, transform them back to data space, and get standard deviation
  beta_var = as.numeric(read.table(paste0(data_path, data_prefix, ".fph.var.txt"))[sel_geno_IX,2:1025])
  beta_var_dataS = as.vector(matrix(data=beta_var, nr = 1, nc = 1024)%*%as.matrix(W2mat_1024))
  beta_sd_dataS = sqrt(beta_var_dataS)

  ## Visualize estimated effect size in the data space
  ymin_beta = min(beta_dataS - 3*beta_sd_dataS) - abs(min(beta_dataS - 3*beta_sd_dataS))*0.0000000001
  ymax_beta = max(beta_dataS + 3*beta_sd_dataS) + abs(max(beta_dataS + 3*beta_sd_dataS))*0.0000000001

  p <- effect_size_plot(y_min = ymin_beta
                        , y_max = ymax_beta
                        , beta_mean = beta_dataS
                        , beta_sd = beta_sd_dataS
                        , x_range = c(1,1024)
                        , plot_title = plot_title)

  return(list(
    beta_dataS = beta_dataS
    ,beta_sd_dataS = beta_sd_dataS
    ,p = p$p
    ,col_posi = p$col_posi
  ))
}

# data_path - where data is saved
# dataset - prefix used to save the output as per WaveQTL algo
# waveqtl_dataset - the prefix of the outputs from WaveQTL analysis (no HMT) on the same thing - needed for the scaling coef.
# Wmat_1024 - the IDWT matrix
# geno_select - effect size of which element in the log? (usually corresponds to which SNP in the various SNPs)
with_hmt_effect_size <- function(data_path, dataset, waveqtl_dataset, Wmat_1024, geno_select = 1
                                 ,num_sims = 1000, plot_title = "Posterior mean +/-3 posterior standard deviation"){

  a_1 <- as.matrix(read.table(paste0(data_path,dataset,".fph.pp.txt")))[geno_select,]
  b_11 <- as.matrix(read.table(paste0(data_path,dataset,".fph.pp_joint_11.txt")))[geno_select,]
  b_10 <- as.matrix(read.table(paste0(data_path,dataset,".fph.pp_joint_10.txt")))[geno_select,]
  # dim(a_1);dim(b_11);dim(b_10);
  # Just take the 1023 numeric values (excl first one as it's the scaling coefficient), from cols 3:1025. Also, we exp() our values as our software returned everything in logs.
  # Keep in mind that the first value of b_11, b_10 is just a placeholder - as the top element of the tree has no parent.
  a_1 <- exp(as.numeric(a_1[3:1025]))
  b_11 <- exp(as.numeric(b_11[3:1025]))
  b_10 <- exp(as.numeric(b_10[3:1025]))

  waveqtl_phi <- as.matrix(read.table(paste0(data_path,waveqtl_dataset,".fph.phi.txt")))[geno_select,]
  waveqtl_phi <- as.numeric(waveqtl_phi[2])

  # Now, let's simulate a sequence of gammas:
  gamma_seq <- numeric()
  post_prob_seq <- numeric()
  set.seed(10)
  rand_seq <- runif(1024)

  # Scaling coefficient
  gamma_seq[1] <- ifelse(rand_seq[1] < waveqtl_phi, 1, 0)
  post_prob_seq[1] <- waveqtl_phi

  # Head of tree
  gamma_seq[2] <- ifelse(rand_seq[2] < a_1[1], 1, 0)
  post_prob_seq[2] <- a_1[1]

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
    post_prob_seq[i+1] <- post_prob

    gamma_seq[i+1] <- ifelse(rand_seq[i+1] < post_prob, 1, 0)

  }

  # Load mean, var outputs from HMT
  mean1 <- as.matrix(read.table(paste0(data_path,dataset,".fph.mean1.txt")))[geno_select,]
  var1 <- as.matrix(read.table(paste0(data_path,dataset,".fph.var1.txt")))[geno_select,]

  # Load mean, var outputs from WaveQTL
  mean1_waveqtl <- as.matrix(read.table(paste0(data_path,waveqtl_dataset,".fph.mean1.txt")))[geno_select,2]
  var1_waveqtl <- as.matrix(read.table(paste0(data_path,waveqtl_dataset,".fph.var1.txt")))[geno_select,2]

  # Append top coeff to the 1023 numeric values (excl first one as it's the scaling coefficient), from cols 3:1025 from HMT
  mean1 <- c(as.numeric(mean1_waveqtl), as.numeric(mean1[3:1025]))
  var1 <- c(as.numeric(var1_waveqtl), as.numeric(var1[3:1025]))

  ## back out a, b (from t-dist) parameters:
  t_nu <- 70
  t_a <- mean1
  t_b <- var1*(t_nu-2)/t_nu
  num_pheno <- length(mean1)

  t_sample <- stats::rt(n = num_pheno, df = t_nu)
  t_sample_3p <- t_a+(sqrt(t_b)*t_sample)

  ## Simulate beta, being either 3-param t-dist, or 0
  beta_seq <- rep(0,num_pheno)
  gamma_1_indx <- which(gamma_seq == 1)
  beta_seq[gamma_1_indx] <- t_sample_3p[gamma_1_indx]

  ### '-ve' is taken to represent biological definition of effects (akin to a base level)
  beta_dataS = as.vector(-matrix(data=beta_seq, nr = 1, nc = 1024)%*%as.matrix(Wmat_1024))
  # plot(beta_dataS, main = "Simulation 1", type = "l")
  # abline(h = 0, col = "red")

  ### Now, run simulations
  num_samples <- num_sims
  set.seed(10)
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
    t_sample <- stats::rt(n = num_pheno, df = t_nu)
    t_sample_3p <- t_a+(sqrt(t_b)*t_sample)

    ## Simulate beta, being either 3-param t-dist, or 0
    beta_seq <- rep(0,num_pheno)
    gamma_1_indx <- which(gamma_seq == 1)
    beta_seq[gamma_1_indx] <- t_sample_3p[gamma_1_indx]

    # Transform into data space
    beta_data_samples[j,] = as.vector(-matrix(data=beta_seq, nr = 1, nc = num_pheno)%*%as.matrix(Wmat_1024))
  }

  sample_mean <- apply(beta_data_samples,MARGIN = 2,mean)
  sample_sd <- apply(beta_data_samples,MARGIN = 2,sd)

  ymin_beta = min(sample_mean - 3*sample_sd) - abs(min(sample_mean - 3*sample_sd))*0.0000000001
  ymax_beta = max(sample_mean + 3*sample_sd) + abs(max(sample_mean + 3*sample_sd))*0.0000000001

  p <- effect_size_plot(y_min = ymin_beta
                        , y_max = ymax_beta
                        , beta_mean = sample_mean
                        , beta_sd = sample_sd
                        , x_range = c(1,1024)
                        , plot_title = plot_title)

  return(list(
    beta_dataS = sample_mean
    ,beta_sd_dataS = sample_sd
    ,p = p$p
    ,col_posi = p$col_posi
  ))

}
