## Functions required for the paper simulations.

# These are functions which have been adapted (usually only very lightly) from
# the original thesis work in:
# - code/sim3_functions.R
# - code/sim4_functions.R


# Get parent of given index in tree ---------------------------------------

# Function which grabs parent indices for particular index on a tree, relative
# to the index that we denote to the tree's root.

# from 'code/sim2_script.R'
get_parent_indices <- function(indx, tree_root_indx = 1){
  return_indices <- indx
  return_indices <- (return_indices - tree_root_indx + 1) %/% 2

  # root of tree doesn't return any index
  return_indices[which(indx == tree_root_indx)] <- NA_integer_
  return(return_indices)
}

# Extract effect size in data space ---------------------------------------

# Given some simulated mean and SDs of the beta (effect size) and some plotting
# parameters, create a plot in the WaveQTL style.

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


# Generate HMT effect size ------------------------------------------------

# Function which runs a simulation to generate the effect size with HMT.
# Runs simulations, and also returns effect size plots.

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
  # set.seed(10)
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
  # set.seed(10)
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


# Read in data and generate effect size -----------------------------------
read_in_gen_eff_size_for_paper <- function(geno_select = 11
                                           ,pheno_data_file
                                           ,library_read_depth_file
                                           ,covariates_file
                                           ,wmat_matrix_file = "~/Cpp/WaveQTL/data/DWT/Wmat_1024"
                                           ,data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
                                           ,dataset = "tree_tie_noQT"
                                           ,waveqtl_dataset = "test.no.QT"){
  # Read in data ------------------------------------------------------------

  pheno.dat = as.matrix(read.table(pheno_data_file))
  dim(pheno.dat)
  #[1]   70 1024

  ### Is this useful at all? For later on
  ## read library read depth
  library.read.depth = scan(library_read_depth_file)
  length(library.read.depth)

  ## read Covariates
  Covariates = as.matrix(read.table(covariates_file))

  ## Read DWT matrix
  Wmat_1024 = read.table(wmat_matrix_file,as.is = TRUE)
  W2mat_1024 = Wmat_1024*Wmat_1024

  # ## Read in SNPs
  # geno_data = read.table("~/Cpp/WaveQTL/data/dsQTL/chr17.10160989.10162012.2kb.cis.geno",as.is = TRUE)
  # geno_data = geno_data[11,4:73]
  #
  # # Group based on midpoint of the data
  # med_data <- median(as.numeric(as.vector(geno_data[1,])))
  # group_data <- as.numeric(as.numeric(as.vector(geno_data[1,])) >= med_data)
  # write.table(t(c("blah","A","A",group_data)), file= paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/", "sim3.cis.geno"), row.names=FALSE, col.names = FALSE, quote=FALSE)

  # # Generate own SNPs data
  # # set.seed(10)
  # group_data <- rbinom(n = 70,size = 1,prob = 0.5)
  # write.table(t(c("blah","A","A",group_data)), file= paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/", "sim3.cis.geno"), row.names=FALSE, col.names = FALSE, quote=FALSE)

  # Summarise pheno.data ----------------------------------------------------

  # Count summation
  seq_sum <- apply(pheno.dat,MARGIN = 2,sum)
  # Count average
  seq_avg <- apply(pheno.dat,MARGIN = 2,mean)

  # WaveQTL effect size -----------------------------------------------------

  effect_size <- with_hmt_effect_size(data_path = data_path
                                      ,dataset = dataset
                                      ,waveqtl_dataset = waveqtl_dataset
                                      ,Wmat_1024 = Wmat_1024
                                      ,geno_select = geno_select
                                      ,plot_title = "Posterior mean +/3 posterior standard deviation")
  return(list(
    pheno.dat = pheno.dat
    ,library.read.depth = library.read.depth
    ,Covariates = Covariates
    ,Wmat_1024 = Wmat_1024
    ,W2mat_1024 = W2mat_1024
    # ,group_data = group_data
    ,seq_sum = seq_sum
    ,seq_avg = seq_avg
    ,effect_size = effect_size
  ))
}

# Wavelet cleaning wrapper script -----------------------------------------

# All the functionality of WaveQTL_preprocess_example.R from original WaveQTL repository.
wavelet_cleaning_wrapper_function_nonRMD <- function(
  pheno.dat
  , output.path
  , meanR.thresh = 2
  , library.read.depth
  , Covariates
  , no.QT = TRUE){

  ## read functions for WaveQTL preprocess
  ### Exactly the same file from the WaveQTL repository
  # source("~/WaveQTL/R/WaveQTL_preprocess_funcs.R")

  # data.path = "../code/WaveQTL/data/dsQTL/"

  ## set seed
  # set.seed(1)

  # ## read library read depth
  # library.read.depth = scan(paste0(data.path, "library.read.depth.dat"))
  #
  # ## read Covariates
  # Covariates = as.matrix(read.table(paste0(data.path, "PC4.dat")))

  res = WaveQTL_preprocess(Data = pheno.dat, library.read.depth=library.read.depth , Covariates = Covariates, meanR.thresh = meanR.thresh, no.QT)

  # No quantile transform
  if(no.QT){
    ## for effect size estimation, we need WCs without QT.
    write.table(res$WCs, file= paste0(output.path, "WCs.no.QT.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
  }else{
    # Else, WCs with QT
    write.table(res$WCs, file= paste0(output.path, "WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
  }

  cat(res$filtered.WCs, file = paste0(output.path, "use.txt"))

  ## produce group information and save it as a file
  group.info = generate_Group(dim(res$WCs)[2])
  cat(group.info, file = paste0(output.path, "group.txt"))

}



# Test simulation function ------------------------------------------------

### V2
### Amended to include the likelihood of the scaling coefficient

run_sim4_v2 <- function(
  sequencing_sums = seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 32 # must be one of 8, 16, 32, 64
  , effect_interval
  , effect_size_data
  , no.QT = FALSE # turn off/on QT data
  , over_disp_param = 70000
  , Wmat_1024 = Wmat_1024
  , W2mat_1024 = W2mat_1024
  , library.read.depth = library.read.depth
  , Covariates = Covariates
  # , group_data = group_data
  , effect_size = 10 # 'square' effect size, if effect_multiple null
  , effect_multiple = NULL # turn on/off the square effect size. NULL = off, number = data-driven effect size scaling constant
  , num_trials = NULL # If not null, uniform number of trials across all bases used. Else...
  , trials_multiple = 1 # if above is null, then scale up the data-driven number of trials by a constant
  , number_sims = 50
  , verbose = T
  , outputAlias # Output prefix for this simulation run
  , WaveQTL_directory # Directory where WaveQTL compilation lives
){

  # Exact same pre-process file as in WaveQTL repo
  source(paste0(WaveQTL_directory,"/R/WaveQTL_preprocess_funcs.R"))

  # Assign folders/paths for simulation -------------------------------------

  output_folder <- paste0(WaveQTL_directory, "/test/dsQTL/sims/",outputAlias,"/")
  if(!dir.exists(output_folder)){
    dir.create(output_folder, recursive = T)
  }

  ## Add in at start
  output_folder_null <- paste0(output_folder,"null/")
  output_folder_alt <- paste0(output_folder,"alt/")
  if(!dir.exists(output_folder_null)){
    dir.create(output_folder_null, recursive = T)
  }
  if(!dir.exists(output_folder_alt)){
    dir.create(output_folder_alt, recursive = T)
  }

  # Create SNP file ---------------------------------------------------------

  # Generate own SNPs data
  group_data <- rbinom(n = 70,size = 1,prob = 0.5)
  write.table(t(c("blah","A","A",group_data)), file= paste0(WaveQTL_directory, "/data/dsQTL/",outputAlias,".cis.geno"), row.names=FALSE, col.names = FALSE, quote=FALSE)

  # Convert effect size to ratio --------------------------------------------
  if(verbose){print("Generating effect size params...")}

  if(is.null(effect_multiple)){
    # If effect_multiple is NULL we use effect_size.
    # This generates the same effect size across all the selected bases.

    ### Here's a ridiculous effect size
    effect_ratio <- rep(effect_size, num_bases)
    effect_ratio[sequencing_sums == 0] <- 1
  }else{
    # If effect_multiple is NOT NULL we use effect_multiple.
    # This takes the absolute of the effect size detected from the data
    # and scales it up by a multiple (as per the thesis simulation methodology)

    ### Here's a sensible effect size
    #### HACK: currently doing absolute of effect size to prevent negatives...
    effect_ratio <- 1 + (num_indivs*abs(effect_size_data$beta_dataS)*effect_multiple/sequencing_sums)
    effect_ratio[sequencing_sums == 0] <- 1
  }

  # Convert ratio to beta-binomial param ------------------------------------

  p1_vector <- rep(1/num_indivs,num_bases)
  p2_vector <- rep(1/num_indivs,num_bases)

  # Add in effects, where required
  p1_vector[effect_interval] <- 2/num_indivs * (1/(1 + effect_ratio[effect_interval]))
  p2_vector[effect_interval] <- 2/num_indivs * (effect_ratio[effect_interval]/(1 + effect_ratio[effect_interval]))

  p1_alpha <- over_disp_param*p1_vector
  p1_beta <- over_disp_param - p1_alpha
  p2_alpha <- over_disp_param*p2_vector
  p2_beta <- over_disp_param - p2_alpha

  # Generate null and alt datasets ------------------------------------------

  if(verbose){print("Generating null and alt datasets...")}

  if(!is.null(num_trials)){
    # Use a generic number, uniform across all bases
    trials_vect <- rep(num_trials,num_bases)
  }else{
    # Use the sequencing counts
    trials_vect <- as.numeric(as.vector(ceiling(sequencing_sums)))
    # Scale them up, as desired
    trials_vect <- trials_vect*trials_multiple
  }


  # Simulate many times -----------------------------------------------------

  # Set up likelihood results vector
  null_waveqtl_lhood <- c()
  alt_waveqtl_lhood <- c()
  null_waveqtl_hmt_lhood <- c()
  alt_waveqtl_hmt_lhood <- c()

  for(lhood_vect_it in 1:number_sims){

    print(paste("Simulation:",lhood_vect_it))

    # set.seed(6)
    # Null
    null_data <- matrix(nrow = num_indivs,ncol = num_bases)
    for(i in 1:num_indivs){
      # Everyone uses same beta
      null_beta <- rbeta(n = num_bases
                         ,shape1 = over_disp_param*1/num_indivs
                         ,shape2 = over_disp_param - (over_disp_param*1/num_indivs))
      null_data[i,] <- rbinom(n = num_bases
                              , size = trials_vect
                              , prob = null_beta)
    }

    # Alt
    alt_data <- matrix(nrow = num_indivs,ncol = num_bases)
    for(i in 1:num_indivs){
      if(group_data[i] == 0){
        alt_beta <- rbeta(n = num_bases
                          ,shape1 = p1_alpha
                          ,shape2 = p1_beta)
      }else{
        alt_beta <- rbeta(n = num_bases
                          ,shape1 = p2_alpha
                          ,shape2 = p2_beta)
      }
      alt_data[i, ] <- rbinom(n = num_bases
                              ,size = trials_vect
                              ,prob = alt_beta)
    }

    # Clean datasets ----------------------------------------------------------

    if(verbose){print("Cleaning datasets...")}
    wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = null_data
                                             ,output.path = output_folder_null
                                             ,library.read.depth = library.read.depth
                                             ,Covariates = Covariates
                                             ,no.QT = no.QT)
    wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = alt_data
                                             ,output.path = output_folder_alt
                                             ,library.read.depth = library.read.depth
                                             ,Covariates = Covariates
                                             ,no.QT = no.QT)

    # Execute WaveQTL and WAveQTL_HMT scripts ---------------------------------

    if(no.QT){
      setwd(paste0(WaveQTL_directory,"/test/dsQTL/"))
      #### Run null dataset
      if(verbose){print("Executing null datasets...")}

      if(verbose){print("Executing null datasets...")}
      # No HMT
      # system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.no.QT.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null -f ",num_bases," -fph 1")
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.no.QT.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      # system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.no.QT.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null_HMT -f ",num_bases," -hmt 1")
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.no.QT.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)
      #### Run alt dataset
      if(verbose){print("Executing alt datasets...")}
      # No HMT
      # system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.no.QT.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt -f ",num_bases," -fph 1")
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.no.QT.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      # system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.no.QT.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt_HMT -f ",num_bases," -hmt 1")
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.no.QT.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)

    }else{
      setwd(paste0(WaveQTL_directory,"/test/dsQTL/"))
      #### Run null dataset
      if(verbose){print("Executing null datasets...")}
      # No HMT
      # system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null -f ",num_bases," -fph 1")
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      # system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null_HMT -f ",num_bases," -hmt 1")
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)
      #### Run alt dataset
      if(verbose){print("Executing alt datasets...")}
      # No HMT
      # system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt -f ",num_bases," -fph 1")
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      # system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt_HMT -f ",num_bases," -hmt 1")
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)

    }


    # Analysis ----------------------------------------------------------------
    # This part is different to last time.

    output_path = paste0(WaveQTL_directory, "/test/dsQTL/output/")
    null_data_prefix = paste0(outputAlias,"_null")
    alt_data_prefix = paste0(outputAlias,"_alt")

    # Note that output is in log10 space (keep it all in log 10 space for)

    ### No HMT
    # Null - WaveQTL
    null_waveqtl_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(output_path,null_data_prefix,".fph.logLR.txt"))[2])
    # Alt - WaveQTL
    alt_waveqtl_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(output_path,alt_data_prefix,".fph.logLR.txt"))[2])

    ### HMT
    # Null - WaveQTL_HMT
    # Grab scaling coefficient logBF, plus scaling coefficient pi. (from non-HMT analysis - will remain unchanged in HMT)
    sc_pi_null <- as.numeric(read.table(paste0(output_path,null_data_prefix,".fph.pi.txt"))[2])
    # Convert out of log10
    sc_BF_null <- 10^(as.numeric(read.table(paste0(output_path,null_data_prefix,".fph.logLR.txt"))[3]))
    # Calc lhood and convert into natural log (which is the base of the logLR)
    sc_logL_null <- log(sc_BF_null*sc_pi_null + (1-sc_pi_null))
    null_waveqtl_hmt_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(output_path,null_data_prefix,"_HMT.fph.logLR.txt"))[2])
    # Adjust with scaling coeff
    null_waveqtl_hmt_lhood[lhood_vect_it] <- null_waveqtl_hmt_lhood[lhood_vect_it] + sc_logL_null

    # Alt - WaveQTL_HMT
    # Grab scaling coefficient logBF, plus scaling coefficient pi. (from non-HMT analysis - will remain unchanged in HMT)
    sc_pi_alt <- as.numeric(read.table(paste0(output_path,alt_data_prefix,".fph.pi.txt"))[2])
    # Convert out of log10
    sc_BF_alt <- 10^(as.numeric(read.table(paste0(output_path,alt_data_prefix,".fph.logLR.txt"))[3]))
    # Calc lhood and convert into natural log (which is the base of the logLR)
    sc_logL_alt <- log(sc_BF_alt*sc_pi_alt + (1-sc_pi_alt))
    alt_waveqtl_hmt_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(output_path,alt_data_prefix,"_HMT.fph.logLR.txt"))[2])
    # Adjust with scaling coeff
    alt_waveqtl_hmt_lhood[lhood_vect_it] <- alt_waveqtl_hmt_lhood[lhood_vect_it] + sc_logL_alt
  }

  # List of parameters
  params_list <- list(
    num_indivs = num_indivs
    , num_bases = num_bases
    , effect_length = effect_length
    , effect_interval = effect_interval
    , no.QT = no.QT
    , over_disp_param = over_disp_param
    , group_data = group_data
    , effect_size = effect_size
    , effect_multiple = effect_multiple
    , num_trials = num_trials
    , trials_multiple = trials_multiple
  )


  # Should make the output aliases and directories customisable -- again, will do in time.
  return(list(
    # effect_ratio = effect_ratio
    # ,effect_interval = effect_interval
    # ,p1_vector = p1_vector
    # ,p2_vector = p2_vector
    # ,null_data = null_data
    # ,alt_data = alt_data
    # ,null_effect_plot = null_effect
    # ,alt_effect_plot = alt_effect
    null_waveqtl_lhood = null_waveqtl_lhood
    ,alt_waveqtl_lhood = alt_waveqtl_lhood
    ,null_waveqtl_hmt_lhood = null_waveqtl_hmt_lhood
    ,alt_waveqtl_hmt_lhood = alt_waveqtl_hmt_lhood
    ,params_list = params_list
  ))

}


