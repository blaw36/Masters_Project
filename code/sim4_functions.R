
# Preliminaries -----------------------------------------------------------

# input_data_path = "~/Cpp/WaveQTL_HMT/data/dsQTL/"
# data_path <- "~/Cpp/WaveQTL_HMT/test/dsQTL/output/" # change this when you port it all over to the Masters Git repo
# dataset <- "tree_tie_noQT"
# waveqtl_data_path <- "~/Cpp/WaveQTL_HMT/test/dsQTL/output/WaveQTL/"
# waveqtl_dataset <- "test.no.QT"
# # geno_select <- 11 # the one used in the demo

# library(rmutil) # for beta-binomial distribution
# library(dplyr)
# source("code/sim3_functions.R")


# Read in data and generate effect size -----------------------------------
read_in_gen_eff_size <- function(geno_select = 11
                                 ,input_data_path = "~/Cpp/WaveQTL_HMT/data/dsQTL/"
                                 ,data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
                                 ,dataset = "tree_tie_noQT"
                                 ,waveqtl_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/WaveQTL/"
                                 ,waveqtl_dataset = "test.no.QT"){
  # Read in data ------------------------------------------------------------

  pheno.dat = as.matrix(read.table(paste0(input_data_path, "chr17.10160989.10162012.pheno.dat")))
  dim(pheno.dat)
  #[1]   70 1024

  ### Is this useful at all? For later on
  ## read library read depth
  library.read.depth = scan(paste0(input_data_path, "library.read.depth.dat"))
  length(library.read.depth)

  ## read Covariates
  Covariates = as.matrix(read.table(paste0(input_data_path, "PC4.dat")))

  ## Read DWT matrix
  Wmat_1024 = read.table("~/Cpp/WaveQTL/data/DWT/Wmat_1024",as.is = TRUE)
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
                                      ,waveqtl_dataset = paste0("WaveQTL/",waveqtl_dataset)
                                      ,Wmat_1024 = Wmat_1024
                                      ,geno_select = geno_select
                                      ,plot_title = "Posterior mean +/3 posterior standard deviaion - SNP 11")
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



# Test sim function -------------------------------------------------------

run_sim4 <- function(
  sequencing_sums = seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 32 # must be one of 8, 16, 32, 64
  , effect_interval
  , effect_size_data
  , use_qt_data = FALSE # Will turn QT on and off, in time. Haven't integrated yet
  , over_disp_param = 70000
  , Wmat_1024 = Wmat_1024
  , W2mat_1024 = W2mat_1024
  , library.read.depth = library.read.depth
  , Covariates = Covariates
  # , group_data = group_data
  , effect_size = 10 # 'square' effect size, if effect_multiple null
  , effect_multiple = NULL # turn on/off the square effect size. NULL = off, number = data-driven effect size scaling constant
  , num_trials = NULL # uniform number of trials across all bases, if not null. else...
  , trials_multiple = 1 # if above is null, then scale up the data-driven number of trials by a constant
  , number_sims = 50
  , verbose = T
  , rMarkdownMode = T # TRUE if running in RMD, else False
  , outputAlias = "sim3"
){


  # Create SNP file ---------------------------------------------------------

  # Generate own SNPs data
  group_data <- rbinom(n = 70,size = 1,prob = 0.5)
  write.table(t(c("blah","A","A",group_data)), file= paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/",outputAlias,".cis.geno"), row.names=FALSE, col.names = FALSE, quote=FALSE)

  # Convert effect size to ratio --------------------------------------------
  if(verbose){print("Generating effect size params...")}

  if(is.null(effect_multiple)){
    ### Here's a ridulous effect size
    effect_ratio <- rep(effect_size, num_bases)
    effect_ratio[sequencing_sums == 0] <- 1
  }else{
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


    # # Plotting difference of datasets -----------------------------------------
    # null_data_avg <- apply(null_data,2,sum)
    # alt_data_avg <- apply(alt_data,2,sum)
    #
    # g0_mean_n <- apply(null_data[which(group_data == 0),],2,mean)
    # g1_mean_n <- apply(null_data[which(group_data == 1),],2,mean)
    # g0_mean_a <- apply(alt_data[which(group_data == 0),],2,mean)
    # g1_mean_a <- apply(alt_data[which(group_data == 1),],2,mean)
    #
    # y_min <- min(min(g0_mean_n - g1_mean_n),min(g0_mean_a - g1_mean_a))
    # y_max <- max(max(g0_mean_n - g1_mean_n),max(g0_mean_a - g1_mean_a))
    # plt_rng_y_2 <- c(min(alt_data_avg - null_data_avg) * 0.9999999999999, max(alt_data_avg - null_data_avg) * 1.000000000000001)
    #
    # pdf(NULL)
    # par(mfrow=c(1,1))
    # plot(1,1,type="n"
    #      , xlab = "Base location"
    #      , ylab = "simulated avg counts"
    #      , ylim=c(y_min, y_max)
    #      , xlim=c(1, 1024)
    #      , main ="Simulated NULL data - g0 mean (no effect) minus g1 mean"
    #      , axes=FALSE)
    # axis(2)
    # axis(1, at = c(1,seq(128,1024,128)))
    # if(length(effect_interval) > 0){
    #   for(j in 1:length(effect_interval)){
    #     polygon(c(effect_interval[j]-0.5, effect_interval[j]-0.5, effect_interval[j]+0.5, effect_interval[j]+0.5), c(plt_rng_y_2[1], plt_rng_y_2[2], plt_rng_y_2[1], plt_rng_y_2[2]), col ="pink", border = NA)
    #   }
    # }
    # lines(g0_mean_n - g1_mean_n, col = "red")
    # # legend("topleft", legend=c("p1", "p2"),
    # #        col=c("red", "green"), lty=c(1,1), cex=0.8,
    # #        box.lty=0)
    # box()
    # null_effect <- recordPlot()
    # invisible(dev.off())
    #
    # pdf(NULL)
    # par(mfrow=c(1,1))
    # plot(1,1,type="n"
    #      , xlab = "Base location"
    #      , ylab = "simulated avg counts"
    #      , ylim=c(y_min, y_max)
    #      , xlim=c(1, 1024)
    #      , main ="Simulated ALT data - g0 (no effect) mean minus g1 mean"
    #      , axes=FALSE)
    # axis(2)
    # axis(1, at = c(1,seq(128,1024,128)))
    # if(length(effect_interval) > 0){
    #   for(j in 1:length(effect_interval)){
    #     polygon(c(effect_interval[j]-0.5, effect_interval[j]-0.5, effect_interval[j]+0.5, effect_interval[j]+0.5), c(plt_rng_y_2[1], plt_rng_y_2[2], plt_rng_y_2[1], plt_rng_y_2[2]), col ="pink", border = NA)
    #   }
    # }
    # lines(g0_mean_a - g1_mean_a, col = "red")
    # # legend("topleft", legend=c("p1", "p2"),
    # #        col=c("red", "green"), lty=c(1,1), cex=0.8,
    # #        box.lty=0)
    # box()
    # alt_effect <- recordPlot()
    # invisible(dev.off())


    # Clean datasets ----------------------------------------------------------

    if(verbose){print("Cleaning datasets...")}

    if(rMarkdownMode){
      wavelet_cleaning_wrapper_function(pheno.dat = null_data
                                        ,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/null_data/")
                                        ,library.read.depth = library.read.depth
                                        ,Covariates = Covariates)
      wavelet_cleaning_wrapper_function(pheno.dat = alt_data
                                        ,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/alt_data/")
                                        ,library.read.depth = library.read.depth
                                        ,Covariates = Covariates)
    }else{
      wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = null_data
                                               ,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/null_data/")
                                               ,library.read.depth = library.read.depth
                                               ,Covariates = Covariates)
      wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = alt_data
                                               ,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/alt_data/")
                                               ,library.read.depth = library.read.depth
                                               ,Covariates = Covariates)
    }



    # Execute WaveQTL and WAveQTL_HMT scripts ---------------------------------

    if(!use_qt_data){
      setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
      #### Run null dataset
      if(verbose){print("Executing null datasets...")}
      # No HMT
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/length_",effect_length,"/null_data/WCs.no.QT.txt -u sims/length_",effect_length,"/null_data/use.txt -o ",outputAlias,"_null -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/length_",effect_length,"/null_data/WCs.no.QT.txt -u sims/length_",effect_length,"/null_data/use.txt -o ",outputAlias,"_null_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)
      #### Run alt dataset
      if(verbose){print("Executing alt datasets...")}
      # No HMT
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/length_",effect_length,"/alt_data/WCs.no.QT.txt -u sims/length_",effect_length,"/alt_data/use.txt -o ",outputAlias,"_alt -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/length_",effect_length,"/alt_data/WCs.no.QT.txt -u sims/length_",effect_length,"/alt_data/use.txt -o ",outputAlias,"_alt_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)

      if(rMarkdownMode){
        setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/analysis")
      }else{
        setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git")
      }

    }else{
      setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
      #### Run null dataset
      if(verbose){print("Executing null datasets...")}
      # No HMT
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/length_",effect_length,"/null_data/WCs.txt -u sims/length_",effect_length,"/null_data/use.txt -o ",outputAlias,"_null -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/length_",effect_length,"/null_data/WCs.txt -u sims/length_",effect_length,"/null_data/use.txt -o ",outputAlias,"_null_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)
      #### Run alt dataset
      if(verbose){print("Executing alt datasets...")}
      # No HMT
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/length_",effect_length,"/alt_data/WCs.txt -u sims/length_",effect_length,"/alt_data/use.txt -o ",outputAlias,"_alt -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/length_",effect_length,"/alt_data/WCs.txt -u sims/length_",effect_length,"/alt_data/use.txt -o ",outputAlias,"_alt_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)

      if(rMarkdownMode){
        setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/analysis")
      }else{
        setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git")
      }
    }


    # Analysis ----------------------------------------------------------------
    # This part is different to last time.

    null_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
    null_data_prefix = paste0(outputAlias,"_null")

    alt_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
    alt_data_prefix = paste0(outputAlias,"_alt")


    # Note that output is in log10 space (keep it all in log 10 space for)
    # Null - WaveQTL
    null_waveqtl_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(null_data_path,null_data_prefix,".fph.logLR.txt"))[2])
    # Alt - WaveQTL
    alt_waveqtl_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(alt_data_path,alt_data_prefix,".fph.logLR.txt"))[2])
    # Null - WaveQTL_HMT
    null_waveqtl_hmt_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(null_data_path,null_data_prefix,"_HMT.fph.logLR.txt"))[2])
    # Alt - WaveQTL_HMT
    alt_waveqtl_hmt_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(alt_data_path,alt_data_prefix,"_HMT.fph.logLR.txt"))[2])
  }

  # List of parameters
  params_list <- list(
    num_indivs = num_indivs
    , num_bases = num_bases
    , effect_length = effect_length
    , effect_interval = effect_interval
    , use_qt_data = use_qt_data
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


# Test sim function -------------------------------------------------------
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
  , num_trials = NULL # uniform number of trials across all bases, if not null. else...
  , trials_multiple = 1 # if above is null, then scale up the data-driven number of trials by a constant
  , number_sims = 50
  , verbose = T
  , rMarkdownMode = T # TRUE if running in RMD, else False
  , outputAlias = "sim3"
){

  # Assign folders/paths for simulation -------------------------------------

  output_folder <- paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/",outputAlias)
  if(!dir.exists(output_folder)){
    dir.create(output_folder)
  }

  ## Add in at start
  output_folder_null <- paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/",outputAlias,"/null/")
  output_folder_alt <- paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/",outputAlias,"/alt/")
  if(!dir.exists(output_folder_null)){
    dir.create(output_folder_null)
  }
  if(!dir.exists(output_folder_alt)){
    dir.create(output_folder_alt)
  }

  # Create SNP file ---------------------------------------------------------

  # Generate own SNPs data
  group_data <- rbinom(n = 70,size = 1,prob = 0.5)
  write.table(t(c("blah","A","A",group_data)), file= paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/",outputAlias,".cis.geno"), row.names=FALSE, col.names = FALSE, quote=FALSE)

  # Convert effect size to ratio --------------------------------------------
  if(verbose){print("Generating effect size params...")}

  if(is.null(effect_multiple)){
    ### Here's a ridulous effect size
    effect_ratio <- rep(effect_size, num_bases)
    effect_ratio[sequencing_sums == 0] <- 1
  }else{
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
    if(rMarkdownMode){
      wavelet_cleaning_wrapper_function(pheno.dat = null_data
                                        ,output.path = output_folder_null
                                        ,library.read.depth = library.read.depth
                                        ,Covariates = Covariates
                                        ,no.QT = no.QT)
      wavelet_cleaning_wrapper_function(pheno.dat = alt_data
                                        ,output.path = output_folder_alt
                                        ,library.read.depth = library.read.depth
                                        ,Covariates = Covariates
                                        ,no.QT = no.QT)
    }else{
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
    }

    # Execute WaveQTL and WAveQTL_HMT scripts ---------------------------------

    if(no.QT){
      setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
      #### Run null dataset
      if(verbose){print("Executing null datasets...")}

      if(verbose){print("Executing null datasets...")}
      # No HMT
      system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.no.QT.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.no.QT.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)
      #### Run alt dataset
      if(verbose){print("Executing alt datasets...")}
      # No HMT
      system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.no.QT.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.no.QT.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)

      if(rMarkdownMode){
        setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/analysis")
      }else{
        setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git")
      }

    }else{
      setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
      #### Run null dataset
      if(verbose){print("Executing null datasets...")}
      # No HMT
      system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/null/WCs.txt -u sims/",outputAlias,"/null/use.txt -o ",outputAlias,"_null_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)
      #### Run alt dataset
      if(verbose){print("Executing alt datasets...")}
      # No HMT
      system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt -f ",num_bases," -fph 1")
             ,show.output.on.console = F)
      # HMT
      system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/",outputAlias,".cis.geno -p sims/",outputAlias,"/alt/WCs.txt -u sims/",outputAlias,"/alt/use.txt -o ",outputAlias,"_alt_HMT -f ",num_bases," -hmt 1")
             ,show.output.on.console = F)

      if(rMarkdownMode){
        setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/analysis")
      }else{
        setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git")
      }

    }


    # Analysis ----------------------------------------------------------------
    # This part is different to last time.

    null_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
    null_data_prefix = paste0(outputAlias,"_null")

    alt_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
    alt_data_prefix = paste0(outputAlias,"_alt")


    # Note that output is in log10 space (keep it all in log 10 space for)

    ### No HMT
    # Null - WaveQTL
    null_waveqtl_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(null_data_path,null_data_prefix,".fph.logLR.txt"))[2])
    # Alt - WaveQTL
    alt_waveqtl_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(alt_data_path,alt_data_prefix,".fph.logLR.txt"))[2])

    ### HMT
    # Null - WaveQTL_HMT
    # Grab scaling coefficient logBF, plus scaling coefficient pi. (from non-HMT analysis - will remain unchanged in HMT)
    sc_pi_null <- as.numeric(read.table(paste0(null_data_path,null_data_prefix,".fph.pi.txt"))[2])
    # Convert out of log10
    sc_BF_null <- 10^(as.numeric(read.table(paste0(null_data_path,null_data_prefix,".fph.logLR.txt"))[3]))
    # Calc lhood and convert into natural log (which is the base of the logLR)
    sc_logL_null <- log(sc_BF_null*sc_pi_null + (1-sc_pi_null))
    null_waveqtl_hmt_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(null_data_path,null_data_prefix,"_HMT.fph.logLR.txt"))[2])
    # Adjust with scaling coeff
    null_waveqtl_hmt_lhood[lhood_vect_it] <- null_waveqtl_hmt_lhood[lhood_vect_it] + sc_logL_null

    # Alt - WaveQTL_HMT
    # Grab scaling coefficient logBF, plus scaling coefficient pi. (from non-HMT analysis - will remain unchanged in HMT)
    sc_pi_alt <- as.numeric(read.table(paste0(alt_data_path,alt_data_prefix,".fph.pi.txt"))[2])
    # Convert out of log10
    sc_BF_alt <- 10^(as.numeric(read.table(paste0(alt_data_path,alt_data_prefix,".fph.logLR.txt"))[3]))
    # Calc lhood and convert into natural log (which is the base of the logLR)
    sc_logL_alt <- log(sc_BF_alt*sc_pi_alt + (1-sc_pi_alt))
    alt_waveqtl_hmt_lhood[lhood_vect_it] <- as.numeric(read.table(paste0(alt_data_path,alt_data_prefix,"_HMT.fph.logLR.txt"))[2])
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


waveqtl_diags <- function(sims_list,num_sims, plot_indiv_curves = T){

  preds_nohmt <- matrix(c(sims_list$null_waveqtl_lhood
                        ,sims_list$alt_waveqtl_lhood
                        ,rep(0,num_sims),rep(1,num_sims))
                      ,nrow = 2*num_sims,ncol = 2,byrow = F)
  preds_hmt <- matrix(c(sims_list$null_waveqtl_hmt_lhood
                      ,sims_list$alt_waveqtl_hmt_lhood
                      ,rep(0,num_sims),rep(1,num_sims))
                    ,nrow = 2*num_sims,ncol = 2,byrow = F)

  # Plot ROC
  ## Example with ROCR
  perf_nohmt <- performance(prediction(preds_nohmt[,1],preds_nohmt[,2]),measure = "auc")
  perf_hmt <- performance(prediction(preds_hmt[,1],preds_hmt[,2]),measure = "auc")

  if(plot_indiv_curves){
    plot(performance(prediction(preds_nohmt[,1],preds_nohmt[,2]),"tpr","fpr"),colorize=TRUE)
    plot(performance(prediction(preds_hmt[,1],preds_hmt[,2]),"tpr","fpr"),colorize=TRUE)
  }

  ## DIY
  n = dim(preds_hmt)[1]

  O1 = order(preds_nohmt[,1], decreasing =TRUE)
  O2 = order(preds_hmt[,1], decreasing =TRUE)

  C1  = c(0,cumsum(preds_nohmt[O1,2])) / sum(preds_nohmt[,2]) #True positives proportions or sensitivity
  C2  = c(0,cumsum(preds_hmt[O2,2])) / sum(preds_hmt[,2])

  FP1 = c(0,cumsum(1-preds_nohmt[O1,2])) / (n-sum(preds_nohmt[,2])) #false positives proportions based on Model 1.
  FP2 = c(0,cumsum(1-preds_hmt[O2,2])) / (n-sum(preds_hmt[,2]))

  plot(FP1,C1, type="l", lwd=1, col="black", xlab="False Positive rate", ylab="True Positive Rate")
  lines(FP2,C2, type="l", lwd=1, col="red")
  legend("bottomright",
         c(paste0("WaveQTL: auc - ",round(attr(perf_nohmt,"y.value")[[1]],4))
           ,paste0("WaveQTL_HMT: auc - ",round(attr(perf_hmt,"y.value")[[1]],4))),
         cex = 0.8,
         lty=c(1,1),
         lwd=c(1,1),
         col=c("black","red"))
  roc_plot <- recordPlot()

  return(list(
    perf_nohmt = perf_nohmt
    ,perf_hmt = perf_hmt
    ,roc_plot = roc_plot
  ))

}
