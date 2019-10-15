
# Preliminaries -----------------------------------------------------------

input_data_path = "~/Cpp/WaveQTL_HMT/data/dsQTL/"
data_path <- "~/Cpp/WaveQTL_HMT/test/dsQTL/output/" # change this when you port it all over to the Masters Git repo
dataset <- "tree_tie_noQT"
waveqtl_data_path <- "~/Cpp/WaveQTL_HMT/test/dsQTL/output/WaveQTL/"
waveqtl_dataset <- "test.no.QT"
geno_select <- 11 # the one used in the demo

library(rmutil) # for beta-binomial distribution
library(dplyr)
source("code/sim3_functions.R")

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

# Generate own SNPs data
# set.seed(10)
group_data <- rbinom(n = 70,size = 1,prob = 0.5)
write.table(t(c("blah","A","A",group_data)), file= paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/", "sim3.cis.geno"), row.names=FALSE, col.names = FALSE, quote=FALSE)

# Summarise pheno.data ----------------------------------------------------

# Count summation
seq_sum <- apply(pheno.dat,MARGIN = 2,sum)
# Count average
seq_avg <- apply(pheno.dat,MARGIN = 2,mean)

# WaveQTL effect size -----------------------------------------------------

waveqtl_hmt_geno11 <- with_hmt_effect_size(data_path = data_path
                                           ,dataset = dataset
                                           ,waveqtl_dataset = paste0("WaveQTL/",waveqtl_dataset)
                                           ,Wmat_1024 = Wmat_1024
                                           ,geno_select = 11
                                           ,plot_title = "Posterior mean +/3 posterior standard deviaion - SNP 11")


# Sim function ------------------------------------------------------------
### This is a production one
run_sim3 <- function(
  sequencing_sums = seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 50 # must be one of 8, 16, 32, 64
  , effect_interval
  , effect_size_data
  , use_qt_data = FALSE # Will turn QT on and off, in time. Haven't integrated yet
  , over_disp_param = 1/70
  , Wmat_1024 = Wmat_1024
  , W2mat_1024 = W2mat_1024
  , library.read.depth = library.read.depth
  , Covariates = Covariates
  , group_data = group_data
){

  # Convert effect size to ratio --------------------------------------------
  print("Generating effect size params...")
  effect_ratio <- 1 + (num_indivs*effect_size_data$beta_dataS/sequencing_sums)
  effect_ratio[sequencing_sums == 0] <- 1

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

  print("Generating null and alt datasets...")

  # Null
  null_data <- matrix(nrow = 70,ncol = 1024)
  for(j in 1:1024){
    null_data[,j] <- rmutil::rbetabinom(n = 70, size = as.numeric(as.vector(ceiling(sequencing_sums[j])))
                                        , m = 1/70
                                        , s = over_disp_param)
  }

  # Alt
  # For alt dataset, create a 70 X 1024 matrix, based on group membership, assigning p1 or p2, respectively
  param_mtx <- matrix(nrow = 70,ncol = 1024)
  n <- 1
  for(i in group_data){
    if(i == 1){
      param_mtx[n,] <- (p2_alpha/(p2_alpha+p2_beta))
    }else{
      param_mtx[n,] <- (p1_alpha/(p1_alpha+p1_beta))
    }
    n <- n + 1
  }

  alt_data <- matrix(nrow = 70,ncol = 1024)
  for(i in 1:70){
    for(j in 1:1024){
      alt_data[i,j] <- rmutil::rbetabinom(n = 1, size = as.numeric(as.vector(ceiling(sequencing_sums[j])))
                                          , m = param_mtx[i,j]
                                          , s = over_disp_param)
    }
  }
  # # Null
  # # Alternatively, do column-wise
  # null_data <- matrix(nrow = num_indivs,ncol = num_bases)
  # for(i in 1:num_bases){
  #   null_data[,i] <- rmutil::rbetabinom(n = num_indivs, size = ceiling(sequencing_sums[i])
  #                                          , m = (p1_alpha/(p1_alpha+p1_beta))[i]
  #                                          , s = (p1_alpha+p1_beta)[i])
  # }
  #
  # # Alt
  # # Params should be p2 throughout
  # alt_data <- matrix(nrow = num_indivs,ncol = num_bases)
  # for(i in 1:num_bases){
  #   alt_data[,i] <- rmutil::rbetabinom(n = num_indivs, size = ceiling(sequencing_sums[i])
  #                                         , m = (p2_alpha/(p2_alpha+p2_beta))[i]
  #                                         , s = (p2_alpha+p2_beta)[i])
  # }


  # Clean datasets ----------------------------------------------------------

  print("Cleaning datasets...")
  wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = null_data
                                           ,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/null_data/")
                                           ,library.read.depth = library.read.depth
                                           ,Covariates = Covariates)
  wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = alt_data
                                           ,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/alt_data/")
                                           ,library.read.depth = library.read.depth
                                           ,Covariates = Covariates)


  # Execute WaveQTL and WAveQTL_HMT scripts ---------------------------------

  if(!use_qt_data){
    setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
    #### Run null dataset
    print("Executing null datasets...")
    # No HMT
    system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/null_data/WCs.no.QT.txt -u sims/length_",effect_length,"/null_data/use.txt -o sim3_noQT_null -f ",num_bases," -fph 1")
           ,show.output.on.console = F)
    # HMT
    system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/null_data/WCs.no.QT.txt -u sims/length_",effect_length,"/null_data/use.txt -o sim3_noQT_null_HMT -f ",num_bases," -hmt 1")
           ,show.output.on.console = F)
    #### Run alt dataset
    print("Executing alt datasets...")
    # No HMT
    system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/alt_data/WCs.no.QT.txt -u sims/length_",effect_length,"/alt_data/use.txt -o sim3_noQT_alt -f ",num_bases," -fph 1")
           ,show.output.on.console = F)
    # HMT
    system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/alt_data/WCs.no.QT.txt -u sims/length_",effect_length,"/alt_data/use.txt -o sim3_noQT_alt_HMT -f ",num_bases," -hmt 1")
           ,show.output.on.console = F)
    setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
  }else{
    setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
    #### Run null dataset
    print("Executing null datasets...")
    # No HMT
    system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/null_data/WCs.txt -u sims/length_",effect_length,"/null_data/use.txt -o sim3_QT_null -f ",num_bases," -fph 1")
           ,show.output.on.console = F)
    # HMT
    system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/null_data/WCs.txt -u sims/length_",effect_length,"/null_data/use.txt -o sim3_QT_null_HMT -f ",num_bases," -hmt 1")
           ,show.output.on.console = F)
    #### Run alt dataset
    print("Executing alt datasets...")
    # No HMT
    system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/alt_data/WCs.txt -u sims/length_",effect_length,"/alt_data/use.txt -o sim3_QT_alt -f ",num_bases," -fph 1")
           ,show.output.on.console = F)
    # HMT
    system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/alt_data/WCs.txt -u sims/length_",effect_length,"/alt_data/use.txt -o sim3_QT_alt_HMT -f ",num_bases," -hmt 1")
           ,show.output.on.console = F)
    setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
  }


  # Analysis ----------------------------------------------------------------

  ### Analysis - no HMT
  print("Analysing effect size - no HMT...")
  ##### Null
  null_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
  null_data_prefix = "sim3_noQT_null"
  null <- no_hmt_effect_size(data_path = null_data_path
                             ,data_prefix = null_data_prefix
                             ,Wmat_1024 = Wmat_1024
                             ,W2mat_1024 = W2mat_1024
                             ,sel_geno_IX = 1
                             ,plot_title = "Posterior mean +/-3 posterior standard deviation - null")

  ##### Alt
  alt_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
  alt_data_prefix = "sim3_noQT_alt"
  alt <- no_hmt_effect_size(data_path = alt_data_path
                            ,data_prefix = alt_data_prefix
                            ,Wmat_1024 = Wmat_1024
                            ,W2mat_1024 = W2mat_1024
                            ,sel_geno_IX = 1
                            ,plot_title = "Posterior mean +/-3 posterior standard deviation - alt")

  ### Analysis - with HMT
  print("Analysing effect size - with HMT...")
  ##### Null
  null_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
  null_data_prefix = "sim3_noQT_null"
  null_hmt <- with_hmt_effect_size(data_path = null_data_path
                                   ,dataset = paste0(null_data_prefix,"_HMT")
                                   ,waveqtl_dataset = null_data_prefix
                                   ,Wmat_1024 = Wmat_1024
                                   ,geno_select = 1
                                   ,plot_title = "Posterior mean +/-3 posterior standard deviation - null - HMT")

  ##### Alt
  alt_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
  alt_data_prefix = "sim3_noQT_alt"
  alt_hmt <- with_hmt_effect_size(data_path = alt_data_path
                                  ,dataset = paste0(alt_data_prefix,"_HMT")
                                  ,waveqtl_dataset = alt_data_prefix
                                  ,Wmat_1024 = Wmat_1024
                                  ,geno_select = 1
                                  ,plot_title = "Posterior mean +/-3 posterior standard deviation - alt - HMT")


  # Should make the output aliases and directories customisable -- again, will do in time.
  return(list(
    effect_ratio = effect_ratio
    ,effect_interval = effect_interval
    ,p1_vector = p1_vector
    ,p2_vector = p2_vector
    ,null_data = null_data
    ,alt_data = alt_data
    ,null_analysis = null
    ,alt_analysis = alt
    ,null_hmt_analysis = null_hmt
    ,alt_hmt_analysis = alt_hmt
  ))

}


# Test sim function -------------------------------------------------------

### Has some ridiculous customisations, for TESTS only
run_sim3_test <- function(
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
  , group_data = group_data
  , effect_size = 10
  , effect_multiple = NULL
  , num_trials = NULL
  , trials_multiple = 1
){

  # Convert effect size to ratio --------------------------------------------
  print("Generating effect size params...")
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

  print("Generating null and alt datasets...")

  if(!is.null(num_trials)){
    # Use a generic number, uniform across all bases
    trials_vect <- rep(num_trials,num_bases)
  }else{
    # Use the sequencing counts
    trials_vect <- as.numeric(as.vector(ceiling(sequencing_sums)))
    # Scale them up, as desired
    trials_vect <- trials_vect*trials_multiple
  }

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


  # Plotting difference of datasets -----------------------------------------
  null_data_avg <- apply(null_data,2,sum)
  alt_data_avg <- apply(alt_data,2,sum)

  g0_mean_n <- apply(null_data[which(group_data == 0),],2,mean)
  g1_mean_n <- apply(null_data[which(group_data == 1),],2,mean)
  g0_mean_a <- apply(alt_data[which(group_data == 0),],2,mean)
  g1_mean_a <- apply(alt_data[which(group_data == 1),],2,mean)

  y_min <- min(min(g0_mean_n - g1_mean_n),min(g0_mean_a - g1_mean_a))
  y_max <- max(max(g0_mean_n - g1_mean_n),max(g0_mean_a - g1_mean_a))
  plt_rng_y_2 <- c(min(alt_data_avg - null_data_avg) * 0.9999999999999, max(alt_data_avg - null_data_avg) * 1.000000000000001)

  par(mfrow=c(1,1))
  plot(1,1,type="n"
       , xlab = "Base location"
       , ylab = "simulated avg counts"
       , ylim=c(y_min, y_max)
       , xlim=c(1, 1024)
       , main ="Simulated NULL data - g0 mean (no effect) minus g1 mean"
       , axes=FALSE)
  axis(2)
  axis(1, at = c(1,seq(128,1024,128)))
  if(length(effect_interval) > 0){
    for(j in 1:length(effect_interval)){
      polygon(c(effect_interval[j]-0.5, effect_interval[j]-0.5, effect_interval[j]+0.5, effect_interval[j]+0.5), c(plt_rng_y_2[1], plt_rng_y_2[2], plt_rng_y_2[1], plt_rng_y_2[2]), col ="pink", border = NA)
    }
  }
  lines(g0_mean_n - g1_mean_n, col = "red")
  # legend("topleft", legend=c("p1", "p2"),
  #        col=c("red", "green"), lty=c(1,1), cex=0.8,
  #        box.lty=0)
  box()
  null_effect <- recordPlot()

  par(mfrow=c(1,1))
  plot(1,1,type="n"
       , xlab = "Base location"
       , ylab = "simulated avg counts"
       , ylim=c(y_min, y_max)
       , xlim=c(1, 1024)
       , main ="Simulated ALT data - g0 (no effect) mean minus g1 mean"
       , axes=FALSE)
  axis(2)
  axis(1, at = c(1,seq(128,1024,128)))
  if(length(effect_interval) > 0){
    for(j in 1:length(effect_interval)){
      polygon(c(effect_interval[j]-0.5, effect_interval[j]-0.5, effect_interval[j]+0.5, effect_interval[j]+0.5), c(plt_rng_y_2[1], plt_rng_y_2[2], plt_rng_y_2[1], plt_rng_y_2[2]), col ="pink", border = NA)
    }
  }
  lines(g0_mean_a - g1_mean_a, col = "red")
  # legend("topleft", legend=c("p1", "p2"),
  #        col=c("red", "green"), lty=c(1,1), cex=0.8,
  #        box.lty=0)
  box()
  alt_effect <- recordPlot()


  # Clean datasets ----------------------------------------------------------

  print("Cleaning datasets...")
  # wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = null_data
  #                                          ,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/null_data/")
  #                                          ,library.read.depth = library.read.depth
  #                                          ,Covariates = Covariates)
  wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = alt_data
                                           ,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/alt_data/")
                                           ,library.read.depth = library.read.depth
                                           ,Covariates = Covariates
                                           ,no.QT = !use_qt_data)


  # Execute WaveQTL and WAveQTL_HMT scripts ---------------------------------

  if(!use_qt_data){
    setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
    # #### Run null dataset
    # print("Executing null datasets...")
    # # No HMT
    # system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/null_data/WCs.no.QT.txt -u sims/length_",effect_length,"/null_data/use.txt -o sim3_null -f ",num_bases," -fph 1")
    #        ,show.output.on.console = F)
    # # HMT
    # system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/null_data/WCs.no.QT.txt -u sims/length_",effect_length,"/null_data/use.txt -o sim3_null_HMT -f ",num_bases," -hmt 1")
    #        ,show.output.on.console = F)
    #### Run alt dataset
    print("Executing alt datasets...")
    # No HMT
    system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/alt_data/WCs.no.QT.txt -u sims/length_",effect_length,"/alt_data/use.txt -o sim3_alt -f ",num_bases," -fph 1")
           ,show.output.on.console = F)
    # HMT
    system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/alt_data/WCs.no.QT.txt -u sims/length_",effect_length,"/alt_data/use.txt -o sim3_alt_HMT -f ",num_bases," -hmt 1")
           ,show.output.on.console = F)
    setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
  }else{
    setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
    # #### Run null dataset
    # print("Executing null datasets...")
    # # No HMT
    # system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/null_data/WCs.txt -u sims/length_",effect_length,"/null_data/use.txt -o sim3_null -f ",num_bases," -fph 1")
    #        ,show.output.on.console = F)
    # # HMT
    # system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/null_data/WCs.txt -u sims/length_",effect_length,"/null_data/use.txt -o sim3_null_HMT -f ",num_bases," -hmt 1")
    #        ,show.output.on.console = F)
    #### Run alt dataset
    print("Executing alt datasets...")
    # No HMT
    system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/alt_data/WCs.txt -u sims/length_",effect_length,"/alt_data/use.txt -o sim3_alt -f ",num_bases," -fph 1")
           ,show.output.on.console = F)
    # HMT
    system(paste0("../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/sim3.cis.geno -p sims/length_",effect_length,"/alt_data/WCs.txt -u sims/length_",effect_length,"/alt_data/use.txt -o sim3_alt_HMT -f ",num_bases," -hmt 1")
           ,show.output.on.console = F)
    setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
  }


  # Analysis ----------------------------------------------------------------

  ### Analysis - no HMT
  print("Analysing effect size - no HMT...")
  # ##### Null
  # null_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
  # null_data_prefix = "sim3_null"
  # null <- no_hmt_effect_size(data_path = null_data_path
  #                            ,data_prefix = null_data_prefix
  #                            ,Wmat_1024 = Wmat_1024
  #                            ,W2mat_1024 = W2mat_1024
  #                            ,sel_geno_IX = 1
  #                            ,plot_title = "Posterior mean +/-3 posterior standard deviation - null")

  ##### Alt
  alt_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
  alt_data_prefix = "sim3_alt"
  alt <- no_hmt_effect_size(data_path = alt_data_path
                            ,data_prefix = alt_data_prefix
                            ,Wmat_1024 = Wmat_1024
                            ,W2mat_1024 = W2mat_1024
                            ,sel_geno_IX = 1
                            ,plot_title = "Posterior mean +/-3 posterior standard deviation - alt")

  ### Analysis - with HMT
  print("Analysing effect size - with HMT...")
  # ##### Null
  # null_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
  # # null_data_prefix = "sim3_noQT_null"
  # null_hmt <- with_hmt_effect_size(data_path = null_data_path
  #                                  ,dataset = paste0(null_data_prefix,"_HMT")
  #                                  ,waveqtl_dataset = null_data_prefix
  #                                  ,Wmat_1024 = Wmat_1024
  #                                  ,geno_select = 1
  #                                  ,plot_title = "Posterior mean +/-3 posterior standard deviation - null - HMT")

  ##### Alt
  alt_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
  # alt_data_prefix = "sim3_noQT_alt"
  alt_hmt <- with_hmt_effect_size(data_path = alt_data_path
                                  ,dataset = paste0(alt_data_prefix,"_HMT")
                                  ,waveqtl_dataset = alt_data_prefix
                                  ,Wmat_1024 = Wmat_1024
                                  ,geno_select = 1
                                  ,plot_title = "Posterior mean +/-3 posterior standard deviation - alt - HMT")


  # Zoomed plots ------------------------------------------------------------
  # Determine graph boundaries
  # No HMT
  sample_mean <- alt$beta_dataS
  sample_sd <- alt$beta_sd_dataS
  ymin_beta_noHMT = min(sample_mean - 3*sample_sd) - abs(min(sample_mean - 3*sample_sd))*0.0000000001
  ymax_beta_noHMT = max(sample_mean + 3*sample_sd) + abs(max(sample_mean + 3*sample_sd))*0.0000000001
  # HMT
  sample_mean <- alt_hmt$beta_dataS
  sample_sd <- alt_hmt$beta_sd_dataS
  ymin_beta_HMT = min(sample_mean - 3*sample_sd) - abs(min(sample_mean - 3*sample_sd))*0.0000000001
  ymax_beta_HMT = max(sample_mean + 3*sample_sd) + abs(max(sample_mean + 3*sample_sd))*0.0000000001

  ymin_beta = min(ymin_beta_HMT,ymin_beta_noHMT)
  ymax_beta = min(ymax_beta_HMT,ymax_beta_noHMT)

  # No HMT
  p_alt_zoom <- effect_size_plot(
    y_min = ymin_beta
    , y_max = ymax_beta
    , beta_mean = alt$beta_dataS
    , beta_sd = alt$beta_sd_dataS
    , x_range = c(effect_interval[1], rev(effect_interval)[1])
    , x_ticks = seq(effect_interval[1],rev(effect_interval)[1],by = 10)
    , plot_title = "Posterior mean +/-3 posterior standard deviation - alt - no-HMT")

  # HMT
  p_alt_hmt_zoom <- effect_size_plot(
    y_min = ymin_beta
    , y_max = ymax_beta
    , beta_mean = alt_hmt$beta_dataS
    , beta_sd = alt_hmt$beta_sd_dataS
    , x_range = c(effect_interval[1], rev(effect_interval)[1])
    , x_ticks = seq(effect_interval[1],rev(effect_interval)[1],by = 10)
    , plot_title = "Posterior mean +/-3 posterior standard deviation - alt - HMT")


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
    effect_ratio = effect_ratio
    ,effect_interval = effect_interval
    ,p1_vector = p1_vector
    ,p2_vector = p2_vector
    # ,null_data = null_data
    ,alt_data = alt_data
    # ,null_analysis = null
    ,alt_analysis = alt
    # ,null_hmt_analysis = null_hmt
    ,alt_hmt_analysis = alt_hmt
    ,null_effect_plot = null_effect
    ,alt_effect_plot = alt_effect
    ,p_alt_zoom = p_alt_zoom
    ,p_alt_hmt_zoom = p_alt_hmt_zoom
    ,params_list = params_list
  ))

}

