
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

run_sim3 <- function(
  sequencing_sums = seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 50 # must be one of 5, 10, 50
  , effect_size_data
  , use_qt_data = FALSE # Will turn QT on and off, in time. Haven't integrated yet
  , over_disp_param = 70
  , Wmat_1024 = Wmat_1024
  , W2mat_1024 = W2mat_1024
  , library.read.depth = library.read.depth
  , Covariates = Covariates
){

  # Pick effect bucket ------------------------------------------------------
  print("Picking effect bucket...\n")
  int_table <- summarise_effect_intervals(effect_size_data$col_posi)

  set.seed(10)
  effect_interval <- effect_length_picker(int_table,effect_length)


  # Convert effect size to ratio --------------------------------------------
  print("Generating effect size params...\n")
  effect_ratio <- 1+(num_indivs*effect_size_data$beta_dataS/(seq_sum + 1))


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

  print("Generating null and alt datasets...\n")
  # Null
  # Alternatively, do column-wise
  null_data <- matrix(nrow = num_indivs,ncol = num_bases)
  for(i in 1:num_bases){
    null_data[,i] <- rmutil::rbetabinom(n = num_indivs, size = ceiling(seq_sum[i])
                                           , m = (p1_alpha/(p1_alpha+p1_beta))[i]
                                           , s = (p1_alpha+p1_beta)[i])
  }

  # Alt
  # Params should be p2 throughout
  alt_data <- matrix(nrow = num_indivs,ncol = num_bases)
  for(i in 1:num_bases){
    alt_data[,i] <- rmutil::rbetabinom(n = num_indivs, size = ceiling(seq_sum[i])
                                          , m = (p2_alpha/(p2_alpha+p2_beta))[i]
                                          , s = (p2_alpha+p2_beta)[i])
  }


  # Clean datasets ----------------------------------------------------------

  print("Cleaning datasets...\n")
  wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = null_data
                                           ,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/null_data/")
                                           ,library.read.depth = library.read.depth
                                           ,Covariates = Covariates)
  wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = alt_data
                                           ,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/alt_data/")
                                           ,library.read.depth = library.read.depth
                                           ,Covariates = Covariates)


  # Execute WaveQTL and WAveQTL_HMT scripts ---------------------------------

  setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
  #### Run null dataset
  print("Executing null datasets...\n")
  # No HMT
  system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/test_chr17.10161485.cis.geno -p sims/length_",effect_length,"/null_data/WCs.no.QT.txt -u sims/length_",effect_length,"/null_data/use.txt -o sim3_noQT_null -f ",num_bases," -fph 1")
         ,show.output.on.console = F)
  # HMT
  system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/test_chr17.10161485.cis.geno -p sims/length_",effect_length,"/null_data/WCs.no.QT.txt -u sims/length_",effect_length,"/null_data/use.txt -o sim3_noQT_null_HMT -f ",num_bases," -hmt 1")
         ,show.output.on.console = F)
  #### Run alt dataset
  print("Executing alt datasets...\n")
  # No HMT
  system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/test_chr17.10161485.cis.geno -p sims/length_",effect_length,"/alt_data/WCs.no.QT.txt -u sims/length_",effect_length,"/alt_data/use.txt -o sim3_noQT_alt -f ",num_bases," -fph 1")
         ,show.output.on.console = F)
  # HMT
  system(paste0("../../WaveQTL -gmode 1 -g ../../data/dsQTL/test_chr17.10161485.cis.geno -p sims/length_",effect_length,"/alt_data/WCs.no.QT.txt -u sims/length_",effect_length,"/alt_data/use.txt -o sim3_noQT_alt_HMT -f ",num_bases," -hmt 1")
         ,show.output.on.console = F)
  setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")


  # Analysis ----------------------------------------------------------------

  ### Analysis - no HMT
  print("Analysing effect size - no HMT...\n")
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
  print("Analysing effect size - with HMT...\n")
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


# Run ---------------------------------------------------------------------

sim3_50 <- run_sim3(sequencing_sums = seq_sum
                 , num_indivs = 70
                 , num_bases = 1024
                 , effect_length = 50 # must be one of 5, 10, 50
                 , effect_size_data = waveqtl_hmt_geno11
                 , use_qt_data = FALSE
                 , Wmat_1024 = Wmat_1024
                 , W2mat_1024 = W2mat_1024
                 , library.read.depth = library.read.depth
                 , Covariates = Covariates)

sim3_10 <- run_sim3(sequencing_sums = seq_sum
                 , num_indivs = 70
                 , num_bases = 1024
                 , effect_length = 10
                 , effect_size_data = waveqtl_hmt_geno11
                 , use_qt_data = FALSE
                 , Wmat_1024 = Wmat_1024
                 , W2mat_1024 = W2mat_1024
                 , library.read.depth = library.read.depth
                 , Covariates = Covariates)

sim3_5 <- run_sim3(sequencing_sums = seq_sum
                 , num_indivs = 70
                 , num_bases = 1024
                 , effect_length = 5 # must be one of 5, 10, 50
                 , effect_size_data = waveqtl_hmt_geno11
                 , use_qt_data = FALSE
                 , Wmat_1024 = Wmat_1024
                 , W2mat_1024 = W2mat_1024
                 , library.read.depth = library.read.depth
                 , Covariates = Covariates)


sim3_50$null_analysis$p
sim3_50$null_hmt_analysis$p
sim3_50$alt_analysis$p
sim3_50$alt_hmt_analysis$p

sim3_10$null_analysis$p
sim3_10$null_hmt_analysis$p
sim3_10$alt_analysis$p
sim3_10$alt_hmt_analysis$p

sim3_5$null_analysis$p
sim3_5$null_hmt_analysis$p
sim3_5$alt_analysis$p
sim3_5$alt_hmt_analysis$p
