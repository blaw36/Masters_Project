### Script runs:
# For a single effect length, tries to find area where the effects are detectable
# and where they fail to become detectable (by HMT)

# Initialisation ----------------------------------------------------------

library(dplyr)
library(ROCR)

# Update this to source sim4_1_4_functions.R, whereever you've placed the script!
source("code/20210106_paper_prep/20211128_chr12_loc171_sims/sim4_1_4_functions.R")

## Reads in phenotype, PCA, library read data, and generates effect size
# through running WaveQTL_HMT on the particular phenotype and genotype data.

effect_size_and_data <- read_in_gen_eff_size_for_paper(geno_select = 19
                                                       ,pheno_data_file = "~/Cpp/20211214_Wqtl_Spartan/data/pheno_data/chr12.pheno.dat.171"
                                                       ,library_read_depth_file = "~/Cpp/20211214_Wqtl_Spartan/data/dsQTL/library.read.depth.dat" # Same as file in WaveQTL repo
                                                       ,covariates_file = "~/Cpp/20211214_Wqtl_Spartan/data/dsQTL/PC4.dat" # Same as file in WaveQTL repo
                                                       ,wmat_matrix_file = "~/Cpp/20211214_Wqtl_Spartan/data/DWT/Wmat_1024" # Same as file in WaveQTL repo
                                                       ,data_path = "~/Cpp/20211214_Wqtl_Spartan/test/dsQTL/output/" # Assumes that both HMT and WQtl outputs in same directory
                                                       ,dataset = "chr.12.171.hmt" # HMT output file prefix
                                                       ,waveqtl_dataset = "chr.12.171.nohmt" # WQtl output file prefix
                                                       )

# Pick effect sizes:
effect_interval_4 <- 519:522
effect_interval_8 <- 519:526
effect_interval_16 <- 450:465
effect_interval_32 <- 437:468
effect_interval_64 <- 417:480

effect_interval_24 <- 442:465
effect_interval_40 <- 433:472
effect_interval_48 <- 500:547
effect_interval_56 <- 417:472


# Parameters --------------------------------------------------------------

## Overdispersion
# Currently at 210 (from thesis). We've also tried 70.
od <- 210

## Effect strength
# We controlled effect strength by multiplying the effect
# size (obtained above) by some number.

effect_multiplier <- 1e8*seq(0.4,1.4,by = 0.1)

## Number of simulations
num_sims <- 200

# Effect length
effect_length <- 8

## Simulation output prefixes and filenames
sim_label <- "paper_p1"
today_dt <- format(Sys.Date(),"%Y%m%d")
savename <- "paper_nogrp_od210"
output_path <- "~/"

# Thesis parameters, for reference ----------------------------------------

## For reference, in our thesis, we tried the following combinations
# of effect lengths and strengths on the WQtl data:

## Sites where we simulated effects
# effect_interval_8 <- 519:526
# effect_interval_16 <- 450:465
# effect_interval_32 <- 437:468
# effect_interval_64 <- 417:480

## Weaker effects as effect length increased
# 8: 1e8*seq(0.4,1.4,by = 0.1)
# 16: 4e7*seq(0.4,1.4,by = 0.1)
# 32: 3e7*seq(0.4,1.4,by = 0.1)
# 64: 1.5e7*seq(0.4,1.4,by = 0.1)

# All overdispersions set to 210.


# Run simulations ---------------------------------------------------------

# Looping the simulations for various effect sizes
for(len in effect_length){

  print(len)
  tmp <- vector("list",length(effect_multiplier))
  n <- 1

  # Looping the simulations for various effect strengths
  for(ef_size in effect_multiplier){
    print(ef_size)
    set.seed(1)
    tmp[[n]] <- run_sim4_v2(
      sequencing_sums = effect_size_and_data$seq_sum
      , num_indivs = 70
      , num_bases = 1024
      , effect_length = len
      , effect_interval = get(paste0("effect_interval_",len))
      , effect_size_data = effect_size_and_data$effect_size
      , no.QT = FALSE
      , over_disp_param = od
      , Wmat_1024 = effect_size_and_data$Wmat_1024
      , W2mat_1024 = effect_size_and_data$W2mat_1024
      , library.read.depth = effect_size_and_data$library.read.depth
      , Covariates = effect_size_and_data$Covariates
      , effect_multiple = ef_size
      , trials_multiple = 10
      , number_sims = num_sims
      , verbose = F
      , outputAlias = sim_label
      , WaveQTL_directory = "~/Cpp/20211214_Wqtl_Spartan/"
    )
    n <- n+1
  }

  assign(x = paste0("l",len,"_res")
         ,value = tmp)
  output_filename <- paste0(today_dt,"_l",len,"_",savename,".RDS")
  saveRDS(tmp, paste0(output_path,output_filename), compress = T)

}
