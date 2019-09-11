# 500 sim batches


# Stepped effect size -----------------------------------------------------

number_sims <- 500
od <- 70

set.seed(1)
t_8_s3 <- run_sim4(
  sequencing_sums = effect_size_and_data$seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 8
  , effect_interval = effect_interval_8
  , effect_size_data = effect_size_and_data$effect_size
  , use_qt_data = FALSE
  , over_disp_param = od
  , Wmat_1024 = effect_size_and_data$Wmat_1024
  , W2mat_1024 = effect_size_and_data$W2mat_1024
  , library.read.depth = effect_size_and_data$library.read.depth
  , Covariates = effect_size_and_data$Covariates
  , effect_multiple = 1.5e8
  , trials_multiple = 10
  , number_sims = number_sims
  , verbose = F
  , outputAlias = "sim4"
  , rMarkdownMode = F
)
t_8_s3_diags <- waveqtl_diags(t_8_s3,number_sims)

set.seed(1)
t_16_s3 <- run_sim4(
  sequencing_sums = effect_size_and_data$seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 16
  , effect_interval = effect_interval_16
  , effect_size_data = effect_size_and_data$effect_size
  , use_qt_data = FALSE
  , over_disp_param = od
  , Wmat_1024 = effect_size_and_data$Wmat_1024
  , W2mat_1024 = effect_size_and_data$W2mat_1024
  , library.read.depth = effect_size_and_data$library.read.depth
  , Covariates = effect_size_and_data$Covariates
  , effect_multiple = 8e7
  , trials_multiple = 10
  , number_sims = number_sims
  , verbose = F
  , outputAlias = "sim4"
  , rMarkdownMode = F
)
t_16_s3_diags <- waveqtl_diags(t_16_s3,number_sims)

set.seed(1)
t_32_s3 <- run_sim4(
  sequencing_sums = effect_size_and_data$seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 32
  , effect_interval = effect_interval_32
  , effect_size_data = effect_size_and_data$effect_size
  , use_qt_data = FALSE
  , over_disp_param = od
  , Wmat_1024 = effect_size_and_data$Wmat_1024
  , W2mat_1024 = effect_size_and_data$W2mat_1024
  , library.read.depth = effect_size_and_data$library.read.depth
  , Covariates = effect_size_and_data$Covariates
  , effect_multiple = 4e7
  , trials_multiple = 10
  , number_sims = number_sims
  , verbose = F
  , outputAlias = "sim4"
  , rMarkdownMode = F
)
t_32_s3_diags <- waveqtl_diags(t_32_s3,number_sims)

set.seed(1)
t_64_s3 <- run_sim4(
  sequencing_sums = effect_size_and_data$seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 64
  , effect_interval = effect_interval_64
  , effect_size_data = effect_size_and_data$effect_size
  , use_qt_data = FALSE
  , over_disp_param = od
  , Wmat_1024 = effect_size_and_data$Wmat_1024
  , W2mat_1024 = effect_size_and_data$W2mat_1024
  , library.read.depth = effect_size_and_data$library.read.depth
  , Covariates = effect_size_and_data$Covariates
  , effect_multiple = 2e7
  , trials_multiple = 10
  , number_sims = number_sims
  , verbose = F
  , outputAlias = "sim4"
  , rMarkdownMode = F
)
t_64_s3_diags <- waveqtl_diags(t_64_s3,number_sims)

saveRDS( t_8_s3, "../data/20190908_sim4_l8_s3_stepped_500sims.RDS",compress = T)
saveRDS(t_16_s3,"../data/20190908_sim4_l16_s3_stepped_500sims.RDS",compress = T)
saveRDS(t_32_s3,"../data/20190908_sim4_l32_s3_stepped_500sims.RDS",compress = T)
saveRDS(t_64_s3,"../data/20190908_sim4_l64_s3_stepped_500sims.RDS",compress = T)


# Same effect size --------------------------------------------------------

# Re-run: 8,16,32

efMult <- 1.2e8
num_sums <- 500


set.seed(1)
t_8_s2 <- run_sim4(
  sequencing_sums = effect_size_and_data$seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 8
  , effect_interval = effect_interval_8
  , effect_size_data = effect_size_and_data$effect_size
  , use_qt_data = FALSE
  , over_disp_param = 70
  , Wmat_1024 = effect_size_and_data$Wmat_1024
  , W2mat_1024 = effect_size_and_data$W2mat_1024
  , library.read.depth = effect_size_and_data$library.read.depth
  , Covariates = effect_size_and_data$Covariates
  , effect_multiple = efMult
  , trials_multiple = 10
  , number_sims = num_sums
  , verbose = F
  , outputAlias = "sim3"
  , rMarkdownMode = F
)
t_8_s2_diags <- waveqtl_diags(t_8_s2,num_sums)

set.seed(1)
t_16_s2 <- run_sim4(
  sequencing_sums = effect_size_and_data$seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 16
  , effect_interval = effect_interval_16
  , effect_size_data = effect_size_and_data$effect_size
  , use_qt_data = FALSE
  , over_disp_param = 70
  , Wmat_1024 = effect_size_and_data$Wmat_1024
  , W2mat_1024 = effect_size_and_data$W2mat_1024
  , library.read.depth = effect_size_and_data$library.read.depth
  , Covariates = effect_size_and_data$Covariates
  , effect_multiple = efMult
  , trials_multiple = 10
  , number_sims = num_sums
  , verbose = F
  , outputAlias = "sim3"
  , rMarkdownMode = F
)
t_16_s2_diags <- waveqtl_diags(t_16_s2,num_sums)

set.seed(1)
t_32_s2 <- run_sim4(
  sequencing_sums = effect_size_and_data$seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 32
  , effect_interval = effect_interval_32
  , effect_size_data = effect_size_and_data$effect_size
  , use_qt_data = FALSE
  , over_disp_param = 70
  , Wmat_1024 = effect_size_and_data$Wmat_1024
  , W2mat_1024 = effect_size_and_data$W2mat_1024
  , library.read.depth = effect_size_and_data$library.read.depth
  , Covariates = effect_size_and_data$Covariates
  , effect_multiple = efMult
  , trials_multiple = 10
  , number_sims = num_sums
  , verbose = F
  , outputAlias = "sim3"
  , rMarkdownMode = F
)
t_32_s2_diags <- waveqtl_diags(t_32_s2,num_sums)

set.seed(1)
t_64_s2 <- run_sim4(
  sequencing_sums = effect_size_and_data$seq_sum
  , num_indivs = 70
  , num_bases = 1024
  , effect_length = 64
  , effect_interval = effect_interval_64
  , effect_size_data = effect_size_and_data$effect_size
  , use_qt_data = FALSE
  , over_disp_param = 70
  , Wmat_1024 = effect_size_and_data$Wmat_1024
  , W2mat_1024 = effect_size_and_data$W2mat_1024
  , library.read.depth = effect_size_and_data$library.read.depth
  , Covariates = effect_size_and_data$Covariates
  , effect_multiple = efMult
  , trials_multiple = 10
  , number_sims = num_sums
  , verbose = F
  , outputAlias = "sim3"
  , rMarkdownMode = F
)
t_64_s2_diags <- waveqtl_diags(t_64_s2,num_sums)

saveRDS(t_8_s2,"../data/20190908_sim4_l8_s2_same_500sims.RDS", compress = T)
saveRDS(t_16_s2,"../data/20190908_sim4_l16_s2_same_500sims.RDS", compress = T)
saveRDS(t_32_s2,"../data/20190908_sim4_l32_s2_same_500sims.RDS", compress = T)
saveRDS(t_64_s2,"../data/20190908_sim4_l64_s2_same_500sims.RDS", compress = T)
