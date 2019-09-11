number_sims <- 200
od <- 70
# mult <- 0.9

t_8_500_list <- list()
t_16_500_list <- list()
t_32_500_list <- list()
t_64_500_list <- list()

n <- 1
for(mult in seq(0.6,1.4,0.2)){

  set.seed(1)
  t_8_500_list[[n]] <- run_sim4(
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
    , effect_multiple = 1.5e8*mult
    , trials_multiple = 10
    , number_sims = number_sims
    , verbose = F
    , outputAlias = "sim4"
    , rMarkdownMode = F
  )

  set.seed(1)
  t_16_500_list[[n]] <- run_sim4(
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
    , effect_multiple = 8e7*mult
    , trials_multiple = 10
    , number_sims = number_sims
    , verbose = F
    , outputAlias = "sim4"
    , rMarkdownMode = F
  )

  set.seed(1)
  t_32_500_list[[n]] <- run_sim4(
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
    , effect_multiple = 4e7*mult
    , trials_multiple = 10
    , number_sims = number_sims
    , verbose = F
    , outputAlias = "sim3"
    , rMarkdownMode = F
  )

  set.seed(1)
  t_64_500_list[[n]] <- run_sim4(
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
    , effect_multiple = 2e7*mult
    , trials_multiple = 10
    , number_sims = number_sims
    , verbose = F
    , outputAlias = "sim3"
    , rMarkdownMode = F
  )

  n <- n+1
}

saveRDS(t_8_500_list,"data/20190911_sim4_l8_manyEfs_200sims.RDS",compress = T)
saveRDS(t_16_500_list,"data/20190911_sim4_l16_manyEfs_200sims.RDS",compress = T)
saveRDS(t_32_500_list,"data/20190911_sim4_l32_manyEfs_200sims.RDS",compress = T)
saveRDS(t_64_500_list,"data/20190911_sim4_l64_manyEfs_200sims.RDS",compress = T)
