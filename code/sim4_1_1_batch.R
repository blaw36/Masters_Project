# Running some really large batch jobs

# Computer 1 --------------------------------------------------------------

results_list <- list()

efMult <- 1e8
num_sums <- 50
for(i in 1:20){
  set.seed(i)
  results_list[[i]] <- run_sim4(
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
}

a <- lapply(results_list,waveqtl_diags,num_sums)
unlist(lapply(a,function(x){attr(x[["perf_nohmt"]],"y.values")}))
unlist(lapply(a,function(x){attr(x[["perf_hmt"]],"y.values")}))
a[[1]]$roc_plot
a[[2]]$roc_plot
a[[3]]$roc_plot
a[[4]]$roc_plot
a[[5]]$roc_plot
a[[6]]$roc_plot
a[[7]]$roc_plot
a[[8]]$roc_plot
a[[9]]$roc_plot
a[[10]]$roc_plot
a[[11]]$roc_plot
a[[12]]$roc_plot
a[[13]]$roc_plot
a[[14]]$roc_plot
a[[15]]$roc_plot
a[[16]]$roc_plot
a[[17]]$roc_plot
a[[18]]$roc_plot
a[[19]]$roc_plot
a[[20]]$roc_plot
all_null <- unlist(lapply(results_list,function(x){x[[1]]}))
all_alt <- unlist(lapply(results_list,function(x){x[[2]]}))
all_null_hmt <- unlist(lapply(results_list,function(x){x[[3]]}))
all_alt_hmt <- unlist(lapply(results_list,function(x){x[[4]]}))
all_diag <- waveqtl_diags(sims_list = list(null_waveqtl_lhood = all_null
                                           ,null_waveqtl_hmt_lhood = all_null_hmt
                                           ,alt_waveqtl_lhood = all_alt
                                           ,alt_waveqtl_hmt_lhood = all_alt_hmt)
                          ,num_sims = 1000)
saveRDS(results_list,"data/20190908_sim4_l8_batch.RDS", compress=T)

results_list_64 <- list()
for(i in 1:20){
  set.seed(i)
  results_list_64[[i]] <- run_sim4(
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
    , effect_multiple = 2e7
    , trials_multiple = 10
    , number_sims = num_sums
    , verbose = F
    , outputAlias = "sim3"
    , rMarkdownMode = F
  )
}

a <- lapply(results_list_64,waveqtl_diags,num_sums)
unlist(lapply(a,function(x){attr(x[["perf_nohmt"]],"y.values")}))
unlist(lapply(a,function(x){attr(x[["perf_hmt"]],"y.values")}))
a[[1]]$roc_plot
a[[2]]$roc_plot
a[[3]]$roc_plot
a[[4]]$roc_plot
a[[5]]$roc_plot
a[[6]]$roc_plot
a[[7]]$roc_plot
a[[8]]$roc_plot
a[[9]]$roc_plot
a[[10]]$roc_plot
a[[11]]$roc_plot
a[[12]]$roc_plot
a[[13]]$roc_plot
a[[14]]$roc_plot
a[[15]]$roc_plot
a[[16]]$roc_plot
a[[17]]$roc_plot
a[[18]]$roc_plot
a[[19]]$roc_plot
a[[20]]$roc_plot
all_null <- unlist(lapply(results_list_64,function(x){x[[1]]}))
all_alt <- unlist(lapply(results_list_64,function(x){x[[2]]}))
all_null_hmt <- unlist(lapply(results_list_64,function(x){x[[3]]}))
all_alt_hmt <- unlist(lapply(results_list_64,function(x){x[[4]]}))
all_diag <- waveqtl_diags(sims_list = list(null_waveqtl_lhood = all_null
                                           ,null_waveqtl_hmt_lhood = all_null_hmt
                                           ,alt_waveqtl_lhood = all_alt
                                           ,alt_waveqtl_hmt_lhood = all_alt_hmt)
                          ,num_sims = 1000)
saveRDS(results_list_64,"data/20190908_sim4_l64_batch.RDS", compress=T)

# Ideal cases
# Effect size 1.4

results_list_8_eSize1.4 <- list()

num_sums <- 50
for(i in 1:20){
  set.seed(i)
  results_list_8_eSize1.4[[i]] <- run_sim4(
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
    , effect_size = 1.4
    , trials_multiple = 10
    , number_sims = num_sums
    , verbose = F
    , outputAlias = "sim3"
    , rMarkdownMode = F
  )
}

a <- lapply(results_list_8_eSize1.4,waveqtl_diags,num_sums)
unlist(lapply(a,function(x){attr(x[["perf_nohmt"]],"y.values")}))
unlist(lapply(a,function(x){attr(x[["perf_hmt"]],"y.values")}))
a[[1]]$roc_plot
a[[2]]$roc_plot
a[[3]]$roc_plot
a[[4]]$roc_plot
a[[5]]$roc_plot
a[[6]]$roc_plot
a[[7]]$roc_plot
a[[8]]$roc_plot
a[[9]]$roc_plot
a[[10]]$roc_plot
a[[11]]$roc_plot
a[[12]]$roc_plot
a[[13]]$roc_plot
a[[14]]$roc_plot
a[[15]]$roc_plot
a[[16]]$roc_plot
a[[17]]$roc_plot
a[[18]]$roc_plot
a[[19]]$roc_plot
a[[20]]$roc_plot
all_null <- unlist(lapply(results_list_8_eSize1.4,function(x){x[[1]]}))
all_alt <- unlist(lapply(results_list_8_eSize1.4,function(x){x[[2]]}))
all_null_hmt <- unlist(lapply(results_list_8_eSize1.4,function(x){x[[3]]}))
all_alt_hmt <- unlist(lapply(results_list_8_eSize1.4,function(x){x[[4]]}))
all_diag <- waveqtl_diags(sims_list = list(null_waveqtl_lhood = all_null
                                           ,null_waveqtl_hmt_lhood = all_null_hmt
                                           ,alt_waveqtl_lhood = all_alt
                                           ,alt_waveqtl_hmt_lhood = all_alt_hmt)
                          ,num_sims = 1000)
saveRDS(results_list_8_eSize1.4,"data/20190908_sim4_l8_eSize1-4_batch.RDS", compress=T)

# DIDN'T RUN
# results_list_64_eSize1.4 <- list()
# for(i in 1:20){
#   set.seed(i)
#   results_list_64_eSize1.4[[i]] <- run_sim4(
#     sequencing_sums = effect_size_and_data$seq_sum
#     , num_indivs = 70
#     , num_bases = 1024
#     , effect_length = 64
#     , effect_interval = effect_interval_64
#     , effect_size_data = effect_size_and_data$effect_size
#     , use_qt_data = FALSE
#     , over_disp_param = 70
#     , Wmat_1024 = effect_size_and_data$Wmat_1024
#     , W2mat_1024 = effect_size_and_data$W2mat_1024
#     , library.read.depth = effect_size_and_data$library.read.depth
#     , Covariates = effect_size_and_data$Covariates
#     , effect_size = 1.4
#     , trials_multiple = 10
#     , number_sims = num_sums
#     , verbose = F
#     , outputAlias = "sim3"
#     , rMarkdownMode = F
#   )
# }



# Computer 2 --------------------------------------------------------------

results_list <- list()

num_sums <- 50
for(i in 1:20){
  set.seed(i)
  results_list[[i]] <- run_sim4(
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
    , effect_multiple = 8e7
    , trials_multiple = 10
    , number_sims = number_sims
    , verbose = F
    , outputAlias = "sim4"
    , rMarkdownMode = F
  )
}

a <- lapply(results_list,waveqtl_diags,num_sums)
unlist(lapply(a,function(x){attr(x[["perf_nohmt"]],"y.values")}))
unlist(lapply(a,function(x){attr(x[["perf_hmt"]],"y.values")}))
a[[1]]$roc_plot
a[[2]]$roc_plot
a[[3]]$roc_plot
a[[4]]$roc_plot
a[[5]]$roc_plot
a[[6]]$roc_plot
a[[7]]$roc_plot
a[[8]]$roc_plot
a[[9]]$roc_plot
a[[10]]$roc_plot
a[[11]]$roc_plot
a[[12]]$roc_plot
a[[13]]$roc_plot
a[[14]]$roc_plot
a[[15]]$roc_plot
a[[16]]$roc_plot
a[[17]]$roc_plot
a[[18]]$roc_plot
a[[19]]$roc_plot
a[[20]]$roc_plot
all_null <- unlist(lapply(results_list,function(x){x[[1]]}))
all_alt <- unlist(lapply(results_list,function(x){x[[2]]}))
all_null_hmt <- unlist(lapply(results_list,function(x){x[[3]]}))
all_alt_hmt <- unlist(lapply(results_list,function(x){x[[4]]}))
all_diag <- waveqtl_diags(sims_list = list(null_waveqtl_lhood = all_null
                                           ,null_waveqtl_hmt_lhood = all_null_hmt
                                           ,alt_waveqtl_lhood = all_alt
                                           ,alt_waveqtl_hmt_lhood = all_alt_hmt)
                          ,num_sims = 1000)
saveRDS(results_list,"data/20190908_sim4_l16_batch.RDS", compress=T)

results_list_32 <- list()

num_sums <- 50
for(i in 1:20){
  set.seed(i)
  results_list_32[[i]] <- run_sim4(
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
    , effect_multiple = 4e7
    , trials_multiple = 10
    , number_sims = number_sims
    , verbose = F
    , outputAlias = "sim4"
    , rMarkdownMode = F
  )
}
a <- lapply(results_list_32,waveqtl_diags,num_sums)
unlist(lapply(a,function(x){attr(x[["perf_nohmt"]],"y.values")}))
unlist(lapply(a,function(x){attr(x[["perf_hmt"]],"y.values")}))
a[[1]]$roc_plot
a[[2]]$roc_plot
a[[3]]$roc_plot
a[[4]]$roc_plot
a[[5]]$roc_plot
a[[6]]$roc_plot
a[[7]]$roc_plot
a[[8]]$roc_plot
a[[9]]$roc_plot
a[[10]]$roc_plot
a[[11]]$roc_plot
a[[12]]$roc_plot
a[[13]]$roc_plot
a[[14]]$roc_plot
a[[15]]$roc_plot
a[[16]]$roc_plot
a[[17]]$roc_plot
a[[18]]$roc_plot
a[[19]]$roc_plot
a[[20]]$roc_plot
all_null <- unlist(lapply(results_list_32,function(x){x[[1]]}))
all_alt <- unlist(lapply(results_list_32,function(x){x[[2]]}))
all_null_hmt <- unlist(lapply(results_list_32,function(x){x[[3]]}))
all_alt_hmt <- unlist(lapply(results_list_32,function(x){x[[4]]}))
all_diag <- waveqtl_diags(sims_list = list(null_waveqtl_lhood = all_null
                                           ,null_waveqtl_hmt_lhood = all_null_hmt
                                           ,alt_waveqtl_lhood = all_alt
                                           ,alt_waveqtl_hmt_lhood = all_alt_hmt)
                          ,num_sims = 1000)
saveRDS(results_list_32,"data/20190908_sim4_l32_batch.RDS", compress=T)

# Ideal cases
# Effect size 1.4
results_list_16_eSize1.4 <- list()

num_sums <- 50
for(i in 1:20){
  set.seed(i)
  results_list_16_eSize1.4[[i]] <- run_sim4(
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
    , effect_size = 1.4
    , trials_multiple = 10
    , number_sims = number_sims
    , verbose = F
    , outputAlias = "sim4"
    , rMarkdownMode = F
  )
}
a <- lapply(results_list_16_eSize1.4,waveqtl_diags,num_sums)
unlist(lapply(a,function(x){attr(x[["perf_nohmt"]],"y.values")}))
unlist(lapply(a,function(x){attr(x[["perf_hmt"]],"y.values")}))
a[[1]]$roc_plot
a[[2]]$roc_plot
a[[3]]$roc_plot
a[[4]]$roc_plot
a[[5]]$roc_plot
a[[6]]$roc_plot
a[[7]]$roc_plot
a[[8]]$roc_plot
a[[9]]$roc_plot
a[[10]]$roc_plot
a[[11]]$roc_plot
a[[12]]$roc_plot
a[[13]]$roc_plot
a[[14]]$roc_plot
a[[15]]$roc_plot
a[[16]]$roc_plot
a[[17]]$roc_plot
a[[18]]$roc_plot
a[[19]]$roc_plot
a[[20]]$roc_plot
all_null <- unlist(lapply(results_list_16_eSize1.4,function(x){x[[1]]}))
all_alt <- unlist(lapply(results_list_16_eSize1.4,function(x){x[[2]]}))
all_null_hmt <- unlist(lapply(results_list_16_eSize1.4,function(x){x[[3]]}))
all_alt_hmt <- unlist(lapply(results_list_16_eSize1.4,function(x){x[[4]]}))
all_diag <- waveqtl_diags(sims_list = list(null_waveqtl_lhood = all_null
                                           ,null_waveqtl_hmt_lhood = all_null_hmt
                                           ,alt_waveqtl_lhood = all_alt
                                           ,alt_waveqtl_hmt_lhood = all_alt_hmt)
                          ,num_sims = 1000)
saveRDS(results_list_16_eSize1.4,"data/20190908_sim4_l16_eSize1-4__batch.RDS", compress=T)

### NOT RUN
# results_list_32_eSize1.4 <- list()
# num_sums <- 50
# for(i in 1:20){
#   set.seed(i)
#   results_list_32_eSize1.4[[i]] <- run_sim4(
#     sequencing_sums = effect_size_and_data$seq_sum
#     , num_indivs = 70
#     , num_bases = 1024
#     , effect_length = 32
#     , effect_interval = effect_interval_32
#     , effect_size_data = effect_size_and_data$effect_size
#     , use_qt_data = FALSE
#     , over_disp_param = 70
#     , Wmat_1024 = effect_size_and_data$Wmat_1024
#     , W2mat_1024 = effect_size_and_data$W2mat_1024
#     , library.read.depth = effect_size_and_data$library.read.depth
#     , Covariates = effect_size_and_data$Covariates
#     , effect_size = 1.4
#     , trials_multiple = 10
#     , number_sims = number_sims
#     , verbose = F
#     , outputAlias = "sim4"
#     , rMarkdownMode = F
#   )
# }
