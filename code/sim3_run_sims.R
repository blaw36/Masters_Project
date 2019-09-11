# ## Experiments on sim3.
# rm(list = ls());gc();cat("\014");
# source("code/sim3_wflow.R")
# int_table <- summarise_effect_intervals(waveqtl_hmt_geno11$col_posi)
#
# # Length 50 ---------------------------------------------------------------
#
# set.seed(10)
# effect_length = 50
# effect_interval_50 <- effect_length_picker(int_table,effect_length)
#
# # Actual effect
#
# effect_size_constant <- 10^seq(0,5)
# overdisp <- 70^seq(-5,3)
#
# sim3_50_actual_eff_res <- list()
# n <- 1
#
# for(i in effect_size_constant){
#   for(j in overdisp){
#     set.seed(1)
#     cat(paste("effect_size_const: ",i,"...\n"))
#     cat(paste("over-dispersion: ",j,"...\n"))
#     mod_waveqtl_hmt_geno11 <- waveqtl_hmt_geno11
#     mod_waveqtl_hmt_geno11$beta_dataS <- mod_waveqtl_hmt_geno11$beta_dataS*i
#     sim3_50_actual_eff_res[[n]] <- run_sim3(sequencing_sums = seq_sum
#                                             , num_indivs = 70
#                                             , num_bases = 1024
#                                             , effect_length = effect_length # must be one of 8,16,32,64
#                                             , effect_interval = effect_interval_50
#                                             , effect_size_data = mod_waveqtl_hmt_geno11
#                                             , use_qt_data = FALSE
#                                             , over_disp_param = j
#                                             , Wmat_1024 = Wmat_1024
#                                             , W2mat_1024 = W2mat_1024
#                                             , library.read.depth = library.read.depth
#                                             , Covariates = Covariates
#                                             , group_data = group_data)
#     sim3_50_actual_eff_res[[n]]$effect_size_constant <- i
#     sim3_50_actual_eff_res[[n]]$over_disp <- j
#     n <- n+1
#   }
# }
#
# saveRDS(sim3_50_actual_eff_res,"~/Cpp/WaveQTL_HMT/sim3_results/sim3_50_actual_eff_res.RDS", compress = T)
#
# sim3_50_actual_eff_res[[1]]$alt_analysis$p
# sim3_50_actual_eff_res[[1]]$alt_hmt_analysis$p
# sim3_50_actual_eff_res[[2]]$alt_analysis$p
# sim3_50_actual_eff_res[[2]]$alt_hmt_analysis$p
# sim3_50_actual_eff_res[[3]]$alt_analysis$p
# sim3_50_actual_eff_res[[3]]$alt_hmt_analysis$p
# sim3_50_actual_eff_res[[4]]$alt_analysis$p
# sim3_50_actual_eff_res[[4]]$alt_hmt_analysis$p
# sim3_50_actual_eff_res[[5]]$alt_analysis$p
# sim3_50_actual_eff_res[[5]]$alt_hmt_analysis$p
# sim3_50_actual_eff_res[[6]]$alt_analysis$p
# sim3_50_actual_eff_res[[6]]$alt_hmt_analysis$p
# sim3_50_actual_eff_res[[7]]$alt_analysis$p
# sim3_50_actual_eff_res[[7]]$alt_hmt_analysis$p
# sim3_50_actual_eff_res[[8]]$alt_analysis$p
# sim3_50_actual_eff_res[[8]]$alt_hmt_analysis$p
# sim3_50_actual_eff_res[[9]]$alt_analysis$p
# sim3_50_actual_eff_res[[9]]$alt_hmt_analysis$p
#
# # Plot all alts
# # for(i in 1:9){ # for a scale, by overdisp
# for(i in seq(4,54,9)){ # for overdisp, by scale
#   both_plot_data <- sim3_50_actual_eff_res[[i]]
#   y_min <- min(min(both_plot_data$alt_analysis$beta_dataS - 3*both_plot_data$alt_analysis$beta_dataS)
#                ,min(both_plot_data$alt_hmt_analysis$beta_dataS - 3*both_plot_data$alt_hmt_analysis$beta_dataS))
#   y_max <- max(max(both_plot_data$alt_analysis$beta_dataS + 3*both_plot_data$alt_analysis$beta_dataS)
#                ,max(both_plot_data$alt_hmt_analysis$beta_dataS + 3*both_plot_data$alt_hmt_analysis$beta_dataS))
#   effect_size_plot(y_min = y_min
#                    ,y_max = y_max
#                    ,beta_mean = both_plot_data$alt_analysis$beta_dataS
#                    ,beta_sd = both_plot_data$alt_analysis$beta_sd_dataS
#                    ,plot_title = paste0("overdisp: ", both_plot_data$over_disp,", effect-scale: ",both_plot_data$effect_size_constant," - no hmt"))
#   effect_size_plot(y_min = y_min
#                    ,y_max = y_max
#                    ,beta_mean = both_plot_data$alt_hmt_analysis$beta_dataS
#                    ,beta_sd = both_plot_data$alt_hmt_analysis$beta_sd_dataS
#                    ,plot_title = paste0("overdisp: ", both_plot_data$over_disp,", effect-scale: ",both_plot_data$effect_size_constant," - hmt"))
#
# }
#
#
# # 1/0 effect (either in region, or not)
#
# effect_size_constant <- 10^seq(0,5)
# overdisp <- 70^seq(-5,3)
#
# sim3_50_bin_eff_res <- list()
# n <- 1
#
# for(i in effect_size_constant){
#   for(j in overdisp){
#     set.seed(1)
#     cat(paste("effect_size_const: ",i,"...\n"))
#     cat(paste("over-dispersion: ",j,"...\n"))
#     mod_waveqtl_hmt_geno11 <- waveqtl_hmt_geno11
#     mod_waveqtl_hmt_geno11$beta_dataS[mod_waveqtl_hmt_geno11$col_posi] <- 1*i
#     mod_waveqtl_hmt_geno11$beta_dataS[!mod_waveqtl_hmt_geno11$col_posi] <- 0*i
#     sim3_50_bin_eff_res[[n]] <- run_sim3(sequencing_sums = seq_sum
#                                             , num_indivs = 70
#                                             , num_bases = 1024
#                                             , effect_length = effect_length # must be one of 8,16,32,64
#                                             , effect_interval = effect_interval_50
#                                             , effect_size_data = mod_waveqtl_hmt_geno11
#                                             , use_qt_data = FALSE
#                                             , over_disp_param = j
#                                             , Wmat_1024 = Wmat_1024
#                                             , W2mat_1024 = W2mat_1024
#                                             , library.read.depth = library.read.depth
#                                             , Covariates = Covariates
#                                             , group_data = group_data)
#     sim3_50_bin_eff_res[[n]]$effect_size_constant <- i
#     sim3_50_bin_eff_res[[n]]$over_disp <- j
#     n <- n+1
#   }
# }
# saveRDS(sim3_50_bin_eff_res,"~/Cpp/WaveQTL_HMT/sim3_results/sim3_50_bin_eff_res.RDS", compress = T)
#
# # Plot all alts
# # for(i in 1:9){ # for a scale, by overdisp
# for(i in seq(4,54,9)){ # for overdisp, by scale
#   both_plot_data <- sim3_50_bin_eff_res[[i]]
#   y_min <- min(min(both_plot_data$alt_analysis$beta_dataS - 3*both_plot_data$alt_analysis$beta_dataS)
#                ,min(both_plot_data$alt_hmt_analysis$beta_dataS - 3*both_plot_data$alt_hmt_analysis$beta_dataS))
#   y_max <- max(max(both_plot_data$alt_analysis$beta_dataS + 3*both_plot_data$alt_analysis$beta_dataS)
#                ,max(both_plot_data$alt_hmt_analysis$beta_dataS + 3*both_plot_data$alt_hmt_analysis$beta_dataS))
#   effect_size_plot(y_min = y_min
#                    ,y_max = y_max
#                    ,beta_mean = both_plot_data$alt_analysis$beta_dataS
#                    ,beta_sd = both_plot_data$alt_analysis$beta_sd_dataS
#                    ,plot_title = paste0("overdisp: ", both_plot_data$over_disp,", effect-scale: ",both_plot_data$effect_size_constant," - no hmt"))
#   effect_size_plot(y_min = y_min
#                    ,y_max = y_max
#                    ,beta_mean = both_plot_data$alt_hmt_analysis$beta_dataS
#                    ,beta_sd = both_plot_data$alt_hmt_analysis$beta_sd_dataS
#                    ,plot_title = paste0("overdisp: ", both_plot_data$over_disp,", effect-scale: ",both_plot_data$effect_size_constant," - hmt"))
#
# }
#
# # scale up the sum_seqs, with 1/0 effect
#
# effect_size_constant <- 10^seq(0,5)
# overdisp <- 70^seq(-5,3)
#
# sim3_50_sum_seq_scale <- list()
# n <- 1
#
# for(i in effect_size_constant){
#   for(j in overdisp){
#     set.seed(1)
#     cat(paste("effect_size_const: ",i,"...\n"))
#     cat(paste("over-dispersion: ",j,"...\n"))
#     mod_waveqtl_hmt_geno11 <- waveqtl_hmt_geno11
#     mod_waveqtl_hmt_geno11$beta_dataS[mod_waveqtl_hmt_geno11$col_posi] <- 1*i
#     mod_waveqtl_hmt_geno11$beta_dataS[!mod_waveqtl_hmt_geno11$col_posi] <- 0*i
#     sim3_50_sum_seq_scale[[n]] <- run_sim3(sequencing_sums = 10*seq_sum
#                                          , num_indivs = 70
#                                          , num_bases = 1024
#                                          , effect_length = effect_length # must be one of 8,16,32,64
#                                          , effect_interval = effect_interval_50
#                                          , effect_size_data = mod_waveqtl_hmt_geno11
#                                          , use_qt_data = FALSE
#                                          , over_disp_param = j
#                                          , Wmat_1024 = Wmat_1024
#                                          , W2mat_1024 = W2mat_1024
#                                          , library.read.depth = library.read.depth
#                                          , Covariates = Covariates
#                                          , group_data = group_data)
#     sim3_50_sum_seq_scale[[n]]$effect_size_constant <- i
#     sim3_50_sum_seq_scale[[n]]$over_disp <- j
#     n <- n+1
#   }
# }
# saveRDS(sim3_50_sum_seq_scale,"~/Cpp/WaveQTL_HMT/sim3_results/sim3_50_sum_seq_scale.RDS", compress = T)
#
# # Plot all alts
# # for(i in 1:9){ # for a scale, by overdisp
# for(i in seq(4,54,9)){ # for overdisp, by scale
#   both_plot_data <- sim3_50_sum_seq_scale[[i]]
#   y_min <- min(min(both_plot_data$alt_analysis$beta_dataS - 3*both_plot_data$alt_analysis$beta_dataS)
#                ,min(both_plot_data$alt_hmt_analysis$beta_dataS - 3*both_plot_data$alt_hmt_analysis$beta_dataS))
#   y_max <- max(max(both_plot_data$alt_analysis$beta_dataS + 3*both_plot_data$alt_analysis$beta_dataS)
#                ,max(both_plot_data$alt_hmt_analysis$beta_dataS + 3*both_plot_data$alt_hmt_analysis$beta_dataS))
#   effect_size_plot(y_min = y_min
#                    ,y_max = y_max
#                    ,beta_mean = both_plot_data$alt_analysis$beta_dataS
#                    ,beta_sd = both_plot_data$alt_analysis$beta_sd_dataS
#                    ,plot_title = paste0("overdisp: ", both_plot_data$over_disp,", effect-scale: ",both_plot_data$effect_size_constant," - no hmt"))
#   effect_size_plot(y_min = y_min
#                    ,y_max = y_max
#                    ,beta_mean = both_plot_data$alt_hmt_analysis$beta_dataS
#                    ,beta_sd = both_plot_data$alt_hmt_analysis$beta_sd_dataS
#                    ,plot_title = paste0("overdisp: ", both_plot_data$over_disp,", effect-scale: ",both_plot_data$effect_size_constant," - hmt"))
#
# }
#
#
# # Length 32 ---------------------------------------------------------------
#
#
#
# # Length 16 ---------------------------------------------------------------
#
#
# set.seed(10)
# effect_length = 16
# effect_interval_16 <- effect_length_picker(int_table,effect_length)
#
# # Actual effect
#
# effect_size_constant <- 10^seq(0,5)
# overdisp <- 70^seq(-5,3)
#
# sim3_16_actual_eff_res <- list()
# n <- 1
#
# for(i in effect_size_constant){
#   for(j in overdisp){
#     set.seed(1)
#     cat(paste("effect_size_const: ",i,"...\n"))
#     cat(paste("over-dispersion: ",j,"...\n"))
#     mod_waveqtl_hmt_geno11 <- waveqtl_hmt_geno11
#     mod_waveqtl_hmt_geno11$beta_dataS <- mod_waveqtl_hmt_geno11$beta_dataS*i
#     sim3_16_actual_eff_res[[n]] <- run_sim3(sequencing_sums = seq_sum
#                                             , num_indivs = 70
#                                             , num_bases = 1024
#                                             , effect_length = effect_length # must be one of 8,16,32,64
#                                             , effect_interval = effect_interval_16
#                                             , effect_size_data = mod_waveqtl_hmt_geno11
#                                             , use_qt_data = FALSE
#                                             , over_disp_param = j
#                                             , Wmat_1024 = Wmat_1024
#                                             , W2mat_1024 = W2mat_1024
#                                             , library.read.depth = library.read.depth
#                                             , Covariates = Covariates
#                                             , group_data = group_data)
#     sim3_16_actual_eff_res[[n]]$effect_size_constant <- i
#     sim3_16_actual_eff_res[[n]]$over_disp <- j
#     n <- n+1
#   }
# }
#
# saveRDS(sim3_16_actual_eff_res,"~/Cpp/WaveQTL_HMT/sim3_results/sim3_16_actual_eff_res.RDS", compress = T)
#
# sim3_16_actual_eff_res[[1]]$alt_analysis$p
# sim3_16_actual_eff_res[[1]]$alt_hmt_analysis$p
# sim3_16_actual_eff_res[[2]]$alt_analysis$p
# sim3_16_actual_eff_res[[2]]$alt_hmt_analysis$p
# sim3_16_actual_eff_res[[3]]$alt_analysis$p
# sim3_16_actual_eff_res[[3]]$alt_hmt_analysis$p
# sim3_16_actual_eff_res[[4]]$alt_analysis$p
# sim3_16_actual_eff_res[[4]]$alt_hmt_analysis$p
# sim3_16_actual_eff_res[[5]]$alt_analysis$p
# sim3_16_actual_eff_res[[5]]$alt_hmt_analysis$p
# sim3_16_actual_eff_res[[6]]$alt_analysis$p
# sim3_16_actual_eff_res[[6]]$alt_hmt_analysis$p
# sim3_16_actual_eff_res[[7]]$alt_analysis$p
# sim3_16_actual_eff_res[[7]]$alt_hmt_analysis$p
# sim3_16_actual_eff_res[[8]]$alt_analysis$p
# sim3_16_actual_eff_res[[8]]$alt_hmt_analysis$p
# sim3_16_actual_eff_res[[9]]$alt_analysis$p
# sim3_16_actual_eff_res[[9]]$alt_hmt_analysis$p
#
# # Plot all alts
# # for(i in 1:9){ # for a scale, by overdisp
# for(i in seq(4,54,9)){ # for overdisp, by scale
#   both_plot_data <- sim3_16_actual_eff_res[[i]]
#   y_min <- min(min(both_plot_data$alt_analysis$beta_dataS - 3*both_plot_data$alt_analysis$beta_dataS)
#                ,min(both_plot_data$alt_hmt_analysis$beta_dataS - 3*both_plot_data$alt_hmt_analysis$beta_dataS))
#   y_max <- max(max(both_plot_data$alt_analysis$beta_dataS + 3*both_plot_data$alt_analysis$beta_dataS)
#                ,max(both_plot_data$alt_hmt_analysis$beta_dataS + 3*both_plot_data$alt_hmt_analysis$beta_dataS))
#   effect_size_plot(y_min = y_min
#                    ,y_max = y_max
#                    ,beta_mean = both_plot_data$alt_analysis$beta_dataS
#                    ,beta_sd = both_plot_data$alt_analysis$beta_sd_dataS
#                    ,plot_title = paste0("overdisp: ", both_plot_data$over_disp,", effect-scale: ",both_plot_data$effect_size_constant," - no hmt"))
#   effect_size_plot(y_min = y_min
#                    ,y_max = y_max
#                    ,beta_mean = both_plot_data$alt_hmt_analysis$beta_dataS
#                    ,beta_sd = both_plot_data$alt_hmt_analysis$beta_sd_dataS
#                    ,plot_title = paste0("overdisp: ", both_plot_data$over_disp,", effect-scale: ",both_plot_data$effect_size_constant," - hmt"))
#
# }
#
# # Effect scale not much effect. Library read depth governs this a lot more.
# # Overdispersion, or effect shape will probably do more
# # Effect length also.
#
# # Length 8 ----------------------------------------------------------------
#
#
#
