## Experiments on sim3.
# These are using the 'test' function to generate ridiculous effect sizes under controlled conditions

rm(list = ls());gc();cat("\014");
source("code/sim3_wflow.R")
int_table <- summarise_effect_intervals(waveqtl_hmt_geno11$col_posi)

# Generate effect sizes
set.seed(10)
effect_interval_8 <- effect_length_picker(int_table,8)
effect_interval_16 <- effect_length_picker(int_table,16)
effect_interval_32 <- effect_length_picker(int_table,32)
effect_interval_64 <- effect_length_picker(int_table,64)

# Experiments -------------------------------------------------------------

# Length 8 ----------------------------------------------------------------
effect_size_vect <- seq(1.6,2,0.2)
l8_res_list <- list()
n <- 1
for(i in 1:length(effect_size_vect)){
# for(i in 1:1){#length(effect_size_vect)){
    set.seed(10)
    l8_res_list[[n]] <- run_sim3_test(
      sequencing_sums = seq_sum
      , num_indivs = 70
      , num_bases = 1024
      , effect_length = 8
      , effect_interval = effect_interval_8
      , effect_size_data = waveqtl_hmt_geno11
      , use_qt_data = FALSE
      , over_disp_param = 700
      , Wmat_1024 = Wmat_1024
      , W2mat_1024 = W2mat_1024
      , library.read.depth = library.read.depth
      , Covariates = Covariates
      , group_data = group_data
      , effect_size = effect_size_vect[i]
      , num_trials = 100
      , trials_multiple = 1
    )
    n <- n+1
}

# # Good cases:
# # Effect between 1.6, 2
# # Overdispersion: 700
# # Trials: 100
# saveRDS(l8_res_list,"data/sim3_length8.RDS",compress = T)

# Now, over-dispersion between 1, 700
over_disp <- c(1,100,300,500,700)
l8_res_list <- list()
n <- 1
for(i in 1:length(over_disp)){
  # for(i in 1:1){#length(effect_size_vect)){
  set.seed(10)
  l8_res_list[[n]] <- run_sim3_test(
    sequencing_sums = seq_sum
    , num_indivs = 70
    , num_bases = 1024
    , effect_length = 8
    , effect_interval = effect_interval_8
    , effect_size_data = waveqtl_hmt_geno11
    , use_qt_data = FALSE
    , over_disp_param = over_disp[i]
    , Wmat_1024 = Wmat_1024
    , W2mat_1024 = W2mat_1024
    , library.read.depth = library.read.depth
    , Covariates = Covariates
    , group_data = group_data
    , effect_size = 1.8
    , num_trials = 100
    , trials_multiple = 1
  )
  n <- n+1
}
# # Another good case, (1.8, 300, 100)
# saveRDS(l8_res_list[[3]],"data/sim3_length8_good2.RDS",compress = T)

## Would like to show, around any inflection point,
# Over, say, 20 runs, how reproducibly accurate is the HMT result over the non-HMT?
# Involves finding sweet spots at all 3/4 signal lengths

l8_res_list <- list()
n <- 1
set.seed(10)
for(i in 1:20){
  l8_res_list[[n]] <- run_sim3_test(
    sequencing_sums = seq_sum
    , num_indivs = 70
    , num_bases = 1024
    , effect_length = 8
    , effect_interval = effect_interval_8
    , effect_size_data = waveqtl_hmt_geno11
    , use_qt_data = FALSE
    , over_disp_param = 300
    , Wmat_1024 = Wmat_1024
    , W2mat_1024 = W2mat_1024
    , library.read.depth = library.read.depth
    , Covariates = Covariates
    , group_data = group_data
    , effect_size = 1.8
    , num_trials = 100
    , trials_multiple = 1
  )
  n <- n+1
}

# saveRDS(l8_res_list,"data/sim3_length8_od300_eff1.8_tri100.RDS",compress = T)
alt_analysis_findings <- lapply(l8_res_list,function(x){x[["alt_analysis"]][["col_posi"]]})
alt_analysis_hmt_findings <- lapply(l8_res_list,function(x){x[["alt_hmt_analysis"]][["col_posi"]]})

# How many in the interval?
alt_analysis_findings_tp <- unlist(lapply(lapply(alt_analysis_findings,function(x){intersect(x,effect_interval_8)}), length))
alt_analysis_hmt_findings_tp <- unlist(lapply(lapply(alt_analysis_hmt_findings,function(x){intersect(x,effect_interval_8)}), length))
# How many not in the interval?
alt_analysis_findings_fp <- unlist(lapply(lapply(alt_analysis_findings,function(x){setdiff(x,effect_interval_8)}), length))
alt_analysis_hmt_findings_fp <- unlist(lapply(lapply(alt_analysis_hmt_findings,function(x){setdiff(x,effect_interval_8)}), length))

# Plot?
plot(0,0,xlim = c(1,20), ylim = c(min(min(alt_analysis_findings_tp),min(alt_analysis_hmt_findings_tp))
                                  ,max(max(alt_analysis_findings_tp),max(alt_analysis_hmt_findings_tp))), main = "True positives")
lines(x = seq(1,20,1), alt_analysis_findings_tp,col = "red")
lines(x = seq(1,20,1), alt_analysis_hmt_findings_tp,col = "blue")
legend("topright",legend = c("noHMT","HMT"),col = c("red","blue"),lty = rep(1,2))

plot(0,0,xlim = c(1,20), ylim = c(min(min(alt_analysis_findings_fp),min(alt_analysis_hmt_findings_fp))
                                  ,max(max(alt_analysis_findings_fp),max(alt_analysis_hmt_findings_fp))), main = "False positives")
lines(x = seq(1,20,1), alt_analysis_findings_fp,col = "red")
lines(x = seq(1,20,1), alt_analysis_hmt_findings_fp,col = "blue")
legend("topright",legend = c("noHMT","HMT"),col = c("red","blue"),lty = rep(1,2))


# Length 16 ---------------------------------------------------------------
##### Try something similar for length 16?

l16_res_list <- list()
n <- 1
set.seed(10)
for(i in 1:20){
  l16_res_list[[n]] <- run_sim3_test(
    sequencing_sums = seq_sum
    , num_indivs = 70
    , num_bases = 1024
    , effect_length = 16
    , effect_interval = effect_interval_16
    , effect_size_data = waveqtl_hmt_geno11
    , use_qt_data = FALSE
    , over_disp_param = 300
    , Wmat_1024 = Wmat_1024
    , W2mat_1024 = W2mat_1024
    , library.read.depth = library.read.depth
    , Covariates = Covariates
    , group_data = group_data
    , effect_size = 1.4
    , num_trials = 100
    , trials_multiple = 1
  )
  n <- n+1
}

saveRDS(l16_res_list,"data/sim3_length16_od300_eff1.8_tri100.RDS",compress = T)
alt_analysis_findings <- lapply(l16_res_list,function(x){x[["alt_analysis"]][["col_posi"]]})
alt_analysis_hmt_findings <- lapply(l16_res_list,function(x){x[["alt_hmt_analysis"]][["col_posi"]]})

# How many in the interval?
alt_analysis_findings_tp <- unlist(lapply(lapply(alt_analysis_findings,function(x){intersect(x,effect_interval_16)}), length))
alt_analysis_hmt_findings_tp <- unlist(lapply(lapply(alt_analysis_hmt_findings,function(x){intersect(x,effect_interval_16)}), length))
# How many not in the interval?
alt_analysis_findings_fp <- unlist(lapply(lapply(alt_analysis_findings,function(x){setdiff(x,effect_interval_16)}), length))
alt_analysis_hmt_findings_fp <- unlist(lapply(lapply(alt_analysis_hmt_findings,function(x){setdiff(x,effect_interval_16)}), length))

# Plot?
plot(0,0,xlim = c(1,20), ylim = c(min(min(alt_analysis_findings_tp),min(alt_analysis_hmt_findings_tp))
                                  ,max(max(alt_analysis_findings_tp),max(alt_analysis_hmt_findings_tp))), main = "True positives")
lines(x = seq(1,20,1), alt_analysis_findings_tp,col = "red")
lines(x = seq(1,20,1), alt_analysis_hmt_findings_tp,col = "blue")
legend("topright",legend = c("noHMT","HMT"),col = c("red","blue"),lty = rep(1,2))

plot(0,0,xlim = c(1,20), ylim = c(min(min(alt_analysis_findings_fp),min(alt_analysis_hmt_findings_fp))
                                  ,max(max(alt_analysis_findings_fp),max(alt_analysis_hmt_findings_fp))), main = "False positives")
lines(x = seq(1,20,1), alt_analysis_findings_fp,col = "red")
lines(x = seq(1,20,1), alt_analysis_hmt_findings_fp,col = "blue")
legend("topright",legend = c("noHMT","HMT"),col = c("red","blue"),lty = rep(1,2))


# Higher over-disp --------------------------------------------------------

### Try it again on a more ideal overdispersion of 1000?

# Dropping assumptions ----------------------------------------------------

# Then, drop some of the assumptions:
  # The main one is seq_sum. Scale it by a constant.

l8_res_list <- list()
n <- 1
set.seed(10)
for(i in 1:20){
  l8_res_list[[n]] <- run_sim3_test(
    sequencing_sums = seq_sum
    , num_indivs = 70
    , num_bases = 1024
    , effect_length = 8
    , effect_interval = effect_interval_8
    , effect_size_data = waveqtl_hmt_geno11
    , use_qt_data = FALSE
    , over_disp_param = 1000
    , Wmat_1024 = Wmat_1024
    , W2mat_1024 = W2mat_1024
    , library.read.depth = library.read.depth
    , Covariates = Covariates
    , group_data = group_data
    , effect_size = 2
    , trials_multiple = 5
  )
  n <- n+1
}

l8_res_list_100 <- list()
n <- 1
set.seed(10)
for(i in 1:20){
  l8_res_list_100[[n]] <- run_sim3_test(
    sequencing_sums = seq_sum
    , num_indivs = 70
    , num_bases = 1024
    , effect_length = 8
    , effect_interval = effect_interval_8
    , effect_size_data = waveqtl_hmt_geno11
    , use_qt_data = FALSE
    , over_disp_param = 1000
    , Wmat_1024 = Wmat_1024
    , W2mat_1024 = W2mat_1024
    , library.read.depth = library.read.depth
    , Covariates = Covariates
    , group_data = group_data
    , effect_size = 2
    , trials_multiple = 10
  )
  n <- n+1
}

l16_res_list_10 <- list()
n <- 1
set.seed(10)
for(i in 1:20){
  l16_res_list_10[[n]] <- run_sim3_test(
    sequencing_sums = seq_sum
    , num_indivs = 70
    , num_bases = 1024
    , effect_length = 16
    , effect_interval = effect_interval_16
    , effect_size_data = waveqtl_hmt_geno11
    , use_qt_data = FALSE
    , over_disp_param = 1000
    , Wmat_1024 = Wmat_1024
    , W2mat_1024 = W2mat_1024
    , library.read.depth = library.read.depth
    , Covariates = Covariates
    , group_data = group_data
    , effect_size = 1.4
    , trials_multiple = 10
  )
  n <- n+1
}

l32_res_list_10 <- list()
n <- 1
set.seed(10)
for(i in 1:20){
  l32_res_list_10[[n]] <- run_sim3_test(
    sequencing_sums = seq_sum
    , num_indivs = 70
    , num_bases = 1024
    , effect_length = 32
    , effect_interval = effect_interval_32
    , effect_size_data = waveqtl_hmt_geno11
    , use_qt_data = FALSE
    , over_disp_param = 1000
    , Wmat_1024 = Wmat_1024
    , W2mat_1024 = W2mat_1024
    , library.read.depth = library.read.depth
    , Covariates = Covariates
    , group_data = group_data
    , effect_size = 1.4
    , trials_multiple = 10
  )
  n <- n+1
}

l64_res_list_10 <- list()
n <- 1
set.seed(10)
for(i in 1:20){
  l64_res_list_10[[n]] <- run_sim3_test(
    sequencing_sums = seq_sum
    , num_indivs = 70
    , num_bases = 1024
    , effect_length = 64
    , effect_interval = effect_interval_64
    , effect_size_data = waveqtl_hmt_geno11
    , use_qt_data = FALSE
    , over_disp_param = 1000
    , Wmat_1024 = Wmat_1024
    , W2mat_1024 = W2mat_1024
    , library.read.depth = library.read.depth
    , Covariates = Covariates
    , group_data = group_data
    , effect_size = 1.4
    , trials_multiple = 10
  )
  n <- n+1
}

extract_tps_fps <- function(output_list, effect_int){
  alt_analysis_findings <- lapply(output_list,function(x){x[["alt_analysis"]][["col_posi"]]})
  alt_analysis_hmt_findings <- lapply(output_list,function(x){x[["alt_hmt_analysis"]][["col_posi"]]})

  # How many in the interval?
  alt_analysis_findings_tp <- unlist(lapply(lapply(alt_analysis_findings,function(x){intersect(x,effect_int)}), length))
  alt_analysis_hmt_findings_tp <- unlist(lapply(lapply(alt_analysis_hmt_findings,function(x){intersect(x,effect_int)}), length))
  # How many not in the interval?
  alt_analysis_findings_fp <- unlist(lapply(lapply(alt_analysis_findings,function(x){setdiff(x,effect_int)}), length))
  alt_analysis_hmt_findings_fp <- unlist(lapply(lapply(alt_analysis_hmt_findings,function(x){setdiff(x,effect_int)}), length))
  return(list(
    alt_analysis_findings_tp = alt_analysis_findings_tp
    ,alt_analysis_hmt_findings_tp = alt_analysis_hmt_findings_tp
    ,alt_analysis_findings_fp = alt_analysis_findings_fp
    ,alt_analysis_hmt_findings_fp = alt_analysis_hmt_findings_fp
  ))
}

l8_diags <- extract_tps_fps(l8_res_list_100, effect_interval_8)
l16_diags <- extract_tps_fps(l16_res_list_10, effect_interval_16)
l32_diags <- extract_tps_fps(l32_res_list_10, effect_interval_32)
l64_diags <- extract_tps_fps(l64_res_list_10, effect_interval_64)

l8_diags_df <- data.frame(l8_diags,length = 8)
l16_diags_df <- data.frame(l16_diags,length = 16)
l32_diags_df <- data.frame(l32_diags,length = 32)
l64_diags_df <- data.frame(l64_diags,length = 64)

library(reshape2)
summary_df <- rbind.data.frame(l8_diags_df,l16_diags_df)
summary_df <- rbind.data.frame(summary_df,l32_diags_df)
summary_df <- rbind.data.frame(summary_df,l64_diags_df)
names(summary_df) <- c("TruePos - no HMT","TruePos - HMT","FalsePos - no HMT","FalsePos - HMT","Length")

# Plot?
library(ggplot2)
bplot <- ggplot() +
  geom_boxplot(data = melt(summary_df,id.vars = "Length")
               ,aes(x = factor(Length), y = value, colour = factor(variable))) +
  ggtitle("20 simulations on each length, Over-disp: 1000, Effect: 2/1.4/1.4/1.4, TrialsX: 10") +
  xlab("Effect length") +
  ylab("Counts") +
  labs(colour = "Count")
ggsave(filename = "images/20190822_sim3_bplot.png", plot = bplot)

saveRDS(l8_res_list_100,"data/20190822_sim3_l8_od1000_ef2X_tri10X.RDS",compress=T)
saveRDS(l16_res_list_10,"data/20190822_sim3_l16_od1000_ef1-4X_tri10X.RDS",compress=T)
saveRDS(l32_res_list_10,"data/20190822_sim3_l32_od1000_ef1-4X_tri10X.RDS",compress=T)
saveRDS(l64_res_list_10,"data/20190822_sim3_l64_od1000_ef1-4X_tri10X.RDS",compress=T)



# Drop another assumption:
# Effect size is a multiple of the original (rather than just a square effect)

# 8 different effect windows, 5 trials each

l8_res_list <- list()
n <- 1
set.seed(10)
for(i in 1:8){
  effect_interval_8 <- effect_length_picker(int_table,8)
  for(j in 1:5){
    l8_res_list[[n]] <- run_sim3_test(
      sequencing_sums = seq_sum
      , num_indivs = 70
      , num_bases = 1024
      , effect_length = 8
      , effect_interval = effect_interval_8
      , effect_size_data = waveqtl_hmt_geno11
      , use_qt_data = FALSE
      , over_disp_param = 1000
      , Wmat_1024 = Wmat_1024
      , W2mat_1024 = W2mat_1024
      , library.read.depth = library.read.depth
      , Covariates = Covariates
      , group_data = group_data
      , effect_multiple = 1.5e8
      , trials_multiple = 10
    )
    n <- n+1
  }
}
saveRDS(l8_res_list,"data/20190823_sim3_l8mix_od1000_ef1-5e8X_tri10X.RDS", compress = T)
