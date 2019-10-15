### Script runs:
# For a single effect length, tries to find area where the effects are detectable
# and where they fail to become detectable (by HMT)

# Initialisation ----------------------------------------------------------

input_data_path = "~/Cpp/WaveQTL_HMT/data/dsQTL/"
data_path <- "~/Cpp/WaveQTL_HMT/test/dsQTL/output/" # change this when you port it all over to the Masters Git repo
dataset <- "tree_tie_noQT"
waveqtl_data_path <- "~/Cpp/WaveQTL_HMT/test/dsQTL/output/WaveQTL/"
waveqtl_dataset <- "test.no.QT"
geno_select <- 11 # the one used in the demo

library(dplyr)
library(ROCR)
source("code/sim3_functions.R")
source("code/sim4_functions.R")

effect_size_and_data <- read_in_gen_eff_size(geno_select = 11)

# Pick effect sizes:
effect_interval_8 <- 519:526
effect_interval_16 <- 450:465
effect_interval_32 <- 437:468
effect_interval_64 <- 417:480

effect_interval_24 <- 442:465
effect_interval_40 <- 433:472
effect_interval_48 <- 500:547
effect_interval_56 <- 417:472

# Params ------------------------------------------------------------------

# ef_mult <- seq(1e8,2e8,2e7)
# num_sims <- 50
# od <- 1/70
# od <- c(1/70,1,30,70,700,700000)
od <- c(1,seq(10,100,by = 20))


# ef_mult_proc1 <- 1.6e8*seq(0.4,1.4,by = 0.1)
ef_mult_proc1 <- 1.6e8
ef_mult_proc2 <- 8e7*seq(0.4,1.4,by = 0.1)
# Try: seq(8e6,2.88e8,by = 2e7)
num_sims <- 200

proc1_label <- "sim3"
proc2_label <- "sim4"

# proc1_lengths <- seq(8,32,8)
# proc2_lengths <- seq(40,64,8)
proc1_lengths <- 8
proc2_lengths <- 16

# 8: 1.6e8
# 16: 8e7
# 32: 4e7
# 48: 3e7
# 64: 2e7

today_dt <- format(Sys.Date(),"%Y%m%d")
savename <- "v2_multiOD"

# Computer 1 --------------------------------------------------------------

for(len in proc1_lengths){
  print(len)
  # tmp <- list()
  tmp <- vector("list",length(od))
  n <- 1
  for(od_spec in od){
    print(od_spec)
    set.seed(1)
    tmp[[n]] <- run_sim4_v2(
      sequencing_sums = effect_size_and_data$seq_sum
      , num_indivs = 70
      , num_bases = 1024
      , effect_length = len
      , effect_interval = get(paste0("effect_interval_",len))
      , effect_size_data = effect_size_and_data$effect_size
      , no.QT = FALSE
      , over_disp_param = od_spec
      , Wmat_1024 = effect_size_and_data$Wmat_1024
      , W2mat_1024 = effect_size_and_data$W2mat_1024
      , library.read.depth = effect_size_and_data$library.read.depth
      , Covariates = effect_size_and_data$Covariates
      , effect_multiple = ef_mult_proc1
      , trials_multiple = 10
      , number_sims = num_sims
      , verbose = F
      , rMarkdownMode = F
      , outputAlias = proc1_label
    )
    n <- n+1
  }
  assign(x = paste0("l",len,"_res")
         ,value = tmp)
  saveRDS(tmp, paste0("data/",today_dt,"_l",len,"_",savename,".RDS"),compress = T)
}

# Then, try 700000 (almost binomial sampling) where
# effect size is half of what it is now


# Diags -------------------------------------------------------------------

l8_diag <- lapply(l8_res,waveqtl_diags,num_sims = num_sims)
dt_8 <- sapply(l8_diag,function(x){
  return(
    c(unlist(attr(x[['perf_nohmt']],"y.values"))
      ,unlist(attr(x[['perf_hmt']],"y.values"))
    )
  )})

ef_size_8 <- sapply(l8_res, function(x){x[["params_list"]][["effect_multiple"]]})

# Indiv effect lnegths ----------------------------------------------------

# for(i in c(8,16,32,64)){
for(i in c(8)){
  dt_effects <- get(paste0("dt_",i))
  eff_sizes <- get(paste0("ef_size_",i))
  y_range <- c(min(dt_effects)*0.95,min(max(dt_effects)*1.05,1))
  # plot.new()
  # pdf(paste0("images/ch3_l",i,"_multieff_auROC.pdf")
  #     , width = 10, height = 6, pointsize = 14)
  plot(x = od
       ,y = dt_effects[1,]
       ,col = "red", type = "o", pch = 20, cex = 0.6
       ,xlab = "Effect Size"
       ,ylab = "auROC"
       ,ylim = y_range
       ,xaxt = "n"
       ,main = paste0("Effect Length ",i))
  axis(side = 1
       ,at = od
       ,labels = od
       ,tick = od)
  lines(x = od
        ,y = dt_effects[2,], col = "blue"
        ,type = "o", pch = 20, cex = 0.6)
  legend("bottomright",legend = c("no-HMT","HMT")
         ,col=c("red","blue"),lty=rep(1,2),pch = c(20,20), cex = c(1.2,1.2))
  # dev.off()
}


# Indiv lengths - ROCs ----------------------------------------------------

extract_roc_plot_data <- function(sims_list,num_sims){
  preds_nohmt <- matrix(c(sims_list$null_waveqtl_lhood
                          ,sims_list$alt_waveqtl_lhood
                          ,rep(0,num_sims),rep(1,num_sims))
                        ,nrow = 2*num_sims,ncol = 2,byrow = F)
  preds_hmt <- matrix(c(sims_list$null_waveqtl_hmt_lhood
                        ,sims_list$alt_waveqtl_hmt_lhood
                        ,rep(0,num_sims),rep(1,num_sims))
                      ,nrow = 2*num_sims,ncol = 2,byrow = F)

  perf_nohmt <- performance(prediction(preds_nohmt[,1],preds_nohmt[,2]),measure = "auc")
  perf_hmt <- performance(prediction(preds_hmt[,1],preds_hmt[,2]),measure = "auc")

  n = dim(preds_hmt)[1]

  O1 = order(preds_nohmt[,1], decreasing =TRUE)
  O2 = order(preds_hmt[,1], decreasing =TRUE)

  C1  = c(0,cumsum(preds_nohmt[O1,2])) / sum(preds_nohmt[,2]) #True positives proportions or sensitivity
  C2  = c(0,cumsum(preds_hmt[O2,2])) / sum(preds_hmt[,2])

  FP1 = c(0,cumsum(1-preds_nohmt[O1,2])) / (n-sum(preds_nohmt[,2])) #false positives proportions based on Model 1.
  FP2 = c(0,cumsum(1-preds_hmt[O2,2])) / (n-sum(preds_hmt[,2]))

  return(list(
    FP1 = FP1
    ,C1 = C1
    ,FP2 = FP2
    ,C2 = C2
    ,perf_nohmt = unlist(attr(perf_nohmt,"y.values"))
    ,perf_hmt = unlist(attr(perf_hmt,"y.values"))
  ))
}

l8_roc_data <- lapply(l8_res,extract_roc_plot_data,num_sims = num_sims)

for(i in c(8)){
# for(i in c(8,16,32,64)){
  # for(i in c(8,16)){
  roc_data <- get(paste0("l",i,"_roc_data"))
  # indices_to_use <- c(1,seq(3,length(roc_data),by = 3))
  indices_to_use <- 1:length(roc_data)
  # eff_sizes <- get(paste0("ef_size_",i))
  od_sizes <- od

  # plot.new()
  # pdf(paste0("images/ch3_l",i,"_multieff_ROC.pdf")
  #     , width = 15, height = 9, pointsize = 20)

  par(xpd = T, mar = par()$mar + c(0,0,0,9))
  plot(roc_data[[1]]$FP1,roc_data[[1]]$C1, type="l", lwd=1.2, col="black", lty = 3, xlab="False Positive rate", ylab="True Positive Rate"
       , main = paste0("Effect length ",i))
  lines(roc_data[[1]]$FP2,roc_data[[1]]$C2, lty = 1, lwd=1.2, col = "black")

  n <- 2
  for(j in indices_to_use[2:length(indices_to_use)]){
    lines(roc_data[[j]]$FP1,roc_data[[j]]$C1, lty = 3, lwd=1.2, col = n)
    lines(roc_data[[j]]$FP2,roc_data[[j]]$C2, lty = 1, lwd=1.2, col = n)
    n <- n+1
  }
  lines(x = seq(0,1,length.out = 10000)
        ,y= seq(0,1,length.out = 10000),col = "grey60")

  legend_string <- c()
  for(k in indices_to_use){
    legend_string <- c(legend_string,
                       paste0("no-HMT: ",od_sizes[[k]],": ",round(roc_data[[k]]$perf_nohmt,4)))
    legend_string <- c(legend_string,
                       paste0("HMT: ",od_sizes[[k]],": ",round(roc_data[[k]]$perf_hmt,4)))
  }
  legend(1.05,1
         ,legend = legend_string
         ,cex = 0.8
         ,bg = "transparent"
         ,bty = "n"
         ,lty=rep(c(3,1),length(indices_to_use))
         ,lwd=rep(c(1.2,1.2),length(indices_to_use))
         ,col=rep(c(1:length(indices_to_use)),each = 2))
  par(mar=c(5, 4, 4, 2) + 0.1)
  # dev.off()
}
