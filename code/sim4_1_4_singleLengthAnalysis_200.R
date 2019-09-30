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
od <- 70
ef_mult_proc1 <- 1.6e8*seq(0.4,1.4,by = 0.1)
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
savename <- "v2_multiEff_tieg15"

# Computer 1 --------------------------------------------------------------

for(len in proc1_lengths){
  print(len)
  # tmp <- list()
  tmp <- vector("list",length(ef_mult_proc1))
  n <- 1
  for(ef_size in ef_mult_proc1){
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
      , rMarkdownMode = F
      , outputAlias = proc1_label
    )
    n <- n+1
  }
  assign(x = paste0("l",len,"_res")
         ,value = tmp)
  saveRDS(tmp, paste0("data/",today_dt,"_l",len,"_",savename,".RDS"),compress = T)
}

### Effect length 8
# Effect sizes: 1.6e8*seq(1.5,1.8,by = 0.1)
# savename <- "v2_multiEff_biggerEff"

# More variance: 40
# savename <- "v2_multiEff_od40"

# Computer 2 --------------------------------------------------------------

for(len in proc2_lengths){
  print(len)
  # tmp <- list()
  tmp <- vector("list",length(ef_mult_proc2))
  n <- 1
  for(ef_size in ef_mult_proc2){
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
      , rMarkdownMode = F
      , outputAlias = proc2_label
    )
    n <- n+1
  }
  assign(x = paste0("l",len,"_res")
         ,value = tmp)
  saveRDS(tmp, paste0("data/",today_dt,"_l",len,"_",savename,".RDS"),compress = T)
}

# Same length, different size ---------------------------------------------

l8_diag <- lapply(l8_res,waveqtl_diags,num_sims = 200)
l8_diag[[1]]
l8_diag[[2]]
l8_diag[[3]]
l8_diag[[4]]
l8_diag[[5]]
l8_diag[[6]]
l8_diag[[7]]
l8_diag[[8]]
l8_diag[[9]]
l8_diag[[10]]
l8_diag[[11]]

dt <- sapply(l8_diag,function(x){
  return(
    c(unlist(attr(x[['perf_nohmt']],"y.values"))
      ,unlist(attr(x[['perf_hmt']],"y.values"))
    )
  )})

y_range <- c(min(dt)*0.95,min(max(dt)*1.05,1))
plot(x = ef_mult_proc1
     ,y = dt[1,]
     ,col = "red", type = "o"
     ,xlab = "Effect Size"
     ,ylab = "auROC"
     ,ylim = y_range
     ,xaxt = "n"
     ,main = "Effect Length 8")
axis(side = 1
     ,at = ef_mult_proc1
     ,labels = ef_mult_proc1
     ,tick = ef_mult_proc1)
lines(x = ef_mult_proc1
      ,y = dt[2,], col = "blue"
      ,type = "o")
