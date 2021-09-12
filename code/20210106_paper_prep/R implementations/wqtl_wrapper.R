rm(list = ls()); gc();

source("~/Cpp/20210817_Wqtl/R/R_Wqtl_Implementations/wqtl_hmt_func.R")
source("~/Cpp/20210817_Wqtl/R/R_Wqtl_Implementations/wqtl_func.R")
library(data.table)


# HMT, no sc --------------------------------------------------------------

a <- as.matrix(read.table("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/chr.2.696.hmt.fph.logLR.txt"))
logBFs = log(10)*as.numeric(a[6,4:1026]) # Exclude scaling coeff
groups = floor(log2((1:length(logBFs))))+1

c2_s696_hmt <- wqtl_hmt(
  logBFs = logBFs,
  groups = groups,
  tying_groups = c(1,2,4,8,16,32,64,128,256,512), # NULL for no tying
  iterations_max = 1000,
  conv_tol = 0.005,
  sumlog_use = T,
  init_pi = 0.5,
  init_eps_11 = 0.5,
  init_eps_10 = 0.5
) # LogL: 0.5865738

## Try and incorporate scale 0 for WQtl
a <- as.matrix(read.table("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/chr.2.696.hmt.fph.logLR.txt"))
logBFs = log(10)*as.numeric(a[6,3:1026]) # Include scaling coeff


c2_s696_nohmt <- wqtl(
  logBFs = logBFs,
  iterations_max = 1000,
  conv_tol = 0.005,
  init_pi = 0.5
) # 6.51
## This is VERY close. The Pi's suggest we're slightly off.
a_nohmt <- as.matrix(read.table("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/chr.2.696.nohmt.fph.logLR.txt"))
a_nohmt_pi <- as.matrix(read.table("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/chr.2.696.nohmt.fph.pi.txt"))

# Comparing LogLs
as.numeric(a_nohmt[6,2])
sum(last(c2_s696_nohmt$logl_list),na.rm = T)

# Comparing Pis
round(as.numeric(a_nohmt_pi[6,-1]),5)
round(last(c2_s696_nohmt$tmp_pi),5)

# Comparing SC LogL from WQtl output vs from the R script
# Using R's pi
sc_pi = last(c2_s696_nohmt$tmp_pi)[1]
sc_bf = 10^as.numeric(a_nohmt[6,3])
log(sc_pi*sc_bf + (1-sc_pi))
last(c2_s696_nohmt$logl_list)[1]

# Using WQtl's pi
sc_pi = as.numeric(a_nohmt_pi[6,2])
sc_bf = 10^as.numeric(a_nohmt[6,3])
log(sc_pi*sc_bf + (1-sc_pi))
last(c2_s696_nohmt$logl_list)[1]
# We get the same calculation both times. The R output logL a little higher (rounding?)

## Therefore, the HMT Likelihood should be:
c2_s696_hmt$logLR + log(sc_pi*sc_bf + (1-sc_pi)) # 6.628281, OR
c2_s696_hmt$logLR + last(c2_s696_nohmt$logl_list)[1] # 6.628275


# Check SC LogL from algo vs outputs --------------------------------------

chr_sites <- list.files("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/", pattern = "*.nohmt.fph.logLR.txt")
chr_sites_pi <- list.files("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/", pattern = "*.nohmt.fph.pi.txt")
for(i in 1:length(chr_sites)){
# for(i in 5:5){
  print(chr_sites[i])
  bfs <- as.matrix(read.table(paste0("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/", chr_sites[i])))
  pis <- as.matrix(read.table(paste0("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/", chr_sites_pi[i])))
  # Pick highest LogLR (and min row in case of tiebreak) to do
  chr_max = data.table(bfs[,1:2])[,min(which(V2==max(V2)))]
  logBFs = log(10)*as.numeric(bfs[chr_max,3:1026]) # Include scaling coeff
  
  wqtl_nohmt <- wqtl(
    logBFs = logBFs,
    iterations_max = 1000,
    conv_tol = 0.005,
    init_pi = 0.5
  )

  cat("C++ logL:", as.numeric(bfs[chr_max,2]),"\n")
  cat("R LogL:",sum(last(wqtl_nohmt$logl_list),na.rm = T),"\n")
  
  cat("C++ pi:", round(as.numeric(pis[chr_max,-1]),5),"\n")
  cat("R pi:", round(last(wqtl_nohmt$tmp_pi),5),"\n")
  
  sc_pi = as.numeric(pis[chr_max,2])
  sc_bf = 10^(as.numeric(bfs[chr_max,3]))
  cat("C++ SC LogL: ",log(sc_pi*sc_bf + (1-sc_pi)),"\n")
  cat("R SC LogL: ",last(wqtl_nohmt$logl_list)[1],"\n")
}


# Check WQtl vs WQtl_HMT with SC ------------------------------------------

source("~/Cpp/20210817_Wqtl/R/R_Wqtl_Implementations/wqtl_hmt_func_v2.R")

chr_sites <- list.files("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/", pattern = "*.nohmt.fph.logLR.txt")
chr_sites_pi <- list.files("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/", pattern = "*.nohmt.fph.pi.txt")
for(i in 1:length(chr_sites)){
  # for(i in 5:5){
  print(chr_sites[i])
  bfs <- as.matrix(read.table(paste0("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/", chr_sites[i])))
  pis <- as.matrix(read.table(paste0("~/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/", chr_sites_pi[i])))
  # Pick highest LogLR (and min row in case of tiebreak) to do
  chr_max = data.table(bfs[,1:2])[,min(which(V2==max(V2)))]
  logBFs = log(10)*as.numeric(bfs[chr_max,3:1026]) # Include scaling coeff
  
  # WQtl
  wqtl_nohmt <- wqtl(
    logBFs = logBFs,
    iterations_max = 1000,
    conv_tol = 0.005,
    init_pi = 0.5
  )
  
  # HMT, no SC
  logBFs_no_sc = log(10)*as.numeric(bfs[chr_max,4:1026]) # Exclude scaling coeff
  groups = floor(log2((1:length(logBFs_no_sc))))+1
  wqtl_hmt_no_sc <- wqtl_hmt(
    logBFs = logBFs_no_sc,
    groups = groups,
    tying_groups = c(1,2,4,8,16,32,64,128,256,512), # NULL for no tying
    iterations_max = 1000,
    conv_tol = 0.005,
    sumlog_use = T,
    init_pi = 0.5,
    init_eps_11 = 0.5,
    init_eps_10 = 0.5
  )
  
  # HMT
  wqtl_hmt_with_sc <- wqtl_hmt_v2(
    logBFs = logBFs,
    iterations_max = 1000,
    conv_tol = 0.005,
    sumlog_use = T,
    init_pi = 0.5,
    init_eps_11 = 0.5,
    init_eps_10 = 0.5
  )
  
  cat("WQTL Diagnostics: ---------------------------\n")
  cat("C++ logL:", as.numeric(bfs[chr_max,2]),"\n")
  cat("R LogL:",sum(last(wqtl_nohmt$logl_list),na.rm = T),"\n")
  
  cat("C++ pi:", round(as.numeric(pis[chr_max,-1]),5),"\n")
  cat("R pi:", round(last(wqtl_nohmt$tmp_pi),5),"\n")
  
  sc_pi = as.numeric(pis[chr_max,2])
  sc_bf = 10^(as.numeric(bfs[chr_max,3]))
  cat("C++ SC LogL: ",log(sc_pi*sc_bf + (1-sc_pi)),"\n")
  cat("R SC LogL: ",last(wqtl_nohmt$logl_list)[1],"\n")
  cat("R SC logL from HMT: ", last(wqtl_hmt_with_sc$logl_list), "\n")
  
  cat("HMT Diagnostics: ---------------------------\n")
  cat("HMT with SC LogL: ", wqtl_hmt_with_sc$logLR, "\n")
  cat("HMT without SC LogL: ", wqtl_hmt_no_sc$logLR, "\n")
  cat("HMT without SC LogL plus SC LogL: ", wqtl_hmt_no_sc$logLR + last(wqtl_hmt_with_sc$logl_list), "\n")
}
