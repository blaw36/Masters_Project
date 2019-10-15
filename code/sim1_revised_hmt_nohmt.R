source("code/sim1_revised_hmt_nohmt_functions.R")

# Run sims ----------------------------------------------------------------

# Read in Wmat
Wmat_1024 = read.table("~/Cpp/WaveQTL_HMT/data/DWT/Wmat_1024",as.is = TRUE)
W2mat_1024 = Wmat_1024*Wmat_1024

# WaveQTL
system.time({
  wqtl_samples <- draw_samples_no_hmt(
    data_path = "~/Cpp/WaveQTL/test/dsQTL/output/"
    ,dataset = "test.no.QT"
    ,geno_select = 11
    ,wmat = Wmat_1024
    ,seed_no = 10
    ,num_samples = 1000
    ,num_pheno = 1024
    ,num_indiv = 70
  )
})

# WaveQTL_HMT
system.time({
  wqtl_hmt_samples <- draw_samples_hmt(
    data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/"
    ,dataset = "tree_tie_noQT"
    ,waveqtl_data_path = "~/Cpp/WaveQTL_HMT/test/dsQTL/output/WaveQTL/"
    ,waveqtl_dataset = "test.no.QT"
    ,geno_select = 11
    ,wmat = Wmat_1024
    ,seed_no = 10
    ,num_samples = 1000
    ,num_pheno = 1024
    ,num_indiv = 70
  )
})

### Plots

# +- 3 SD
wqtl_mean <- apply(wqtl_samples,MARGIN = 2,mean)
wqtl_sd <- apply(wqtl_samples,MARGIN = 2,sd)
wqtl_lwr <- apply(wqtl_samples,MARGIN = 2,quantile,probs = 0.005)
wqtl_upr <- apply(wqtl_samples,MARGIN = 2,quantile,probs = 0.995)

wqtl_hmt_mean <- apply(wqtl_hmt_samples,MARGIN = 2,mean)
wqtl_hmt_sd <- apply(wqtl_hmt_samples,MARGIN = 2,sd)
wqtl_hmt_lwr <- apply(wqtl_hmt_samples,MARGIN = 2,quantile,probs = 0.005)
wqtl_hmt_upr <- apply(wqtl_hmt_samples,MARGIN = 2,quantile,probs = 0.995)

# Wqtl
ymin_beta = min(wqtl_mean - 3*wqtl_sd) - abs(min(wqtl_mean - 3*wqtl_sd))*0.0000000001
ymax_beta = max(wqtl_mean + 3*wqtl_sd) + abs(max(wqtl_mean + 3*wqtl_sd))*0.0000000001

beta_l = wqtl_mean - 3*wqtl_sd
beta_r = wqtl_mean + 3*wqtl_sd

wh_l = which(beta_l > 0)
wh_r = which(beta_r < 0)
high_wh = sort(unique(union(wh_l, wh_r)))

xval = 1:1024
col_posi = xval[high_wh]

par(mar = c(2,4,4,2))
plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024),main ="Wqtl - +- 3 sd method", axes=FALSE)
axis(2)
if(length(col_posi) > 0){
  for(j in 1:length(col_posi)){
    polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin_beta-2, ymax_beta+2, ymax_beta+2, ymin_beta-2), col ="pink", border = NA)
  }
}

abline(h = 0, col = "red")
points(xval, wqtl_mean, col = "blue", type="l")
points(xval, beta_l, col = "skyblue", type="l")
points(xval, beta_r, col = "skyblue", type="l")
box()

# Wqtl_HMT
ymin_beta = min(wqtl_hmt_mean - 3*wqtl_hmt_sd) - abs(min(wqtl_hmt_mean - 3*wqtl_hmt_sd))*0.0000000001
ymax_beta = max(wqtl_hmt_mean + 3*wqtl_hmt_sd) + abs(max(wqtl_hmt_mean + 3*wqtl_hmt_sd))*0.0000000001

beta_l = wqtl_hmt_mean - 3*wqtl_hmt_sd
beta_r = wqtl_hmt_mean + 3*wqtl_hmt_sd

wh_l = which(beta_l > 0)
wh_r = which(beta_r < 0)
high_wh = sort(unique(union(wh_l, wh_r)))

xval = 1:1024
col_posi = xval[high_wh]

par(mar = c(2,4,4,2))
plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024),main ="Wqtl HMT - +- 3 sd method", axes=FALSE)
axis(2)
if(length(col_posi) > 0){
  for(j in 1:length(col_posi)){
    polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin_beta-2, ymax_beta+2, ymax_beta+2, ymin_beta-2), col ="pink", border = NA)
  }
}

abline(h = 0, col = "red")
points(xval, wqtl_hmt_mean, col = "blue", type="l")
points(xval, beta_l, col = "skyblue", type="l")
points(xval, beta_r, col = "skyblue", type="l")
box()


# Posterior sampling
wqtl_mean <- apply(wqtl_samples,MARGIN = 2,mean)
wqtl_sd <- apply(wqtl_samples,MARGIN = 2,sd)
wqtl_lwr <- apply(wqtl_samples,MARGIN = 2,quantile,probs = 0.005)
wqtl_upr <- apply(wqtl_samples,MARGIN = 2,quantile,probs = 0.995)

wqtl_hmt_mean <- apply(wqtl_hmt_samples,MARGIN = 2,mean)
wqtl_hmt_sd <- apply(wqtl_hmt_samples,MARGIN = 2,sd)
wqtl_hmt_lwr <- apply(wqtl_hmt_samples,MARGIN = 2,quantile,probs = 0.005)
wqtl_hmt_upr <- apply(wqtl_hmt_samples,MARGIN = 2,quantile,probs = 0.995)

# Wqtl
ymin_beta = min(wqtl_mean - 3*wqtl_sd) - abs(min(wqtl_mean - 3*wqtl_sd))*0.0000000001
ymax_beta = max(wqtl_mean + 3*wqtl_sd) + abs(max(wqtl_mean + 3*wqtl_sd))*0.0000000001

beta_l = wqtl_mean - 3*wqtl_sd
beta_r = wqtl_mean + 3*wqtl_sd

wh_l = which(beta_l > 0)
wh_r = which(beta_r < 0)
high_wh = sort(unique(union(wh_l, wh_r)))

xval = 1:1024
col_posi = xval[high_wh]

par(mar = c(2,4,4,2))
plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024),main ="Wqtl - posterior sampling method", axes=FALSE)
axis(2)
if(length(col_posi) > 0){
  for(j in 1:length(col_posi)){
    polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin_beta-2, ymax_beta+2, ymax_beta+2, ymin_beta-2), col ="pink", border = NA)
  }
}

abline(h = 0, col = "red")
points(xval, wqtl_mean, col = "blue", type="l")
points(xval, beta_l, col = "skyblue", type="l")
points(xval, beta_r, col = "skyblue", type="l")
box()

# Wqtl_HMT
ymin_beta = min(wqtl_hmt_mean - 3*wqtl_hmt_sd) - abs(min(wqtl_hmt_mean - 3*wqtl_hmt_sd))*0.0000000001
ymax_beta = max(wqtl_hmt_mean + 3*wqtl_hmt_sd) + abs(max(wqtl_hmt_mean + 3*wqtl_hmt_sd))*0.0000000001

beta_l = wqtl_hmt_mean - 3*wqtl_hmt_sd
beta_r = wqtl_hmt_mean + 3*wqtl_hmt_sd

wh_l = which(beta_l > 0)
wh_r = which(beta_r < 0)
high_wh = sort(unique(union(wh_l, wh_r)))

xval = 1:1024
col_posi = xval[high_wh]

par(mar = c(2,4,4,2))
plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024),main ="Wqtl HMT - posterior sampling method", axes=FALSE)
axis(2)
if(length(col_posi) > 0){
  for(j in 1:length(col_posi)){
    polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin_beta-2, ymax_beta+2, ymax_beta+2, ymin_beta-2), col ="pink", border = NA)
  }
}

abline(h = 0, col = "red")
points(xval, wqtl_hmt_mean, col = "blue", type="l")
points(xval, beta_l, col = "skyblue", type="l")
points(xval, beta_r, col = "skyblue", type="l")
box()

# Closed form -------------------------------------------------------------

# WaveQTL
system.time({
  wqtl_samples_closed <- closed_form_no_hmt(
    data_path = "~/Cpp/WaveQTL/test/dsQTL/output/"
    ,dataset = "test.no.QT"
    ,geno_select = 11
    ,wmat = Wmat_1024
    ,wmat2 = W2mat_1024
    ,num_pheno = 1024
  )
})

# WaveQTL_HMT

### Plots

# +- 3 SD

# Wqtl
ymin_beta = min(wqtl_samples_closed$beta_mean_dataS - 3*wqtl_samples_closed$beta_sd_dataS) - abs(min(wqtl_samples_closed$beta_mean_dataS - 3*wqtl_samples_closed$beta_sd_dataS))*0.0000000001
ymax_beta = max(wqtl_samples_closed$beta_mean_dataS + 3*wqtl_samples_closed$beta_sd_dataS) + abs(max(wqtl_samples_closed$beta_mean_dataS + 3*wqtl_samples_closed$beta_sd_dataS))*0.0000000001

beta_l = wqtl_samples_closed$beta_mean_dataS - 3*wqtl_samples_closed$beta_sd_dataS
beta_r = wqtl_samples_closed$beta_mean_dataS + 3*wqtl_samples_closed$beta_sd_dataS

wh_l = which(beta_l > 0)
wh_r = which(beta_r < 0)
high_wh = sort(unique(union(wh_l, wh_r)))

xval = 1:1024
col_posi = xval[high_wh]

par(mar = c(2,4,4,2))
plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024),main ="Wqtl - +- 3 sd method", axes=FALSE)
axis(2)
if(length(col_posi) > 0){
  for(j in 1:length(col_posi)){
    polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin_beta-2, ymax_beta+2, ymax_beta+2, ymin_beta-2), col ="pink", border = NA)
  }
}

abline(h = 0, col = "red")
points(xval, wqtl_samples_closed$beta_mean_dataS, col = "blue", type="l")
points(xval, beta_l, col = "skyblue", type="l")
points(xval, beta_r, col = "skyblue", type="l")
box()

# Wqtl HMT
