# Function which does effect size plots:
# Inputs: No QT mean and var outputs, WMat, QT pval output (to get strongest SNP index)
# Outputs: plot (written to disk), constitutents (col_posi, beta_l, beta_r, beta_dataS)

effect_size_plot <- function(
  site
  ,chrIX
  ,wqtl_qt_output_path # ie. "/data/projects/punim0614/brendan/WaveQTL_HMT_wperm/20210627_run4_smller_rndm_batch/output/"
  ,wqtl_no_qt_output_path # may be same place as above; where noQT outputs are stored
  ,plot_output_path # ie. "/data/gpfs/projects/punim0614/brendan/WaveQTL_HMT_wperm/20210627_run4_smller_rndm_batch/plots/"
  ,snp_to_use = 'nearest' # other options; 'hmt', 'nohmt'
){

  # source("code/sim1_revised_hmt_nohmt_functions.R")
  source("/home/bklaw/WaveQTL_HMT_wperm/R/sim1_revised_hmt_nohmt_functions.R")

  ## Read DWT matrix
  Wmat_1024 = read.table("/home/bklaw/WaveQTL/WaveQTL-master/data/DWT/Wmat_1024",as.is = TRUE)
  W2mat_1024 = Wmat_1024*Wmat_1024

  # Determine which SNPs to plot
  if(snp_to_use == 'nearest'){
    ## Grab index of SNP with strongest effect by reading in pval file
    nearby_snp_nohmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".nohmt.fph.pval.txt"), stringsAsFactors = F)$V1[1]
    ## Grab index of SNP with strongest effect by reading in pval file
    nearby_snp_hmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".hmt.fph.pval.txt"), stringsAsFactors = F)$V1[1]
  }else if(snp_to_use == 'nohmt'){
    nearby_snp_nohmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".nohmt.fph.pval.txt"), stringsAsFactors = F)$V1[1]
    nearby_snp_hmt = nearby_snp_nohmt
  }else if(snp_to_use == 'hmt'){
    nearby_snp_hmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".hmt.fph.pval.txt"), stringsAsFactors = F)$V1[1]
    nearby_snp_nohmt = nearby_snp_hmt
  }

  # WaveQTL samples
  pval_nohmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".nohmt.fph.pval.txt"), stringsAsFactors = F)$V1[4]
  test = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".nohmt.fph.mean.txt"), stringsAsFactors = F) # Read in one of the mean or var files we'll be using here
  sel_geno_IX = which(nearby_snp_nohmt == test[,1]) # test[,1] should relate to the row names, ie the SNP names; we want to grab the row of the strongest SNP.
  wqtl_samples <- draw_samples_no_hmt(
    data_path = wqtl_no_qt_output_path
    ,dataset = paste0("chr.",chrIX,".",site,".nohmt.no.QT")
    ,geno_select = sel_geno_IX
    ,wmat = Wmat_1024
    ,seed_no = 10
    ,num_samples = 1000
    ,num_pheno = 1024
    ,num_indiv = 70
  )

  # WaveQTL_HMT samples
  pval_hmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".hmt.fph.pval.txt"), stringsAsFactors = F)$V1[4]
  test = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".hmt.fph.mean1.txt"), stringsAsFactors = F) # Read in one of the mean or var files we'll be using here
  sel_geno_IX = which(nearby_snp_hmt == test[,1]) # test[,1] should relate to the row names, ie the SNP names; we want to grab the row of the strongest SNP.
  wqtl_hmt_samples <- draw_samples_hmt(
    data_path = wqtl_no_qt_output_path
    ,dataset = paste0("chr.",chrIX,".",site,".hmt.no.QT")
    ,waveqtl_data_path = wqtl_no_qt_output_path
    ,waveqtl_dataset = paste0("chr.",chrIX,".",site,".nohmt.no.QT")
    ,geno_select = sel_geno_IX
    ,wmat = Wmat_1024
    ,seed_no = 10
    ,num_samples = 1000
    ,num_pheno = 1024
    ,num_indiv = 70
  )

  wqtl_mean <- apply(wqtl_samples,MARGIN = 2,mean)
  wqtl_sd <- apply(wqtl_samples,MARGIN = 2,sd)

  wqtl_hmt_mean <- apply(wqtl_hmt_samples,MARGIN = 2,mean)
  wqtl_hmt_sd <- apply(wqtl_hmt_samples,MARGIN = 2,sd)

  ymin_beta_wqtl = min(wqtl_mean - 3*wqtl_sd) - abs(min(wqtl_mean - 3*wqtl_sd))*0.0000000001
  ymax_beta_wqtl = max(wqtl_mean + 3*wqtl_sd) + abs(max(wqtl_mean + 3*wqtl_sd))*0.0000000001
  ymin_beta_wqtlhmt = min(wqtl_hmt_mean - 3*wqtl_hmt_sd) - abs(min(wqtl_hmt_mean - 3*wqtl_hmt_sd))*0.0000000001
  ymax_beta_wqtlhmt = max(wqtl_hmt_mean + 3*wqtl_hmt_sd) + abs(max(wqtl_hmt_mean + 3*wqtl_hmt_sd))*0.0000000001

  # To keep the axes common
  ymin_beta = min(ymin_beta_wqtl,ymin_beta_wqtlhmt)
  ymax_beta = max(ymax_beta_wqtl,ymax_beta_wqtlhmt)


  # WaveQTL Plot

  # Wqtl
  # ymin_beta = min(wqtl_mean - 3*wqtl_sd) - abs(min(wqtl_mean - 3*wqtl_sd))*0.0000000001
  # ymax_beta = max(wqtl_mean + 3*wqtl_sd) + abs(max(wqtl_mean + 3*wqtl_sd))*0.0000000001

  beta_l = wqtl_mean - 3*wqtl_sd
  beta_r = wqtl_mean + 3*wqtl_sd

  wh_l = which(beta_l > 0)
  wh_r = which(beta_r < 0)
  high_wh = sort(unique(union(wh_l, wh_r)))

  xval = 1:1024
  wqtl_col_posi = xval[high_wh]

  # Plot titles
  if(snp_to_use == 'nearest'){
    nohmt_title = paste0("No HMT - Chr: ",chrIX,", Site: ",site,", Nearby SNP: ",nearby_snp_nohmt,", pval: ",round(as.numeric(pval_nohmt),4))
    hmt_title = paste0("HMT - Chr: ",chrIX,", Site: ",site,", Nearby SNP: ",nearby_snp_hmt,", pval: ",round(as.numeric(pval_hmt),4))
  }else if(snp_to_use == 'nohmt'){
    nohmt_title = paste0("No HMT - Chr: ",chrIX,", Site: ",site,", Nearby SNP: ",nearby_snp_nohmt,", pval: ",round(as.numeric(pval_nohmt),4))
    hmt_title = paste0("HMT - Chr: ",chrIX,", Site: ",site,", Nearby SNP: ",nearby_snp_nohmt)
  }else if(snp_to_use == 'hmt'){
    hmt_title = paste0("HMT - Chr: ",chrIX,", Site: ",site,", Nearby SNP: ",nearby_snp_hmt,", pval: ",round(as.numeric(pval_hmt),4))
    nohmt_title = paste0("No HMT - Chr: ",chrIX,", Site: ",site,", Nearby SNP: ",nearby_snp_hmt)
  }

  # No HMT plot
  plot.new()
  pdf(paste0(plot_output_path,"chr.",chrIX,".",site,".nohmt.effect.plot.PDF"), width = 8, height = 3)
  par(mar = c(2,4,4,2))
  plot(1,1,type="n", xlab = "Base location", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024)
       ,main = nohmt_title, axes=FALSE)

  axis(2)
  axis(1,at = c(1,seq(256,1024,by = 256)),tick = c(1,seq(256,1024,by = 256)))
  if(length(wqtl_col_posi) > 0){
    for(j in 1:length(wqtl_col_posi)){
      polygon(c(wqtl_col_posi[j]-0.5, wqtl_col_posi[j]-0.5, wqtl_col_posi[j]+0.5, wqtl_col_posi[j]+0.5), c(ymin_beta-2, ymax_beta+2, ymax_beta+2, ymin_beta-2), col ="pink", border = NA)
    }
  }

  abline(h = 0, col = "red")
  points(xval, wqtl_mean, col = "blue", type="l")
  points(xval, beta_l, col = "skyblue", type="l")
  points(xval, beta_r, col = "skyblue", type="l")
  box()

  dev.off()

  # WaveQTL_HMT plot

  # Wqtl_HMT
  # ymin_beta = min(wqtl_hmt_mean - 3*wqtl_hmt_sd) - abs(min(wqtl_hmt_mean - 3*wqtl_hmt_sd))*0.0000000001
  # ymax_beta = max(wqtl_hmt_mean + 3*wqtl_hmt_sd) + abs(max(wqtl_hmt_mean + 3*wqtl_hmt_sd))*0.0000000001

  beta_l = wqtl_hmt_mean - 3*wqtl_hmt_sd
  beta_r = wqtl_hmt_mean + 3*wqtl_hmt_sd

  wh_l = which(beta_l > 0)
  wh_r = which(beta_r < 0)
  high_wh = sort(unique(union(wh_l, wh_r)))

  xval = 1:1024
  wqtl_hmt_col_posi = xval[high_wh]

  plot.new()
  pdf(paste0(plot_output_path,"chr.",chrIX,".",site,".hmt.effect.plot.PDF"), width = 8, height = 3)
  par(mar = c(2,4,4,2))
  plot(1,1,type="n", xlab = "Base location", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024)
       ,main = hmt_title, axes=FALSE)

  axis(2)
  axis(1,at = c(1,seq(256,1024,by = 256)),tick = c(1,seq(256,1024,by = 256)))
  if(length(wqtl_hmt_col_posi) > 0){
    for(j in 1:length(wqtl_hmt_col_posi)){
      polygon(c(wqtl_hmt_col_posi[j]-0.5, wqtl_hmt_col_posi[j]-0.5, wqtl_hmt_col_posi[j]+0.5, wqtl_hmt_col_posi[j]+0.5), c(ymin_beta-2, ymax_beta+2, ymax_beta+2, ymin_beta-2), col ="pink", border = NA)
    }
  }

  abline(h = 0, col = "red")
  points(xval, wqtl_hmt_mean, col = "blue", type="l")
  points(xval, beta_l, col = "skyblue", type="l")
  points(xval, beta_r, col = "skyblue", type="l")
  box()

  dev.off()

  return(
    list(wqtl_col_posi = wqtl_col_posi
         ,wqtl_hmt_col_posi = wqtl_hmt_col_posi
         ,wqtl_mean = wqtl_mean
         ,wqtl_sd = wqtl_sd
         ,wqtl_hmt_mean = wqtl_hmt_mean
         ,wqtl_hmt_sd = wqtl_hmt_sd
    )
  )
}
