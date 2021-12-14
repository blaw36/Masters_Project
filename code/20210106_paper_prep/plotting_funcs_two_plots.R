## Some plotting helper functions for paper

# Data space plot function adapted from 'C:\Users\brend\Dropbox\Uni Stuff - Masters\Research Project\Masters_Project_Git\images\ch1_images.R'
# Inputs: geno and pheno file, snp_no, and where to output the plot
# Effect size plot functions adapted from 'C:\Users\brend\Dropbox\Uni Stuff - Masters\Research Project\Masters_Project_Git\code\effect_size_plot_spartan.R'

## Combined plot
# This function will NOT generate no.QT WCs
# This will take spartan outputs, as inputs (ie. the plotting data RDS)
# This is written as if all the outputs have been moved from Spartan onto local PC for ad-hoc analysis.

# This differs from plotting_funcs.R in terms of how it handles mismatched strong SNPs, but doesn't allow for specific SNP input.
# This function will create two copies of plots for each site/chr pair; one using the strongest SNP from WaveQTL analysis, the other using
# strongest SNP from HMT analysis.

# It also takes in qvalue data from Spartan. Assumed this is in an RDS file in the root directory.
combined_plot <- function(
  geno_path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/"
  ,pheno_path = "/home/bklaw/paper_data_clean/analysis/simulation/sample_size/simulation_manydsQTL_v1/data/"
  ,site
  ,chrIX
  ,wqtl_qt_output_path = "output/"
  ,wqtl_no_qt_output_path = "output/"
  ,plot_output_path = "plots/"
  ,qval_file_name = "all_qval_data.RDS"
  # ,snp_to_use = 'nearest' # other options; 'hmt', 'nohmt'
){

  # source("C:/Users/brend/Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/code/sim1_revised_hmt_nohmt_functions.R")
  source("/home/bklaw/WaveQTL_HMT_wperm/R/sim1_revised_hmt_nohmt_functions.R")

  ## Read DWT matrix
  # Wmat_1024 = read.table("C:/Users/brend/Documents/Cpp/WaveQTL/data/DWT/Wmat_1024",as.is = TRUE)
  Wmat_1024 = read.table("/home/bklaw/WaveQTL/WaveQTL-master/data/DWT/Wmat_1024",as.is = TRUE)

  W2mat_1024 = Wmat_1024*Wmat_1024

  # Strongest SNP - WaveQTL
  nearby_snp_nohmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".nohmt.fph.pval.txt"), stringsAsFactors = F)$V1[1]
  test = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".nohmt.fph.logLR.txt"), stringsAsFactors = F) # Read in one of the mean or var files we'll be using here
  sel_geno_IX = which(nearby_snp_nohmt == test[,1]) # test[,1] should relate to the row names, ie the SNP names; we want to grab the row of the strongest SNP.

  # Strongest SNP - HMT
  nearby_snp_hmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".hmt.fph.pval.txt"), stringsAsFactors = F)$V1[1]
  testhmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".hmt.fph.logLR.txt"), stringsAsFactors = F) # Read in one of the mean or var files we'll be using here
  sel_geno_IX_hmt = which(nearby_snp_hmt == testhmt[,1]) # test[,1] should relate to the row names, ie the SNP names; we want to grab the row of the strongest SNP.

  if(sel_geno_IX == sel_geno_IX_hmt){
    same_snp_mode = TRUE
  }else{
    same_snp_mode = FALSE
  }

  # WaveQTL strongest SNP plot. Mention this at top of plot (WaveQTL strongest SNP: 'SNP name')
  # WAveQTL pval, qval only

  # If same, include all details

  # WaveQTL strongest SNP plot ----------------------------------------------


  # WaveQTL samples
  logLR_nohmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".nohmt.fph.logLR.txt"), stringsAsFactors = F)$V2[sel_geno_IX]
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
  logLR_hmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".hmt.fph.logLR.txt"), stringsAsFactors = F)$V2[sel_geno_IX]
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

  library(data.table)

  # Plot titles

  pval_nohmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".nohmt.fph.pval.txt"), stringsAsFactors = F)$V1[4]
  qval_data = setDT(readRDS(qval_file_name))
  site_idx = site
  qval_nohmt = qval_data[chr == chrIX & site == site_idx, nohmt_qval]

  if(same_snp_mode){
    # SNPs are same
    pval_hmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".hmt.fph.pval.txt"), stringsAsFactors = F)$V1[4]
    qval_hmt = qval_data[chr == chrIX & site == site_idx, hmt_qval]

    nohmt_title = paste0("WaveQTL, pval: ",round(as.numeric(pval_nohmt),4), ", qval: ", round(as.numeric(qval_nohmt),4), ", logLR: ", round(as.numeric(logLR_nohmt),4))
    hmt_title = paste0("HMT, pval: ",round(as.numeric(pval_hmt),4), ", qval: ", round(as.numeric(qval_hmt),4), ", logLR: ", round(as.numeric(logLR_hmt),4))
  }else{
    # Strongest SNPs not the same, can't retrieve p-val or q-val data for HMT and SNP isn't strongest for HMT
    nohmt_title = paste0("WaveQTL, pval: ",round(as.numeric(pval_nohmt),4), ", qval: ", round(as.numeric(qval_nohmt),4), ", logLR: ", round(as.numeric(logLR_nohmt),4))
    hmt_title = paste0("HMT, logLR: ", round(as.numeric(logLR_hmt),4))
  }

  ## Data space plot data
  library(reshape2)

  # Load genotype (SNP) data
  # eg_geno <- as.matrix(read.table(paste0(data.path,"chr17.10160989.10162012.2kb.cis.geno")))
  eg_geno <- as.matrix(read.table(paste0(geno_path,"chr",chrIX,".",site,".geno")))

  eg_geno.2 = data.frame(t(eg_geno[,c(4:73)]))
  names(eg_geno.2) <- eg_geno[,1]

  lapply(eg_geno.2,table)

  # Read in phenotype data
  # pheno.dat = as.matrix(read.table(paste0(data.path,"chr17.10160989.10162012.pheno.dat")))
  pheno.dat = as.matrix(read.table(paste0(pheno_path,"chr",chrIX,".pheno.dat.",site)))
  pheno.dat.2=data.frame(t(pheno.dat))
  pheno.dat.2 = melt(pheno.dat.2)
  pheno.dat.2$loc = rep(1:1024,70)


  # Bucketed up based on SNP value ----
  # Read in genotype, round to nearest bucket (0,1,2)
  specific_snp <- as.numeric(eg_geno[sel_geno_IX,][4:73])
  specific_snp_bucketed <- round(specific_snp)
  # Count num of unique buckets
  unique_buckets = unique(specific_snp_bucketed)
  num_unique_buckets = length(unique_buckets)
  mean_pheno_list = list()
  for(i in 1:num_unique_buckets){
    bucket_val = unique_buckets[i]
    bucket_pheno_data = pheno.dat[specific_snp_bucketed == bucket_val, ]
    if(is.null(dim(bucket_pheno_data))){
      mean_pheno_list[[i]] = data.frame(mean = bucket_pheno_data, loc = 1:1024)
    }else{
      mean_pheno_list[[i]] = data.frame(mean = apply(bucket_pheno_data,2,mean), loc = 1:1024)
    }
  }
  # mean.pheno_0 = data.frame(mean = apply(pheno.dat[specific_snp_bucketed == 0, ],2,mean), loc = 1:1024)
  # mean.pheno_1 = data.frame(mean = apply(pheno.dat[specific_snp_bucketed == 1, ],2,mean), loc = 1:1024)
  # mean.pheno_2 = data.frame(mean = apply(pheno.dat[specific_snp_bucketed == 2, ],2,mean), loc = 1:1024)
  # y_bounds <- c(0,max(max(mean.pheno_0$mean,mean.pheno_1$mean),mean.pheno_2$mean))*1.05
  y_bounds <- c(0,max(unlist(lapply(mean_pheno_list,function(x){max(x$mean)}))))*1.05

  data_space_title = paste0("WaveQTL strongest SNP - Chr: ",chrIX,", Site: ",site,", SNP: ",nearby_snp_nohmt)

  # Plot effect size plots together

  # Data space plot
  plot.new()
  pdf(paste0(plot_output_path,"chr.",chrIX,".",site,".wqtl.snp.all.effect.plot.PDF"), width = 8, height = 8)
  par(mar = c(2,4,4,2))
  par(mfrow=c(3,1))

  plot(1,1
       ,type = "l"
       ,xlab = "Base location"
       ,ylab = "Normalised count"
       ,main = data_space_title
       ,xaxt = "n"
       ,xlim = c(0,dim(pheno.dat)[2])
       ,ylim = y_bounds)
  axis(1, at = seq(0,1024,128),labels = seq(0,1024,128))
  for(i in rev(1:num_unique_buckets)){
    lines(mean_pheno_list[[i]]$mean, col = i)
  }
  # lines(mean.pheno_2$mean, col = "red")
  # lines(mean.pheno_1$mean, col = "green")
  # lines(mean.pheno_0$mean, col = "black")
  # legend(x = 750, y = 0.7, title = "Genotype value",legend = unique_buckets
  #        ,cex=0.8,col=rev(1:num_unique_buckets),lty = rep(1,num_unique_buckets)
  #        ,box.lty=0)
  legend(x = "topright", title = "Genotype value",legend = unique_buckets, bg = "transparent"
         ,cex=0.8,col=rev(1:num_unique_buckets),lty = rep(1,num_unique_buckets)
         ,box.lty=0)

  # WQtl plot pre-amble
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

  # WQtl plot
  beta_l = wqtl_mean - 3*wqtl_sd
  beta_r = wqtl_mean + 3*wqtl_sd

  wh_l = which(beta_l > 0)
  wh_r = which(beta_r < 0)
  high_wh = sort(unique(union(wh_l, wh_r)))

  xval = 1:1024
  wqtl_col_posi = xval[high_wh]

  plot(1,2,type="n", xlab = "Base location", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024)
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

  # HMT Plot

  beta_l_hmt = wqtl_hmt_mean - 3*wqtl_hmt_sd
  beta_r_hmt = wqtl_hmt_mean + 3*wqtl_hmt_sd

  wh_l_hmt = which(beta_l_hmt > 0)
  wh_r_hmt = which(beta_r_hmt < 0)
  high_wh_hmt = sort(unique(union(wh_l_hmt, wh_r_hmt)))

  xval = 1:1024
  wqtl_hmt_col_posi = xval[high_wh_hmt]

  plot(3,1,type="n", xlab = "Base location", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024)
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
  points(xval, beta_l_hmt, col = "skyblue", type="l")
  points(xval, beta_r_hmt, col = "skyblue", type="l")
  box()

  dev.off()

  # If SNP is same, just repeat the plot
  if(same_snp_mode){

    data_space_title = paste0("HMT strongest SNP - Chr: ",chrIX,", Site: ",site,", SNP: ",nearby_snp_nohmt)

    # Plot effect size plots together

    # Data space plot
    plot.new()
    pdf(paste0(plot_output_path,"chr.",chrIX,".",site,".hmt.snp.all.effect.plot.PDF"), width = 8, height = 8)
    par(mar = c(2,4,4,2))
    par(mfrow=c(3,1))

    plot(1,1
         ,type = "l"
         ,xlab = "Base location"
         ,ylab = "Normalised count"
         ,main = data_space_title
         ,xaxt = "n"
         ,xlim = c(0,dim(pheno.dat)[2])
         ,ylim = y_bounds)
    axis(1, at = seq(0,1024,128),labels = seq(0,1024,128))
    for(i in rev(1:num_unique_buckets)){
      lines(mean_pheno_list[[i]]$mean, col = i)
    }
    # lines(mean.pheno_2$mean, col = "red")
    # lines(mean.pheno_1$mean, col = "green")
    # lines(mean.pheno_0$mean, col = "black")
    # legend(x = 750, y = 0.7, title = "Genotype value",legend = unique_buckets
    #        ,cex=0.8,col=rev(1:num_unique_buckets),lty = rep(1,num_unique_buckets)
    #        ,box.lty=0)
    legend(x = "topright", title = "Genotype value",legend = unique_buckets, bg = "transparent"
           ,cex=0.8,col=rev(1:num_unique_buckets),lty = rep(1,num_unique_buckets)
           ,box.lty=0)

    # WQtl plot
    plot(1,2,type="n", xlab = "Base location", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024)
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

    # HMT Plot
    plot(3,1,type="n", xlab = "Base location", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024)
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
    points(xval, beta_l_hmt, col = "skyblue", type="l")
    points(xval, beta_r_hmt, col = "skyblue", type="l")
    box()

    dev.off()

  }else{

    # HMT strongest SNP plot --------------------------------------------------
    # Only in case where HMT strongest SNP not equal to WQtl strongest SNP

    # WaveQTL samples
    logLR_nohmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".nohmt.fph.logLR.txt"), stringsAsFactors = F)$V2[sel_geno_IX_hmt]
    wqtl_samples <- draw_samples_no_hmt(
      data_path = wqtl_no_qt_output_path
      ,dataset = paste0("chr.",chrIX,".",site,".nohmt.no.QT")
      ,geno_select = sel_geno_IX_hmt
      ,wmat = Wmat_1024
      ,seed_no = 10
      ,num_samples = 1000
      ,num_pheno = 1024
      ,num_indiv = 70
    )

    # WaveQTL_HMT samples
    logLR_hmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".hmt.fph.logLR.txt"), stringsAsFactors = F)$V2[sel_geno_IX_hmt]
    wqtl_hmt_samples <- draw_samples_hmt(
      data_path = wqtl_no_qt_output_path
      ,dataset = paste0("chr.",chrIX,".",site,".hmt.no.QT")
      ,waveqtl_data_path = wqtl_no_qt_output_path
      ,waveqtl_dataset = paste0("chr.",chrIX,".",site,".nohmt.no.QT")
      ,geno_select = sel_geno_IX_hmt
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

    beta_l = wqtl_mean - 3*wqtl_sd
    beta_r = wqtl_mean + 3*wqtl_sd

    wh_l = which(beta_l > 0)
    wh_r = which(beta_r < 0)
    high_wh = sort(unique(union(wh_l, wh_r)))

    xval = 1:1024
    wqtl_col_posi = xval[high_wh]

    # library(data.table)

    # Plot titles
    pval_hmt = read.table(paste0(wqtl_qt_output_path,"chr.",chrIX,".",site,".hmt.fph.pval.txt"), stringsAsFactors = F)$V1[4]
    qval_hmt = qval_data[chr == chrIX & site == site_idx, hmt_qval]
    nohmt_title = paste0("WaveQTL, logLR: ", round(as.numeric(logLR_nohmt),4))
    hmt_title = paste0("HMT, pval: ",round(as.numeric(pval_hmt),4), ", qval: ", round(as.numeric(qval_hmt),4), ", logLR: ", round(as.numeric(logLR_hmt),4))

    ## Data space plot data
    # library(reshape2)

    # Load genotype (SNP) data
    # eg_geno <- as.matrix(read.table(paste0(data.path,"chr17.10160989.10162012.2kb.cis.geno")))
    eg_geno <- as.matrix(read.table(paste0(geno_path,"chr",chrIX,".",site,".geno")))

    eg_geno.2 = data.frame(t(eg_geno[,c(4:73)]))
    names(eg_geno.2) <- eg_geno[,1]

    lapply(eg_geno.2,table)

    # Read in phenotype data
    # pheno.dat = as.matrix(read.table(paste0(data.path,"chr17.10160989.10162012.pheno.dat")))
    pheno.dat = as.matrix(read.table(paste0(pheno_path,"chr",chrIX,".pheno.dat.",site)))
    pheno.dat.2=data.frame(t(pheno.dat))
    pheno.dat.2 = melt(pheno.dat.2)
    pheno.dat.2$loc = rep(1:1024,70)


    # Bucketed up based on SNP value ----
    # Read in genotype, round to nearest bucket (0,1,2)
    specific_snp <- as.numeric(eg_geno[sel_geno_IX_hmt,][4:73])
    specific_snp_bucketed <- round(specific_snp)
    # Count num of unique buckets
    unique_buckets = unique(specific_snp_bucketed)
    num_unique_buckets = length(unique_buckets)
    mean_pheno_list = list()
    for(i in 1:num_unique_buckets){
      bucket_val = unique_buckets[i]
      bucket_pheno_data = pheno.dat[specific_snp_bucketed == bucket_val, ]
      if(is.null(dim(bucket_pheno_data))){
        mean_pheno_list[[i]] = data.frame(mean = bucket_pheno_data, loc = 1:1024)
      }else{
        mean_pheno_list[[i]] = data.frame(mean = apply(bucket_pheno_data,2,mean), loc = 1:1024)
      }
    }
    # mean.pheno_0 = data.frame(mean = apply(pheno.dat[specific_snp_bucketed == 0, ],2,mean), loc = 1:1024)
    # mean.pheno_1 = data.frame(mean = apply(pheno.dat[specific_snp_bucketed == 1, ],2,mean), loc = 1:1024)
    # mean.pheno_2 = data.frame(mean = apply(pheno.dat[specific_snp_bucketed == 2, ],2,mean), loc = 1:1024)
    # y_bounds <- c(0,max(max(mean.pheno_0$mean,mean.pheno_1$mean),mean.pheno_2$mean))*1.05
    y_bounds <- c(0,max(unlist(lapply(mean_pheno_list,function(x){max(x$mean)}))))*1.05

    data_space_title = paste0("HMT strongest SNP - Chr: ",chrIX,", Site: ",site,", SNP: ",nearby_snp_hmt)

    # Plot effect size plots together

    plot.new()
    pdf(paste0(plot_output_path,"chr.",chrIX,".",site,".hmt.snp.all.effect.plot.PDF"), width = 8, height = 8)
    par(mar = c(2,4,4,2))
    par(mfrow=c(3,1))

    plot(1,1
         ,type = "l"
         ,xlab = "Base location"
         ,ylab = "Normalised count"
         ,main = data_space_title
         ,xaxt = "n"
         ,xlim = c(0,dim(pheno.dat)[2])
         ,ylim = y_bounds)
    axis(1, at = seq(0,1024,128),labels = seq(0,1024,128))
    for(i in rev(1:num_unique_buckets)){
      lines(mean_pheno_list[[i]]$mean, col = i)
    }
    # lines(mean.pheno_2$mean, col = "red")
    # lines(mean.pheno_1$mean, col = "green")
    # lines(mean.pheno_0$mean, col = "black")
    # legend(x = 750, y = 0.7, title = "Genotype value",legend = unique_buckets
    #        ,cex=0.8,col=rev(1:num_unique_buckets),lty = rep(1,num_unique_buckets)
    #        ,box.lty=0)
    legend(x = "topright", title = "Genotype value",legend = unique_buckets, bg = "transparent"
           ,cex=0.8,col=rev(1:num_unique_buckets),lty = rep(1,num_unique_buckets)
           ,box.lty=0)

    plot(1,2,type="n", xlab = "Base location", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024)
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

    beta_l = wqtl_hmt_mean - 3*wqtl_hmt_sd
    beta_r = wqtl_hmt_mean + 3*wqtl_hmt_sd

    wh_l = which(beta_l > 0)
    wh_r = which(beta_r < 0)
    high_wh = sort(unique(union(wh_l, wh_r)))

    xval = 1:1024
    wqtl_hmt_col_posi = xval[high_wh]

    plot(3,1,type="n", xlab = "Base location", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024)
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
  }

}

# Local version
# # Same SNPS
# combined_plot(
#   geno_path = "C:/Users/brend/Documents/Cpp/20210817_Wqtl/data/geno_data/"
#   ,pheno_path = "C:/Users/brend/Documents/Cpp/20210817_Wqtl/data/pheno_data/"
#   ,site = 845
#   ,chrIX = 4
#   ,wqtl_qt_output_path = "C:/Users/brend/Documents/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/"
#   ,wqtl_no_qt_output_path = "C:/Users/brend/Documents/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/"
#   ,plot_output_path = "C:/Users/brend/Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/analysis/20210918_paperprep_run5/effect_size/"
#   ,qval_file_name = "C:/Users/brend/Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/analysis/20210918_paperprep_run5/all_qval_data.RDS"
# )
# # Different SNPs
# combined_plot(
#   geno_path = "C:/Users/brend/Documents/Cpp/20210817_Wqtl/data/geno_data/"
#   ,pheno_path = "C:/Users/brend/Documents/Cpp/20210817_Wqtl/data/pheno_data/"
#   ,site = 825
#   ,chrIX = 1
#   ,wqtl_qt_output_path = "C:/Users/brend/Documents/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/"
#   ,wqtl_no_qt_output_path = "C:/Users/brend/Documents/Cpp/20210817_Wqtl/test/dsQTL/output_from_spartan/"
#   ,plot_output_path = "C:/Users/brend/Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/analysis/20210705_paperprep_run4/effect_size/"
#   ,qval_file_name = "C:/Users/brend/Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/analysis/20210918_paperprep_run5/all_qval_data.RDS"
# )

# Spartan version
combined_plot(site = 23, chrIX = 1)
combined_plot(site = 111, chrIX = 1)
combined_plot(site = 453, chrIX = 1)
combined_plot(site = 1396, chrIX = 1)
combined_plot(site = 1559, chrIX = 1)
combined_plot(site = 1914, chrIX = 1)
combined_plot(site = 2109, chrIX = 1)
combined_plot(site = 2391, chrIX = 1)
combined_plot(site = 2562, chrIX = 1)
combined_plot(site = 3213, chrIX = 1)
combined_plot(site = 3352, chrIX = 1)
combined_plot(site = 3506, chrIX = 1)
combined_plot(site = 4063, chrIX = 1)
combined_plot(site = 4280, chrIX = 1)
combined_plot(site = 4705, chrIX = 1)
combined_plot(site = 268, chrIX = 2)
combined_plot(site = 694, chrIX = 2)
combined_plot(site = 808, chrIX = 2)
combined_plot(site = 1066, chrIX = 2)
combined_plot(site = 1333, chrIX = 2)
combined_plot(site = 1617, chrIX = 2)
combined_plot(site = 1736, chrIX = 2)
combined_plot(site = 1791, chrIX = 2)
combined_plot(site = 1918, chrIX = 2)
combined_plot(site = 2105, chrIX = 2)
combined_plot(site = 2376, chrIX = 2)
combined_plot(site = 2623, chrIX = 2)
combined_plot(site = 2796, chrIX = 2)
combined_plot(site = 2941, chrIX = 2)
combined_plot(site = 3082, chrIX = 2)
combined_plot(site = 3327, chrIX = 2)
combined_plot(site = 113, chrIX = 3)
combined_plot(site = 1159, chrIX = 3)
combined_plot(site = 1327, chrIX = 3)
combined_plot(site = 1579, chrIX = 3)
combined_plot(site = 2242, chrIX = 3)
combined_plot(site = 2398, chrIX = 3)
combined_plot(site = 10, chrIX = 4)
combined_plot(site = 587, chrIX = 4)
combined_plot(site = 792, chrIX = 4)
combined_plot(site = 1294, chrIX = 4)
combined_plot(site = 138, chrIX = 5)
combined_plot(site = 206, chrIX = 5)
combined_plot(site = 247, chrIX = 5)
combined_plot(site = 540, chrIX = 5)
combined_plot(site = 686, chrIX = 5)
combined_plot(site = 689, chrIX = 5)
combined_plot(site = 900, chrIX = 5)
combined_plot(site = 994, chrIX = 5)
combined_plot(site = 1117, chrIX = 5)
combined_plot(site = 1397, chrIX = 5)
combined_plot(site = 1952, chrIX = 5)
combined_plot(site = 2023, chrIX = 5)
combined_plot(site = 175, chrIX = 6)
combined_plot(site = 458, chrIX = 6)
combined_plot(site = 662, chrIX = 6)
combined_plot(site = 897, chrIX = 6)
combined_plot(site = 1238, chrIX = 6)
combined_plot(site = 89, chrIX = 7)
combined_plot(site = 554, chrIX = 7)
combined_plot(site = 1684, chrIX = 7)
combined_plot(site = 1765, chrIX = 7)
combined_plot(site = 2217, chrIX = 7)
combined_plot(site = 37, chrIX = 8)
combined_plot(site = 145, chrIX = 8)
combined_plot(site = 179, chrIX = 8)
combined_plot(site = 302, chrIX = 8)
combined_plot(site = 520, chrIX = 8)
combined_plot(site = 935, chrIX = 8)
combined_plot(site = 965, chrIX = 8)
combined_plot(site = 1204, chrIX = 8)
combined_plot(site = 1565, chrIX = 8)
combined_plot(site = 1749, chrIX = 8)
combined_plot(site = 350, chrIX = 9)
combined_plot(site = 478, chrIX = 9)
combined_plot(site = 951, chrIX = 9)
combined_plot(site = 1131, chrIX = 9)
combined_plot(site = 1477, chrIX = 9)
combined_plot(site = 1568, chrIX = 9)
combined_plot(site = 1580, chrIX = 9)
combined_plot(site = 1648, chrIX = 9)
combined_plot(site = 1658, chrIX = 9)
combined_plot(site = 2453, chrIX = 9)
combined_plot(site = 18, chrIX = 10)
combined_plot(site = 239, chrIX = 10)
combined_plot(site = 392, chrIX = 10)
combined_plot(site = 456, chrIX = 10)
combined_plot(site = 1089, chrIX = 10)
combined_plot(site = 1161, chrIX = 10)
combined_plot(site = 1168, chrIX = 10)
combined_plot(site = 1221, chrIX = 10)
combined_plot(site = 1452, chrIX = 10)
combined_plot(site = 2081, chrIX = 10)
combined_plot(site = 231, chrIX = 11)
combined_plot(site = 908, chrIX = 11)
combined_plot(site = 1088, chrIX = 11)
combined_plot(site = 1448, chrIX = 11)
combined_plot(site = 1717, chrIX = 11)
combined_plot(site = 1813, chrIX = 11)
combined_plot(site = 2099, chrIX = 11)
combined_plot(site = 2391, chrIX = 11)
combined_plot(site = 42, chrIX = 12)
combined_plot(site = 67, chrIX = 12)
combined_plot(site = 252, chrIX = 12)
combined_plot(site = 326, chrIX = 12)
combined_plot(site = 327, chrIX = 12)
combined_plot(site = 1354, chrIX = 12)
combined_plot(site = 1504, chrIX = 12)
combined_plot(site = 2284, chrIX = 12)
combined_plot(site = 2368, chrIX = 12)
combined_plot(site = 210, chrIX = 13)
combined_plot(site = 452, chrIX = 13)
combined_plot(site = 498, chrIX = 13)
combined_plot(site = 746, chrIX = 13)
combined_plot(site = 792, chrIX = 13)
combined_plot(site = 797, chrIX = 13)
combined_plot(site = 10, chrIX = 14)
combined_plot(site = 488, chrIX = 14)
combined_plot(site = 694, chrIX = 14)
combined_plot(site = 810, chrIX = 14)
combined_plot(site = 1109, chrIX = 14)
combined_plot(site = 1180, chrIX = 14)
combined_plot(site = 1461, chrIX = 14)
combined_plot(site = 1584, chrIX = 14)
combined_plot(site = 213, chrIX = 15)
combined_plot(site = 856, chrIX = 15)
combined_plot(site = 1163, chrIX = 15)
combined_plot(site = 236, chrIX = 16)
combined_plot(site = 239, chrIX = 16)
combined_plot(site = 1985, chrIX = 16)
combined_plot(site = 2107, chrIX = 16)
combined_plot(site = 566, chrIX = 17)
combined_plot(site = 902, chrIX = 17)
combined_plot(site = 1221, chrIX = 17)
combined_plot(site = 1727, chrIX = 17)
combined_plot(site = 1972, chrIX = 17)
combined_plot(site = 2101, chrIX = 17)
combined_plot(site = 2131, chrIX = 17)
combined_plot(site = 2553, chrIX = 17)
combined_plot(site = 2757, chrIX = 17)
combined_plot(site = 2850, chrIX = 17)
combined_plot(site = 3026, chrIX = 17)
combined_plot(site = 3129, chrIX = 17)
combined_plot(site = 43, chrIX = 18)
combined_plot(site = 319, chrIX = 18)
combined_plot(site = 477, chrIX = 18)
combined_plot(site = 786, chrIX = 18)
combined_plot(site = 121, chrIX = 19)
combined_plot(site = 225, chrIX = 19)
combined_plot(site = 437, chrIX = 19)
combined_plot(site = 750, chrIX = 19)
combined_plot(site = 1060, chrIX = 19)
combined_plot(site = 2379, chrIX = 19)
combined_plot(site = 3114, chrIX = 19)
combined_plot(site = 3125, chrIX = 19)
combined_plot(site = 3322, chrIX = 19)
combined_plot(site = 375, chrIX = 20)
combined_plot(site = 537, chrIX = 20)
combined_plot(site = 932, chrIX = 20)
combined_plot(site = 1087, chrIX = 20)
combined_plot(site = 80, chrIX = 21)
combined_plot(site = 249, chrIX = 21)
combined_plot(site = 131, chrIX = 22)
combined_plot(site = 133, chrIX = 22)
combined_plot(site = 209, chrIX = 22)
combined_plot(site = 478, chrIX = 22)
combined_plot(site = 527, chrIX = 22)
combined_plot(site = 791, chrIX = 22)
combined_plot(site = 1386, chrIX = 22)
combined_plot(site = 1453, chrIX = 22)
combined_plot(site = 1530, chrIX = 22)
