# From https://github.com/heejungshim/multiseq-ms-figures/blob/master/scripts/ATACseq_qvalues_multiseq_WaveQTL_DESeq2.Rmd

# library(BiocManager)
# BiocManager::install("qvalue")
library("qvalue")
library(data.table)

data_path = "~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/analysis/20210918_paperprep_run5/"

# Read p-value data
hmt_data = setDT(
  readRDS(paste0(data_path,"hmt_pvals.RDS"))
)
nohmt_data = setDT(
  readRDS(paste0(data_path,"nohmt_pvals.RDS"))
)

hmt_data[is.na(pval)]
nohmt_data[is.na(pval)]

hmt_pval = hmt_data$pval
nohmt_pval = nohmt_data$pval

## Compute q-value
del.ix = which(is.na(hmt_pval)==TRUE)
numTests = length(hmt_pval)
length(hmt_pval)
## [1] 8827
totalTests = length(hmt_pval) - length(del.ix)
totalTests
## [1] 8827

# qval.hmt = qvalue(hmt_pval[-del.ix])
# qval.nohmt = qvalue(nohmt_pval[-del.ix])
qval.hmt = qvalue(hmt_pval)
qval.nohmt = qvalue(nohmt_pval)
qval.hmt$pi0
qval.nohmt$pi0
## [1] 0.8354974
## [1] 0.8496002

### WaveQTL_HMT vs WaveQTL
path = paste0(data_path,"qval.png")
png(path, units="in", width = 5, height = 5, res =300)
par(mar = c(4, 4, 1, 1))
yval = -log(qval.hmt$qvalues,10)
xval = -log(qval.nohmt$qvalues,10)
sig.cut = -log(0.05, 10)
nonsig.cut = -log(0.5, 10)
tblue = "#00008B80"
tred = "#ff634780"
tblack = "#CCCCCC80"
tgreen = "#458b0080"
col.list = rep(tblack, length(xval))
wh = which((xval > sig.cut) & (yval > sig.cut))
col.list[wh] = tblue
wh = which((xval < nonsig.cut) & (yval > sig.cut))
col.list[wh] = tred
wh = which((xval > sig.cut) & (yval < nonsig.cut))
col.list[wh] = tgreen
plot(xval, yval, pch=46, cex=4, col= col.list, xlab = "-log10(qvalue) from WaveQTL", ylab ="-log10(qvalue) from WaveQTL HMT")
lines(c(min(xval, yval),max(xval, yval)), c(min(xval, yval),max(xval, yval)), lty = 2, col ="darkgrey")
dev.off()

# Blue is BOTH above -log(0.05,10) ~= 1.30
sum(col.list == tblue)
## 119

# Red is HMT > -log(0.05,10) ~= 1.30, and No-HMT < -log(0.5,10) ~= 0.30
sum(col.list == tred)
## 37

# Green is HMT < -log(0.5,10) ~= 0.30 and No-HMT > -log(0.05,10) ~= 1.30
sum(col.list == tgreen)
## 0

wh  = which(col.list == tred)
qval.hmt$qvalues[wh]
qval.hmt$pvalues[wh]
qval.nohmt$qvalues[wh]
qval.nohmt$pvalues[wh]

wh
# [1]  324  385 1001 1082 1433 1498 1737 2060 2197 2255 2270 2481 2902 2948 3349 3720 3888 3901 3937 4702 4775 4793 5103 5112 5166
# [26] 6114 6571 6719 6804 6815 7013 7234 7659 7660 7678 8019 8275

# Note: should be the same as:
# intersect(which(hmt_data[,pval < 0.05]),which(nohmt_data[,pval > 0.5]))

cbind.data.frame(
  hmt_data[wh]
  ,nohmt_data[wh]
)

### WaveQTL_HMT vs WaveQTL p-values
path = paste0(data_path,"pval.png")
png(path, units="in", width = 5, height = 5, res =300)
par(mar = c(4, 4, 1, 1))
yval = -log(hmt_pval,10)
xval = -log(nohmt_pval,10)
sig.cut = -log(0.05, 10)
nonsig.cut = -log(0.5, 10)
tblue = "#00008B80"
tred = "#ff634780"
tblack = "#CCCCCC80"
tgreen = "#458b0080"
col.list = rep(tblack, length(xval))
wh = which((xval > sig.cut) & (yval > sig.cut))
col.list[wh] = tblue
wh = which((xval < nonsig.cut) & (yval > sig.cut))
col.list[wh] = tred
wh = which((xval > sig.cut) & (yval < nonsig.cut))
col.list[wh] = tgreen
plot(xval, yval, pch=46, cex=4, col= col.list, xlab = "-log10(pvalue) from WaveQTL", ylab ="-log10(pvalue) from WaveQTL HMT")
lines(c(min(xval, yval),max(xval, yval)), c(min(xval, yval),max(xval, yval)), lty = 2, col ="darkgrey")
dev.off()

# Blue is BOTH above -log(0.05,10) ~= 1.30
sum(col.list == tblue)
## 835

# Red is HMT > -log(0.05,10) ~= 1.30, and No-HMT < -log(0.5,10) ~= 0.30
sum(col.list == tred)
## 37

# Green is HMT < -log(0.5,10) ~= 0.30 and No-HMT > -log(0.05,10) ~= 1.30
sum(col.list == tgreen)
## 0

cbind.data.frame(
  hmt_data[col.list == tgreen]
  ,nohmt_data[col.list == tgreen]
)

merged_data <- merge(hmt_data[,.(chr,site
                                 ,hmt_strong_snp = strong_snp
                                 ,hmt_perms = perms
                                 ,hmt_above_thold = above_thold
                                 ,hmt_pval = pval)]
                     ,nohmt_data[,.(chr,site
                                    ,nohmt_strong_snp = strong_snp
                                    ,nohmt_perms = perms
                                    ,nohmt_above_thold = above_thold
                                    ,nohmt_pval = pval)]
                     ,by = c("chr","site"))
merged_data[,table(hmt_strong_snp == nohmt_strong_snp, hmt_pval < nohmt_pval)]

merged_data[,"pval_diff" := hmt_pval - nohmt_pval]
merged_data[,"pval_diff_pct" := pval_diff/nohmt_pval]
merged_data[,hist(pval_diff)]

merged_data[,.(min(hmt_pval),max(hmt_pval),min(pval_diff),median(pval_diff),max(pval_diff),mean(pval_diff))
            ,by = .(dplyr::ntile(hmt_pval,10))][order(dplyr)]
merged_data[,median(pval_diff)]

merged_data[,.(median(pval_diff),mean(pval_diff))
            ,by = .(hmt_strong_snp == nohmt_strong_snp)]

merged_data[hmt_strong_snp == nohmt_strong_snp][order(pval_diff)][1:50]
merged_data[hmt_strong_snp == nohmt_strong_snp][order(-pval_diff)][1:50]

merged_data[hmt_strong_snp != nohmt_strong_snp][hmt_pval < 0.05][order(pval_diff)][1:50]
merged_data[hmt_strong_snp != nohmt_strong_snp][nohmt_pval < 0.05][order(-pval_diff)][1:50]
