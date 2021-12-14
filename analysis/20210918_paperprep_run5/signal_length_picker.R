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

# Both strong
both_data = merge(
  hmt_data[,.(chr,site,hmt_snp=strong_snp,hmt_pval=pval)]
  ,nohmt_data[,.(chr,site,nohmt_snp=strong_snp,nohmt_pval=pval)]
  ,by = c("chr","site")
)

# Both strong, randomly shuffled
set.seed(1)
both_data[sample(1:nrow(both_data),size = nrow(both_data))][hmt_snp==nohmt_snp & (hmt_pval < 0.05) & (nohmt_pval < 0.05)][1:20]

# chr site               hmt_snp   hmt_pval             nohmt_snp nohmt_pval
# 1:  14  872        chr14.76563562 0.04291860        chr14.76563562 0.03793070
# 2:  10  323        chr10.22768032 0.01407000        chr10.22768032 0.01268170
# 3:  17 3026        chr17.76683351 0.00009999        chr17.76683351 0.00009999
# 4:   2  694         chr2.42989319 0.00109989         chr2.42989319 0.00109989
# 5:  13  792             rs7338868 0.00009999             rs7338868 0.00009999
# 6:   6 1901        chr6.156759073 0.00569943        chr6.156759073 0.00569943
# 7:  22 1401        chr22.43056674 0.01630880        chr22.43056674 0.01370100
# 8:  10 1500 chr10.indel.101460848 0.01742360 chr10.indel.101460848 0.01724240
# 9:  12  535        chr12.31165389 0.02221110        chr12.31165389 0.02365620
# 10:  12  698        chr12.47637590 0.03023240        chr12.47637590 0.02156710
# 11:  19  225         chr19.1440303 0.00109989         chr19.1440303 0.00109989
# 12:  19 1595        chr19.17854868 0.00239976        chr19.17854868 0.00239976
# 13:   8 1104        chr8.101342698 0.01243410        chr8.101342698 0.01167580
# 14:  16   54          chr16.373624 0.00539946          chr16.373624 0.00479952
# 15:  17 3027        chr17.76682930 0.02443010        chr17.76682930 0.02459470
# 16:   5 1952        chr5.176363320 0.00009999        chr5.176363320 0.00009999
# 17:  19 1060        chr19.10926409 0.00009999        chr19.10926409 0.00009999
# 18:   3  869         chr3.49820010 0.00269973         chr3.49820010 0.00259974
# 19:  15 1341        chr15.81298565 0.01933800        chr15.81298565 0.01934080
# 20:  10 1161        chr10.80819803 0.00069993        chr10.80819803 0.00069993


# HMT strong and WQtl relatively less strong
# hmt_pval > 0.00009999 to avoid those weird cases
# both_data[hmt_snp==nohmt_snp & (hmt_pval < 0.05) & (hmt_pval > 0.00009999) & (nohmt_pval > 0.05)][order(-((nohmt_pval-hmt_pval)/nohmt_pval))]
both_data[hmt_snp==nohmt_snp & (hmt_pval > 0.00009999)][order(-((nohmt_pval-hmt_pval)/nohmt_pval))][1:20]

# chr site              hmt_snp   hmt_pval            nohmt_snp nohmt_pval
# 1:   1 1703        chr1.32309095 0.23004600        chr1.32309095 0.61849600
# 2:  19 2339       chr19.45665523 0.19901800       chr19.45665523 0.51294100
# 3:  22  595       chr22.29359004 0.01080460       chr22.29359004 0.02476160
# 4:   1 2585        chr1.90060230 0.03123930        chr1.90060230 0.06162830
# 5:  18  516       chr18.46978722 0.19127300       chr18.46978722 0.34516900
# 6:   9  884        chr9.99500358 0.02604260        chr9.99500358 0.04631460
# 7:  16 1211       chr16.30878442 0.21868400       chr16.30878442 0.38209400
# 8:  20 1549       chr20.62142221 0.05345080       chr20.62142221 0.09242750
# 9:  14 1461      chr14.103620224 0.00029997      chr14.103620224 0.00049995
# 10:   6  691  chr6.indel.34870101 0.02772790  chr6.indel.34870101 0.04605380
# 11:  16 1042 chr16.indel.28798733 0.33211800 chr16.indel.28798733 0.54117200
# 12:  16 1995       chr16.74026651 0.42988100       chr16.74026651 0.68423400
# 13:  16 1267       chr16.31791941 0.22999400       chr16.31791941 0.35817500
# 14:  19  843        chr19.7455281 0.12885400        chr19.7455281 0.19970600
# 15:   4  520         chr4.7509761 0.19306700         chr4.7509761 0.29632300
# 16:  10  774 chr10.indel.69315364 0.20453100 chr10.indel.69315364 0.31229600
# 17:   1  328         chr1.2506641 0.19760100         chr1.2506641 0.30114600
# 18:   5 1999       chr5.176812691 0.15957200       chr5.176812691 0.24290000
# 19:   9 1806       chr9.134219921 0.24565100       chr9.134219921 0.37259300
# 20:   9 1742       chr9.133257830 0.04025690       chr9.133257830 0.06085030

# statistically significant ones:
# chr site              hmt_snp   hmt_pval            nohmt_snp nohmt_pval
# 3:  22  595       chr22.29359004 0.01080460       chr22.29359004 0.02476160
# 4:   1 2585        chr1.90060230 0.03123930        chr1.90060230 0.06162830
# 6:   9  884        chr9.99500358 0.02604260        chr9.99500358 0.04631460
# 8:  20 1549       chr20.62142221 0.05345080       chr20.62142221 0.09242750
# 9:  14 1461      chr14.103620224 0.00029997      chr14.103620224 0.00049995
# 10:   6  691  chr6.indel.34870101 0.02772790  chr6.indel.34870101 0.04605380
# 20:   9 1742       chr9.133257830 0.04025690       chr9.133257830 0.06085030

### Based on FDR

# q-value data
hmt_qvalue = qvalue(both_data$hmt_pval)
nohmt_qvalue = qvalue(both_data$nohmt_pval)

hmt_qvalue$pi0
nohmt_qvalue$pi0

both_data[,c("hmt_qval","nohmt_qval") :=
            .(hmt_qvalue$qvalues, nohmt_qvalue$qvalues)]

# Where are qval < 0.05?
# Both (ordered by most similar to least similar qvals)
both_data[hmt_qval < 0.05 & nohmt_qval < 0.05][order(abs(hmt_qval-nohmt_qval))]
# HMT only
both_data[hmt_qval < 0.05 & nohmt_qval >= 0.05][order(-nohmt_qval/hmt_qval)]
# Non-HMT only
both_data[hmt_qval >= 0.05 & nohmt_qval < 0.05][order(nohmt_qval/hmt_qval)]

# Let's look at:
# Both
# 1,2391 (similar qvals)
# 22,133 (less similar qvals)

# HMT only
# 4,1294 (less similar qvals)
# 16,1985 (more similar qvals)
