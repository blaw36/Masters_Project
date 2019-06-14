# Date: 5th June 2019.
# Check for similarity between logLRs at bases and overall, between no-HMT and HMT model, on all 24 SNPs.

# It's just the joint logLR we're after, the remaining values are just the BFs at each scale-loc,
# ie the ratio of the posterior prob of y given gamma = 1 over y given gamma = 0

library(ggplot2)
library(dplyr)
# before HMT
logLR_before = as.matrix(read.table("~/Cpp/WaveQTL/test/dsQTL/output/test1.fph.logLR.txt"))
dim(logLR_before)
# after HMT
logLR_after = as.matrix(read.table("~/Cpp/WaveQTL_HMT/test/dsQTL/output/all_snp_hmt.fph.logLR.txt"))
dim(logLR_after)

logLR_before[1:10,1:10]
logLR_after[1:10,1:10]

lhood_compare <- tibble(snp_name = character()
                        ,before_HMT_logLR = character()
                        ,after_HMT_logLR = character()
)
for(i in 1:24){
  before_dat <- logLR_before[i,]
  after_dat <- logLR_after[i,]
  snp_name <- before_dat[1]
  before_lhood <- as.numeric(before_dat[2])
  after_lhood <- as.numeric(after_dat[2])
  lhood_compare <- add_row(lhood_compare
                           ,snp_name = snp_name
                           ,before_HMT_logLR = before_lhood
                           ,after_HMT_logLR = after_lhood)
}
