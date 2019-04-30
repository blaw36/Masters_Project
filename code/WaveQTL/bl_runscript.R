# WaveQTL_manual
# Author: Brendan Law
# Date: 21 March 2019

# First attempt at running through commands in WaveQTL_manual.PDF from 
# HJ's WaveQTL package.
# https://github.com/heejungshim/WaveQTL
# /doc/manual/WaveQTL_manual.pdf

# install.packages("wavethresh")

library(ggplot2)
library(reshape2)

# Section 3: Functional phenotypic... -------------------------------------

setwd("~/Cpp/WaveQTL/R/")
source("WaveQTL_preprocess_funcs.R")

# wc <- wd(c(USAccDeaths)[1:64],filter.number = 1,family = "DaubExPhase")

## Set directory to:
  # Folder with data in Figure 2 of Shim and Stephens (2014) paper
  # This is the "dsQTL data...containing DNase-seq data at chr17.10160989.10162012
  # and genotypes at 24 SNPs in 2kb cis-candidate region on 70 individuals."
data.path = "~/Cpp/WaveQTL/data/dsQTL/"

## The following code has been copied from the manual:

# read functional phenotypic data
pheno.dat = as.matrix(read.table(paste0(data.path,"chr17.10160989.10162012.pheno.dat")))
dim(pheno.dat)
# BL: 70 individuals, on a region of 1024 base pairs (bps).
# BL: 'In our dsQTL analysis of [4], we obtained functional phe-
# notypic data after applying multiple procedures into read counts to avoid potential biases.'

# Lots of zeroes
sum(pheno.dat == 0)/(nrow(pheno.dat)*ncol(pheno.dat))
# [1] 0.9563616
# Mainly 0s, means around 0.01 - 0.03, max is 3.5 for all individuals across all bases
apply(pheno.dat,1,summary)

# Can plot them...
pheno.dat.2=data.frame(t(pheno.dat))
pheno.dat.2 = melt(pheno.dat.2)
pheno.dat.2$loc = rep(1:1024,70)
ggplot(pheno.dat.2) + geom_path(aes(x = loc, y = value, group = factor(variable), alpha = 0.5))
ggplot(melt(t(apply(pheno.dat,2,summary)))) + 
    geom_path(aes(x = as.numeric(Var1), y = value, group = factor(Var2), colour = factor(Var2))) + 
    scale_x_continuous(breaks = seq(0,1024,by = 64)) + xlab("Base location") + 
    guides(colour = guide_legend("Stat"))

# The means are the most interesting
mean.pheno = data.frame(mean = apply(pheno.dat,2,mean), loc = 1:1024)
ggplot(mean.pheno) + 
	geom_line(aes(x = loc, y = mean)) + 
	scale_x_continuous(breaks = seq(0,1024,by = 64)) + 
	xlab("Base location")    

# read library read depth
library.read.depth = scan(paste0(data.path, "library.read.depth.dat"))
length(library.read.depth)
# BL: 'read depth' (??) for each of the 70 individuals

# read Covariates
Covariates = as.matrix(read.table(paste0(data.path, "PC4.dat")))
dim(Covariates)
# BL: are these like 'g' for the 70 individuals? 4 regressor factors?
# BL: Why does the WC transform depend on these covariates (if used) for normalisation?
# BL: Covariates contains multiple covariates to be
# corrected for (here, four principal components we used to control for confounding factors in our
# dsQTL analysis of [4]).
# BL: Oh. This is the part where we believe there are 'other' covariates outside of the y,g we're studying
# and so, as part of pre-processing (but after doing wavelet t/forms), 
# we find the (4 in this case) principal components which capture
# the most variation. We PCA regress these covariates out (to attempt to remove the effects of anything else
# from our analysis), and then we take the residuals afterwards, and start from there.
# Ie we do our WC transforms and wavelet analysis on 'what's left' after significant other covariates (latent,
# in this case, picked up by PCA) have been regressed away, keeping just the residuals. See 'corrected_forCovariates'
# and 'Normalize.WCs' functions for more details.

# Data -> WCs -> Quantile ranking -> Regress out covariates through PCA regression (if provided) ->
# Quantile rank again for standard normality of 'regressed out' -> WCs to begin analysis with.
# Only done for inference, NOT for estimation, otherwise effect size in data space would be difficult
# to obtain and invert with this series of transformations.

## Second part - waveQTL_preprocessing
# 'meanR.thresh' = 2 to filter out low count WCs.
meanR.thresh = 2

# Some under the hood things to know:
  # Uses family '1' of the 'wd' function in 'wd' package as well as
  # filter 'DaubExPhase', both of which combine to make a DWT as desired.
  # 1) Filters out if mean count below the thresh. (0 if filtered out, 1 if retained)
    # This does a mean of the sum of counts of each individual, across various scales and locations.
    # The scales and locations are in line with the 'span' of the WCs
    # ie sum across all locs for mean counts 1 and 2
    # sum across first half for 3, second half for 4, ..., (B-1) to B at scale J.
  # 2) Correct counts for read depth (divide by read depth)
  # 3) Do the wavelet t/form
  # 4) Normalise:
    # a) quantile of rank of WCs
    # b) take residuals after regressing covariates out of WCs
    # c) quantile of rank of residuals
res = WaveQTL_preprocess(Data = pheno.dat
                         , library.read.depth = library.read.depth
                         , Covariates = Covariates
                         , meanR.thresh = meanR.thresh)
### A more in-depth look into the FWT function and the underlying
### 'wd' function for the wavelet transform can be found in my 
### 'Wavelets_in_R.Rmd' file.

### In summary, FWT spits out:
# WCs = an N * T matrix, where, as desired:
# 1st element = sum of all raw data
# 2nd element = difference of the two halves
# etc ...
# last element = difference between Jth and (J-1)th element
# 'difference between halves' = right - left
# All differences and sums are scaled at each level. Once at
# the bottom level, and then each level on top of that uses
# the scaled sums and differences of that below. Hence the top level
# is scaled 'nlevels' times on the raw data.

str(res)

# If you step through the step 4, you can see the progression
# and how the quantile rank makes the tails heavier and dist'n less peaky at 0
# then how the residuals makes the values more normal, and then quantile ranking
# doing the same effect.

# res$WCs:
  # Ordered from low res to high res WC (scaling coefficient)
  # {(0,0),(1,1),(2,1),(2,2),...,(J,1),...,(J,2^(J-1))}

# Save as text files
output.path = "../test/dsQTL/"
write.table(res$WCs, file= paste0(output.path, "WCs.txt"), row.names=FALSE,
            col.names = FALSE, quote=FALSE)
cat(res$filtered.WCs, file = paste0(output.path, "use.txt"))



# Section 4: WaveQTL ------------------------------------------------------

## Genotype data (covariates)
# Here's an example:
# 24 x 73 matrix which shows 24 SNPs for 70 individuals:
# Cols 1 - 3: SNP id, allele types with minor allele first
# Cols 4 - 70: (posterior) mean genotypes (0 - 2) of each individual, 1,...,70
eg_geno <- as.matrix(read.table(paste0(data.path,"chr17.10160989.10162012.2kb.cis.geno")))

# Mean values for each 24 SNPs
plot(apply(eg_geno[,4:73],1,function(x){mean(as.numeric(x))}),type = 'l')

# Can take non-genotype covariate (control/treatment group indicator) for association analysis:
# Arbitrary SNP id and allele coding using ACGT in first three columns.
# User should use "-notsnp" option to disable the minor allele frequency cutoff and
# to use any numerical values as covariates.
# Is it something like this? Eg: first 35 in grp 1, second 35 not in grp 1

eg_non_geno <- matrix(data = c("mbr_grp_1","A","A"
                               ,rep(1,35)
                               ,rep(0,35))
                      ,nrow = 1)
dim(eg_non_geno)

## Phenotype data (response)
# WCs, with 70 rows (in the same order as the 70 columns of genotype data)
# after passing through WaveQTL_preprocess above. 
# Columns represent low - high res/coarse - fine grain/high - low scale WCs.

# Note the supplementary datasets:
  # filtered.WCs: which of the 1024 columns to filter out due to low counts: filter with "-u" option
  # (if provided), which scale-locations the hyperparam 'pi' corresponds to: "-group" option

# The second file comes in the form:
# Vector of multiple positive integers, indicating the start position of each group of WCs in the phenotype file.
# Eg: 1,2,9,17
# 1, 2-8, 9-16, 17-1024
# Defaults to the paper structure: each scale shares same pi.
# Input is to generate_Group, which takes in a vector of scales.
# c(0,1,4,5) says group 0th scale, 1st - 3rd, 4th, 5th -  end scale outputting:
# c(1,2,9,17) (1, 2 - 8, 9 - 16, 17 - 1024)



# Section 5: Estimating effect size ---------------------------------------

# Generating/pre-processing data with NO quantile transform: (no.QT = TRUE)
# setwd("~/WaveQTL/R")
# source("WaveQTL_preprocess_funcs.R")
# data.path = "../data/dsQTL/"
pheno.dat = as.matrix(read.table(paste0(data.path,
                                        "chr17.10160989.10162012.pheno.dat")))
library.read.depth = scan(paste0(data.path, "library.read.depth.dat"))
Covariates = as.matrix(read.table(paste0(data.path, "PC4.dat")))

set.seed(1)
meanR.thresh = 2
res.noQT = WaveQTL_preprocess(Data = pheno.dat, library.read.depth =
                                library.read.depth , Covariates = Covariates, meanR.thresh = meanR.thresh, no.QT = TRUE)
output.path = "../test/dsQTL/"
write.table(res.noQT$WCs, file= paste0(output.path, "WCs.no.QT.txt"), row.names=FALSE,
            col.names = FALSE, quote=FALSE)

## Run 4.2.1 of the guide (testing phenotype vs 24 'near-by' SNPs)
## run stuff through command line - runs the WaveQTL.exe
# Here's some diagnostics...

# 1) Evidence for association between any SNPs
# Here, one snp, with p-val at the bottom (intermediates relate to permutation testing)
pval = as.matrix(read.table("~/Cpp/WaveQTL/test/dsQTL/output/test.fph.pval.txt"))

# 2) LogLR of each of 24 'nearby' SNPs. 
# 1st col = name of SNP
# 2nd col = LogLR for each SNP, ie the ratio of pi vs pi = 0
# 3rd and remaining cols = BFs for each scale,loc, for that SNP
logLR = as.matrix(read.table("~/Cpp/WaveQTL/test/dsQTL/output/test.fph.logLR.txt"))

### ~~~ We will need to change this with HMT! Give 1024 estimates for each SNP, generated by HMT algo
# 3) Max L'hood estimates for pi, for EACH SNP
# (If you were to generalise this model to a linear combo of Betas of different g's each with
# either own pi parameter vectors dictating proportion at each scale, s)
# Note: 24 SNPs, 11 scales
mlePi = as.matrix(read.table("~/Cpp/WaveQTL/test/dsQTL/output/test.fph.pi.txt"))

# 4) Posterior mean of beta effect size, at each s-l, for each SNP
postMean = as.matrix(read.table("~/Cpp/WaveQTL/test/dsQTL/output/test.fph.mean.txt"))

# 5) Posterior variance of effect size beta, at each s-l, for each SNP
postVar = as.matrix(read.table("~/Cpp/WaveQTL/test/dsQTL/output/test.fph.var.txt"))

# Use this script to get_effectSizeInDataSpace.R
source("get_effectSizeinDataSpace.R")

# What this script does:
  # Reads in DWT transform matrix (for the respective data length - 512, 1024 and 2048 pre-calculated)
  # Posterior mean (like in the supp mtl, for each base, sum of wlet x mean across all sl's)
  # Posterior variance (like in the supp mtl, for each base, sum of wlet^2 x mean across all sl's)
  # Looking for either significantly above or significantly below zero:
    # left bounds (data - 3 SDs) > 0
    # right bounds (data + 3 SDs) < 0 
  # Outputs pdf of plot:
    # effect size (y) vs base location (x)
    # pink bars are base location of significant effect (either above or below zero)
    # Three lines: lower and upper bounds, and mean of effect size

# All the above analysis focussed on the 11th SNP, which was the one as identified in 'test.' (multiple SNP analysis
# with permutations) as the one which has the strongest association with the data.
# The above shows effect size of this SNP with the DNase-seq data at the given sites.

