# Script for preparing data for simulations, adapted for use on UniMelb cluster and for data in the UniMelb cluster.
# We will be doing this for 50k sites, well beyond the 578 dsQTLs identified for simulations as below.
# This is so that we can run WaveQTL and WaveQTL-HMT on these sites to identify the benefit of WaveQTL-HMT on sites which are not
# just tailored for/selected for WaveQTL simulation.
# https://github.com/heejungshim/multiscale_analysis/blob/master/analysis/simulation/sample_size/simulation_manydsQTL_v1/prepareData/script/prepare.data.for.simulation.R

#!/usr/bin/env Rscript

## Aim : This file contains Rscript to prepare phenotype and genotype data for 578 dsQTLs for simulations. The 578 dsQTLs are idetinfied by Shim and Stephens 2014 (see its Supplementary Materials for details of those dsQTLs). This information will be used in simulation.
## I modified two scripts: "/mnt/lustre/home/shim/wavelets/revision/code/simulation.explore.578.sites.R", "/mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_footprint/code/better.effect.size.R", and `~/multiscale_analysis/analysis/simulation/sample_size/simulation_578/code/prepare.data.for.simulation.R'
##
## Usage R CMD BATCH --no-save --no-restore "--args ss=$SGE_TASK_ID" prepare.data.for.simulation.R
## See /mnt/lustre/home/shim/multiscale_analysis/analysis/simulation/sample_size/simulation_manydsQTL_v1/com/com.Data.sh
##
##
## Copyright (C) 2014 Heejung Shim
##
## License: GPL3+

## A few things to note.
# Chromosome naming conventions; ChrNum.Region, eg. Chr1.1 is Num 1, Region 1.
# SNP location:
# Chr1.702091 A G 0.018 0 702091
# Ie. Name (with location), minor alleles, minor allele freq, something, then the loc (again)


## ss = 305
# TODO: this errors out for me
# Note that these are runtime commands for when running on the server, perhaps.
args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))

print(paste(Sys.time(),"Setting parameters...", sep = " - "))

# setwd("/mnt/lustre/home/shim/multiscale_analysis")
# # TODO: Figure out if this file exists in the repo?
# # [HJ] It’s this repository here: https://github.com/heejungshim/multiscale_analysis/
# multiscale.analysis.repodir <- scan(".multiscale_analysis.repodir.txt", what=character())
# # TODO: Figure out if this file exists in the repo?
# # [HJ] It’s this repository here: https://github.com/heejungshim/WaveQTL
# WaveQTL.repodir <- scan(".WaveQTL.repodir.txt", what=character())

## BL:
# After a bit of thinking, I think what HJ is saying is that 'multiscale.analysis.repodir' was some dynamically typed path or file
# which is actually just the directory of wherever your 'multiscale_analysis' repo is stored. For me on my computer, it's currently:
# Multiscale: "../../../../Documents/R/multiscale_analysis/"
# WaveQTL: "../../../../Documents/Cpp/WaveQTL/"

# Current location on Spartan
multiscale.analysis.repodir <- "/home/bklaw/paper_data_clean/"
WaveQTL.repodir <- "/home/bklaw/paper_data_clean/"

## set working directory
wd.path = paste0(multiscale.analysis.repodir, "/analysis/simulation/sample_size/simulation_manydsQTL_v1/prepareData/")
setwd(wd.path)


## set path to files
## Path to directory which contain DNase-seq data as hdf5 format,
# hdf5.data.path = "/mnt/lustre/data/internal/genome_db/hg18/dnase/"
hdf5.data.path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/hg18/dnase/"

## Path to mappability information as hdf5 format
# hdf5.mapp.path  = "/mnt/lustre/data/internal/genome_db/hg18/mappability/roger_20bp_mapping_uniqueness.h5"
hdf5.mapp.path  = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/hg18/mappability/roger_20bp_mapping_uniqueness.h5"

## path to directory which contains information on SNPs located region of interest
# geno.info.dir.path = "/mnt/lustre/home/shim/wavelets/data/DNase/geno_01_step1/geno/"
geno.info.dir.path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno/"

## path to directory which contains genotype data
# geno.dir.path = "/mnt/lustre/home/shim/wavelets/data/DNase/geno_01_step1/geno_maf/"
geno.dir.path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/"
# BL: This is a subset of the above’s /geno/, filtered down based on 4th column > 0.05 (minor allele frequency)
# The files here are maf_chr1.1.geno (maf chromosome 1.1)

## path to directory which contains location information on 578 sites
# locus.path = "/mnt/lustre/home/shim/wavelets/data/DNase/region_01_sel_step1/"
locus.path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/region_01_sel_step1/"
# This contains information on 50K randomly selected regions. 578 regions are a subset of 50K regions.
# Don’t really need the info on the 578 regions as that was WaveQTL. We want to go further here (to the 50k).


## path to output directory
output.path = paste0(multiscale.analysis.repodir, "/analysis/simulation/sample_size/simulation_manydsQTL_v1/data/")

## read a list of individual IDs
# It contains YRI 70 individuals’ IDs.
# inds.IDs = scan(paste0(WaveQTL.repodir, "/data/Shim_2014_etc/DNaseI.individuals.oneline.txt"), what="")
inds.IDs = scan(paste0(WaveQTL.repodir, "/data/Shim_2014_etc/DNaseI.individuals.oneline.txt"), what="")

## read functions to read DNase data and preprocess
# TODO: CHANGE THESE PATHS TO POINT TO WHERE THE APPROPRIATE FUNCTIONS ARE ON THE UNIMELB SERVER

# This function (read.DNase.data) is documented in comments_on_mycode_2.docx
# source(paste0(multiscale.analysis.repodir, "/src/R/prepare.DNase.funcs.R"))
source(paste0(multiscale.analysis.repodir, "/src/R/02_prepare_DNase_funcs.R"))

# The above function contains get.counts.h5, which is called in 'read.DNase.data'
# source(paste0(multiscale.analysis.repodir, "/src/R/utils.R"))
source(paste0(multiscale.analysis.repodir, "/src/R/03_utils.R"))

## read information on 578 sites
# data = read.table("/mnt/lustre/home/shim/wavelets/revision/etc/simu.578.sites.txt", header=T)
data = read.table("/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/etc/simu.578.sites.txt", header=T)
##names(data)
##[1] "index"       "chr"         "site"        "genoIX"      "FDR.10.wave"

# We won’t focus on 578 sites.
# “chr” indicates chromosome of the region.
# “site” indicates row number of the region in “chrX.loc” file in /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/region_01_sel_step1/
# “genoIX” indicates the index (in chrX.X.geno file in /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/)
# of the SNP with strongest association – we probably won’t use this information because we are supposed to run for all SNPs (in cis area with MAF > 0.05).


# chr.list = data$chr
# site.list = data$site
# genoIX.list = data$genoIX


## for each dsQTL, let's prepare data
# for(ss in 1:578){

# chr = chr.list[ss]
# site = site.list[ss]
# genoIX = genoIX.list[ss]

## [BL]; we can run the above loop for all the 50k sites, with the data found in:
# /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/region_01_sel_step1
# Site will be row numbers of each file...ie. row 1 of chr1 will be site 1
# This will correspond to chr1.1.geno, for example (which is chr1, site1)
# chr will be 1-22
# genoIX is not required

## Note that we'll use chr18.loc.back, NOT chr18.loc (as using .back sums up to 50k rows)

## read location information
# [HJ] the following code extracts chromosome name (“chrIX”),
# start position of region (“locus.start”), and end position of region (“locus.end”)
# from “chrX.loc” files in /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/region_01_sel_step1/

# TODO: Remove and use the commented out path later. This is hardcoded as a test location file
path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/region_01_sel_step1/test.loc"
# path = paste0(locus.path, "chr", chr, ".loc")

loc_dat = read.table(path, as.is = TRUE)

num_rows = nrow(loc_dat)

print(paste(Sys.time(),"Reading in data...", sep = " - "))
for(site in 1:1){
  #for(site in num_rows){
  print(paste(Sys.time(), paste("Site:",site), sep = " - "))

  chrIX = loc_dat[site,1]
  locus.start = loc_dat[site,2]
  locus.end = loc_dat[site, 3] - 1

  ## path to genotype information (to correct for cutting preference)
  # geno.info.path = paste0(geno.info.dir.path, "maf_chr", chr, ".", site, ".geno")
  geno.info.path = paste0(geno.info.dir.path, "maf_", chrIX, ".", site, ".geno")
  # [HJ] these files (“maf_chrX.X.geno” in /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno/)
  # contain information (e.g, position) on all SNPs in cis-area for a given region.
  # We will use this position information to handle cutting preference inside “read.DNase.data” function.

  print(paste(Sys.time(), "Running read.DNase.data function", sep = " - "))
  ## run function to read DNase data
  # [HJ] this function is implemented in https://github.com/heejungshim/multiscale_analysis/blob/master/src/R/prepare.DNase.funcs.R
  res = read.DNase.data(hdf5.data.path = hdf5.data.path, hdf5.mapp.path = hdf5.mapp.path, geno.info.path = geno.info.path, inds.IDs = inds.IDs, chrIX = chrIX, locus.start = locus.start  , locus.end = locus.end)

  # [HJ] DNase-seq counts from the function “read.DNase.data” (res$DNase.dat) are average from two strands.
  # So they are not necessary counts, so using ceiling function, we made them counts data.
  phenoD = ceiling(res$DNase.dat)

  ## output information to files
  path.output = paste0(output.path, chrIX, ".pheno.dat.", site)
  print(paste(Sys.time(), paste("Writing output to:",path.output), sep = " - "))
  write.table(phenoD, path.output, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# [BL] - not required
# ## get genotype informaiton
# # [HJ] the following code is to obtain genotypes for SNP with strongest association. But we won’t need this information.
# geno.path = paste0(geno.dir.path, "chr", chr, ".", site, ".geno")
# genoF = read.table(geno.path, as.is = TRUE)
# genoD = as.numeric(genoF[genoIX, 4:73])
#
# ## output information to files
# path.output = paste0(output.path, "pheno.dat.", ss)
# write.table(phenoD, path.output, row.names = FALSE, col.names = FALSE, quote = FALSE)
#
# path.output = paste0(output.path, "orig.geno.dat.", ss)
# cat(genoD, file = path.output)
