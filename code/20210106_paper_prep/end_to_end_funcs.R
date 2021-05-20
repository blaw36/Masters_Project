run_end_to_end <- function(
  pheno.dat.file
  ,library.read.depth.file
  ,covariates.file
  ,output.path
  ,output.prefix
  ,grouping.file){

## A script which goes end-to-end for WaveQTL and WaveQTL_HMT

# Other scripts which need to be in the same directory:
# WaveQTL_preprocess_funcs.R ; from WaveQTL repository

## Parameters
num_perms = 1000

## Preamble functions

## read functions for WaveQTL preprocess
source("WaveQTL_preprocess_funcs.R")

wavelet_cleaning_wrapper_function_nonRMD <- function(
  pheno.dat
  , output.path
  , output.prefix
  , meanR.thresh = 2
  , library.read.depth
  , Covariates
  , no.QT = TRUE){

  # data.path = "../code/WaveQTL/data/dsQTL/"

  ## set seed
  set.seed(1)

  res = WaveQTL_preprocess(Data = pheno.dat, library.read.depth=library.read.depth , Covariates = Covariates, meanR.thresh = meanR.thresh, no.QT)

  # No quantile transform
  if(no.QT){
    ## for effect size estimation, we need WCs without QT.
    write.table(res$WCs, file= paste0(output.path, output.prefix, "_WCs.no.QT.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
    }else{
    # Else, WCs with QT
    write.table(res$WCs, file= paste0(output.path, output.prefix, "_WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
  }

  cat(res$filtered.WCs, file = paste0(output.path, output.prefix, "_use.txt"))

  # ## produce group information and save it as a file
  # group.info = generate_Group(dim(res$WCs)[2])
  # cat(group.info, file = paste0(output.path, output.prefix, "_group.txt"))

}

## Geno and pheno pre-processing

pre_process_start <- proc.time()

print(paste(Sys.time(), "Cleaning datasets", sep = " - "))

## read functional data
#pheno.dat = as.matrix(read.table(paste0(data.path, "chr17.10160989.10162012.pheno.dat")))
pheno.dat = as.matrix(read.table(pheno.dat.file))
#dim(pheno.dat)
#[1]   70 1024

## read library read depth
#library.read.depth = scan(paste0(data.path, "library.read.depth.dat"))
library.read.depth = scan(library.read.depth.file)
#length(library.read.depth)
#[1] 70

## read Covariates
#Covariates = as.matrix(read.table(paste0(data.path, "PC4.dat")))
Covariates = as.matrix(read.table(covariates.file))
#dim(Covariates)

wavelet_cleaning_wrapper_function_nonRMD(pheno.dat = pheno.dat
  #,output.path = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/length_",effect_length,"/null_data/")
  ,output.path = output.path
  ,output.prefix = output.prefix
  ,library.read.depth = library.read.depth
  ,Covariates = Covariates
  ,no.QT = FALSE)

pre_process_end <- proc.time()
pre_process_elapsed <- pre_process_end - pre_process_start

## Run WaveQTL for ALL NEARBY SNPs with permutations
waveqtl_start <- proc.time()

print(paste(Sys.time(), "WaveQTL", sep = " - "))

if(no.QT){
  setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
    # No HMT
      system(paste0("~/WaveQTL_HMT/WaveQTL"
         ," -gmode 1"
         ," -g /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno70.txt" # genotype (eg ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno)
         ," -p ", paste0(output.path, output.prefix, "_WCs.no.QT.txt")
         ," -group ", grouping.file # grouping file
         ," -u ", paste0(output.path, output.prefix, "_use.txt")
         ," -o paper_waveqtl_",output.prefix # output name
         ," -f 1024"
         ," -numPerm ", num_perms
         ," -fph 3"
         ," > out.",output.prefix," 2> err.",output.prefix)
      ,show.output.on.console = F)
    setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
    }else{
      setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
    # No HMT
      system(paste0("~/WaveQTL_HMT/WaveQTL"
         ," -gmode 1"
         ," -g /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno70.txt" # genotype (eg ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno)
         ," -p ", paste0(output.path, output.prefix, "_WCs.txt")
         ," -group ", grouping.file # grouping file
         ," -u ", paste0(output.path, output.prefix, "_use.txt")
         ," -o paper_waveqtl_",output.prefix # output name
         ," -f 1024"
         ," -numPerm ", num_perms
         ," -fph 3"
         ," > out.",output.prefix," 2> err.",output.prefix)
      ,show.output.on.console = F)
    setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
  }

  waveqtl_end <- proc.time()
  waveqtl_elapsed <- waveqtl_end - waveqtl_start

## Run WaveQTL_HMT for ALL NEARBY SNPs with permutations
waveqtl_hmt_start <- proc.time()

print(paste(Sys.time(), "WaveQTL_HMT", sep = " - "))

if(no.QT){
  setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
    # HMT
      system(paste0("~/WaveQTL_HMT/WaveQTL"
         ," -gmode 1"
         ," -g /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno70.txt" # genotype (eg ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno)
         ," -p ", paste0(output.path, output.prefix, "_WCs.no.QT.txt")
         ," -group ", grouping.file, # grouping file
         ," -u ", paste0(output.path, output.prefix, "_use.txt"),
         ," -o paper_waveqtl_",output.prefix # output name
         ," -f 1024"
         ," -numPerm ", num_perms
         ," -fph 3"
         ," -hmt 1"
         ," > out.",output.prefix," 2> err.",output.prefix)
      ,show.output.on.console = F)
    setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
    }else{
      setwd("~/Cpp/WaveQTL_HMT/test/dsQTL/")
    # HMT
      system(paste0("~/WaveQTL_HMT/WaveQTL"
         ," -gmode 1"
         ," -g /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno70.txt" # genotype (eg ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno)
         ," -p ", paste0(output.path, output.prefix, "_WCs.txt")
         ," -group ", grouping.file, # grouping file
         ," -u ", paste0(output.path, output.prefix, "_use.txt"),
         ," -o paper_waveqtl_",output.prefix # output name
         ," -f 1024"
         ," -numPerm ", num_perms
         ," -fph 3"
         ," -hmt 1"
         ," > out.",output.prefix," 2> err.",output.prefix)
      ,show.output.on.console = F)
    setwd("~/../Dropbox/Uni Stuff - Masters/Research Project/Masters_Project_Git/")
  }

  waveqtl_hmt_end <- proc.time()
  waveqtl_hmt_elapsed <- waveqtl_hmt_end - waveqtl_hmt_start

}