run_end_to_end <- function(
  pheno.dat.file
  ,library.read.depth.file
  ,covariates.file
  ,output.path
  ,site
  ,chrIX
  ,grouping.file = NULL
  ,num_perms
  ,no.QT = TRUE){

  ## A script which goes end-to-end for WaveQTL and WaveQTL_HMT

  # Other scripts which need to be in the same directory:
  # WaveQTL_preprocess_funcs.R ; from WaveQTL repository

  # ## Parameters
  # num_perms = 1000

  ## Preamble functions

  ## read functions for WaveQTL preprocess
  # source("~/Cpp/WaveQTL_HMT_wperm/R/WaveQTL_preprocess_funcs.R")
  source("/home/bklaw/WaveQTL_HMT_wperm/R/WaveQTL_preprocess_funcs.R")

  output.prefix = paste0("chr.",chrIX,".",site)

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
      cat(res$filtered.WCs, file = paste0(output.path, output.prefix, "_use.no.QT.txt"))
    }else{
      # Else, WCs with QT
      write.table(res$WCs, file= paste0(output.path, output.prefix, "_WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
      cat(res$filtered.WCs, file = paste0(output.path, output.prefix, "_use.txt"))
    }

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
                                           ,output.path = output.path
                                           ,output.prefix = output.prefix
                                           ,library.read.depth = library.read.depth
                                           ,Covariates = Covariates
                                           ,no.QT = no.QT)

  pre_process_end <- proc.time()
  pre_process_elapsed <- pre_process_end - pre_process_start

  ## Run WaveQTL for ALL NEARBY SNPs with permutations
  waveqtl_start <- proc.time()

  print(paste(Sys.time(), "WaveQTL", sep = " - "))

  if(no.QT){
    # No HMT
    system(paste0("../WaveQTL" #"~/Cpp/WaveQTL_HMT_wperm/WaveQTL"
                  ," -gmode 1"
                  ," -g /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno"
                  ," -p ", paste0(output.path, output.prefix, "_WCs.no.QT.txt")
                  #," -group ", grouping.file # grouping file
                  ," -u ", paste0(output.path, output.prefix, "_use.no.QT.txt")
                  ," -o ", output.prefix,".nohmt.no.QT" # output name
                  ," -f 1024"
                  ," -numPerm ", num_perms
                  ," -fph 3")
           )
  }else{
    # No HMT
    system(paste0("../WaveQTL" #"~/Cpp/WaveQTL_HMT_wperm/WaveQTL"
                  ," -gmode 1"
                  ," -g /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno"
                  ," -p ", paste0(output.path, output.prefix, "_WCs.txt")
                  #," -group ", grouping.file # grouping file
                  ," -u ", paste0(output.path, output.prefix, "_use.txt")
                  ," -o ", output.prefix,".nohmt" # output name
                  ," -f 1024"
                  ," -numPerm ", num_perms
                  ," -fph 3")
           )
  }

  waveqtl_end <- proc.time()
  waveqtl_elapsed <- waveqtl_end - waveqtl_start

  ## Run WaveQTL_HMT for ALL NEARBY SNPs with permutations
  waveqtl_hmt_start <- proc.time()

  print(paste(Sys.time(), "WaveQTL_HMT", sep = " - "))

  if(no.QT){
    # HMT
    system(paste0("../WaveQTL" #"~/Cpp/WaveQTL_HMT_wperm/WaveQTL"
                  ," -gmode 1"
                  ," -g /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno"
                  ," -p ", paste0(output.path, output.prefix, "_WCs.no.QT.txt")
                  #," -group ", grouping.file # grouping file
                  ," -u ", paste0(output.path, output.prefix, "_use.no.QT.txt")
                  ," -o ", output.prefix,".hmt.no.QT" # output name
                  ," -f 1024"
                  ," -numPerm ", num_perms
                  ," -fph 3"
                  ," -hmt 1")
           )
  }else{
    # HMT
    system(paste0("../WaveQTL" #"~/Cpp/WaveQTL_HMT_wperm/WaveQTL"
                  ," -gmode 1"
                  ," -g /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno"
                  ," -p ", paste0(output.path, output.prefix, "_WCs.txt")
                  #," -group ", grouping.file # grouping file
                  ," -u ", paste0(output.path, output.prefix, "_use.txt")
                  ," -o ", output.prefix,".hmt" # output name
                  ," -f 1024"
                  ," -numPerm ", num_perms
                  ," -fph 3"
                  ," -hmt 1")
           )
  }

  waveqtl_hmt_end <- proc.time()
  waveqtl_hmt_elapsed <- waveqtl_hmt_end - waveqtl_hmt_start

  return(timings_list = list(
    pre_process_elapsed = pre_process_elapsed
    ,waveqtl_elapsed = waveqtl_elapsed
    ,waveqtl_hmt_elapsed = waveqtl_hmt_elapsed
  ))

}

# # Sample usage
# t <- run_end_to_end(
#   pheno.dat.file = "~/Cpp/WaveQTL_HMT_wperm/data/dsQTL/chr17.10160989.10162012.pheno.dat"
#   ,library.read.depth.file = "~/Cpp/WaveQTL_HMT_wperm/data/dsQTL/library.read.depth.dat"
#   ,covariates.file = "~/Cpp/WaveQTL_HMT_wperm/data/dsQTL/PC4.dat"
#   ,output.path = "~/Cpp/WaveQTL_HMT_wperm/batch_test/"
#   ,output.prefix = "sample_site"
#   ,num_perms = 100
#   ,no.QT = TRUE)


run_end_to_end_no_perm <- function(
  pheno.dat.file
  ,library.read.depth.file
  ,covariates.file
  ,output.path
  ,site
  ,chrIX
  ,grouping.file = NULL
  ,num_perms
  ,no.QT = TRUE){

  ## A script which goes end-to-end for WaveQTL and WaveQTL_HMT

  # Other scripts which need to be in the same directory:
  # WaveQTL_preprocess_funcs.R ; from WaveQTL repository

  # ## Parameters
  # num_perms = 1000

  ## Preamble functions

  ## read functions for WaveQTL preprocess
  # source("~/Cpp/WaveQTL_HMT_wperm/R/WaveQTL_preprocess_funcs.R")
  source("/home/bklaw/WaveQTL_HMT_wperm/R/WaveQTL_preprocess_funcs.R")

  output.prefix = paste0("chr.",chrIX,".",site)

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
      cat(res$filtered.WCs, file = paste0(output.path, output.prefix, "_use.no.QT.txt"))
    }else{
      # Else, WCs with QT
      write.table(res$WCs, file= paste0(output.path, output.prefix, "_WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
      cat(res$filtered.WCs, file = paste0(output.path, output.prefix, "_use.txt"))
    }

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
                                           ,output.path = output.path
                                           ,output.prefix = output.prefix
                                           ,library.read.depth = library.read.depth
                                           ,Covariates = Covariates
                                           ,no.QT = no.QT)

  pre_process_end <- proc.time()
  pre_process_elapsed <- pre_process_end - pre_process_start

  ## Run WaveQTL for ALL NEARBY SNPs with permutations
  waveqtl_start <- proc.time()

  print(paste(Sys.time(), "WaveQTL", sep = " - "))

  if(no.QT){
    # No HMT
    system(paste0("../WaveQTL" #"~/Cpp/WaveQTL_HMT_wperm/WaveQTL"
                  ," -gmode 1"
                  ," -g /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno"
                  ," -p ", paste0(output.path, output.prefix, "_WCs.no.QT.txt")
                  #," -group ", grouping.file # grouping file
                  ," -u ", paste0(output.path, output.prefix, "_use.no.QT.txt")
                  ," -o ", output.prefix,".nohmt.no.QT" # output name
                  ," -f 1024"
                  ," -fph 1")
           )
  }else{
    # No HMT
    system(paste0("../WaveQTL" #"~/Cpp/WaveQTL_HMT_wperm/WaveQTL"
                  ," -gmode 1"
                  ," -g /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno"
                  ," -p ", paste0(output.path, output.prefix, "_WCs.txt")
                  #," -group ", grouping.file # grouping file
                  ," -u ", paste0(output.path, output.prefix, "_use.txt")
                  ," -o ", output.prefix,".nohmt" # output name
                  ," -f 1024"
                  ," -fph 1")
           )
  }

  waveqtl_end <- proc.time()
  waveqtl_elapsed <- waveqtl_end - waveqtl_start

  ## Run WaveQTL_HMT for ALL NEARBY SNPs with permutations
  waveqtl_hmt_start <- proc.time()

  print(paste(Sys.time(), "WaveQTL_HMT", sep = " - "))

  if(no.QT){
    # HMT
    system(paste0("../WaveQTL" #"~/Cpp/WaveQTL_HMT_wperm/WaveQTL"
                  ," -gmode 1"
                  ," -g /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno"
                  ," -p ", paste0(output.path, output.prefix, "_WCs.no.QT.txt")
                  #," -group ", grouping.file # grouping file
                  ," -u ", paste0(output.path, output.prefix, "_use.no.QT.txt")
                  ," -o ", output.prefix,".hmt.no.QT" # output name
                  ," -f 1024"
                  ," -fph 1"
                  ," -hmt 1")
           )
  }else{
    # HMT
    system(paste0("../WaveQTL" #"~/Cpp/WaveQTL_HMT_wperm/WaveQTL"
                  ," -gmode 1"
                  ," -g /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno"
                  ," -p ", paste0(output.path, output.prefix, "_WCs.txt")
                  #," -group ", grouping.file # grouping file
                  ," -u ", paste0(output.path, output.prefix, "_use.txt")
                  ," -o ", output.prefix,".hmt" # output name
                  ," -f 1024"
                  ," -fph 1"
                  ," -hmt 1")
           )
  }

  waveqtl_hmt_end <- proc.time()
  waveqtl_hmt_elapsed <- waveqtl_hmt_end - waveqtl_hmt_start

  return(timings_list = list(
    pre_process_elapsed = pre_process_elapsed
    ,waveqtl_elapsed = waveqtl_elapsed
    ,waveqtl_hmt_elapsed = waveqtl_hmt_elapsed
  ))

}