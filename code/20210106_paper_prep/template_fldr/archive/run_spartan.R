#!/usr/bin/env Rscript

source("/home/bklaw/WaveQTL_HMT_wperm/R/end_to_end_funcs_spartan.R")

print(getwd())

# # SAMPLE
# t <- run_end_to_end(
#   pheno.dat.file = "/home/bklaw/WaveQTL_HMT_wperm/data/dsQTL/chr17.10160989.10162012.pheno.dat"
#   ,library.read.depth.file = "/home/bklaw/WaveQTL_HMT_wperm/data/dsQTL/library.read.depth.dat"
#   ,covariates.file = "/home/bklaw/WaveQTL_HMT_wperm/data/dsQTL/PC4.dat"
#   ,output.path = "/home/bklaw/WaveQTL_HMT_wperm/20210616_test/wcs/"
#   ,output.prefix = "sample_site"
#   ,num_perms = 100
#   ,no.QT = TRUE)

# REAL
# Read in chr/sites
args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
print(args[1])

chr_sites = read.table(paste0("/home/bklaw/paper_data_clean/site_splits/grp1/chunk",args[1],".txt"), header = T)

for(i in 1:nrow(chr_sites)){
  chrIX = chr_sites[i,"chr"]
  site = chr_sites[i,"site"]
  print(paste0("Chromosome: ", chrIX, ", Site: ", site))

  # See if phenotype data exists
  pheno_data_path = paste0("/home/bklaw/paper_data_clean/analysis/simulation/sample_size/simulation_manydsQTL_v1/data/chr", chrIX, ".pheno.dat.", site)

  if(file.exists(pheno_data_path)){
    t <- run_end_to_end(
      pheno.dat.file = pheno_data_path
      ,library.read.depth.file = "/home/bklaw/WaveQTL_HMT_wperm/data/dsQTL/library.read.depth.dat"
      ,covariates.file = "/home/bklaw/WaveQTL_HMT_wperm/data/dsQTL/PC4.dat"
      ,output.path = "/home/bklaw/WaveQTL_HMT_wperm/20210617_run1/wcs/"
      ,site = site
      ,chrIX = chrIX
      #,output.prefix = paste0("chr.",chrIX,".",site)
      ,num_perms = 10000
      ,no.QT = FALSE)

    print(paste0("Pre-process: ", round(t$pre_process_elapsed[3],2)))
    print(paste0("WaveQTL: ", round(t$waveqtl_elapsed[3],2)))
    print(paste0("WaveQTL HMT: ", round(t$waveqtl_hmt_elapsed[3],2)))

  }else{
    print("No data")
  }

}
