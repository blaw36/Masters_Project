#!/usr/bin/env Rscript

#source("/home/bklaw/WaveQTL_HMT_wperm/R/end_to_end_funcs_spartan.R")

print(getwd())
curr_wd = getwd()

# source(paste0(curr_wd,"/end_to_end_funcs_spartan_no_wc.R"))
source(paste0(curr_wd,"/end_to_end_funcs_spartan.R"))

# REAL
# Read in chr/sites
args = (commandArgs(TRUE)) 
eval(parse(text=args[[1]]))
print(args[1])
eval(parse(text=args[[2]]))
print(args[2])

chrIX = args[1]
site = args[2]

print(paste0("Chromosome: ", chrIX, ", Site: ", site))

# See if phenotype data exists
pheno_data_path = paste0("/home/bklaw/paper_data_clean/analysis/simulation/sample_size/simulation_manydsQTL_v1/data/chr", chrIX, ".pheno.dat.", site)

if(file.exists(pheno_data_path)){
    t <- run_end_to_end(
        pheno.dat.file = pheno_data_path
        ,library.read.depth.file = "/home/bklaw/WaveQTL_HMT_wperm/data/dsQTL/library.read.depth.dat"
        ,covariates.file = "/home/bklaw/WaveQTL_HMT_wperm/data/dsQTL/PC4.dat"
        #,output.path = "/home/bklaw/WaveQTL_HMT_wperm/20210617_run1/wcs/"
        ,output.path = paste0(curr_wd,"/wcs/")
        ,site = site
        ,chrIX = chrIX
        #,output.prefix = paste0("chr.",chrIX,".",site)
        ,num_perms = 10000
        ,no.QT = FALSE)

    #print(paste0("Pre-process: ", round(t$pre_process_elapsed[3],2)))
    print(paste0("WaveQTL: ", round(t$waveqtl_elapsed[3],2)))
    print(paste0("WaveQTL HMT: ", round(t$waveqtl_hmt_elapsed[3],2)))

}else{
    print("No data")
}