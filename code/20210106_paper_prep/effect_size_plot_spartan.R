#!/usr/bin/env Rscript
## Run this in the base run folder (NOT the batch scripts folder)

source("/home/bklaw/WaveQTL_HMT_wperm/R/end_to_end_funcs_spartan.R")
source("/home/bklaw/WaveQTL_HMT_wperm/R/effect_size_plot_funcs_spartan.R")

print(getwd())
curr_wd = getwd()

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
  t <- run_end_to_end_no_perm(
    pheno.dat.file = pheno_data_path
    ,library.read.depth.file = "/home/bklaw/WaveQTL_HMT_wperm/data/dsQTL/library.read.depth.dat"
    ,covariates.file = "/home/bklaw/WaveQTL_HMT_wperm/data/dsQTL/PC4.dat"
    ,output.path = paste0(curr_wd,"/wcs/")
    ,site = site
    ,chrIX = chrIX
    ,no.QT = TRUE)

  effect_size_plot_data <- effect_size_plot(
    site = site
    ,chrIX = chrIX
    ,wqtl_qt_output_path = paste0(curr_wd,"/output/")
    ,wqtl_no_qt_output_path = paste0(curr_wd,"/output/")
    ,plot_output_path = paste0(curr_wd,"/plots/")
  )
  saveRDS(effect_size_plot_data, paste0(curr_wd,"/plots/chr.",chrIX,".",site,".plot.data.RDS"), compress = TRUE)

}else{
  print("No data")
}
