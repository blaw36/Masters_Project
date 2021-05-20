# Test script
print(Sys.time())

#print("Installing rhdf5 package...")
#install.packages("BiocManager", repos = "https://cran.ms.unimelb.edu.au/")
#BiocManager::install("rhdf5")

print(paste(Sys.time(),"Setting parameters...", sep = " - "))
wd.path = "/home/bklaw/paper_data_clean/analysis/simulation/sample_size/simulation_manydsQTL_v1/prepareData/"
setwd(wd.path)
source("/home/bklaw/paper_data_clean/src/R/02_prepare_DNase_funcs.R")
source("/home/bklaw/paper_data_clean/src/R/03_utils.R")

path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/region_01_sel_step1/test.loc"
geno.info.dir.path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno/"
hdf5.data.path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/hg18/dnase/"
hdf5.mapp.path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/hg18/mappability/roger_20bp_mapping_uniqueness.h5"
inds.IDs = scan(paste0("/home/bklaw/paper_data_clean/data/Shim_2014_etc/DNaseI.individuals.oneline.txt"), what="")
output.path = "/home/bklaw/paper_data_clean/analysis/simulation/sample_size/simulation_manydsQTL_v1/data/"

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
