#!/usr/bin/env Rscript

## Data space plot script for Spartan
## Run this in the base run folder (NOT the batch scripts folder)

print(getwd())
curr_wd = getwd()

# Read in chr/sites
args = (commandArgs(TRUE))
eval(parse(text=args[[1]]))
print(args[1])
eval(parse(text=args[[2]]))
print(args[2])
#eval(parse(text=args[[3]]))
#print(args[3])

chrIX = args[1]
site = args[2]
snp_name = args[3]
print(paste0("Chromosome: ", chrIX, ", Site: ", site))
print(paste0("Nearby snp: ", snp_name))

data_space_plot <- function(
  geno_data_path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/"
  ,pheno_data_path = "/home/bklaw/paper_data_clean/analysis/simulation/sample_size/simulation_manydsQTL_v1/data/"
  ,site
  ,chrIX
  ,snp_name
  ,plot_output_path = paste0(curr_wd,"/plots_orig_data/")
){
  library(data.table)
  library(reshape2)

  # Load genotype (SNP) data
  eg_geno <- as.matrix(read.table(paste0(geno_data_path,"chr",chrIX,".",site,".geno")))

  # Read in phenotype data
  pheno.dat = as.matrix(read.table(paste0(pheno_data_path,"chr", chrIX, ".pheno.dat.", site)))

  # Bucketed up based on SNP value ----
  # Read in genotype, round to nearest bucket
  snpIX = which(eg_geno[,1] == snp_name)
  if(length(snpIX) == 0){
    stop(paste0("No such SNP ",snp_name," found"))
  }

  # Data from the 70 individuals, bucketed by rounded genotype value
  geno_nearby <- as.numeric(eg_geno[snpIX,][4:73])
  geno_nearby_bucketed <- round(geno_nearby)
  num_buckets = sort(unique(geno_nearby_bucketed))
  bucket_means = list()

  n = 1
  for(i in num_buckets){
    num_ind_in_bucket = sum(geno_nearby_bucketed == i)
    if(num_ind_in_bucket > 1){
      # If num obs > 1
      bucket_means[[n]] = data.frame(mean = apply(pheno.dat[geno_nearby_bucketed == i, ],2,mean), loc = 1:1024)
    }else if(num_ind_in_bucket == 1){
      # If num obs == 1
      bucket_means[[n]] = data.frame(mean = pheno.dat[geno_nearby_bucketed == i, ], loc = 1:1024)
    }
    n = n + 1
  }

  all_means <- data.table::rbindlist(bucket_means)
  y_bounds <- c(0,max(all_means$mean))*1.05

  plot.new()
  # Open file
  pdf(paste0(plot_output_path,"chr.",chrIX,".",site,".data.plot.PDF"), width = 8, height = 3)
  # pdf("images/ch1_seq_mean_70indivs_bucketed.pdf", width = 12, height = 7, pointsize = 20)
  # Create plot
  plot(1,1
       ,type = "l"
       ,xlab = "Base location"
       ,ylab = "Normalised count"
       ,main = paste0("Chr: ",chrIX,", site: ",site,", SNP: ",snp_name)
       ,xaxt = "n"
       ,xlim = c(0,dim(pheno.dat)[2])
       ,ylim = y_bounds)
  axis(1, at = seq(0,1024,128),labels = seq(0,1024,128))

  # Loop through the lines
  for(i in 1:length(bucket_means)){
    lines(bucket_means[[i]]$mean, col = i)
  }

  legend(x = 750, y = 0.7, title = "Genotype value",legend = c(2:0)
         ,cex=0.8
         ,col = 1:length(bucket_means) # ,col=c("red","green","black")
         ,lty = rep(1,3)
         ,box.lty=0)
  # Close the file
  dev.off()
}

data_space_plot(
  site = site
  ,chrIX = chrIX
  ,snp_name = snp_name
)
