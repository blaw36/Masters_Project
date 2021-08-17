library(stringr)
library(data.table)

# Compile all the possible locs
locus.path = "/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/region_01_sel_step1/"
chr = list.files(locus.path,"*.loc")
site_list = data.table(
	chr = integer()
	,site = integer())
for(i in chr){
	path = paste0(locus.path, i)
	loc_dat = read.table(path, as.is = TRUE)
	num_sites = nrow(loc_dat)
	site_list = rbind(site_list,
		data.table(
			chr = rep(i, num_sites)
			,site = 1:num_sites))
}
site_list[,"chr" := gsub("chr","",chr)]
site_list[,"chr" := gsub(".loc","",chr)]

site_list[,"chr" := as.integer(chr)]
site_list[,"site" := as.integer(site)]
setorder(site_list,chr,site)

# Randomise data
set.seed(1)
# site_allocations = split(sample(50000,50000,replace = F), as.factor(1:5))
site_allocations = split(sample(nrow(site_list),nrow(site_list),replace = F), as.factor(1:5))
grp1 = site_list[site_allocations[[1]]]
grp2 = site_list[site_allocations[[2]]]
grp3 = site_list[site_allocations[[3]]]
grp4 = site_list[site_allocations[[4]]]
grp5 = site_list[site_allocations[[5]]]

setorder(grp1, chr, site)
setorder(grp2, chr, site)
setorder(grp3, chr, site)
setorder(grp4, chr, site)
setorder(grp5, chr, site)

# Create the file names
grp1[,"filename" := paste0("chr", chr, ".pheno.dat.", site)]
grp2[,"filename" := paste0("chr", chr, ".pheno.dat.", site)]
grp3[,"filename" := paste0("chr", chr, ".pheno.dat.", site)]
grp4[,"filename" := paste0("chr", chr, ".pheno.dat.", site)]
grp5[,"filename" := paste0("chr", chr, ".pheno.dat.", site)]

save_path = "/home/bklaw/paper_data_clean/site_splits/"
write.table(x = grp1, file = paste0(save_path,"sites_grp1.txt"), row.names = F)
write.table(x = grp2, file = paste0(save_path,"sites_grp2.txt"), row.names = F)
write.table(x = grp3, file = paste0(save_path,"sites_grp3.txt"), row.names = F)
write.table(x = grp4, file = paste0(save_path,"sites_grp4.txt"), row.names = F)
write.table(x = grp5, file = paste0(save_path,"sites_grp5.txt"), row.names = F)