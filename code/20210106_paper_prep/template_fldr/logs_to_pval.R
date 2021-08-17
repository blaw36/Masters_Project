library(data.table)
library(stringr)
path = "output/"

hmt_files = list.files(path, pattern = "\\.hmt\\.")
hmt_files = grep(pattern = "\\.pval\\.txt",x = hmt_files,value = T)

hmt_dat = data.table(
  chr = integer()
  , site = integer()
  , strong_snp = character()
  , perms = integer()
  , above_thold = integer()
  , pval = numeric())

for(i in 1:length(hmt_files)){
  file_name_split = stringr::str_split(hmt_files[i],pattern = "\\.",simplify = T)
  hmt_chr = as.integer(file_name_split[,2])
  hmt_site = as.integer(file_name_split[,3])
  hmt_pval = read.table(paste0(path,hmt_files[i]), stringsAsFactors = F)
  # hmt_dat = rbind(hmt_dat, as.list(hmt_pval$V1))
  hmt_dat = rbind(hmt_dat, list(chr = hmt_chr
                                , site = hmt_site
                                , strong_snp = hmt_pval$V1[1]
                                , perms = as.integer(hmt_pval$V1[2])
                                , above_thold = as.integer(hmt_pval$V1[3])
                                , pval = as.numeric(hmt_pval$V1[4])))
}

nohmt_files = list.files(path, pattern = "\\.nohmt\\.")
nohmt_files = grep(pattern = "\\.pval\\.txt",x = nohmt_files,value = T)

nohmt_dat = data.table(
  chr = integer()
  , site = integer()
  , strong_snp = character()
  , perms = integer()
  , above_thold = integer()
  , pval = numeric())

for(i in 1:length(nohmt_files)){
  file_name_split = stringr::str_split(nohmt_files[i],pattern = "\\.",simplify = T)
  nohmt_chr = as.integer(file_name_split[,2])
  nohmt_site = as.integer(file_name_split[,3])
  nohmt_pval = read.table(paste0(path,nohmt_files[i]), stringsAsFactors = F)
  # nohmt_dat = rbind(hmt_dat, as.list(nohmt_pval$V1))
  nohmt_dat = rbind(nohmt_dat, list(chr = nohmt_chr
                                    , site = nohmt_site
                                    , strong_snp = nohmt_pval$V1[1]
                                    , perms = as.integer(nohmt_pval$V1[2])
                                    , above_thold = as.integer(nohmt_pval$V1[3])
                                    , pval = as.numeric(nohmt_pval$V1[4])))
}

saveRDS(hmt_dat,"hmt_pvals.RDS", compress = T)
saveRDS(nohmt_dat,"nohmt_pvals.RDS", compress = T)
