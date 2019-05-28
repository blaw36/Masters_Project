## Few checks to look for equivalence between the WaveQTL manual workings.

# Read in the files with almost identical naming structure.
# Some won't read - the logs and snpdata, which we're not concerned about anyway.
path <- "Cpp/WaveQTL/test/dsQTL/output/"
files <- list.files(path)
test_files <- grep("test\\.no\\.QT\\."
                   ,grep("test\\.",files,value = T)
                   ,value = T, invert = T)
testNoQT_files <- grep("test\\.no\\.QT\\.",files,value = T)
test1_files <- grep("test1\\.",files,value = T)
test2_files <- grep("test2\\.",files,value = T)

test_list <- list()
for(i in test_files){
  tmp_name = stringr::str_split_fixed(string = i, ".fph.", n = 2)[2]
  test_list[[tmp_name]] = as.matrix(read.table(paste0(path,i)))
}

testNoQT_list <- list()
for(i in testNoQT_files){
  tmp_name = stringr::str_split_fixed(string = i, ".fph.", n = 2)[2]
  testNoQT_list[[tmp_name]] = as.matrix(read.table(paste0(path,i)))
}

test1_list <- list()
for(i in test1_files){
  tmp_name = stringr::str_split_fixed(string = i, ".fph.", n = 2)[2]
  test1_list[[tmp_name]] = as.matrix(read.table(paste0(path,i)))
}

test2_list <- list()
for(i in test2_files){
  tmp_name = stringr::str_split_fixed(string = i, ".fph.", n = 2)[2]
  test2_list[[tmp_name]] = as.matrix(read.table(paste0(path,i)))
}

# Similarity between test (all snp test), test2 (individual snp test)
for(i in names(test_list)){
  print(i)
  print(all.equal(test_list[[i]],test2_list[[i]]))
}

# Difference in the pval file, first one lists only the strongest SNP.
# test2 lists all SNP (tests individually). check for similarity between
# the inidividual SNP record found in both.
all.equal(test_list$pval.txt[,1]
          , test2_list$pval.txt[,which(test2_list$pval.txt[1, ] == "chr17.10161485")])

lapply(test_list, dim)
lapply(testNoQT_list, dim)
lapply(test1_list, dim)
lapply(test2_list, dim)

# Test1 and test2 should be very similar. Test1 is just test2 without permutations, meaning
# that the presence of the p-val file should be the only difference.
for(i in names(test1_list)){
  print(i)
  print(all.equal(test1_list[[i]],test2_list[[i]]))
}

# noQT is test1, but on the noQT file. These might be quite different as one uses quant-t/formed
# WCs, the other doesn't.
for(i in names(test1_list)){
  print(i)
  print(all.equal(test1_list[[i]],testNoQT_list[[i]]))
}
# as expected.
