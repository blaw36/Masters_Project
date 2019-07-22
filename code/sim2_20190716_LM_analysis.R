# Linear modelling ----
y_mtx <- read.table("~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_16wc.txt")
g_seq <- as.numeric(
  read.table("~/Cpp/WaveQTL_HMT/test/dsQTL/DeepdiveExamples/toy_eg1_16wc.cis.geno")[,4:73])

lm_list <- list()
for(i in 1:ncol(y_mtx)){
  tmp <- data.frame(y = y_mtx[[i]]
                    ,g = g_seq)
  lm_list[[i]] <- lm(y ~ 0 + g, data = tmp)
}
