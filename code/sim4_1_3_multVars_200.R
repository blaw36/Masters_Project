
# Computer 1 --------------------------------------------------------------

for(len in proc1_lengths){
  print(len)
  # tmp <- list()
  tmp <- vector("list",length(ef_mult))
  n <- 1
  for(ef_size in ef_mult){
    print(ef_size)
    set.seed(1)
    tmp[[n]] <- run_sim4(
      sequencing_sums = effect_size_and_data$seq_sum
      , num_indivs = 70
      , num_bases = 1024
      , effect_length = len
      , effect_interval = get(paste0("effect_interval_",len))
      , effect_size_data = effect_size_and_data$effect_size
      , use_qt_data = FALSE
      , over_disp_param = od
      , Wmat_1024 = effect_size_and_data$Wmat_1024
      , W2mat_1024 = effect_size_and_data$W2mat_1024
      , library.read.depth = effect_size_and_data$library.read.depth
      , Covariates = effect_size_and_data$Covariates
      , effect_multiple = ef_size
      , trials_multiple = 10
      , number_sims = num_sims
      , verbose = F
      , rMarkdownMode = F
      , outputAlias = proc1_label
    )
    n <- n+1
  }
  assign(paste0("l",len,"_res"),tmp)
}

saveRDS(l8_res, "data/20190913_l8_grid_200sims.RDS",compress = T)
saveRDS(l16_res,"data/20190913_l16_grid_200sims.RDS",compress = T)
saveRDS(l24_res,"data/20190913_l24_grid_200sims.RDS",compress = T)
saveRDS(l32_res,"data/20190913_l32_grid_200sims.RDS",compress = T)


# Computer 2 --------------------------------------------------------------

for(len in proc2_lengths){
  print(len)
  # tmp <- list()
  tmp <- vector("list",length(ef_mult))
  n <- 1
  for(ef_size in ef_mult){
    print(ef_size)
    set.seed(1)
    tmp[[n]] <- run_sim4(
      sequencing_sums = effect_size_and_data$seq_sum
      , num_indivs = 70
      , num_bases = 1024
      , effect_length = len
      , effect_interval = get(paste0("effect_interval_",len))
      , effect_size_data = effect_size_and_data$effect_size
      , use_qt_data = FALSE
      , over_disp_param = od
      , Wmat_1024 = effect_size_and_data$Wmat_1024
      , W2mat_1024 = effect_size_and_data$W2mat_1024
      , library.read.depth = effect_size_and_data$library.read.depth
      , Covariates = effect_size_and_data$Covariates
      , effect_multiple = ef_size
      , trials_multiple = 10
      , number_sims = num_sims
      , verbose = F
      , rMarkdownMode = F
      , outputAlias = proc2_label
    )
    n <- n+1
  }
  assign(paste0("l",len,"_res"),tmp)
}

saveRDS(l40_res,"data/20190913_l40_grid_200sims.RDS",compress = T)
saveRDS(l48_res,"data/20190913_l48_grid_200sims.RDS",compress = T)
saveRDS(l56_res,"data/20190913_l56_grid_200sims.RDS",compress = T)
saveRDS(l64_res,"data/20190913_l64_grid_200sims.RDS",compress = T)


# Analysis ----------------------------------------------------------------
ef_mult <- seq(7e7,1.4e8,by = 1e7)

 l8_res <-  readRDS("data/20190913_l8_grid_200sims.RDS")
l16_res <- readRDS("data/20190913_l16_grid_200sims.RDS")
l24_res <- readRDS("data/20190913_l24_grid_200sims.RDS")
l32_res <- readRDS("data/20190913_l32_grid_200sims.RDS")
l40_res <- readRDS("data/20190913_l40_grid_200sims.RDS")
l48_res <- readRDS("data/20190913_l48_grid_200sims.RDS")
l56_res <- readRDS("data/20190913_l56_grid_200sims.RDS")
l64_res <- readRDS("data/20190913_l64_grid_200sims.RDS")

 l8_p <- lapply( l8_res,waveqtl_diags,num_sims = 200)
l16_p <- lapply(l16_res,waveqtl_diags,num_sims = 200)
l24_p <- lapply(l24_res,waveqtl_diags,num_sims = 200)
l32_p <- lapply(l32_res,waveqtl_diags,num_sims = 200)
l40_p <- lapply(l40_res,waveqtl_diags,num_sims = 200)
l48_p <- lapply(l48_res,waveqtl_diags,num_sims = 200)
l56_p <- lapply(l56_res,waveqtl_diags,num_sims = 200)
l64_p <- lapply(l64_res,waveqtl_diags,num_sims = 200)

l8_results <- sapply(l8_p,function(x){
  return(c(attr(x$perf_nohmt,"y.values")[[1]]
           ,attr(x$perf_hmt,"y.values")[[1]]))})
l16_results <- sapply(l16_p,function(x){
  return(c(attr(x$perf_nohmt,"y.values")[[1]]
           ,attr(x$perf_hmt,"y.values")[[1]]))})
l24_results <- sapply(l24_p,function(x){
  return(c(attr(x$perf_nohmt,"y.values")[[1]]
           ,attr(x$perf_hmt,"y.values")[[1]]))})
l32_results <- sapply(l32_p,function(x){
  return(c(attr(x$perf_nohmt,"y.values")[[1]]
           ,attr(x$perf_hmt,"y.values")[[1]]))})
l40_results <- sapply(l40_p,function(x){
  return(c(attr(x$perf_nohmt,"y.values")[[1]]
           ,attr(x$perf_hmt,"y.values")[[1]]))})
l48_results <- sapply(l48_p,function(x){
  return(c(attr(x$perf_nohmt,"y.values")[[1]]
           ,attr(x$perf_hmt,"y.values")[[1]]))})
l56_results <- sapply(l56_p,function(x){
  return(c(attr(x$perf_nohmt,"y.values")[[1]]
           ,attr(x$perf_hmt,"y.values")[[1]]))})
l64_results <- sapply(l64_p,function(x){
  return(c(attr(x$perf_nohmt,"y.values")[[1]]
           ,attr(x$perf_hmt,"y.values")[[1]]))})

a <- cbind(t(l8_results),ef_mult,rep(8,8))
b <- cbind(t(l16_results),ef_mult,rep(16,8))
c <- cbind(t(l24_results),ef_mult,rep(24,8))
d <- cbind(t(l32_results),ef_mult,rep(32,8))
f <- cbind(t(l40_results),ef_mult,rep(40,8))
g <- cbind(t(l48_results),ef_mult,rep(48,8))
h <- cbind(t(l56_results),ef_mult,rep(56,8))
i <- cbind(t(l64_results),ef_mult,rep(64,8))

e <- as.data.frame(rbind(rbind(rbind(a,b),c),d))
names(e) <- c("noHMT","HMT","efMult","efLength")
e2 <- as.data.frame(rbind(rbind(rbind(f,g),h),i))
names(e2) <- c("noHMT","HMT","efMult","efLength")

library(ggplot2)
library(reshape2)

ggplot(melt(e,id.vars=c("efMult","efLength"))) +
  geom_line(aes(x = efMult, y = value, colour = factor(variable))) +
  geom_point(aes(x = efMult, y = value, colour = factor(variable))) +
  facet_grid(efLength ~ .) +
  scale_x_continuous(breaks = ef_mult, labels = ef_mult)
ggplot(melt(e2,id.vars=c("efMult","efLength"))) +
  geom_line(aes(x = efMult, y = value, colour = factor(variable))) +
  geom_point(aes(x = efMult, y = value, colour = factor(variable))) +
  facet_grid(efLength ~ .) +
  scale_x_continuous(breaks = ef_mult, labels = ef_mult)

e3 <- rbind(e,e2)
ggplot(melt(e3,id.vars=c("efMult","efLength"))) +
  geom_line(aes(x = efMult, y = value, colour = factor(variable))) +
  geom_point(aes(x = efMult, y = value, colour = factor(variable))) +
  facet_grid(efLength ~ .) +
  scale_x_continuous(breaks = ef_mult, labels = ef_mult)

