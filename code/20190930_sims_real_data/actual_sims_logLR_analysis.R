t <- readRDS("../../2019 Sem 2/20191003_own_sims/20191005_newLen-1-2.4_results.RDS")
null <- t[1:578]
null_hmt <- t[(578+1):(2*578)]
alt <- t[((2*578)+1):(3*578)]
alt_hmt <- t[((3*578)+1):(4*578)]

res_list <- list(
  null_waveqtl_lhood = t[1:578]
  ,null_waveqtl_hmt_lhood = t[(578+1):(2*578)]
  ,alt_waveqtl_lhood = t[((2*578)+1):(3*578)]
  ,alt_waveqtl_hmt_lhood = t[((3*578)+1):(4*578)]
)

library(ROCR)
extract_roc_plot_data <- function(sims_list,num_sims){
  preds_nohmt <- matrix(c(sims_list$null_waveqtl_lhood
                          ,sims_list$alt_waveqtl_lhood
                          ,rep(0,num_sims),rep(1,num_sims))
                        ,nrow = 2*num_sims,ncol = 2,byrow = F)
  preds_hmt <- matrix(c(sims_list$null_waveqtl_hmt_lhood
                        ,sims_list$alt_waveqtl_hmt_lhood
                        ,rep(0,num_sims),rep(1,num_sims))
                      ,nrow = 2*num_sims,ncol = 2,byrow = F)

  perf_nohmt <- performance(prediction(preds_nohmt[,1],preds_nohmt[,2]),measure = "auc")
  perf_hmt <- performance(prediction(preds_hmt[,1],preds_hmt[,2]),measure = "auc")

  n = dim(preds_hmt)[1]

  O1 = order(preds_nohmt[,1], decreasing =TRUE)
  O2 = order(preds_hmt[,1], decreasing =TRUE)

  C1  = c(0,cumsum(preds_nohmt[O1,2])) / sum(preds_nohmt[,2]) #True positives proportions or sensitivity
  C2  = c(0,cumsum(preds_hmt[O2,2])) / sum(preds_hmt[,2])

  FP1 = c(0,cumsum(1-preds_nohmt[O1,2])) / (n-sum(preds_nohmt[,2])) #false positives proportions based on Model 1.
  FP2 = c(0,cumsum(1-preds_hmt[O2,2])) / (n-sum(preds_hmt[,2]))

  return(list(
    FP1 = FP1
    ,C1 = C1
    ,FP2 = FP2
    ,C2 = C2
    ,perf_nohmt = unlist(attr(perf_nohmt,"y.values"))
    ,perf_hmt = unlist(attr(perf_hmt,"y.values"))
  ))
}

a <- extract_roc_plot_data(res_list,num_sims = 578)

plot(a$FP1,a$C1, type="l", lwd=1.2, col="black", lty = 3, xlab="False Positive rate", ylab="True Positive Rate"
     , main = paste0("578 simulations"))
lines(a$FP2,a$C2, lty = 1, lwd=1.2, col = "black")
legend("bottomright"
       ,legend = c(paste0("noHMT: ",round(a$perf_nohmt,4))
                   ,paste0("HMT: ",round(a$perf_hmt,4)))
       ,lty = c(3,1)
       ,cex = 0.8
       ,bty = "n")
