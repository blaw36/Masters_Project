## EM algorithm with HMT test script.
# This is an R script version of the tester I made in excel.

# Input: first 7 logBFs of the dataset
# Output: set of parameter estimates from the HMT EM algorithm implementation
rm(list = ls()); gc();
library(data.table)

# # "chr17.10159002"
# logBFs = c(0.350769
#            ,-2.62E-05
#            ,-0.0155632
#            ,-0.0895906
#            ,0.147018
#            ,-0.0964465
#            ,-0.0330294
# )

# "chr17.10161485"
logBFs = c(-0.0134792
           ,9.7995
           ,5.45466
           ,-0.334555
           ,15.2401
           ,8.63809
           ,0.20578
)


groups = c(1, 2, 2, rep(3,4))
iterations_max = 0 #1000
conv_tol = 0.0005
diff <- Inf

# Helper functions with trees ---------------------------------------------

get_parent_indices <- function(indx, tree_root_indx = 1){
  return_indices <- indx
  return_indices <- (return_indices - tree_root_indx + 1) %/% 2

  # root of tree doesn't return any index
  return_indices[which(indx == tree_root_indx)] <- NA_integer_
  return(return_indices)
}

# Returns matrix: ij-th element is ith child of jth index
get_child_indices <- function(indx, tree_root_indx = 1){
  return_indices <- indx
  return_indices <- c((return_indices - tree_root_indx + 1) * 2
                      , ((return_indices - tree_root_indx + 1) * 2) + 1)
  return_indices <- matrix(return_indices,ncol = length(indx), byrow = T)
  return(return_indices)
}

tree_levels <- function(elements){
  return(seq(1:ceiling(log2(elements + 1))))
}

## sumlog function to prevent underflow
## will return the value of log(sum(exp(a + b)))
# sumlog_waveQTL <- function(a,b){
#   if(a > b){
#     res = a + log(1 + exp(b - a))
#   }else{
#     res = b + log(1 + exp(a - b))
#   }
#   return(res)
# }
## This one takes in vectors a and b
sumlog_waveQTL_gen <- function(a,b){
  mtx <- matrix(c(a,b), ncol = 2)
  sorted_mtx <- apply(mtx, MARGIN = 1, sort)
  t_sorted_mtx <- t(sorted_mtx)
  res <- t_sorted_mtx[,2] + log(1 + exp(t_sorted_mtx[,1] - t_sorted_mtx[,2]))
  return(res)
}

# Convert logBFs to BFs ---------------------------------------------------
bfs = exp(logBFs)

n = 0
levels = tree_levels(length(bfs))

b_i <- data.table(`1` = rep(NA_real_, length(bfs))
                  ,`0` = rep(NA_real_, length(bfs)))
b_i_pi <- data.table(`1` = rep(NA_real_, length(bfs))
                     ,`0` = rep(NA_real_, length(bfs)))
b_pi_no_i <- data.table(`1` = rep(NA_real_, length(bfs))
                        ,`0` = rep(NA_real_, length(bfs)))
a_i <- data.table(`1` = rep(NA_real_, length(bfs))
                  ,`0` = rep(NA_real_, length(bfs)))
pp_i <- data.table(`1` = rep(NA_real_, length(bfs))
                   ,`0` = rep(NA_real_, length(bfs)))
pp_j_i <- data.table(`11` = rep(NA_real_, length(bfs))
                     ,`10` = rep(NA_real_, length(bfs))
                     ,`01` = rep(NA_real_, length(bfs))
                     ,`00` = rep(NA_real_, length(bfs)))
a <- data.table(`1` = rep(NA_real_, length(bfs))
                ,`0` = rep(NA_real_, length(bfs)))
b <- data.table(`11` = rep(NA_real_, length(bfs))
                ,`10` = rep(NA_real_, length(bfs))
                ,`01` = rep(NA_real_, length(bfs))
                ,`00` = rep(NA_real_, length(bfs)))

# Parameters
eps <- data.table(`11` = c(NA, rep(0.5, length(bfs) - 1))
                  ,`10` = c(NA, rep(0.5, length(bfs) - 1)))

pi <- 0.5


# Initial E-step ----------------------------------------------------------

cat(paste0("Initial E-step \n"))

# Up step -----------------------------------------------------------------
for(i in rev(levels[-1])){
  elements_to_it <- which(groups == i)

  # Last level only
  if(i == rev(levels[-1])[1]){
    b_i[elements_to_it, "1" := bfs[elements_to_it]]
    b_i[elements_to_it, "0" := 1]
  }

  b_i_pi[elements_to_it, "1" := eps[elements_to_it, `11`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `11`]) * b_i[elements_to_it, `0`]]
  b_i_pi[elements_to_it, "0" := eps[elements_to_it, `10`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `10`]) * b_i[elements_to_it, `0`]]

  all_parents <- get_parent_indices(elements_to_it)
  parents <- unique(all_parents)
  for(j in parents){
    children <- get_child_indices(j)[,1]
    b_i[j, "1" := bfs[j]*prod(b_i_pi[children,`1`])]
    b_i[j, "0" := 1*prod(b_i_pi[children,`0`])]
  }

  b_pi_no_i[elements_to_it, "1" := b_i[all_parents,`1`]/b_i_pi[elements_to_it,`1`]]
  b_pi_no_i[elements_to_it, "0" := b_i[all_parents,`0`]/b_i_pi[elements_to_it,`0`]]

}

# Down step ---------------------------------------------------------------
for(i in levels){
  elements_to_it <- which(groups == i)
  if(i == levels[1]){
    a_i[elements_to_it, c("1","0") := .(pi, 1 - pi)]
  }else{
    all_parents <- get_parent_indices(elements_to_it)
    a_i[elements_to_it, `1` := (eps[elements_to_it, `11`] * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
          (eps[elements_to_it, `10`] * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
    a_i[elements_to_it, `0` := ((1 - eps[elements_to_it, `11`]) * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
          ((1 - eps[elements_to_it, `10`]) * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
  }
}

### DEBUGGING ###
cat(paste0("Number of NaNs in b_i: ",sum(is.nan(c(unlist(b_i))))))
cat(paste0("Max b_i: ",max(c(unlist(b_i)),na.rm = T)))
cat(paste0("Min b_i: ",min(c(unlist(b_i)),na.rm = T)))

cat(paste0("Number of NaNs in b_i_pi: ",sum(is.nan(c(unlist(b_i_pi))))))
cat(paste0("Max b_i_pi: ",max(c(unlist(b_i_pi)),na.rm = T)))
cat(paste0("Min b_i_pi: ",min(c(unlist(b_i_pi)),na.rm = T)))

cat(paste0("Number of NaNs in b_pi_no_i: ",sum(is.nan(c(unlist(b_pi_no_i))))))
cat(paste0("Max b_pi_no_i: ",max(c(unlist(b_pi_no_i)),na.rm = T)))
cat(paste0("Min b_pi_no_i: ",min(c(unlist(b_pi_no_i)),na.rm = T)))

cat(paste0("Number of NaNs in a_i: ",sum(is.nan(c(unlist(a_i))))))
cat(paste0("Max a_i: ",max(c(unlist(a_i)),na.rm = T)))
cat(paste0("Min a_i: ",min(c(unlist(a_i)),na.rm = T)))
### DEBUGGING END ###

# Posteriors --------------------------------------------------------------
denoms <- (b_i[,`1`]*a_i[,`1`]) + (b_i[,`0`]*a_i[,`0`])

pp_i[,`1` := b_i[,`1`]*a_i[,`1`]/denoms]
pp_i[,`0` := b_i[,`0`]*a_i[,`0`]/denoms]

# checks
apply(pp_i,1,sum)

all_parents <- get_parent_indices(seq(1:length(bfs)))
pp_j_i[,`11` := eps[, `11`] * b_i[,`1`] * b_pi_no_i[,`1`] * a_i[all_parents,`1`] / denoms]
pp_j_i[,`10` := eps[, `10`] * b_i[,`1`] * b_pi_no_i[,`0`] * a_i[all_parents,`0`] / denoms]
pp_j_i[,`01` := (1-eps[, `11`]) * b_i[,`0`] * b_pi_no_i[,`1`] * a_i[all_parents,`1`] / denoms]
pp_j_i[,`00` := (1-eps[, `10`]) * b_i[,`0`] * b_pi_no_i[,`0`] * a_i[all_parents,`0`] / denoms]

# checks
apply(pp_j_i,1,sum)

new_logL <- log(b_i[1,`1`]*a_i[1,`1`] + b_i[1,`0`]*a_i[1,`0`])
# cat(paste0("New LogL: ", new_logL,"\n"))

# Loop --------------------------------------------------------------------

while(n <= iterations_max && abs(diff) > conv_tol){

  old_logL <- new_logL
  cat(paste0("Iteration: ", n,"\n"))
  cat(paste0("Old LogL: ", old_logL,"\n"))

  # M-step ------------------------------------------------------------------
  # Updates -----------------------------------------------------------------

  a[,`1` := pp_i[,`1`]]
  a[,`0` := pp_i[,`0`]]

  # checks
  apply(a,1,sum)

  b[,`11` := pp_j_i[,`11`]/a[all_parents,`1`]]
  b[,`10` := pp_j_i[,`10`]/a[all_parents,`0`]]
  b[,`01` := pp_j_i[,`01`]/a[all_parents,`1`]]
  b[,`00` := pp_j_i[,`00`]/a[all_parents,`0`]]

  # checks
  apply(b[,c(1,3)],1,sum)
  apply(b[,c(2,4)],1,sum)

  # Update parameters -------------------------------------------------------
  pi <- a[1, `1`]
  eps <- b[,.(`11`,`10`)]

  # E-step ------------------------------------------------------------------
  # Up step -----------------------------------------------------------------
  for(i in rev(levels[-1])){
    elements_to_it <- which(groups == i)

    # Last level only
    if(i == rev(levels[-1])[1]){
      b_i[elements_to_it, "1" := bfs[elements_to_it]]
      b_i[elements_to_it, "0" := 1]
    }

    b_i_pi[elements_to_it, "1" := eps[elements_to_it, `11`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `11`]) * b_i[elements_to_it, `0`]]
    b_i_pi[elements_to_it, "0" := eps[elements_to_it, `10`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `10`]) * b_i[elements_to_it, `0`]]

    all_parents <- get_parent_indices(elements_to_it)
    parents <- unique(all_parents)
    for(j in parents){
      children <- get_child_indices(j)[,1]
      b_i[j, "1" := bfs[j]*prod(b_i_pi[children,`1`])]
      b_i[j, "0" := 1*prod(b_i_pi[children,`0`])]
    }

    b_pi_no_i[elements_to_it, "1" := b_i[all_parents,`1`]/b_i_pi[elements_to_it,`1`]]
    b_pi_no_i[elements_to_it, "0" := b_i[all_parents,`0`]/b_i_pi[elements_to_it,`0`]]

  }

  # Down step ---------------------------------------------------------------
  for(i in levels){
    elements_to_it <- which(groups == i)
    if(i == levels[1]){
      a_i[elements_to_it, c("1","0") := .(pi, 1 - pi)]
    }else{
      all_parents <- get_parent_indices(elements_to_it)
      a_i[elements_to_it, `1` := (eps[elements_to_it, `11`] * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
            (eps[elements_to_it, `10`] * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
      a_i[elements_to_it, `0` := ((1 - eps[elements_to_it, `11`]) * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
            ((1 - eps[elements_to_it, `10`]) * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
    }
  }


  ### DEBUGGING ###
  cat(paste0("Number of NaNs in b_i: ",sum(is.nan(c(unlist(b_i)))),"\n"))
  cat(paste0("Max b_i: ",max(c(unlist(b_i)),na.rm = T),"\n"))
  cat(paste0("Min b_i: ",min(c(unlist(b_i)),na.rm = T),"\n"))

  cat(paste0("Number of NaNs in b_i_pi: ",sum(is.nan(c(unlist(b_i_pi)))),"\n"))
  cat(paste0("Max b_i_pi: ",max(c(unlist(b_i_pi)),na.rm = T),"\n"))
  cat(paste0("Min b_i_pi: ",min(c(unlist(b_i_pi)),na.rm = T),"\n"))

  cat(paste0("Number of NaNs in b_pi_no_i: ",sum(is.nan(c(unlist(b_pi_no_i)))),"\n"))
  cat(paste0("Max b_pi_no_i: ",max(c(unlist(b_pi_no_i)),na.rm = T),"\n"))
  cat(paste0("Min b_pi_no_i: ",min(c(unlist(b_pi_no_i)),na.rm = T),"\n"))

  cat(paste0("Number of NaNs in a_i: ",sum(is.nan(c(unlist(a_i)))),"\n"))
  cat(paste0("Max a_i: ",max(c(unlist(a_i)),na.rm = T),"\n"))
  cat(paste0("Min a_i: ",min(c(unlist(a_i)),na.rm = T),"\n"))
  ### DEBUGGING END ###

  # Posteriors --------------------------------------------------------------
  denoms <- (b_i[,`1`]*a_i[,`1`]) + (b_i[,`0`]*a_i[,`0`])

  pp_i[,`1` := b_i[,`1`]*a_i[,`1`]/denoms]
  pp_i[,`0` := b_i[,`0`]*a_i[,`0`]/denoms]

  # checks
  apply(pp_i,1,sum)

  all_parents <- get_parent_indices(seq(1:length(bfs)))
  pp_j_i[,`11` := eps[, `11`] * b_i[,`1`] * b_pi_no_i[,`1`] * a_i[all_parents,`1`] / denoms]
  pp_j_i[,`10` := eps[, `10`] * b_i[,`1`] * b_pi_no_i[,`0`] * a_i[all_parents,`0`] / denoms]
  pp_j_i[,`01` := (1-eps[, `11`]) * b_i[,`0`] * b_pi_no_i[,`1`] * a_i[all_parents,`1`] / denoms]
  pp_j_i[,`00` := (1-eps[, `10`]) * b_i[,`0`] * b_pi_no_i[,`0`] * a_i[all_parents,`0`] / denoms]

  # checks
  apply(pp_j_i,1,sum)


  # LogL --------------------------------------------------------------------

  # old_logL <- new_logL
  new_logL <- log(b_i[1,`1`]*a_i[1,`1`] + b_i[1,`0`]*a_i[1,`0`])

  diff <- new_logL - old_logL

  n <- n + 1

}

