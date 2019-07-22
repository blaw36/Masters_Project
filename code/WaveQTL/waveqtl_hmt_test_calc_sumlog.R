## EM algorithm with HMT test script.
# This is an R script version of the tester I made in excel.
# Now this has the sumlog trick to prevent under/overflow.
# Also, it tries to work with logBFs rather than BFs. Converting only
# back to a non-log scale when we deal with probabilities.

# Input: first 7 logBFs of the dataset
# Output: set of parameter estimates from the HMT EM algorithm implementation
rm(list = ls()); gc();
library(data.table)

# You can actually retrieve these values from the logLR files, and multiply by
# log(10) (current values are stored as log of base 10). For example,
# a <- as.matrix(read.table("~/Cpp/WaveQTL_HMT/test/dsQTL/output/all_snp_hmt.fph.logLR.txt"))
# log(10)*as.numeric(a[11,4:1026])

# "chr17.10159002"
# logBFs = c(0.350769
#            ,-2.62E-05
#            ,-0.0155632
#            ,-0.0895906
#            ,0.147018
#            ,-0.0964465
#            ,-0.0330294
# )
# # Omit the scaling coefficient from this script!
# a <- as.matrix(read.table("~/Cpp/WaveQTL_HMT/test/dsQTL/output/all_snp_hmt.fph.logLR.txt"))
# logBFs = log(10)*as.numeric(a[1,4:1026])

# # "chr17.10161485"
# logBFs = c(-0.0134792
#            ,9.7995
#            ,5.45466
#            ,-0.334555
#            ,15.2401
#            ,8.63809
#            ,0.20578
# )

# a <- as.matrix(read.table("~/Cpp/WaveQTL_HMT/test/dsQTL/output/sim_07_03_bad.fph.logLR.txt"))
# logBFs = log(10)*as.numeric(a[1,4:(dim(a)[2])])

# Toy example ----
a <- as.matrix(read.table("~/Cpp/WaveQTL_HMT/test/dsQTL/output/toy_eg1_16wc_noTie_noHmt.fph.logLR.txt"))
logBFs = log(10)*as.numeric(a[1,4:(dim(a)[2])])

# # LogBF alterations
logBFs[between(logBFs,-2,3)] <- 0
logBFs[!between(logBFs,-2,3)] <- 1
# logBFs[between(logBFs,-2,3)] <- -1
# logBFs[!between(logBFs,-2,3)] <- 1

# # Tying alterations
# Tree level tying
tying_groups = c(1,2,4,8)
# Individual tying
# tying_groups = 1:15

# Toy example END ----

groups = floor(log2((1:length(logBFs))))+1

# starting group indicators,
# tying_groups = NULL # for no tying
# tying_groups = c(1,2,4,8,16,32,64,128,256,512)

# tying_groups = c(1,2,64)
# tying_groups = c(1,2)
# tying_groups = c(1,2,4)
# Note this should be c(1,2,3,5) in the c++ code to account for
# scaling coefficient @ 1 always being in its own group.


iterations_max = 1000 #1
conv_tol = 0.0005
diff <- Inf
sumlog_use <- T
init_pi = 0.5
init_eps_11 = 0.5
init_eps_10 = 0.5

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
## This one takes in vectors a and b
sumlog_waveQTL_gen <- function(a,b){
  mtx <- matrix(c(a,b), ncol = 2)
  sorted_mtx <- apply(mtx, MARGIN = 1, sort)
  t_sorted_mtx <- t(sorted_mtx)
  res <- t_sorted_mtx[,2] + log(1 + exp(t_sorted_mtx[,1] - t_sorted_mtx[,2]))
  return(res)
}

# Convert logBFs to BFs ---------------------------------------------------

if(!sumlog_use){
  bfs = exp(logBFs)
}else{
  #### LOG ADD ####
  bfs = logBFs
  ## END LOG ADD ##
}

if(is.null(tying_groups)){
  tying_groups = c(1:length(groups))
}

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
if(!sumlog_use){
  eps <- data.table(`11` = c(NA, rep(init_eps_11, length(bfs) - 1))
                    ,`10` = c(NA, rep(init_eps_10, length(bfs) - 1)))

  pi <- init_pi
}else{
  #### LOG ADD ####
  # eps <- data.table(`11` = c(NA, rep(log(init_eps_11), length(bfs) - 1))
  #                   ,`01` = c(NA, rep(log(1-init_eps_11), length(bfs) - 1))
  #                   ,`10` = c(NA, rep(log(init_eps_10), length(bfs) - 1))
  #                   ,`00` = c(NA, rep(log(1-init_eps_10), length(bfs) - 1)))
  eps <- data.table(`11` = c(NA, rep(log(0.5),2), rep(log(0.1), 4), rep(log(0.75), 8))
                    ,`01` = c(NA, rep(log(0.5),2), rep(log(0.9), 4), rep(log(0.25), 8))
                    ,`10` = c(NA, rep(log(0.5),2), rep(log(0.9), 4), rep(log(0.5), 8))
                    ,`00` = c(NA, rep(log(0.5),2), rep(log(0.1), 4), rep(log(0.5), 8)))
  # With bottom level tarnsition to 0
  # eps <- data.table(`11` = c(NA, rep(log(init_eps_11), 62), rep(log(0.000000001),64))
  #                   ,`01` = c(NA, rep(log(1-init_eps_11), 62), rep(log(1-0.000000001),64))
  #                   ,`10` = c(NA, rep(log(init_eps_10), 62), rep(log(0.000000001),64))
  #                   ,`00` = c(NA, rep(log(1-init_eps_10), 62), rep(log(1-0.000000001),64)))
  pi <- log(init_pi)
  pi_1_minus <- log(init_pi)
  ## END LOG ADD ##
}

# Initial E-step ----------------------------------------------------------

cat(paste0("Initial E-step \n"))

# Up step -----------------------------------------------------------------
for(i in rev(levels[-1])){
  elements_to_it <- which(groups == i)

  # Last level only
  if(i == rev(levels[-1])[1]){
    b_i[elements_to_it, "1" := bfs[elements_to_it]]
    if(!sumlog_use){
      b_i[elements_to_it, "0" := 1]
    }else{
      #### LOG ADD ####
      b_i[elements_to_it, "0" := log(1)]
      ## END LOG ADD ##
    }
  }

  if(!sumlog_use){
    b_i_pi[elements_to_it, "1" := eps[elements_to_it, `11`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `11`]) * b_i[elements_to_it, `0`]]
    b_i_pi[elements_to_it, "0" := eps[elements_to_it, `10`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `10`]) * b_i[elements_to_it, `0`]]
  }else{
    #### LOG ADD ####
    # log(b_i_pi) = log(eps1*b_i_1 + eps0*b_i_0)
    # = sumlog(log(eps1*b_i_1) + log(eps0*b_i_0))
    # = sumlog(log(eps1) + b_i_1, log(eps0) + b_i_0) as eps, b_i already logged.
    b_i_pi[elements_to_it, "1" := sumlog_waveQTL_gen(eps[elements_to_it, `11`] + b_i[elements_to_it, `1`]
                                                     ,eps[elements_to_it, `01`] + b_i[elements_to_it, `0`])]
    b_i_pi[elements_to_it, "0" := sumlog_waveQTL_gen(eps[elements_to_it, `10`] + b_i[elements_to_it, `1`]
                                                     ,eps[elements_to_it, `00`] + b_i[elements_to_it, `0`])]
    ## END LOG ADD ##
  }



  all_parents <- get_parent_indices(elements_to_it)
  parents <- unique(all_parents)
  for(j in parents){
    children <- get_child_indices(j)[,1]
    if(!sumlog_use){
      b_i[j, "1" := bfs[j]*prod(b_i_pi[children,`1`])]
      b_i[j, "0" := 1*prod(b_i_pi[children,`0`])]
    }else{
      #### LOG ADD ####
      b_i[j, "1" := bfs[j] + sum(b_i_pi[children,`1`])]
      b_i[j, "0" := sum(b_i_pi[children,`0`])]
      ## END LOG ADD ##
    }

  }

  if(!sumlog_use){
    b_pi_no_i[elements_to_it, "1" := b_i[all_parents,`1`]/b_i_pi[elements_to_it,`1`]]
    b_pi_no_i[elements_to_it, "0" := b_i[all_parents,`0`]/b_i_pi[elements_to_it,`0`]]
  }else{
    #### LOG ADD ####
    b_pi_no_i[elements_to_it, "1" := b_i[all_parents,`1`] - b_i_pi[elements_to_it,`1`]]
    b_pi_no_i[elements_to_it, "0" := b_i[all_parents,`0`] - b_i_pi[elements_to_it,`0`]]
    ## END LOG ADD ##
  }


}

# Down step ---------------------------------------------------------------
for(i in levels){
  elements_to_it <- which(groups == i)
  if(i == levels[1]){
    if(!sumlog_use){
      a_i[elements_to_it, c("1","0") := .(pi, 1 - pi)]
    }else{
      #### LOG ADD ####
      a_i[elements_to_it, c("1","0") := .(pi, pi_1_minus)]
      ## END LOG ADD ##
    }
  }else{
    all_parents <- get_parent_indices(elements_to_it)
    if(!sumlog_use){
      a_i[elements_to_it, `1` := (eps[elements_to_it, `11`] * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
            (eps[elements_to_it, `10`] * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
      a_i[elements_to_it, `0` := ((1 - eps[elements_to_it, `11`]) * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
            ((1 - eps[elements_to_it, `10`]) * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
    }else{
      #### LOG ADD ####
      a_i[elements_to_it, `1` := sumlog_waveQTL_gen(eps[elements_to_it, `11`] + a_i[all_parents,`1`] + b_pi_no_i[elements_to_it, `1`]
                                                    ,eps[elements_to_it, `10`] + a_i[all_parents,`0`] + b_pi_no_i[elements_to_it, `0`])]
      a_i[elements_to_it, `0` := sumlog_waveQTL_gen(eps[elements_to_it, `01`] + a_i[all_parents,`1`] + b_pi_no_i[elements_to_it, `1`]
                                                    ,eps[elements_to_it, `00`] + a_i[all_parents,`0`] + b_pi_no_i[elements_to_it, `0`])]
      ## END LOG ADD ##
    }
  }
}

# ### DEBUGGING ###
# cat(paste0("Number of NaNs in b_i: ",sum(is.nan(c(unlist(b_i))))))
# cat(paste0("Max b_i: ",max(c(unlist(b_i)),na.rm = T)))
# cat(paste0("Min b_i: ",min(c(unlist(b_i)),na.rm = T)))
#
# cat(paste0("Number of NaNs in b_i_pi: ",sum(is.nan(c(unlist(b_i_pi))))))
# cat(paste0("Max b_i_pi: ",max(c(unlist(b_i_pi)),na.rm = T)))
# cat(paste0("Min b_i_pi: ",min(c(unlist(b_i_pi)),na.rm = T)))
#
# cat(paste0("Number of NaNs in b_pi_no_i: ",sum(is.nan(c(unlist(b_pi_no_i))))))
# cat(paste0("Max b_pi_no_i: ",max(c(unlist(b_pi_no_i)),na.rm = T)))
# cat(paste0("Min b_pi_no_i: ",min(c(unlist(b_pi_no_i)),na.rm = T)))
#
# cat(paste0("Number of NaNs in a_i: ",sum(is.nan(c(unlist(a_i))))))
# cat(paste0("Max a_i: ",max(c(unlist(a_i)),na.rm = T)))
# cat(paste0("Min a_i: ",min(c(unlist(a_i)),na.rm = T)))
# ### DEBUGGING END ###

# Posteriors --------------------------------------------------------------

if(!sumlog_use){
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
}else{
  #### LOG ADD ####
  # Calculate these qtys in log for now.
  denoms <- sumlog_waveQTL_gen(b_i[,`1`] + a_i[,`1`]
                               , b_i[,`0`] + a_i[,`0`])

  pp_i[,`1` := b_i[,`1`] + a_i[,`1`] - denoms]
  pp_i[,`0` := b_i[,`0`] + a_i[,`0`] - denoms]

  # checks
  apply(exp(pp_i),1,sum)

  all_parents <- get_parent_indices(seq(1:length(bfs)))
  pp_j_i[,`11` := eps[, `11`] + b_i[,`1`] + b_pi_no_i[,`1`] + a_i[all_parents,`1`] - denoms]
  pp_j_i[,`10` := eps[, `10`] + b_i[,`1`] + b_pi_no_i[,`0`] + a_i[all_parents,`0`] - denoms]
  pp_j_i[,`01` := eps[, `01`] + b_i[,`0`] + b_pi_no_i[,`1`] + a_i[all_parents,`1`] - denoms]
  pp_j_i[,`00` := eps[, `00`] + b_i[,`0`] + b_pi_no_i[,`0`] + a_i[all_parents,`0`] - denoms]

  # checks
  apply(exp(pp_j_i),1,sum)

  new_logL <- sumlog_waveQTL_gen(b_i[1,`1`] + a_i[1,`1`]
                                 , b_i[1,`0`] + a_i[1,`0`])
  # cat(paste0("New LogL: ", new_logL,"\n"))
  ## END LOG ADD ##
}


# Loop --------------------------------------------------------------------

while(n < iterations_max && abs(diff) > conv_tol){

  old_logL <- new_logL
  cat(paste0("Iteration: ", n,"\n"))
  cat(paste0("Old LogL: ", old_logL,"\n"))

  # M-step ------------------------------------------------------------------
  # Updates -----------------------------------------------------------------

  a[,`1` := pp_i[,`1`]]
  a[,`0` := pp_i[,`0`]]

  # checks
  if(!sumlog_use){
    apply(a,1,sum)
  }else{
    #### LOG ADD ####
    apply(exp(a),1,sum)
    ## END LOG ADD ##
  }

  if(!sumlog_use){
    b[,`11` := pp_j_i[,`11`]/a[all_parents,`1`]]
    b[,`10` := pp_j_i[,`10`]/a[all_parents,`0`]]
    b[,`01` := pp_j_i[,`01`]/a[all_parents,`1`]]
    b[,`00` := pp_j_i[,`00`]/a[all_parents,`0`]]
  }else{
    #### LOG ADD ####
    b[,`11` := pp_j_i[,`11`] - a[all_parents,`1`]]
    b[,`10` := pp_j_i[,`10`] - a[all_parents,`0`]]
    b[,`01` := pp_j_i[,`01`] - a[all_parents,`1`]]
    b[,`00` := pp_j_i[,`00`] - a[all_parents,`0`]]
    ## END LOG ADD ##
  }



  # checks
  if(!sumlog_use){
    apply(b[,c(1,3)],1,sum)
    apply(b[,c(2,4)],1,sum)
  }else{
    #### LOG ADD ####
    apply(exp(b[,c(1,3)]),1,sum)
    apply(exp(b[,c(2,4)]),1,sum)
    ## END LOG ADD ##
  }

  # Update parameters -------------------------------------------------------

  # Set up tying groups
  num_tying_grps = length(tying_groups)

  tying_start = integer()
  tying_end = integer()
  num_pheno_tying_grp = integer()
  for(i in 1:(num_tying_grps - 1)){
    tying_start = c(tying_start,tying_groups[i])
    tying_end = c(tying_end,tying_groups[i + 1] - 1)
    num_pheno_tying_grp = c(num_pheno_tying_grp,(tying_groups[i + 1] - 1) - tying_groups[i] + 1)
  }
  tying_start[num_tying_grps] = tying_groups[num_tying_grps]
  tying_end[num_tying_grps] = length(groups)
  num_pheno_tying_grp[num_tying_grps] = tying_end[num_tying_grps] - tying_start[num_tying_grps] + 1


  if(!sumlog_use){
    # Tying pi
    # Pi only has 1 coefficient - the first tying group
    # pi <- pp_i[tying_start[1]:tying_end[1], mean(`1`)]
    pi <- pp_i[1, mean(`1`)]

    # Tying epsilon
    for(i in 1:num_tying_grps){
      # Need to exclude tree-root: undefined parameter for that element
      start_excl_head = max(2,tying_start[i])
      indices_to_sum <- start_excl_head:tying_end[i]
      parents_indices_to_sum <- get_parent_indices(indices_to_sum)
      eps[indices_to_sum,] <-
        pp_j_i[indices_to_sum , .(`11` = sum(`11`), `10` = sum(`10`), `01` = sum(`01`), `00` = sum(`00`))]/
        pp_i[parents_indices_to_sum , .(`11` = sum(`1`), `10` = sum(`0`), `01` = sum(`1`), `00` = sum(`0`))]
    }

  }else{
    # Tying pi
    # Pi only has 1 coefficient - the first tying group
    indices_to_sum <- tying_start[1]:tying_end[1]

    num_pheno_grp_1 <- length(indices_to_sum)
    # if(num_pheno_grp_1 == 1){
    #   pi <- pp_i[indices_to_sum, (`1`)]
    #   pi_1_minus <- pp_i[indices_to_sum, (`0`)]
    # }else{
    #   pi_1_values <- pp_i[indices_to_sum, (`1`)]
    #   pi_0_values <- pp_i[indices_to_sum, (`0`)]
    #
    #   pi <- pi_1_values[1]
    #   for(i in 1:(num_pheno_grp_1 - 1)){
    #     pi <- sumlog_waveQTL_gen(pi,pi_1_values[i+1])
    #   }
    #   pi <- pi - log(num_pheno_grp_1)
    #
    #   pi_1_minus <- pi_0_values[1]
    #   for(i in 1:(num_pheno_grp_1 - 1)){
    #     pi_1_minus <- sumlog_waveQTL_gen(pi_1_minus,pi_0_values[i+1])
    #   }
    #   pi_1_minus <- pi_1_minus - log(num_pheno_grp_1)
    #
    # }
    pi <- pp_i[1, (`1`)]
    pi_1_minus <- pp_i[1, (`0`)]


    # Tying epsilon - need to exclude tree-root: undefined parameter for that element
    for(i in 1:num_tying_grps){
      # Need to exclude tree-root: undefined parameter for that element
      start_excl_head = max(2,tying_start[i])
      indices_to_sum <- start_excl_head:tying_end[i]
      parents_indices_to_sum <- get_parent_indices(indices_to_sum)

      if(tying_end[i] < start_excl_head){
        next
      }else{
        num_pheno <- length(indices_to_sum)

        if(num_pheno == 1){
          eps[indices_to_sum,] <- b[indices_to_sum, .(`11`, `01`, `10`, `00`)]

        }else{
          pp_values <- pp_i[parents_indices_to_sum,.(`1`,`1`,`0`,`0`)]
          pp_joint_values <- pp_j_i[indices_to_sum, .(`11`,`01`,`10`,`00`)]

          for(j in 1:4){
            eps_logs <- pp_joint_values[1][[j]]
            eps_denom_logs <- pp_values[1][[j]]

            for(k in 1:(num_pheno - 1)){
              eps_logs <- sumlog_waveQTL_gen(eps_logs,pp_joint_values[k + 1][[j]])
              eps_denom_logs <- sumlog_waveQTL_gen(eps_denom_logs,pp_values[k + 1][[j]])
            }
            eps[indices_to_sum][[j]] <- eps_logs - eps_denom_logs
          }

        }
      }
    }

  }

  # E-step ------------------------------------------------------------------
  # Up step -----------------------------------------------------------------
  for(i in rev(levels[-1])){
    elements_to_it <- which(groups == i)

    # Last level only
    if(i == rev(levels[-1])[1]){
      b_i[elements_to_it, "1" := bfs[elements_to_it]]
      if(!sumlog_use){
        b_i[elements_to_it, "0" := 1]
      }else{
        #### LOG ADD ####
        b_i[elements_to_it, "0" := log(1)]
        ## END LOG ADD ##
      }
    }

    if(!sumlog_use){
      b_i_pi[elements_to_it, "1" := eps[elements_to_it, `11`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `11`]) * b_i[elements_to_it, `0`]]
      b_i_pi[elements_to_it, "0" := eps[elements_to_it, `10`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `10`]) * b_i[elements_to_it, `0`]]
    }else{
      #### LOG ADD ####
      # log(b_i_pi) = log(eps1*b_i_1 + eps0*b_i_0)
      # = sumlog(log(eps1*b_i_1) + log(eps0*b_i_0))
      # = sumlog(log(eps1) + b_i_1, log(eps0) + b_i_0) as b_i already logged.
      b_i_pi[elements_to_it, "1" := sumlog_waveQTL_gen(eps[elements_to_it, `11`] + b_i[elements_to_it, `1`]
                                                       ,eps[elements_to_it, `01`] + b_i[elements_to_it, `0`])]
      b_i_pi[elements_to_it, "0" := sumlog_waveQTL_gen(eps[elements_to_it, `10`] + b_i[elements_to_it, `1`]
                                                       ,eps[elements_to_it, `00`] + b_i[elements_to_it, `0`])]
      ## END LOG ADD ##
    }



    all_parents <- get_parent_indices(elements_to_it)
    parents <- unique(all_parents)
    for(j in parents){
      children <- get_child_indices(j)[,1]
      if(!sumlog_use){
        b_i[j, "1" := bfs[j]*prod(b_i_pi[children,`1`])]
        b_i[j, "0" := 1*prod(b_i_pi[children,`0`])]
      }else{
        #### LOG ADD ####
        b_i[j, "1" := bfs[j] + sum(b_i_pi[children,`1`])]
        b_i[j, "0" := sum(b_i_pi[children,`0`])]
        ## END LOG ADD ##
      }

    }

    if(!sumlog_use){
      b_pi_no_i[elements_to_it, "1" := b_i[all_parents,`1`]/b_i_pi[elements_to_it,`1`]]
      b_pi_no_i[elements_to_it, "0" := b_i[all_parents,`0`]/b_i_pi[elements_to_it,`0`]]
    }else{
      #### LOG ADD ####
      b_pi_no_i[elements_to_it, "1" := b_i[all_parents,`1`] - b_i_pi[elements_to_it,`1`]]
      b_pi_no_i[elements_to_it, "0" := b_i[all_parents,`0`] - b_i_pi[elements_to_it,`0`]]
      ## END LOG ADD ##
    }


  }

  # Down step ---------------------------------------------------------------
  for(i in levels){
    elements_to_it <- which(groups == i)
    if(i == levels[1]){
      if(!sumlog_use){
        a_i[elements_to_it, c("1","0") := .(pi, 1 - pi)]
      }else{
        #### LOG ADD ####
        a_i[elements_to_it, c("1","0") := .(pi, pi_1_minus)]
        ## END LOG ADD ##
      }
    }else{
      all_parents <- get_parent_indices(elements_to_it)
      if(!sumlog_use){
        a_i[elements_to_it, `1` := (eps[elements_to_it, `11`] * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
              (eps[elements_to_it, `10`] * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
        a_i[elements_to_it, `0` := ((1 - eps[elements_to_it, `11`]) * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
              ((1 - eps[elements_to_it, `10`]) * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
      }else{
        #### LOG ADD ####
        a_i[elements_to_it, `1` := sumlog_waveQTL_gen(eps[elements_to_it, `11`] + a_i[all_parents,`1`] + b_pi_no_i[elements_to_it, `1`]
                                                      ,eps[elements_to_it, `10`] + a_i[all_parents,`0`] + b_pi_no_i[elements_to_it, `0`])]
        a_i[elements_to_it, `0` := sumlog_waveQTL_gen(eps[elements_to_it, `01`] + a_i[all_parents,`1`] + b_pi_no_i[elements_to_it, `1`]
                                                      ,eps[elements_to_it, `00`] + a_i[all_parents,`0`] + b_pi_no_i[elements_to_it, `0`])]
        ## END LOG ADD ##
      }
    }
  }


  # ### DEBUGGING ###
  # cat(paste0("Number of NaNs in b_i: ",sum(is.nan(c(unlist(b_i)))),"\n"))
  # cat(paste0("Max b_i: ",max(c(unlist(b_i)),na.rm = T),"\n"))
  # cat(paste0("Min b_i: ",min(c(unlist(b_i)),na.rm = T),"\n"))
  #
  # cat(paste0("Number of NaNs in b_i_pi: ",sum(is.nan(c(unlist(b_i_pi)))),"\n"))
  # cat(paste0("Max b_i_pi: ",max(c(unlist(b_i_pi)),na.rm = T),"\n"))
  # cat(paste0("Min b_i_pi: ",min(c(unlist(b_i_pi)),na.rm = T),"\n"))
  #
  # cat(paste0("Number of NaNs in b_pi_no_i: ",sum(is.nan(c(unlist(b_pi_no_i)))),"\n"))
  # cat(paste0("Max b_pi_no_i: ",max(c(unlist(b_pi_no_i)),na.rm = T),"\n"))
  # cat(paste0("Min b_pi_no_i: ",min(c(unlist(b_pi_no_i)),na.rm = T),"\n"))
  #
  # cat(paste0("Number of NaNs in a_i: ",sum(is.nan(c(unlist(a_i)))),"\n"))
  # cat(paste0("Max a_i: ",max(c(unlist(a_i)),na.rm = T),"\n"))
  # cat(paste0("Min a_i: ",min(c(unlist(a_i)),na.rm = T),"\n"))
  # ### DEBUGGING END ###

  # Posteriors --------------------------------------------------------------

  if(!sumlog_use){
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
  }else{
    #### LOG ADD ####
    # Calculate these qtys in log for now.
    denoms <- sumlog_waveQTL_gen(b_i[,`1`] + a_i[,`1`]
                                 , b_i[,`0`] + a_i[,`0`])

    pp_i[,`1` := b_i[,`1`] + a_i[,`1`] - denoms]
    pp_i[,`0` := b_i[,`0`] + a_i[,`0`] - denoms]

    # checks
    apply(exp(pp_i),1,sum)

    all_parents <- get_parent_indices(seq(1:length(bfs)))
    pp_j_i[,`11` := eps[, `11`] + b_i[,`1`] + b_pi_no_i[,`1`] + a_i[all_parents,`1`] - denoms]
    pp_j_i[,`10` := eps[, `10`] + b_i[,`1`] + b_pi_no_i[,`0`] + a_i[all_parents,`0`] - denoms]
    pp_j_i[,`01` := eps[, `01`] + b_i[,`0`] + b_pi_no_i[,`1`] + a_i[all_parents,`1`] - denoms]
    pp_j_i[,`00` := eps[, `00`] + b_i[,`0`] + b_pi_no_i[,`0`] + a_i[all_parents,`0`] - denoms]

    # checks
    apply(exp(pp_j_i),1,sum)

    new_logL <- sumlog_waveQTL_gen(b_i[1,`1`] + a_i[1,`1`]
                                   , b_i[1,`0`] +a_i[1,`0`])
    # cat(paste0("New LogL: ", new_logL,"\n"))
    ## END LOG ADD ##
  }


  diff <- new_logL - old_logL

  n <- n + 1
print(exp(eps[tying_groups,]))
print(exp(pp_i))
print(cbind(round(exp(pp_i$`1`),5),logBFs))
}

cat(paste0("New LogL: ", new_logL,"\n"))
round(exp(eps[tying_groups,]),5)
round(exp(pp_i),5)
print(cbind(round(exp(pp_i$`1`),5),logBFs))
