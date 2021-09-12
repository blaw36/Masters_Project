## EM algorithm with HMT test script.
### This script allows INIT PARAMETER TO BE VARIED AND PLAYED WITH
### This v2 now adds functionality for the scaling coefficient built-in.
### The logBFs and groups inputs must reflect this.

wqtl_hmt_v2 <- function(
  logBFs,
  # groups,
  tying_groups = c(1,c(1,2,4,8,16,32,64,128,256,512)+1),
  iterations_max,
  conv_tol,
  sumlog_use,
  init_pi,
  init_eps_11,
  init_eps_10
){
  
  eps_vals <- list()
  logl_list <- list()
  
  inits_loop=1
  
  # init_try = 0.5
  
  diff <- Inf
  # init_pi = 0.5
  # # if(inits_loop <= 4){
  # if(inits_loop <= 3){
  #   init_eps_11 = init_try
  #   init_eps_10 = init_try
  # }else{
  #   init_eps_11 = init_try
  #   init_eps_10 = 1 - init_try
  # }
  
  logl_list[[inits_loop]] <- 0
  tmp_eps <- list()
  
  
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
    return(1:ceiling(log2(elements + 1)))
  }
  
  ## sumlog function to prevent underflow
  ## This one takes in vectors a and b
  sumlog_waveQTL_gen <- function(a,b){
    mtx <- matrix(c(a,b), ncol = 2)
    sorted_mtx <- apply(mtx, MARGIN = 1, sort)
    res <- sorted_mtx[2,] + log(1 + exp(sorted_mtx[1,] - sorted_mtx[2,]))
    # t_sorted_mtx <- t(sorted_mtx)
    # res <- t_sorted_mtx[,2] + log(1 + exp(t_sorted_mtx[,1] - t_sorted_mtx[,2]))
    return(res)
  }

  group_converter = function(tying_groups,logBFs){
    group_vector = rep(NA,length(logBFs))
    for(i in 1:(length(tying_groups)-1)){
      start <- tying_groups[i]
      end <- tying_groups[i+1]-1
      group_vector[start:end] <- i
    }
    group_vector[(end+1):length(group_vector)] <- i+1
    return(group_vector)
  }
  
  # Convert logBFs to BFs ---------------------------------------------------
  
  if(!sumlog_use){
    bfs = exp(logBFs)
  }else{
    #### LOG ADD ####
    bfs = logBFs
    ## END LOG ADD ##
  }
  
  # Split out BFs and tying_group data for SC and tree ------------------------

  # BFs
  bfs_sc = bfs[1]
  bfs = bfs[-1]

  # Tying group information
  sc_tying_groups = tying_groups[1]
  tying_groups = tying_groups[-1] - 1

  # if(is.null(tying_groups)){
  #   tying_groups = c(1:length(groups))
  # }
  groups = group_converter(tying_groups,bfs)
  levels = uniqueN(groups)
  
  n = 0
  # levels = tree_levels(length(bfs))
  
  # b_i <- data.table(`1` = rep(NA_real_, length(bfs))
  #                   ,`0` = rep(NA_real_, length(bfs)))
  # b_i_pi <- data.table(`1` = rep(NA_real_, length(bfs))
  #                      ,`0` = rep(NA_real_, length(bfs)))
  # b_pi_no_i <- data.table(`1` = rep(NA_real_, length(bfs))
  #                         ,`0` = rep(NA_real_, length(bfs)))
  # a_i <- data.table(`1` = rep(NA_real_, length(bfs))
  #                   ,`0` = rep(NA_real_, length(bfs)))
  # pp_i <- data.table(`1` = rep(NA_real_, length(bfs))
  #                    ,`0` = rep(NA_real_, length(bfs)))
  # pp_j_i <- data.table(`11` = rep(NA_real_, length(bfs))
  #                      ,`10` = rep(NA_real_, length(bfs))
  #                      ,`01` = rep(NA_real_, length(bfs))
  #                      ,`00` = rep(NA_real_, length(bfs)))
  # a <- data.table(`1` = rep(NA_real_, length(bfs))
  #                 ,`0` = rep(NA_real_, length(bfs)))
  # b <- data.table(`11` = rep(NA_real_, length(bfs))
  #                 ,`10` = rep(NA_real_, length(bfs))
  #                 ,`01` = rep(NA_real_, length(bfs))
  #                 ,`00` = rep(NA_real_, length(bfs)))
  
  b_i <- matrix(rep(NA_real_,length(bfs)*2),ncol = 2)
  b_i_pi <- matrix(rep(NA_real_,length(bfs)*2),ncol = 2)
  b_pi_no_i <- matrix(rep(NA_real_,length(bfs)*2),ncol = 2)
  a_i <- matrix(rep(NA_real_,length(bfs)*2),ncol = 2)
  pp_i <- matrix(rep(NA_real_,length(bfs)*2),ncol = 2)
  pp_j_i <- matrix(rep(NA_real_,length(bfs)*4),ncol = 4)
  a <- matrix(rep(NA_real_,length(bfs)*2),ncol = 2)
  b <- matrix(rep(NA_real_,length(bfs)*4),ncol = 4)
  
  # Parameters
  if(!sumlog_use){
    # eps <- data.table(`11` = c(NA, rep(init_eps_11, length(bfs) - 1))
    #                   ,`10` = c(NA, rep(init_eps_10, length(bfs) - 1)))
    eps <- matrix(c(
      c(NA, rep(init_eps_11, length(bfs) - 1))
      ,c(NA, rep(init_eps_10, length(bfs) - 1))
    ),ncol = 2)
    
    pi <- init_pi
  }else{
    #### LOG ADD ####
    # eps <- data.table(`11` = c(NA, rep(log(init_eps_11), length(bfs) - 1))
    #                   ,`01` = c(NA, rep(log(1-init_eps_11), length(bfs) - 1))
    #                   ,`10` = c(NA, rep(log(init_eps_10), length(bfs) - 1))
    #                   ,`00` = c(NA, rep(log(1-init_eps_10), length(bfs) - 1)))
    eps <- matrix(c(
      c(NA, rep(log(init_eps_11), length(bfs) - 1))
      ,c(NA, rep(log(1-init_eps_11), length(bfs) - 1))
      ,c(NA, rep(log(init_eps_10), length(bfs) - 1))
      ,c(NA, rep(log(1-init_eps_10), length(bfs) - 1))
    ),ncol = 4)
    # # With bottom level transition to 0
    # eps <- data.table(`11` = c(NA, rep(log(init_eps_11), 62), rep(log(0.000000001),64))
    #                   ,`01` = c(NA, rep(log(1-init_eps_11), 62), rep(log(1-0.000000001),64))
    #                   ,`10` = c(NA, rep(log(init_eps_10), 62), rep(log(0.000000001),64))
    #                   ,`00` = c(NA, rep(log(1-init_eps_10), 62), rep(log(1-0.000000001),64)))
    pi <- log(init_pi)
    pi_1_minus <- log(1-init_pi)
    ## END LOG ADD ##
  }
  
  # Initial E-step ----------------------------------------------------------
  
  cat(paste0("Initial E-step \n"))
  
  # Up step -----------------------------------------------------------------
  for(i in rev((1:levels)[-1])){
    elements_to_it <- which(groups == i)
    
    # Last level only
    if(i == rev((1:levels)[-1])[1]){
      # b_i[elements_to_it, "1" := bfs[elements_to_it]]
      b_i[elements_to_it,1] = bfs[elements_to_it]
      if(!sumlog_use){
        # b_i[elements_to_it, "0" := 1]
        b_i[elements_to_it,2] = 1
      }else{
        #### LOG ADD ####
        # b_i[elements_to_it, "0" := log(1)]
        b_i[elements_to_it,2] = log(1)
        ## END LOG ADD ##
      }
    }
    
    if(!sumlog_use){
      # b_i_pi[elements_to_it, "1" := eps[elements_to_it, `11`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `11`]) * b_i[elements_to_it, `0`]]
      # b_i_pi[elements_to_it, "0" := eps[elements_to_it, `10`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `10`]) * b_i[elements_to_it, `0`]]
      b_i_pi[elements_to_it,1] =
        eps[elements_to_it, 1] * b_i[elements_to_it, 1] + (1 - eps[elements_to_it, 1]) * b_i[elements_to_it, 2]
      b_i_pi[elements_to_it,2] =
        eps[elements_to_it, 3] * b_i[elements_to_it, 1] + (1 - eps[elements_to_it, 3]) * b_i[elements_to_it, 2]
    }else{
      #### LOG ADD ####
      # log(b_i_pi) = log(eps1*b_i_1 + eps0*b_i_0)
      # = sumlog(log(eps1*b_i_1) + log(eps0*b_i_0))
      # = sumlog(log(eps1) + b_i_1, log(eps0) + b_i_0) as eps, b_i already logged.
      # b_i_pi[elements_to_it, "1" := sumlog_waveQTL_gen(eps[elements_to_it, `11`] + b_i[elements_to_it, `1`]
      #                                                  ,eps[elements_to_it, `01`] + b_i[elements_to_it, `0`])]
      # b_i_pi[elements_to_it, "0" := sumlog_waveQTL_gen(eps[elements_to_it, `10`] + b_i[elements_to_it, `1`]
      #                                                  ,eps[elements_to_it, `00`] + b_i[elements_to_it, `0`])]
      b_i_pi[elements_to_it,1] = sumlog_waveQTL_gen(eps[elements_to_it, 1] + b_i[elements_to_it, 1]
                                                    ,eps[elements_to_it, 2] + b_i[elements_to_it, 2])
      b_i_pi[elements_to_it,2] = sumlog_waveQTL_gen(eps[elements_to_it, 3] + b_i[elements_to_it, 1]
                                                    ,eps[elements_to_it, 4] + b_i[elements_to_it, 2])
      ## END LOG ADD ##
    }
    
    
    
    all_parents <- get_parent_indices(elements_to_it)
    parents <- unique(all_parents)
    for(j in parents){
      children <- get_child_indices(j)[,1]
      if(!sumlog_use){
        # b_i[j, "1" := bfs[j]*prod(b_i_pi[children,`1`])]
        # b_i[j, "0" := 1*prod(b_i_pi[children,`0`])]
        b_i[j, 1] = bfs[j]*prod(b_i_pi[children,1])
        b_i[j, 2] = 1*prod(b_i_pi[children,2])
      }else{
        #### LOG ADD ####
        # b_i[j, "1" := bfs[j] + sum(b_i_pi[children,`1`])]
        # b_i[j, "0" := sum(b_i_pi[children,`0`])]
        b_i[j, 1] = bfs[j] + sum(b_i_pi[children,1])
        b_i[j, 2] = sum(b_i_pi[children,2])
        ## END LOG ADD ##
      }
      
    }
    
    if(!sumlog_use){
      # b_pi_no_i[elements_to_it, "1" := b_i[all_parents,`1`]/b_i_pi[elements_to_it,`1`]]
      # b_pi_no_i[elements_to_it, "0" := b_i[all_parents,`0`]/b_i_pi[elements_to_it,`0`]]
      b_pi_no_i[elements_to_it, 1] = b_i[all_parents,1]/b_i_pi[elements_to_it,1]
      b_pi_no_i[elements_to_it, 2] = b_i[all_parents,2]/b_i_pi[elements_to_it,2]
    }else{
      #### LOG ADD ####
      # b_pi_no_i[elements_to_it, "1" := b_i[all_parents,`1`] - b_i_pi[elements_to_it,`1`]]
      # b_pi_no_i[elements_to_it, "0" := b_i[all_parents,`0`] - b_i_pi[elements_to_it,`0`]]
      b_pi_no_i[elements_to_it, 1] = b_i[all_parents,1] - b_i_pi[elements_to_it,1]
      b_pi_no_i[elements_to_it, 2] = b_i[all_parents,2] - b_i_pi[elements_to_it,2]
      ## END LOG ADD ##
    }
    
    
  }
  
  # Down step ---------------------------------------------------------------
  for(i in 1:levels){
    elements_to_it <- which(groups == i)
    if(i == (1:levels)[1]){
      if(!sumlog_use){
        # a_i[elements_to_it, c("1","0") := .(pi, 1 - pi)]
        a_i[elements_to_it,1] = pi
        a_i[elements_to_it,2] = 1-pi
      }else{
        #### LOG ADD ####
        # a_i[elements_to_it, c("1","0") := .(pi, pi_1_minus)]
        a_i[elements_to_it, 1] = pi
        a_i[elements_to_it, 2] = pi_1_minus
        ## END LOG ADD ##
      }
    }else{
      all_parents <- get_parent_indices(elements_to_it)
      if(!sumlog_use){
        # a_i[elements_to_it, `1` := (eps[elements_to_it, `11`] * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
        #       (eps[elements_to_it, `10`] * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
        # a_i[elements_to_it, `0` := ((1 - eps[elements_to_it, `11`]) * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
        #       ((1 - eps[elements_to_it, `10`]) * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
        a_i[elements_to_it, 1] = (eps[elements_to_it, 1] * a_i[all_parents,1] * b_pi_no_i[elements_to_it, 1]) +
          (eps[elements_to_it, 3] * a_i[all_parents,2] * b_pi_no_i[elements_to_it, 2])
        a_i[elements_to_it, 2] = ((1 - eps[elements_to_it, 1]) * a_i[all_parents,1] * b_pi_no_i[elements_to_it, 1]) +
          ((1 - eps[elements_to_it, 3]) * a_i[all_parents,2] * b_pi_no_i[elements_to_it, 2])
      }else{
        #### LOG ADD ####
        # a_i[elements_to_it, `1` := sumlog_waveQTL_gen(eps[elements_to_it, `11`] + a_i[all_parents,`1`] + b_pi_no_i[elements_to_it, `1`]
        #                                               ,eps[elements_to_it, `10`] + a_i[all_parents,`0`] + b_pi_no_i[elements_to_it, `0`])]
        # a_i[elements_to_it, `0` := sumlog_waveQTL_gen(eps[elements_to_it, `01`] + a_i[all_parents,`1`] + b_pi_no_i[elements_to_it, `1`]
        #                                               ,eps[elements_to_it, `00`] + a_i[all_parents,`0`] + b_pi_no_i[elements_to_it, `0`])]
        a_i[elements_to_it, 1] = sumlog_waveQTL_gen(eps[elements_to_it, 1] + a_i[all_parents, 1] + b_pi_no_i[elements_to_it, 1]
                                                    ,eps[elements_to_it, 3] + a_i[all_parents,2] + b_pi_no_i[elements_to_it, 2])
        a_i[elements_to_it, 2] = sumlog_waveQTL_gen(eps[elements_to_it, 2] + a_i[all_parents, 1] + b_pi_no_i[elements_to_it, 1]
                                                    ,eps[elements_to_it, 4] + a_i[all_parents,2] + b_pi_no_i[elements_to_it, 2])
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
    # denoms <- (b_i[,`1`]*a_i[,`1`]) + (b_i[,`0`]*a_i[,`0`])
    denoms <- (b_i[,1]*a_i[,1]) + (b_i[,2]*a_i[,2])
    
    # pp_i[,`1` := b_i[,`1`]*a_i[,`1`]/denoms]
    # pp_i[,`0` := b_i[,`0`]*a_i[,`0`]/denoms]
    pp_i[,1] = b_i[,1]*a_i[,1]/denoms
    pp_i[,2] = b_i[,2]*a_i[,2]/denoms
    
    # # checks
    # apply(pp_i,1,sum)
    
    all_parents <- get_parent_indices(1:length(bfs))
    # pp_j_i[,`11` := eps[, `11`] * b_i[,`1`] * b_pi_no_i[,`1`] * a_i[all_parents,`1`] / denoms]
    # pp_j_i[,`10` := eps[, `10`] * b_i[,`1`] * b_pi_no_i[,`0`] * a_i[all_parents,`0`] / denoms]
    # pp_j_i[,`01` := (1-eps[, `11`]) * b_i[,`0`] * b_pi_no_i[,`1`] * a_i[all_parents,`1`] / denoms]
    # pp_j_i[,`00` := (1-eps[, `10`]) * b_i[,`0`] * b_pi_no_i[,`0`] * a_i[all_parents,`0`] / denoms]
    pp_j_i[,1] = eps[, 1] * b_i[,1] * b_pi_no_i[,1] * a_i[all_parents,1] / denoms
    pp_j_i[,2] = eps[, 3] * b_i[,1] * b_pi_no_i[,2] * a_i[all_parents,2] / denoms
    pp_j_i[,3] = (1-eps[, 1]) * b_i[,2] * b_pi_no_i[,1] * a_i[all_parents,1] / denoms
    pp_j_i[,4] = (1-eps[, 3]) * b_i[,2] * b_pi_no_i[,2] * a_i[all_parents,2] / denoms
    
    # # checks
    # apply(pp_j_i,1,sum)
    
    # new_logL <- log(b_i[1,`1`]*a_i[1,`1`] + b_i[1,`0`]*a_i[1,`0`])
    new_logL <- log(b_i[1,1]*a_i[1,1] + b_i[1,2]*a_i[1,2])
    # cat(paste0("New LogL: ", new_logL,"\n"))
  }else{
    #### LOG ADD ####
    # Calculate these qtys in log for now.
    # denoms <- sumlog_waveQTL_gen(b_i[,`1`] + a_i[,`1`]
    #                              , b_i[,`0`] + a_i[,`0`])
    denoms <- sumlog_waveQTL_gen(b_i[,1] + a_i[,1]
                                 , b_i[,2] + a_i[,2])
    
    # pp_i[,`1` := b_i[,`1`] + a_i[,`1`] - denoms]
    # pp_i[,`0` := b_i[,`0`] + a_i[,`0`] - denoms]
    pp_i[,1] = b_i[,1] + a_i[,1] - denoms
    pp_i[,2] = b_i[,2] + a_i[,2] - denoms
    
    # # checks
    # apply(exp(pp_i),1,sum)
    
    all_parents <- get_parent_indices(1:length(bfs))
    # pp_j_i[,`11` := eps[, `11`] + b_i[,`1`] + b_pi_no_i[,`1`] + a_i[all_parents,`1`] - denoms]
    # pp_j_i[,`10` := eps[, `10`] + b_i[,`1`] + b_pi_no_i[,`0`] + a_i[all_parents,`0`] - denoms]
    # pp_j_i[,`01` := eps[, `01`] + b_i[,`0`] + b_pi_no_i[,`1`] + a_i[all_parents,`1`] - denoms]
    # pp_j_i[,`00` := eps[, `00`] + b_i[,`0`] + b_pi_no_i[,`0`] + a_i[all_parents,`0`] - denoms]
    pp_j_i[,1] = eps[, 1] + b_i[,1] + b_pi_no_i[,1] + a_i[all_parents,1] - denoms
    pp_j_i[,2] = eps[, 3] + b_i[,1] + b_pi_no_i[,2] + a_i[all_parents,2] - denoms
    pp_j_i[,3] = eps[, 2] + b_i[,2] + b_pi_no_i[,1] + a_i[all_parents,1] - denoms
    pp_j_i[,4] = eps[, 4] + b_i[,2] + b_pi_no_i[,2] + a_i[all_parents,2] - denoms
    
    # # checks
    # apply(exp(pp_j_i),1,sum)
    
    # new_logL <- sumlog_waveQTL_gen(b_i[1,`1`] + a_i[1,`1`]
    #                                , b_i[1,`0`] +a_i[1,`0`])
    new_logL <- sumlog_waveQTL_gen(b_i[1,1] + a_i[1,1]
                                   , b_i[1,2] +a_i[1,2])
    # cat(paste0("New LogL: ", new_logL,"\n"))
    ## END LOG ADD ##
  }
  
  
  # Loop --------------------------------------------------------------------
  
  while(n < iterations_max && abs(diff) > conv_tol){
    
    old_logL <- new_logL
    # cat(paste0("Iteration: ", n,"\n"))
    # cat(paste0("Old LogL: ", old_logL,"\n"))
    cat(paste0("Iteration: ", n, "\tOld LogL:", round(old_logL,4),"\n"))
    
    # M-step ------------------------------------------------------------------
    # Updates -----------------------------------------------------------------
    
    # a[,`1` := pp_i[,`1`]]
    # a[,`0` := pp_i[,`0`]]
    a[,1] = pp_i[,1]
    a[,2] = pp_i[,2]
    
    # # checks
    # if(!sumlog_use){
    #   apply(a,1,sum)
    # }else{
    #   #### LOG ADD ####
    #   apply(exp(a),1,sum)
    #   ## END LOG ADD ##
    # }
    
    if(!sumlog_use){
      # b[,`11` := pp_j_i[,`11`]/a[all_parents,`1`]]
      # b[,`10` := pp_j_i[,`10`]/a[all_parents,`0`]]
      # b[,`01` := pp_j_i[,`01`]/a[all_parents,`1`]]
      # b[,`00` := pp_j_i[,`00`]/a[all_parents,`0`]]
      b[,1] = pp_j_i[,1]/a[all_parents,1]
      b[,2] = pp_j_i[,2]/a[all_parents,2]
      b[,3] = pp_j_i[,3]/a[all_parents,1]
      b[,4] = pp_j_i[,4]/a[all_parents,2]
    }else{
      #### LOG ADD ####
      # b[,`11` := pp_j_i[,`11`] - a[all_parents,`1`]]
      # b[,`10` := pp_j_i[,`10`] - a[all_parents,`0`]]
      # b[,`01` := pp_j_i[,`01`] - a[all_parents,`1`]]
      # b[,`00` := pp_j_i[,`00`] - a[all_parents,`0`]]
      b[,1] = pp_j_i[,1] - a[all_parents,1]
      b[,2] = pp_j_i[,2] - a[all_parents,2]
      b[,3] = pp_j_i[,3] - a[all_parents,1]
      b[,4] = pp_j_i[,4] - a[all_parents,2]
      ## END LOG ADD ##
    }
    
    
    
    # # checks
    # if(!sumlog_use){
    #   apply(b[,c(1,3)],1,sum)
    #   apply(b[,c(2,4)],1,sum)
    # }else{
    #   #### LOG ADD ####
    #   apply(exp(b[,c(1,3)]),1,sum)
    #   apply(exp(b[,c(2,4)]),1,sum)
    #   ## END LOG ADD ##
    # }
    
    # Update parameters -------------------------------------------------------
    
    # Set up tying groups
    # num_tying_grps = length(tying_groups)
    
    tying_start = integer()
    tying_end = integer()
    num_pheno_tying_grp = integer()
    for(i in 1:(levels - 1)){
      tying_start = c(tying_start,tying_groups[i])
      tying_end = c(tying_end,tying_groups[i + 1] - 1)
      num_pheno_tying_grp = c(num_pheno_tying_grp,(tying_groups[i + 1] - 1) - tying_groups[i] + 1)
    }
    tying_start[levels] = tying_groups[levels]
    tying_end[levels] = length(groups)
    num_pheno_tying_grp[levels] = tying_end[levels] - tying_start[levels] + 1
    
    
    if(!sumlog_use){
      # Tying pi
      # Pi only has 1 coefficient - the first tying group
      # pi <- pp_i[tying_start[1]:tying_end[1], mean(`1`)]
      # pi <- pp_i[1, mean(`1`)]
      pi <- mean(pp_i[1, 1])
      
      # Tying epsilon
      for(i in 1:levels){
        # Need to exclude tree-root: undefined parameter for that element
        start_excl_head = max(2,tying_start[i])
        indices_to_sum <- start_excl_head:tying_end[i]
        parents_indices_to_sum <- get_parent_indices(indices_to_sum)
        # eps[indices_to_sum,] <-
        #   pp_j_i[indices_to_sum , .(`11` = sum(`11`), `10` = sum(`10`), `01` = sum(`01`), `00` = sum(`00`))]/
        #   pp_i[parents_indices_to_sum , .(`11` = sum(`1`), `10` = sum(`0`), `01` = sum(`1`), `00` = sum(`0`))]
        eps[indices_to_sum,] <-
          apply(pp_j_i[indices_to_sum],2,sum)/c(apply(pp_i[parents_indices_to_sum],2,sum),apply(pp_i[parents_indices_to_sum],2,sum))
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
      # pi <- pp_i[1, (`1`)]
      # pi_1_minus <- pp_i[1, (`0`)]
      pi <- pp_i[1, 1]
      pi_1_minus <- pp_i[1, 2]
      
      
      # Tying epsilon - need to exclude tree-root: undefined parameter for that element
      for(i in 1:levels){
        # Need to exclude tree-root: undefined parameter for that element
        start_excl_head = max(2,tying_start[i])
        indices_to_sum <- start_excl_head:tying_end[i]
        parents_indices_to_sum <- get_parent_indices(indices_to_sum)
        
        if(tying_end[i] < start_excl_head){
          next
        }else{
          num_pheno <- length(indices_to_sum)
          
          if(num_pheno == 1){
            # eps[indices_to_sum,] <- b[indices_to_sum, .(`11`, `01`, `10`, `00`)]
            eps[indices_to_sum,] <- b[indices_to_sum, c(1,3,2,4)]
            
          }else{
            # pp_values <- pp_i[parents_indices_to_sum,.(`1`,`1`,`0`,`0`)]
            # pp_joint_values <- pp_j_i[indices_to_sum, .(`11`,`01`,`10`,`00`)]
            pp_values <- pp_i[parents_indices_to_sum, c(1,1,2,2)]
            pp_joint_values <- pp_j_i[indices_to_sum, c(1,3,2,4)]
            
            for(j in 1:4){
              # eps_logs <- pp_joint_values[1][[j]]
              # eps_denom_logs <- pp_values[1][[j]]
              eps_logs <- pp_joint_values[1,j]
              eps_denom_logs <- pp_values[1,j]
              
              for(k in 1:(num_pheno - 1)){
                # eps_logs <- sumlog_waveQTL_gen(eps_logs,pp_joint_values[k + 1][[j]])
                # eps_denom_logs <- sumlog_waveQTL_gen(eps_denom_logs,pp_values[k + 1][[j]])
                eps_logs <- sumlog_waveQTL_gen(eps_logs,pp_joint_values[k + 1,j])
                eps_denom_logs <- sumlog_waveQTL_gen(eps_denom_logs,pp_values[k + 1,j])
              }
              # eps[indices_to_sum][[j]] <- eps_logs - eps_denom_logs
              eps[indices_to_sum,j] <- eps_logs - eps_denom_logs
            }
            
          }
        }
      }
      
    }
    
    # E-step ------------------------------------------------------------------
    # Up step -----------------------------------------------------------------
    for(i in rev((1:levels)[-1])){
      elements_to_it <- which(groups == i)
      
      # Last level only
      if(i == rev((1:levels)[-1])[1]){
        # b_i[elements_to_it, "1" := bfs[elements_to_it]]
        b_i[elements_to_it, 1] = bfs[elements_to_it]
        if(!sumlog_use){
          # b_i[elements_to_it, "0" := 1]
          b_i[elements_to_it, 2] = 1
        }else{
          #### LOG ADD ####
          # b_i[elements_to_it, "0" := log(1)]
          b_i[elements_to_it, 2] = log(1)
          ## END LOG ADD ##
        }
      }
      
      if(!sumlog_use){
        # b_i_pi[elements_to_it, "1" := eps[elements_to_it, `11`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `11`]) * b_i[elements_to_it, `0`]]
        # b_i_pi[elements_to_it, "0" := eps[elements_to_it, `10`] * b_i[elements_to_it, `1`] + (1 - eps[elements_to_it, `10`]) * b_i[elements_to_it, `0`]]
        b_i_pi[elements_to_it, 1] = eps[elements_to_it, 1] * b_i[elements_to_it, 1] + (1 - eps[elements_to_it, 1]) * b_i[elements_to_it, 2]
        b_i_pi[elements_to_it, 2] = eps[elements_to_it, 3] * b_i[elements_to_it, 1] + (1 - eps[elements_to_it, 3]) * b_i[elements_to_it, 2]
      }else{
        #### LOG ADD ####
        # log(b_i_pi) = log(eps1*b_i_1 + eps0*b_i_0)
        # = sumlog(log(eps1*b_i_1) + log(eps0*b_i_0))
        # = sumlog(log(eps1) + b_i_1, log(eps0) + b_i_0) as b_i already logged.
        # b_i_pi[elements_to_it, "1" := sumlog_waveQTL_gen(eps[elements_to_it, `11`] + b_i[elements_to_it, `1`]
        #                                                  ,eps[elements_to_it, `01`] + b_i[elements_to_it, `0`])]
        # b_i_pi[elements_to_it, "0" := sumlog_waveQTL_gen(eps[elements_to_it, `10`] + b_i[elements_to_it, `1`]
        #                                                  ,eps[elements_to_it, `00`] + b_i[elements_to_it, `0`])]
        b_i_pi[elements_to_it, 1] = sumlog_waveQTL_gen(eps[elements_to_it, 1] + b_i[elements_to_it, 1]
                                                       ,eps[elements_to_it, 2] + b_i[elements_to_it, 2])
        b_i_pi[elements_to_it, 2] = sumlog_waveQTL_gen(eps[elements_to_it, 3] + b_i[elements_to_it, 1]
                                                       ,eps[elements_to_it, 4] + b_i[elements_to_it, 2])
        ## END LOG ADD ##
      }
      
      
      
      all_parents <- get_parent_indices(elements_to_it)
      parents <- unique(all_parents)
      for(j in parents){
        children <- get_child_indices(j)[,1]
        if(!sumlog_use){
          # b_i[j, "1" := bfs[j]*prod(b_i_pi[children,`1`])]
          # b_i[j, "0" := 1*prod(b_i_pi[children,`0`])]
          b_i[j, 1] = bfs[j]*prod(b_i_pi[children,1])
          b_i[j, 2] = 1*prod(b_i_pi[children,2])
        }else{
          #### LOG ADD ####
          # b_i[j, "1" := bfs[j] + sum(b_i_pi[children,`1`])]
          # b_i[j, "0" := sum(b_i_pi[children,`0`])]
          b_i[j, 1] = bfs[j] + sum(b_i_pi[children,1])
          b_i[j, 2] = sum(b_i_pi[children,2])
          ## END LOG ADD ##
        }
        
      }
      
      if(!sumlog_use){
        # b_pi_no_i[elements_to_it, "1" := b_i[all_parents,`1`]/b_i_pi[elements_to_it,`1`]]
        # b_pi_no_i[elements_to_it, "0" := b_i[all_parents,`0`]/b_i_pi[elements_to_it,`0`]]
        b_pi_no_i[elements_to_it, 1] = b_i[all_parents,1]/b_i_pi[elements_to_it,1]
        b_pi_no_i[elements_to_it, 2] = b_i[all_parents,2]/b_i_pi[elements_to_it,2]
      }else{
        #### LOG ADD ####
        # b_pi_no_i[elements_to_it, "1" := b_i[all_parents,`1`] - b_i_pi[elements_to_it,`1`]]
        # b_pi_no_i[elements_to_it, "0" := b_i[all_parents,`0`] - b_i_pi[elements_to_it,`0`]]
        b_pi_no_i[elements_to_it, 1] = b_i[all_parents,1] - b_i_pi[elements_to_it,1]
        b_pi_no_i[elements_to_it, 2] = b_i[all_parents,2] - b_i_pi[elements_to_it,2]
        ## END LOG ADD ##
      }
      
      
    }
    
    # Down step ---------------------------------------------------------------
    for(i in 1:levels){
      elements_to_it <- which(groups == i)
      if(i == (1:levels)[1]){
        if(!sumlog_use){
          # a_i[elements_to_it, c("1","0") := .(pi, 1 - pi)]
          a_i[elements_to_it, 1] = pi
          a_i[elements_to_it, 2] = 1 - pi
        }else{
          #### LOG ADD ####
          # a_i[elements_to_it, c("1","0") := .(pi, pi_1_minus)]
          a_i[elements_to_it, 1] = pi
          a_i[elements_to_it, 2] = pi_1_minus
          ## END LOG ADD ##
        }
      }else{
        all_parents <- get_parent_indices(elements_to_it)
        if(!sumlog_use){
          # a_i[elements_to_it, `1` := (eps[elements_to_it, `11`] * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
          #       (eps[elements_to_it, `10`] * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
          # a_i[elements_to_it, `0` := ((1 - eps[elements_to_it, `11`]) * a_i[all_parents,`1`] * b_pi_no_i[elements_to_it, `1`]) +
          #       ((1 - eps[elements_to_it, `10`]) * a_i[all_parents,`0`] * b_pi_no_i[elements_to_it, `0`])]
          a_i[elements_to_it, 1] = (eps[elements_to_it, 1] * a_i[all_parents,1] * b_pi_no_i[elements_to_it, 1]) +
            (eps[elements_to_it, 3] * a_i[all_parents,2] * b_pi_no_i[elements_to_it, 2])
          a_i[elements_to_it, 2] = ((1 - eps[elements_to_it, 1]) * a_i[all_parents,1] * b_pi_no_i[elements_to_it, 1]) +
            ((1 - eps[elements_to_it, 3]) * a_i[all_parents,2] * b_pi_no_i[elements_to_it, 2])
        }else{
          #### LOG ADD ####
          # a_i[elements_to_it, `1` := sumlog_waveQTL_gen(eps[elements_to_it, `11`] + a_i[all_parents,`1`] + b_pi_no_i[elements_to_it, `1`]
          #                                               ,eps[elements_to_it, `10`] + a_i[all_parents,`0`] + b_pi_no_i[elements_to_it, `0`])]
          # a_i[elements_to_it, `0` := sumlog_waveQTL_gen(eps[elements_to_it, `01`] + a_i[all_parents,`1`] + b_pi_no_i[elements_to_it, `1`]
          #                                               ,eps[elements_to_it, `00`] + a_i[all_parents,`0`] + b_pi_no_i[elements_to_it, `0`])]
          a_i[elements_to_it, 1] = sumlog_waveQTL_gen(eps[elements_to_it, 1] + a_i[all_parents,1] + b_pi_no_i[elements_to_it, 1]
                                                      ,eps[elements_to_it, 3] + a_i[all_parents,2] + b_pi_no_i[elements_to_it, 2])
          a_i[elements_to_it, 2] = sumlog_waveQTL_gen(eps[elements_to_it, 2] + a_i[all_parents,1] + b_pi_no_i[elements_to_it, 1]
                                                      ,eps[elements_to_it, 4] + a_i[all_parents,2] + b_pi_no_i[elements_to_it, 2])
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
      # denoms <- (b_i[,`1`]*a_i[,`1`]) + (b_i[,`0`]*a_i[,`0`])
      denoms <- (b_i[,1]*a_i[,1]) + (b_i[,2]*a_i[,2])
      
      # pp_i[,`1` := b_i[,`1`]*a_i[,`1`]/denoms]
      # pp_i[,`0` := b_i[,`0`]*a_i[,`0`]/denoms]
      pp_i[,1] = b_i[,1]*a_i[,1]/denoms
      pp_i[,2] = b_i[,2]*a_i[,2]/denoms
      
      # # checks
      # apply(pp_i,1,sum)
      
      all_parents <- get_parent_indices(1:length(bfs))
      # pp_j_i[,`11` := eps[, `11`] * b_i[,`1`] * b_pi_no_i[,`1`] * a_i[all_parents,`1`] / denoms]
      # pp_j_i[,`10` := eps[, `10`] * b_i[,`1`] * b_pi_no_i[,`0`] * a_i[all_parents,`0`] / denoms]
      # pp_j_i[,`01` := (1-eps[, `11`]) * b_i[,`0`] * b_pi_no_i[,`1`] * a_i[all_parents,`1`] / denoms]
      # pp_j_i[,`00` := (1-eps[, `10`]) * b_i[,`0`] * b_pi_no_i[,`0`] * a_i[all_parents,`0`] / denoms]
      pp_j_i[,1] = eps[, 1] * b_i[,1] * b_pi_no_i[,1] * a_i[all_parents,1] / denoms
      pp_j_i[,2] = eps[, 3] * b_i[,1] * b_pi_no_i[,2] * a_i[all_parents,2] / denoms
      pp_j_i[,3] = (1-eps[, 1]) * b_i[,2] * b_pi_no_i[,1] * a_i[all_parents,1] / denoms
      pp_j_i[,4] = (1-eps[, 3]) * b_i[,2] * b_pi_no_i[,2] * a_i[all_parents,2] / denoms
      
      # # checks
      # apply(pp_j_i,1,sum)
      
      # new_logL <- log(b_i[1,`1`]*a_i[1,`1`] + b_i[1,`0`]*a_i[1,`0`])
      new_logL <- log(b_i[1,1]*a_i[1,1] + b_i[1,2]*a_i[1,2])
      # cat(paste0("New LogL: ", new_logL,"\n"))
    }else{
      #### LOG ADD ####
      # Calculate these qtys in log for now.
      # denoms <- sumlog_waveQTL_gen(b_i[,`1`] + a_i[,`1`]
      #                              , b_i[,`0`] + a_i[,`0`])
      denoms <- sumlog_waveQTL_gen(b_i[,1] + a_i[,1]
                                   , b_i[,2] + a_i[,2])
      
      # pp_i[,`1` := b_i[,`1`] + a_i[,`1`] - denoms]
      # pp_i[,`0` := b_i[,`0`] + a_i[,`0`] - denoms]
      pp_i[,1] = b_i[,1] + a_i[,1] - denoms
      pp_i[,2] = b_i[,2] + a_i[,2] - denoms
      
      # # checks
      # apply(exp(pp_i),1,sum)
      
      all_parents <- get_parent_indices(1:length(bfs))
      # pp_j_i[,`11` := eps[, `11`] + b_i[,`1`] + b_pi_no_i[,`1`] + a_i[all_parents,`1`] - denoms]
      # pp_j_i[,`10` := eps[, `10`] + b_i[,`1`] + b_pi_no_i[,`0`] + a_i[all_parents,`0`] - denoms]
      # pp_j_i[,`01` := eps[, `01`] + b_i[,`0`] + b_pi_no_i[,`1`] + a_i[all_parents,`1`] - denoms]
      # pp_j_i[,`00` := eps[, `00`] + b_i[,`0`] + b_pi_no_i[,`0`] + a_i[all_parents,`0`] - denoms]
      pp_j_i[,1] = eps[, 1] + b_i[,1] + b_pi_no_i[,1] + a_i[all_parents,1] - denoms
      pp_j_i[,2] = eps[, 3] + b_i[,1] + b_pi_no_i[,2] + a_i[all_parents,2] - denoms
      pp_j_i[,3] = eps[, 2] + b_i[,2] + b_pi_no_i[,1] + a_i[all_parents,1] - denoms
      pp_j_i[,4] = eps[, 4] + b_i[,2] + b_pi_no_i[,2] + a_i[all_parents,2] - denoms
      
      # # checks
      # apply(exp(pp_j_i),1,sum)
      
      # new_logL <- sumlog_waveQTL_gen(b_i[1,`1`] + a_i[1,`1`]
      #                                , b_i[1,`0`] +a_i[1,`0`])
      new_logL <- sumlog_waveQTL_gen(b_i[1,1] + a_i[1,1]
                                     , b_i[1,2] +a_i[1,2])
      # cat(paste0("New LogL: ", new_logL,"\n"))
      ## END LOG ADD ##
    }
    
    
    diff <- new_logL - old_logL
    
    n <- n + 1
    # print(exp(eps[tying_groups,]))
    logl_list[[inits_loop]] <- c(logl_list[[inits_loop]], new_logL)
    tmp_eps[[n]] <- exp(eps[tying_groups,])
  }
  eps_vals[[inits_loop]] <- tmp_eps
  inits_loop=inits_loop+1
  cat(paste0("New LogL: ", new_logL,"\n"))

  # WQtl step for scaling coefficient ---------------------------------------
  
  if(!sumlog_use){
    stop("Scaling coefficient not made for sumlog_use == False. Set sumlog_use == True")
  }
  
  logLR_sc = 0
  # elements_to_it <- which(groups == 1)
  n_obllikli = 0
  pp_sl = 0
  pi_sc = log(init_pi)
  pi_1_minus_sc = log(1-init_pi)
  logL_sc_list = 0
  pi_sc_list = exp(pi_sc)

  # For each location in the group
  # for(loc in elements_to_it){
    
  # }
  log_pi_bf = bfs_sc + pi_sc
  log_den = sumlog_waveQTL_gen(log_pi_bf, pi_1_minus_sc)
  pp_sl = pp_sl + exp(log_pi_bf - log_den)
  n_obllikli = n_obllikli + log_den  
  o_obllikli = n_obllikli

  n = 0
  diff <- Inf
  while(n < iterations_max){

    cat(paste0("Scaling coeff iteration: ", n, "\tOld LogL:", round(o_obllikli,4),"\n"))
    pi_sc = log(pp_sl)
    pi_1_minus_sc = log(1 - pp_sl)
    pi_sc_list <- c(pi_sc_list,exp(pi_sc))
    # logpi  = log(pi_s[g])
    # log1pi = log(1-pi)
    n_obllikli = 0
    pp_sl = 0

    # for(loc in elements_to_it){
    # }
    log_pi_bf = bfs_sc + pi_sc
    log_den = sumlog_waveQTL_gen(log_pi_bf, pi_1_minus_sc)
    pp_sl = pp_sl + exp(log_pi_bf - log_den)
    n_obllikli = n_obllikli + log_den

    diff = n_obllikli - o_obllikli
    logL_sc_list <- c(logL_sc_list, n_obllikli)
    
    if(diff < conv_tol){
      pi_sc <- log(pp_sl)
      pi_sc_list <- c(pi_sc_list,exp(pi_sc))
      break
    }else{
      o_obllikli = n_obllikli
      n = n + 1
    }

  }

  logLR_sc = logLR_sc + n_obllikli
  cat(paste0("LogL after scaling coeff: ", round(logLR_sc,4),"\n"))

  # Augment outputs to consider scaling coef ----------------------------------
  
  cat(paste0("SC LogL: ", round(logLR_sc,5) ,"\n"))
  cat(paste0("HMT LogL: ", round(new_logL,5) ,"\n"))
  cat(paste0("Total LogL: ", round(logLR_sc + new_logL,5) ,"\n"))

  return(
    list(
      eps_list = eps_vals
      ,pi_list = c(exp(pi_sc), exp(pi))
      ,logl_list = append(logl_list, logL_sc_list)
      ,logLR = logLR_sc + new_logL
    )
  )
  
}