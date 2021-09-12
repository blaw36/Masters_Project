## EM algorithm with HMT test script.
### This script allows INIT PARAMETER TO BE VARIED AND PLAYED WITH

wqtl <- function(
  logBFs,
  # groups,
  tying_groups = c(1,c(1,2,4,8,16,32,64,128,256,512)+1),
  iterations_max,
  conv_tol,
  # sumlog_use, # we will only allow sumlog method here as per c++
  init_pi
){
  
  eps_vals <- list()
  logl_list <- list()
  
  inits_loop=1
  
  # init_try = 0.5
  
  logl_list[[inits_loop]] <- 0
  tmp_pi <- list()
  
  
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
  
  bfs = logBFs
  
  # if(is.null(tying_groups)){
  #   tying_groups = c(1:length(groups))
  # }
  
  # levels = tree_levels(length(bfs))
  groups = group_converter(tying_groups,logBFs)
  levels = uniqueN(groups)

  pp_sl <- matrix(rep(NA_real_,levels*2),ncol = 2) # It's kept at a group level in c++ code
  n_obllikli <- rep(NA_real_,levels)
  o_obllikli <- rep(NA_real_,levels)
  logLR <- rep(NA_real_,levels)
    
  # Parameters
  pi_s <- rep(init_pi,levels)

  # Loop over each group ----------------------------------------------------

  logLR = 0
  total_iter_ticker = 1
  for(g in 1:levels){
    
    elements_to_it <- which(groups == g)

    n_obllikli[g] = 0
    pp_sl[g,1] = 0

    # For each location in the group
    for(loc in elements_to_it){
      log_pi_bf = bfs[loc] + log(pi_s[g])
      log_den = sumlog_waveQTL_gen(log_pi_bf, log(1-pi_s[g]))
      pp_sl[g,1] = pp_sl[g,1] + exp(log_pi_bf - log_den)
      n_obllikli[g] = n_obllikli[g] + log_den
    }
    o_obllikli[g] = n_obllikli[g]

    n = 0
    diff <- Inf
    # while(n < iterations_max && abs(diff) > conv_tol){
    while(n < iterations_max){

      cat(paste0("Group: ",g,"\tIteration: ", n, "\tOld LogL:", round(o_obllikli[g],4),"\n"))
      pi_s[g] = pp_sl[g,1]/length(elements_to_it)
      tmp_pi[[total_iter_ticker]] <- pi_s
      # logpi  = log(pi_s[g])
      # log1pi = log(1-pi)
      n_obllikli[g] = 0
      pp_sl[g,1] = 0

      for(loc in elements_to_it){
        log_pi_bf = bfs[loc] + log(pi_s[g])
        log_den = sumlog_waveQTL_gen(log_pi_bf, log(1-pi_s[g]))
        pp_sl[g,1] = pp_sl[g,1] + exp(log_pi_bf - log_den)
        n_obllikli[g] = n_obllikli[g] + log_den
      }

      diff = n_obllikli[g] - o_obllikli[g]
      logl_list[[total_iter_ticker]] <- n_obllikli
      total_iter_ticker = total_iter_ticker + 1
      
      if(diff < conv_tol){
        pi_s[g] <- pp_sl[g,1]/length(elements_to_it)
        tmp_pi[[total_iter_ticker]] <- pi_s
        break
      }else{
        o_obllikli[g] = n_obllikli[g]
        n = n + 1
      }

    }

    logLR = logLR + n_obllikli[g];
    cat(paste0("LogL after ",g," groups: ", round(logLR,4),"\n"))

  }

  return(
    list(
      tmp_pi = tmp_pi
      ,logl_list = logl_list
      ,logLR = logLR
      )
    )

}