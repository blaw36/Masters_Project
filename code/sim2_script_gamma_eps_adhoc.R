### Check the simulation
# Epsilon was generated by:
# Eps_11
c(rep(0.9,(num_tying_grps-1)),0)
# Eps_10
c(rep(0.1,(num_tying_grps-1)),0)

# Gamma simul
gamma_seq_tree <- gamma_seq[-1]
length(gamma_seq_tree)

gamma_and_parent <- matrix(c(gamma_seq_tree[get_parent_indices(1:1023)],gamma_seq_tree)
                           , nrow = 2, ncol = 1023
                           , byrow = T)
gamma_and_parent_dt <- as.data.table(t(gamma_and_parent))
setnames(gamma_and_parent_dt,names(gamma_and_parent_dt),c("Parent","Child"))
gamma_and_parent_dt[,"Transition" := paste0(Child,Parent)]
gamma_and_parent_dt[,"TreeLvl" := floor(log2(.I))]

# Exclude the root
gamma_and_parent_dt <- copy(gamma_and_parent_dt[-1])
gamma_and_parent_stats <- dcast.data.table(
  gamma_and_parent_dt[,.N,by = .(Transition,TreeLvl)]
  ,formula = TreeLvl ~ Transition
  ,value.var = "N"
  ,fill = 0
)
gamma_and_parent_stats[,"Total" := apply(.SD,1,sum),.SDcols = 2:5]


# Top levels will be noisy when no tying (not a large neough sample to get, say. 90% transition probabilities out of, say, 4 transitions)

# Simulated data
gamma_and_parent_stats[,c("sim_eps_11","sim_eps_10") :=
                         .(`11`/(`11`+`01`)
                           ,`10`/(`10`+`00`))]
gamma_and_parent_stats[,c("param_eps_11","param_eps_10") :=
                         .(c(rep(0.9,8),0)
                           ,c(rep(0.1,8),0))]

# Out of a dataset which has 2^i-1 rows...
extract_tied_params <- function(data){
  return_data <- c()
  for(i in 0:floor(log2(length(data)))){
    return_data <- c(return_data,data[2^i])
  }
  return(return_data)
}

# Eps results
gamma_and_parent_stats[,c("act_eps_11","act_eps_10") :=
                         .(extract_tied_params(eps_11_file)[-1]
                           ,extract_tied_params(eps_10_file)[-1])]

## This is an issue, there is a disconnect between:
  # - Parameterised epsilon parameters
  # - Simulated (sample) epsilons based on the simulated gammas
  # - Actual (algorithm) retrieved epsilons
