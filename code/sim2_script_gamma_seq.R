n_ind = 70
# n_pheno = 1024
# tying_grp = c(1,2,3,5,9,17,33,65,129,257,513)
n_pheno = 64
tying_grp = c(1,2,33)
param_pi_00 = 1
param_pi_11 = 0
coeff_mu = 0
coeff_beta = 2
param_gi_prob = 0.4
param_sigma_beta = 0.5
num_sims = 100
seed = 20

# Round to nearest 100th, for simplicity
num_tying_grps <- length(tying_grp)
set.seed(seed)

# Generate gamma sequence -------------------------------------------------

#### ~~~ Default method
# Eps_11 (elements 2 -> 1023 of tree)
grped_eps_11 <- round(runif(n=num_tying_grps)*100)/100
param_eps_11 <- numeric()
for(i in 1:(num_tying_grps - 1)){
  indx_start <- tying_grp[i]
  indx_end <- tying_grp[i+1]
  param_eps_11[indx_start:indx_end] <- grped_eps_11[i]
}
param_eps_11[tying_grp[num_tying_grps]:n_pheno] <- grped_eps_11[num_tying_grps]
param_eps_11 <- param_eps_11[-(1:2)]

# Eps_10 (elements 2 -> 1023 of tree)
grped_eps_10 <- round(runif(n=num_tying_grps)*100)/100
param_eps_10 <- numeric()
for(i in 1:(num_tying_grps - 1)){
  indx_start <- tying_grp[i]
  indx_end <- tying_grp[i+1]
  param_eps_10[indx_start:indx_end] <- grped_eps_10[i]
}
param_eps_10[tying_grp[num_tying_grps]:n_pheno] <- grped_eps_10[num_tying_grps]
param_eps_10 <- param_eps_10[-(1:2)]

#### ~~~ New method



# Experiment around with gamma sequence generation ------------------------

gamma_seq <- numeric()
rand_seq <- runif(n_pheno)

# Scaling coefficient
gamma_seq[1] <- ifelse(rand_seq[1] < param_pi_00, 1, 0)

# Head of tree
gamma_seq[2] <- ifelse(rand_seq[2] < param_pi_11, 1, 0)

# i is the index of tree, where i = 1 is the head of the tree
# Using this notation because that's how 'get_parent_indices' has been written
for(i in 2:(n_pheno - 1)){
  indx <- i
  parent_indx <- get_parent_indices(indx)
  parent_gamma <- gamma_seq[parent_indx + 1]

  if(parent_gamma == 1){
    sl_prob <- param_eps_11[indx - 1]
  }else if(parent_gamma == 0){
    sl_prob <- param_eps_10[indx - 1]
  }

  gamma_seq[i+1] <- ifelse(rand_seq[i+1] < sl_prob, 1, 0)
}

# Sort each tree level so that all 0s on one side, all 1s on the other
num_tree_levels <- log2(n_pheno)
for(i in 1:(num_tree_levels-1)){
  grp_start <- 2^i + 1
  grp_finish <- 2^(i+1)
  gamma_seq[grp_start:grp_finish] <- sort(gamma_seq[grp_start:grp_finish])
}
# final grp - make all 0s
grp_start <- tying_grp[num_tying_grps]
grp_finish <- n_pheno
gamma_seq[grp_start:grp_finish] <- 0


### Step 2: Generate beta,mu,eps,g,y ----
beta_seq <- rep(0,n_pheno)
beta_seq[which(gamma_seq == 1)] <- coeff_beta

mu_seq <- rep(0,n_pheno)
mu_seq[which(gamma_seq == 1)] <- coeff_mu
mu_mtx <- matrix(rep(mu_seq,70),nrow = n_ind,ncol = n_pheno,byrow = T)

eps_seq <- rnorm(n_pheno*n_ind,mean = 0,sd = sqrt(param_sigma_beta))
eps_mtx <- matrix(eps_seq,nrow = n_ind,ncol = n_pheno,byrow = T)

g_seq <- rbinom(n_ind,size = 2,prob = param_gi_prob)
# cat(c("sim2","A","A",g_seq), file = paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/sim2.cis.geno"))
beta_mtx <- g_seq %*% t(beta_seq)

y_mtx <- mu_mtx + beta_mtx + eps_mtx
# write.table(y_mtx, file= paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sim2_WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)
# cat(rep(1,n_pheno), file = paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/use_all.txt"))

### stop experiment



# Plots -------------------------------------------------------------------

vector_centeriser <- function(vect){
  in_between_zeros <- length(vect) - 1

  res_vect <- c(0,vect[1],0)
  if(in_between_zeros > 0){
    for(i in 1:in_between_zeros){
      res_vect <- c(res_vect,vect[i + 1],0)
    }
  }
  return(res_vect)
}

# Kind of a start to visualising tree coefficients - without the scaling coeff (11 rows too large for plot)
gamma_seq_plot <- gamma_seq[-1]
par(mfrow = c(log2(n_pheno),1),mar = c(2,2,2,2))
plot(vector_centeriser(gamma_seq_plot[1]),type = "h",ylab = "",axes = F
     ,main = "Gamma value, by tree level (0 or 1)")
for(i in 1:(log2(n_pheno)-1)){
  plot(vector_centeriser(gamma_seq_plot[(2^i):(2^(i+1) - 1)]),type = "h",ylab = "",axes = F)
}

y_mtx
beta_mtx

# y means, by group
y_mtx_w_group <- cbind(y_mtx,g_seq)
y_mtx_means_g0 <- apply(y_mtx[g_seq == 0,],2,mean)
y_mtx_means_g1 <- apply(y_mtx[g_seq == 1,],2,mean)
y_mtx_means_g2 <- apply(y_mtx[g_seq == 2,],2,mean)

y_axis_bounds <- c(min(min(y_mtx_means_g0),min(y_mtx_means_g1),min(y_mtx_means_g2))*1.1
                   ,max(max(y_mtx_means_g0),max(y_mtx_means_g1),max(y_mtx_means_g2))*1.1)
xval = 1:n_pheno
graphics.off()
plot(1,1, type = "n", xlab = "Wavelet scale-loc", ylab = "Mean WC value", main = "Mean WCs by group", xaxt = "n"
     ,xlim = c(0,n_pheno)
     ,ylim = y_axis_bounds)
axis(1,at = seq(0,n_pheno,by = 16),labels = seq(0,n_pheno,by = 16),las = 2)
lines(y_mtx_means_g0,col="red")
lines(y_mtx_means_g1,col="blue")
lines(y_mtx_means_g2,col="green")
par(xpd=TRUE)
legend(1,6,legend=c("g=0", "g=1" ,"g=2")
       # ,"topleft"
       ,col=c("red", "blue" ,"green")
       ,lty = 1
       ,box.lty = 0
       ,cex = 0.5)


# Y, beta means
y_mtx_means <- apply(y_mtx,MARGIN = 2,FUN = mean)
beta_mtx_means <- apply(beta_mtx,MARGIN = 2,FUN = mean)

plot(y_mtx_means,type="l", main = "WC means, by scale-loc")
plot(beta_mtx_means,type="l", main = "beta means, by scale-loc")

# Y means, by level
y_mtx_means_plot <- y_mtx_means[-1]
y_min <- min(y_mtx_means_plot)
y_max <- max(y_mtx_means_plot)
par(mfrow = c(log2(n_pheno),1)
    ,mar = c(2,2,2,2))
plot(vector_centeriser(y_mtx_means_plot[1]),type = "h",ylab = "",axes = F
     ,main = "WC means by tree level",ylim = c(y_min,y_max))
for(i in 1:(log2(n_pheno)-1)){
  plot(vector_centeriser(y_mtx_means_plot[(2^i):(2^(i+1) - 1)])
       ,type = "h",ylab = "",axes = F,ylim = c(y_min,y_max))
}

# Beta means, by level
beta_mtx_means_plot <- beta_mtx_means[-1]
beta_y_min <- min(beta_mtx_means_plot)
beta_y_max <- max(beta_mtx_means_plot)
par(mfrow = c(log2(n_pheno),1)
    ,mar = c(2,2,2,2))
plot(vector_centeriser(beta_mtx_means_plot[1]),type = "h",ylab = "",axes = F
     ,main = "beta means by tree level", ylim = c(beta_y_min, beta_y_max))
for(i in 1:(log2(n_pheno)-1)){
  plot(vector_centeriser(beta_mtx_means_plot[(2^i):(2^(i+1) - 1)]),type = "h",ylab = "",axes = F
       , ylim = c(beta_y_min, beta_y_max))
}
