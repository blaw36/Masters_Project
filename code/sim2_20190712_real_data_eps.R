library(ggplot2)
library(reshape2)

p_tying_grp = c(1,2,3,5,9,17,33,65,129,513)

eps_10 <- exp(as.matrix(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/tree_tie_noQT.fph.eps_10.txt"))[,-1]))[,p_tying_grp]
eps_11 <- exp(as.matrix(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/tree_tie_noQT.fph.eps_11.txt"))[,-1]))[,p_tying_grp]

pp <- exp(as.matrix(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/tree_tie_noQT.fph.pp.txt"))[,-1]))[,p_tying_grp]
pp_joint_01 <- exp(as.matrix(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/tree_tie_noQT.fph.pp_joint_01.txt"))[,-1]))[,p_tying_grp]
pp_joint_10 <- exp(as.matrix(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/tree_tie_noQT.fph.pp_joint_10.txt"))[,-1]))[,p_tying_grp]
pp_joint_11 <- exp(as.matrix(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/tree_tie_noQT.fph.pp_joint_11.txt"))[,-1]))[,p_tying_grp]

# Extract params from each tying group

eps_11_df <- data.frame(eps_11)
names(eps_11_df) <- c(paste0("g",1:10))
eps_11_df$geno_num = row.names(eps_11_df)
ggplot(melt(eps_11_df,id.vars = "geno_num")) +
  geom_boxplot(aes(x = factor(variable), y = value)) +
  geom_point(position=position_dodge(width=0.75)
             ,aes(x = factor(variable), y = value, group = factor(variable))) +
  xlab("Tying group (top of tree to bottom)") +
  ylab("Probability") +
  ggtitle("Distribution of Epsilon_11 parameters across real data")


eps_10_df <- data.frame(eps_10)
names(eps_10_df) <- c(paste0("g",1:10))
eps_10_df$geno_num = row.names(eps_10_df)
ggplot(melt(eps_10_df,id.vars = "geno_num")) +
  geom_boxplot(aes(x = factor(variable), y = value)) +
  geom_point(position=position_dodge(width=0.75)
             ,aes(x = factor(variable), y = value, group = factor(variable))) +
  xlab("Tying group (top of tree to bottom)") +
  ylab("Probability") +
  ggtitle("Distribution of Epsilon_10 parameters across real data")

ggplot(melt(eps_11_df[11,],id.vars = "geno_num")) +
  geom_line(position=position_dodge(width=0.75)
             ,aes(x = factor(variable), y = value, group = factor(geno_num))) +
  geom_point(aes(x = factor(variable), y = value)) +
  xlab("Tying group (top of tree to bottom)") +
  ylab("Probability") +
  ggtitle("Distribution of Epsilon_11 parameters across real data")

ggplot(melt(eps_10_df[11,],id.vars = "geno_num")) +
  geom_line(position=position_dodge(width=0.75)
             ,aes(x = factor(variable), y = value, group = factor(geno_num))) +
  geom_point(aes(x = factor(variable), y = value)) +
  xlab("Tying group (top of tree to bottom)") +
  ylab("Probability") +
  ggtitle("Distribution of Epsilon_10 parameters across real data")

# readLines(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/use.txt"))
# plot(as.matrix(read.table(paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/chr17.10160989.10162012.pheno.dat")))[1,],type = "l")
# plot(apply(as.matrix(read.table(paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/chr17.10160989.10162012.pheno.dat"))),2,mean),type = "l")
# which(apply(as.matrix(read.table(paste0("~/Cpp/WaveQTL_HMT/data/dsQTL/chr17.10160989.10162012.pheno.dat"))),2,mean) >= 2)
