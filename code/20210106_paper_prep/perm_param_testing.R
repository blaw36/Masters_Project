library(ggplot2)
library(reshape2)

setwd("C:/Users/brend/Documents/Cpp/WaveQTL_HMT_wperm/test/dsQTL/")
# Try on 100 permutations


# Tests on 100 perms ------------------------------------------------------

## Current run
# (niter = 1000, tol = 0.0005)
start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210601_hmt_curr -f 1024 -numPerm 100 -fph 3 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
curr_run <- end_time - start_time
# user  system elapsed
# 0.02    0.03  186.94

### Change some of the parameters for the PERMUTATION runs
## Reduce # iters to 100
# (niter = 100, tol = 0.0005)
start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210601_hmt_less_iter -f 1024 -numPerm 100 -fph 3 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
less_iter <- end_time - start_time
# user  system elapsed
# 0.05    0.03  120.0

## Increase tol to 0.005
# (niter = 1000, tol = 0.005)
start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210601_hmt_more_tol -f 1024 -numPerm 100 -fph 3 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
more_tol <- end_time - start_time
# user  system elapsed
# 0.05    0.07   53.07

## Reduce # iters to 100 AND increase tol to 0.005
# (niter = 100, tol = 0.005)
start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210601_hmt_both -f 1024 -numPerm 100 -fph 3 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
both <- end_time - start_time
# user  system elapsed
# 0.02    0.00   51.11

## Increase # iters to 1000 AND decrease tol to 0.00005
# (niter = 100, tol = 0.005)
start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210601_hmt_precise -f 1024 -numPerm 100 -fph 3 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
precise <- end_time - start_time
# user  system elapsed
# 0.31    0.16  633.53


### Tests with custom initialised values

# niter = 100, tol = 0.05
start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210602_hmt_100p_low_custominit -f 1024 -numPerm 100 -fph 3 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
lowprec_custominit <- end_time - start_time
# user  system elapsed
# 0.43    0.13   13.06

start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210602_hmt_100p_custominit -f 1024 -numPerm 100 -fph 3 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
custominit <- end_time - start_time
# user  system elapsed
# 0.28    0.20  211.17

# (niter = 1000, tol = 0.05)
start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210602_hmt_100p_low2_custominit -f 1024 -numPerm 100 -fph 3 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
low2_custominit <- end_time - start_time
# user  system elapsed
# 0.11    0.14   10.77

# Histogram of permuted logLRs
curr_perms <- read.table("output/20210601_hmt_curr.fph.perm.logLR.txt")$V1
less_iter_perms <- read.table("output/20210601_hmt_less_iter.fph.perm.logLR.txt")$V1
more_tol_perms <- read.table("output/20210601_hmt_more_tol.fph.perm.logLR.txt")$V1
both_perms <- read.table("output/20210601_hmt_both.fph.perm.logLR.txt")$V1
precise_perms <- read.table("output/20210601_hmt_precise.fph.perm.logLR.txt")$V1
lowprec_custominit_perms <- read.table("output/20210602_hmt_100p_low_custominit.fph.perm.logLR.txt")$V1
custominit_perms <- read.table("output/20210602_hmt_100p_custominit.fph.perm.logLR.txt")$V1
lowprec2_custominit_perms <- read.table("output/20210602_hmt_100p_low2_custominit.fph.perm.logLR.txt")$V1

summary(curr_perms)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1908  1.2325  2.0376  2.3601  2.9491  8.7890
summary(less_iter_perms)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0865  1.1638  1.9905  2.3005  2.8988  8.7522
summary(more_tol_perms)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -0.09127  1.04524  1.86923  2.17174  2.83195  8.64260
summary(both_perms)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -0.09127  1.04524  1.86923  2.17174  2.83195  8.64260
summary(precise_perms)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2799  1.3092  2.0858  2.4247  3.0233  8.8713
summary(lowprec_custominit_perms)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -0.113   0.543   1.381   1.726   2.386   8.005
summary(custominit_perms)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2552  1.2472  2.0494  2.3747  2.9711  8.8042
summary(lowprec2_custominit_perms)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -0.113   0.543   1.381   1.726   2.386   8.005

hist(curr_perms)
hist(less_iter_perms)
hist(more_tol_perms)
hist(both_perms)
hist(precise_perms)
hist(lowprec_custominit_perms)
hist(custominit_perms)
hist(lowprec2_custominit_perms)

perms_df <- data.frame(perms = c(curr_perms,less_iter_perms,more_tol_perms,both_perms,precise_perms,lowprec_custominit_perms,custominit_perms,lowprec2_custominit_perms))
perms_df$label = c(rep("curr",100)
                   ,rep("less_iter",100)
                   ,rep("more_tol",100)
                   ,rep("both",100)
                   ,rep("precise",100)
                   ,rep("lowprec_custominit",100)
                   ,rep("custominit",100)
                   ,rep("lowprec2_custominit",100))
ggplot(perms_df, aes(x = perms, fill = label)) +
  geom_histogram()


# Tests on non-permuted data ----------------------------------------------

## How does the precision of the EM actually change if we change the params?
# This is for the first run.

# Current run (niter = 1000, tol = 0.0005)
start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210601_hmt_norm_prec -f 1024 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
norm_prec <- end_time - start_time
# user  system elapsed
# 0.00    0.00    2.55

# (niter = 10000, tol = 0.00005)
start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210601_hmt_high_prec -f 1024 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
high_prec <- end_time - start_time
# user  system elapsed
# 0.00    0.01    7.14

# (niter = 100, tol = 0.005)
start_time <- proc.time()
system(command = "../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210601_hmt_low_prec -f 1024 -hmt 1"
       ,show.output.on.console = T)
end_time <- proc.time()
low_prec <- end_time - start_time
# user  system elapsed
# 0.02    0.02    1.33

## logLR of each SNP under the three methods
log_norm_prec <- read.table("output/20210601_hmt_norm_prec.fph.logLR.txt")
log_high_prec <- read.table("output/20210601_hmt_high_prec.fph.logLR.txt")
log_low_prec <- read.table("output/20210601_hmt_low_prec.fph.logLR.txt")

log_norm_prec[,1:2]
log_high_prec[,1:2]
log_low_prec[,1:2]

log_values <- log_norm_prec[,1:2]
# names(log_values) <- c("snp","norm_prec")
log_values <- cbind(log_values
                    ,log_high_prec[,2])
log_values <- cbind(log_values
                    ,log_low_prec[,2])
names(log_values) <- c("snp","norm_prec","high_prec","low_prec")

ggplot(melt(log_values
            ,id.vars = "snp")
       ,aes(x = snp, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(hjust = 1, angle = 60))

# Scaled as % diff vs normal...
log_values$norm_prec_idx = ((log_values$norm_prec - log_values$norm_prec)/log_values$norm_prec)
log_values$high_prec_idx = ((log_values$high_prec - log_values$norm_prec)/log_values$norm_prec)
log_values$low_prec_idx = ((log_values$low_prec - log_values$norm_prec)/log_values$norm_prec)

ggplot(melt(log_values[,c(1,5:7)]
            ,id.vars = "snp")
       ,aes(x = snp, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(hjust = 1, angle = 60))

ggplot(melt(log_values[-24,c(1,5:7)]
            ,id.vars = "snp")
       ,aes(x = snp, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(hjust = 1, angle = 60))


# Tests on 10k perms ------------------------------------------------------

# 12:37 AM - 1:50 AM (~ 1h, 12m)
# (niter = 100, tol = 0.005)
"../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210602_hmt_10k_p -f 1024 -numPerm 10000 -fph 3 -hmt 1"

# 12:06 PM - 12:1 PM (~ 6m)
# (niter = 1000, tol = 0.05) AND custom initialisations
"../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o 20210602_hmt_10kp_low2_custominit -f 1024 -numPerm 10000 -fph 3 -hmt 1"

lowprec_10k <- read.table("output/20210602_hmt_10k_p.fph.perm.logLR.txt")$V1
lowprec2_custominit_10k <- read.table("output/20210602_hmt_10kp_low2_custominit.fph.perm.logLR.txt")$V1
summary(lowprec_10k)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -0.1282  1.1581  1.8739  2.2832  2.8985 24.3197
summary(lowprec2_custominit_10k)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -0.1417  0.5680  1.3768  1.7669  2.3939 23.7961

perms_10k_df <- data.frame(perms = c(lowprec_10k,lowprec2_custominit_10k))
perms_10k_df$label = c(rep("lowprec_10k",10000)
                   ,rep("lowprec2_custominit_10k",10000))
ggplot(perms_10k_df, aes(x = perms, fill = label, col = label, alpha = 0.1)) +
  geom_histogram(position = "identity")
