file_prefix <- "toy_eg1_16wc_noTie_hmt"
logLR <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.logLR.txt"))[3:18]))
tree_plot(logLR,yaxis_lims = c(min(logLR),max(logLR)),plot_title = "logLR")
pi_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pi.txt"))[2]))
eps_11_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.eps_11.txt"))[2:17]))
eps_10_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.eps_10.txt"))[2:17]))
eps_11_file[c(1,2,3,5,9)]
eps_10_file[c(1,2,3,5,9)]
pp_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pp.txt"))[2:17]))
pp_joint_01_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pp_joint_01.txt"))[2:17]))
pp_joint_10_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pp_joint_10.txt"))[2:17]))
pp_joint_11_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pp_joint_11.txt"))[2:17]))
round(pp_file,5)
round(pp_joint_01_file,5)
round(pp_joint_10_file,5)
round(pp_joint_11_file,5)


file_prefix <- "toy_eg1_16wc_wTie_hmt"
pi_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pi.txt"))[2]))
eps_11_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.eps_11.txt"))[2:17]))
eps_10_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.eps_10.txt"))[2:17]))
eps_11_file[c(1,2,3)]
eps_10_file[c(1,2,3)]
pp_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pp.txt"))[2:17]))
pp_joint_01_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pp_joint_01.txt"))[2:17]))
pp_joint_10_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pp_joint_10.txt"))[2:17]))
pp_joint_11_file <- exp(as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pp_joint_11.txt"))[2:17]))

file_prefix <- "toy_eg1_16wc_noTie_noHmt"
pi_file <- as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pi.txt"))[2:6])

file_prefix <- "toy_eg1_16wc_wTie_noHmt"
pi_file <- as.numeric(read.table(paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/output/",file_prefix,".fph.pi.txt"))[2:4])
