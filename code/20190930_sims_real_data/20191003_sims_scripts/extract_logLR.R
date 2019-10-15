null_lhoods <- c()
null_hmt_lhoods <- c()
alt_lhoods <- c()
alt_hmt_lhoods <- c()

output_path <- "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v4/"

for(i in 1:578){
  # Read in null l-hood file
  null_lhoods_file <- as.numeric(read.table(paste0(output_path,"null/output/fullread.70ind.",i,".fph.logLR.txt"))[2:3])

  # NULL
  null_lhoods[i] <- null_lhoods_file[1]

  # NULL-HMT
  # Grab scaling coefficient logBF, plus scaling coefficient pi. (from non-HMT analysis - will remain unchanged in HMT)
  sc_pi_null <- as.numeric(read.table(paste0(output_path,"null/output/fullread.70ind.",i,".fph.pi.txt"))[2])
  # Convert out of log10
  sc_BF_null <- 10^(null_lhoods_file[2])
  # Calc lhood and convert back to log10
  sc_logL_null <- log(sc_BF_null*sc_pi_null + (1-sc_pi_null),base = 10)
  null_hmt_lhoods[i] <- as.numeric(read.table(paste0(output_path,"null/output/fullread.70ind.",i,".fph.logLR.txt"))[2]) + sc_logL_null
}

for(i in 1:578){
  # Read in alt l-hood file
  alt_lhoods_file <- as.numeric(read.table(paste0(output_path,"alt/output/fullread.70ind.",i,".fph.logLR.txt"))[2:3])

  # alt
  alt_lhoods[i] <- alt_lhoods_file[1]

  # alt-HMT
  # Grab scaling coefficient logBF, plus scaling coefficient pi. (from non-HMT analysis - will remain unchanged in HMT)
  sc_pi_alt <- as.numeric(read.table(paste0(output_path,"alt/output/fullread.70ind.",i,".fph.pi.txt"))[2])
  # Convert out of log10
  sc_BF_alt <- 10^(alt_lhoods_file[2])
  # Calc lhood and convert back to log10
  sc_logL_alt <- log(sc_BF_alt*sc_pi_alt + (1-sc_pi_alt),base = 10)
  alt_hmt_lhoods[i] <- as.numeric(read.table(paste0("/home/bklaw/actual_sims_output/alt-hmt/output/fullread.70ind.",i,".fph.logLR.txt"))[2]) + sc_logL_alt
}

all_lhoods <- c(null_lhoods
                ,null_hmt_lhoods
                ,alt_lhoods
                ,alt_hmt_lhoods)

today_dt <- format(Sys.Date(),"%Y%m%d")
saveRDS(all_lhoods,paste0("/home/bklaw/actual_sims_output/",today_dt,"logLR_results.RDS"),compress = T)
