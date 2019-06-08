# wavelet checker

file_prefix = "tie_test1"
path = "Cpp/WaveQTL_HMT/test/dsQTL/output/"

files_to_chk <- list.files(path)
files_to_chk <- grep(paste0("^",file_prefix,".") ,files_to_chk, value = T)

cpp_output <- list()
for(i in files_to_chk){
  tmp_name = stringr::str_split_fixed(string = i, file_prefix, n = 2)[2]
  if(grepl("log.txt",tmp_name)){
    cpp_output[[tmp_name]] = readLines(paste0(path,i))
  }else{
    cpp_output[[tmp_name]] = as.matrix(read.table(paste0(path,i)))
  }
}

if(sumlog_use){
  all.equal(round(a_i[,`1`],5),as.numeric(cpp_output$.fph.alpha_sl_1.txt[-(1:2)]))
  all.equal(round(a_i[,`0`],5),as.numeric(cpp_output$.fph.alpha_sl_0.txt[-(1:2)]))

  all.equal(round(b_i[,`1`],5),as.numeric(cpp_output$.fph.beta_sl_1.txt[-(1:2)]))
  all.equal(round(b_i[,`0`],5),as.numeric(cpp_output$.fph.beta_sl_0.txt[-(1:2)]))

  all.equal(round(b_i_pi[-1,`1`],5),as.numeric(cpp_output$.fph.beta_sl_psl_1.txt[-(1:3)]))
  all.equal(round(b_i_pi[-1,`0`],5),as.numeric(cpp_output$.fph.beta_sl_psl_0.txt[-(1:3)]))

  all.equal(round(b_pi_no_i[-1,`1`],5),as.numeric(cpp_output$.fph.beta_psl_no_sl_1.txt[-(1:3)]))
  all.equal(round(b_pi_no_i[-1,`0`],5),as.numeric(cpp_output$.fph.beta_psl_no_sl_0.txt[-(1:3)]))

  all.equal(round(pp_i[,`1`],5),as.numeric(cpp_output$.fph.pp.txt[-(1:2)]))
  all.equal(round(pp_j_i[-1,`11`],5),as.numeric(cpp_output$.fph.pp_joint_11.txt[-(1:3)]))
  all.equal(round(pp_j_i[-1,`10`],5),as.numeric(cpp_output$.fph.pp_joint_10.txt[-(1:3)]))
  all.equal(round(pp_j_i[-1,`01`],5),as.numeric(cpp_output$.fph.pp_joint_01.txt[-(1:3)]))

  all.equal(round(pi,5),as.numeric(cpp_output$.fph.pi.txt[-(1)]))
  all.equal(round(eps[-1,`11`],5),as.numeric(cpp_output$.fph.eps_11.txt[-(1:3)]))
  all.equal(round(eps[-1,`10`],5),as.numeric(cpp_output$.fph.eps_10.txt[-(1:3)]))

  # Final logL, plus the log_10(BFs) for each other
  all.equal(round(c(new_logL,bfs/log(10)),5),as.numeric(cpp_output$.fph.logLR.txt[-c(1,3)]))

}else{
  all.equal(round(exp(a_i[,`1`]),5),as.numeric(cpp_output$.fph.alpha_sl_1.txt[-(1:2)]))
  all.equal(round(exp(a_i[,`0`]),5),as.numeric(cpp_output$.fph.alpha_sl_0.txt[-(1:2)]))

  all.equal(round(exp(b_i[,`1`]),5),as.numeric(cpp_output$.fph.beta_sl_1.txt[-(1:2)]))
  all.equal(round(exp(b_i[,`0`]),5),as.numeric(cpp_output$.fph.beta_sl_0.txt[-(1:2)]))

  all.equal(round(exp(b_i_pi[-1,`1`]),5),as.numeric(cpp_output$.fph.beta_sl_psl_1.txt[-(1:3)]))
  all.equal(round(exp(b_i_pi[-1,`0`]),5),as.numeric(cpp_output$.fph.beta_sl_psl_0.txt[-(1:3)]))

  all.equal(round(exp(b_pi_no_i[-1,`1`]),5),as.numeric(cpp_output$.fph.beta_psl_no_sl_1.txt[-(1:3)]))
  all.equal(round(exp(b_pi_no_i[-1,`0`]),5),as.numeric(cpp_output$.fph.beta_psl_no_sl_0.txt[-(1:3)]))

  all.equal(round(exp(pp_i[,`1`]),5),as.numeric(cpp_output$.fph.pp.txt[-(1:2)]))
  all.equal(round(exp(pp_j_i[-1,`11`]),5),as.numeric(cpp_output$.fph.pp_joint_11.txt[-(1:3)]))
  all.equal(round(exp(pp_j_i[-1,`10`]),5),as.numeric(cpp_output$.fph.pp_joint_10.txt[-(1:3)]))
  all.equal(round(exp(pp_j_i[-1,`01`]),5),as.numeric(cpp_output$.fph.pp_joint_01.txt[-(1:3)]))

  all.equal(round(pi,5),as.numeric(cpp_output$.fph.pi.txt[-(1)]))
  all.equal(round(eps[-1,`11`],5),as.numeric(cpp_output$.fph.eps_11.txt[-(1:3)]))
  all.equal(round(eps[-1,`10`],5),as.numeric(cpp_output$.fph.eps_10.txt[-(1:3)]))

  # Final logL, plus the log_10(BFs) for each other
  all.equal(round(c(new_logL,bfs/log(10)),5),as.numeric(cpp_output$.fph.logLR.txt[-c(1,3)]))
}
