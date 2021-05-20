### ALT


# ALT - WaveQTL
# setwd("~/actual_sims_output/alt/")
#for(task_id in 1:2){
  for(task_id in 1:578){
  print(paste0("Running segment ",task_id))
  system(paste0("~/WaveQTL_HMT/WaveQTL"
         ," -gmode 1"
         ," -g /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno70.txt"
         ," -p /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/alt/wave/fullread.70ind.data/DNase.",task_id,".txt"
         ," -group ~/WaveQTL/WaveQTL-master/test/dsQTL/group.txt"
         ," -u /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/alt/wave/fullread.70ind.data/use.",task_id,".txt"
         ," -o fullread.70ind.",task_id
         ," -f 1024 -fph 1"
         ," > out.",task_id," 2> err.",task_id))
}

# ALT R4 - WaveQTL
# setwd("~/actual_sims_output/alt/")
# for(task_id in 1:2){
  for(task_id in 1:578){
  print(paste0("Running segment ",task_id))
  system(paste0("~/WaveQTL_HMT/WaveQTL"
         ," -gmode 1"
         ," -g /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno4.txt"
         ," -p /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/alt/wave/fullread.4ind.data/DNase.",task_id,".txt"
         ," -group ~/WaveQTL/WaveQTL-master/test/dsQTL/group.txt"
         ," -u /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/alt/wave/fullread.4ind.data/use.",task_id,".txt"
         ," -o fullread.4ind.",task_id
         ," -f 1024 -fph 1"
         ," > out.",task_id," 2> err.",task_id))
}


### ALT HMT
# ALT - WaveQTL-HMT
setwd("~/actual_sims_output/alt-hmt/")
#for(task_id in 1:2){
  for(task_id in 1:578){
  print(paste0("Running segment ",task_id))
  system(paste0("~/WaveQTL_HMT/WaveQTL"
         ," -gmode 1"
         ," -g /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno70.txt"
         ," -p /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/alt/wave/fullread.70ind.data/DNase.",task_id,".txt"
         ," -group ~/WaveQTL_HMT/test/dsQTL/g15_1024.txt"
         ," -u /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/alt/wave/fullread.70ind.data/use.",task_id,".txt"
         ," -o fullread.70ind.",task_id
         ," -f 1024 -hmt 1"
         ," > out.",task_id," 2> err.",task_id))
}


# ALT R4 - WaveQTL-HMT
#setwd("~/actual_sims_output/alt-hmt/")
#for(task_id in 1:2){
  for(task_id in 1:578){
  print(paste0("Running segment ",task_id))
  system(paste0("~/WaveQTL_HMT/WaveQTL"
         ," -gmode 1"
         ," -g /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno4.txt"
         ," -p /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/alt/wave/fullread.4ind.data/DNase.",task_id,".txt"
         ," -group ~/WaveQTL_HMT/test/dsQTL/g15_1024.txt"
         ," -u /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/alt/wave/fullread.4ind.data/use.",task_id,".txt"
         ," -o fullread.4ind.",task_id
         ," -f 1024 -hmt 1"
         ," > out.",task_id," 2> err.",task_id))
}

../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o test_1kperm_hmt -f 1024 -numPerm 1000 -fph 3 -hmt 1
../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o test_1kperm_no_hmt -f 1024 -numPerm 1000 -fph 3 -hmt 1