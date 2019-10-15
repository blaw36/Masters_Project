# NULL - WaveQTL
# setwd("~/actual_sims_output/null/")
# for(task_id in 1:2){
  for(task_id in 1:578){
  print(paste0("Running segment ",task_id))
  system(paste0("~/WaveQTL_HMT/WaveQTL"
         ," -gmode 1"
         ," -g /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno4.txt"
         ," -p /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v1/null/wave/fullread.4ind.data/DNase.",task_id,".txt"
         ," -group ~/WaveQTL/WaveQTL-master/test/dsQTL/group.txt"
         ," -u /data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v1/null/wave/fullread.4ind.data/use.",task_id,".txt"
         ," -o fullread.4ind.",task_id
         ," -f 1024 -fph 1"
         ," > out.",task_id," 2> err.",task_id))
}
