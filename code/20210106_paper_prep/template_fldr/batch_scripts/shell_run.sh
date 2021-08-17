#!/bin/bash

for i in `seq 1 100`
do

 let grp1=(i*2)-1
 let grp2=(i*2)
 
 sbatch grp_${grp1}_${grp2}.sh

done
