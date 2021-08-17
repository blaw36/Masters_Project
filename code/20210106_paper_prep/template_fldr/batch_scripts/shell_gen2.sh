#!/bin/bash

orig_dir=`pwd`

cd ..
home_dir=`pwd`

cd $orig_dir

for i in `seq 1 100`
do

 let grp1=(i*2)-1
 let grp2=(i*2)

 touch grp_${grp1}_${grp2}.sh
 
cat << EOF > grp_${grp1}_${grp2}.sh
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename


#SBATCH --job-name=g_${grp1}_${grp2}                # A single job name for the array
#SBATCH --cpus-per-task=1                       # Number of cores
#SBATCH --ntasks=1                       # All cores on one machine
#SBATCH --mem 80000                 # Memory request (4Gb)
#SBATCH -t 0-24:00 #Maximum execution time (D-HH:MM)

#SBATCH -p mig                     # Partition
#SBATCH --account=punim0614

#SBATCH --mail-user=bklaw@student.unimelb.edu.au   # email address
#SBATCH --mail-type=END         # only send email if END
#SBATCH --mail-type=FAIL        # only send email if FAIL
#SBATCH -D $home_dir 

#SBATCH -o logs/grp_${grp1}_${grp2}.out        # Standard output
#SBATCH -e logs/grp_${grp1}_${grp2}.err        # Standard error

# LOAD MODULES
module load 'r/3.6.0'

# INSERT CODE
Rscript run_spartan.R $grp1
Rscript run_spartan.R $grp2
EOF

done
