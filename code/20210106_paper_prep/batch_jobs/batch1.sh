#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename


#SBATCH --job-name=data_prep1                # A single job name for the array
#SBATCH --cpus-per-task=1                       # Number of cores
#SBATCH --ntasks=1                       # All cores on one machine
#SBATCH --mem 80000                 # Memory request (4Gb)
#SBATCH -t 0-03:00                  # Maximum execution time (D-HH:MM)

#SBATCH -p mig                     # Partition
#SBATCH --account=punim0614

#SBATCH --mail-user=bklaw@student.unimelb.edu.au   # email address
#SBATCH --mail-type=END         # only send email if END
#SBATCH --mail-type=FAIL        # only send email if FAIL

#SBATCH -o /home/bklaw/paper_data_clean/batch_logs/test1.out        # Standard output
#SBATCH -e /home/bklaw/paper_data_clean/batch_logs/test1.err        # Standard error



# LOAD MODULES
module load 'r/3.6.0'

# INSERT CODE
Rscript /home/bklaw/paper_data_clean/01_prepare_data_for_simulation_SERVER_EDIT.R 1
