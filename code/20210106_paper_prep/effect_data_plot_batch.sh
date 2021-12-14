#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename


#SBATCH --job-name=effect_data_plots                # A single job name for the array
#SBATCH --cpus-per-task=1                       # Number of cores
#SBATCH --ntasks=1                       # All cores on one machine
#SBATCH --mem 80000                 # Memory request (4Gb)
#SBATCH -t 0-24:00 #Maximum execution time (D-HH:MM)

#SBATCH -p mig                     # Partition
#SBATCH --account=punim0614

#SBATCH --mail-user=bklaw@student.unimelb.edu.au   # email address
#SBATCH --mail-type=END         # only send email if END
#SBATCH --mail-type=FAIL        # only send email if FAIL
#SBATCH -D /data/gpfs/projects/punim0614/brendan/WaveQTL_HMT_wperm/20210918_run5_scale_coeff

#SBATCH -o effect_data_plots.out        # Standard output
#SBATCH -e effect_data_plots.err        # Standard error

# LOAD MODULES
module load 'r/3.6.0'

# INSERT CODE
Rscript plotting_funcs_two_plots.R