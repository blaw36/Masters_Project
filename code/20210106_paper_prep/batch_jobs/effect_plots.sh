#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename


#SBATCH --job-name=g_1                # A single job name for the array
#SBATCH --cpus-per-task=1                       # Number of cores
#SBATCH --ntasks=1                       # All cores on one machine
#SBATCH --mem 80000                 # Memory request (4Gb)
#SBATCH -t 0-24:00 #Maximum execution time (D-HH:MM)

#SBATCH -p mig                     # Partition
#SBATCH --account=punim0614

#SBATCH --mail-user=bklaw@student.unimelb.edu.au   # email address
#SBATCH --mail-type=END         # only send email if END
#SBATCH --mail-type=FAIL        # only send email if FAIL
#SBATCH -D /data/gpfs/projects/punim0614/brendan/WaveQTL_HMT_wperm/20210627_run4_smller_rndm_batch 

#SBATCH -o effect_plots.out        # Standard output
#SBATCH -e effect_plots.err        # Standard error

# LOAD MODULES
module load 'r/3.6.0'

# INSERT CODE
Rscript effect_size_plot_spartan.R 1 2722
Rscript effect_size_plot_spartan.R 1 4596
Rscript effect_size_plot_spartan.R 1 825
Rscript effect_size_plot_spartan.R 2 696
Rscript effect_size_plot_spartan.R 20 1364
Rscript effect_size_plot_spartan.R 21 281
Rscript effect_size_plot_spartan.R 7 784
Rscript effect_size_plot_spartan.R 8 225
Rscript effect_size_plot_spartan.R 4 845
Rscript effect_size_plot_spartan.R 5 428
Rscript effect_size_plot_spartan.R 1 825
