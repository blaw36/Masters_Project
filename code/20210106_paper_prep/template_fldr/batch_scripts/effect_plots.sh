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
#SBATCH -D /data/gpfs/projects/punim0614/brendan/WaveQTL_HMT_wperm/20210918_run5_scale_coeff

#SBATCH -o effect_plots.out        # Standard output
#SBATCH -e effect_plots.err        # Standard error

# LOAD MODULES
module load 'r/3.6.0'

# INSERT CODE
Rscript effect_size_plot_spartan.R 1 23
Rscript effect_size_plot_spartan.R 1 111
Rscript effect_size_plot_spartan.R 1 453
Rscript effect_size_plot_spartan.R 1 1396
Rscript effect_size_plot_spartan.R 1 1559
Rscript effect_size_plot_spartan.R 1 1914
Rscript effect_size_plot_spartan.R 1 2109
Rscript effect_size_plot_spartan.R 1 2391
Rscript effect_size_plot_spartan.R 1 2562
Rscript effect_size_plot_spartan.R 1 3213
Rscript effect_size_plot_spartan.R 1 3352
Rscript effect_size_plot_spartan.R 1 3506
Rscript effect_size_plot_spartan.R 1 4063
Rscript effect_size_plot_spartan.R 1 4280
Rscript effect_size_plot_spartan.R 1 4705
Rscript effect_size_plot_spartan.R 2 268
Rscript effect_size_plot_spartan.R 2 694
Rscript effect_size_plot_spartan.R 2 808
Rscript effect_size_plot_spartan.R 2 1066
Rscript effect_size_plot_spartan.R 2 1333
Rscript effect_size_plot_spartan.R 2 1617
Rscript effect_size_plot_spartan.R 2 1736
Rscript effect_size_plot_spartan.R 2 1791
Rscript effect_size_plot_spartan.R 2 1918
Rscript effect_size_plot_spartan.R 2 2105
Rscript effect_size_plot_spartan.R 2 2376
Rscript effect_size_plot_spartan.R 2 2623
Rscript effect_size_plot_spartan.R 2 2796
Rscript effect_size_plot_spartan.R 2 2941
Rscript effect_size_plot_spartan.R 2 3082
Rscript effect_size_plot_spartan.R 2 3327
Rscript effect_size_plot_spartan.R 3 113
Rscript effect_size_plot_spartan.R 3 1159
Rscript effect_size_plot_spartan.R 3 1327
Rscript effect_size_plot_spartan.R 3 1579
Rscript effect_size_plot_spartan.R 3 2242
Rscript effect_size_plot_spartan.R 3 2398
Rscript effect_size_plot_spartan.R 4 10
Rscript effect_size_plot_spartan.R 4 587
Rscript effect_size_plot_spartan.R 4 792
Rscript effect_size_plot_spartan.R 4 1294
Rscript effect_size_plot_spartan.R 5 138
Rscript effect_size_plot_spartan.R 5 206
Rscript effect_size_plot_spartan.R 5 247
Rscript effect_size_plot_spartan.R 5 540
Rscript effect_size_plot_spartan.R 5 686
Rscript effect_size_plot_spartan.R 5 689
Rscript effect_size_plot_spartan.R 5 900
Rscript effect_size_plot_spartan.R 5 994
Rscript effect_size_plot_spartan.R 5 1117
Rscript effect_size_plot_spartan.R 5 1397
Rscript effect_size_plot_spartan.R 5 1952
Rscript effect_size_plot_spartan.R 5 2023
Rscript effect_size_plot_spartan.R 6 175
Rscript effect_size_plot_spartan.R 6 458
Rscript effect_size_plot_spartan.R 6 662
Rscript effect_size_plot_spartan.R 6 897
Rscript effect_size_plot_spartan.R 6 1238
Rscript effect_size_plot_spartan.R 7 89
Rscript effect_size_plot_spartan.R 7 554
Rscript effect_size_plot_spartan.R 7 1684
Rscript effect_size_plot_spartan.R 7 1765
Rscript effect_size_plot_spartan.R 7 2217
Rscript effect_size_plot_spartan.R 8 37
Rscript effect_size_plot_spartan.R 8 145
Rscript effect_size_plot_spartan.R 8 179
Rscript effect_size_plot_spartan.R 8 302
Rscript effect_size_plot_spartan.R 8 520
Rscript effect_size_plot_spartan.R 8 935
Rscript effect_size_plot_spartan.R 8 965
Rscript effect_size_plot_spartan.R 8 1204
Rscript effect_size_plot_spartan.R 8 1565
Rscript effect_size_plot_spartan.R 8 1749
Rscript effect_size_plot_spartan.R 9 350
Rscript effect_size_plot_spartan.R 9 478
Rscript effect_size_plot_spartan.R 9 951
Rscript effect_size_plot_spartan.R 9 1131
Rscript effect_size_plot_spartan.R 9 1477
Rscript effect_size_plot_spartan.R 9 1568
Rscript effect_size_plot_spartan.R 9 1580
Rscript effect_size_plot_spartan.R 9 1648
Rscript effect_size_plot_spartan.R 9 1658
Rscript effect_size_plot_spartan.R 9 2453
Rscript effect_size_plot_spartan.R 10 18
Rscript effect_size_plot_spartan.R 10 239
Rscript effect_size_plot_spartan.R 10 392
Rscript effect_size_plot_spartan.R 10 456
Rscript effect_size_plot_spartan.R 10 1089
Rscript effect_size_plot_spartan.R 10 1161
Rscript effect_size_plot_spartan.R 10 1168
Rscript effect_size_plot_spartan.R 10 1221
Rscript effect_size_plot_spartan.R 10 1452
Rscript effect_size_plot_spartan.R 10 2081
Rscript effect_size_plot_spartan.R 11 231
Rscript effect_size_plot_spartan.R 11 908
Rscript effect_size_plot_spartan.R 11 1088
Rscript effect_size_plot_spartan.R 11 1448
Rscript effect_size_plot_spartan.R 11 1717
Rscript effect_size_plot_spartan.R 11 1813
Rscript effect_size_plot_spartan.R 11 2099
Rscript effect_size_plot_spartan.R 11 2391
Rscript effect_size_plot_spartan.R 12 42
Rscript effect_size_plot_spartan.R 12 67
Rscript effect_size_plot_spartan.R 12 252
Rscript effect_size_plot_spartan.R 12 326
Rscript effect_size_plot_spartan.R 12 327
Rscript effect_size_plot_spartan.R 12 1354
Rscript effect_size_plot_spartan.R 12 1504
Rscript effect_size_plot_spartan.R 12 2284
Rscript effect_size_plot_spartan.R 12 2368
Rscript effect_size_plot_spartan.R 13 210
Rscript effect_size_plot_spartan.R 13 452
Rscript effect_size_plot_spartan.R 13 498
Rscript effect_size_plot_spartan.R 13 746
Rscript effect_size_plot_spartan.R 13 792
Rscript effect_size_plot_spartan.R 13 797
Rscript effect_size_plot_spartan.R 14 10
Rscript effect_size_plot_spartan.R 14 488
Rscript effect_size_plot_spartan.R 14 694
Rscript effect_size_plot_spartan.R 14 810
Rscript effect_size_plot_spartan.R 14 1109
Rscript effect_size_plot_spartan.R 14 1180
Rscript effect_size_plot_spartan.R 14 1461
Rscript effect_size_plot_spartan.R 14 1584
Rscript effect_size_plot_spartan.R 15 213
Rscript effect_size_plot_spartan.R 15 856
Rscript effect_size_plot_spartan.R 15 1163
Rscript effect_size_plot_spartan.R 16 236
Rscript effect_size_plot_spartan.R 16 239
Rscript effect_size_plot_spartan.R 16 1985
Rscript effect_size_plot_spartan.R 16 2107
Rscript effect_size_plot_spartan.R 17 566
Rscript effect_size_plot_spartan.R 17 902
Rscript effect_size_plot_spartan.R 17 1221
Rscript effect_size_plot_spartan.R 17 1727
Rscript effect_size_plot_spartan.R 17 1972
Rscript effect_size_plot_spartan.R 17 2101
Rscript effect_size_plot_spartan.R 17 2131
Rscript effect_size_plot_spartan.R 17 2553
Rscript effect_size_plot_spartan.R 17 2757
Rscript effect_size_plot_spartan.R 17 2850
Rscript effect_size_plot_spartan.R 17 3026
Rscript effect_size_plot_spartan.R 17 3129
Rscript effect_size_plot_spartan.R 18 43
Rscript effect_size_plot_spartan.R 18 319
Rscript effect_size_plot_spartan.R 18 477
Rscript effect_size_plot_spartan.R 18 786
Rscript effect_size_plot_spartan.R 19 121
Rscript effect_size_plot_spartan.R 19 225
Rscript effect_size_plot_spartan.R 19 437
Rscript effect_size_plot_spartan.R 19 750
Rscript effect_size_plot_spartan.R 19 1060
Rscript effect_size_plot_spartan.R 19 2379
Rscript effect_size_plot_spartan.R 19 3114
Rscript effect_size_plot_spartan.R 19 3125
Rscript effect_size_plot_spartan.R 19 3322
Rscript effect_size_plot_spartan.R 20 375
Rscript effect_size_plot_spartan.R 20 537
Rscript effect_size_plot_spartan.R 20 932
Rscript effect_size_plot_spartan.R 20 1087
Rscript effect_size_plot_spartan.R 21 80
Rscript effect_size_plot_spartan.R 21 249
Rscript effect_size_plot_spartan.R 22 131
Rscript effect_size_plot_spartan.R 22 133
Rscript effect_size_plot_spartan.R 22 209
Rscript effect_size_plot_spartan.R 22 478
Rscript effect_size_plot_spartan.R 22 527
Rscript effect_size_plot_spartan.R 22 791
Rscript effect_size_plot_spartan.R 22 1386
Rscript effect_size_plot_spartan.R 22 1453
Rscript effect_size_plot_spartan.R 22 1530