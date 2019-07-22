## WaveQTL_HMT running calls:

- All 24 SNPs on all 1024 phenotypes, using the pre-processing steps as per the WaveQTL software manual.
- Will use its "use/do not use" vector for WCs
- Will use quantile-transformed for likelihood ratio testing
- Will use non-quantile transformed for measuring effect size in data space.

### Tying
Several tying techniques:

### With quantile transform (for likelihood ratio tests)
#### Group 1: Top (root) of tree on its own level

1) Tree-level tying (uses default group file, which is based on tree levels)
```
../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o tree_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

2) Tree levels 2-3 tied together
```
../../WaveQTL -gmode 1 -group g23_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g23_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

3) Tree levels 2-4 tied together
```
../../WaveQTL -gmode 1 -group g24_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g24_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

4) Tree levels 2-5 tied together
```
../../WaveQTL -gmode 1 -group g25_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g25_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

#### Group 2: Top (root) of tree grouped with the first tying level

5) Tree levels 1-3 tied together
```
../../WaveQTL -gmode 1 -group g13_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g13_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

6) Tree levels 1-4 tied together
```
../../WaveQTL -gmode 1 -group g14_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g14_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

7) Tree levels 1-5 tied together
```
../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g15_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

### Without quantile transform (for effect size measurement)
1) Tree-level tying (uses default group file, which is based on tree levels)
```
../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.no.QT.txt -u use.txt -o tree_tie_noQT -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

### Deepdiving for simulations
../../WaveQTL -gmode 1 -group DeepdiveExamples/deepdive_grouping.txt -g DeepdiveExamples/g_07_03_good.cis.geno -p DeepdiveExamples/WCs_07_03_good.txt -u DeepdiveExamples/use_all.txt -o sim_07_03_good -f 128 -hmt 1 > DeepdiveExamples/sim_07_03_good_out.txt 2> DeepdiveExamples/sim_07_03_good_err.txt

../../WaveQTL -gmode 1 -group DeepdiveExamples/deepdive_grouping.txt -g DeepdiveExamples/g_07_03_bad.cis.geno -p DeepdiveExamples/WCs_07_03_bad.txt -u DeepdiveExamples/use_all.txt -o sim_07_03_bad -f 128 -hmt 1 > DeepdiveExamples/sim_07_03_bad_out.txt 2> DeepdiveExamples/sim_07_03_bad_err.txt

../../WaveQTL -gmode 1 -group DeepdiveExamples/deepdive_grouping.txt -g DeepdiveExamples/g_09_01_good_1.cis.geno -p DeepdiveExamples/WCs_09_01_good_1.txt -u DeepdiveExamples/use_all.txt -o sim_09_01_good_1 -f 128 -hmt 1 > DeepdiveExamples/sim_09_01_good_1_out.txt 2> DeepdiveExamples/sim_09_01_good_1_err.txt

../../WaveQTL -gmode 1 -group DeepdiveExamples/deepdive_grouping.txt -g DeepdiveExamples/g_09_01_good_2.cis.geno -p DeepdiveExamples/WCs_09_01_good_2.txt -u DeepdiveExamples/use_all.txt -o sim_09_01_good_2 -f 128 -hmt 1 > DeepdiveExamples/sim_09_01_good_2_out.txt 2> DeepdiveExamples/sim_09_01_good_2_err.txt

../../WaveQTL -gmode 1 -group DeepdiveExamples/deepdive_grouping.txt -g DeepdiveExamples/g_09_01_bad_1.cis.geno -p DeepdiveExamples/WCs_09_01_bad_1.txt -u DeepdiveExamples/use_all.txt -o sim_09_01_bad_1 -f 128 -hmt 1 > DeepdiveExamples/sim_09_01_bad_1_out.txt 2> DeepdiveExamples/sim_09_01_bad_1_err.txt

../../WaveQTL -gmode 1 -group DeepdiveExamples/deepdive_grouping.txt -g DeepdiveExamples/g_09_01_bad_2.cis.geno -p DeepdiveExamples/WCs_09_01_bad_2.txt -u DeepdiveExamples/use_all.txt -o sim_09_01_bad_2 -f 128 -hmt 1 > DeepdiveExamples/sim_09_01_bad_2_out.txt 2> DeepdiveExamples/sim_09_01_bad_2_err.txt

../../WaveQTL -gmode 1 -group DeepdiveExamples/deepdive_grouping.txt -g DeepdiveExamples/g_09_01_bad_nan.cis.geno -p DeepdiveExamples/WCs_09_01_bad_nan.txt -u DeepdiveExamples/use_all.txt -o sim_09_01_bad_nan -f 128 -hmt 1 > DeepdiveExamples/sim_09_01_bad_nan_out.txt 2> DeepdiveExamples/sim_09_01_bad_nan_err.txt

../../WaveQTL -gmode 1 -group DeepdiveExamples/deepdive_grouping.txt -g DeepdiveExamples/beta_sens_07_03_b1.cis.geno -p DeepdiveExamples/beta_sens_07_03_b1.txt -u DeepdiveExamples/use_all.txt -o beta_sens_07_03_b1 -f 128 -hmt 1 > DeepdiveExamples/beta_sens_07_03_b1_out.txt 2> DeepdiveExamples/beta_sens_07_03_b1_nan_err.txt

../../WaveQTL -gmode 1 -group DeepdiveExamples/deepdive_grouping.txt -g DeepdiveExamples/beta_sens_07_03_b2.cis.geno -p DeepdiveExamples/beta_sens_07_03_b2.txt -u DeepdiveExamples/use_all.txt -o beta_sens_07_03_b2 -f 128 -hmt 1 > DeepdiveExamples/beta_sens_07_03_b2_out.txt 2> DeepdiveExamples/beta_sens_07_03_b2_nan_err.txt

<!-- TOY EXAMPLE -->
<!-- WITH HMT -->
<!-- No tying -->
../../WaveQTL -gmode 1 -g DeepdiveExamples/toy_eg1_16wc.cis.geno -p DeepdiveExamples/toy_eg1_16wc.txt -u use_all_16.txt -o toy_eg1_16wc_noTie_hmt -f 16 -hmt 1 > DeepdiveExamples/toy_eg1_16wc_out.txt 2> DeepdiveExamples/toy_eg1_16wc_err.txt
<!-- Tie all -->
../../WaveQTL -gmode 1 -group tree_grp_8_g2.txt -g DeepdiveExamples/toy_eg1_16wc.cis.geno -p DeepdiveExamples/toy_eg1_16wc.txt -u use_all_16.txt -o toy_eg1_16wc_wTie_hmt -f 16 -hmt 1 > DeepdiveExamples/toy_eg1_16wc_out.txt 2> DeepdiveExamples/toy_eg1_16wc_err.txt
<!-- WITHOUT HMT -->
<!-- No tying -->
../../WaveQTL -gmode 1 -g DeepdiveExamples/toy_eg1_16wc.cis.geno -p DeepdiveExamples/toy_eg1_16wc.txt -u use_all_16.txt -o toy_eg1_16wc_noTie_noHmt -f 16 -fph 1 > DeepdiveExamples/toy_eg1_16wc_out.txt 2> DeepdiveExamples/toy_eg1_16wc_err.txt
<!-- Tie all -->
../../WaveQTL -gmode 1 -group tree_grp_8_g2.txt -g DeepdiveExamples/toy_eg1_16wc.cis.geno -p DeepdiveExamples/toy_eg1_16wc.txt -u use_all_16.txt -o toy_eg1_16wc_wTie_noHmt -f 16 -fph 1 > DeepdiveExamples/toy_eg1_16wc_out.txt 2> DeepdiveExamples/toy_eg1_16wc_err.txt

<!-- 32-tree toy -->
../../WaveQTL -gmode 1 -group tree_grp_8_g2.txt -g DeepdiveExamples/toy_eg1_16wc.cis.geno -p DeepdiveExamples/toy_eg1_32wc.txt -u use_all_32.txt -o toy_eg1_32wc_wTie_hmt -f 32 -hmt 1 > DeepdiveExamples/toy_eg1_32wc_out.txt 2> DeepdiveExamples/toy_eg1_32wc_err.txt
<!-- 1024-tree toy -->
../../WaveQTL -gmode 1 -group tree_grp_8_g2.txt -g DeepdiveExamples/toy_eg1_1024wc.cis.geno -p DeepdiveExamples/toy_eg1_1024wc.txt -u use_all_1024.txt -o toy_eg1_1024wc_wTie_hmt -f 1024 -hmt 1 > DeepdiveExamples/toy_eg1_1024wc_out.txt 2> DeepdiveExamples/toy_eg1_1024wc_err.txt