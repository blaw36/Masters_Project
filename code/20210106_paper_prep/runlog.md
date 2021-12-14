## 2nd June 2021

Objective is to get the batch run going. 

### Current settings:_

__no HMT__
_non-permuted_
niter = 1000, tol = 0.0005

_permuted_
niter = 1000, tol = 0.05

__HMT__
_Non-permuted_
niter = 1000, tol = 0.0005
Initialisations
double pi_1 = log(0.5);
double pi_0 = log(0.5);

double eps_11 = log(0.5);
double eps_01 = log(0.5);
double eps_10 = log(0.5);
double eps_00 = log(0.5);

_Permuted_
niter = 1000, tol = 0.05
Initialisations
double pi_1 = log(0.1);
double pi_0 = log(0.9);

double eps_11 = log(0.5);
double eps_01 = log(0.5);
double eps_10 = log(0.1);
double eps_00 = log(0.9);

10,000 sites:
Run in chunks of 50, email Yupei about configs and allowable compute times.

High-level loop structure:
non-HMT
for each SNP...
	for each group...
		EM algorithm, until diff in logLR < tol
		break
	end
end

HMT
for each SNP...
	EM algorithm, until diff in logLR < tol
	break
end

We want to:
- run the R script which formats site data into waveqtl compatible data
- do runtime (see Shim 2014 4.5.2 for details) for each SITE (so we can subsample by site) for:
	- normalisation
	- WaveQTL
	- WaveQTL_HMT
- all nearby SNPs WITH permutations (10000 perms, fph -3)
- random sampling (into 5?) groups of 10,000 so we can run in 10,000's at a time until we're satisfied.
- be careful of null p-values/NA returns
- List of p-values -> FDR curve
- histogram of p-values for each method. 
- Under null should be uniform. 
- Under alternative should have a peak near 0 p-values
- If ok, then use q-value R-package to make the calculations

Implementing permuation tests in the CODE
- WaveQTL and WaveQTL_HMT should yield one p-val (should be just as small) with the SNP with the strongest association (should be the same)
- run for 10000 random sites

Things to worry about later:
	- can we change epsilon inits for each of the bottom groups?
	- how?
	- do not use default waveQTL groupings, use our own. Need to tune this.

Paths to (50k) sites:
"/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/region_01_sel_step1/"

Paths to (phenotype) data (70 x 1024): (only for 47024 sites)
paste0("/home/bklaw/paper_data_clean/analysis/simulation/sample_size/simulation_manydsQTL_v1/data/", chrIX, ".pheno.dat.", site)

Paths to (genotype/nearby SNP) data (xx x 70):
/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno/ (HAS MAF? is bigger)
/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/ (USE THIS ONE) (HAS NO MAF filter applied? is smaller, subsets of the above data)

## 14th June 2021
Created the file which splits 50k sites into 5 randomly selected groups of 10k sites. Done this from the ORIGINAL site list; note that not all 50k sites have data, only around 47024 have data to read in.

Paths to site splits: 
'/home/bklaw/paper_data_clean/site_splits'

Script is here: '~/paper_data_clean/site_splits.R'

## 16th June 2021
Changes made to src/Makefile on Spartan so it can compile. Specifically:

`CFLAGS += -I/user/local/include -L/user/local/lib` has become `CFLAGS += -L$(EBROOTGSL)/lib -I$(EBROOTGSL)/include` and
`-static` removed from `fp: fp.o $(OBJS); $(CC) -static $(CFLAGS) $(OBJS) fp.o $(STATICLIBS) $(LIBS) -o WaveQTL`  as the version of gcc on Spartan doesn't do static linking. This allows us to compile on Spartan with:
- `module load gcc/8.3.0`
- `module load gsl/2.5`

Managed to get the end to end R script onto spartan and get it working. HAve a template for running en masse now. Just need to adapt the parameters, really.

## 17th June 2021
Ran the script end-to-end (phenotype -> WC transform -> WaveQTL -> output)!
Need to do:
	- Scale it out to run over all 200 groups (of 50) [DONE]
	- Figure out how to run batch scripts from a separate folder but changing the working directory to start where we'd like it to be ('/home/bklaw/WaveQTL_HMT_wperm/20210616_test') [DONE]
	- Recompile with updated feedback from the HPC guys [DONE]
	- Tune parameters?
### Qs
- library read depth for each chr/site? yes, all the same [DONE]
- covariates are same for all chr/site? yes, all the same [DONE]
- use quantile transform [DONE]

To make it EVEN EASIER we should pass through the current run folder directory into the R script as a command line arg in the bash script! [TODO]

Results are stored in the folder '20210617_run1'

## 19th June 2021
Took the outputs and processed q-values as per HJ's scripts. Histograms and results don't look good at all. Let HJ know that we'll need to chat further.

Actually realised it's because I used the same genotype file for every chr-site. Changed R script to reflect this, and re-run. Changed the script to use site and chrIX as inputs so that the proper genotype data can be used.

## 20th June 2021
Looked at the results; still not great. I'm going to try and change the permuted settings to
`niter = 1000, tol = 0.005`, rather than `tol = 0.05`, because we can afford some run-time.

Findings in 20210620_run2

It didn't run properly; a lot of things stopped because of data issues. I think my home directory has run out of space, instead I need to save stuff here: `cd /data/gpfs/projects/punim0614`.

Check home directory: `check_home_usage`
Check project directory (see dir above): `check_project_usage`

Everything will now be in `/data/gpfs/projects/punim0614/brendan/WaveQTL_HMT_wperm`

Also decided to do a parallel run using WaveQTL (original, to verify) in `/data/gpfs/projects/punim0614/brendan/WaveQTL`

## 24th June 2021
There was some experimentation with the original WaveQTL to see what thresholds we'd need to adjust to get our p-val histogram to look like it should. The conclusion was that as long as the epsilon/thresholds were set the same for both permutation AND the original dataset, we get p-values looking pretty good. Runs were saved in `/data/gpfs/projects/punim0614/brendan/WaveQTL` and analysis results in `C:\Users\brend\Dropbox\Uni Stuff - Masters\Research Project\Masters_Project_Git\analysis\` Settings and folder names and corresponding results in git are:

- `20210620_verify_wqtl`; running original WaveQTL with original settings (epsilon = 0.0005 for both). Results in git at `\20210620_verify_wqtl`
- `20210622_run2_thold`; running original WaveQTL with epsilon = 0.005 for permutation. Results in git at `\20210622_wqtl_run2_param` 
- `20210624_run3_thold_both`; running original WaveQTL with epsilon = 0.005 for both datasets. Results in git at `\20210624_wqtl_run3_thold_both`

As a result, in `/data/gpfs/projects/punim0614/brendan/WaveQTL_HMT_wperm/20210624_run3`, we're going to try running on a WaveQTL whereby epsilon = tol = 0.005 for both original datasets and in permutation mode.

We've moved the old binary (with slightly amended settings as of 20th June) to a `bin/` folder, and called it `WaveQTL_20210620`

- p-value plots
- WaveQTL effect size diagram using the example cases

Good results; making the epsilon = tol = 0.005 for both the datasets helped out dramatically. Our p-value distributions look a lot nicer now.

## 4th July 2021

1) Find example where HMT is much better than no-HMT.
	- (1) scatter plot of -log10(pvalue) from HMT vs -log10(pvalue) from no-HMT 
	- (2) scatter plot of -log10(qvalue) from HMT vs -log10(qvalue) from no-HMT?
Example of (2) is Fig 2B, 2D here: https://github.com/heejungshim/multiseq_Shim_et_al/blob/main/manuscript/multiseq.pdf). Code: https://github.com/heejungshim/multiseq-ms-figures/blob/master/scripts/ATACseq_qvalues_multiseq_WaveQTL_DESeq2.Rmd

The idea is to find an example where HMT is much better than no-HMT (we will use an example in the paper; HMT's pvalue or qvalue much smaller) and another example where no-HMT is much better than HMT (we will use it to improve HMT because intuitively HMT should be better than no-HMT; no-HMT's pvalue or qvalue much smaller).

2) Find the example site where much better, and then run through effect size calculations.

Did 1), still need to do 2).

## 13th July 2021
Configured the 'logs_to_pval.R' script so that it now includes the chr and site number of each log file. Changed 'chr' field to 'strong_snp' being the strongest nearby SNP with signal.

Effect size plots to do:
- HMT strong and non-HMT weak: (1, 2722), (10,1880) (all seem to have similar results)
- non-HMT strong and HMT weak: 
	- all cases are because different nearby SNP identified
	- Multiple sites identified because potentially very different types of results
	- (1,4596), (1,825), (2,696), (20,1364), (21,281), (7,784), (8,225)
	- Note that (7,784) is the only situation where strong_snps are the same between HMT, non-HMT

					HMT_pval < Non-HMT_pval
        					FALSE TRUE
Same strong-SNP?  	FALSE	  981  614
  					TRUE	 2589 4643

A few others to check:
- Same sites:
	- HMT strong and non-HMT weak: (4,845)
- Different sites:
	- HMT strong and non-HMT weak: (5,428); HMT: chr5.7422839; noHMT: chr5.7423920
	- Non-HMT strong and HMT weak: (1,825); HMT: chr1.11991133; noHMT: chr1.indel.11988665

Began work on end-to-end script. Hope that it will input:
- chr
- site
- closest SNP name (need to match on the name)
- whether HMT or not
and output:
- the plot
- some of its constitutents (col_posi, beta_l, beta_r, beta_dataS)

The steps are:
1. Re-run with noQT, but also no permutations; just get the outputs as required. This will need us to replicate end_to_end but specially just to output the parts required for effect size plotting.
2. Read in mean and var outputs
3. Output as per get_effectSizeinDataSpace.R

Note that we'll need to Wmat matrices as per WaveQTL.

## 18th July 2021
Finished the scripts to do the effect size plots. Put them on spartan and ran some effect size plots in '/plots'. Sent some plots off to HJ. Findings:

- HMT strong and non-HMT weak seem to have 0 effect size in 2 of them (maybe representative of the 37) so not sure what 'signal' HMT has picked up
- HMT weak and non-HMT strong; largest differences are situations where both algorithms picked up different 'strongest SNPs'
- Situations where they've picked up the same 'strongest SNPs' show less difference, but still many instances where non-HMT pval < HMT pval

## 22nd July 2021
Plot out effect size - does it match with the data? Data figures are coloured based on the groups of SNPs/alleles-thingys related to each individual. Use the strongest SNP which we're creating the effect size plots for.
For category 2), we've identified the SNPs which SHOULD have signal from WaveQTL. So the key question here is why didn't HMT pick up that signal from that SNP that WaveQTL found? The things HMT found weren't meant to have signal anyway (weren't picked up by WaveQTL).

## 3rd August 2021
What i've done:
0) Re-run the plots under investigation on the same SNPs (WaveQTL)
	- Change effect_size_plot_funcs_spartan.R to have arguments for whether to use nearest snp from nohmt, hmt, or each run's nearest. Arguments will be 'nohmt', 'hmt', 'nearest'. Default will be 'nearest' (current behaviour).
	- Created a copy of `effect_size_plot_spartan.R` called `effect_size_plot_spartan_sameSNP.R`. This will call effect_size_plot with argument 'nohmt'
	- Will save in a folder called /plots_nohmt_snp. (Other option will save in a folder called /plots_hmt_snp).
1) Data plots
	- Data plots; in thesis, used cis.geno (not the noMAF), vs pheno.dat (70 individuals x counts at 1024 sites)
	- Grab the rounded SNP value for each individual, join to the pheno data, and then plot. (all as per thesis)
	- data_space_plot_spartan.R
	- use the data we ran with (ie. ref the directory)

Questions for HJ:
	- Data plots; in thesis, used cis.geno (not the noMAF) is that ok?

there is 'effect_plots.sh' or something like that in 'batch_scripts'. Consolidate the two batch scripts in the home dir (show examples of both) and then move it into templates folder.

Rscript --no-save --no-restore --verbose effect_size_plot_spartan_nohmt_snp.R 1 825 > effect_plots_nohmt_snp.Rout 2> effect_plots_nohmt_snp.Rerr
Rscript --no-save --no-restore --verbose data_space_plot_spartan.R 1 825 chr1.indel.11988665 > data_space_plot.Rout 2> data_space_plot.Rerr

nohmtLR = read.table("/data/gpfs/projects/punim0614/brendan/WaveQTL_HMT_wperm/20210627_run4_smller_rndm_batch/output/chr.1.2722.nohmt.fph.logLR.txt")
hmtLR = read.table("/data/gpfs/projects/punim0614/brendan/WaveQTL_HMT_wperm/20210627_run4_smller_rndm_batch/output/chr.1.2722.hmt.fph.logLR.txt")
hmtLR_perms = read.table("/data/gpfs/projects/punim0614/brendan/WaveQTL_HMT_wperm/20210627_run4_smller_rndm_batch/output/chr.1.2722.hmt.fph.perm.logLR.txt")



## 5th August 2021
Re-run the effect size plots with the data we ran our algo on (maf/nomaf?)

Issue 1: p-values calculated weirdly for this site. The permutations would suggest that the p-values should be high.
permuted likelihood ratios that leads to the p-values?
With chr1, site2722, why did the permutations calculate so weirdly?
[Print out logLR for nohmt permutations]

Issue 2: HMT doesn't do as well as noHMT; misses the bit in the middle
Chr 2, Site 696 (HMT poor and no-HMT good); why does HMT not do as well/misses the bit in the middle?
[logLR and logLR at different scales
Pi param estimates at different scales]

Issue 3: HMT poor, no-HMT good, but HMT should be 'stronger'. It has more 'pink regions' than no-HMT, but the p-value is higher.
Chr 20, Site 1364 (HMT poor and no-HMT good but should be HMT stronger?); HMT has more 'pink regions' than no-HMT but p-val is higher.
[effect size in original algo space (with quantile transform); Make effect size WITH quantile transform (ie. original data which we ran algo on) for everything]

## 17th August 2021
1) Check which pheno data we ran our algo on.
	In our WaveQTL runscripts, we used the following genotype data:
	`/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno"`
	In data size plot, we used:
	`as.matrix(read.table(paste0("/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno")))`
	They're the same, no further work required. Phew!

2) In order to debug properly, we should move data onto local computer:
(chr,site) properties
	- To run WaveQTL, we need:
		+ pheno_data_path = paste0("/home/bklaw/paper_data_clean/analysis/simulation/sample_size/simulation_manydsQTL_v1/data/chr", chrIX, ".pheno.dat.", site) [DONE]
		+ current version WaveQTL_HMT as per Spartan. [DONE]
		+ /data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno" [DONE]
		+ WCs (both QT and noQT) [DONE]
		+ use file [DONE]
	- To analyse WaveQTL, we need:
		- all the outputs [DONE]
	- To run the effect size plots, we need geno and pheno data, both from above

Local R script versions:
	- end_to_end_funcs_spartan_no_wc.R [DONE]
	- "/home/bklaw/WaveQTL_HMT_wperm/R/WaveQTL_preprocess_funcs.R" [DONE]
	- all the '.R' scripts in a folder [DONE]
	- "/home/bklaw/WaveQTL_HMT_wperm/R/end_to_end_funcs_spartan.R" [DONE]
	- "/home/bklaw/WaveQTL_HMT_wperm/R/effect_size_plot_funcs_spartan.R" [DONE]

## 18th August 2021
Outstanding work:
- Adjust WaveQTL_HMT so it prints out LogLR at each permutation with noHMT (.perm.logLR.txt) [DONE]
- effect size in original algo space (with quantile transform); Make effect size WITH quantile transform (ie. original data which we ran algo on) for everything [DONE]
- do logLR at different scales, Pi param estimates at different scales

Here's how to run locally using the current file structure (currently no `--group` option used). Run from `/test/dsQTL` folder
### HMT
`../../WaveQTL -gmode 1 -g ../../data/geno_data/chr1.825.geno -p chr.1.825_WCs.txt -u chr.1.825_use.txt -o chr.1.825.hmt -f 1024 -numPerm 10000 -fph 3 -hmt 1`

### No HMT
`../../WaveQTL -gmode 1 -g ../../data/geno_data/chr1.825.geno -p chr.1.825_WCs.txt -u chr.1.825_use.txt -o chr.1.825.nohmt -f 1024 -numPerm 10000 -fph 3`

For no QT, change the `-p`, `-u` options and probably `-o` to signal you've run something without quantile transformed inputs.

../../WaveQTL -gmode 1 -g ../../data/geno_data/chr1.825.geno -p chr.1.825_WCs.txt -u chr.1.825_use.txt -o chr.1.825.nohmt -f 1024 -numPerm 10000 -fph 3

__Debug issue 1) The p-values calculated weirdly for (1,2722)__
Issue here is that both display the same thing, but for HMT pval is really small (~0), and for non-HMT, it's really big ~1.
- I think there's clearly something going on here.
- noHMT has logLR of 0, and permutations of 0
- HMT has logLR of 0, and permutations of 4.75675e-13, which are all _just_ above 0, and therefore it runs for the full number of permutations.
- It doesn't make sense as it's fine returning a 0 logLR in the non-permutation part with exactly the same parameters. Why is it returning not quite 0 for the permutation bit? Need to check that logLRs are being treated the same in the first bit as the second.
- noHMT outputs maxLR and perms as: max_logLR: 0, max_logLR_perm: 0
- HMT outputs maxLR and perms as: max_logLR: 5.74651e-13, max_logLR_perm: 4.59244e-13, and the reason why it continues to run all these iterations is because 5.74 > 4.59. Somehow in the non-perm, it gets 5.74, but in all the subsequent perm runs, it only gets 4.59.
- Not sure still.

../../WaveQTL -gmode 1 -g ../../data/geno_data/chr1.2722.geno -p chr.1.2722_WCs.txt -u chr.1.2722_use.txt -o chr.1.2722.nohmt.sml.eps -f 1024 -numPerm 500 -fph 3 > chr.1.2722.nohmt.out
../../WaveQTL -gmode 1 -g ../../data/geno_data/chr1.2722.geno -p chr.1.2722_WCs.txt -u chr.1.2722_use.txt -o chr.1.2722.hmt.sml.eps -f 1024 -numPerm 500 -fph 3 -hmt 1 > chr.1.2722.hmt.out

## 19th August 2021
Goal: understand what the code is doing and outputting regarding parameters text files in both models.
Understand WaveQTL_HMT again, how model works, how it runs in code.
VScode and get gdb for c++!!!

## 31st August 2021
Verified that the R version of WaveQTL_HMT still works. It works on (2,696) on the first SNP. Also verified it on the highest logLR SNP for this (chr,site) as it's the one plotted and which shows the problem space. I've also done a few checks to verify:

- logLR output has logLR, and logBFs; all the logBFs are the same between hmt and nohmt outputs. Other things which should be the same between hmt and nohmt outputs are the mean1 and var1 (as they all come from the bf_uni function). Verified this is true.
- eps, pi are the hyperparams for Wqtl_HMT, and pp and pp_joint are the posterior probabilities of gamma, gamma joint respectively (required to evaluate eps and pi)

To get a better idea of what's going on, it'd be advisable to see if:
- We can speed up the R implementation of WQtl_HMT
- Implement standard WQtl in R, so we can do nice, easy comparisons between the two

## 1st September 2021

So we have a couple of issues we need to address in the R script:
- In BOTH we seemingly have the 'tying_groups' disconnected to the 'groups' in that they don't reflect each other. Groups is calculated before tying_groups, which doesn't really make sense? Groups needs to reflect tying_groups.
- Need to figure out why the WQtl one isn't iterating properly between groups.

## 2nd September 2021
- It's perhaps the scaling coefficient;
	- There needs to be a separate pi hyperparameter for the scaling coefficient; it doesn't appear we have this (please check)
	- The scaling coefficient needs to be a part of the logLR calculation; it also doesn't appear that we have this (please check)
- The scaling coefficient is as per WaveQTL; there are no epsilons/transitions to/from it, it only has a pi (probability of underlying state)

In our plotting functions, we compensate for this; HMT grabs some information from the output of WQtl;
```
  # Grab relevant quantities from WaveQTL (scaling coefficient)
  waveqtl_phi <- as.numeric(as.matrix(read.table(paste0(waveqtl_data_path,waveqtl_dataset,".fph.phi.txt")))[geno_select,2])

```
In this case, it grabs the phi param to get the gamma. I think we need to figure out, for the logLR, how we:
- Update our pi to incorporate the probability gamma = 1 for the scaling coeff
- In WQtl, incorporate the WC into the logLR stat
- Do the same thing for HMT

## 4th September 2021
- I think we've successfully ported the WQtl to R. Need to verify that we're there because the Pi's are a little off.
- Indeed it seems that we've got it right and that adding the scaling coeff will do the trick.
- Our LogLR for HMT is higher in this case (without the WC)! Now we need to figure out how to do the estimation with scaling coeff.

## 11th September 2021
Looked through `C:\Users\brend\Dropbox\Uni Stuff - Masters\Research Project\Masters_Project_Git\code\sim4_functions.R`, and specifically the `run_sim4_v2` function, where (as part of our thesis work), we made adjustments to the LogL to handle the scaling coefficient in the HMT case. We'll need to replicate this logic in our code, but also add the method to calculate LogL using WQtl to the HMT code, rather than relying on post-processing and adding the input ourselves (so we can run it end to end in C++). The method we used:

- HMT LogL is as per the output
- LogL of scaling coefficient is computed using the following:
	+ Grab the 1st pi (pi of scaling coeff) from WQtl output. Note that this output is NOT in a log scale, so can be used raw.
	+ Grab the BF corresponding to the scaling coefficient (ie. the first number in the logLR file after the likelihood), and do 10 to the power of it. This converts the BF (which is expressed as log with base 10) to a BF _NOT_ a logBF.
	+ Formula is that LogL of scaling coefficient is `log(BF*pi + (1-pi))`. In this case, we've done a log with base 10 in our R function.
	+ We need to verify that the logL is actually saved as log with base 10 in the 'logLR.txt' output. __I don't think it is! That is, the logLR is calculated with natural log and saved in the file as a natural log. There's an issue here. Everything is done in natural log, so why have we chosen to do our SC logL in base 10 when adding? We should be doing it in natural log!__
	+ We then add the two together (logLR is additive as LR is a product).
	+ Everything here has been done correctly, but we shouldn't do log base 10.

- If we want to incorporate this into C++, we need to:
	+ Get WQtl in R printing out the same results as per WQtl in C++ [Done, up to like, the 4th decimal place!]
	+ Verify that, for the scaling coeff, WQtl gives the same LogLR value as that if we were to compute using the method above from the outputs. [Also done, also up to like, the 4th decimal place!]
	+ Test by incorporating the WQtl process for scaling coeff in R (iterate the EM over scaling coeff) [DONE]
	+ Verify that it gives the same output as using the above procedure where we add WQtl's outputs as a 'post-process'. [DONE]
	+ If so, need to integrate into C++ code

## 12th September 2021
Figured out that I think the easiest way to incorporate this into C++ code is:
- We run current HMT as it is
- We run the WQtl at the end. We'll need to run a special process for it to pull out the required BF, do all the required things.
- Calculate:
	+ pi
	+ pp
	+ logL
- And then make sure we print it at the appropriate place in the text document. Ie. when we're printing stuff out to our document, we slot in the pi so that:
	+ pi is now two columns wide
	+ logL is the addition of the two
	+ pp ... we'll have to figure out where it goes (what currently gets printed? It looks like we just print a 0 where the pp for the scaling coefficient should be, so just slot it in there.)
This was the pattern we used in the R code which was successful. Do the same way in C++


## 18th September 2021
Trying to integrate scaling coef EM in c++.

```bash
../../WaveQTL -gmode 1 -g ../../data/geno_data/chr4.845.geno -p chr.4.845_WCs.txt -u chr.4.845_use.txt -o chr.4.845.sc -f 1024 -fph 1 -hmt 1

../../WaveQTL -gmode 1 -g ../../data/geno_data/chr4.845.geno -p chr.4.845_WCs.txt -u chr.4.845_use.txt -o chr.4.845.sc.wperm -f 1024 -numPerm 1 -fph 3 -hmt 1
```

- Appears we still haven't gotten rid of the '0 logLR, not significant under noHMT but really significant under HMT' bug. It still appears to be there (see (1,2722)). It's some weird sort of underflow character rounding thing which causes:
	+ the logLR in the HMT setting to be seen as marginally higher than 0 (small positive value) and permute 10000 times (as all the permutations are 0)
	+ the logLR in the non-HMT setting to be seen as pretty much 0, and cut permutations off at 100
	+ even though both logLRs are pretty much zero (when rounded) as are their permuted LogLRs.

- Done! Just need to get it up on the server now.
- Done! Kicked it off. 20210918_run5...

## 20th September 2021
Raced to look at the results (took 24h). A little better than last time.
- Effect_plots.sh isn't quite working properly:
	+ Had to chmod 755 for it to run (why? Did it have different permissions to shell_gen2.sh?)
	+ Couldn't find the effect_size_plot_spartan.R for some reason, even though the file is there.
	+ Changed the starting working dir by which it runs but it still didn't work. Investigate...

## 23rd September 2021
1. Real data [in progress]
- Two groups:
	+ One is both (HMT and WQtl) strong; this signal should be broader
	+ One is HMT strong vs WQtl relatively less strong; this signal should be narrower
- Try and do something like average # of locs significant between the two
- Would want to show there's a big difference between HMT and WQtl at shorter lengths

2. Simulated data [Done (preliminary)]
- Try effect length 8 simulation (as per thesis) to get a grasp over whether new scaling coefficient (different log scale) makes a difference
- Need to figure out where:
	+ `l8 <- readRDS("data/20191010_l8_v3_multiEff_tieg15_od210.RDS")` in `ch3_auroc_images.R` is created. Once this is done, it's straight forward.
	+ Ok, found it. It's in `code/sim4_1_4_singleLengthAnalysis_200.R`
	+ Re-do it with the updated WQtl (so re-run the WQtl HMT as required) given the new log transformed logLRs
	+ Here were some of the details from that run:
	```
	input_data_path = "~/Cpp/WaveQTL_HMT/data/dsQTL/"
	data_path <- "~/Cpp/WaveQTL_HMT/test/dsQTL/output/" # change this when you port it all over to the Masters Git repo
	dataset <- "tree_tie_noQT" # need to replicate the 'tree_tie_noQT' run
	waveqtl_data_path <- "~/Cpp/WaveQTL_HMT/test/dsQTL/output/WaveQTL/"
	waveqtl_dataset <- "test.no.QT"
	geno_select <- 11 # the one used in the demo
	```
	+ Firstly, tried to reproduce the analysis for 100% strength, length 8. [DONE, reproduced (up to the simulation randomness)]
	+ Now, time to dig into the code and figure out what change I need to make to account for the scaling coefficient properly. Turns out all that's required is we change the R code to make sure it converts the scaling coeff's logL to natural log so that it's compatible with the tree's logLR.
	+ I've re-run the analysis for 100% strength, length 8 with this updated R code. I manually use the code from `ch3_auroc_images.R` to use this data and plot some ROC curves, and compare this (and auROC) with that from my thesis.
	+ Outputs saved in /analysis/20211006_paperprep_sim
	+ The impact seems to be minimal. It's worth just running it for all though, just to see.

3. Grouping [Done]
- Try grouping the top 5-6 levels (excl scaling coefficient) as one group (all others as individual groups). Kick that off, see if it makes a difference.
`../../WaveQTL -gmode 1 -group grouping1.txt -g ../../data/geno_data/chr1.825.geno -p chr.1.825_WCs.txt -u chr.1.825_use.txt -o chr.1.825.sc.grp -f 1024 -fph 1 -hmt 1`
- Results look more or less the same
- The histogram looks funny; a sharp peak near '1' for HMT histogram. Is there a bug in the code (argh)?
- Something isn't right, we'll handle it later.

4. Look into the 'red group'; the group where no effect, but HMT says significant and WQtl says not significant

5. Follow-up with Sean about Spartan access on staff email account; lawb@unimelb.edu.au (username: lawb) [Done]

## 7th October 2021
1. Simulations based on chr12:6264339-6265362, rather than the current loc-site we're using at the moment.
	+ chr12:6264339-6265362 is chr12, loc 171 [FOUND, DONE]
	+ Reproduce the effect size plot using WaveQTL, will show that we're on the right track [DONE. Verified is same.]
	+ Reproduce data space shape, should be similar to the top plot [DONE, is similar]
	+ Find effect size location and try and simulate higher or lower effect size
2. Use FDR instead! In all locations; fdr rather than p-value (even in effect size plots, etc)
	+ Can we find the one or two examples? Use a broad brush first and try and look for effect sizes
	+ amend the script to do it vertically; data, nohmt, hmt and p-val, fdr, logLR, real data shape
	+ fdr < 0.05 all examples, print them all.
	+ We want to get a feel to see if it's worth pursuing this avenue of analysis

## 26th October 2021
- Creating a new plotting function so we can have data plot, hmt and no-hmt effect plot all in one.
- Then we can stack (as required) for easier reference
- File in C:\Users\brend\Dropbox\Uni Stuff - Masters\Research Project\Masters_Project_Git\code\20210106_paper_prep\plotting_funcs.R	
To do:
- Data space plots for all [done, just need to move to server and align with server files. test on the chr12.171]
- Add the WaveQTL and WaveQTL_HMT funcitons into that plotting_funcs thing [done]
- Combine so we have three plots on top of each other in a function too. [done]
- Make SNP an input so we can control that too. Make sure we add in q-value [done]

## 10th November 2021
1. Make SNP num an input so we can control that too. Make sure we add in q-value on the plot. This involves adding qvalues to the p-value data, and then having this function read that p-value data and extract the q-value, and plotting it. [DONE]
2. Try new plotting function on server to do some of the other sites, as per signal_strength_picker.R: [DONE]
- Both
- 1,2391 (similar qvals)
- 22,133 (less similar qvals)
- HMT only
- 4,1294 (less similar qvals)
- 16,1985 (more similar qvals)
How did I do it? Copied the script over to Spartan and ran the .R script through the Rscript command. This requires me to edit the R script and 'queue' up multiple function calls over all the sites and chrs I want to plot. Plots saved in run 5's analysis directory.
3. Start simulation work. What do we need? Which file did we use last time to do simulations?

### 15th November 2021
- Two sets of plots; one uses HMT strongest plot and the other uses non-HMT strongest plot. Otherwise one plot. In fact, do this all the time. [DONE]
- 'No HMT' -> 'WaveQTL' [DONE]
- For all FDR < 0.05; how many cases HMT find [DONE]
- Then start on the simulations, and then we may be all good.

## 16th November 2021
Did all the pairs of plots. Use adobe to merge in batches of 100:
Start to 5.1117
5.1397 to 11.2099
11.2391 to 19.437
19.750 to end
11.1448; an example of two separate plots; did we do it correctly? (pg 103-104). [DONE, seems reasonable]

## 28th November 2021
Re-doing simulations on new data. The idea here is that in my thesis, we did simulations based off the counts from the data in the WaveQTL package, which is chr17:10160989-10162012 as per page 14 of WQtl paper. This corresponds to Chr 17, Site 570. 

Now what we want to do is re-run the simulations, but based off the counts from another chr-site which shows more 'narrower' effects. This is chr12:6264339-6265362, or chr12, loc 171. It seems the data space has more sudden and 'localised' peaks for this chr-site. So now we'll re-create the simulation script which was originally in `code/sim4_1_4_singleLengthAnalysis_200.R`. I've copied this script into code/20210106_paper_prep/20211128_chr12_loc171_sims

- Will adapt this script to run in the 'Masters_Project_Git' R Project. [DONE]
- Try effect length 8 simulation first (as per thesis) to get a grasp over the simulation parameters.
- Current code had all the data in the `~/Cpp/WaveQTL_HMT/data/dsQTL` directory. I'm going to try and replicate one in `~/Cpp/20210817_Wqtl/data/dsQTL`. Here's what I need:
	+ Need to make the `read_in_gen_eff_size` function a little more generic to be able to read in the phenotype data I want. [DONE]
	+ HMT dataset used in original script was `WaveQTL_HMT/test/dsQTL/output/tree_tie_noQT` in WaveQTL_HMT folder. Its command was ` ../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.no.QT.txt -u use.txt -o tree_tie_noQT -f 1024 -hmt 1 `. Ie. default grouping, that chr-site's geno, used the no QT wavelet coefficients, with hmt mode.
	+ WQtl dataset used in original script was `WaveQTL_HMT/test/dsQTL/output/WaveQTL/test.no.QT.___`. Its command was `../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.no.QT.txt -u use.txt -o test.no.QT -f 1024 -fph 1` which is the same as above, just on WQtl, NOT HMT mode.
	+ We need:
		+ Chr12-Loc171 phenotype data [DONE]
		+ Outputs from running Chr12-Loc171 on HMT mode, with no grouping, on the No QT WCs [DONE]
		+ Output from running Chr12-Loc171 on WQtl mode, with no grouping, on the No QT WCs [DONE]
- May need to generalise the `read_in_gen_effect_size` function to take in other phenotype file names. (Path is configurable, file name currently isn't) [DONE]
- Change the file savename
- `run_sim4_v2` function in `sim4_functions.R` is the big driver of the work here. Lots of things to change:
	+ Pass input to change all sims output location from `paste0("~/Cpp/WaveQTL_HMT/test/dsQTL/sims/",outputAlias)` to `paste0("~/Cpp/20210817_Wqtl/test/dsQTL/sims/",outputAlias)` [DONE]
	+ Unfortunately, not genericised at the moment. For expedience, make a copy into my sims folder of just what we need. [DONE]
	+ Update the system commands to run WaveQTL and HMT from the current locations, with the new syntax (if required)
	+ Also, can we make this thing run any faster?

### Previous run details from re-running thesis simulations on updated code with scaling coefficient
+ Here were some of the details from that run:
```
input_data_path = "~/Cpp/WaveQTL_HMT/data/dsQTL/"
data_path <- "~/Cpp/WaveQTL_HMT/test/dsQTL/output/" # change this when you port it all over to the Masters Git repo
dataset <- "tree_tie_noQT" # need to replicate the 'tree_tie_noQT' run
waveqtl_data_path <- "~/Cpp/WaveQTL_HMT/test/dsQTL/output/WaveQTL/"
waveqtl_dataset <- "test.no.QT"
geno_select <- 11 # the one used in the demo
```

### 13th December 2021
- In the original thesis runs, I used 'WaveQTL_HMT/test/dsQTL/g15_1024.txt', which grouped all levels 1-5 together as follows: `1 2 33 65 129 257 513` when doing the simulations. Note that these simulations were run off output files which DIDN'T have this grouping; the no QT output files were run with default groupings.
- Send HJ the scripts and a play by play commentary (well commented code) about what these scripts actually do.
- If it's easy to do, port to Spartan
- How to login, submit jobs through slurm etc