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
permuted likelihood ratios that leads to the p-values?
With chr1, site2722, why did the permutations calculate so weirdly?
Print out logLR for nohmt permutations

Debug case: 
Chr 2, Site 696 (HMT poor and no-HMT good)
Chr 20, Site 1364 (HMT poor and no-HMT good but should be HMT stronger?)
effect size in original algo space (with quantile transform); Make effect size WITH quantile transform (ie. original data which we ran algo on) for everything
logLR and logLR at different scales
Pi param estimates at different scales

## 17th August 2021
1) Check which pheno data we ran our algo on.
	In our WaveQTL runscripts, we used the following genotype data:
	`/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno"`
	In data size plot, we used:
	`as.matrix(read.table(paste0("/data/gpfs/projects/punim0614/shared/shared_data/internal/multi.scale/WaveQTL/DNase/geno_01_step1/geno_maf/chr",chrIX,".",site,".geno")))`
	They're the same, no further work required. Phew!

Outstanding work: