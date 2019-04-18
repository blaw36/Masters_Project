## Summary of the C++ code for WaveQTL

The bit we need to change is the implementation of the EM algorithm. This is a brief summary of what the code does, the variables it creates, what those variables should 'look' like and what the script 'does'. The following lines have been chosen because the only cases we're interested in analysing are:

- single_snp_functional_phenotype (functional phenotype), and more specifically
- mode = 1 (not 2 or 3), as we're not interested in finding association between multiple SNPs to a phenotype at once, just one SNP, and therefore, not interested in doing permutation tests for statistical significance for a group of SNPs vs a phenotype at once.

Relevant extracts, from model.cpp
1) Lines 7925 - 8350
2) Lines 8871 - 8887

Also, __think of number of phenotypes__ as equivalent to __number of wavelet coefficients (ie, a scale-loc pair)__

### 1) Lines 7925 - 8350

The aim of this is to run several loops to populate some variables with inputs from the input window, input files and command line option. Then it will run the EM algorithm by iterating over each phenotype (WC sl), by looping, for all individuals, over each _phenotype groupings_, which should all share the same WC = 0 proportion (pi-s).

Order of happenings seems to be:
1) Figure out which phenotypes are to be used, which aren't to be used (filtered out)
2) Assign outfiles to print and append to
3) For each SNP but across all individuals, we run a loop...
  a) get the genotype data for that SNP (across all individuals) and put into 2nd column of the design matrix, gMat
  b) for WCs which are used, write the logBF, (posterior) mean and variance of the beta distributions (for gamma equals 1)
  c) for WCs which are not used, set all to 0
  d) run the EM algorithm: details below
  e) print out to output
  f) summarise all the final parameters (phi)
  g) calculate data space estimations for mean, variance,


The EM algorithm is (for each SNP/covariate/genotype):
1) For each phenotype group (ie. scale level),
	a) Initialise logL to 0, pi to 0.5
	b) Calculate the logLR quantity, posterior denominator, and sum to get poterior and logL of all scale-locs in that group
	c) Now, loop and iterate, altering pi, logPi and log1Pi, as well as resetting logL to 0 so it can be re-calculated. All these quantities represent the E and M steps, and are summed over all scale-locs in that group to get updated quantities.
2) When all done, calculate the final pi for that group, as well as adding the logLR for all WCs in that group to the total logLR.
3) Whole logLR is then composed of all logLRs across all phenotype groups (and therefore, incorporating all individuals and WCs)

- _niter_: max iterations
- _epsilon_: convergence tolerance
- _delta_: for permutation test (not required for us)
- _col_: = 2 in this case, probably just 2 as one mu, one beta (for one SNP at a time)
- _inv_va_: = 0.0. Goes into 'bf_uni' function. Not sure what it does...
- _inv_vd_: = 0.0. Goes into 'bf_uni' function. Not sure what it does... 
- _logbf_: (real) is the log of bayes factor, for a given scale-loc
- _nPH\_use_ (int): is the NUMBER of phenotypes to use
- _nPH\_no_use_ (int): is the NUMBER of phonotypes not to use (should be total - nPH_use)
- _pi\_list_: vector of the pi probabilities we want to calculuate (should be equal to the number of different groups/same scales we want to calculate params for). _in HMT, we'd want to change this to be the length of number of phenos_
- _group\_end_: vector of end positions for each 'scale-group'. _Scale-groups won't really be a thing in HMT so we'd want to get rid of these_
- _nPH\_GP_: vector of number of phenotypes in each group. _again, in HMT, this will be 1 phenotype, for each group_.
- _logLR\_list_: It's a vector of (log) likelihood ratios between the phenotype and each genotype/SNP, for each of how many SNPs (covariates) we are making inferences on. It's used more for permutation testing. _For our purposes, this will really be a 1x1 vector (a scalar)_

#### Seemingly global variables:
- nPH: Assuming it's number of phenotypes <==> number of WCs for each indiv
- col: = 1 + m\_df, which in this case means col = 2.
- nCohort: Assuming this is nIndiv (n), if nPanel = 0. I'm hoping in all our cases, nPanel is always 0. Then nCohort and nIndiv are the same. Otherwise, need to check this.
- nIndiv: Defined in 'model.h' as nCohort + nPanel. What's nPanel?
- vv_phval: WCs, at scale-locs 1,...,nPH, for each individual, 1,...,nIndiv
- nLoci: Assuming it's number of SNPs (covariates)

#### gsl matrices and vectors defined within:
- logBFs: nPH sized vector of log BFs of the evidence supporting gamma not being 0 for each scale-loc
- phMat: n x nPH sized matrix of WCs for each individual, at each scale-loc
- gMat: n x 2 sized matrix - the design matrix for each linear regression; a vector of 1s, and vector of genotype (SNP) values for that SNP 
- ph: an n sized vector of all individuals' WCs at a given scale-loc, to be used in bf_uni, and ONLY if that WC has been marked as one to be used.
- mean1: nPH sized vector of means of the beta coefficient at each scale-loc
- var1:nPH sized vector of variances of the beta coefficient at each scale-loc

#### During the computation (EM algo)
- _use\_pheno\_t_: vector (up to size nPH) of indices of phenotypes to be used (ie the indices of WCs which are marked as '1' in the text file)
- _use\_pheno\_f_: vector (size nPH - size of the use\_pheno\_t vector) of indices of phenotypes not to be used (indices marked as '0')
- _O_obllikli_: old observed log-likelihood ratio (at the start of each iteration, the obs logL from the end of the last iteration)
- _N_obllikli_: new observed log-likelihood ratio (post updating, at the end of each iteration)
- _logLR_: a scalar of the log likelihood ratio a given SNP vs the WCs. It's calculated by summing all the log likelihoods of all the scale-locations, across all individuals
- _pi_: estimate of pi for that scale-loc; sum of posterior probs, divided by the number of scale-locs in that group (pi\_s)
- _numP_: number of scale-locs in that group
- _st_: starting index of the group
- _logpi_: the log pi value at the start of an iteration (like pi\_t). Initialises as log(0.5)
- _log1pi_: the log (1-pi) value (complementary value of logpi)
- _logden_: the log of the denominator used for the posterior (E-step) calulations.
- _logPiBF_: the log-likelihood ratio for each scale-loc
- _pp_: 'posterior probability', which would be the probability of gamma = 1, given data, hyperparams. Sums over all scales in the grouping, to give the numerator of the (M-step) calculation for next pi.
- _diff_: difference in log-likelihood ratio from start to end of iteration for each group. used as stopping criterion

#### Need to figure out what these functions do (if have time):
- bf_uni


### Other things to note
- Most of the mathematical functions are in 'fpmath.cpp'
- Most of the input handling from command line are in 'control.cpp'
- The command line interface should be where main is, which is 'fp.cpp'
- Most of the structs and classes are defined in 'model.h'
- In 'control.cpp' is the definition of when to use 'single_snp_functional_phenotype', as well as its mode, as well as when to use its alternatives. Specifically, in lines 923-935:

```
if(mph > 0) 
{
	pMD->single_snp_multi_phenotype(mph); 
	return; 
}

//--- wavelets start --//
if(fph > 0) 
{
    pMD->read_extra_information_for_functional_phenotype(); 
    pMD->single_snp_functional_phenotype(fph, numPerm, nullcheck);
    return; 
}
```
 where, from the command-line input (lines 191-193): 
```
hcom["-mph"] = com_multi_ph; 
//--- wavelets start --//
hcom["-fph"] = com_functional_ph;
```
- Note that 'fph' and 'mph' both default to 0 unless otherwise initialised through command-line.
- The idea may be, for example, to create 'single_snp_multi_phenotype_hmt' as another function to be used instead, when '-hmt 1' is invoked as an option (in place of -fph >0 or -mph > 0)

#### Main differences req'd with HMT
* No groups, each scale-loc (1,...,nPH) will have its own calculations and probabilities
* Need to follow backward-forward to derive all the betas, all the alphas, then all the gamma\_sl's and joint gamma\_sl,gamma\_p(sl)'s
* Then calculate pi's (for first coefficient), and epsilon's (for all other relevant parent-child relationship), based on the values
* May be easiest to work it backwards
* Not sure how the 'logPiBF' equivalent will look like. ie. need to figure out a logLR ratio equivalent for HMT version when s-l not indp't? But it probably will be the same. The y | gamma_sl are still conditionally independent. The only difference is how we split it up/'de-marginalise' it given that prob gamma | pi is no longer just pi or (1 - pi) but a big long chain of probabilities.

### 2) Lines 8871 - 8887
This is purely just clearing memory and removing vectors/matrices from the environment, etc.