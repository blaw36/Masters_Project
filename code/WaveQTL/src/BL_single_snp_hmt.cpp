// New code for the single_snp_functional_phenotype, now with HMT. Note:
// No modes (no permutation testing will be performed).
// Hence, only the default computation of the EM algorithm will be performed.
// Adapted from single_snp_functional in 'model.cpp'.

//// ~~~~~ The idea is to integrate this into model.cpp (and model.h) later.

// Things which won't change:
// 1) Bayes' Factor calculation should be the same.
// 2) Which WCs to filter out, which to use
// 3) Choice of priors for variances

// Major changes:
// 1) Backward-forward to calculate alphas and betas, from 'top-to-bottom'
// and then from 'bottom-to-top' levels. How does this work across individuals?
// 2) Anything regarding phenotype groups, and numbers in each group will change.
// Each case will be treated individually, and we'll need to loop through each s,l
// according to backward/forward algorithm.
// 3) Reducing the pi_list vector to now be of size 1, as we only need pi_{1,1}, for 
// state = 1. (intiialise in downward algo)
// 4) Inclusion of an eps_list ("epsilon list") vector to be of size (s,l), to capture
// epsilon values for each s,l. There needs to be three of these to capture each of the 
// 3 state combinations required to fully parameterise. eps_list_11, eps_list_01, eps_list_10.
// To be initialised at 0.25, each.
// 5) pp_joint for the joint posterior between pheno and parent, size (s,l) - 1. Each phenotype
// needs to have a pp_joint, and each only has 1 parent. However, (1,1) doesn't have a parent,
// so only (s,l) - 1 values here.
// 6) Up-down algo needs:
  // beta_sl, for each s,l, for states 0, and 1
  // beta_sl_psl, for each s,l, for states 0, and 1
  // beta_psl, for each s,l, for states 0, and 1
  // beta_psl_no_sl, for each s,l, for states 0, and 1
  // alpha_sl, for each s,l, for states 0, and 1
// 7) group_start must be set to a default (), such that group sizes and start/ends
// represent the standard wavelet transform scales; [0, 1, 2, 4, 8, ..., 512, 1024].
// Keeping this structure will allow us to run the algorithm between scale levels.


//-- Pre-amble --//
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include "model.h"
#include "fpmath.h"
#include "control.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_permute.h"
#include "gsl/gsl_cdf.h"
using namespace std;

#if defined (MPI_ENABLED)
#include "mpi.h"
#endif
//-- End Pre-amble --//


// void ModelnData::single_snp_functional_phenotype_HMT(int mode, int numPerm, int nullcheck)
void ModelnData::single_snp_functional_phenotype_HMT(int nullcheck)
{


  int niter = 1000;
  double tol = 0.0005;
  
  //-- WaveQTL.1.1 start --//
  double delta = 0.01;  
  //-- WaveQTL.1.1 end --//

  m_df = 1;
  int col =  1 + m_df; 
  real inv_va = 0.0; 
  real inv_vd = 0.0;
  real logbf;
  int nPH_use;
  int nPH_no_use;
  
  vector<double> pi_list; // log of a prior prob for each phenotype being associated with genotype

  //--- BL_HMT start ---//
  // log of a prior joint (phenotype-parent) prob for that phenotype being associated with genotype
  // 3 different state combinations to fully parameterise
  vector<double> eps_11_list;
  vector<double> eps_10_list;
  vector<double> eps_01_list;
  //--- BL_HMT end ---//

  vector<int> group_end; // end position for each group
  vector<int> nPH_GP;    // number of phenotypes in each group
  vector<double> logLR_list;  // logLR for each SNP (for permutation test)



  // which phenotypes will be used or not
  for(int i = 0; i < nPH; i++){
    if(use_pheno[i] == 1){
      use_pheno_t.push_back(i);
    }else{
      use_pheno_f.push_back(i);
    }
  }
  
  nPH_use = (int) use_pheno_t.size();
  nPH_no_use = (int) use_pheno_f.size(); 



  //set group end positions
  //set number of phenotypes in each group
  //--- BL_HMT start ---//
  // Must ensure this is kept at default, and not changeable for the HMT version.
  //--- BL_HMT end ---//
  int numG = (int)group_start.size();
  group_end.resize(0);
  nPH_GP.resize(0);
  for(int i = 1; i < numG; i++){
    group_end.push_back(group_start[i]-1);
    nPH_GP.push_back(group_start[i] - group_start[i-1]);
  }
  group_end.push_back(nPH);
  nPH_GP.push_back(nPH - group_end[numG-2]);



  if(vsigma_a.size() == 0){
    fplog << "## BIMBAM: Use default priors on additive effects" << endl; 
    vsigma_a.clear(); 
    vsigma_a.push_back(0.05); 
    vsigma_a.push_back(0.1); 
    vsigma_a.push_back(0.2);   
    vsigma_a.push_back(0.4);
  }else{
    fplog << "## BIMBAM: Use user specified priors on additive effects" << endl; 
  }
  
  if(vsigma_d.size() > 0)
    vsigma_d.clear(); 


  gsl_vector * logBFs = gsl_vector_alloc(nPH);
  gsl_matrix * phMat = gsl_matrix_alloc(nCohort, nPH);   
  gsl_matrix * gMat = gsl_matrix_alloc(nCohort, col); 
  gsl_vector * ph = gsl_vector_alloc(nCohort); 


  int ni = 0; 
  for (int i = 0; i < nIndiv; i++){
    if(pIndiv[i]->GetisPanel()) continue; 
    for (int p = 0; p < nPH; p++){
      gsl_matrix_set(phMat, ni, p, vv_phval.at(p).at(i));
    }
    gsl_matrix_set(gMat, ni, 0, 1);
    ni++;
  }   




  fstream outfile; 
  string sfn("output/");
  sfn.append(fnOutput);
  sfn.append(".fph.logLR.txt");
  outfile.open(sfn.c_str(), ios::out);
  if(!outfile.is_open()) {
    cout << "can't open file ... " << endl;  
    exit(0); 
  }

  
  fstream outfile_pi; 
  string sfn_pi("output/");
  sfn_pi.append(fnOutput);
  sfn_pi.append(".fph.pi.txt");
  outfile_pi.open(sfn_pi.c_str(), ios::out);
  if(!outfile_pi.is_open()) {
    cout << "can't open file ... " << endl;  
    exit(0); 
  }

  //--- BL_HMT start ---//
  fstream outfile_eps_11; 
  string sfn_eps_11("output/");
  sfn_eps_11.append(fnOutput);
  sfn_eps_11.append(".fph.eps.txt");
  outfile_eps_11.open(sfn_eps_11.c_str(), ios::out);
  if(!outfile_eps_11.is_open()) {
    cout << "can't open file ... " << endl;  
    exit(0); 
  }

  fstream outfile_eps_10; 
  string sfn_eps_10("output/");
  sfn_eps_10.append(fnOutput);
  sfn_eps_10.append(".fph.eps.txt");
  outfile_eps_10.open(sfn_eps_10.c_str(), ios::out);
  if(!outfile_eps_10.is_open()) {
    cout << "can't open file ... " << endl;  
    exit(0); 
  }

  fstream outfile_eps_01; 
  string sfn_eps_01("output/");
  sfn_eps_01.append(fnOutput);
  sfn_eps_01.append(".fph.eps.txt");
  outfile_eps_01.open(sfn_eps_01.c_str(), ios::out);
  if(!outfile_eps_01.is_open()) {
    cout << "can't open file ... " << endl;  
    exit(0); 
  }    
  //--- BL_HMT end ---//


  fstream outfile_mean; 
  string sfn_mean("output/");
  sfn_mean.append(fnOutput);
  sfn_mean.append(".fph.mean.txt");
  outfile_mean.open(sfn_mean.c_str(), ios::out);
  if(!outfile_mean.is_open()) {
    cout << "can't open file ... " << endl;  
    exit(0); 
  }


  fstream outfile_var; 
  string sfn_var("output/");
  sfn_var.append(fnOutput);
  sfn_var.append(".fph.var.txt");
  outfile_var.open(sfn_var.c_str(), ios::out);
  if(!outfile_var.is_open()) {
    cout << "can't open file ... " << endl;  
    exit(0); 
  }


  gsl_vector * mean1 = gsl_vector_alloc(nPH);
  gsl_vector * var1 = gsl_vector_alloc(nPH);

  //--- wavelets_v1.3 end ---//


  //--- wavelets_v2_2 start ---//
  
  fstream outfile_mean1; 
  string sfn_mean1("output/");
  sfn_mean1.append(fnOutput);
  sfn_mean1.append(".fph.mean1.txt");
  outfile_mean1.open(sfn_mean1.c_str(), ios::out);
  if(!outfile_mean1.is_open()) {
    cout << "can't open file ... " << endl;  
    exit(0); 
  }

  fstream outfile_var1; 
  string sfn_var1("output/");
  sfn_var1.append(fnOutput);
  sfn_var1.append(".fph.var1.txt");
  outfile_var1.open(sfn_var1.c_str(), ios::out);
  if(!outfile_var1.is_open()) {
    cout << "can't open file ... " << endl;  
    exit(0); 
  }

  fstream outfile_phi; 
  string sfn_phi("output/");
  sfn_phi.append(fnOutput);
  sfn_phi.append(".fph.phi.txt");
  outfile_phi.open(sfn_phi.c_str(), ios::out);
  if(!outfile_phi.is_open()) {
    cout << "can't open file ... " << endl;  
    exit(0); 
  }

  //--- wavelets_v2_2 end ---//
  
  logLR_list.resize(0);
  // Calculate for each genotype!!
  for(int g = 0; g < nLoci; g++){

    // assign each genotype 
    ni = 0;
    for (int i = 0; i < nIndiv; i++){
      if(pIndiv[i]->GetisPanel()) continue; 
      gsl_matrix_set(gMat, ni, 1, pIndiv[i]->get_snpmgt(g));
      ni++; 
    }	

    //--- wavelets_v1.3 start ---//
    double * res = new double[3]; 
    //--- wavelets_v1.3 end ---//

    for(int i =0; i < nPH_use; i++){
      int p = use_pheno_t[i];
      gsl_matrix_get_col(ph, phMat, p);
      //--- wavelets_v1.3 start ---//
      bf_uni(inv_va, inv_vd, nCohort, 1, gMat, ph, res);
      gsl_vector_set(logBFs, p, res[0]);
      gsl_vector_set(mean1, p, res[1]);
      gsl_vector_set(var1, p, res[2]);
      //--- wavelets_v1.3 end ---//
    }

    //--- wavelets_v1.3 start ---//
    delete[] res; 
    //--- wavelets_v1.3 end ---//


    for(int i =0; i < nPH_no_use; i++){
      int p = use_pheno_f[i];  
      gsl_vector_set(logBFs, p, 0.0);
      //--- wavelets_v1.3 start ---//
      gsl_vector_set(mean1, p, 0.0);
      gsl_vector_set(var1, p, 0.0);
      //--- wavelets_v1.3 end ---//
    }



  //--- BL_HMT start ---//

  /**********************/
  //     EM algorithm    //
  /***********************/

    double O_obllikli;
    double N_obllikli;
    double logLR;
    pi_list.resize(0);
    eps_11_list.resize(0);
    eps_10_list.resize(0);
    eps_01_list.resize(0);

  // With HMT, need to do the EM algorithm for all groups (scales) in one step.
  // Can no longer do separately as they are all linked together.
  // ~~~ TO CONFIRM: nCohort same as nIndiv?

  // ---- Initialisation of parameters ---- //

  // Alpha and beta quantities
    gsl_vector * beta_sl_0 = gsl_vector_alloc(nPH);
    gsl_vector * beta_sl_1 = gsl_vector_alloc(nPH);

    gsl_vector * beta_sl_psl_1 = gsl_vector_alloc(nPH);
    gsl_vector * beta_sl_psl_0 = gsl_vector_alloc(nPH);

    gsl_vector * beta_psl_no_sl_1 = gsl_vector_alloc(nPH);
    gsl_vector * beta_psl_no_sl_0 = gsl_vector_alloc(nPH);

    gsl_vector * alpha_sl_1 = gsl_vector_alloc(nPH);
    gsl_vector * alpha_sl_0 = gsl_vector_alloc(nPH);

    gsl_vector * pp = gsl_vector_alloc(nPH);
    gsl_vector * pp_joint_11 = gsl_vector_alloc(nPH);
    gsl_vector * pp_joint_10 = gsl_vector_alloc(nPH);
    gsl_vector * pp_joint_01 = gsl_vector_alloc(nPH);

  // Vector of epsilon quantities
    gsl_vector * eps_11_vect = gsl_vector_alloc(nPH);
    gsl_vector * eps_10_vect = gsl_vector_alloc(nPH);
    gsl_vector * eps_01_vect = gsl_vector_alloc(nPH);

  // ---- Initial E-step ---- //
  // 'Zero-th' up-down algorithm using all individuals

  // This pi now only represents the 1 node at the coarsest scale (1,1)
  // Probably easier to have pi not logged this time (?)
  // double logpi = log(0.5);
  // double log1pi = logpi;
    double pi_1 = 0.5;
    double pi_0 = 1 - pi_1;


  // Probably (?) easier to also have epsilons not logged
  // double logeps_11 = log(0.25);
  // double logeps_01 = log(0.25);
  // double logeps_10 = log(0.25);
  // double logeps_00 = log(0.25);
    double eps_11 = 0.25;
    double eps_01 = 0.25;
    double eps_10 = 0.25;
    double eps_00 = 1 - eps_11 - eps_01 - eps_10;

  // Define some variables
    double logBF, logParentBF;

    double b_sl_1, b_sl_0;
  // The beta of the adjacent child (which shares the same parent)
    double b_sl_adj_1, b_sl_adj_0;

    double b_sl_psl_1, b_sl_psl_0;
    double b_psl_1, b_psl_0;
    double b_psl_no_sl_1, b_psl_no_sl_0;

    double a_sl_1, a_sl_0;
    double a_psl_1, a_psl_0;

    double denom;
    double pp_sl;
    double pp_psl;
    double pp_joint_sl_11, pp_joint_sl_10, pp_joint_sl_01;

  // Up step - don't do this step at coarsest scale
    for(int gi = (numG - 1); gi > 0; gi--){

      int numP = nPH_GP[gi];
      int st = group_start[gi] - 1;
      int parent_st = group_start[gi - 1] - 1;

    // Step 0 at finest scale only.
      if(gi == (numG - 1)){
        for(int i = 0; i < numP; i++){

          logBF = gsl_vector_get(logBFs,st+i);

          gsl_vector_set(beta_sl_1, st + i, exp(logBF));
          gsl_vector_set(beta_sl_0, st + i, 1);
        }
      }

      for(int i = 0; i < numP; i++){

        int parent_indx = wc/2; // integer division
        logParentBF = gsl_vector_get(logBFs,parent_st+parent_indx);

        b_sl_1 = gsl_vector_get(beta_sl_1, st + i);
        b_sl_0 = gsl_vector_get(beta_sl_0, st + i);

      // If even, grab one to the right, else one to the left
        if(wc % 2 == 0){
          b_sl_adj_1 = gsl_vector_get(beta_sl_1, st + i + 1);
          b_sl_adj_0 = gsl_vector_get(beta_sl_0, st + i + 1);
        }else{
          b_sl_adj_1 = gsl_vector_get(beta_sl_1, st + i - 1);
          b_sl_adj_0 = gsl_vector_get(beta_sl_0, st + i - 1);
        }

        b_sl_psl_1 = b_sl_1*eps_11 + b_sl_0*eps_01;
        b_sl_psl_0 = b_sl_0*eps_00 + b_sl_1*eps_10;

        b_psl_1 = exp(logParentBF) * b_sl_adj_1 * b_sl_1;
        b_psl_0 = 1 * b_sl_adj_0 * b_sl_0;

        b_psl_no_sl_1 = b_psl_1/b_sl_psl_1;
        b_psl_no_sl_0 = b_psl_0/b_sl_psl_0;

        gsl_vector_set(beta_sl_psl_1, st + i, b_sl_psl_1);
        gsl_vector_set(beta_sl_psl_0, st + i, b_sl_psl_0);

      // Both the children will try and update the parents' beta, so
      // only do this once.
        if(i % 2 == 0){
          gsl_vector_set(beta_sl_1, parent_st + parent_indx, b_psl_1);
          gsl_vector_set(beta_sl_0, parent_st + parent_indx, b_psl_0);  
        }

        gsl_vector_set(beta_psl_no_sl_1, st + i, b_psl_no_sl_1);
        gsl_vector_set(beta_psl_no_sl_0, st + i, b_psl_no_sl_0);

      }

    }

  // Down step
  // Start at coarsest scale, then go until finest - 1.
    for(int gi = 0; gi < (numG - 1); gi++){

      int numP = nPH_GP[gi];
      int st = group_start[gi] - 1;
      int parent_st = group_start[gi - 1] - 1;

    // Step 0 at finest scale only.
      if(gi == 0){

        gsl_vector_set(alpha_sl_1, st, pi_1);
        gsl_vector_set(alpha_sl_0, st, pi_0);

      }else{

        for(int i = 0; i < numP; i++){

          b_psl_no_sl_1 = gsl_vector_get(beta_sl_psl_1, st + i);
          b_psl_no_sl_0 = gsl_vector_get(beta_sl_psl_0, st + i);

          int parent_indx = i/2;
          a_psl_1 = gsl_vector_get(alpha_sl_1, parent_st + parent_indx);
          a_psl_0 = gsl_vector_get(alpha_sl_0, parent_st + parent_indx);

          a_sl_1 = eps_11*a_psl_1*b_psl_no_sl_1 + eps_10*a_psl_0*b_psl_no_sl_0;
          a_sl_0 = eps_01*a_psl_1*b_psl_no_sl_1 + eps_00*a_psl_0*b_psl_no_sl_0;
          gsl_vector_set(alpha_sl_1, st + i, a_sl_1);
          gsl_vector_set(alpha_sl_0, st + i, a_sl_0);

        }
      }
      
    } 

  // ---- Calculate the E-step quantities (posterior probabilities) ---- //
    for(int wc = 0; wc < nPH; wc++){

      // Posterior marginal of gamma_sl
      b_sl_1 = gsl_vector_get(beta_sl_1, wc);
      b_sl_0 = gsl_vector_get(beta_sl_0, wc);
      a_sl_1 = gsl_vector_get(alpha_sl_1, wc);
      a_sl_0 = gsl_vector_get(alpha_sl_0, wc);

      denom = b_sl_1*a_sl_1 + b_sl_0*a_sl_0;
      pp_sl = (b_sl_1*a_sl_1)/denom;
      gsl_vector_set(pp, wc, pp_sl);

    // Joint marginal NOT DEFINED for (1,1) - set to 0 for now.
      if(wc > 0){
      // Posterior joint of gamma_sl, gamma_psl
      // This is the easiest way I can think of (for now) to get the parents' index.
      // Use the 'divide' by 2 trick, on a tree with starting index of 1 (not 0)
        int parent_wc = ((wc + 1) / 2) - 1;
        a_psl_1 = gsl_vector_get(alpha_sl_1, parent_wc);
        a_psl_0 = gsl_vector_get(alpha_sl_0, parent_wc);
        b_psl_no_sl_1 = gsl_vector_get(beta_psl_no_sl_1, wc);
        b_psl_no_sl_0 = gsl_vector_get(beta_psl_no_sl_0, wc);       

        pp_joint_sl_11 = (b_sl_1*eps_11*a_psl_1*b_psl_no_sl_1)/denom;
        pp_joint_sl_10 = (b_sl_1*eps_10*a_psl_0*b_psl_no_sl_0)/denom;
        pp_joint_sl_01 = (b_sl_0*eps_01*a_psl_1*b_psl_no_sl_1)/denom;

        gsl_vector_set(pp_joint_11, wc, pp_joint_sl_11);
        gsl_vector_set(pp_joint_10, wc, pp_joint_sl_10);
        gsl_vector_set(pp_joint_01, wc, pp_joint_sl_01);  
      }else{
        gsl_vector_set(pp_joint_11, wc, 0);
        gsl_vector_set(pp_joint_10, wc, 0);
        gsl_vector_set(pp_joint_01, wc, 0);
      }

    }

  // ---- Calculate initial observed log-likelihood ---- //
  // Observed log-likelihood calculated by marginalising beta_11 over the possible
  // states of gamma_11.

    double b_11_1 = gsl_vector_get(beta_sl_1, 0);
    double b_11_0 = gsl_vector_get(beta_sl_0, 0);

    N_obllikli = log(b_11_1*pi_1 + b_11_0*pi_0);
    O_obllikli = N_obllikli;


  // ---- Iterate through first M-step, and subsequent E-steps ---- //
    for(int iter = 0; iter < niter; iter++){

    // ---- M-step ---- //
    // ~~~ TOFIX: Aggregating several parameters to get 'grouped' parameter estimates
    // for pi and epsilons.
    // ie A may be average of A's from several nodes near the top,
    // B may be average of all B/A's from all locations in that scale.

    // Pi param, for (1,1)
      pi_1 = gsl_vector_get(pp, 0);
      pi_0 = 1 - pi_1;

    // Eps param, for each s,l apart from (1,1)
      for(int i = 0; i < numP; i++){
      // Joint marginal NOT DEFINED for (1,1) - set to 0 for now.
        if(i > 0){

          pp_joint_sl_11 = gsl_vector_get(pp_joint_11, i);
          pp_joint_sl_10 = gsl_vector_get(pp_joint_10, i);
          pp_joint_sl_01 = gsl_vector_get(pp_joint_01, i);
          pp_joint_sl_00 = 1 - pp_joint_sl_11 - pp_joint_sl_10 - pp_joint_sl_01;

          int parent_wc = ((i + 1) / 2) - 1;
          pp_psl = gsl_vector_get(pp, parent_wc);

          eps_11 = pp_joint_sl_11/pp_psl;
          eps_01 = pp_joint_sl_01/pp_psl;
          eps_10 = pp_joint_sl_10/(1-pp_psl);

          gsl_vector_set(eps_11_vect, i, eps_11);
          gsl_vector_set(eps_01_vect, i, eps_01);
          gsl_vector_set(eps_10_vect, i, eps_10);

        }else{

          gsl_vector_set(eps_11_vect, i, 0);
          gsl_vector_set(eps_01_vect, i, 0);
          gsl_vector_set(eps_10_vect, i, 0);

        }
      }

    // Reset obs log-likelihood for new iteration
      N_obllikli = 0;

    // Up step - don't do this step at coarsest scale
      for(int gi = (numG - 1); gi > 0; gi--){

        int numP = nPH_GP[gi];
        int st = group_start[gi] - 1;
        int parent_st = group_start[gi - 1] - 1;

      // Step 0 at finest scale only.
      // Does this remain the same with each iteration?
        if(gi == (numG - 1)){
          for(int i = 0; i < numP; i++){

            logBF = gsl_vector_get(logBFs,st + i);

            gsl_vector_set(beta_sl_1, st + i, exp(logBF));
            gsl_vector_set(beta_sl_0, st + i, 1);
          }
        }

        for(int i = 0; i < numP; i++){

          int parent_indx = i/2; // integer division
          logParentBF = gsl_vector_get(logBFs,parent_st + parent_indx);

          eps_11 = gsl_vector_get(eps_11_vect, st + i);
          eps_01 = gsl_vector_get(eps_01_vect, st + i);
          eps_10 = gsl_vector_get(eps_10_vect, st + i);
          eps_00 = 1 - eps_11 - eps_01 - eps_10;

          b_sl_1 = gsl_vector_get(beta_sl_1, st + i);
          b_sl_0 = gsl_vector_get(beta_sl_0, st + i);

        // If even, grab one to the right, else one to the left
          if(i % 2 == 0){
            b_sl_adj_1 = gsl_vector_get(beta_sl_1, st + i + 1);
            b_sl_adj_0 = gsl_vector_get(beta_sl_0, st + i + 1);
          }else{
            b_sl_adj_1 = gsl_vector_get(beta_sl_1, st + i - 1);
            b_sl_adj_0 = gsl_vector_get(beta_sl_0, st + i - 1);
          }

          b_sl_psl_1 = b_sl_1*eps_11 + b_sl_0*eps_01;
          b_sl_psl_0 = b_sl_0*eps_00 + b_sl_1*eps_10;

          b_psl_1 = exp(logParentBF) * b_sl_adj_1 * b_sl_1;
          b_psl_0 = 1 * b_sl_adj_0 * b_sl_0;

          b_psl_no_sl_1 = b_psl_1/b_sl_psl_1;
          b_psl_no_sl_0 = b_psl_0/b_sl_psl_0;

          gsl_vector_set(beta_sl_psl_1, st + i, b_sl_psl_1);
          gsl_vector_set(beta_sl_psl_0, st + i, b_sl_psl_0);

        // Both the children will try and update the parents' beta, so
        // only do this once.
          if(wc % 2 == 0){
            gsl_vector_set(beta_sl_1, parent_st + parent_indx, b_psl_1);
            gsl_vector_set(beta_sl_0, parent_st + parent_indx, b_psl_0);  
          }

          gsl_vector_set(beta_psl_no_sl_1, st + i, b_psl_no_sl_1);
          gsl_vector_set(beta_psl_no_sl_0, st + i, b_psl_no_sl_0);

        }

      }

    // Down step
    // Start at coarsest scale, then go until finest - 1.
      for(int gi = 0; gi < (numG - 1); gi++){

        int numP = nPH_GP[gi];
        int st = group_start[gi] - 1;
        int parent_st = group_start[gi - 1] - 1;

      // Step 0 at finest scale only.
        if(gi == 0){

          gsl_vector_set(alpha_sl_1, st, pi_1);
          gsl_vector_set(alpha_sl_0, st, pi_0);

        }else{

          for(int i = 0; i < numP; i++){

            b_psl_no_sl_1 = gsl_vector_get(beta_sl_psl_1, st + i);
            b_psl_no_sl_0 = gsl_vector_get(beta_sl_psl_0, st + i);

            eps_11 = gsl_vector_get(eps_11_vect, st + i);
            eps_01 = gsl_vector_get(eps_01_vect, st + i);
            eps_10 = gsl_vector_get(eps_10_vect, st + i);
            eps_00 = 1 - eps_11 - eps_01 - eps_10;

            int parent_indx = i/2;
            a_psl_1 = gsl_vector_get(alpha_sl_1, parent_st + parent_indx);
            a_psl_0 = gsl_vector_get(alpha_sl_0, parent_st + parent_indx);

            a_sl_1 = eps_11*a_psl_1*b_psl_no_sl_1 + eps_10*a_psl_0*b_psl_no_sl_0;
            a_sl_0 = eps_01*a_psl_1*b_psl_no_sl_1 + eps_00*a_psl_0*b_psl_no_sl_0;
            gsl_vector_set(alpha_sl_1, st + i, a_sl_1);
            gsl_vector_set(alpha_sl_0, st + i, a_sl_0);

          }
        }

      } 

    // ---- Calculate the E-step quantities (posterior probabilities) ---- //
      for(int wc = 0; wc < nPH; wc++){

        // Posterior marginal of gamma_sl
        b_sl_1 = gsl_vector_get(beta_sl_1, wc);
        b_sl_0 = gsl_vector_get(beta_sl_0, wc);
        a_sl_1 = gsl_vector_get(alpha_sl_1, wc);
        a_sl_0 = gsl_vector_get(alpha_sl_0, wc);

        denom = b_sl_1*a_sl_1 + b_sl_0*a_sl_0;
        pp_sl = (b_sl_1*a_sl_1)/denom;
        gsl_vector_set(pp, wc, pp_sl);

      // Joint marginal NOT DEFINED for (1,1) - set to 0 for now.
        if(wc > 0){
        // Posterior joint of gamma_sl, gamma_psl
        // This is the easiest way I can think of (for now) to get the parents' index.
        // Use the 'divide' by 2 trick, on a tree with starting index of 1 (not 0)
          int parent_wc = ((wc + 1) / 2) - 1;
          a_psl_1 = gsl_vector_get(alpha_sl_1, parent_wc);
          a_psl_0 = gsl_vector_get(alpha_sl_0, parent_wc);
          b_psl_no_sl_1 = gsl_vector_get(beta_psl_no_sl_1, wc);
          b_psl_no_sl_0 = gsl_vector_get(beta_psl_no_sl_0, wc);       

          eps_11 = gsl_vector_get(eps_11_vect, st + i);
          eps_01 = gsl_vector_get(eps_01_vect, st + i);
          eps_10 = gsl_vector_get(eps_10_vect, st + i);
          eps_00 = 1 - eps_11 - eps_01 - eps_10;

          pp_joint_sl_11 = (b_sl_1*eps_11*a_psl_1*b_psl_no_sl_1)/denom;
          pp_joint_sl_10 = (b_sl_1*eps_10*a_psl_0*b_psl_no_sl_0)/denom;
          pp_joint_sl_01 = (b_sl_0*eps_01*a_psl_1*b_psl_no_sl_1)/denom;

          gsl_vector_set(pp_joint_11, wc, pp_joint_sl_11);
          gsl_vector_set(pp_joint_10, wc, pp_joint_sl_10);
          gsl_vector_set(pp_joint_01, wc, pp_joint_sl_01);  
        }else{
          gsl_vector_set(pp_joint_11, wc, 0);
          gsl_vector_set(pp_joint_10, wc, 0);
          gsl_vector_set(pp_joint_01, wc, 0);
        }

      }

      // ---- Calculate initial observed log-likelihood ---- //
      // Observed log-likelihood calculated by marginalising beta_11 over the possible
      // states of gamma_11.

      double b_11_1 = gsl_vector_get(beta_sl_1, 0);
      double b_11_0 = gsl_vector_get(beta_sl_0, 0);

      N_obllikli = log(b_11_1*pi_1 + b_11_0*pi_0);
      
      diff = N_obllikli - O_obllikli;

      if(diff < tol){
        break;
      }else{
        O_obllikli = N_obllikli;    
      }
    }

  // ---- Print things to outfiles ---- //

  // Only one of these values for the whole algorithm
    logLR = N_obllikli;
    pi_list.push_back(pi_1);
    logLR_list.push_back(logLR);

    for(int p = 0; p < nPH; p++){
      if(p == 0){
        eps_11_list.push_back(0);
        eps_10_list.push_back(0);
        eps_01_list.push_back(0);
      }else{
        eps_11 = gsl_vector_get(eps_11_vect, p);
        eps_01 = gsl_vector_get(eps_01_vect, p);
        eps_10 = gsl_vector_get(eps_10_vect, p);

        eps_11_list.push_back(eps_11);
        eps_10_list.push_back(eps_01);
        eps_01_list.push_back(eps_10);
      }
    }

  // logLR and logBF    
    outfile << vsRsnum.at(g) << " "; 
    char buf[100]; 

    if(logLR < 1e-5){
      sprintf(buf, "%.5f ", logLR); 
    }else{
      sprintf(buf, "%+.5f ", logLR); 
    }
    outfile << buf;

    for(int p = 0; p < nPH; p++){

      logbf = gsl_vector_get(logBFs, p)/log(10);

      if(logbf < 1e-5){
        sprintf(buf, "%.5f ", logbf);
      }else{
        sprintf(buf, "%+.5f ", logbf);  
      }
      outfile << buf;
    }
    outfile << endl;

  // pi
    outfile_pi << vsRsnum.at(g) << " ";

    if(pi_1 < 1e-5){
      sprintf(buf, "%.5f ", pi_1);
    }else{
      sprintf(buf, "%+.5f ", pi_1); 
    }
    outfile_pi << buf;
    outfile_pi << endl; 

  // epsilon
    outfile_eps_11 << vsRsnum.at(g) << " ";
    outfile_eps_10 << vsRsnum.at(g) << " ";
    outfile_eps_01 << vsRsnum.at(g) << " ";
    for(int p = 0; p < nPH; p++){
      eps_11 = eps_11_list[p];
      if(eps_11 < 1e-5){
        sprintf(buf, "%.5f ", eps_11);  
      }else{
        sprintf(buf, "%+.5f ", eps_11);   
      }
      outfile_eps_11 << buf;

      eps_10 = eps_10_list[p];
      if(eps_10 < 1e-5){
        sprintf(buf, "%.5f ", eps_10);  
      }else{
        sprintf(buf, "%+.5f ", eps_10);   
      }
      outfile_eps_10 << buf;

      eps_01 = eps_01_list[p];
      if(eps_01 < 1e-5){
        sprintf(buf, "%.5f ", eps_01);  
      }else{
        sprintf(buf, "%+.5f ", eps_01);   
      }
      outfile_eps_01 << buf;
    }
    outfile_eps_11 << endl;
    outfile_eps_10 << endl;
    outfile_eps_01 << endl;

  // ---- Compile estimates of phi, means and vars ---- //

  // phi, mean, and var
    double phi, mean_out, var_out, mean1_out, var1_out;
    outfile_mean << vsRsnum.at(g) << " ";  
    outfile_var << vsRsnum.at(g) << " ";
    //-- wavelets_v2_2 start --//
    outfile_mean1 << vsRsnum.at(g) << " ";  
    outfile_var1 << vsRsnum.at(g) << " ";
    outfile_phi << vsRsnum.at(g) << " ";
    //-- wavelets_v2_2 end --//

  // ~~~ TOFIX; Means, vars, likely to change with simulation.
    for(int p = 0, p < nPH, p++){
      if(use_pheno[p] == 1){
        phi = pp[p];
        mean1_out = gsl_vector_get(mean1,p);
        var1_out = gsl_vector_get(var1,p);
        mean_out = phi*mean1_out;
        var_out = phi*(var1_out + mean1_out*mean1_out*(1-phi));
        outfile_mean << mean_out << " ";
        outfile_var << var_out << " ";
      //-- wavelets_v2_2 start --//
        outfile_mean1 << mean1_out << " ";
        outfile_var1 << var1_out << " ";
        outfile_phi << phi << " ";
      //-- wavelets_v2_2 end --//
      }else{
        sprintf(buf, "%.5f ", 0.0);
        outfile_mean << buf;
        outfile_var << buf;
      //-- wavelets_v2_2 start --//
        outfile_mean1 << buf;
        outfile_var1 << buf;
        outfile_phi << buf;
      //-- wavelets_v2_2 end --//
      }
    }
  //--- BL_HMT end ---/

    outfile_mean << endl;
    outfile_var << endl;
    //-- wavelets_v2_2 start --//
    outfile_mean1 << endl;
    outfile_var1 << endl;
    outfile_phi << endl; 
   //-- wavelets_v2_2 end --//


  }




  outfile.close(); 
  cout << sfn << " has been created." << endl;


  outfile_pi.close(); 
  cout << sfn_pi << " has been created." << endl;


  //--- BL_HMT start ---//
  outfile_eps_11.close(); 
  cout << sfn_eps_11 << " has been created." << endl;

  outfile_eps_10.close(); 
  cout << sfn_eps_10 << " has been created." << endl;

  outfile_eps_01.close(); 
  cout << sfn_eps_01 << " has been created." << endl;

  outfile_eps_00.close(); 
  cout << sfn_eps_00 << " has been created." << endl;
  //--- BL_HMT end ---//


  //--- wavelets_v1.3 start ---//


  outfile_mean.close(); 
  cout << sfn_mean << " has been created." << endl;


  outfile_var.close(); 
  cout << sfn_var << " has been created." << endl;


  //-- wavelets_v2_2 start --//
  outfile_mean1.close(); 
  cout << sfn_mean1 << " has been created." << endl;

  outfile_var1.close(); 
  cout << sfn_var1 << " has been created." << endl;

  outfile_phi.close(); 
  cout << sfn_phi << " has been created." << endl;

  pi_list.resize(0);
  //--- BL_HMT start ---//
  eps_11_list.resize(0);
  eps_10_list.resize(0);
  eps_01_list.resize(0);
  //--- BL_HMT end ---//
  nPH_GP.resize(0);
  group_end.resize(0);


  group_start.resize(0);
  use_pheno.resize(0);
  use_pheno_t.resize(0);
  use_pheno_f.resize(0);

  gsl_vector_free(ph); 
  gsl_vector_free(logBFs); 
  gsl_matrix_free(gMat); 
  gsl_matrix_free(phMat);


  logLR_list.resize(0);

  //--- BL_HMT start ---//
  gsl_vector_free(beta_sl_0);
  gsl_vector_free(beta_sl_1);
  gsl_vector_free(beta_sl_psl_1);
  gsl_vector_free(beta_sl_psl_0);
  gsl_vector_free(beta_psl_no_sl_1);
  gsl_vector_free(beta_psl_no_sl_0);
  gsl_vector_free(alpha_sl_1);
  gsl_vector_free(alpha_sl_0);
  gsl_vector_free(pp);
  gsl_vector_free(pp_joint_11);
  gsl_vector_free(pp_joint_10);
  gsl_vector_free(pp_joint_01);
  gsl_vector_free(eps_11_vect);
  gsl_vector_free(eps_10_vect);
  gsl_vector_free(eps_01_vect);
  //--- BL_HMT end ---//

}