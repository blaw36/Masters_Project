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
  vector<double> eps_list_11;
  vector<double> eps_list_10;
  vector<double> eps_list_01;
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
  fstream outfile_eps; 
  string sfn_eps("output/");
  sfn_eps.append(fnOutput);
  sfn_eps.append(".fph.eps.txt");
  outfile_eps.open(sfn_eps.c_str(), ios::out);
  if(!outfile_eps.is_open()) {
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



    /**********************/
    //     EM algorithm    //
    /***********************/

    double O_obllikli;
    double N_obllikli;
    pi_list.resize(0);

    double logLR = 0;
    double pi;
    for(int gi = 0; gi < numG; gi++){

      // EM algorithm for each group separately 


      N_obllikli = 0;
      int numP = nPH_GP[gi];
      int st = group_start[gi] - 1;
      double logpi = log(0.5);
      double log1pi = logpi;
      
      double logden;
      double logPiBF;
      double pp;
      double diff;

      pp = 0;
      for(int i = 0; i < numP; i++){
       logPiBF = gsl_vector_get(logBFs,st+i) + logpi;
       logden = sumlog(logPiBF, log1pi);
       pp += exp(logPiBF - logden);
       N_obllikli += logden;
     } 
     O_obllikli = N_obllikli;


     for(int iter = 0; iter < niter; iter++){

       pi = pp/(double)numP;
       logpi  = log(pi);
       log1pi = log(1-pi);
       N_obllikli = 0;

       pp = 0;
       for(int i = 0; i < numP; i++){
         logPiBF = gsl_vector_get(logBFs,st+i) + logpi;
         logden = sumlog(logPiBF, log1pi);
         pp += exp(logPiBF - logden);
         N_obllikli += logden;
       }

       diff = N_obllikli - O_obllikli;

       if(diff < tol){
         break;
       }else{
         O_obllikli = N_obllikli;		
       }
     }

      //-- WaveQTL.1.1 start --//
     if(nullcheck == 1){
       if(N_obllikli < 0){
         N_obllikli = 0;
         pp = 0;
       }
     }
      //-- WaveQTL.1.1 end --//

     pi_list.push_back(pp/(double)numP);
     logLR += N_obllikli;

   }

   logLR_list.push_back(logLR);



    // logLR and logBF    
   outfile << vsRsnum.at(g) << " "; 
   char buf[100]; 

   if(logLR < 1e-5) 
    sprintf(buf, "%.5f ", logLR); 
  else 
    sprintf(buf, "%+.5f ", logLR); 
  outfile << buf;
  
  for(int p =0; p < nPH; p++){

    logbf = gsl_vector_get(logBFs, p)/log(10);

    if(logbf < 1e-5) 
     sprintf(buf, "%.5f ", logbf); 
   else 
     sprintf(buf, "%+.5f ", logbf); 
   outfile << buf;
 }
 outfile << endl;



    // pi

 outfile_pi << vsRsnum.at(g) << " "; 
 
 for(int p =0; p < numG; p++){

  pi = pi_list[p];
  if(pi < 1e-5) 
   sprintf(buf, "%.5f ", pi); 
 else 
   sprintf(buf, "%+.5f ", pi); 
 outfile_pi << buf;
}
outfile_pi << endl;


    // phi, mean, and var
double phi, mean_out, var_out, mean1_out, var1_out;
outfile_mean << vsRsnum.at(g) << " ";  
outfile_var << vsRsnum.at(g) << " ";
    //-- wavelets_v2_2 start --//
outfile_mean1 << vsRsnum.at(g) << " ";  
outfile_var1 << vsRsnum.at(g) << " ";
outfile_phi << vsRsnum.at(g) << " ";
    //-- wavelets_v2_2 end --//

int p;
for(int gi =0; gi < numG; gi++){
  pi = pi_list[gi]; 
  int numP = nPH_GP[gi];
  int st = group_start[gi] - 1;
  for(int i = 0; i < numP; i++){
   p = st+i;
   if(use_pheno[p] == 1){
     double logpi  = log(pi);
     double log1pi = log(1-pi);
     double logPiBF = gsl_vector_get(logBFs,p) + logpi;
     double logden = sumlog(logPiBF, log1pi);
     phi = exp(logPiBF - logden);
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
}

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
outfile_eps.close(); 
cout << sfn_eps << " has been created." << endl;
//--- BL_HMT start ---//


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

//// ~~~~~ Remember to free up any memory for any new variables created.

logLR_list.resize(0);

}