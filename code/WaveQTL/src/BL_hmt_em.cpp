//--- BL_HMT start ---//
// Just a working page for now. Integrate into BL_single_snp_hmt.cpp when done.

// Only need pi for the coarsest scale (gi == 0)
// Due to the nature of the Up-down algorithm, start at the finest scale (gi == numG)
    /**********************/
    //     EM algorithm    //
    /***********************/

    double O_obllikli;
    double N_obllikli;
    pi_list.resize(0);

    double logLR = 0;
    double pi;

    // With HMT, need to do the EM algorithm for all groups (scales) in one step.
    // Can no longer do separately as they are all linked together.
    // ~~~ TO CONFIRM: What does the obl log-l/hood look like in this version?
    // ~~~ TO CONFIRM: How are we handling multiple individuals? Are they separate, or same trees?
    // ~~~ TO CONFIRM: nCohort same as nIndiv?
    // Would each individuals' tree posterior probabilities be multiplicative (indp't)?

    // ---- Initialisation of parameters ---- //

	// Alpha and beta quantities
    gsl_matrix * beta_sl_1 = gsl_matrix_alloc(nIndiv, nPH)
    gsl_matrix * beta_sl_0 = gsl_matrix_alloc(nIndiv, nPH)

    gsl_matrix * beta_sl_psl_1 = gsl_matrix_alloc(nIndiv, nPH)
    gsl_matrix * beta_sl_psl_0 = gsl_matrix_alloc(nIndiv, nPH)

    gsl_matrix * beta_psl_no_sl_1 = gsl_matrix_alloc(nIndiv, nPH)
    gsl_matrix * beta_psl_no_sl_0 = gsl_matrix_alloc(nIndiv, nPH)

    gsl_matrix * alpha_sl_1 = gsl_matrix_alloc(nIndiv, nPH)
    gsl_matrix * alpha_sl_0 = gsl_matrix_alloc(nIndiv, nPH)

    gsl_matrix * pp_mtx = gsl_matrix_alloc(nIndiv, nPH)
    gsl_matrix * pp_joint_11_mtx = gsl_matrix_alloc(nIndiv, nPH)
    gsl_matrix * pp_joint_10_mtx = gsl_matrix_alloc(nIndiv, nPH)
    gsl_matrix * pp_joint_01_mtx = gsl_matrix_alloc(nIndiv, nPH)

    // ---- Initial E-step ---- //
    // 'Zero-th' up-down algorithm for each individual
    for(int ind = 0, ind < nIndiv, ind++){

		// ~~~ This pi actually only needed for coarsest scale (init of alpha at 1,1)
		double logpi = log(0.5);
		double log1pi = logpi;

		// Probably (?) easier to also have epsilons as logs
		double logeps_11 = log(0.25);
		double logeps_01 = log(0.25);
		double logeps_10 = log(0.25);
		double logeps_00 = log(0.25);

    	// Up step
    	for(int gi = (numG - 1); gi > 0; gi--){

			int numP = nPH_GP[gi];
			int st = group_start[gi] - 1;
			int parent_st = group_start[gi - 1] - 1;

			// Step 0 at finest scale only.
			if(gi == (numG - 1)){
				for(int wc = 0; wc < numP; wc++){

			    	double logBF;
			    	logBF = gsl_vector_get(logBFs,st+wc);

			    	double b_sl_1 = exp(logBF);
			    	double b_sl_0 = 1;

					gsl_matrix_set(beta_sl_1, ind, st+wc, b_sl_1);
					gsl_matrix_set(beta_sl_0, ind, st+wc, b_sl_0);
				}
			}

			for(int wc = 0; wc < numP; wc++){

		    	double b_sl_1;
		    	double b_sl_0;
		    	// The beta of the adjacent child (which shares the same parent)
				double b_sl_adj_1;
		    	double b_sl_adj_0;

		    	double b_sl_psl_1;
		    	double b_sl_psl_0;
		    	double b_psl_1;
		    	double b_psl_0;
		    	double b_psl_no_sl_1;
		    	double b_psl_no_sl_0;

		    	double logParentBF;
		    	int parent_indx = wc/2;
		    	logParentBF = gsl_vector_get(logBFs,parent_st+parent_indx);


				b_sl_1 = gsl_matrix_get(beta_sl_1, ind, st+wc);
				b_sl_0 = gsl_matrix_get(beta_sl_0, ind, st+wc);

				// If even, grab one to the right, else one to the left
				if(wc % 2 == 0){
					b_sl_adj_1 = gsl_matrix_get(beta_sl_1, ind, st+wc + 1);
					b_sl_adj_0 = gsl_matrix_get(beta_sl_0, ind, st+wc + 1);
				}else{
					b_sl_adj_1 = gsl_matrix_get(beta_sl_1, ind, st+wc - 1);
					b_sl_adj_0 = gsl_matrix_get(beta_sl_0, ind, st+wc - 1);
				}

				b_sl_psl_1 = b_sl_1*exp(logeps_11) + b_sl_1*exp(logeps_01);
				b_sl_psl_0 = b_sl_0*exp(logeps_00) + b_sl_0*exp(logeps_10);

				b_psl_1 = exp(logParentBF) * b_sl_adj_1 * b_sl_1;
				b_psl_0 = exp(logParentBF) * b_sl_adj_0 * b_sl_0;

				b_psl_no_sl_1 = b_psl_1/b_sl_psl_1;
				b_psl_no_sl_0 = b_psl_0/b_sl_psl_0;

				gsl_matrix_set(beta_sl_psl_1, ind, st+wc, b_sl_psl_1);
				gsl_matrix_set(beta_sl_psl_0, ind, st+wc, b_sl_psl_0);
				
				// Both the children will try and update the parents' beta, so
				// only do this once.
				if(wc % 2 == 0){
					gsl_matrix_set(beta_sl_1, ind, st+wc, b_psl_1);
					gsl_matrix_set(beta_sl_0, ind, st+wc, b_psl_0);	
				}

				gsl_matrix_set(beta_psl_no_sl_1, ind, st+wc, b_psl_no_sl_1);
				gsl_matrix_set(beta_psl_no_sl_0, ind, st+wc, b_psl_no_sl_0);

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

		    	double a_sl_1 = exp(logpi);
	    		double a_sl_0 = exp(logpi);

				gsl_matrix_set(alpha_sl_1, ind, st, a_sl_1);
				gsl_matrix_set(alpha_sl_0, ind, st, a_sl_0);
			}else{

				for(int wc = 0; wc < numP; wc++){

			    	double b_psl_no_sl_1;
			    	double b_psl_no_sl_0;
			    	b_psl_no_sl_1 = gsl_matrix_get(beta_sl_psl_1, ind, st+wc);
			    	b_psl_no_sl_0 = gsl_matrix_get(beta_sl_psl_0, ind, st+wc);

			    	int parent_indx = wc/2;
			    	double a_psl_1 = gsl_matrix_get(alpha_sl_1, ind, st+parent_indx);
			    	double a_psl_0 = gsl_matrix_get(alpha_sl_0, ind, parent_st+parent_indx);

			    	double a_sl_1 = exp(logeps_11)*a_psl_1*b_psl_no_sl_1 + exp(logeps_10)*a_psl_0*b_psl_no_sl_0;
			    	double a_sl_0 = exp(logeps_01)*a_psl_1*b_psl_no_sl_1 + exp(logeps_00)*a_psl_0*b_psl_no_sl_0;
			    	gsl_matrix_set(alpha_sl_1, ind, st+wc, a_sl_1);
			    	gsl_matrix_set(alpha_sl_0, ind, st+wc, a_sl_0);

	    		}
			}
    		
    	}

		// ---- Calculate the E-step quantities (posterior probabilities) ---- //
	    for(int wc = 0; wc < nPH; wc++){

	    	/// ~~~ FIX: joint matrix NOT DEFINED for (1,1).
		    
	    	double b_sl_1 = gsl_matrix_get;
	    	double b_sl_0;
	    	double a_sl_1;
	    	double a_sl_0;

	    	double a_psl_1;
	    	double a_psl_0;
	    	double b_psl_1;
	    	double b_psl_0;	    	

	    	double denom = b_sl_1*a_sl_1 + b_sl_0*a_sl_0;
	    	double pp = (b_sl_1*a_sl_1)/denom;
	    	double pp_joint_11 = (a_psl_1*b_psl_1*b_sl_1*exp(logeps_11))/denom;
	    	double pp_joint_10 = (a_psl_0*b_psl_0*b_sl_1*exp(logeps_10))/denom;
	    	double pp_joint_01 = (a_psl_1*b_psl_1*b_sl_0*exp(logeps_01))/denom;

			gsl_matrix_set(pp_mtx, ind, wc, pp);
			gsl_matrix_set(pp_joint_11_mtx, ind, wc, pp_joint_11);
			gsl_matrix_set(pp_joint_10_mtx, ind, wc, pp_joint_10);
			gsl_matrix_set(pp_joint_01_mtx, ind, wc, pp_joint_01);		    
			
	    }

    }
	

    for(int gi = 0; gi < numG; g++){
    	N_obllikli = 0;
    	for(int gi = 0; gi < (numG - 1); gi++){

    	}
    }




	// ---- Iterate through first M-step, and subsequent E-steps ---- //

	// ---- Print things to outfiles ---- //

	// ---- Compile estimates of phi, means and vars ---- //

    for(int gi = 0; gi < numG; gi++){

      // EM algorithm for each group separately 


      
      
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

//--- BL_HMT end ---/