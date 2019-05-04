// Just a working page for now. Integrate into BL_single_snp_hmt.cpp when done.

// Only need pi for the coarsest scale (gi == 0)
// Due to the nature of the Up-down algorithm, start at the finest scale (gi == numG)

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

	// double pi;

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