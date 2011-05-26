#include "host_hmc.h"

using namespace std;

//CP: molecular dynamics update for the gauge momenta:
//p_out = p_in - eps/2 force(u_in, phi)
//it is assumed that the force term has already been computed. then one only has real-vectors and this is essentially adding one vector to another...
hmc_error md_update_gauge_momenta(hmc_float eps, hmc_gauge_momentum * p_inout, hmc_gauge_momentum * force_in){
	for(int i = 0; i<GAUGEMOMENTASIZE; i++){
		p_inout[i] -= eps*force_in[i];
	}
	return HMC_SUCCESS;
}

//molecular dynamics update for the gaugefield:
//u_out = exp(i eps p_in) u_in
hmc_error md_update_gaugefield(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gaugefield * u_inout){
	int index;
	for(int t = 0; t<NTIME; t++){
		for(int pos = 0; pos < VOLSPACE; pos++){
			for(int mu = 0; mu<NDIM; mu++){
			hmc_su3matrix tmp;
			hmc_su3matrix tmp2;
			index= get_global_link_pos(mu, pos, t);
			// an su3 algebra element has NC*NC-1 = 8 hmc_float entries
			// &(p_in[index*8]) should point to the right position for the pos-th element of the long gaugemomentum vector p_in
			build_su3matrix_by_exponentiation(&(p_in[index*8]), &tmp2, eps);
			get_su3matrix(&tmp, u_inout, pos, t, mu);
			accumulate_su3matrix_prod( &tmp2, &tmp);
			put_su3matrix(u_inout, &tmp2, pos, t, mu);
	}}}
	return HMC_SUCCESS;
}

#ifdef _FERMIONS_
//phi = Q+ chi
hmc_error md_update_spinorfield(hmc_spinor_field * in, hmc_spinor_field * out, hmc_gaugefield * field, inputparameters * parameters){
	//TODO check again if it is M or Mdagger here
	M(parameters, in, field, out);
	gamma_5_psi(out);
	
	return HMC_SUCCESS;
}
#endif

// beta * sum_links sum_nu>mu ( 3 - Tr Re Plaquette )
//CP: since one is only interested in differences of s_gauge, the constant part can be left out!!
hmc_float s_gauge(hmc_gaugefield * field, hmc_float beta){
	/** @TODO CP: implement saving of plaquette measurement (and possibly t_plaq and s_plaq and also polyakov-loop??) */
	hmc_float plaq=0;

	//CP: obvious method
// 	for(int t=0;t<NTIME;t++) {
// 		for(int n=0;n<VOLSPACE;n++) {
// 			for(int mu=0; mu<NDIM; mu++) {
// 				for(int nu=0;nu<mu; nu++) {
// 					hmc_su3matrix prod;
// 					local_plaquette(field, &prod, n, t, mu, nu );
// 					hmc_float tmpfloat = - trace_su3matrix(&prod).re;
// 					plaq += tmpfloat;
// 				}}}}
// 	return beta/3.*plaq;
	
	//CP: alternative method: use already existing plaquette-functions
	hmc_float t_plaq;
	hmc_float s_plaq;
	plaq = plaquette(field, &t_plaq, &s_plaq);
	//plaq is normalized by factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC), so one has to divide by it again
	hmc_float factor = 2.0/static_cast<hmc_float>(VOL4D*NDIM*(NDIM-1)*NC);
	return beta/6./factor*( - plaq);
	
}

#ifdef _FERMIONS_
// sum_links phi*_i (M^+M)_ij^-1 phi_j
// here it is assumed that the rhs has already been computed... (because usually it has been during the leapfrog..)
hmc_complex s_fermion(hmc_spinor_field * phi, hmc_spinor_field * QplusQminusphi){
	return scalar_product(phi, QplusQminusphi);
}
#endif /* _FERMIONS_ */

//S_gauge + S_fermion + S_gaugemomenta
hmc_complex hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_gauge_momentum * p
												#ifdef _FERMIONS_
												, hmc_spinor_field * phi, hmc_spinor_field * phi_inv
												#endif
												){
	hmc_complex result;
	hmc_complex tmp;
	(result) = {0.,0.};
	(result).re += s_gauge(field, beta);
#ifdef _FERMIONS_
	tmp  = s_fermion(phi, phi_inv);
	complexaccumulate(&result, &tmp);
#endif
	//s_gm = 1/2*squarenorm(Pl)
	hmc_float s_gm;
	gaugemomenta_squarenorm(p, &s_gm);
	result.re += 0.5*s_gm;
	
	return result;
}

hmc_error gauge_force(inputparameters * parameters, hmc_gaugefield * field, hmc_gauge_momentum * out){
	hmc_float beta = (*parameters).get_beta();
	int globalpos;
	hmc_3x3matrix V;
	hmc_3x3matrix tmp;
	hmc_su3matrix U;
	
	//Gauge force is factor*Im(i Tr(T_i U V))
	//   with T_i being the SU3-Generator in i-th direction and V the staplematrix
	//   and the factor being -beta/3. (for standard Wilson-action)	
	for(int t = 0; t < NTIME; t++){
		for(int n = 0; n < VOLSPACE; n++){
			for(int mu = 0; mu < NDIM; mu++){
				globalpos = get_global_link_pos(mu, n, t);

				calc_staple(field, &V, n, t, mu);
				get_su3matrix(&U, field, n, t, mu);

				/** @TODO CP: this is not valid for REC12 */
				//when not using Reconstruct12 staplematrix and 3x3matrix are just the same!
				multiply_su3matrices (&tmp, &U, &V);

				hmc_algebraelement out_tmp;
				//CP: way one with function-calls
				//iterate through the different directions for i of which there are NC*NC-1 = 8
	 			//in the function the gen_index runs from 1 to 8 (including 8) !! -> i +1 
	 			//out is a vector of length NDIM*VOL4D*(NC*NC-1) = NDIM*VOL4D*8
// 				hmc_3x3matrix tmp2;
// 				for(int i = 0; i<8; i++){
// 					multiply_generator_3x3matrix (&tmp2, i+1, &tmp);
// 					trace_3x3matrix (&trace, &tmp2);
// 					//CP: Minus???
// 					(out_tmp)[i] = (trace.im);
// 				}
				//CP: hardcoded way: (like in tmlqcd)
				(out_tmp)[0]  = ( -(tmp)[1][0].im - (tmp)[0][1].im);
				(out_tmp)[1] = (+(tmp)[1][0].re-(tmp)[0][1].re);
				(out_tmp)[2] = (-(tmp)[0][0].im+(tmp)[1][1].im);
				(out_tmp)[3] = (-(tmp)[2][0].im-(tmp)[0][2].im);
				(out_tmp)[4] = (+(tmp)[2][0].re-(tmp)[0][2].re);
				(out_tmp)[5] = (-(tmp)[2][1].im-(tmp)[1][2].im);
				(out_tmp)[6] = (+(tmp)[2][1].re-(tmp)[1][2].re);
				(out_tmp)[7] = (-(tmp)[0][0].im-(tmp)[1][1].im + 2.0*(tmp)[2][2].im)*0.577350269189625;

				hmc_float factor = -beta/3.;
				for(int i = 0; i<8; i++){
  					out[globalpos*8 + i] = factor*out_tmp[i];
				}
			}
		}
	}
	
	return HMC_SUCCESS;
}

#ifdef _FERMIONS_

//CP: fermion_force = (gamma_5 Y)^dagger iT_i
//	it is assumed that the results can be added to out!!
hmc_error fermion_force(inputparameters * parameters, hmc_gaugefield * field, hmc_spinor_field * Y, hmc_spinor_field * X, hmc_gauge_momentum * out){

  hmc_su3matrix U_up, U_down;
	hmc_3x3matrix v1, v2;
  hmc_su3vector psia,psib,phia,phib;
	hmc_full_spinor y, plus;
	int nup, ndown;
	hmc_algebraelement out_tmp;
	int global_link_pos;
	int global_link_pos_down;
	hmc_float factor;
	
	//main loop
	for(int t = 0; t<NTIME; t++){
		for(int n = 0; n<VOLSPACE; n++){
			get_spinor_from_field(Y, y, n, t);
			///////////////////////////////////
			// Calculate gamma_5 y
			///////////////////////////////////
			gamma_5_spinor(y);

			//go through the different directions
			///////////////////////////////////
			// mu = +0
			///////////////////////////////////
			nup = (t+1)%NTIME;
			global_link_pos = get_global_link_pos(0, n, t);
			get_spinor_from_field(X, plus, n, nup);
			
			get_su3matrix(&U_up, field, n, t, 0);
		
			//psi = (1-gamma_mu)plus
			spinproj_gamma0_a(plus, psia, -hmc_one_f);
			spinproj_gamma0_b(plus, psib, -hmc_one_f);
      
			//phi = (1-gamma_mu)y
			spinproj_gamma0_a(y, phia, -hmc_one_f);
			spinproj_gamma0_b(y, phib, -hmc_one_f);

			// v1 = Tr(phi*psi_dagger)
			tr_v_times_u_dagger(phia, psia, phib, psib, v1);

			//U*v1 = U*(phi_a)
			/** @todo what about REC12 and this call??*/
			multiply_3x3matrix (&v2, &U_up, &v1);
		
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			tr_lambda_u(v2, out_tmp);
			//what is the factor here??
			factor = 1.;
			for(int i = 0; i<8; i++){
  					out[global_link_pos*8 + i] += factor*out_tmp[i];
			}

			///////////////////////////////////
			//mu = -0
			///////////////////////////////////
			ndown = (t+NTIME-1)%NTIME;
			global_link_pos_down = get_global_link_pos(0, n, ndown);
			get_spinor_from_field(X, plus, n, ndown);
			
			get_su3matrix(&U_down, field, n, ndown, 0);
		
			//psi = (1+gamma_mu)plus
			spinproj_gamma0_a(plus, psia, hmc_one_f);
			spinproj_gamma0_b(plus, psib, hmc_one_f);
      
			//phi = (1+gamma_mu)y
			spinproj_gamma0_a(y, phia, hmc_one_f);
			spinproj_gamma0_b(y, phib, hmc_one_f);

			//CP: here is the difference with regard to +mu-direction: psi and phi interchanged!!
			// v1 = Tr(psi*phi_dagger)
			tr_v_times_u_dagger(psia, phia, psib, phib, v1);

			//U*v1 = U*(phi_a)
			/** @todo what about REC12 and this call??*/
			multiply_3x3matrix (&v2, &U_down, &v1);
		
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			
			//this must become v1 if the above is included again
			tr_lambda_u(v2, out_tmp);
			//what is the factor here??
			factor = (*parameters).get_kappa();
			for(int i = 0; i<8; i++){
  					out[global_link_pos_down*8 + i] += factor*out_tmp[i];
			}
		}}
	return HMC_SUCCESS;
}
#endif

//CP: this essentially calculates a hmc_gauge_momentum vector
//CP: if fermions are used, here is the point where the inversion has to be performed
hmc_error force(inputparameters * parameters, hmc_gaugefield * field
#ifdef _FERMIONS_
	, hmc_spinor_field * phi, hmc_spinor_field * phi_inv
#endif
	, hmc_gauge_momentum * out){
	cout << "\t\tstart calculating the force..." << endl;
	//CP: make sure that the output field is set to zero
	set_zero_gaugemomenta(out);
	//add contributions
	cout << "\t\tcalc gauge_force..." << endl;
	gauge_force(parameters, field, out);
#ifdef _FERMIONS_
	cout << "\t\tinvert fermion field..." << endl;
	//the algorithm needs two spinor-fields
	hmc_spinor_field* X = new hmc_spinor_field[SPINORFIELDSIZE];
	//CP: to begin with, consider only the cg-solver
	//source is at 0
	int k = 0;
	int use_cg = TRUE;
	//CP: at the moment, use_eo = 0 so that even-odd is not used!!!!!
	
	//debugging
	int err = 0;
	
	if(use_cg){
		if(!use_eo){
			//the inversion calculates Y = (QplusQminus)^-1 phi = phi_inv
			hmc_spinor_field b[SPINORFIELDSIZE];
			create_point_source(parameters,k,0,0,b);
			cout << "\t\t\tstart solver" << endl;
			err = solver(parameters, phi, b, field, use_cg, phi_inv);
		}
		else{
			hmc_eoprec_spinor_field be[EOPREC_SPINORFIELDSIZE];
			hmc_eoprec_spinor_field bo[EOPREC_SPINORFIELDSIZE];
			
			create_point_source_eoprec(parameters, k,0,0, field, be,bo);
			solver_eoprec(parameters, phi, be, bo, field, use_cg, phi_inv);
		}
		if (err != HMC_SUCCESS) cout << "\t\tsolver did not solve!!" << endl;
		cout << "\t\t\tcalc X" << endl;
		//X = Qminus Y = Qminus phi_inv 
		Qminus(parameters, phi_inv, field, X);
	}
	else{
		//here, one has first to invert (with BiCGStab) Qplus phi = X and then invert Qminus X => Qminus^-1 Qplus^-1 phi = (QplusQminus)^-1 phi = Y = phi_inv
	}

/** @todo control the fields here again!!! */
	fermion_force(parameters, field, phi_inv, X, out);
	
	delete [] X;
#endif
	return HMC_SUCCESS;
}

hmc_error metropolis(hmc_float rndnumber, hmc_float beta, 
										 #ifdef _FERMIONS_
										 hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi, 
										 #endif
										 hmc_gaugefield * field,	hmc_gauge_momentum * p, hmc_gaugefield * new_field, hmc_gauge_momentum * new_p){
	// takes:
	//		phi and beta as constant
	//		new/old versions of gaugefield and of momenta
	//		and a random number
	//		if it has to be, performs the change old->new, and returns true if there are no failures.
	hmc_complex h_old = hamiltonian(field, beta, p
																	#ifdef _FERMIONS_
																	, phi, MdaggerMphi
																	#endif
																	);
	hmc_complex h_new = hamiltonian(new_field, beta, new_p
																	#ifdef _FERMIONS_
																	, phi, MdaggerMphi
																	#endif
																	);
																	
																	
	if(h_old.im > projectioneps){
		printf("\n\tError: imaginary part in H_OLD [in function: metropolis(...)].\n");
		return HMC_COMPLEX_HAMILTONIANERROR;
	}
	if(h_new.im > projectioneps){
		printf("\n\tError: imaginary part in H_NEW [in function: metropolis(...)].\n");
		return HMC_COMPLEX_HAMILTONIANERROR;
	}
	/** @todo CP:  export h_diff */
	hmc_float h_diff = h_old.re - h_new.re;
	hmc_float compare_prob;
	if(h_diff<0){
		compare_prob = exp(h_diff);
	}else{
		compare_prob = 1.0;
	}
	// SL: the following can be tuned, whether it is more costly to draw always the rnd number even when compare_prob=1
	//     and whether the "if compare_prob==1" costs more or less than always evaluating the exp ...
	if(rndnumber <= compare_prob){
		// perform the change nonprimed->primed !
		copy_gaugefield(new_field, field);
		copy_gaugemomenta(new_p, p);
		// SL: this works as long as p and field are pointers to the *original* memory locations!
		cout << "new configuration accepted" << endl;
	}
	else{
		cout << "new configuration rejected" << endl;
	}
	return HMC_SUCCESS;
}

//it is assumed that gaugefield and gaugemomentum have been set to the old ones already
hmc_error leapfrog(inputparameters * parameters, 
									 #ifdef _FERMIONS_
									 hmc_spinor_field * phi, hmc_spinor_field * phi_inv, 
									 #endif
									 hmc_gaugefield * u_out, hmc_gauge_momentum * p_out	){
	// CP: it operates directly on the fields p_out and u_out
	int steps = (*parameters).get_integrationsteps1() ;	
	hmc_float stepsize = ((*parameters).get_tau()) /((hmc_float) steps);
	int k;
	hmc_float stepsize_half = 0.5*stepsize;
	
	hmc_gauge_momentum* force_tmp = new hmc_gauge_momentum[GAUGEMOMENTASIZE];

	//initial step
	cout << "\tinitial step:" << endl;
	force(parameters, u_out ,
		#ifdef _FERMIONS_
		phi, phi_inv, 
		#endif
		force_tmp);
	md_update_gauge_momenta(stepsize_half, p_out, force_tmp);
	
	//intermediate steps
	if(steps > 1) cout << "\tperform " << steps << " intermediate steps " << endl;
	for(k = 1; k<steps; k++){
		md_update_gaugefield(stepsize, p_out, u_out);
		force(parameters, u_out ,
			#ifdef _FERMIONS_
			phi, phi_inv, 
			#endif
			force_tmp);
		md_update_gauge_momenta(stepsize, p_out, force_tmp);
	}
	
	//final step
	cout << "\tfinal step" << endl;
	md_update_gaugefield(stepsize, p_out, u_out);
	force(parameters, u_out ,
		#ifdef _FERMIONS_
		phi, phi_inv, 
		#endif
		force_tmp);
	md_update_gauge_momenta(stepsize_half, p_out,force_tmp); 
	
	delete [] force_tmp;
	
	cout << "\tfinished leapfrog" << endl;
	return HMC_SUCCESS;
}

