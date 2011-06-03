/** @file
 * Device code for the heatbath update
 */

//opencl_update_heatbath.cl

Matrix3x3 calc_staple(__global hmc_ocl_gaugefield* field, const int pos, const int t, const int mu_in)
{
	Matrixsu3 prod;
	Matrixsu3 prod2;
	Matrixsu3 tmp;
	Matrix3x3 staple;
	int nu, newpos, newt;

 	staple = zero_matrix3x3();
	
	//iterate through the three directions other than mu
 	for(int i = 1; i<NDIM; i++) {

		prod = zero_matrixsu3();
		
		nu = (mu_in + i)%NDIM;
		//first staple
		//u_nu(x+mu)
		if(mu_in==0) {
			newt = (t+1)%NTIME;
			tmp = get_matrixsu3(field,pos,newt,nu);

		} else {
			tmp = get_matrixsu3(field,get_neighbor(pos,mu_in),t,nu);
		}
		prod = copy_matrixsu3(tmp);
		
		//adjoint(u_mu(x+nu))
		if(nu==0) {
			newt = (t+1)%NTIME;
			tmp = get_matrixsu3(field,pos,newt,mu_in);
		} else {
			tmp = get_matrixsu3(field,get_neighbor(pos,nu),t,mu_in);
		}
		tmp = adjoint_matrixsu3(tmp);
		
		prod = multiply_matrixsu3(prod, tmp);
		
		//adjoint(u_nu(x))
		tmp = get_matrixsu3(field,pos,t,nu);
		tmp = adjoint_matrixsu3(tmp);
		prod = multiply_matrixsu3 (prod, tmp);
		
		//second staple
		//adjoint (u_nu(x+mu-nu))
		//newpos is "pos-nu" (spatial)
		newpos = get_lower_neighbor(pos, nu);
		if(mu_in==0) {
			newt = (t+1)%NTIME;
			tmp = get_matrixsu3(field,newpos,newt,nu);
		} else if (nu == 0) {
			newt = (t-1+NTIME)%NTIME;
			tmp = get_matrixsu3(field,get_neighbor(pos,mu_in),newt,nu);
		} else {
			tmp = get_matrixsu3(field,get_neighbor(newpos,mu_in),t,nu);
		}
		prod2 = adjoint_matrixsu3(tmp);
		//adjoint(u_mu(x-nu))
		if(mu_in==0) {
			tmp = get_matrixsu3(field,newpos,t,mu_in);
		} else if (nu == 0) {
			newt = (t-1+NTIME)%NTIME;
			tmp = get_matrixsu3(field,pos,newt,mu_in);
		} else {
			tmp = get_matrixsu3(field,newpos,t,mu_in);
		}
		tmp = adjoint_matrixsu3(tmp);
		prod2 = multiply_matrixsu3(prod2, tmp);
		//adjoint(u_nu(x-nu))
		if(mu_in==0) {
			tmp = get_matrixsu3(field,newpos,t,nu);
		} else if (nu == 0) {
			newt = (t-1+NTIME)%NTIME;
			tmp = get_matrixsu3(field,pos,newt,nu);
		} else {
			tmp = get_matrixsu3(field,newpos,t,nu);
		}
		prod2 = multiply_matrixsu3(prod2, tmp);

		
		Matrix3x3 dummy;
		dummy = matrix_su3to3x3 (prod);
		staple = add_matrix3x3 (staple, dummy );
		dummy = matrix_su3to3x3 (prod2);
		staple = add_matrix3x3 (staple, dummy );
	}
	
	return staple;
}

void SU2Update(__private hmc_float dst [su2_entries], const hmc_float alpha, __global hmc_ocl_ran * rnd)
{
	hmc_float delta;
	hmc_float a0 ;
	hmc_float eta ;
	do {
		delta = -log(ocl_new_ran(rnd))/alpha*pow(cos(2. * PI * ocl_new_ran(rnd)), 2.) -log(ocl_new_ran(rnd))/alpha;
		a0 = 1.-delta;
		eta = ocl_new_ran(rnd);
	} while ( (1.-0.5*delta) < eta*eta);
	hmc_float phi = 2.*PI*ocl_new_ran(rnd);
	hmc_float theta = asin(2.*ocl_new_ran(rnd) - 1.);
	dst[0] = a0;
	dst[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
	dst[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
	dst[3] = sqrt(1.-a0 * a0)*sin(theta);
}

void inline perform_heatbath(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd, int pos, int t, int id)
{
	Matrixsu3 U;
	Matrix3x3 W;
	Matrix3x3 staplematrix;
	int order[3];
	hmc_complex w [su2_entries];
	hmc_float w_pauli[su2_entries];
	hmc_float k;
	hmc_float r_pauli[su2_entries];
	hmc_float beta_new;

	random_1_2_3(order, &rnd[id]);
	
	U = get_matrixsu3(gaugefield, pos, t, mu);
	
	staplematrix = calc_staple(gaugefield, pos, t, mu);
	
	for(int i=0; i<NC; i++) {
	  	
		W = matrix_su3to3x3 (U);
		W = multiply_matrix3x3 (W, staplematrix);
		
 		reduction(w, W, order[i]);

		w_pauli[0] = 0.5*(w[0].re + w[3].re);
		w_pauli[1] = 0.5*(w[1].im + w[2].im);
		w_pauli[2] = 0.5*(w[1].re - w[2].re);
		w_pauli[3] = 0.5*(w[0].im - w[3].im);
		k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

		beta_new =  2.*beta / NC*k;
		SU2Update(r_pauli, beta_new, &rnd[id]);

		//dispensable?
		/*
		w[0].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
		w[0].im = (w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
		w[1].re = (w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
		w[1].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
		w[2].re = -(w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
		w[2].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
		w[3].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
		w[3].im = -(w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
		*/

		//old:
		w_pauli[0] = w_pauli[0]/k;
		w_pauli[1] = -w_pauli[1]/k;
		w_pauli[2] = -w_pauli[2]/k;
		w_pauli[3] = -w_pauli[3]/k;

		hmc_float su2_tmp[su2_entries];
		su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
		su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
		su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
		su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
		r_pauli[0] = su2_tmp[0];
		r_pauli[1] = su2_tmp[1];
		r_pauli[2] = su2_tmp[2];
		r_pauli[3] = su2_tmp[3];

		//go back to a su2 matrix in standard basis
		w[0].re = r_pauli[0];
		w[0].im = r_pauli[3];
		w[1].re = r_pauli[2];
		w[1].im = r_pauli[1];
		w[2].re = -r_pauli[2];
		w[2].im = r_pauli[1];
		w[3].re = r_pauli[0];
		w[3].im = -r_pauli[3];

		Matrixsu3 extW;
		extW = extend (order[i], w);
		
		//Matthias
		//für order[i]=2, 3 gibt es endliche ergebniss
		//für order[i]=1 unendliche
		//Die Funktion extend macht aber das was sie sollte...
		//
		//extW = extend (3, w);

		//Matthias
		//Schreibt die Elemente w aus
		// !!! wird dieser Teil auskommentiert, gehen die Einträge von U gegen unendlich
		//
		//printf("order %i: \n", order[i]);
 		//printf("%f %f \t %f %f \t %f %f \t %f %f \n", w[0].re, w[0].im, w[1].re, w[1].im, w[2].re, w[2].im, w[3].re, w[3].im);

		U = multiply_matrixsu3 (extW, U);
	}
	
	project_su3(U);

	put_matrixsu3(gaugefield, U, pos, t, mu);

	//Matthias
	//Schreibt die erzeugte Matrix aus
// 	printf("%f %f \t %f %f \t %f %f \n",U.e00.re, U.e00.im, U.e01.re, U.e01.im, U.e02.re, U.e02.im);
// 	printf("%f %f \t %f %f \t %f %f \n",U.e10.re, U.e10.im, U.e11.re, U.e11.im, U.e12.re, U.e12.im);
// 	printf("%f %f \t %f %f \t %f %f \n",U.e20.re, U.e20.im, U.e21.re, U.e21.im, U.e22.re, U.e22.im);
// 	printf("\n");
	

	//Matthias
	//Überprüft, ob die erzeugte Matrix unitär ist
	//Ja, falls trace.re = 3.0 und trace.im = 0.0
// 	Matrixsu3 blubb;
// 	Matrixsu3 adjU;
// 	adjU = adjoint_matrixsu3(U);
// 	blubb = multiply_matrixsu3 (U, adjU);
// 	hmc_complex trace;
// 	trace = trace_matrixsu3 (blubb);
// 	printf (" U * adj(u) %f \n", trace.re);
// 	printf (" U * adj(u) %f \n", trace.im);
	
}



__kernel void heatbath_even(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd)
{
	int t, pos, id, id_tmp, size;
	id_tmp = get_global_id(0);
	size = get_global_size(0);
	for(id = id_tmp; id<VOLSPACE*NTIME/2; id+=size) {
		get_even_site(id, &pos, &t);
		perform_heatbath(gaugefield, beta, mu, rnd, pos, t, id_tmp);
	}
	return;
}

__kernel void heatbath_odd(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd)
{
	int t, pos, id, id_tmp, size;
	id_tmp = get_global_id(0);
	size = get_global_size(0);
	for(id = id_tmp; id<VOLSPACE*NTIME/2; id+=size) {
		get_odd_site(id, &pos, &t);
		perform_heatbath(gaugefield, beta, mu, rnd, pos, t, id_tmp);
	}
	return;
}

void inline perform_overrelaxing(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd, int pos, int t, int id)
{

	Matrixsu3 U;
	Matrix3x3 W;
	Matrix3x3 staplematrix;

	hmc_complex w [su2_entries];
	hmc_float w_pauli[su2_entries];
	hmc_float k;
	int order[3];

	random_1_2_3(order, &rnd[id]);
	U = get_matrixsu3(gaugefield, pos, t, mu);

//why?
	project_su3(U);

	staplematrix = calc_staple(gaugefield, pos, t, mu);

	Matrixsu3 extW;

	for(int i=0; i<NC; i++) {
		W = matrix_su3to3x3 (U);
		W = multiply_matrix3x3 (W, staplematrix);
	  
		reduction(w, W, order[i]);

		w_pauli[0] = 0.5*(w[0].re + w[3].re);
		w_pauli[1] = 0.5*(w[1].im + w[2].im);
		w_pauli[2] = 0.5*(w[1].re - w[2].re);
		w_pauli[3] = 0.5*(w[0].im - w[3].im);
		k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

		w[0].re = (w_pauli[0]*w_pauli[0] - w_pauli[1]*w_pauli[1] - w_pauli[2]*w_pauli[2] - w_pauli[3]*w_pauli[3])/k/k;
		w[0].im = (-2.*w_pauli[0]*w_pauli[3])/k/k;
		w[1].re = (-2.*w_pauli[0]*w_pauli[2])/k/k;
		w[1].im = (-2.*w_pauli[0]*w_pauli[1])/k/k;
		w[2].re = (2.*w_pauli[0]*w_pauli[2])/k/k;
		w[2].im = (-2.*w_pauli[0]*w_pauli[1])/k/k;
		w[3].re = (w_pauli[0]*w_pauli[0] - w_pauli[1]*w_pauli[1] - w_pauli[2]*w_pauli[2] - w_pauli[3]*w_pauli[3])/k/k;
		w[3].im = (2.*w_pauli[0]*w_pauli[3])/k/k;

		extW = extend (order[i], w);
		U = multiply_matrixsu3(extW, U);
	}
	put_matrixsu3(gaugefield, U, pos, t, mu);

	return;
}

__kernel void overrelax_even(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd)
{
	int t, pos, id, id_tmp, size;
	id_tmp = get_global_id(0);
	size = get_global_size(0);
	for(id = id_tmp; id<VOLSPACE*NTIME/2; id+=size) {
		get_even_site(id, &pos, &t);
		perform_overrelaxing(gaugefield, beta, mu, rnd, pos, t, id_tmp);
	}
	return;
}

__kernel void overrelax_odd(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd)
{
	int t, pos, id, id_tmp, size;
	id_tmp = get_global_id(0);
	size = get_global_size(0);
	for(id = id_tmp; id<VOLSPACE*NTIME/2; id+=size) {
		get_odd_site(id, &pos, &t);
		perform_overrelaxing(gaugefield, beta, mu, rnd, pos, t, id_tmp);
	}
	return;
}

