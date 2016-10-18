/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
  Calculation of S_det = Log(Det[(1+T_ee)^2])
*/

hmc_float log_det_squared_qr(Matrix6x6 a)
{
    //map Matrix6x6 struct to C array6x6
    const int rows = 6; const int cols = 6;
    hmc_complex T[6][6];
    //0.row
    T[0][0].re = a.e00.re; T[0][0].im = a.e00.im;
    T[0][1].re = a.e01.re; T[0][1].im = a.e01.im;
    T[0][2].re = a.e02.re; T[0][2].im = a.e02.im;
    T[0][3].re = a.e03.re; T[0][3].im = a.e03.im;
    T[0][4].re = a.e04.re; T[0][4].im = a.e04.im;
    T[0][5].re = a.e05.re; T[0][5].im = a.e05.im;
    //1.row
    T[1][0].re = a.e10.re; T[1][0].im = a.e10.im;
    T[1][1].re = a.e11.re; T[1][1].im = a.e11.im;
    T[1][2].re = a.e12.re; T[1][2].im = a.e12.im;
    T[1][3].re = a.e13.re; T[1][3].im = a.e13.im;
    T[1][4].re = a.e14.re; T[1][4].im = a.e14.im;
    T[1][5].re = a.e15.re; T[1][5].im = a.e15.im;
    //2.row
    T[2][0].re = a.e20.re; T[2][0].im = a.e20.im;
    T[2][1].re = a.e21.re; T[2][1].im = a.e21.im;
    T[2][2].re = a.e22.re; T[2][2].im = a.e22.im;
    T[2][3].re = a.e23.re; T[2][3].im = a.e23.im;
    T[2][4].re = a.e24.re; T[2][4].im = a.e24.im;
    T[2][5].re = a.e25.re; T[2][5].im = a.e25.im;
    //3.row
    T[3][0].re = a.e30.re; T[3][0].im = a.e30.im;
    T[3][1].re = a.e31.re; T[3][1].im = a.e31.im;
    T[3][2].re = a.e32.re; T[3][2].im = a.e32.im;
    T[3][3].re = a.e33.re; T[3][3].im = a.e33.im;
    T[3][4].re = a.e34.re; T[3][4].im = a.e34.im;
    T[3][5].re = a.e35.re; T[3][5].im = a.e35.im;
    //4.row
    T[4][0].re = a.e40.re; T[4][0].im = a.e40.im;
    T[4][1].re = a.e41.re; T[4][1].im = a.e41.im;
    T[4][2].re = a.e42.re; T[4][2].im = a.e42.im;
    T[4][3].re = a.e43.re; T[4][3].im = a.e43.im;
    T[4][4].re = a.e44.re; T[4][4].im = a.e44.im;
    T[4][5].re = a.e45.re; T[4][5].im = a.e45.im;
    //5.row
    T[5][0].re = a.e50.re; T[5][0].im = a.e50.im;
    T[5][1].re = a.e51.re; T[5][1].im = a.e51.im;
    T[5][2].re = a.e52.re; T[5][2].im = a.e52.im;
    T[5][3].re = a.e53.re; T[5][3].im = a.e53.im;
    T[5][4].re = a.e54.re; T[5][4].im = a.e54.im;
    T[5][5].re = a.e55.re; T[5][5].im = a.e55.im;

	//compute the determinant of a 6x6 Matrix via QR-decomposition
	//build R_1,...,R_n-1,T via vectors u_k
    hmc_complex u[6];
    hmc_complex N1[6][6];
    hmc_complex R[6-1][6][6] = {0.};

	//square input matrix
    for(unsigned int i=0; i<rows; ++i){
        for(unsigned int j=0; j<cols; ++j){
            N1[i][j] = hmc_complex_zero;
            for(unsigned int l=0; l<rows; ++l){
                N1[i][j] = complex_add(N1[i][j], complex_mult(T[i][l], T[l][j]));
            }}}
    for(unsigned int i=0; i<rows; ++i){
        for(unsigned int j=0; j<cols; ++j){
            T[i][j] = N1[i][j];}}
    
    for(int k=0; k<(cols-1); ++k){
        //build u_k and norm^2 of u_k
        hmc_float norm_u_squared = 0.;
        for(int l=0; l<rows; ++l){ //vector u according to Lüscher-OpenQCD doc
            if(l<k){u[l] = hmc_complex_zero;}
            else if(l==k){hmc_float r = 0.0;
                for(unsigned int j=k; j<rows; ++j){r = r + complex_abs_value(T[j][k]) * complex_abs_value(T[j][k]);}
                if(T[l][k].re == 0 && T[l][k].im == 0){u[l].re = - sqrt(r);} // 0/abs(0)!=1
                else{u[l] = complex_add(T[k][k], complex_mult(complex_divid(T[k][k], convert_float_to_complex(complex_abs_value(T[k][k]))), convertfloattocomplex(sqrt(r))));}}
            else{u[l] = T[l][k];}
            //printf("(%f,%f)", u[l]);
            norm_u_squared = norm_u_squared + complex_abs_value(u[l]) * complex_abs_value(u[l]);}

    // build up R_k by u_k und perform R_k * R_k-1*...*R_1*T
    for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
            R[k][m][n] = complex_mult(convert_float_to_complex(-2./norm_u_squared), complex_mult(u[m], complex_conj(u[n])));
            if(m==n){R[k][m][n] = complex_add(hmc_complex_one, R[k][m][n]);}}}
    for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
            N1[m][n] = hmc_complex_zero;
            for(unsigned int l=0; l<rows; ++l){
                N1[m][n] = complex_add(N1[m][n], complex_mult(R[k][m][l], T[l][n]));
            }}}
    for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
            T[m][n] = N1[m][n];}}
    }//end for loop over k-->now we have R_n-1,...,R_1,T=upper triangular
	/*for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
	  		printf("(%f,%f)", T[m][n].re, T[m][n].im);}
		printf("\n");}*/
	
	hmc_complex det = hmc_complex_one;
	//compute determinant of T
	for(unsigned int i=0; i<rows; ++i){
		//printf("%f", T[i][i]);
		//printf("\n");
		det = complex_mult(det, T[i][i]);
	//printf("(%f,%f)", det.re, det.im);
	//printf("\n");
	}
	//printf("(%f,%f)", det.re, det.im);
	//printf("\n");
	det.re *= -1; //^=(-1)^k with k number of householder reflections
	/* logarithm */
    hmc_float result = log(det.re);
	//printf("%f", result);
	return result;
}

hmc_float log_det_squared(Matrix6x6 a)
{
    //map Matrix6x6 struct to C array6x6
    const int n = 5;
    hmc_complex T[6][6];
    //0.row
    T[0][0].re = a.e00.re; T[0][0].im = a.e00.im;
    T[0][1].re = a.e01.re; T[0][1].im = a.e01.im;
    T[0][2].re = a.e02.re; T[0][2].im = a.e02.im;
    T[0][3].re = a.e03.re; T[0][3].im = a.e03.im;
    T[0][4].re = a.e04.re; T[0][4].im = a.e04.im;
    T[0][5].re = a.e05.re; T[0][5].im = a.e05.im;
    //1.row
    T[1][0].re = a.e10.re; T[1][0].im = a.e10.im;
    T[1][1].re = a.e11.re; T[1][1].im = a.e11.im;
    T[1][2].re = a.e12.re; T[1][2].im = a.e12.im;
    T[1][3].re = a.e13.re; T[1][3].im = a.e13.im;
    T[1][4].re = a.e14.re; T[1][4].im = a.e14.im;
    T[1][5].re = a.e15.re; T[1][5].im = a.e15.im;
    //2.row
    T[2][0].re = a.e20.re; T[2][0].im = a.e20.im;
    T[2][1].re = a.e21.re; T[2][1].im = a.e21.im;
    T[2][2].re = a.e22.re; T[2][2].im = a.e22.im;
    T[2][3].re = a.e23.re; T[2][3].im = a.e23.im;
    T[2][4].re = a.e24.re; T[2][4].im = a.e24.im;
    T[2][5].re = a.e25.re; T[2][5].im = a.e25.im;
    //3.row
    T[3][0].re = a.e30.re; T[3][0].im = a.e30.im;
    T[3][1].re = a.e31.re; T[3][1].im = a.e31.im;
    T[3][2].re = a.e32.re; T[3][2].im = a.e32.im;
    T[3][3].re = a.e33.re; T[3][3].im = a.e33.im;
    T[3][4].re = a.e34.re; T[3][4].im = a.e34.im;
    T[3][5].re = a.e35.re; T[3][5].im = a.e35.im;
    //4.row
    T[4][0].re = a.e40.re; T[4][0].im = a.e40.im;
    T[4][1].re = a.e41.re; T[4][1].im = a.e41.im;
    T[4][2].re = a.e42.re; T[4][2].im = a.e42.im;
    T[4][3].re = a.e43.re; T[4][3].im = a.e43.im;
    T[4][4].re = a.e44.re; T[4][4].im = a.e44.im;
    T[4][5].re = a.e45.re; T[4][5].im = a.e45.im;
    //5.row
    T[5][0].re = a.e50.re; T[5][0].im = a.e50.im;
    T[5][1].re = a.e51.re; T[5][1].im = a.e51.im;
    T[5][2].re = a.e52.re; T[5][2].im = a.e52.im;
    T[5][3].re = a.e53.re; T[5][3].im = a.e53.im;
    T[5][4].re = a.e54.re; T[5][4].im = a.e54.im;
    T[5][5].re = a.e55.re; T[5][5].im = a.e55.im;

    //compute the determinant of a 6x6 Matrix
    hmc_float s;
    hmc_complex sigma, det, q, z;
    hmc_complex p[6];
    det = hmc_complex_one;
    hmc_complex N1[6][6];

	//square input matrix
    for(unsigned int i=0; i<=n; ++i){
        for(unsigned int j=0; j<=n; ++j){
            N1[i][j] = hmc_complex_zero;
            for(unsigned int l=0; l<=n; ++l){
                N1[i][j] = complex_add(N1[i][j], complex_mult(T[i][l], T[l][j]));
            }}}
    for(unsigned int i=0; i<=n; ++i){
        for(unsigned int j=0; j<=n; ++j){
            T[i][j] = N1[i][j];}}

    for(unsigned int k=0; k<n; ++k) {
        s = 0.;
        for(unsigned int j=k+1; j<=n; ++j) {
            s += complex_abs_value(T[j][k]) * complex_abs_value(T[j][k]);
        }
        s = sqrt(1. + s / (complex_abs_value(T[k][k]) * complex_abs_value(T[k][k])));
        sigma = complex_mult(convert_float_to_complex(s), T[k][k]);
        
        /* determinant */
        det = complex_mult(det, sigma);
        q = complex_mult(sigma, complex_conj(T[k][k]));
        
        T[k][k] = complex_add(T[k][k], sigma);
        p[k] = complex_mult(sigma, complex_conj(T[k][k]));
        
        /* reflect all columns to the right */
        for(unsigned int j=k+1; j<=n; j++) {
            z = hmc_complex_zero;
            for(unsigned int i=k; i<=n; i++) {
                z = complex_add(z, complex_mult(complex_conj(T[i][k]), T[i][j]));
            }
            z = complex_divid(z, p[k]);
            for(unsigned int i=k; i<=n; i++) {
                T[i][j] = complex_sub(T[i][j], complex_mult(z, T[i][k]));
            }
        }
    }
    sigma = T[n][n];
    
    /* determinant */
    det = complex_mult(det, sigma);
	/* logarithm */
    hmc_float result = log(det.re);
    return result;
}


hmc_float log_det_matrix6x6_for_site(__global Matrixsu3StorageType const * const restrict field, st_idx const pos, hmc_float kappa_in, hmc_float csw)
{
	//calculate clover_explicit blocks and and call for every block log_det_squared and add results --> input = gaugefield
    Matrix6x6 UpperLeft = clover_eoprec_unified_local_upper_left_block(field, pos, kappa_in, csw);
	Matrix6x6 LowerRight = clover_eoprec_unified_local_lower_right_block(field, pos, kappa_in, csw);
	hmc_float tmp1 = log_det_squared(UpperLeft);
	hmc_float tmp2 = log_det_squared(LowerRight);
	hmc_float res_tmp = tmp1 + tmp2;
	//printf("%f", res_tmp);
	return res_tmp;
}

__kernel void clover_eo_log_det(__global Matrixsu3StorageType const * const restrict field, __global hmc_float * res, hmc_float kappa_in, hmc_float csw)
{
	hmc_float tmp = 0.;
    printf("kernel is called\n");
    PARALLEL_FOR(id_local, EOPREC_SPINORFIELDSIZE_LOCAL) {
        st_idx pos = get_even_st_idx_local(id_local);
        tmp -= log_det_matrix6x6_for_site(field, pos, kappa_in, csw);
        //printf("%f", log_det_matrix6x6_for_site(in, pos));
        //printf("\n");
    }
    //printf("%f", tmp);
	*res = tmp;
}


