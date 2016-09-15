/* This function performs the inversion a complex 6x6 matrix
   via Householder-Triangularization(see OpenQCD documentation "Implementation of the Dirac Operator" section 5
 */

Matrix6x6 inverse_6x6_via_Householder_triangularization(Matrix6x6 a)
{
    Matrix6x6 out;
    //map Matrix6x6 struct to C array6x6
    const int rows = 6; const int cols = 6;
    hmc_complex T[6][6];
    /*T = {{a.e00.re + I * a.e00.im, a.e01.re + I * a.e01.im, a.e02.re + I * a.e02.im, a.e03.re + I * a.e03.im, a.e04.re + I * a.e04.im, a.e05.re + I * a.e05.im},
        {a.e10.re + I * a.e10.im, a.e11.re + I * a.e11.im, a.e12.re + I * a.e12.im, a.e13.re + I * a.e13.im, a.e14.re + I * a.e14.im, a.e15.re + I * a.e05.im},
        {a.e20.re + I * a.e20.im, a.e21.re + I * a.e21.im, a.e22.re + I * a.e22.im, a.e23.re + I * a.e23.im, a.e24.re + I * a.e24.im, a.e25.re + I * a.e25.im},
        {a.e30.re + I * a.e30.im, a.e31.re + I * a.e31.im, a.e32.re + I * a.e32.im, a.e33.re + I * a.e33.im, a.e34.re + I * a.e34.im, a.e35.re + I * a.e35.im},
        {a.e40.re + I * a.e40.im, a.e41.re + I * a.e41.im, a.e42.re + I * a.e42.im, a.e43.re + I * a.e43.im, a.e44.re + I * a.e44.im, a.e45.re + I * a.e45.im},
        {a.e50.re + I * a.e50.im, a.e51.re + I * a.e51.im, a.e52.re + I * a.e52.im, a.e53.re + I * a.e53.im, a.e54.re + I * a.e54.im, a.e55.re + I * a.e55.im}};*/
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
    
    //build R_1,...,R_n-1,T via vectors u_k
    hmc_complex u[6];
    hmc_complex N1[6][6];
    hmc_complex R[6-1][6][6] = {0.};
    
    for(unsigned int k=0; k<(cols-1); ++k){
        //build u_k and norm^2 of u_k
        hmc_float norm_u_squared = 0;
        for(unsigned int l=0; l<rows; ++l){ //vector u according to Lüscher-OpenQCD doc
            if(l<k){u[l] = hmc_complex_zero;}
            else if(l==k){hmc_float r = 0.0;
                for(unsigned int j=k; j<rows; ++j){r = r + complex_abs_value(T[j][k]) * complex_abs_value(T[j][k]);}
                if(T[l][k].re ==0 && T[l][k].im == 0){u[l].re = - sqrt(r);} // 0/abs(0)!=1
                else{u[l] = complex_sub(T[k][k], complex_mult(complex_divid(T[k][k], convert_float_to_complex(complex_abs_value(T[k][k]))), convertfloattocomplex(sqrt(r))));}}
            else{u[l] = T[l][k];}
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
                N1[m][n] = complex_add(N1[m][n], complex_add(R[k][m][l], T[l][n]));
            }}}
    for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
            T[m][n] = N1[m][n];}}
    }//end for loop over k-->now we have R_n-1,...,R_1,T=upper triangular


    //inversion of T=S^(-1) according to OpenQCD doc equation (5.8)ff.
    //note: inverse of upper triangular matrix is upper triangular aswell
    hmc_complex S[6][6] = {0.};
    for(unsigned int j=0; j<cols; j=j+1){ //diagonal part: S_jj = T_jj^(-1)
        S[j][j] = complex_divid(hmc_complex_one, T[j][j]);
    }
    for(int k=rows-1; k>0; k=k-1){ //
        for(int i=k-1; i>=0; i=i-1){
            hmc_complex r = hmc_complex_zero;
            for(int j=i+1; j<=k; j=j+1){
                r = complex_add(r, complex_mult(T[i][j], S[j][k]));}
            S[i][k] = complexmult(hmc_complex_minusone, complex_mult(S[i][i], r));}
    }


    //R_prod = R_n-1,...,R_1
    hmc_complex R_prod[6][6] = {0.};
    for(unsigned int m=0; m<rows; ++m){ //initialise as R_1
        for(unsigned int n=0; n<cols; ++n){
            R_prod[m][n] = R[0][m][n];}}
    
    for(unsigned int k=1; k<(cols-1); k=k+1){ //multiply R_k's
        for(unsigned int m=0; m<rows; ++m){
            for(unsigned int n=0; n<cols; ++n){
                N1[m][n] = hmc_complex_zero;
                for(unsigned int l=0; l<rows; ++l){
                    N1[m][n] = complex_add(N1[m][n], complex_add(R[k][m][l], R_prod[l][n]));}}}
        for(unsigned int m=0; m<rows; ++m){
            for(unsigned int n=0; n<cols; ++n){
                R_prod[m][n] = N1[m][n];}}
    }

    // result = r = T^(-1) = S * R_prod
    hmc_complex r[6][6];
    for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
            r[m][n] = hmc_complex_zero;
            for(unsigned int l=0; l<rows; ++l){
                r[m][n] = complex_add(r[m][n], complex_mult(S[m][l], R_prod[l][n]));}}
    }


    //map C array6x6 to Matrix6x6 struct
    //0. row
    out.e00.re = r[0][0].re; out.e00.im = r[0][0].im;
    out.e01.re = r[0][1].re; out.e01.im = r[0][1].im;
    out.e02.re = r[0][2].re; out.e02.im = r[0][2].im;
    out.e03.re = r[0][3].re; out.e03.im = r[0][3].im;
    out.e04.re = r[0][4].re; out.e04.im = r[0][4].im;
    out.e05.re = r[0][5].re; out.e05.im = r[0][5].im;
    //1. row
    out.e10.re = r[1][0].re; out.e10.im = r[1][0].im;
    out.e11.re = r[1][1].re; out.e11.im = r[1][1].im;
    out.e12.re = r[1][2].re; out.e12.im = r[1][2].im;
    out.e13.re = r[1][3].re; out.e13.im = r[1][3].im;
    out.e14.re = r[1][4].re; out.e14.im = r[1][4].im;
    out.e15.re = r[1][5].re; out.e15.im = r[1][5].im;
    //2. row
    out.e20.re = r[2][0].re; out.e20.im = r[2][0].im;
    out.e21.re = r[2][1].re; out.e21.im = r[2][1].im;
    out.e22.re = r[2][2].re; out.e22.im = r[2][2].im;
    out.e23.re = r[2][3].re; out.e23.im = r[2][3].im;
    out.e24.re = r[2][4].re; out.e24.im = r[2][4].im;
    out.e25.re = r[2][5].re; out.e25.im = r[2][5].im;
    //3. row
    out.e30.re = r[3][0].re; out.e30.im = r[3][0].im;
    out.e31.re = r[3][1].re; out.e31.im = r[3][1].im;
    out.e32.re = r[3][2].re; out.e32.im = r[3][2].im;
    out.e33.re = r[3][3].re; out.e33.im = r[3][3].im;
    out.e34.re = r[3][4].re; out.e34.im = r[3][4].im;
    out.e35.re = r[3][5].re; out.e35.im = r[3][5].im;
    //4. row
    out.e40.re = r[4][0].re; out.e40.im = r[4][0].im;
    out.e41.re = r[4][1].re; out.e41.im = r[4][1].im;
    out.e42.re = r[4][2].re; out.e42.im = r[4][2].im;
    out.e43.re = r[4][3].re; out.e43.im = r[4][3].im;
    out.e44.re = r[4][4].re; out.e44.im = r[4][4].im;
    out.e45.re = r[4][5].re; out.e45.im = r[4][5].im;
    //5. row
    out.e50.re = r[5][0].re; out.e50.im = r[5][0].im;
    out.e51.re = r[5][1].re; out.e51.im = r[5][1].im;
    out.e52.re = r[5][2].re; out.e52.im = r[5][2].im;
    out.e53.re = r[5][3].re; out.e53.im = r[5][3].im;
    out.e54.re = r[5][4].re; out.e54.im = r[5][4].im;
    out.e55.re = r[5][5].re; out.e55.im = r[5][5].im;
    return out;
}

void clover_eo_inverse_explizit_upper_left_for_site(__global Matrix6x6StorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, hmc_float kappa_in, hmc_float csw, st_idx const pos)
{/*
    Matrix6x6 out_tmp, tmp;
    
    tmp = clover_eoprec_unified_local_upper_left_block(field, pos, csw);
    out_tmp = inverse_6x6_via_Householder_triangularization(tmp);
    put6x6(out, get_site_idx(pos), out_tmp);
*/}

void clover_eo_inverse_explizit_lower_right_for_site(__global Matrix6x6StorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, hmc_float kappa_in, hmc_float csw, st_idx const pos)
{/*
    Matrix6x6 out_tmp, tmp;
    
    tmp = clover_eoprec_unified_local_lower_right_block(field, pos, csw);
    out_tmp = inverse_6x6_via_Householder_triangularization(tmp);
    
    put6x6(out, get_site_idx(pos), out_tmp);
*/}

__kernel void clover_eo_inverse_explizit_upper_left(__global Matrix6x6StorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, hmc_float kappa_in, hmc_float csw)
{/*
    PARALLEL_FOR(id_local, SPINORFIELDSIZE_LOCAL) {
        st_idx pos = get_st_idx_from_site_idx(id_local);
        clover_eo_inverse_explizit_upper_left_for_site(out, field, kappa_in, csw, pos);
    }
*/}

__kernel void clover_eo_inverse_explizit_lower_right(__global Matrix6x6StorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, hmc_float kappa_in, hmc_float csw)
{/*
    PARALLEL_FOR(id_local, SPINORFIELDSIZE_LOCAL) {
        st_idx pos = get_st_idx_from_site_idx(id_local);
        clover_eo_inverse_explizit_lower_right_for_site(out, field, kappa_in, csw, pos);
    }
*/}


//clover Matrix (1+T)={{B_plus,0},{0,B_minus}} is blockdiagonal, therefore (1+T)^(-1)={{A_plus,0},{0,A_minus}} also
void clover_eo_inverse_for_site(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, hmc_float kappa_in, hmc_float csw, st_idx const pos)
{
    spinor out_tmp, phi;
    site_idx pos_eo = get_eo_site_idx_from_st_idx(pos);
    phi = getSpinor_eo(in, pos_eo);
    halfspinor tmp1, tmp2;
    tmp1.e0 = phi.e0;
    tmp1.e1 = phi.e1;
    tmp2.e0 = phi.e2;
    tmp2.e1 = phi.e3;
    
    //get 6x6 blocks
    Matrix6x6 B_plus = clover_eoprec_unified_local_upper_left_block(field, pos, csw);
    Matrix6x6 B_minus = clover_eoprec_unified_local_lower_right_block(field, pos, csw);
    
    //inversion of 6x6 blocks
    Matrix6x6 A_plus = inverse_6x6_via_Householder_triangularization(B_plus);
    Matrix6x6 A_minus = inverse_6x6_via_Householder_triangularization(B_minus);
    
    //Inverse * spinor
    tmp1 = matrix6x6_times_halfspinor(A_plus, tmp1);
    tmp2 = matrix6x6_times_halfspinor(A_minus, tmp1);
    
    out_tmp.e0 = tmp1.e0;
    out_tmp.e1 = tmp1.e1;
    out_tmp.e2 = tmp2.e0;
    out_tmp.e3 = tmp2.e1;
    
    putSpinor_eo(out, get_eo_site_idx_from_st_idx(pos), out_tmp);
}

__kernel void clover_eo_inverse(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in, hmc_float csw)
{/*
    PARALLEL_FOR(id_local, EOPREC_SPINORFIELDSIZE_LOCAL) {
        st_idx pos = (evenodd == EVEN) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);
        clover_eo_inverse_for_site(in, out, field, kappa_in, csw, pos);
    }
*/}
