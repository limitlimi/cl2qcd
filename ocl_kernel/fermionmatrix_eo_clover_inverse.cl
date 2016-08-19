void clover_eo_inverse_for_site(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, hmc_float kappa_in, hmc_float csw, st_idx const pos)
{
    halfspinor tmp1, tmp2;
    tmp1.e1 = in.e1;
    tmp1.e2 = in.e2;
    tmp2.e1 = in.e3;
    tmp2.e2 = in.e4;
    spinor out_tmp;
    
    Matrix6x6 B_plus = clover_eoprec_unified_local_upper_left_block(field, pos, csw);
    Matrix6x6 B_minus = clover_eoprec_unified_local_lower_right_block(field, pos, csw);
    Matrix6x6 A_plus = inverse_6x6_via_Householder_triangularization(B_plus);
    Matrix6x6 A_minus = inverse_6x6_via_Householder_triangularization(B_minus);
    tmp1 = matrix6x6_times_halfspinor(A_plus, tmp1);
    tmp2 = matrix6x6_times_halfspinor(A_minus, tmp1);
    
    out_tmp.e1 = tmp1.e1;
    out_tmp.e2 = tmp1.e2;
    out_tmp.e3 = tmp2.e1;
    out_tmp.e4 = tmp2.e2;
    
    putSpinor_eo(out, get_eo_site_idx_from_st_idx(pos), out_tmp);
}

__kernel void clover_eo_inverse(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in, hmc_float csw)
{
    PARALLEL_FOR(id_local, EOPREC_SPINORFIELDSIZE_LOCAL) {
        st_idx pos = (evenodd == EVEN) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);
        clover_eo_inverse_for_site(in, out, field, kappa_in, csw, pos);
    }
}

/* This function performs the inversion a complex 6x6 matrix
 via Householder-Triangularization(see OpenQCD documentation "Implementation of the Dirac Operator" section 5
 */
Matrix6x6 inverse_6x6_via_Householder_triangularization(Matrix6x6 a)
{
    Matrix6x6 out;
    //map Matrix6x6 struct to C array6x6
    const int rows = 6; const int cols = 6;
    double complex T[rows][cols];
    T = {{a.e00.re + I * a.e00.im, a.e01.re + I * a.e01.im, a.e02.re + I * a.e02.im, a.e03.re + I * a.e03.im, a.e04.re + I * a.e04.im, a.e05.re + I * a.e05.im},
        {a.e10.re + I * a.e10.im, a.e11.re + I * a.e11.im, a.e12.re + I * a.e12.im, a.e13.re + I * a.e13.im, a.e14.re + I * a.e14.im, a.e15.re + I * a.e05.im},
        {a.e20.re + I * a.e20.im, a.e21.re + I * a.e21.im, a.e22.re + I * a.e22.im, a.e23.re + I * a.e23.im, a.e24.re + I * a.e24.im, a.e25.re + I * a.e25.im},
        {a.e30.re + I * a.e30.im, a.e31.re + I * a.e31.im, a.e32.re + I * a.e32.im, a.e33.re + I * a.e33.im, a.e34.re + I * a.e34.im, a.e35.re + I * a.e35.im},
        {a.e40.re + I * a.e40.im, a.e41.re + I * a.e41.im, a.e42.re + I * a.e42.im, a.e43.re + I * a.e43.im, a.e44.re + I * a.e44.im, a.e45.re + I * a.e45.im},
        {a.e50.re + I * a.e50.im, a.e51.re + I * a.e51.im, a.e52.re + I * a.e52.im, a.e53.re + I * a.e53.im, a.e54.re + I * a.e54.im, a.e55.re + I * a.e55.im}};
    
    //build R_1,...,R_n-1,T via vectors u_k
    double complex u[rows];
    double complex N1[rows][cols];
    double complex R[cols-1][rows][cols] = {0.};
    
    for(unsigned int k=0; k<(cols-1); ++k){
        //build u_k and norm^2 of u_k
        double norm_u_squared = 0;
        for(unsigned int l=0; l<rows; ++l){ //vector u according to Lüscher-OpenQCD doc
            if(l<k){u[l] = 0.;}
            else if(l==k){double r = 0.0;
                for(unsigned int j=k; j<rows; ++j){r = r + cabs(T[j][k]) * cabs(T[j][k]);}
                if(T[l][k] ==0){u[l] = - sqrt(r);} // 0/abs(0)!=1
                else{u[l] = T[k][k] - T[k][k]/cabs(T[k][k])*sqrt(r);}}
            else{u[l] = T[l][k];}
            norm_u_squared = norm_u_squared + cabs(u[l])*cabs(u[l]);}
    // build up R_k by u_k und perform R_k * R_k-1*...*R_1*T
    for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
            R[k][m][n] = -2.0/norm_u_squared * u[m] * conj(u[n]);
            if(m==n){R[k][m][n] = 1 + R[k][m][n];}}}
    for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
            N1[m][n] = 0.;
            for(unsigned int l=0; l<rows; ++l){
                N1[m][n] = N1[m][n] + R[k][m][l] * T[l][n];
            }}}
    for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
            T[m][n] = N1[m][n];}}
    }//end for loop over k-->now we have R_n-1,...,R_1,T=upper triangular


    //inversion of T=S^(-1) according to OpenQCD doc equation (5.8)ff.
    //note: inverse of upper triangular matrix is upper triangular aswell
    double complex S[rows][cols] = {0.};
    for(unsigned int j=0; j<cols; j=j+1){ //diagonal part: S_jj = T_jj^(-1)
        S[j][j] = 1.0 / T[j][j];
    }
    for(int k=rows-1; k>0; k=k-1){ //
        for(int i=k-1; i>=0; i=i-1){
            double complex r = 0.;
            for(int j=i+1; j<=k; j=j+1){
                r = r + T[i][j] * S[j][k];
            }
            S[i][k] = - S[i][i] * r;
        }
    }


    //R_prod = R_n-1,...,R_1
    double complex R_prod[rows][cols] = {0.};
    for(unsigned int m=0; m<rows; ++m){ //initialise as R_1
    for(unsigned int n=0; n<cols; ++n){
        R_prod[m][n] = R[0][m][n];}}
    for(unsigned int k=1; k<(cols-1); k=k+1){ //multiply R_k's
        for(unsigned int m=0; m<rows; ++m){
            for(unsigned int n=0; n<cols; ++n){
                N1[m][n] = 0.;
                for(unsigned int l=0; l<rows; ++l){
                    N1[m][n] = N1[m][n] + R[k][m][l] * R_prod[l][n];
                }}}
        for(unsigned int m=0; m<rows; ++m){
            for(unsigned int n=0; n<cols; ++n){
                R_prod[m][n] = N1[m][n];}}
    }

    // result = r = T^(-1) = S * R_prod
    double complex r[rows][cols];
    for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
            r[m][n] = 0.;
            for(unsigned int l=0; l<rows; ++l){
                r[m][n] = r[m][n] + S[m][l] * R_prod[l][n];}}
    }


    //map C array6x6 to Matrix6x6 struct
    //0. row
    out.e00.re = creal(r[0][0]); out.e00.im = cimag(r[0][0]);
    out.e01.re = creal(r[0][1]); out.e01.im = cimag(r[0][1]);
    out.e02.re = creal(r[0][2]); out.e02.im = cimag(r[0][2]);
    out.e03.re = creal(r[0][3]); out.e03.im = cimag(r[0][3]);
    out.e04.re = creal(r[0][4]); out.e04.im = cimag(r[0][4]);
    out.e05.re = creal(r[0][5]); out.e05.im = cimag(r[0][5]);
    //1. row
    out.e10.re = creal(r[1][0]); out.e10.im = cimag(r[1][0]);
    out.e11.re = creal(r[1][1]); out.e11.im = cimag(r[1][1]);
    out.e12.re = creal(r[1][2]); out.e12.im = cimag(r[1][2]);
    out.e13.re = creal(r[1][3]); out.e13.im = cimag(r[1][3]);
    out.e14.re = creal(r[1][4]); out.e14.im = cimag(r[1][4]);
    out.e15.re = creal(r[1][5]); out.e15.im = cimag(r[1][5]);
    //2. row
    out.e20.re = creal(r[2][0]); out.e20.im = cimag(r[2][0]);
    out.e21.re = creal(r[2][1]); out.e21.im = cimag(r[2][1]);
    out.e22.re = creal(r[2][2]); out.e22.im = cimag(r[2][2]);
    out.e23.re = creal(r[2][3]); out.e23.im = cimag(r[2][3]);
    out.e24.re = creal(r[2][4]); out.e24.im = cimag(r[2][4]);
    out.e25.re = creal(r[2][5]); out.e25.im = cimag(r[2][5]);
    //3. row
    out.e30.re = creal(r[3][0]); out.e30.im = cimag(r[3][0]);
    out.e31.re = creal(r[3][1]); out.e31.im = cimag(r[3][1]);
    out.e32.re = creal(r[3][2]); out.e32.im = cimag(r[3][2]);
    out.e33.re = creal(r[3][3]); out.e33.im = cimag(r[3][3]);
    out.e34.re = creal(r[3][4]); out.e34.im = cimag(r[3][4]);
    out.e35.re = creal(r[3][5]); out.e35.im = cimag(r[3][5]);
    //4. row
    out.e40.re = creal(r[4][0]); out.e40.im = cimag(r[4][0]);
    out.e41.re = creal(r[4][1]); out.e41.im = cimag(r[4][1]);
    out.e42.re = creal(r[4][2]); out.e42.im = cimag(r[4][2]);
    out.e43.re = creal(r[4][3]); out.e43.im = cimag(r[4][3]);
    out.e44.re = creal(r[4][4]); out.e44.im = cimag(r[4][4]);
    out.e45.re = creal(r[4][5]); out.e45.im = cimag(r[4][5]);
    //5. row
    out.e50.re = creal(r[5][0]); out.e50.im = cimag(r[5][0]);
    out.e51.re = creal(r[5][1]); out.e51.im = cimag(r[5][1]);
    out.e52.re = creal(r[5][2]); out.e52.im = cimag(r[5][2]);
    out.e53.re = creal(r[5][3]); out.e53.im = cimag(r[5][3]);
    out.e54.re = creal(r[5][4]); out.e54.im = cimag(r[5][4]);
    out.e55.re = creal(r[5][5]); out.e55.im = cimag(r[5][5]);
    return out;
}
