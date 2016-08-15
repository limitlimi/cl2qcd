/* This function performs locally the inversion a complex 6x6 matrix
   via Householder-Triangu 
 */

Matrix6x6 inverse_6x6_via_Householder_triangularization(__global Matrix6x6StorageType  const * const restrict field, const st_idx idx_arg)
{
    Matrix6x6 out;
    a = getMatrix6x6(field, idx_arg);
    //map Matrix6x6 struct to C array6x6
    const int rows = 6; const int cols = 6;
    double complex T[rows][cols] ={ 0.};
    
    T[0][0]  = a.e00.re + I * a.e00.im
    
    
    
    out.e00.r = creal(r[0][0]); out.e00.cimag(r[0][0]);
    return out;
}

Inversion()
{
    // initialize matrix which should be inverted
    const int rows = 3;
    const int cols = 3;
    double complex T[rows][cols] ={ 0.};
    T[0][0]=2.0*I; T[0][1]=1; T[0][2]=1.0*I;
    T[1][0]=1.0*I; T[1][1]=1.0*I; T[1][2]=1.0;
    T[2][0]=1.0; T[2][1]=1.0*I; T[2][2]=1.0*I;
    /*for(unsigned int i=0; i<rows; ++i){ \\print T
     for(unsigned int j=0; j<cols; ++j){
     printf("%.2f + %.2fi\t", creal(T[i][j]), cimag(T[i][j]));}
     printf("\n");}*/
    
    
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
                u[l] = T[k][k] - T[k][k]/cabs(T[k][k])*sqrt(r);}
            else{u[l] = T[l][k];}
            norm_u_squared = norm_u_squared + cabs(u[l])*cabs(u[l]);}
        /*printf("%f + %fi\n", creal(u[l]), cimag(u[l]));}
         printf("%f\n", norm_u_squared);*/
        /*double r = 0.0; //vector u according to Harrach script
         for(unsigned int j=k; j<rows; ++j){r = r + cabs(T[j][k]) * cabs(T[j][k]);}
         r = sqrt(r);
         for(unsigned int l=0; l<rows; ++l){
         if(l<k){u[l] = 0.;}
         else{u[l] = T[l][k]/r * cabs(T[k][k])/T[k][k];}
         if(l==k){u[l] = u[l] + 1;}
         norm_u_squared = norm_u_squared + cabs(u[l])*cabs(u[l]);
         printf("%f + %fi\n", creal(u[l]), cimag(u[l]));}
         printf("%f\n", norm_u_squared);*/
        
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
        /*for(unsigned int i=0; i<rows; ++i){ //print R_k
         for(unsigned int j=0; j<cols; ++j){
         printf("%.2f + %.2fi\t", creal(R[k][i][j]), cimag(R[k][i][j]));}
         printf("\n");}
         for(unsigned int i=0; i<rows; ++i){ //print R_k-1*...*R_1*A
         for(unsigned int j=0; j<cols; ++j){
         printf("%.2f + %.2fi\t", creal(T[i][j]), cimag(T[i][j]));}
         printf("\n");}*/
    }//end for loop over k-->now we have R_n-1,...,R_1,T=upper triangular
    
    
    //inversion of T=S^(-1) according to OpenQCD doc equation (5.8)ff.
    //note: inverse of upper triangular matrix is upper triangular aswell
    double complex S[rows][cols] = {0.};
    /*printf("T=\n");
     for(unsigned int i=0; i<rows; ++i){ //print R_k-1*...*R_1*A=upper triangular
     for(unsigned int j=0; j<cols; ++j){
     printf("%.2f + %.2fi\t", creal(T[i][j]), cimag(T[i][j]));}
     printf("\n");}*/
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
    /*printf("S=\n");
     for(unsigned int i=0; i<rows; ++i){ //print T^(-1)
     for(unsigned int j=0; j<cols; ++j){
     printf("%.2f + %.2fi\t", creal(S[i][j]), cimag(S[i][j]));}
     printf("\n");}*/
    
    
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
    
    
    printf("R_prod=\n");
    for(unsigned int i=0; i<rows; ++i){
        for(unsigned int j=0; j<cols; ++j){
            printf("%.2f + %.2fi\t", creal(R_prod[i][j]), cimag(R_prod[i][j]));}
        printf("\n");}
    
    // result = T^(-1) = S * R_prod
    double complex result[rows][cols];
    for(unsigned int m=0; m<rows; ++m){
        for(unsigned int n=0; n<cols; ++n){
            result[m][n] = 0.;
            for(unsigned int l=0; l<rows; ++l){
                result[m][n] = result[m][n] + S[m][l] * R_prod[l][n];
            }}}
    printf("Result=\n");
    for(unsigned int i=0; i<rows; ++i){
        for(unsigned int j=0; j<cols; ++j){
            printf("%.2f + %.2fi\t", creal(result[i][j]), cimag(result[i][j]));}
        printf("\n");}
}