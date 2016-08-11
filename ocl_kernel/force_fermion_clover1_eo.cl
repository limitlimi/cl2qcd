/*
 calculate the even-odd-preconditioned force of the clover term with eoprec spinorfields (cf. paper of Jansen, Liu: Implementation of Symanzik's improvement program for simulations of dynamical wilson fermions in lattice QCD, URL:http://arxiv.org/abs/hep-lat/9603008)
*/


__kernel void fermion_force_clover1_eo_0(__global const Matrixsu3StorageType * const restrict field, __global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Z, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);

        Matrixsu3 U;
        Matrix3x3 v1, v2, tmp;
        su3vec psia, psib, phia, phib;
        spinor y, x, tmp;
        ae out_tmp;
        //this is used to save the BC-conditions...
        hmc_complex bc_tmp;
        int dir;
        int global_link_pos;
        int n = pos.space;
        int t = pos.time;
        int nn_eo;
        
        if(evenodd == ODD) {
            x = getSpinor_eo(X, get_eo_site_idx_from_st_idx(pos));
        } else if(evenodd == EVEN) {
            pos_odd = get_odd_st_idx(id_mem);
            //x_even = -(1+T_ee)^{-1} M_eo x_odd
            tmp = getSpinor_eo(X, get_eo_site_idx_from_st_idx(pos_odd));
            // x = tmp;
        }
        
        y = getSpinor_eo(Y, get_eo_site_idx_from_st_idx(pos));
        global_link_pos = get_link_pos(dir, n, t);
        ///////////////////////////////////
        // Calculate gamma_5 y
        ///////////////////////////////////
        y = gamma5_local(y);
        
        //go through the different directions
        ///////////////////////////////////
        // mu = 0
        ///////////////////////////////////
        dir1 = 0;
        //the 2 here comes from Tr(lambda_ij) = 2delta_ij
        hmc_float c_0_hat = 1/(1 + 64 * kappa_in * kappa_in);
        bc_tmp.re = 2./8. * kappa_in * c_0_hat * csw * TEMPORAL_RE;
        bc_tmp.im = 2./8. * kappa_in * c_0_hat * csw * TEMPORAL_IM;

        //if chemical potential is activated, U has to be multiplied by appropiate factor

        
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        
        v1 = identity_matrix3x3();
        
        //diagrams
        v2 = diagram1a_up(...);
        v1 = add_matrix3x3(matrix_su3to3x3(v2), v1);
        
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        
        
        v1 = multiply_matrix3x3_by_complex(v1, bc_tmp);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}

//calculates the Dirac-Trace of the matrix resulting from multiplying X*Y^dagger =  x1*y1^dagger + x2*y2^dagger + x3*y3^dagger + x4*y4^dagger , where x1,x2,x3,x4,y1,y2,y3,y4 are SU(3)-vectors. The result is a 3x3-matrix
Matrix3x3 tr_dirac_x_times_y_dagger(su3vec x1, su3vec x2, su3vec x3, su3vec x4, su3vec y1, su3vec y2, su3vec y3, su3vec y4)
{
    Matrix3x3 tmp;
    tmp.e00.re = (x1).e0.re * (y1).e0.re + (x1).e0.im * (y1).e0.im + (x2).e0.re * (y2).e0.re + (x2).e0.im * (y2).e0.im + (x3).e0.re * (y3).e0.re + (x3).e0.im * (y3).e0.im + (x4).e0.re * (y4).e0.re + (x4).e0.im * (y4).e0.im;
    tmp.e00.im = (x1).e0.im * (y1).e0.re - (x1).e0.re * (y1).e0.im + (x2).e0.im * (y2).e0.re - (x2).e0.re * (y2).e0.im + (x3).e0.im * (y3).e0.re - (x3).e0.re * (y3).e0.im + (x4).e0.im * (y4).e0.re - (x4).e0.re * (y4).e0.im;
    tmp.e01.re = (x1).e0.re * (y1).e1.re + (x1).e0.im * (y1).e1.im + (x2).e0.re * (y2).e1.re + (x2).e0.im * (y2).e1.im + (x3).e0.re * (y3).e1.re + (x3).e0.im * (y3).e1.im + (x4).e0.re * (y4).e1.re + (x4).e0.im * (y4).e1.im;
    tmp.e01.im = (x1).e0.im * (y1).e1.re - (x1).e0.re * (y1).e1.im + (x2).e0.im * (y2).e1.re - (x2).e0.re * (y2).e1.im + (x3).e0.im * (y3).e1.re - (x3).e0.re * (y3).e1.im + (x4).e0.im * (y4).e1.re - (x4).e0.re * (y4).e1.im;
    tmp.e02.re = (x1).e0.re * (y1).e2.re + (x1).e0.im * (y1).e2.im + (x2).e0.re * (y2).e2.re + (x2).e0.im * (y2).e2.im + (x3).e0.re * (y3).e2.re + (x3).e0.im * (y3).e2.im + (x4).e0.re * (y4).e2.re + (x4).e0.im * (y4).e2.im;
    tmp.e02.im = (x1).e0.im * (y1).e2.re - (x1).e0.re * (y1).e2.im + (x2).e0.im * (y2).e2.re - (x2).e0.re * (y2).e2.im + (x3).e0.im * (y3).e2.re - (x3).e0.re * (y3).e2.im + (x4).e0.im * (y4).e2.re - (x4).e0.re * (y4).e2.im;
    tmp.e10.re = (x1).e1.re * (y1).e0.re + (x1).e1.im * (y1).e0.im + (x2).e1.re * (y2).e0.re + (x2).e1.im * (y2).e0.im + (x3).e1.re * (y3).e0.re + (x3).e1.im * (y3).e0.im + (x4).e1.re * (y4).e0.re + (x4).e1.im * (y4).e0.im;
    tmp.e10.im = (x1).e1.im * (y1).e0.re - (x1).e1.re * (y1).e0.im + (x2).e1.im * (y2).e0.re - (x2).e1.re * (y2).e0.im + (x3).e1.im * (y3).e0.re - (x3).e1.re * (y3).e0.im + (x4).e1.im * (y4).e0.re - (x4).e1.re * (y4).e0.im;
    tmp.e11.re = (x1).e1.re * (y1).e1.re + (x1).e1.im * (y1).e1.im + (x2).e1.re * (y2).e1.re + (x2).e1.im * (y2).e1.im + (x3).e1.re * (y3).e1.re + (x3).e1.im * (y3).e1.im + (x4).e1.re * (y4).e1.re + (x4).e1.im * (y4).e1.im;
    tmp.e11.im = (x1).e1.im * (y1).e1.re - (x1).e1.re * (y1).e1.im + (x2).e1.im * (y2).e1.re - (x2).e1.re * (y2).e1.im + (x3).e1.im * (y3).e1.re - (x3).e1.re * (y3).e1.im + (x4).e1.im * (y4).e1.re - (x4).e1.re * (y4).e1.im;
    tmp.e12.re = (x1).e1.re * (y1).e2.re + (x1).e1.im * (y1).e2.im + (x2).e1.re * (y2).e2.re + (x2).e1.im * (y2).e2.im + (x3).e1.re * (y3).e2.re + (x3).e1.im * (y3).e2.im + (x4).e1.re * (y4).e2.re + (x4).e1.im * (y4).e2.im;
    tmp.e12.im = (x1).e1.im * (y1).e2.re - (x1).e1.re * (y1).e2.im + (x2).e1.im * (y2).e2.re - (x2).e1.re * (y2).e2.im + (x3).e1.im * (y3).e2.re - (x3).e1.re * (y3).e2.im + (x4).e1.im * (y4).e2.re - (x4).e1.re * (y4).e2.im;
    tmp.e20.re = (x1).e2.re * (y1).e0.re + (x1).e2.im * (y1).e0.im + (x2).e2.re * (y2).e0.re + (x2).e2.im * (y2).e0.im + (x3).e2.re * (y3).e0.re + (x3).e2.im * (y3).e0.im + (x4).e2.re * (y4).e0.re + (x4).e2.im * (y4).e0.im;
    tmp.e20.im = (x1).e2.im * (y1).e0.re - (x1).e2.re * (y1).e0.im + (x2).e2.im * (y2).e0.re - (x2).e2.re * (y2).e0.im + (x3).e2.im * (y3).e0.re - (x3).e2.re * (y3).e0.im + (x4).e2.im * (y4).e0.re - (x4).e2.re * (y4).e0.im;
    tmp.e21.re = (x1).e2.re * (y1).e1.re + (x1).e2.im * (y1).e1.im + (x2).e2.re * (y2).e1.re + (x2).e2.im * (y2).e1.im + (x3).e2.re * (y3).e1.re + (x3).e2.im * (y3).e1.im + (x4).e2.re * (y4).e1.re + (x4).e2.im * (y4).e1.im;
    tmp.e21.im = (x1).e2.im * (y1).e1.re - (x1).e2.re * (y1).e1.im + (x2).e2.im * (y2).e1.re - (x2).e2.re * (y2).e1.im + (x3).e2.im * (y3).e1.re - (x3).e2.re * (y3).e1.im + (x4).e2.im * (y4).e1.re - (x4).e2.re * (y4).e1.im;
    tmp.e22.re = (x1).e2.re * (y1).e2.re + (x1).e2.im * (y1).e2.im + (x2).e2.re * (y2).e2.re + (x2).e2.im * (y2).e2.im + (x3).e2.re * (y3).e2.re + (x3).e2.im * (y3).e2.im + (x4).e2.re * (y4).e2.re + (x4).e2.im * (y4).e2.im;
    tmp.e22.im = (x1).e2.im * (y1).e2.re - (x1).e2.re * (y1).e2.im + (x2).e2.im * (y2).e2.re - (x2).e2.re * (y2).e2.im + (x3).e2.im * (y3).e2.re - (x3).e2.re * (y3).e2.im + (x4).e2.im * (y4).e2.re - (x4).e2.re * (y4).e2.im;
    return tmp;
}

//input eo spinor X and Y where one is multiplied by gamma_5*sigma_{dir1,dir2}, Faktor i??
Matrixsu3 diagram1a_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Z, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
    U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
    hmc_complex cpi_tmp = {COSCPI, SINCPI};
    U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3(out, U);
    
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // square(x+nu)
    idx_neigh = get_neighbor_from_st_idx(pos, dir2);
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh)
    out = multiply_matrixsu3(out, square(Z, Y, X, idx_neigh_eo));
    
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    return out;
}

Matrixsu3 diagram1a_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Z, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // square(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh)
    out = multiply_matrixsu3(out, square(Z, Y, X, idx_neigh_eo));
    
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3(out, U);
    
    out = multiply_matrixsu3_by_real (out, -1.)
    return out;
}


Matrixsu3 diagram1b_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Z, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3(out, U);
    
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // square(x)
    idx_arg_eo = get_eo_site_idx_from_st_idx(idx_arg)
    out = multiply_matrixsu3(out, square(Z, Y, X, idx_arg_eo));
    
    return out;
}

Matrixsu3 diagram1b_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Z, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3(out, U);
    
    // square(x)
    idx_arg_eo = get_eo_site_idx_from_st_idx(idx_arg)
    out = multiply_matrixsu3(out, square(Z, Y, X, idx_arg_eo));
    
    out = multiply_matrixsu3_by_real (out, -1.)
    return out;
}


Matrixsu3 diagram1c_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Z, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // square(x+mu)
    idx_neigh = get_neighbor_from_st_idx(pos, dir1);
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh)
    out = multiply_matrixsu3(out, square(Z, Y, X, idx_neigh_eo));
    
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3(out, U);
    
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    return out;
}

Matrixsu3 diagram1c_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Z, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // square(x+mu)
    idx_neigh = get_neighbor_from_st_idx(pos, dir1);
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh)
    out = multiply_matrixsu3(out, square(Z, Y, X, idx_neigh_eo));
    
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3(out, U);
    
    out = multiply_matrixsu3_by_real (out, -1.)
    return out;
}



Matrixsu3 diagram1d_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Z, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3(out, U);

    // square(x+mu+nu)
    idx_neigh1 = get_neighbor_from_st_idx(idx_arg, dir2); // x+nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x+nu)+mu
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh)
    out = multiply_matrixsu3(out, square(Z, Y, X, idx_neigh_eo));
    
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    return out;
}

Matrixsu3 diagram1d_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Z, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);
    
    // square(x-nu+mu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh)
    out = multiply_matrixsu3(out, square(Z, Y, X, idx_neigh_eo));
    
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3_dagger(out, U);

    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
    }
    out = multiply_matrixsu3(out, U);
    
    out = multiply_matrixsu3_by_real (out, -1.)
    return out;
}



//matrixsu3 which is given by Dirac-Trace(gamma_5*sigma_{mu,nu} Y*X^dagger + gamma_5*sigma_{mu,nu} X*Y^dagger)
//cf. equation (22) in Jansen Liu Paper
Matrixsu3 square(__global const spinorStorageType * const restrict Z, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X, const site_idx pos_square)
{
    Matrixsu3 out;
    //if statement for X_even, Y_even
    //sigma_mu,nu*gamma_5*X
    return out;
}
                      
            



__kernel void fermion_force_clover1_eo_1(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);

    }
}

__kernel void fermion_force_clover1_eo_2(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);

    }
}

__kernel void fermion_force_clover1_eo_3(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);

    }
}