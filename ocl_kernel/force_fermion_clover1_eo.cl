/*
 calculate the even-odd-preconditioned force of the clover term with eoprec spinorfields (cf. paper of Jansen, Liu: Implementation of Symanzik's improvement program for simulations of dynamical wilson fermions in lattice QCD, URL:http://arxiv.org/abs/hep-lat/9603008)
*/

//spinorfield input X1/Y1 ^= X_even/Y_even = -(1+T_ee)^(-1) M_eo X_o/Y_o
__kernel void fermion_force_clover1_eo_0(__global const Matrixsu3StorageType * const restrict field, __global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);

        Matrix3x3 v1, v2;
        ae out_tmp;
        //this is used to save the BC-conditions...
        hmc_complex bc_tmp;
        int dir1, dir2;
        int global_link_pos;
        int n = pos.space;
        int t = pos.time;

        //go through the different directions
        ///////////////////////////////////
        // mu = 0
        ///////////////////////////////////
        dir1 = 0;
        global_link_pos = get_link_pos(dir1, n, t);
        
        //the 2 here comes from Tr(lambda_ij) = 2delta_ij
        hmc_float c_0_hat = 1/(1 + 64 * kappa_in * kappa_in);
        bc_tmp.re = 2./8. * kappa_in * c_0_hat * csw * TEMPORAL_RE;
        bc_tmp.im = 2./8. * kappa_in * c_0_hat * csw * TEMPORAL_IM;
        
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_complex(v1, bc_tmp);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}

Matrix3x3 add_up_diagrams(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1 , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrix3x3 out, tmp;
    out = zero_matrix3x3();
    
    tmp = diagram1a_up(field, X, X1, Y, Y1, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(v2), v1);
    tmp = diagram1a_down(field, X, X1, Y, Y1, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(v2), v1);
    tmp = diagram1b_up(field, X, X1, Y, Y1, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(v2), v1);
    tmp = diagram1b_down(field, X, X1, Y, Y1, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(v2), v1);
    tmp = diagram1c_up(field, X, X1, Y, Y1, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(v2), v1);
    tmp = diagram1c_down(field, X, X1, Y, Y1, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(v2), v1);
    tmp = diagram1d_up(field, X, X1, Y, Y1, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(v2), v1);
    tmp = diagram1d_down(field, X, X1, Y, Y1, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(v2), v1);
    
    return out;
}

//input eo spinor X and Y where one is multiplied by gamma_5*sigma_{dir1,dir2}, Faktor i??
Matrixsu3 diagram1a_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1 , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
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
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, X1, Y, Y1, idx_neigh_eo, dir1, dir2, evenodd));
    
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

Matrixsu3 diagram1a_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1 , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
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
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, X1, Y, Y1, idx_neigh_eo, dir1, dir2, evenodd));
    
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


Matrixsu3 diagram1b_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1 , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
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
    idx_arg_eo = get_eo_site_idx_from_st_idx(idx_arg);
    out = multiply_matrixsu3(out, square(X, X1, Y, Y1, idx_arg_eo, dir1, dir2, evenodd));
    
    return out;
}

Matrixsu3 diagram1b_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1 , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
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
    idx_arg_eo = get_eo_site_idx_from_st_idx(idx_arg);
    out = multiply_matrixsu3(out, square(X, X1, Y, Y1, idx_arg_eo, dir1, dir2, evenodd));
    
    out = multiply_matrixsu3_by_real (out, -1.)
    return out;
}


Matrixsu3 diagram1c_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1 , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // square(x+mu)
    idx_neigh = get_neighbor_from_st_idx(pos, dir1);
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, X1, Y, Y1, idx_neigh_eo, dir1, dir2, evenodd));
    
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

Matrixsu3 diagram1c_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1 , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // square(x+mu)
    idx_neigh = get_neighbor_from_st_idx(pos, dir1);
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, X1, Y, Y1, idx_neigh_eo, dir1, dir2, evenodd));
    
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



Matrixsu3 diagram1d_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1 , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
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
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, X1, Y, Y1, idx_neigh_eo, dir1, dir2, evenodd));
    
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

Matrixsu3 diagram1d_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1 , const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
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
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, X1, Y, Y1, idx_neigh_eo, dir1, dir2, evenodd));
    
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    if(dir1 == 0){
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
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



//matrixsu3 which is given by i * Trace_Dirac[gamma_5*sigma_{mu,nu} Y*X^dagger + gamma_5*sigma_{mu,nu} X*Y^dagger]
//cf. equation (22) in Jansen Liu Paper
Matrixsu3 square(__global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict X1, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict Y1, const site_idx pos_square, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 out;
    Matrix3x3 tmp1, tmp2;
    spinor x, y;
    if(evenodd == ODD) {
        x = getSpinor_eo(X, get_eo_site_idx_from_st_idx(pos_square));
        y = getSpinor_eo(Y, get_eo_site_idx_from_st_idx(pos_square));
    }
    else if(evenodd == EVEN) {
        pos_odd = get_odd_st_idx(pos_square);//??????
        //x_even/y_even = -(1+T_ee)^{-1} M_eo x_odd/y_odd
        x = getSpinor_eo(X1, get_eo_site_idx_from_st_idx(pos_odd));
        y = getSpinor_eo(Y1, get_eo_site_idx_from_st_idx(pos_odd));
    }
    //calculate gamma_5 * simga_{mu,nu} * x
    x = gamma5_local(x);
    x = sigma_mu_nu_times_spinor(x, dir1, dir2);
    
    tmp1 = tr_dirac_x_times_y_dagger(y.e0, y.e1, y.e2, y.e3, x.e0, x.e1, x.e2, x.e3);
    tmp2 = tr_dirac_x_times_y_dagger(x.e0, x.e1, x.e2, x.e3, y.e0, y.e1, y.e2, y.e3);
    
    tmp1 = add_matrix3x3(tmp1, tmp2);
    
    out = matrix_3x3tosu3(tmp1);
    out = multiply_matrixsu3_by_complex(out, {0., 1.});
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

        Matrix3x3 v1, v2;
        ae out_tmp;
        //this is used to save the BC-conditions...
        hmc_complex bc_tmp;
        int dir1, dir2;
        int global_link_pos;
        int n = pos.space;
        int t = pos.time;
        
        //go through the different directions
        ///////////////////////////////////
        // mu = 1
        ///////////////////////////////////
        dir1 = 1;
        global_link_pos = get_link_pos(dir, n, t);
        
        //the 2 here comes from Tr(lambda_ij) = 2delta_ij
        hmc_float c_0_hat = 1/(1 + 64 * kappa_in * kappa_in);
        bc_tmp.re = 2./8. * kappa_in * c_0_hat * csw * TEMPORAL_RE;
        bc_tmp.im = 2./8. * kappa_in * c_0_hat * csw * TEMPORAL_IM;
        
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 0
        dir2 = 0;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_complex(v1, bc_tmp);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
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

        Matrix3x3 v1, v2;
        ae out_tmp;
        //this is used to save the BC-conditions...
        hmc_complex bc_tmp;
        int dir1, dir2;
        int global_link_pos;
        int n = pos.space;
        int t = pos.time;
        
        //go through the different directions
        ///////////////////////////////////
        // mu = 2
        ///////////////////////////////////
        dir1 = 2;
        global_link_pos = get_link_pos(dir, n, t);
        
        //the 2 here comes from Tr(lambda_ij) = 2delta_ij
        hmc_float c_0_hat = 1/(1 + 64 * kappa_in * kappa_in);
        bc_tmp.re = 2./8. * kappa_in * c_0_hat * csw * TEMPORAL_RE;
        bc_tmp.im = 2./8. * kappa_in * c_0_hat * csw * TEMPORAL_IM;
        
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 0
        dir2 = 0;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_complex(v1, bc_tmp);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
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

        Matrix3x3 v1, v2;
        ae out_tmp;
        //this is used to save the BC-conditions...
        hmc_complex bc_tmp;
        int dir1, dir2;
        int global_link_pos;
        int n = pos.space;
        int t = pos.time;
        
        //go through the different directions
        ///////////////////////////////////
        // mu = 3
        ///////////////////////////////////
        dir1 = 3;
        global_link_pos = get_link_pos(dir, n, t);
        
        //the 2 here comes from Tr(lambda_ij) = 2delta_ij
        hmc_float c_0_hat = 1/(1 + 64 * kappa_in * kappa_in);
        bc_tmp.re = 2./8. * kappa_in * c_0_hat * csw * TEMPORAL_RE;
        bc_tmp.im = 2./8. * kappa_in * c_0_hat * csw * TEMPORAL_IM;
        
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 0
        dir2 = 0;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        v2 = add_up_diagrams(field,X,X1,Y,Y1,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_complex(v1, bc_tmp);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}