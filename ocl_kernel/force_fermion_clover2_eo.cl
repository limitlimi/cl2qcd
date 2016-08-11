//
//calculate the even-odd-preconditioned force of the clover term
//


__kernel void fermion_force_clover2_eo_0(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
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
        global_link_pos = get_link_pos(dir, n, t);

        //go through the different directions
        ///////////////////////////////////
        // mu = 0
        ///////////////////////////////////
        dir1 = 0;
        //the 2 here comes from Tr(lambda_ij) = 2delta_ij
        //add factor c_sw*kappa/4
        bc_tmp.re = 2./4. * kappa_in * csw * TEMPORAL_RE;
        bc_tmp.im = 2./4. * kappa_in * csw * TEMPORAL_IM;
        
        //if chemical potential is activated, U has to be multiplied by appropiate factor
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
        
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        //diagram1a_up
        v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        
        
    }
}

Matrixsu3 diagram2a_up(__global const Matrixsu3StorageType * const restrict field, const st_index pos, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 out;
    int n = pos.space;
    int t = pos.time;
    U = get_matrixsu3(field, n, t, dir);
    return out;
}


//matrixsu3 which is given by Dirac-Trace(i*sigma_{mu,nu}*(1+T_ee)^{-1})
//cf. equation (22) in Jansen Liu Paper
Matrixsu3 triangle(__global const Matrixsu3StorageType * const restrict one_plus_T_ee_inverse)
{
    Matrixsu3 out;
    
    return out;
}



__kernel void fermion_force_clover2_eo_1(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);
        
    }
}

__kernel void fermion_force_clover2_eo_2(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);
        
    }
}

__kernel void fermion_force_clover2_eo_3(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);
        
    }
}