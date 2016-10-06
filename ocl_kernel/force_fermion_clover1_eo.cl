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
 calculate the even-odd-preconditioned force of the clover term with eoprec spinorfields (cf. paper of Jansen, Liu: Implementation of Symanzik's improvement program for simulations of dynamical wilson fermions in lattice QCD, URL:http://arxiv.org/abs/hep-lat/9603008)
*/

//matrixsu3 which is given by i * Trace_Dirac[gamma_5*sigma_{mu,nu}* Y*X^dagger + gamma_5*sigma_{mu,nu}* X*Y^dagger]
//cf. equation (22) in Jansen Liu Paper
Matrixsu3 square(__global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, const site_idx pos_square, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 out;
    Matrix3x3 tmp1, tmp2;
    spinor x, y, u, v;
    x = getSpinor_eo(X, pos_square);
    y = getSpinor_eo(Y, pos_square);
    
    //calculate u/v = gamma_5 * simga_{mu,nu} * x/y
    u = sigma_mu_nu_times_spinor(x, dir1, dir2);
    u = gamma5_local(u);
    v = sigma_mu_nu_times_spinor(y, dir1, dir2);
    v = gamma5_local(v);
    
    tmp1 = tr_dirac_x_times_y_dagger(v.e0, v.e1, v.e2, v.e3, x.e0, x.e1, x.e2, x.e3);
    tmp2 = tr_dirac_x_times_y_dagger(u.e0, u.e1, u.e2, u.e3, y.e0, y.e1, y.e2, y.e3);
    
    tmp1 = add_matrix3x3(tmp1, tmp2);
    
    out = matrix_3x3tosu3(tmp1);
    out = multiply_matrixsu3_by_complex(out, hmc_complex_i);
    return out;
}


Matrixsu3 diagram1a_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh;
    site_idx idx_neigh_eo;
    
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);
    
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // square(x+nu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, Y, idx_neigh_eo, dir1, dir2, evenodd));
    
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    out = multiply_matrixsu3_dagger(out, U);
    return out;
}

Matrixsu3 diagram1a_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // square(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, Y, idx_neigh_eo, dir1, dir2, evenodd));
    
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);
    
    out = multiply_matrixsu3_by_real (out, -1.);
    return out;
}


Matrixsu3 diagram1b_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh;
    //site_idx idx_neigh_eo;
    
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);
    
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    out = multiply_matrixsu3_dagger(out, U);
    
    // square(x)
    site_idx idx_arg_eo = get_eo_site_idx_from_st_idx(idx_arg);
    out = multiply_matrixsu3(out, square(X, Y, idx_arg_eo, dir1, dir2, evenodd));
    return out;
}

Matrixsu3 diagram1b_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    //site_idx idx_neigh_eo;
    
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);
    
    // square(x)
    site_idx idx_arg_eo = get_eo_site_idx_from_st_idx(idx_arg);
    out = multiply_matrixsu3(out, square(X, Y, idx_arg_eo, dir1, dir2, evenodd));
    
    out = multiply_matrixsu3_by_real (out, -1.);
    return out;
}


Matrixsu3 diagram1c_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh;
    site_idx idx_neigh_eo;
    
    // square(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, Y, idx_neigh_eo, dir1, dir2, evenodd));
    
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);
    
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    out = multiply_matrixsu3_dagger(out, U);
    return out;
}

Matrixsu3 diagram1c_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // square(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, Y, idx_neigh_eo, dir1, dir2, evenodd));
    
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);
    
    out = multiply_matrixsu3_by_real (out, -1.);
    return out;
}



Matrixsu3 diagram1d_up(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);

    // square(x+mu+nu)
    idx_neigh1 = get_neighbor_from_st_idx(idx_arg, dir2); // x+nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x+nu)+mu
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, Y, idx_neigh_eo, dir1, dir2, evenodd));
    
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    out = multiply_matrixsu3_dagger(out, U);
    return out;
}

Matrixsu3 diagram1d_down(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    site_idx idx_neigh_eo;
    
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // square(x-nu+mu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    idx_neigh_eo = get_eo_site_idx_from_st_idx(idx_neigh);
    out = multiply_matrixsu3(out, square(X, Y, idx_neigh_eo, dir1, dir2, evenodd));
    
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);

    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);
    
    out = multiply_matrixsu3_by_real (out, -1.);
    return out;
}

Matrix3x3 add_up_diagrams(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
	Matrixsu3 tmp;
    Matrix3x3 out;
    out = zero_matrix3x3();
    
    tmp = diagram1a_up(field, X, Y, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(tmp), out);
    tmp = diagram1a_down(field, X, Y, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(tmp), out);
    tmp = diagram1b_up(field, X, Y, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(tmp), out);
    tmp = diagram1b_down(field, X, Y, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(tmp), out);
    tmp = diagram1c_up(field, X, Y, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(tmp), out);
    tmp = diagram1c_down(field, X, Y, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(tmp), out);
    tmp = diagram1d_up(field, X, Y, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(tmp), out);
    tmp = diagram1d_down(field, X, Y, idx_arg, dir1, dir2, evenodd);
    out = add_matrix3x3(matrix_su3to3x3(tmp), out);
    
    return out;
}


__kernel void fermion_force_clover1_eo_0(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);

        Matrix3x3 v1, v2;
        ae out_tmp;
        //note: no factor for boundary conditions and chemical potential because clover term is diagonal in the lattice points
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
	    hmc_float factor = 0.125 * kappa_in * c_0_hat * csw;
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_real(v1, factor);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}


__kernel void fermion_force_clover1_eo_1(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);

        Matrix3x3 v1, v2;
        ae out_tmp;
        //note: no factor for boundary conditions and chemical potential because clover term is diagonal in the lattice points
        int dir1, dir2;
        int global_link_pos;
        int n = pos.space;
        int t = pos.time;
        
        //go through the different directions
        ///////////////////////////////////
        // mu = 1
        ///////////////////////////////////
        dir1 = 1;
        global_link_pos = get_link_pos(dir1, n, t);
        
        //the 2 here comes from Tr(lambda_ij) = 2delta_ij
        hmc_float c_0_hat = 1/(1 + 64 * kappa_in * kappa_in);
	    hmc_float factor = 0.125 * kappa_in * c_0_hat * csw;
        
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 0
        dir2 = 0;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_real(v1, factor);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}

__kernel void fermion_force_clover1_eo_2(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);

        Matrix3x3 v1, v2;
        ae out_tmp;
        //note: no factor for boundary conditions and chemical potential because clover term is diagonal in the lattice points
        int dir1, dir2;
        int global_link_pos;
        int n = pos.space;
        int t = pos.time;
        
        //go through the different directions
        ///////////////////////////////////
        // mu = 2
        ///////////////////////////////////
        dir1 = 2;
        global_link_pos = get_link_pos(dir1, n, t);
        
        //the 2 here comes from Tr(lambda_ij) = 2delta_ij
        hmc_float c_0_hat = 1/(1 + 64 * kappa_in * kappa_in);
	    hmc_float factor = 0.125 * kappa_in * c_0_hat * csw;
        
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 0
        dir2 = 0;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_real(v1, factor);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}

__kernel void fermion_force_clover1_eo_3(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict X, __global const spinorStorageType * const restrict Y, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
        // caculate (pos,time) out of id_mem depending on evenodd
        // as we are positioning only on even or odd site we can update up- and downwards link without the danger of overwriting each other
        st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);

        Matrix3x3 v1, v2;
        ae out_tmp;
        //note: no factor for boundary conditions and chemical potential because clover term is diagonal in the lattice points
        int dir1, dir2;
        int global_link_pos;
        int n = pos.space;
        int t = pos.time;
        
        //go through the different directions
        ///////////////////////////////////
        // mu = 3
        ///////////////////////////////////
        dir1 = 3;
        global_link_pos = get_link_pos(dir1, n, t);
        
        //the 2 here comes from Tr(lambda_ij) = 2delta_ij
        hmc_float c_0_hat = 1/(1 + 64 * kappa_in * kappa_in);
	    hmc_float factor = 0.125 * kappa_in * c_0_hat * csw;
        
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 0
        dir2 = 0;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        v2 = add_up_diagrams(field,X,Y,pos,dir1,dir2,evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_real(v1, factor);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}
