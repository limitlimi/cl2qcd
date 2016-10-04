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
  calculate the even-odd-preconditioned force of the clover term
*/

//matrixsu3 which is given by Dirac-Trace(i*sigma_{mu,nu}*(1+T_ee)^{-1})
//cf. equation (22) in Jansen Liu Paper
Matrixsu3 triangle(__global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, st_index pos_triangle, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 out;
    Matrix3x3 tmp1;
    Matrix6x6 B_plus, B_minus;
    B_plus = get6x6(C, get_site_idx(pos_triangle)); //upper left block of (1+T_ee)^(-1) = {{B_11, B_12},{B_21, B_22}}
    B_minus = get6x6(D, get_site_idx(pos_triangle)); //lower right block of (1+T_ee)^(-1) = {{B_33, B_34},{B_43, B_44}}
    if(dir1 ==0){
        if(dir2 == 1){ //B_12 + B_21 - B_34 - B_43
            tmp1 = get_3x3_block_upperright(B_plus);
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerleft(B_plus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_upperright(B_minus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerleft(B_minus));
        }
        else if(dir2 == 2){ //i * ( B_12 - B_21 - B_34 + B_43 )
            tmp1 = get_3x3_block_upperright(B_plus);
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerleft(B_plus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_upperright(B_minus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerleft(B_minus));
            tmp1 = multiply_matrix3x3_by_complex(tmp1, hmc_complex_i);
        }
        else if(dir2 == 3){//B_11 - B_22 - B_33 + B_44
            tmp1 = get_3x3_block_upperleft(B_plus);
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerright(B_plus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_upperleft(B_minus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerright(B_minus));
        }
    }
    else if(dir1 == 1){
        if(dir2 == 0){ //-B_12 - B_21 + B_34 + B_43 = (-1) * (dir1=0,dir2=1)
            tmp1 = get_3x3_block_upperright(B_plus);
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerleft(B_plus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_upperright(B_minus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerleft(B_minus));
            tmp1 = multiply_matrix3x3_by_real(tmp1 , -1.);
        }
        else if(dir2 == 2){ //-B_11 + B_22 - B_33 + B_44
            tmp1 = multiply_matrix3x3_by_real(get_3x3_block_upperleft(B_plus), -1.);
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerright(B_plus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_upperleft(B_minus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerright(B_minus));
        }
        else if(dir2 == 3){ //i * ( B_12 - B_21 + B_34 - B_43)
            tmp1 = get_3x3_block_upperright(B_plus);
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerleft(B_plus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_upperright(B_minus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerleft(B_minus));
            tmp1 = multiply_matrix3x3_by_complex(tmp1, hmc_complex_i);
        }
    }
    else if(dir1 == 2){
        if(dir2 == 0){ //-i * ( B_12 - B_21 - B_34 + B_43 ) = (-1) * (dir1=0,dir2=2)
            tmp1 = get_3x3_block_upperright(B_plus);
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerleft(B_plus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_upperright(B_minus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerleft(B_minus));
            tmp1 = multiply_matrix3x3_by_complex(tmp1, hmc_complex_minusi);
        }
        else if(dir2 == 1){ //B_11 - B_22 + B_33 - B_44
            tmp1 = get_3x3_block_upperleft(B_plus);
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerright(B_plus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_upperleft(B_minus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerright(B_minus));
        }
        else if(dir2 == 3){ //(-1) * (B_12 + B_21 + B_34 + B_43)
            tmp1 = get_3x3_block_upperright(B_plus);
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerleft(B_plus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_upperright(B_minus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerleft(B_minus));
            tmp1 = multiply_matrix3x3_by_real(tmp1, -1.);
        }
    }
    else if(dir1 == 3){
        if(dir2 == 0){//-B_11 + B_22 + B_33 - B_44
            tmp1 = multiply_matrix3x3_by_real(get_3x3_block_upperleft(B_plus), -1.);
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerright(B_plus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_upperleft(B_minus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerright(B_minus));
        }
        else if(dir2 == 1){//-i * ( B_12 - B_21 + B_34 - B_43) = (-1) * (dir1=1,dir2=3)
            tmp1 = get_3x3_block_upperright(B_plus);
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerleft(B_plus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_upperright(B_minus));
            tmp1 = subtract_matrix3x3(tmp1, get_3x3_block_lowerleft(B_minus));
            tmp1 = multiply_matrix3x3_by_complex(tmp1, hmc_complex_minusi);
        }
        else if(dir2 == 2){//B_12 + B_21 + B_34 + B_43
            tmp1 = get_3x3_block_upperright(B_plus);
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerleft(B_plus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_upperright(B_minus));
            tmp1 = add_matrix3x3(tmp1, get_3x3_block_lowerleft(B_minus));
        }
    }
    
    out = matrix_3x3tosu3(tmp1);
    out = multiply_matrixsu3_by_complex(out, hmc_complex_i);
    return out;
}

Matrixsu3 diagram2a_up(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh;
    
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);
    
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);

    // triangle(x+nu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    out = multiply_matrixsu3(out, triangle(C, D, idx_neigh, dir1, dir2));
    
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    out = multiply_matrixsu3_dagger(out, U);
    
    return out;
}

Matrixsu3 diagram2a_down(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // triangle(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    out = multiply_matrixsu3(out, triangle(C, D, idx_neigh, dir1, dir2));
    
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);

    out = multiply_matrixsu3_by_real (out, -1.);
    return out;
}


Matrixsu3 diagram2b_up(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh;
    
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
    
    // triangle(x)
    out = multiply_matrixsu3(out, triangle(C, D, idx_arg, dir1, dir2));
   
    return out;
}

Matrixsu3 diagram2b_down(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;

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
    
    // triangle(x)
    out = multiply_matrixsu3(out, triangle(C, D, idx_arg, dir1, dir2));
 
    out = multiply_matrixsu3_by_real (out, -1.);
    return out;
}


Matrixsu3 diagram2c_up(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh;
    
    // triangle(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    out = multiply_matrixsu3(out, triangle(C, D, idx_neigh, dir1, dir2));
    
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

Matrixsu3 diagram2c_down(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    
    // triangle(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    out = multiply_matrixsu3(out, triangle(C, D, idx_neigh, dir1, dir2));
    
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


Matrixsu3 diagram2d_up(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3(out, U);
    
    // triangle(x+mu+nu)
    idx_neigh1 = get_neighbor_from_st_idx(idx_arg, dir2); // x+nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x+nu)+mu
    out = multiply_matrixsu3(out, triangle(C, D, idx_neigh, dir1, dir2));
    
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    out = multiply_matrixsu3_dagger(out, U);
    
    return out;
}

Matrixsu3 diagram2d_down(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    Matrixsu3 U, out;
    out = unit_matrixsu3();
    st_idx idx_neigh, idx_neigh1;
    
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    out = multiply_matrixsu3_dagger(out, U);
    
    // triangle(x-nu+mu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    out = multiply_matrixsu3(out, triangle(C, D, idx_neigh, dir1, dir2));
    
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

Matrix3x3 add_up_diagrams2(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, const st_index idx_arg, const dir_idx dir1, const dir_idx dir2, int evenodd)
{
    Matrix3x3 out;
    Matrixsu3 tmp;
    out = zero_matrix3x3();
    
    if(evenodd == ODD){ //add up only diagrams (a) and (c)
        tmp = diagram2a_up(field, C, D, idx_arg, dir1, dir2);
        out = add_matrix3x3(matrix_su3to3x3(tmp), out);
        tmp = diagram2a_down(field, C, D, idx_arg, dir1, dir2);
        out = add_matrix3x3(matrix_su3to3x3(tmp), out);
        tmp = diagram2c_up(field, C, D, idx_arg, dir1, dir2);
        out = add_matrix3x3(matrix_su3to3x3(tmp), out);
        tmp = diagram2c_down(field, C, D, idx_arg, dir1, dir2);
        out = add_matrix3x3(matrix_su3to3x3(tmp), out);
    }
    else if(evenodd == EVEN){ //add up only diagrams (b) and (d)
        tmp = diagram2b_up(field, C, D, idx_arg, dir1, dir2);
        out = add_matrix3x3(matrix_su3to3x3(tmp), out);
        tmp = diagram2b_down(field, C, D, idx_arg, dir1, dir2);
        out = add_matrix3x3(matrix_su3to3x3(tmp), out);
        tmp = diagram2d_up(field, C, D, idx_arg, dir1, dir2);
        out = add_matrix3x3(matrix_su3to3x3(tmp), out);
        tmp = diagram2d_down(field, C, D, idx_arg, dir1, dir2);
        out = add_matrix3x3(matrix_su3to3x3(tmp), out);
    }
    
    return out;
}

//input C/D = upperright/lowerleft-block of (1+T)^(-1)
__kernel void fermion_force_clover2_eo_0(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
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
        //add factor 2*c_sw*kappa/4
	hmc_float factor = 0.5 * kappa_in * csw;
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_real(v1, factor);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}

__kernel void fermion_force_clover2_eo_1(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
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
        // mu = 1
        ///////////////////////////////////
        dir1 = 1;
        global_link_pos = get_link_pos(dir1, n, t);
        
        
        //the 2 here comes from Tr(lambda_ij) = 2delta_ij
        //add factor 2*c_sw*kappa/4
	hmc_float factor = 0.5 * kappa_in * csw;
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 0
        dir2 = 0;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_real(v1, factor);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}

__kernel void fermion_force_clover2_eo_2(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
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
        //add factor 2*c_sw*kappa/4
	hmc_float factor = 0.5 * kappa_in * csw;
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 0
        dir2 = 0;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 3
        dir2 = 3;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_real(v1, factor);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}

__kernel void fermion_force_clover2_eo_3(__global const Matrixsu3StorageType * const restrict field, __global const Matrix6x6StorageType * const restrict C, __global const Matrix6x6StorageType * const restrict D, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in, hmc_float csw)
{
    // must include HALO, as we are updating neighbouring sites
    // -> not all local sites will fully updated if we don't calculate on halo indices, too
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
        //add factor 2*c_sw*kappa/4
	hmc_float factor = 0.5 * kappa_in * csw;
        
        v1 = zero_matrix3x3();
        //add up diagrams for all nu unequal mu
        /////////////////////////////////
        // nu = 0
        dir2 = 0;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 1
        dir2 = 1;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        /////////////////////////////////
        // nu = 2
        dir2 = 2;
        v2 = add_up_diagrams2(field, C, D, pos, dir1, dir2, evenodd);
        v1 = add_matrix3x3(v1, v2);
        
        v1 = multiply_matrix3x3_by_real(v1, factor);
        out_tmp = tr_lambda_u(v1);
        update_gaugemomentum(out_tmp, 1., global_link_pos, out);
    }
}
