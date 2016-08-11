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
 clover-term fermion matrix 1+T for eoprec spinorfields (c.f. paper of Jansen, Liu: Implementation of Symanzik's improvement program for simulations of dynamical wilson fermions in lattice QCD, URL:http://arxiv.org/abs/hep-lat/9603008)
 T=i/2*c_sw*kappa*sigma_{\mu \nu}*F_{\mu \nu}*delta_xy, with F=field strength tensor on the lattice
sigma_{\mu \nu} = [gamma_\mu, gamma_\nu]
*/

spinor clover_eoprec_unified_local(__global const spinorStorageType * const restrict in, __global Matrixsu3StorageType  const * const restrict field, const st_idx idx_arg, const dir_idx dir, hmc_float kappa_in, hmc_float csw)
{
    //this is used to save the idx of the neighbors
    st_idx idx_neigh;
    
    spinor out_tmp, plus;
    //site_idx nn_eo;
    dir_idx dir2;
    su3vec psi0, psi1, psi2, psi3, temp;
    Matrixsu3 U;
    //this is used to save the BC-conditions...
    hmc_complex bc_tmp = (dir == TDIR) ? (get_neighb) {//add c_sw here? /// factor 2 because of (1+T)~2*delta_xy // factor i/2 ??
        kappa_in * TEMPORAL_RE, kappa_in * TEMPORAL_IM
    } :
    (hmc_complex) {
        kappa_in * SPATIAL_RE, kappa_in * SPATIAL_IM
    };
    out_tmp = set_spinor_zero();

    //transform normal indices to eoprec index
    pos_eo = get_eo_site_idx_from_st_idx(idx_arg);
    plus = getSpinor_eo(in, pos_eo);
    U = getSU3(field, get_link_idx(dir, idx_arg));
    if(dir == TDIR) {
        /////////////////////////////////
        // nu = 0
        dir2 = TDIR;
        /////////////////////////////////
        //Calculate (1+sigma_00) * plus
        //with 1+sigma_00:
        //| 1  0  0  0  |       | plus.e0 |
        //| 0  1  0  0  | psi = | plus.e1 |
        //| 0  0  1  0  |       | plus.e2 |
        //| 0  0  0  1  |       | plus.e3 |
        /////////////////////////////////
        psi0 = plus.e0;
        psi1 = plus.e1;
        psi2 = plus.e2;
        psi3 = plus.e3;
        //calculate (1+field-strength-tensor) * su3vec
        tmp = field_strength_tensor_times_su3vec(psi0, field, idx_arg, dir, dir2);
        out_tmp.e0 = su3vec_acc(psi0, tmp);
        tmp = field_strength_tensor_times_su3vec(psi1, field, idx_arg, dir, dir2);
        out_tmp.e1 = su3vec_acc(psi1, tmp);
        tmp = field_strength_tensor_times_su3vec(psi2, field, idx_arg, dir, dir2);
        out_tmp.e2 = su3vec_acc(psi2, tmp);
        tmp = field_strength_tensor_times_su3vec(psi3, field, idx_arg, dir, dir2);
        out_tmp.e3 = su3vec_acc(psi3, tmp);
        
        /////////////////////////////////
        // nu = 1
        dir2 = XDIR;
        /////////////////////////////////
        //Calculate (1+sigma_01) * plus
        //with sigma_01:
        //| 1     -2*i  0    0  |       | 2*(-i)*plus.e1+plus.e0 |
        //| -2*i  1     0    0  | psi = | 2*(-i)*plus.e0+plus.e1 |
        //| 0     0     1    2*i|       |    2*i*plus.e3+plus.e2 |
        //| 0     0     2*i  1  |       |    2*i*plus.e2+plus.e3 |
        /////////////////////////////////
        psi0 = su3vec_times_complex(plus.e1,{0,-2});
        psi1 = su3vec_times_complex(plus.e0,{0,-2});
        psi2 = su3vec_times_complex(plus.e3,{0,2});
        psi3 = su3vec_times_complex(plus.e2,{0,2});
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate (1+field-strength-tensor) * su3vec
        tmp = field_strength_tensor_times_su3vec(psi0, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi0, tmp);
        out_tmp.e0 = su3vec_acc(temp, out_tmp.e0);
        tmp = field_strength_tensor_times_su3vec(psi1, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi1, tmp);
        out_tmp.e1 = su3vec_acc(temp, out_tmp.e1);
        tmp = field_strength_tensor_times_su3vec(psi2, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi2, tmp);
        out_tmp.e2 = su3vec_acc(temp, out_tmp.e2);
        tmp = field_strength_tensor_times_su3vec(psi3, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi3, tmp);
        out_tmp.e3 = su3vec_acc(temp, out_tmp.e3);
        
        /////////////////////////////////
        // nu = 2
        dir2 = YDIR;
        /////////////////////////////////
        //Calculate (1+sigma_02) * plus
        //with 1+sigma_02:
        //| 1   -2  0   0  |       | -2*plus.e1+plus.e0 |
        //| 2   1   0   0  | psi = |  2*plus.e0+plus.e1 |
        //| 0   0   1   2  |       |  2*plus.e3+plus.e2 |
        //| 0   0   -2  1  |       | -2*plus.e2+plus.e3 |
        /////////////////////////////////
        psi0 = su3vec_times_real(plus.e1,-2);
        psi1 = su3vec_times_real(plus.e0,2);
        psi2 = su3vec_times_real(plus.e3,2);
        psi3 = su3vec_times_real(plus.e2,-2);
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate (1+field-strength-tensor) * su3vec
        tmp = field_strength_tensor_times_su3vec(psi0, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi0, tmp);
        out_tmp.e0 = su3vec_acc(temp, out_tmp.e0);
        tmp = field_strength_tensor_times_su3vec(psi1, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi1, tmp);
        out_tmp.e1 = su3vec_acc(temp, out_tmp.e1);
        tmp = field_strength_tensor_times_su3vec(psi2, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi2, tmp);
        out_tmp.e2 = su3vec_acc(temp, out_tmp.e2);
        tmp = field_strength_tensor_times_su3vec(psi3, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi3, tmp);
        out_tmp.e3 = su3vec_acc(temp, out_tmp.e3);
        
        
        /////////////////////////////////
        // nu = 3
        dir = ZDIR;
        /////////////////////////////////
        //Calculate 1+sigma_03 * plus
        //with 1+sigma_03:
        //| -2*i+1  0       0      0   |       | 2*(-i)*plus.e0+plus.e0 |
        //| 0       2*i+1   0      0   | psi = |    2*i*plus.e1+plus.e1 |
        //| 0       0       2*i+1  0   |       |    2*i*plus.e2+plus.e2 |
        //| 0       0       0      -2*i|       | 2*(-i)*plus.e3+plus.e3 |
        /////////////////////////////////
        psi0 = su3vec_times_complex(plus.e0,{0,-2});
        psi1 = su3vec_times_complex(plus.e1,{0,2});
        psi2 = su3vec_times_complex(plus.e2,{0,2});
        psi3 = su3vec_times_complex(plus.e3,{0,-2});
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate (1+field-strength-tensor) * su3vec
        tmp = field_strength_tensor_times_su3vec(psi0, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi0, tmp);
        out_tmp.e0 = su3vec_acc(temp, out_tmp.e0);
        tmp = field_strength_tensor_times_su3vec(psi1, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi1, tmp);
        out_tmp.e1 = su3vec_acc(temp, out_tmp.e1);
        tmp = field_strength_tensor_times_su3vec(psi2, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi2, tmp);
        out_tmp.e2 = su3vec_acc(temp, out_tmp.e2);
        tmp = field_strength_tensor_times_su3vec(psi3, field, idx_arg, dir, dir2);
        tmp = su3vec_acc(psi3, tmp);
        out_tmp.e3 = su3vec_acc(temp, out_tmp.e3);
        
    }
    if(dir == XDIR) {
        /////////////////////////////////
        // nu = 0
        /////////////////////////////////
        //Calculate (1+sigma_10) * plus
        //with 1+sigma_10:
        //| 1   2*i 0    0   |       |     2*i*plus.e1+plus.e0  |
        //| 2*i 1   0    0   | psi = |     2*i*plus.e0+plus.e1  |
        //| 0   0   1    -2*i|       |  2*(-i)*plus.e3+plus.e2  |
        //| 0   0   -2*i 1   |       |  2*(-i)*plus.e2+plus.e3  |
        /////////////////////////////////
        psi0 = su3vec_times_complex(plus.e1,{0,2});
        psi1 = su3vec_times_complex(plus.e0,{0,2});
        psi2 = su3vec_times_complex(plus.e3,{0,-2});
        psi3 = su3vec_times_complex(plus.e2,{0,-2});
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate field-strength-tensor * spinor
        
        
        /////////////////////////////////
        // nu = 1
        /////////////////////////////////
        //Calculate (1+sigma_11) * plus
        //with 1+sigma_11:
        //| 1  0  0  0  |       | plus.e0 |
        //| 0  1  0  0  | psi = | plus.e1 |
        //| 0  0  1  0  |       | plus.e2 |
        //| 0  0  0  1  |       | plus.e3 |
        /////////////////////////////////
        psi0 = plus.e0;
        psi1 = plus.e1;
        psi2 = plus.e2;
        psi3 = plus.e3;
        //calculate field-strength-tensor * spinor
        
        
        /////////////////////////////////
        // nu = 2
        /////////////////////////////////
        //Calculate sigma_12 * plus
        //with sigma_12:
        //| 2*i+1   0       0      0      |       |    2*i*plus.e0+plus.e0 |
        //| 0       -2*i+1  0      0      | psi = | 2*(-i)*plus.e1+plus.e1 |
        //| 0       0       2*i+1  0      |       |    2*i*plus.e2+plus.e2 |
        //| 0       0       0      -2*i+1 |       | 2*(-i)*plus.e3+plus.e3 |
        /////////////////////////////////
        psi0 = su3vec_times_complex(plus.e0,{0,2});
        psi1 = su3vec_times_complex(plus.e1,{0,-2});
        psi2 = su3vec_times_complex(plus.e2,{0,2});
        psi3 = su3vec_times_complex(plus.e3,{0,-2});
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate field-strength-tensor * spinor
        
        /////////////////////////////////
        // nu = 3
        /////////////////////////////////
        //Calculate (1+sigma_13) * plus
        //with 1+sigma_13:
        //| 1   -2  0   0  |       | -2*plus.e1+plus.e0 |
        //| 2   1   0   0  | psi = |  2*plus.e0+plus.e1 |
        //| 0   0   1   -2 |       | -2*plus.e3+plus.e2 |
        //| 0   0   2   1  |       |  2*plus.e2+plus.e3 |
        /////////////////////////////////
        psi0 = su3vec_times_real(plus.e1,-2);
        psi1 = su3vec_times_real(plus.e0,2);
        psi2 = su3vec_times_real(plus.e3,-2);
        psi3 = su3vec_times_real(plus.e2,2);
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate field-strength-tensor * spinor
        
    }
    if(dir == YDIR) {
        /////////////////////////////////
        // nu = 0
        /////////////////////////////////
        //Calculate (1+sigma_20) * plus
        //with 1+sigma_20:
        //| 1   2   0   0  |       |  2*plus.e1+plus.e0 |
        //| -2  1   0   0  | psi = | -2*plus.e0+plus.e1 |
        //| 0   0   1   -2 |       | -2*plus.e3+plus.e2 |
        //| 0   0   2   1  |       |  2*plus.e2+plus.e3 |
        /////////////////////////////////
        psi0 = su3vec_times_real(plus.e1,2);
        psi1 = su3vec_times_real(plus.e0,-2);
        psi2 = su3vec_times_real(plus.e3,-2);
        psi3 = su3vec_times_real(plus.e2,2);
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate field-strength-tensor * spinor
        
        /////////////////////////////////
        // nu = 1
        /////////////////////////////////
        //Calculate (1+sigma_21) * plus
        //with 1+sigma_21:
        //| -2*i+1  0       0       0     |       | 2*(-i)*plus.e0+plus.e0 |
        //| 0       2*i+1   0       0     | psi = |    2*i*plus.e1+plus.e1 |
        //| 0       0       -2*i+1  0     |       | 2*(-i)*plus.e2+plus.e2 |
        //| 0       0       0       2*i+1 |       |    2*i*plus.e3+plus.e3 |
        /////////////////////////////////
        psi0 = su3vec_times_complex(plus.e0,{0,2});
        psi1 = su3vec_times_complex(plus.e1,{0,-2});
        psi2 = su3vec_times_complex(plus.e2,{0,2});
        psi3 = su3vec_times_complex(plus.e3,{0,-2});
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate field-strength-tensor * spinor
        
        
        /////////////////////////////////
        // nu = 2
        /////////////////////////////////
        //Calculate (1+sigma_22) * plus
        //with 1+sigma_22:
        //| 1  0  0  0  |       | plus.e0 |
        //| 0  1  0  0  | psi = | plus.e1 |
        //| 0  0  1  0  |       | plus.e2 |
        //| 0  0  0  1  |       | plus.e3 |
        /////////////////////////////////
        psi0 = plus.e0;
        psi1 = plus.e1;
        psi2 = plus.e2;
        psi3 = plus.e3;
        //calculate field-strength-tensor * spinor
        
        
        /////////////////////////////////
        // nu = 3
        /////////////////////////////////
        //Calculate (1+sigma_23) * plus
        //with 1+sigma_23:
        //| 1   2*i 0    0   |       |  2*i*plus.e1+plus.e0  |
        //| 2*i 1   0    0   | psi = |  2*i*plus.e0+plus.e1  |
        //| 0   0   1    2*i |       |  2*i*plus.e3+plus.e2  |
        //| 0   0   2*i  1   |       |  2*i*plus.e2+plus.e3  |
        /////////////////////////////////
        psi0 = su3vec_times_complex(plus.e1,{0,2});
        psi1 = su3vec_times_complex(plus.e0,{0,2});
        psi2 = su3vec_times_complex(plus.e3,{0,2});
        psi3 = su3vec_times_complex(plus.e2,{0,2});
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate field-strength-tensor * spinor
        
        
    }
    if(dir == ZDIR) {
        /////////////////////////////////
        // nu = 0
        /////////////////////////////////
        //Calculate (1+sigma_30) * plus
        //with 1+sigma_30:
        //| 2*i+1   0       0       0     |       |    2*i*plus.e0+plus.e0 |
        //| 0       -2*i+1  0       0     | psi = | 2*(-i)*plus.e1+plus.e1 |
        //| 0       0       -2*i+1  0     |       | 2*(-i)*plus.e2+plus.e2 |
        //| 0       0       0       2*i+1 |       |    2*i*plus.e3+plus.e3 |
        /////////////////////////////////
        psi0 = su3vec_times_complex(plus.e0,{0,2});
        psi1 = su3vec_times_complex(plus.e1,{0,-2});
        psi2 = su3vec_times_complex(plus.e2,{0,-2});
        psi3 = su3vec_times_complex(plus.e3,{0,2});
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate field-strength-tensor * spinor
        
        
        /////////////////////////////////
        // nu = 1
        /////////////////////////////////
        //Calculate (1+sigma_31) * plus
        //with 1+sigma_31:
        //| 1   2   0    0  |       |  2*plus.e1+plus.e0 |
        //| -2  1   0    0  | psi = | -2*plus.e0+plus.e1 |
        //| 0   0   1    2  |       |  2*plus.e3+plus.e2 |
        //| 0   0   -2   1  |       | -2*plus.e2+plus.e3 |
        /////////////////////////////////
        psi0 = su3vec_times_real(plus.e1,2);
        psi1 = su3vec_times_real(plus.e0,-2);
        psi2 = su3vec_times_real(plus.e3,2);
        psi3 = su3vec_times_real(plus.e2,-2);
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate field-strength-tensor * spinor
     
        
        /////////////////////////////////
        // nu = 2
        /////////////////////////////////
        //Calculate (1+sigma_32) * plus
        //with 1+sigma_32:
        //| 1    -2*i 0     0    |       |  -2*i*plus.e1+plus.e0  |
        //| -2*i 1    0     0    | psi = |  -2*i*plus.e0+plus.e1  |
        //| 0    0    1     -2*i |       |  -2*i*plus.e3+plus.e2  |
        //| 0    0    -2*i  1    |       |  -2*i*plus.e2+plus.e3  |
        /////////////////////////////////
        psi0 = su3vec_times_complex(plus.e1,{0,-2});
        psi1 = su3vec_times_complex(plus.e0,{0,-2});
        psi2 = su3vec_times_complex(plus.e3,{0,-2});
        psi3 = su3vec_times_complex(plus.e2,{0,-2});
        psi0 = su3vec_acc(psi0,plus.e0);
        psi1 = su3vec_acc(psi1,plus.e1);
        psi2 = su3vec_acc(psi2,plus.e2);
        psi3 = su3vec_acc(psi3,plus.e3);
        //calculate field-strength-tensor * spinor
        
        
        /////////////////////////////////
        // nu = 3
        /////////////////////////////////
        //Calculate (1+sigma_3) * plus
        //with 1+sigma_33:
        //| 1  0  0  0  |       | plus.e0 |
        //| 0  1  0  0  | psi = | plus.e1 |
        //| 0  0  1  0  |       | plus.e2 |
        //| 0  0  0  1  |       | plus.e3 |
        /////////////////////////////////
        psi0 = plus.e0;
        psi1 = plus.e1;
        psi2 = plus.e2;
        psi3 = plus.e3;
        //calculate field-strength-tensor * spinor
        
    }
    out_tmp = spinor_times_complex(out_temp, bc_temp);
    return out_tmp;
}

void field_strength_tensor_times_su3vec(su3vec inout, __global Matrixsu3StorageType  const * const restrict field, const st_idx idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    //calculation of the lattice-field-strength-tensor according to the paper of Jansen and Liu, equation (6)
    //it consits of a sum of 8 terms which are products of link variables
    //psi1 = 1.term * inout
    //dir1 = mu, dir2 = nu
    st_idx idx_neigh, idx_neigh1;
    Matrixsu3 U;
    //////////////////////
    // 1.term = U_mu(x) * U_nu(x+mu) * U_mu(x+nu)^dagger * U_nu(x)^dagger
    su3vec psi1 = inout;
    //////////////////////
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    psi1 = su3matrix_dagger_times_su3vec(U, psi1);
    //////////////////////
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi1 = su3matrix_dagger_times_su3vec(U, psi1);
    //////////////////////
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi1 = su3matrix_times_su3vec(U, psi1);
    //////////////////////
    // U_mu(x)
    U = getSU3(field, get_link_idx(dir1, idx_arg));
    psi1 = su3matrix_times_su3vec(U, psi1);
    
    //////////////////////
    // 2.term = U_nu(x) * U_mu(x+nu-mu)^dagger * U_nu(x-nu)^dagger * U_mu(x-mu)
    su3vec psi2 = inout;
    //////////////////////
    // U_mu(x-mu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi2 = su3matrix_times_su3vec(U, psi2);
    //////////////////////
    // U_nu(x_nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi2 = su3matrix_dagger_times_su3vec(U, psi2);
    //////////////////////
    // U_mu(x+nu-mu)^dagger
    idx_neigh1 = get_neighbor_from_st_idx(idx_arg, dir2); // x+nu
    idx_neigh = get_lower_neighbor_from_st_idx(idx_neigh1, dir1); //(x+nu)-mu
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi2 = su3matrix_dagger_times_su3vec(U, psi2);
    //////////////////////
    // U_nu(x)
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    psi2 = su3matrix_times_su3vec(U, psi2);
    
    //////////////////////
    // 3.term = U_mu(x-mu)^dagger * U_nu(x-mu-nu)^dagger * U_mu(x-mu-nu) * U_nu(x-nu)
    su3vec psi3 = inout;
    //////////////////////
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi3 = su3matrix_times_su3vec(U, psi3);
    //////////////////////
    // U_mu(x-mu-nu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir1); // x-mu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi3 = su3matrix_times_su3vec(U, psi3);
    //////////////////////
    // U_nu(x-mu-nu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir1); // x-mu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi3 = su3matrix_dagger_times_su3vec(U, psi3);
    //////////////////////
    // U_mu(x-mu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi3 = su3matrix_dagger_times_su3vec(U, psi3);
    
    ///////////////////////
    // 4.term = U_nu(x-nu)^dagger * U_mu(x-nu) * U_nu(x-nu+mu) * U_mu(x)^dagger
    su3vec psi4 = inout;
    // U_mu(x)^dagger
    U = getSU3(field, get_link_idx(dir1, idx_arg));
    psi4 = su3matrix_dagger_times_su3vec(U, psi4);
    ////////////////////////
    // U_nu(x-nu+mu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi4 = su3matrix_times_su3vec(U, psi4);
    ////////////////////////
    // U_mu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi4 = su3matrix_times_su3vec(U, psi4);
    ////////////////////////
    // U_nu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi4 = su3matrix_dagger_times_su3vec(u, psi4);
    
    ////////////////////////
    // 5.term = U_nu(x) * U_mu(x+nu) * U_nu(x+mu)^dagger * U_mu(x)^dagger
    su3vec psi5 = inout;
    // U_mu(x)^dagger
    U = getSU3(field, get_link_idx(dir1, idx_arg));
    psi5 = su3matrix_dagger_times_su3vec(U, psi5);
    ////////////////////////
    // U_nu(x+mu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi5 = su3matrix_dagger_times_su3vec(U, psi5);
    ////////////////////////
    // U_mu(x+nu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi5 = su3matrix_times_su3vec(U, psi5);
    /////////////////////////
    // U_nu(x)
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    psi5 = su3matrix_times_su3vec(U, psi5);
    
    //////////////////////////
    // 6.term = U_mu(x-mu)^dagger * U_nu(x-mu) * U_mu(x+nu-mu) * U_nu(x)^dagger
    su3vec psi6 = inout;
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    psi6 = su3matrix_dagger_times_su3vec(U, psi6);
    //////////////////////////
    // U_mu(x+nu-mu)
    idx_neigh1 = get_neighbor_from_st_idx(idx_arg, dir2); // x+nu
    idx_neigh = get_lower_neighbor_from_st_idx(idx_neigh1, dir1); //(x+nu)-mu
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi6 = su3matrix_times_su3vec(U, psi6);
    /////////////////////////
    // U_nu(x-mu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi6 = su3matrix_times_su3vec(U, psi6);
    // U_mu(x-mu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi6 = su3matrix_dagger_times_su3vec(U, psi6);
    
    ////////////////////////////
    // 7.term = U_nu(x-nu)^dagger * U_mu(x-nu-mu)^dagger * U_nu(x-nu-mu) * U_mu(x-mu)
    su3vec psi7 = inout;
    // U_mu(x-mu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi7 = su3matrix_times_su3vec(U, psi7);
    // U_nu(x-nu-mu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir1); // x-mu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi7 = su3matrix_times_su3vec(U, psi7);
    // U_mu(x-nu-mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir1); // x-mu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi7 = su3matrix_dagger_times_su3vec(U, psi7);
    // U_nu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg,dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi7 = su3matrix_dagger_times_su3vec(U, psi7);
    
    ///////////////////////////
    // 8.term = U_mu(x) * U_nu(x-nu+mu)^dagger * U_mu(x-nu)^dagger * U_nu(x-nu)
    su3vec psi8 = inout;
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi8 = su3matrix_times_su3vec(u, psi8);
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    psi8 = su3matrix_dagger_times_su3vec(U, psi8);
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    psi8 = su3matrix_dagger_times_su3vec(U, psi8);
    // U_mu(x)
    U = getSU3(field, get_link_idx(dir1, idx_arg));
    psi8 = su3matrix_times_su3vec(U, psi8);
    
    
    //add psi1,...,psi8 and multiply by factor 1/8
    hmc_float factor = 1/8;
    inout = psi1;
    inout = su3vec_acc(inout, psi2);
    inout = su3vec_acc(inout, psi3);
    inout = su3vec_acc(inout, psi4);
    inout = su3vec_acc(inout, psi5);
    inout = su3vec_acc(inout, psi6);
    inout = su3vec_acc(inout, psi7);
    inout = su3vec_acc(inout, psi8);
    inout = su3vec_times_real(inout, factor);
}

