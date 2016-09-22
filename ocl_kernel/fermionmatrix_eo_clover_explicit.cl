/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra, 
 * Max Theilig
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

/* The clover Matrix (1+T) is a complex 12x12 Matrix for every lattice site
   T = c_sw sum_{mu,nu=0}^3 i/4 sigma_{mu,nu} F_{mu,nu} with F = lattice-field-strength-tensor
   because of sigma_{mu,nu} (1 + T) is blockdiagonal and one can store it as two complex 6x6 matrices
   the implementation is done according to OpenQCD documentation section 4.1
 */

Matrix3x3 field_strength_tensor(__global Matrixsu3StorageType  const * const restrict field, const st_idx idx_arg, const dir_idx dir1, const dir_idx dir2)
{
    //calculation of the lattice-field-strength-tensor according to the paper of Jansen and Liu, equation (6)
    //it consits of a sum of 8 terms which are products of link variables
    //dir1 = mu, dir2 = nu
    st_idx idx_neigh, idx_neigh1;
    Matrixsu3 U, tmp;
    Matrix3x3 out = zero_matrix3x3();
    
    
    //////////////////////
    // 1.term = U_mu(x) * U_nu(x+mu) * U_mu(x+nu)^dagger * U_nu(x)^dagger
    //////////////////////
    // U_mu(x)
    U = getSU3(field, get_link_idx(dir1, idx_arg));
    tmp = U;
    //////////////////////
    // U_nu(x+mu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    //////////////////////
    // U_mu(x+nu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    //////////////////////
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    /////////////////////
    out = add_matrix3x3(out, matrix_su3to3x3(tmp));
    
    
    //////////////////////
    // 2.term = U_nu(x) * U_mu(x+nu-mu)^dagger * U_nu(x-nu)^dagger * U_mu(x-mu)
    //////////////////////
    // U_nu(x)
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    tmp = U;
    //////////////////////
    // U_mu(x+nu-mu)^dagger
    idx_neigh1 = get_neighbor_from_st_idx(idx_arg, dir2); // x+nu
    idx_neigh = get_lower_neighbor_from_st_idx(idx_neigh1, dir1); //(x+nu)-mu
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    //////////////////////
    // U_nu(x-mu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    //////////////////////
    // U_mu(x-mu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    /////////////////////
    out = add_matrix3x3(out, matrix_su3to3x3(tmp));

    
    //////////////////////
    // 3.term = U_mu(x-mu)^dagger * U_nu(x-mu-nu)^dagger * U_mu(x-mu-nu) * U_nu(x-nu)
    //////////////////////
    // U_mu(x-mu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = adjoint_matrixsu3(U);
    //////////////////////
    // U_nu(x-mu-nu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir1); // x-mu
    idx_neigh = get_lower_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    //////////////////////
    // U_mu(x-mu-nu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir1); // x-mu
    idx_neigh = get_lower_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    //////////////////////
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    /////////////////////
    out = add_matrix3x3(out, matrix_su3to3x3(tmp));
    
    
    ///////////////////////
    // 4.term = U_nu(x-nu)^dagger * U_mu(x-nu) * U_nu(x-nu+mu) * U_mu(x)^dagger
    ////////////////////////
    // U_nu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = adjoint_matrixsu3(U);
    ////////////////////////
    // U_mu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    ////////////////////////
    // U_nu(x-nu+mu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    ////////////////////////
    // U_mu(x)^dagger
    U = getSU3(field, get_link_idx(dir1, idx_arg));
    tmp = multiply_matrixsu3_dagger(tmp ,U);
    /////////////////////
    out = add_matrix3x3(out, matrix_su3to3x3(tmp));


    ////////////////////////
    // 5.term = U_nu(x) * U_mu(x+nu) * U_nu(x+mu)^dagger * U_mu(x)^dagger
    /////////////////////////
    // U_nu(x)
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    tmp = U;
    ////////////////////////
    // U_mu(x+nu)
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    ////////////////////////
    // U_nu(x+mu)^dagger
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    ////////////////////////
    // U_mu(x)^dagger
    U = getSU3(field, get_link_idx(dir1, idx_arg));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    /////////////////////
    out = subtract_matrix3x3(out, matrix_su3to3x3(tmp));
    
    
    //////////////////////////
    // 6.term = U_mu(x-mu)^dagger * U_nu(x-mu) * U_mu(x+nu-mu) * U_nu(x)^dagger
    /////////////////////////
    // U_mu(x-mu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = adjoint_matrixsu3(U);
    /////////////////////////
    // U_nu(x-mu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    //////////////////////////
    // U_mu(x+nu-mu)
    idx_neigh1 = get_neighbor_from_st_idx(idx_arg, dir2); // x+nu
    idx_neigh = get_lower_neighbor_from_st_idx(idx_neigh1, dir1); //(x+nu)-mu
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    /////////////////////////
    // U_nu(x)^dagger
    U = getSU3(field, get_link_idx(dir2, idx_arg));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    /////////////////////
    out = subtract_matrix3x3(out, matrix_su3to3x3(tmp));
    
    
    ////////////////////////////
    // 7.term = U_nu(x-nu)^dagger * U_mu(x-nu-mu)^dagger * U_nu(x-nu-mu) * U_mu(x-mu)
    /////////////////////////
    // U_nu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg,dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = adjoint_matrixsu3(U);
    /////////////////////////
    // U_mu(x-nu-mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir1); // x-mu
    idx_neigh = get_lower_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    /////////////////////////
    // U_nu(x-nu-mu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir1); // x-mu
    idx_neigh = get_lower_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    ////////////////////////
    // U_mu(x-mu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    /////////////////////
    out = subtract_matrix3x3(out, matrix_su3to3x3(tmp));

    
    ///////////////////////////
    // 8.term = U_mu(x) * U_nu(x-nu+mu)^dagger * U_mu(x-nu)^dagger * U_nu(x-nu)
    /////////////////////////
    // U_mu(x)
    U = getSU3(field, get_link_idx(dir1, idx_arg));
    tmp = U;
    /////////////////////////
    // U_nu(x-nu+mu)^dagger
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir2); // x-nu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir1); // (x-nu)+mu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    /////////////////////////
    // U_mu(x-nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    /////////////////////////
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    /////////////////////
    out = subtract_matrix3x3(out, matrix_su3to3x3(tmp));

    return out;
}

/*the upper-left block is given by c_sw * i/16 * sum{k=1}^3 sigma_k (E_k - B_k)
  with sigma_k Pauli matrices, E_k = 8 * F_{0,k}, B_k = sum_{j,l=1}^3 4 * epsilon_{k,l,j} F_{l,j}
 */
Matrix6x6 clover_eoprec_unified_local_upper_left_block(__global Matrixsu3StorageType  const * const restrict field, const st_idx idx_arg, hmc_float csw)
{
    Matrix6x6 out = zero_matrix6x6();
    Matrix3x3 EB1, EB2, EB3, E, B, tmp, tmp1;
    //no bc_tmp needed
    hmc_complex factor = {0., 1./16. * csw};
    //this is used to save the BC-conditions...
    /*hmc_complex bc_tmp = (dir == TDIR) ? (hmc_complex) {
        1./16. * csw * TEMPORAL_RE, 1./16. * csw * TEMPORAL_IM
    } :
    (hmc_complex) {
        1./16. * csw * SPATIAL_RE, 1./16. * csw * SPATIAL_IM
    };*/
    
    //the matrix consits out of 4 3x3 blocks which are calculated now
    
    //E1 = 8 * F_01
    tmp = field_strength_tensor(field, idx_arg, 0, 1);
    E = multiply_matrix3x3_by_real(tmp, 8.);
    //B1 = 4 * (F_23 - F_32)
    tmp = field_strength_tensor(field, idx_arg, 2, 3);
    tmp1 = field_strength_tensor(field, idx_arg, 3, 2);
    B = subtract_matrix3x3(tmp, tmp1);
    B = multiply_matrix3x3_by_real(tmp, 4.);
    //EB1 = E1 - B1
    EB1 = subtract_matrix3x3(E, B);
    EB1 = multiply_matrix3x3_by_complex(EB1, factor);
    
    //E2 = 8 * F_02
    tmp = field_strength_tensor(field, idx_arg, 0, 2);
    E = multiply_matrix3x3_by_real(tmp, 8.);
    //B2 = 4 * (F_31 - F_13)
    tmp = field_strength_tensor(field, idx_arg, 3, 1);
    tmp1 = field_strength_tensor(field, idx_arg, 1, 3);
    B = subtract_matrix3x3(tmp, tmp1);
    B = multiply_matrix3x3_by_real(tmp, 4.);
    //EB2 = E2 - B2
    EB2 = subtract_matrix3x3(E, B);
    EB2 = multiply_matrix3x3_by_complex(EB2, factor);
    
    //E3 = 8 * F_03
    tmp = field_strength_tensor(field, idx_arg, 0, 3);
    E = multiply_matrix3x3_by_real(tmp, 8.);
    //B3 = 4 * (F_12 - F_21)
    tmp = field_strength_tensor(field, idx_arg, 1, 2);
    tmp1 = field_strength_tensor(field, idx_arg, 2, 1);
    B = subtract_matrix3x3(tmp, tmp1);
    B = multiply_matrix3x3_by_real(tmp, 4.);
    //EB3 = E3 - B3
    EB3 = subtract_matrix3x3(E, B);
    EB3 = multiply_matrix3x3_by_complex(EB3, factor);
    
    //upper-left block = 1 + EB3
    tmp = add_matrix3x3(identity_matrix3x3(), EB3);
    out = put_3x3block_matrix6x6_upperleft(out, tmp);
    
    //upper-right block = EB1 - i * EB2
    tmp = multiply_matrix3x3_by_complex(EB2, hmc_complex_i);
    out = put_3x3block_matrix6x6_upperright(out, subtract_matrix3x3(EB1, tmp));
    
    //lower-left block = EB1 + i * EB2
    tmp = multiply_matrix3x3_by_complex(EB2, hmc_complex_i);
    out = put_3x3block_matrix6x6_lowerleft(out, add_matrix3x3(EB1, tmp));
    
    //lower-right block = 1 - EB3
    tmp = subtract_matrix3x3(identity_matrix3x3(), EB3);
    out = put_3x3block_matrix6x6_lowerright(out, tmp);
    
    return out;
}

/*the lower-right block is given by c_sw * i/16 * sum{k=1}^3 (-)sigma_k (E_k + B_k)
 with sigma_k Pauli matrices, E_k = 8 * F_{0,k}, B_k = sum_{j,l=1}^3 4 * epsilon_{k,l,j} F_{l,j}
 */
Matrix6x6 clover_eoprec_unified_local_lower_right_block(__global Matrixsu3StorageType  const * const restrict field, const st_idx idx_arg, hmc_float csw)
{
    Matrix6x6 out = zero_matrix6x6();
    Matrix3x3 EB1, EB2, EB3, E, B, tmp, tmp1;
    //no bc_tmp needed
    hmc_complex factor = {0., 1./16. * csw};
    //this is used to save the BC-conditions...
    /*hmc_complex bc_tmp = (dir == TDIR) ? (hmc_complex) {
        1./16. * csw * TEMPORAL_RE, 1./16. * csw * TEMPORAL_IM
    } :
    (hmc_complex) {
        1./16. * csw * SPATIAL_RE, 1./16. * csw * SPATIAL_IM
    };*/
    
    //the matrix consits out of 4 3x3 blocks which are calculated now(cf. equation 4.1 in openQCD documentation)

    //E1 = 8 * F_01
    tmp = field_strength_tensor(field, idx_arg, 0, 1);
    E = multiply_matrix3x3_by_real(tmp, 8.);
    //B1 = 4 * (F_23 - F_32)
    tmp = field_strength_tensor(field, idx_arg, 2, 3);
    tmp1 = field_strength_tensor(field, idx_arg, 3, 2);
    B = subtract_matrix3x3(tmp, tmp1);
    B = multiply_matrix3x3_by_real(tmp, 4.);
    //EB1 = E1 + B1
    EB1 = add_matrix3x3(E, B);
    EB1 = multiply_matrix3x3_by_complex(EB1, factor);
    
    //E2 = 8 * F_02
    tmp = field_strength_tensor(field, idx_arg, 0, 2);
    E = multiply_matrix3x3_by_real(tmp, 8.);
    //B2 = 4 * (F_31 - F_13)
    tmp = field_strength_tensor(field, idx_arg, 3, 1);
    tmp1 = field_strength_tensor(field, idx_arg, 1, 3);
    B = subtract_matrix3x3(tmp, tmp1);
    B = multiply_matrix3x3_by_real(tmp, 4.);
    //EB2 = E2 + B2
    EB2 = add_matrix3x3(E, B);
    EB2 = multiply_matrix3x3_by_complex(EB2, factor);
    
    //E3 = 8 * F_03
    tmp = field_strength_tensor(field, idx_arg, 0, 3);
    E = multiply_matrix3x3_by_real(tmp, 8.);
    //B3 = 4 * (F_12 - F_21)
    tmp = field_strength_tensor(field, idx_arg, 1, 2);
    tmp1 = field_strength_tensor(field, idx_arg, 2, 1);
    B = subtract_matrix3x3(tmp, tmp1);
    B = multiply_matrix3x3_by_real(tmp, 4.);
    //EB3 = E3 + B3
    EB3 = add_matrix3x3(E, B);
    EB3 = multiply_matrix3x3_by_complex(EB3, factor);
    
    //upper-left block = 1 + EB3
    tmp = add_matrix3x3(identity_matrix3x3(), EB3);
    out = put_3x3block_matrix6x6_upperleft(out, tmp);
    
    //upper-right block = EB1 - i * EB2
    tmp = multiply_matrix3x3_by_complex(EB2, hmc_complex_i);
    out = put_3x3block_matrix6x6_upperright(out, subtract_matrix3x3(EB1, tmp));
    
    //lower-left block = EB1 + i * EB2
    tmp = multiply_matrix3x3_by_complex(EB2, hmc_complex_i);
    out = put_3x3block_matrix6x6_lowerleft(out, add_matrix3x3(EB1, tmp));
    
    //lower-right block = 1 - EB3
    tmp = subtract_matrix3x3(identity_matrix3x3(), EB3);
    out = put_3x3block_matrix6x6_lowerright(out, tmp);

    return out;
}
