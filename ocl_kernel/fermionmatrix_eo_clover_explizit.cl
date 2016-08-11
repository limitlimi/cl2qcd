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
    // U_nu(x_nu)^dagger
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
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
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    //////////////////////
    // U_mu(x-mu-nu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir1); // x-mu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3(U, tmp);
    //////////////////////
    // U_nu(x-nu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir2);
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3(U, tmp);
    /////////////////////
    out = add_matrix3x3(out, matrix_su3to3x3(tmp));
    
    
    ///////////////////////
    // 4.term = U_nu(x-nu)^dagger * U_mu(x-nu) * U_nu(x-nu+mu) * U_mu(x)^dagger
    su3vec psi4 = inout;
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
    out = add_matrix3x3(out, matrix_su3to3x3(tmp));
    
    
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
    out = add_matrix3x3(out, matrix_su3to3x3(tmp));
    
    
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
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3_dagger(tmp, U);
    /////////////////////////
    // U_nu(x-nu-mu)
    idx_neigh1 = get_lower_neighbor_from_st_idx(idx_arg, dir1); // x-mu
    idx_neigh = get_neighbor_from_st_idx(idx_neigh1, dir2); // (x-mu)-nu
    U = getSU3(field, get_link_idx(dir2, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    ////////////////////////
    // U_mu(x-mu)
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir1);
    U = getSU3(field, get_link_idx(dir1, idx_neigh));
    tmp = multiply_matrixsu3(tmp, U);
    /////////////////////
    out = add_matrix3x3(out, matrix_su3to3x3(tmp));

    
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
    out = add_matrix3x3(out, matrix_su3to3x3(tmp));
    
    return out;
}