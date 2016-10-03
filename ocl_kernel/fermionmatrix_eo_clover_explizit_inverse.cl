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


void clover_eo_inverse_explicit_upper_left_for_site(__global Matrix6x6StorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, hmc_float kappa_in, hmc_float csw, st_idx const pos)
{
    Matrix6x6 out_tmp, tmp;
    
    tmp = clover_eoprec_unified_local_upper_left_block(field, pos, kappa_in, csw);
    out_tmp = inverse_6x6_via_Householder_triangularization(tmp);    

    put6x6(out, get_site_idx(pos), out_tmp);
}

void clover_eo_inverse_explicit_lower_right_for_site(__global Matrix6x6StorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, hmc_float kappa_in, hmc_float csw, st_idx const pos)
{
    Matrix6x6 out_tmp, tmp;
    
    tmp = clover_eoprec_unified_local_lower_right_block(field, pos, kappa_in, csw);
    out_tmp = inverse_6x6_via_Householder_triangularization(tmp);
    
    put6x6(out, get_site_idx(pos), out_tmp);
}

__kernel void clover_eo_inverse_explicit_upper_left(__global Matrix6x6StorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, hmc_float kappa_in, hmc_float csw)
{
    PARALLEL_FOR(id_local, SPINORFIELDSIZE_LOCAL) {
        st_idx pos = get_st_idx_from_site_idx(id_local);
        clover_eo_inverse_explicit_upper_left_for_site(out, field, kappa_in, csw, pos);
    }
}

__kernel void clover_eo_inverse_explicit_lower_right(__global Matrix6x6StorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, hmc_float kappa_in, hmc_float csw)
{
    PARALLEL_FOR(id_local, SPINORFIELDSIZE_LOCAL) {
        st_idx pos = get_st_idx_from_site_idx(id_local);
        clover_eo_inverse_explicit_lower_right_for_site(out, field, kappa_in, csw, pos);
    }
}
