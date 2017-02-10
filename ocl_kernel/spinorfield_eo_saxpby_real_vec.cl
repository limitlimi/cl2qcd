/*
 * Copyright 2017 Christopher Czaban
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

// Description of variables of saxpby:
//  - x: The first input spinor field (one spinor per site => vector of VOL4D)
//  - y: The second input spinor field (one spinor per site => vector of VOL4D)
//  - alpha: Vector of constants
//  - beta: Vector of constants
//  - index_alpha: Constant alpha to be used
//  - index_beta: Constant beta to be used
//  - out: The output spinor field: alpha*x+beta*y (site by site)

__kernel void saxpby_eoprec_real_vec(__global const spinorStorageType * const x, __global const spinorStorageType * const y, __global const hmc_float * const alpha, __global hmc_float * beta, const int index_alpha, const int index_beta, __global spinorStorageType * const out)
{
	const int id = get_global_id(0);
	const int global_size = get_global_size(0);
	
	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		const spinor x_tmp = getSpinor_eo(x, id_mem);
		const spinor x_tmp_tmp = real_multiply_spinor(x_tmp, alpha[index_alpha]);
		const spinor y_tmp = getSpinor_eo(y, id_mem);
		const spinor y_tmp_tmp = real_multiply_spinor(y_tmp, beta[index_beta]);
		const result_tmp = spinor_acc(x_tmp_tmp, y_tmp_tmp);
		putSpinor_eo(out, id_mem, result_tmp);
	}
}