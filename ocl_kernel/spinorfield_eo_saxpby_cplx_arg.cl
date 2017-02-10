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
//  - x: The first input spinor field (an spinor per each site => vector of VOL4D
//       components that are spinor varibles)
//  - y: The second input spinor field (an spinor per each site => vector of VOL4D
//       components that are spinor varibles)
//  - alpha: The complex number by which x has to be multiplied
//  - beta:  The complex number by which y has to be multiplied
//  - out: The output spinor field: alpha*x+beta*y (site by site)



__kernel void saxpby_eoprec_cplx_arg(__global const spinorStorageType * const x, __global const spinorStorageType * const y, const hmc_float alpha_re, const hmc_float alpha_im, const hmc_float beta_re, const hmc_float beta_im, __global spinorStorageType * const out)
{

	const int id = get_global_id(0);
	const int global_size = get_global_size(0);

	const hmc_complex alpha = (hmc_complex) {alpha_re, alpha_im};
	const hmc_complex beta = (hmc_complex) { beta_re,  beta_im};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		const spinor x_tmp = getSpinor_eo(x,id_mem);
		const spinor x_tmp_tmp = spinor_times_complex(x_tmp, alpha);
		const spinor y_tmp = getSpinor_eo(y,id_mem);
		const spinor y_tmp_tmp = spinor_times_complex(y_tmp, beta);
		const spinor x_y_sum = spinor_acc(y_tmp_tmp, x_tmp_tmp);
		putSpinor_eo(out, id_mem, x_y_sum);
	}	
}
