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



//alpha*x +y
//CC: defined with a plus!!!
__kernel void saxpy_eoprec_real_arg(__global const spinorStorageType * const x, __global const spinorStorageType * const y, const hmc_float alpha, __global spinorStorageType * const out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		spinor x_tmp = getSpinor_eo(x, id_mem);
		spinor y_tmp = getSpinor_eo(y, id_mem);
		x_tmp = real_multiply_spinor(x_tmp, alpha);
		x_tmp = spinor_acc(y_tmp, x_tmp);
		putSpinor_eo(out, id_mem, x_tmp);
	}
}
