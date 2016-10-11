/*
 * Copyright 2016, Alessandro Sciarra, Tim Breitenfelder
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

/**
 @file fermion-observables
*/

// Correlator is given by:
// C(t)= - 64 * (-1)^t sum_{vec{x}} sum_{c,d} |[D^(-1)_f (vec{x}|0)]_{c,d}|^2
// We drop the factor of -64*(-1)^t since it can be neglected for the final mass extraction 

__kernel void correlator_staggered_ps(__global hmc_float * const restrict out, __global const staggeredStorageType * const restrict phi) //<- Put one additional field
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	
	//suppose that there are NTIME threads (one for each entry of the correlator) divided by two because of using even odd precondition	
	
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size ) 
	{
		hmc_float summedSquarenorms = 0.;
		uint3 coord;
		
		for(coord.x = 0; coord.x < NSPACE; coord.x += 1 ) 
		{
			for(coord.y = 0; coord.y < NSPACE; coord.y += 1 ) 
			{
				for(coord.z = 0; coord.z < NSPACE; coord.z += 1 ) 
				{
					//
					// Here understand if field is even or odd and get temporalField accordingly from phiEven or from phiOdd
					// ATTENTION: use get_su3vec_from_field_eo function (see spinorfield_staggered_eo.cl file)
					//
					
					int nspace = get_nspace(coord);
					su3vec temporalField = phi[get_n_eoprec(nspace, id_tmp)];
					
					// taking squarenorm:
					summedSquarenorms += su3vec_squarenorm(temporalField);
				}
			}
		}
		
		//consider taking into account an normalization factor like in wilson case:
		//hmc_float fac = NSPACE * NSPACE * NSPACE;
		//out[NTIME_OFFSET + id_tmp] += 2. * KAPPA * 2. * KAPPA * correlator / fac;
		
		out[NTIME_OFFSET + id_tmp] += -64 * summedSquarenorms; //<-neglect factor -64 but put division by NSPACE^3
		
	}

	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}
