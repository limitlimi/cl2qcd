// Description of variables of this kernel:
//  - out: The output staggered field: gaussian complex su3vec (site by site)
//  - rngStates: The state of the generator (rngStateStorageType is useful to produce
//               random number --> see random.cl)

__kernel void set_gaussian_spinorfield_stagg(__global su3vec * const restrict out, __global rngStateStorageType * const restrict rngStates)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
		
	// complex to store the gaussian numbers drawn
	hmc_complex tmp;
	// sigma has to be 0.5 here (for the explanation of why, see the documentation
	// of the function set_gaussian_spinorfield_device in spinors_staggered.hpp)
	hmc_float sigma = 0.5;
	
	su3vec out_tmp;

#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif

	prng_state rnd;
	prng_loadState(&rnd, rngStates);

	for(int id_local = id; id_local < SPINORFIELDSIZE_LOCAL; id_local += global_size) {
		/** @todo this must be done more efficient */
		st_index pos = (id_local < SPINORFIELDSIZE_LOCAL / 2) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local - (SPINORFIELDSIZE_LOCAL / 2));

		//There are NC=3 complex elements in the su3vec
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e0.re = tmp.re;
		out_tmp.e0.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e1.re = tmp.re;
		out_tmp.e1.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e2.re = tmp.re;
		out_tmp.e2.im = tmp.im;

		// multiply by sigma because gaussianNormalPair generates a couple
		// of real gaussian number distributed with variance 1 instead of 0.5 (see the 
		// documentation of the function set_gaussian_spinorfield_device in spinors_staggered.hpp)
		// Remark that here the variable sigma is the VARIANCE even if it is called sigma!!!
		
		out_tmp = su3vec_times_real(out_tmp, sqrt(sigma));

		put_su3vec_to_field(out_tmp, out, pos.space, pos.time);
	}

	prng_storeState(rngStates, &rnd);
}
