/** @todo rewrite gaussianComplexVector to handle structs??*/
__kernel void generate_gaussian_gaugemomenta(__global ae * const restrict out, __global rngStateStorageType * const restrict rngStates)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	hmc_ocl_ran rnd = loadRngState(rngStates);

	hmc_complex tmp;
#ifdef _SAME_RND_NUMBERS_
       if(id>0) return;
       global_size = 1;
#endif
	for(int id_tmp = id; id_tmp < GAUGEMOMENTASIZE; id_tmp += global_size) {
		//CP: THERE ARE 8 ELEMENTS IN AE
		tmp = gaussianNormalPair(&rnd);
		out[id_tmp].e0 = tmp.re;
		out[id_tmp].e1 =  tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out[id_tmp].e2 =  tmp.re;
		out[id_tmp].e3 =  tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out[id_tmp].e4 =  tmp.re;
		out[id_tmp].e5 =  tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out[id_tmp].e6 =  tmp.re;
		out[id_tmp].e7 =  tmp.im;
	}

	storeRngState(rngStates, rnd);
}
