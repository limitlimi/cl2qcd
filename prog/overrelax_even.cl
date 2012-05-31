/** @file
 * Device code for the heatbath overrelaxation
 */

__kernel void overrelax_even(__global Matrixsu3StorageType * const restrict gaugefield, const int mu, __global rngStateStorageType * const restrict rngStates)
{
	prng_state rnd;
	prng_loadState(&rnd, rngStates);

	PARALLEL_FOR(id, VOLSPACE * NTIME / 2) {
		st_index pos = get_even_site(id);
		perform_overrelaxing(gaugefield, mu, &rnd, pos.space, pos.time);
	}

	prng_storeState(rngStates, &rnd);
}
