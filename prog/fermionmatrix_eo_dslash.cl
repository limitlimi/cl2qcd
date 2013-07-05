//Here, two cases are possible:
//evenodd = ODD or EVEN
//	ODD corresponds to the D_oe case: Dslash acts on even indices (the "x+mu" in the formulae) and the function
//	saves the outcoming spinor with an odd index.
//	EVEN is then D_eo.

void dslash_eo_for_site(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in, st_idx const pos)
{
	spinor out_tmp = set_spinor_zero();
	spinor out_tmp2;

	//calc dslash (this includes mutliplication with kappa)

	out_tmp2 = dslash_eoprec_unified_local(in, field, pos, TDIR, kappa_in);
	out_tmp = spinor_dim(out_tmp, out_tmp2);
	out_tmp2 = dslash_eoprec_unified_local(in, field, pos, XDIR, kappa_in);
	out_tmp = spinor_dim(out_tmp, out_tmp2);
	out_tmp2 = dslash_eoprec_unified_local(in, field, pos, YDIR, kappa_in);
	out_tmp = spinor_dim(out_tmp, out_tmp2);
	out_tmp2 = dslash_eoprec_unified_local(in, field, pos, ZDIR, kappa_in);
	out_tmp = spinor_dim(out_tmp, out_tmp2);

	putSpinor_eo(out, get_eo_site_idx_from_st_idx(pos), out_tmp);
}

__kernel void dslash_eo(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in)
{
	PARALLEL_FOR(id_local, EOPREC_SPINORFIELDSIZE_LOCAL) {
		st_idx pos = (evenodd == ODD) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);
		dslash_eo_for_site(in, out, field, evenodd, kappa_in, pos);
	}
}

#define REQD_HALO_WIDTH 1

__kernel void dslash_eo_inner(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in)
{
	PARALLEL_FOR(id_local, EOPREC_SPINORFIELDSIZE_LOCAL) {
		st_idx pos = (evenodd == ODD) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);
		if(pos.time >= REQD_HALO_WIDTH && pos.time < (NTIME_LOCAL - REQD_HALO_WIDTH)) {
			dslash_eo_for_site(in, out, field, evenodd, kappa_in, pos);
		}
	}
}

__kernel void dslash_eo_boundary(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in)
{
	// TODO only loop over boundary sites instead of masking everybody else out
	PARALLEL_FOR(id_local, EOPREC_SPINORFIELDSIZE_LOCAL) {
		st_idx pos = (evenodd == ODD) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);
		if(pos.time < REQD_HALO_WIDTH || pos.time >= (NTIME_LOCAL - REQD_HALO_WIDTH)) {
			dslash_eo_for_site(in, out, field, evenodd, kappa_in, pos);
		}
	}
}
