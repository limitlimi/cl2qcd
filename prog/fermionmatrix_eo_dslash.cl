//Here, two cases are possible:
//evenodd = ODD or EVEN
//	ODD corresponds to the D_oe case: Dslash acts on even indices (the "x+mu" in the formulae) and the function
//	saves the outcoming spinor with an odd index.
//	EVEN is then D_eo.
#ifdef _USEGPU_
__attribute__((reqd_work_group_size(128, 1, 1)))
#endif
__kernel void dslash_eoprec(__global const hmc_float * const restrict in, __global hmc_float * const restrict out, __global const hmc_float * const restrict field, const int evenodd)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {
		st_index pos = (evenodd == ODD) ? get_even_site(id_tmp) : get_odd_site(id_tmp);

		spinor out_tmp = set_spinor_zero();
		spinor out_tmp2;

		//calc dslash (this includes mutliplication with kappa)

		out_tmp2 = dslash_eoprec_local_0(in, field, pos.space, pos.time);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_local_1(in, field, pos.space, pos.time);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_local_2(in, field, pos.space, pos.time);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_local_3(in, field, pos.space, pos.time);
		out_tmp = spinor_dim(out_tmp, out_tmp2);

		putSpinorSOA_eo(out, id_tmp, out_tmp);
	}
}

__kernel void convertSpinorfieldToSOA_eo(__global hmc_float * const restrict out, __global const spinor * const restrict in)
{
	for(uint i = get_global_id(0); i < EOPREC_SPINORFIELDSIZE; i += get_global_size(0)) {
		putSpinorSOA_eo(out, i, in[i]);
	}
}
__kernel void convertSpinorfieldFromSOA_eo(__global spinor * const restrict out, __global const hmc_float * const restrict in)
{
	for(uint i = get_global_id(0); i < EOPREC_SPINORFIELDSIZE; i += get_global_size(0)) {
		out[i] = getSpinorSOA_eo(in, i);
	}
}

__kernel void convertGaugefieldToSOA(__global hmc_float * const restrict out, __global const Matrixsu3 * const restrict in)
{
	for(uint i = get_global_id(0); i < NDIM * VOL4D; i += get_global_size(0)) {
		putSU3SOA(out, i, in[i]);
	}
}
__kernel void convertGaugefieldFromSOA(__global Matrixsu3 * const restrict out, __global const hmc_float * const restrict in)
{
	for(uint i = get_global_id(0); i < NDIM * VOL4D; i += get_global_size(0)) {
		out[i] = getSU3SOA(in, i);
	}
}
