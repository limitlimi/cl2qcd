__kernel void set_zero_spinorfield_eoprec( __global hmc_complex * const restrict x )
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {
		putSpinorSOA_eo(x, id_tmp, set_spinor_zero());
	}
	return;
}
