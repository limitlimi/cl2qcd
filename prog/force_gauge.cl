__kernel void gauge_force(__global const Matrixsu3StorageType * const restrict field, __global aeStorageType * const restrict out)
{
	//Gauge force is factor*Im(i Tr(T_i U V))
	//   with T_i being the SU3-Generator in i-th direction and V the staplematrix
	//   and the factor being - beta / NC (for standard Wilson-action) and -c0 * beta / NC (for tlSym)
	PARALLEL_FOR(id_tmp, VOL4D * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		const int pos_tmp = id_tmp % VOL4D;
		const int dir     = id_tmp / VOL4D;

		const st_index pos = (pos_tmp % 2 == 0) ? get_even_site(pos_tmp / 2) : get_odd_site(pos_tmp / 2);

		Matrix3x3 V = calc_staple(field, pos.space, pos.time, dir);
		Matrixsu3 U = get_matrixsu3(field, pos.space, pos.time, dir);
		V = multiply_matrix3x3 (matrix_su3to3x3(U), V);
		ae out_tmp = tr_lambda_u(V);

		hmc_float factor = -BETA / 3.;
#ifdef _USE_RECT_
		factor = factor * C0;
#endif
		int global_link_pos = get_global_link_pos(dir, pos.space, pos.time);
		update_gaugemomentum(out_tmp, factor , global_link_pos, out);
	}
}
