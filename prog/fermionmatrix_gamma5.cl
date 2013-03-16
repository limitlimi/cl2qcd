__kernel void gamma5(__global spinor * const restrict inout)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;

	for(int id_local = id; id_local < SPINORFIELDSIZE_LOCAL; id_local += global_size) {
		/** @todo this must be done more efficient */
		st_index pos = (id_local % 2 == 0) ? get_even_st_idx_local(id_local / 2) : get_odd_st_idx_local(id_local / 2);

		out_tmp = get_spinor_from_field(inout, pos.space, pos.time);
		out_tmp = gamma5_local(out_tmp);
		put_spinor_to_field(out_tmp, inout, pos.space, pos.time);
	}
}
