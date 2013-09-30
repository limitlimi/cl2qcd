/** @file
* Staggered (local) D_KS operator that is called from within kernels \b with EVEN-ODD preconditioning
* \internal (for example, see fermionmatrix_staggered.cl)
*/

/**
  * This kernel is nothing but the local D_KS working on a particular link in a specific direction, 
  * taking into account even-odd preconditioning. Indeed, here there are no conceptual aspects to
  * be implemented to take into account eo-prec. The difference between this kernel and D_KS_local
  * is that here the coordinates of the neighbors have to be transformed into an eoprec index and then
  * the functions get_su3vec_from_field_eo must be used.
  * \internal (see spinorfield_staggered_eo.cl for these 2 functions) \endinternal  
  * The expression of D_KS for a specific (couple of) site(s) and a specific direction is
  *  \f[
     (D_{KS})_{n,m,\mu}=\frac{1}{2} \eta_\mu(n)\Bigl[U_\mu(n)\,\delta_{n+\hat\mu,m} - U^\dag_\mu(n-\hat\mu)\,\delta_{n-\hat\mu,m}\Bigr]
     \f]
  * This function returns the value of the field (D_KS*in) at the site idx_arg [so in the
  * function only values of the field "in" in the nextneighbour of idx_arg will be needed]:
  *  \f[
     \bigl[(D_{KS})_\mu\cdot \text{\texttt{in}}\bigr]_n=\frac{1}{2}\eta_\mu(n) \Bigl[U_\mu(n) \cdot\text{\texttt{in}}_{n+\hat\mu} - U^\dag_\mu(n-\hat\mu)\cdot\text{\texttt{in}}_{n-\hat\mu}\Bigr]
    \f]
  *
  * The variables passed to the kernel are:
  *  @param in The input staggered field on half lattice (either on even or odd sites)
  *  @param field The links configuration
  *  @param idx_arg The superindex of the site where the output field is returned
  *  @param dir The direction in which D_KS works
  * @return The field in the idx_arg site is returned. Remark that if the "in" staggered
  *         field is on even sites then idx_arg will be an odd site and vice-versa.
  * 
  * @note The staggered phases are included in this function with the help of the function
  *       get_modified_stagg_phase @internal(see operations_staggered.cl)@endinternal.
  * \par
  * @note Since we chose to impose boundary conditions modifying staggered phases at the end of the
  *       lattice in each direction, then we make each staggered phase appear EXACTELY
  *       next to the link, and if the link is dagger, we take the complex coniugate
  *       of the staggered phase (that would be complex in general due to the modification).
  *       Next to the link means that it must be calculated in the same site where the
  *       link is considered.
  * 
  * @todo If a chemical potential is introduced, this kernel has to be modified!
  */
su3vec D_KS_eo_local(__global const staggeredStorageType * const restrict in, __global const Matrixsu3StorageType * const restrict field, const st_idx idx_arg, const dir_idx dir)
{
	//this is used to save the idx of the neighbors
	st_idx idx_neigh;
	//this is used to transform the idx of the neighbors, that is a superindex (between 0 to VOL4D-1),
	//to a superindex of type even-odd (between 0 and VOL4D/2-1)
	site_idx nn_eo;
	//this are used for the calculation
	su3vec out_tmp, plus, chi;
	Matrixsu3 U;
	//this is used to take into account the staggered phase and the BC-conditions...
	hmc_complex eta_mod;
	hmc_float eta;
	
	out_tmp = set_su3vec_zero();
	
	//go through the different directions
	///////////////////////////////////
	// mu = +dir
	///////////////////////////////////
	idx_neigh = get_neighbor_from_st_idx(idx_arg, dir);
	//transform normal indices to eoprec index
	nn_eo = get_eo_site_idx_from_st_idx(idx_neigh);
	plus = get_su3vec_from_field_eo(in, nn_eo);
	U = getSU3(field, get_link_idx(dir, idx_arg));
	//Thanks to the variables passed to this kernel I can write all directions in few lines:
	//chi=U*plus
	chi=su3matrix_times_su3vec(U,plus);
	coord_spatial coord = get_coord_spatial(idx_arg.space);
	if(coord.x == (NSPACE-1) || coord.y == (NSPACE-1) || coord.z == (NSPACE-1) ||
	   idx_arg.time == (NTIME_GLOBAL-1)){
		eta_mod = get_modified_stagg_phase(idx_arg.space, idx_arg.time, dir);
		eta_mod.re *= 0.5; //to take into account the factor at the beginning of D_KS
		chi = su3vec_times_complex(chi, eta_mod);
	}else{
		if(dir != XDIR)
		  eta = 0.5 * get_staggered_phase(idx_arg.space, dir);
		else
		  eta = 0.5;
		chi = su3vec_times_real(chi, eta);
	}
	out_tmp=su3vec_acc(out_tmp,chi);
	
	///////////////////////////////////
	// mu = -dir
	///////////////////////////////////
	idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir);
	//transform normal indices to eoprec index
	nn_eo = get_eo_site_idx_from_st_idx(idx_neigh);
	plus = get_su3vec_from_field_eo(in, nn_eo);
	U = getSU3(field, get_link_idx(dir, idx_neigh));
	//Thanks to the variables passed to this kernel I can write all directions in few lines:
	//chi=((0.5*eta)^conjugated * U^dagger) * plus
	chi=su3matrix_dagger_times_su3vec(U,plus);
	coord = get_coord_spatial(idx_neigh.space);
	if(coord.x == (NSPACE-1) || coord.y == (NSPACE-1) || coord.z == (NSPACE-1) ||
	   idx_neigh.time == (NTIME_GLOBAL-1)){
		eta_mod = get_modified_stagg_phase(idx_neigh.space, idx_neigh.time, dir);
		eta_mod.re *= 0.5; //to take into account the factor at the beginning of D_KS
		chi = su3vec_times_complex_conj(chi, eta_mod);
	}else{
		if(dir != XDIR)
		  eta = 0.5 * get_staggered_phase(idx_arg.space, dir);
		else
		  eta = 0.5;
		chi = su3vec_times_real(chi, eta);
	}
	out_tmp=su3vec_dim(out_tmp,chi);

	return out_tmp;
}

