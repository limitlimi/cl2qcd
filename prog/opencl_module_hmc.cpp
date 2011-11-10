#include "opencl_module_hmc.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Hmc::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Fermions::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters()->get_beta() << " -DGAUGEMOMENTASIZE=" << get_parameters()->get_gaugemomentasize();
	return;
}


void Opencl_Module_Hmc::fill_buffers()
{

	Opencl_Module_Fermions::fill_buffers();
	///@todo CP: some of the above buffers are not used and can be deleted again!! especially in the eoprec-case
	
	int spinorfield_size = get_parameters()->get_sf_buf_size();
	int eoprec_spinorfield_size = get_parameters()->get_eo_sf_buf_size();
	int gaugemomentum_size = get_parameters()->get_gm_buf_size();
	int gaugefield_size = get_parameters()->get_gf_buf_size();
	int float_size = sizeof(hmc_float);
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;

	//init mem-objects

	logger.trace() << "Create buffer for HMC...";
	clmem_force = create_rw_buffer(gaugemomentum_size);
	if(get_parameters()->get_use_eo() == true){
		clmem_phi_inv_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_phi_eoprec = create_rw_buffer(eoprec_spinorfield_size);
	} 
	else{
		clmem_phi = create_rw_buffer(spinorfield_size);
	}
	//this is always used for the force-calculation
	clmem_phi_inv = create_rw_buffer(spinorfield_size);
	clmem_new_u = create_rw_buffer(gaugefield_size);
	clmem_p = create_rw_buffer(gaugemomentum_size);
	clmem_new_p = create_rw_buffer(gaugemomentum_size);
	clmem_energy_init = create_rw_buffer(float_size);
	clmem_p2 = create_rw_buffer(float_size);
	clmem_new_p2 = create_rw_buffer(float_size);
	clmem_s_fermion = create_rw_buffer(float_size);

	return;
}

void Opencl_Module_Hmc::fill_kernels()
{
	Opencl_Module_Fermions::fill_kernels();

	basic_hmc_code = basic_fermion_code << "types_hmc.h";

	//init kernels for HMC
	if(get_parameters()->get_use_eo() == true){
		generate_gaussian_spinorfield_eoprec = createKernel("generate_gaussian_spinorfield_eoprec") << basic_hmc_code << "random.cl" << "spinorfield_eo_gaussian.cl";
		fermion_force_eoprec = createKernel("fermion_force_eoprec") << basic_hmc_code << "operations_gaugemomentum.cl" << "fermionmatrix.cl" << "force_fermion_eo.cl";
	}
	else{
		generate_gaussian_spinorfield = createKernel("generate_gaussian_spinorfield") << basic_hmc_code << "random.cl" << "spinorfield_gaussian.cl";
	}
	fermion_force = createKernel("fermion_force") << basic_hmc_code << "operations_gaugemomentum.cl" << "fermionmatrix.cl" << "force_fermion.cl";
	set_zero_gaugemomentum = createKernel("set_zero_gaugemomentum") << basic_hmc_code << "gaugemomentum_zero.cl";
	generate_gaussian_gaugemomenta = createKernel("generate_gaussian_gaugemomenta") << basic_hmc_code << "random.cl" << "gaugemomentum_gaussian.cl";
	md_update_gaugefield = createKernel("md_update_gaugefield") << basic_hmc_code << "md_update_gaugefield.cl";
	md_update_gaugemomenta = createKernel("md_update_gaugemomenta") << basic_hmc_code << "operations_gaugemomentum.cl" << "md_update_gaugemomenta.cl";
	gauge_force = createKernel("gauge_force") << basic_hmc_code << "operations_gaugemomentum.cl" << "force_gauge.cl";
	
	if(get_parameters()->get_use_smearing() == true) {
		stout_smear_fermion_force = createKernel("stout_smear_fermion_force") << basic_hmc_code << "force_fermion_stout_smear.cl";
	}
	gaugemomentum_squarenorm = createKernel("gaugemomentum_squarenorm") << basic_hmc_code << "operations_gaugemomentum.cl" << "gaugemomentum_squarenorm.cl";

	return;
}

void Opencl_Module_Hmc::clear_kernels()
{
	Opencl_Module_Fermions::clear_kernels();

	cl_uint clerr = CL_SUCCESS;

	logger.debug() << "release HMC-kernels.." ;
	if(get_parameters()->get_use_eo() == true){
		clerr = clReleaseKernel(generate_gaussian_spinorfield_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(fermion_force_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	else{
		clerr = clReleaseKernel(generate_gaussian_spinorfield);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(fermion_force);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	clerr = clReleaseKernel(generate_gaussian_gaugemomenta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(md_update_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(md_update_gaugemomenta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(gauge_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(set_zero_gaugemomentum);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(get_parameters()->get_use_smearing() == true) {
		clerr = clReleaseKernel(stout_smear_fermion_force);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}

	return;
}

void Opencl_Module_Hmc::clear_buffers()
{
	Opencl_Module_Fermions::clear_buffers();

	cl_uint clerr = CL_SUCCESS;

	logger.debug() << "release HMC-variables.." ;
	clerr = clReleaseMemObject(clmem_energy_init);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_p2);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_new_p2);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_new_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	if(get_parameters()->get_use_eo() == true){
		clerr = clReleaseMemObject(clmem_phi_inv_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_phi_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
	else{
		clerr = clReleaseMemObject(clmem_phi_inv);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_phi);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
	return;
}

void Opencl_Module_Hmc::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
	Opencl_Module_Fermions::get_work_sizes(kernel, dev_type, ls, gs, num_groups);

	return;
}

////////////////////////////////////////////////////
//Access to members

cl_mem Opencl_Module_Hmc::get_clmem_p()
{
	return clmem_p;
}

cl_mem Opencl_Module_Hmc::get_clmem_new_p()
{
	return clmem_new_p;
}

cl_mem Opencl_Module_Hmc::get_clmem_new_u()
{
	return clmem_new_u;
}

cl_mem Opencl_Module_Hmc::get_clmem_phi()
{
	return clmem_phi;
}

cl_mem Opencl_Module_Hmc::get_clmem_phi_eoprec()
{
	return clmem_phi_eoprec;
}

#ifdef _PROFILING_
usetimer* Opencl_Module_Hmc::get_timer(char * in)
{
	usetimer *noop = NULL;
	noop = Opencl_Module_Fermions::get_timer(in);
	if(noop != NULL) return noop;

	if (strcmp(in, "generate_gaussian_spinorfield") == 0) {
		return &this->timer_generate_gaussian_spinorfield;
	}
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0) {
		return &this->timer_generate_gaussian_gaugemomenta;
	}
	if (strcmp(in, "md_update_gaugefield") == 0) {
		return &this->timer_md_update_gaugefield;
	}
	if (strcmp(in, "md_update_gaugemomenta") == 0) {
		return &this->timer_md_update_gaugemomenta;
	}
	if (strcmp(in, "gauge_force") == 0) {
		return &this->timer_gauge_force;
	}
	if (strcmp(in, "fermion_force") == 0) {
		return &this->timer_fermion_force;
	}
	if (strcmp(in, "fermion_force_eoprec") == 0) {
		return &this->timer_fermion_force_eoprec;
	}
	if (strcmp(in, "set_zero_gaugemomentum") == 0) {
		return &this->timer_set_zero_gaugemomentum;
	}
	if (strcmp(in, "gaugemomentum_squarenorm") == 0) {
		return &this->timer_gaugemomentum_squarenorm;
	}
	if (strcmp(in, "stout_smear_fermion_force") == 0) {
		return &this->timer_stout_smear_fermion_force;
	}
	//if the kernelname has not matched, return NULL
	else {
		return NULL;
	}
}

int Opencl_Module_Hmc::get_read_write_size(char * in, inputparameters * parameters)
{
	int result = Opencl_Module_Fermions::get_read_write_size(in, parameters);
	if (result != 0) return result;
	//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	int R = (*parameters).get_mat_size();
	int S;
	if((*parameters).get_use_eo() == 1)
		S = get_parameters()->get_eoprec_spinorfieldsize();
	else
		S = get_parameters()->get_spinorfieldsize();
	//this is the same as in the function above
	if (strcmp(in, "generate_gaussian_spinorfield") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "generate_gaussian_spinorfield_eoprec") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "md_update_gaugefield") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "md_update_gaugemomenta") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "gauge_force") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "fermion_force") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "fermion_force_eoprec") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "set_zero_gaugemomentum;") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "gaugemomentum_squarenorm") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "stout_smear_fermion_force") == 0) {
		return 10000000000000000000;
	}
	return 0;
}

int Opencl_Module_Hmc::get_flop_size(const char * in, inputparameters * parameters)
{
	int result = Opencl_Module_Fermions::get_flop_size(in, parameters);
	if (result != 0) return result;
	int S;
	if((*parameters).get_use_eo() == 1)
		S = get_parameters()->get_eoprec_spinorfieldsize();
	else
		S = get_parameters()->get_spinorfieldsize();
	//this is the same as in the function above
	if (strcmp(in, "generate_gaussian_spinorfield") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "generate_gaussian_spinorfield_eoprec") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "md_update_gaugefield") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "md_update_gaugemomenta") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "gauge_force") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "fermion_force") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "fermion_force_eoprec") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "set_zero_gaugemomentum;") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "gaugemomentum_squarenorm") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "stout_smear_fermion_force") == 0) {
		return 10000000000000000000;
	}
	return 0;
}

void Opencl_Module_Hmc::print_profiling(std::string filename, int number)
{
	Opencl_Module_Fermions::print_profiling(filename, number);
	char * kernelName;
	kernelName = "generate_gaussian_spinorfield";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "generate_gaussian_spinorfield_eoprec";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "generate_gaussian_gaugemomenta";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "md_update_gaugefield";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "md_update_gaugemomenta";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "gauge_force";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "fermion_force";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "fermion_force_eoprec";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "set_zero_gaugemomentum";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "gaugemomentum_squarenorm";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "stout_smear_fermion_force";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
}
#endif

////////////////////////////////////////////////////
//Methods needed for the HMC-algorithm

void Opencl_Module_Hmc::generate_gaussian_gaugemomenta_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(generate_gaussian_gaugemomenta, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_new_p
	int clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 0, sizeof(cl_mem), &clmem_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 1, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( generate_gaussian_gaugemomenta , gs2, ls2);

	if(logger.beDebug()){
	  cl_mem force_tmp = create_rw_buffer(sizeof(hmc_float));
	  hmc_float resid;
	  this->set_float_to_gaugemomentum_squarenorm_device(clmem_p, force_tmp);
	  get_buffer_from_device(force_tmp, &resid, sizeof(hmc_float));
	  logger.debug() <<  "\tgaussian gaugemomenta:\t" << resid;
	  int clerr = clReleaseMemObject(force_tmp);
	  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	  if(resid != resid){
	    throw Print_Error_Message("calculation of gaussian gm gave nan! Aborting...", __FILE__, __LINE__);
	  }
	}

}

void Opencl_Module_Hmc::generate_spinorfield_gaussian(){
	if(get_parameters()->get_use_eo()== true){
		this->generate_gaussian_spinorfield_eoprec_device();
	}
	else{
		this->generate_gaussian_spinorfield_device();
	}
	return;
}

void Opencl_Module_Hmc::generate_gaussian_spinorfield_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(generate_gaussian_spinorfield, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_phi_inv, which can be done since the gaussian field is only needed in the beginning
	int clerr = clSetKernelArg(generate_gaussian_spinorfield, 0, sizeof(cl_mem), &clmem_phi_inv);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_spinorfield, 1, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(generate_gaussian_spinorfield  , gs2, ls2);

	if(logger.beDebug()){
	  cl_mem force_tmp = create_rw_buffer(sizeof(hmc_float));
	  hmc_float resid;
	  this->set_float_to_global_squarenorm_device(clmem_phi_inv, force_tmp);
	  get_buffer_from_device(force_tmp, &resid, sizeof(hmc_float));
	  logger.debug() <<  "\tinit gaussian spinorfield:\t" << resid;
	  int clerr = clReleaseMemObject(force_tmp);
	  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	  if(resid != resid){
	    throw Print_Error_Message("calculation of gaussian spinorfield gave nan! Aborting...", __FILE__, __LINE__);
	  }
	}

}

void Opencl_Module_Hmc::generate_gaussian_spinorfield_eoprec_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(generate_gaussian_spinorfield_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_phi_inv_eoporec, which can be done since the gaussian field is only needed in the beginning
	int clerr = clSetKernelArg(generate_gaussian_spinorfield_eoprec, 0, sizeof(cl_mem), &clmem_phi_inv_eoprec);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_spinorfield_eoprec, 1, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(generate_gaussian_spinorfield_eoprec  , gs2, ls2);

	if(logger.beDebug()){
	  cl_mem force_tmp = create_rw_buffer(sizeof(hmc_float));
	  hmc_float resid;
	  this->set_float_to_global_squarenorm_eoprec_device(clmem_phi_inv_eoprec, force_tmp);
	  get_buffer_from_device(force_tmp, &resid, sizeof(hmc_float));
	  logger.debug() <<  "\tinit gaussian spinorfield:\t" << resid;
	  int clerr = clReleaseMemObject(force_tmp);
	  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	  if(resid != resid){
	    throw Print_Error_Message("calculation of gaussian spinorfield gave nan! Aborting...", __FILE__, __LINE__);
	  }
	}

}

void Opencl_Module_Hmc::md_update_spinorfield()
{
	//suppose the initial gaussian field is saved in clmem_phi_inv (see above).
	//  then the "phi" = Dpsi from the algorithm is stored in clmem_phi
	//  which then has to be the source of the inversion
	if(get_parameters()->get_use_eo() == true){
	  Opencl_Module_Fermions::Qplus_eoprec (clmem_phi_inv_eoprec, clmem_phi_eoprec , *get_gaugefield());
		if(logger.beDebug()) print_info_inv_field(clmem_phi_eoprec, true, "\tinit field after update ");
	}
	else{
		Opencl_Module_Fermions::Qplus(clmem_phi_inv, clmem_phi , *get_gaugefield());
		if(logger.beDebug()) print_info_inv_field(clmem_phi, false, "\tinit field after update ");
	}
}

void Opencl_Module_Hmc::calc_fermion_force(usetimer * solvertimer)
{
	int converged = -1;
	/**
	 * @NOTE The force is up to this point always calculated using "non-eoprec" spinorfields.
	 * This is done this way since one would not save any operations using eoprec-fields, but would need
	 * another kernel instead.
	 * Therefore, the eoprec-fields are converted back to the normal format. 
	 * However, using an eoprec-force holds the possibility of saving memory, which would be relevant
	 * on e.g. a GPU.
	 * @NOTE A dummy-kernel and corresponding calling function is already implemented if one wishes to 
	 * 	change this one day.
	 */
	if(get_parameters()->get_use_eo() == true){
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters()->get_use_cg() == true) { 
			/**
			* The first inversion calculates
			* X_even = phi = (Qplusminus_eoprec)^-1 psi
			* out of
			* Qplusminus_eoprec phi_even = psi
			*/
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			/**
			  * Trial solution for the spinorfield
			  */
			set_eoprec_spinorfield_cold_device(get_clmem_inout_eoprec());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eoprec(), true, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::cg_eoprec(QplusQminus_eoprec_call, this->get_clmem_inout_eoprec(), this->get_clmem_phi_eoprec(), this->clmem_new_u, get_parameters()->get_force_prec());
			if (converged < 0){
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			}
			else logger.debug() << "\t\t\tsolver solved in "<< converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eoprec(), true, "\tinv. field after inversion ");

			/**
			 * Y_even is now just
			 * 	Y_even = (Qminus_eoprec) X_even = (Qminus_eoprec) (Qplusminus_eoprec)^-1 psi =
			 * 		= (Qplus_eoprec)^-1 psi
			 */
			Opencl_Module_Fermions::Qminus_eoprec(this->get_clmem_inout_eoprec(), clmem_phi_inv_eoprec, this->clmem_new_u);
		} 
		else {
			///@todo if wanted, solvertimer has to be used here..
			logger.debug() << "\t\tcalc fermion force ingredients using bicgstab with eoprec.";
			/**
			* The first inversion calculates
			* Y_even = phi = (Qplus_eoprec)^-1 psi
			* out of
			* Qplus_eoprec phi = psi
			* This is also the energy of the final field!
			*/
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			/**
			  * Trial solution for the spinorfield
			  */
			set_zero_spinorfield_eoprec_device(get_clmem_inout_eoprec());
			gamma5_eoprec_device(get_clmem_inout_eoprec());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eoprec(), true, "\tinv. field before inversion ");
			if(logger.beDebug()) print_info_inv_field(get_clmem_phi_eoprec(), true, "\tsource before inversion ");
			converged = Opencl_Module_Fermions::bicgstab_eoprec(Qplus_eoprec_call, this->get_clmem_inout_eoprec(), this->get_clmem_phi_eoprec(), this->clmem_new_u, get_parameters()->get_force_prec());
			if (converged < 0){
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			}
			else logger.debug() << "\t\t\tsolver solved in "<< converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eoprec(), true, "\tinv. field after inversion ");
			
			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout_eoprec(), clmem_phi_inv_eoprec, get_parameters()->get_eo_sf_buf_size());

			/**
			 * Now, one has to calculate
			 * X = (Qminus_eoprec)^-1 Y = (Qminus_eoprec)^-1 (Qplus_eoprec)^-1 psi = (QplusQminus_eoprec)^-1 psi ??
			 * out of
			 * Qminus_eoprec clmem_inout_eoprec = clmem_phi_inv_eoprec
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv_eoprec (see also above)!!
			 */

			//copy former solution to clmem_source
			copy_buffer_on_device(get_clmem_inout_eoprec(), get_clmem_source_even(), get_parameters()->get_eo_sf_buf_size());
			logger.debug() << "\t\t\tstart solver";

			//this sets clmem_inout cold as trial-solution
			set_eoprec_spinorfield_cold_device(get_clmem_inout_eoprec());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eoprec(), true, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::bicgstab_eoprec(Qminus_eoprec_call, get_clmem_inout_eoprec(), get_clmem_source_even(), clmem_new_u, get_parameters()->get_force_prec());
			if (converged < 0){
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			}
			else logger.debug() << "\t\t\tsolver solved in "<< converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eoprec(), true, "\tinv. field after inversion ");
		}
		/**
		 * At this point, one has calculated X_even and Y_odd.
		 * If one has a fermionmatrix
		 * 	M = R + D
		 * these are:
		 * 	X_odd = -R(-mu)_inv D X_even
		 * 	Y_odd = -R(mu)_inv D Y_even
		 */
		
		///@NOTE the following calculations could also go in a new function for convenience
		//calculate X_odd
		//therefore, clmem_tmp_eoprec_1 is used as intermediate state. The result is saved in clmem_inout, since
		//	this is used as a default in the force-function.
		if(get_parameters()->get_fermact() == WILSON){
			dslash_eoprec_device(get_clmem_inout_eoprec(), get_clmem_tmp_eoprec_1(), clmem_new_u, ODD);
			sax_eoprec_device(get_clmem_tmp_eoprec_1(), get_clmem_minusone(), get_clmem_tmp_eoprec_1());
		}
		else if(get_parameters()->get_fermact() == TWISTEDMASS){
			dslash_eoprec_device(get_clmem_inout_eoprec(), get_clmem_tmp_eoprec_1(), clmem_new_u, ODD);
			M_tm_inverse_sitediagonal_minus_device(get_clmem_tmp_eoprec_1(), get_clmem_tmp_eoprec_2());
			sax_eoprec_device(get_clmem_tmp_eoprec_2(), get_clmem_minusone(), get_clmem_tmp_eoprec_1());
		}
		this->convert_from_eoprec_device(get_clmem_inout_eoprec(), get_clmem_tmp_eoprec_1(), get_clmem_inout());
		
		//calculate Y_odd
		//therefore, clmem_tmp_eoprec_1 is used as intermediate state. The result is saved in clmem_phi_inv, since
		//	this is used as a default in the force-function.
		if(get_parameters()->get_fermact() == WILSON){
			dslash_eoprec_device(clmem_phi_inv_eoprec, get_clmem_tmp_eoprec_1(), clmem_new_u, ODD);
			sax_eoprec_device(get_clmem_tmp_eoprec_1(), get_clmem_minusone(), get_clmem_tmp_eoprec_1());
		}
		else if(get_parameters()->get_fermact() == TWISTEDMASS){
			dslash_eoprec_device(clmem_phi_inv_eoprec, get_clmem_tmp_eoprec_1(), clmem_new_u, ODD);
			M_tm_inverse_sitediagonal_minus_device(get_clmem_tmp_eoprec_1(), get_clmem_tmp_eoprec_2());
			sax_eoprec_device(get_clmem_tmp_eoprec_2(), get_clmem_minusone(), get_clmem_tmp_eoprec_1());
		}
		this->convert_from_eoprec_device(clmem_phi_inv_eoprec, get_clmem_tmp_eoprec_1(), clmem_phi_inv);
	}
	else{
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield 
		if(get_parameters()->get_use_cg() == true) { 
			/**
			* The first inversion calculates
			* X = phi = (Qplusminus)^-1 psi
			* out of
			* Qplusminus phi = psi
			*/
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			/**
			  * Trial solution for the spinorfield
			  */
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = Opencl_Module_Fermions::cg(QplusQminus_call, this->get_clmem_inout(), this->get_clmem_phi(), this->clmem_new_u, get_parameters()->get_force_prec());
			if (converged < 0){
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			}
			else logger.debug() << "\t\t\tsolver solved in "<< converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");

			/**
			 * Y is now just
			 * 	Y = (Qminus) X = (Qminus) (Qplusminus)^-1 psi =
			 * 		= (Qplus)^-1 psi
			 */
			Opencl_Module_Fermions::Qminus(this->get_clmem_inout(), clmem_phi_inv, this->clmem_new_u);
			
		} else  {
			logger.debug() << "\t\tcalc fermion force ingredients using bicgstab without eoprec";

			/**
			* The first inversion calculates
			* Y = phi = (Qplus)^-1 psi
			* out of
			* Qplus phi = psi
			* This is also the energy of the final field!
			*/
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			/**
			  * Trial solution for the spinorfield
			  */
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = Opencl_Module_Fermions::bicgstab(Qplus_call, this->get_clmem_inout(), this->get_clmem_phi(), this->clmem_new_u, get_parameters()->get_force_prec());
			if (converged < 0){
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			}
			else logger.debug() << "\t\t\tsolver solved in "<< converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout(), clmem_phi_inv, sizeof(spinor) * get_parameters()->get_spinorfieldsize());

			/**
			 * Now, one has to calculate
			 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
			 * out of
			 * Qminus clmem_inout = clmem_phi_inv
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv (see also above)!!
			 */

			//copy former solution to clmem_source
			copy_buffer_on_device(get_clmem_inout(), get_clmem_source(), sizeof(spinor) * get_parameters()->get_spinorfieldsize());
			logger.debug() << "\t\t\tstart solver";

			//this sets clmem_inout cold as trial-solution
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::bicgstab(Qminus_call, get_clmem_inout(), get_clmem_source(), clmem_new_u, get_parameters()->get_force_prec());
			if (converged < 0){
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			}
			else logger.debug() << "\t\t\tsolver solved in "<< converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
		}
	}
	logger.debug() << "\t\tcalc fermion_force...";
	//CP: this always calls fermion_force(Y,X) with Y = clmem_phi_inv, X = clmem_inout
	fermion_force_device();
}

void Opencl_Module_Hmc::calc_gauge_force()
{
	logger.debug() << "\t\tcalc gauge_force...";
	gauge_force_device();
}

hmc_float Opencl_Module_Hmc::calc_s_fermion(){
	logger.debug() << "calc final fermion energy...";
	//this function essentially performs the same steps as in the force-calculation, but with higher precision.
	//	therefore, comments are deleted here...
	//	Furthermore, in the bicgstab-case, the second inversions are not needed
	bool converged = false;
	if(get_parameters()->get_use_eo() == true){
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters()->get_use_cg() == true) { 
			logger.debug() << "\t\t\tstart solver";

			set_eoprec_spinorfield_cold_device(get_clmem_inout_eoprec());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eoprec(), true, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::cg_eoprec(QplusQminus_eoprec_call, this->get_clmem_inout_eoprec(), this->get_clmem_phi_eoprec(), this->clmem_new_u, get_parameters()->get_solver_prec());
			if (converged < 0){
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			}
			else logger.debug() << "\t\t\tsolver solved in "<< converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eoprec(), true, "\tinv. field after inversion ");

			Opencl_Module_Fermions::Qminus_eoprec(this->get_clmem_inout_eoprec(), clmem_phi_inv_eoprec, this->clmem_new_u);
		} 
		else {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			set_zero_spinorfield_eoprec_device(get_clmem_inout_eoprec());
			gamma5_eoprec_device(get_clmem_inout_eoprec());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eoprec(), true, "\tinv. field before inversion ");
			if(logger.beDebug()) print_info_inv_field(get_clmem_phi_eoprec(), true, "\tsource before inversion ");
			converged = Opencl_Module_Fermions::bicgstab_eoprec(Qplus_eoprec_call, this->get_clmem_inout_eoprec(), this->get_clmem_phi_eoprec(), this->clmem_new_u, get_parameters()->get_solver_prec());
			if (converged < 0){
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			}
			else logger.debug() << "\t\t\tsolver solved in "<< converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eoprec(), true, "\tinv. field after inversion ");
			
			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout_eoprec(), clmem_phi_inv_eoprec, get_parameters()->get_eo_sf_buf_size());

		}
	}
	else{
		if(get_parameters()->get_use_cg() == true) { 
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::cg(QplusQminus_call, this->get_clmem_inout(), this->get_clmem_phi(), this->clmem_new_u, get_parameters()->get_solver_prec());
			if (converged < 0){
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			}
			else logger.debug() << "\t\t\tsolver solved in "<< converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");

			Opencl_Module_Fermions::Qminus(this->get_clmem_inout(), clmem_phi_inv, this->clmem_new_u);
			
		} else  {

			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::bicgstab(Qplus_call, this->get_clmem_inout(), this->get_clmem_phi(), this->clmem_new_u, get_parameters()->get_solver_prec());
			if (converged < 0){
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			}
			else logger.debug() << "\t\t\tsolver solved in "<< converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout(), clmem_phi_inv, sizeof(spinor) * get_parameters()->get_spinorfieldsize());
		}
	}
	
	if(get_parameters()->get_use_eo() == true){
		set_float_to_global_squarenorm_eoprec_device(clmem_phi_inv_eoprec, clmem_s_fermion);
	}
	else{
		set_float_to_global_squarenorm_device(clmem_phi_inv, clmem_s_fermion);
	}
	hmc_float tmp;
	get_buffer_from_device(clmem_s_fermion, &tmp, sizeof(hmc_float));
	return tmp;
}

hmc_observables Opencl_Module_Hmc::metropolis(hmc_float rnd, hmc_float beta)
{
	//Calc Hamiltonian
	logger.debug() << "Calculate Hamiltonian";

	//Gauge-Part
	hmc_float tplaq, splaq, plaq;
	hmc_float tplaq_new, splaq_new, plaq_new;
	hmc_complex poly;
	hmc_complex poly_new;
	//In this call, the observables are calculated already with appropiate Weighting factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC)
	Opencl_Module::gaugeobservables(*get_gaugefield(), &plaq,  &tplaq, &splaq, &poly);
	Opencl_Module::gaugeobservables(clmem_new_u, &plaq_new,  &tplaq_new, &splaq_new, &poly_new);
	//plaq has to be divided by the norm-factor and multiplied by NC
	//  (because this is in the defintion of the gauge action and not in the normalization) to get s_gauge
	hmc_float factor = 2.0 / static_cast<hmc_float>(parameters->get_vol4d() * NDIM * (NDIM - 1) );
	/** NOTE: the minus here is introduced to fit tmlqcd!!! */
	hmc_float deltaH = -(plaq - plaq_new) * beta / factor;

	logger.debug() << "\tS_gauge(old field) = " << setprecision(10) << plaq << "\t" << plaq* beta  / factor;
	logger.debug() << "\tS_gauge(new field) = " << setprecision(10) << plaq_new << "\t" << plaq_new* beta / factor;
	logger.debug() << "\tdeltaS_gauge = " << setprecision(10) << deltaH;

	//Gaugemomentum-Part
	hmc_float p2, new_p2;
	set_float_to_gaugemomentum_squarenorm_device(clmem_p, clmem_p2);
	set_float_to_gaugemomentum_squarenorm_device(clmem_new_p, clmem_new_p2);
	Opencl_Module_Hmc::get_buffer_from_device(clmem_p2, &p2, sizeof(hmc_float));
	Opencl_Module_Hmc::get_buffer_from_device(clmem_new_p2, &new_p2, sizeof(hmc_float));
	//the energy is half the squarenorm
	deltaH += 0.5 * (p2 - new_p2);

	logger.debug() << "\tS_gaugemom(old field) = " << setprecision(10) << 0.5 * p2;
	logger.debug() << "\tS_gaugemom(new field) = " << setprecision(10) << 0.5 * new_p2;
	logger.debug() << "\tdeltaS_gaugemom = " << setprecision(10) << 0.5 * (p2 - new_p2);

	//Fermion-Part:
	hmc_float spinor_energy_init, s_fermion;
	//initial energy has been computed in the beginning...
	Opencl_Module_Hmc::get_buffer_from_device(clmem_energy_init, &spinor_energy_init, sizeof(hmc_float));
	// sum_links phi*_i (M^+M)_ij^-1 phi_j
	s_fermion = calc_s_fermion();
	deltaH += spinor_energy_init - s_fermion;

	logger.debug() << "\tS_ferm(old field) = " << setprecision(10) <<  spinor_energy_init;
	logger.debug() << "\tS_ferm(new field) = " << setprecision(10) << s_fermion;
	logger.debug() << "\tdeltaS_ferm = " << spinor_energy_init - s_fermion;

	//Metropolis-Part
	hmc_float compare_prob;
	if(deltaH < 0) {
		compare_prob = exp(deltaH);
	} else {
		compare_prob = 1.0;
	}
	logger.debug() << "\tdeltaH = " << deltaH << "\tAcc-Prop = " << compare_prob;
	hmc_observables tmp;
	if(rnd <= compare_prob) {
		tmp.accept = 1;
		tmp.plaq = plaq_new;
		tmp.tplaq = tplaq_new;
		tmp.splaq = splaq_new;
		tmp.poly = poly_new;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
	} else {
		tmp.accept = 0;
		tmp.plaq = plaq;
		tmp.tplaq = tplaq;
		tmp.splaq = splaq;
		tmp.poly = poly;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
	}

	return tmp;
}

void Opencl_Module_Hmc::calc_spinorfield_init_energy()
{
	//Suppose the initial spinorfield is saved in phi_inv
	//  it is created in generate_gaussian_spinorfield_device
	if(get_parameters()->get_use_eo() == true){
		Opencl_Module_Fermions::set_float_to_global_squarenorm_eoprec_device(clmem_phi_inv_eoprec, clmem_energy_init);
	}
	else{
		Opencl_Module_Fermions::set_float_to_global_squarenorm_device(clmem_phi_inv, clmem_energy_init);
	}
}

void Opencl_Module_Hmc::md_update_gaugemomentum_device(hmc_float eps)
{
	//__kernel void md_update_gaugemomenta(hmc_float eps, __global ae * p_inout, __global ae* force_in){
	hmc_float tmp = eps;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(md_update_gaugemomenta, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(md_update_gaugemomenta, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(md_update_gaugemomenta, 1, sizeof(cl_mem), &clmem_new_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(md_update_gaugemomenta, 2, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(md_update_gaugemomenta , gs2, ls2);
}

void Opencl_Module_Hmc::md_update_gaugefield_device(hmc_float eps)
{
	// __kernel void md_update_gaugefield(hmc_float eps, __global ae * p_in, __global ocl_s_gaugefield * u_inout){
	hmc_float tmp = eps;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(md_update_gaugefield, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_force
	int clerr = clSetKernelArg(md_update_gaugefield, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(md_update_gaugefield, 1, sizeof(cl_mem), &clmem_new_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(md_update_gaugefield, 2, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( md_update_gaugefield , gs2, ls2);
}

void Opencl_Module_Hmc::set_zero_clmem_force_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(md_update_gaugemomenta, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_force
	int clerr = clSetKernelArg(set_zero_gaugemomentum, 0, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( set_zero_gaugemomentum , gs2, ls2);
}

void Opencl_Module_Hmc::gauge_force_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gauge_force, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gauge_force, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(gauge_force, 1, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( gauge_force , gs2, ls2);
}

void Opencl_Module_Hmc::fermion_force_device()
{
	//fermion_force(field, Y, X, out);
  cl_mem tmp = get_clmem_inout();
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(fermion_force, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(fermion_force, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 1, sizeof(cl_mem), &clmem_phi_inv);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 2, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 3, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( fermion_force , gs2, ls2);

	if(logger.beDebug()){
	  cl_mem force_tmp = create_rw_buffer(sizeof(hmc_float));
	  hmc_float resid;
	  this->set_float_to_gaugemomentum_squarenorm_device(clmem_force, force_tmp);
	  get_buffer_from_device(force_tmp, &resid, sizeof(hmc_float));
	  logger.debug() <<  "\tforce:\t" << resid;
	  int clerr = clReleaseMemObject(force_tmp);
	  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	  if(resid != resid){
	    throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
	  }
	}
}

void Opencl_Module_Hmc::fermion_force_eoprec_device(cl_mem Y, cl_mem X, int evenodd)
{
	//fermion_force(field, Y, X, out);
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(fermion_force_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(fermion_force_eoprec, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eoprec, 1, sizeof(cl_mem), &Y);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eoprec, 2, sizeof(cl_mem), &X);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eoprec, 3, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eoprec, 3, sizeof(int), &evenodd);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	enqueueKernel( fermion_force_eoprec , gs2, ls2);


}

void Opencl_Module_Hmc::stout_smeared_fermion_force_device(cl_mem * gf_intermediate)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(stout_smear_fermion_force, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
}

void Opencl_Module_Hmc::set_float_to_gaugemomentum_squarenorm_device(cl_mem clmem_in, cl_mem clmem_out)
{
	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gaugemomentum_squarenorm, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	int clerr = clSetKernelArg(gaugemomentum_squarenorm, 0, sizeof(cl_mem), &clmem_in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

//  /** @todo add reduction */
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 1, sizeof(cl_mem), &clmem_out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(gaugemomentum_squarenorm  , gs2, ls2);
}
