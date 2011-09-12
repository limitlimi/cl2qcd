#include "opencl_module_correlator.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Correlator::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Spinors::fill_collect_options(collect_options);
	//CP: give kappa and its negative value
	hmc_float kappa_tmp = get_parameters()->get_kappa();
	*collect_options << " -DKAPPA=" << kappa_tmp;
	*collect_options << " -DMKAPPA=" << -kappa_tmp;
	
	return;
}


void Opencl_Module_Correlator::fill_buffers()
{

	Opencl_Module_Spinors::fill_buffers();
	logger.trace() << "init buffer for correlators...";
	int clerr = CL_SUCCESS;

	int spinorfield_size = sizeof(spinor) * SPINORFIELDSIZE;
	int num_sources;
	if( get_parameters()->get_use_pointsource() == true){
		num_sources = 12;
	}
	else{
		num_sources = get_parameters()->get_num_sources();
	}
	clmem_corr = create_rw_buffer(spinorfield_size*num_sources);
	clmem_source = create_rw_buffer(spinorfield_size);

	return;
}

void Opencl_Module_Correlator::fill_kernels()
{
	Opencl_Module_Spinors::fill_kernels();
	basic_correlator_code = basic_fermion_code;
	logger.debug() << "Create correlator kernels...";
	
	if(get_parameters()->get_use_pointsource() == true)
		create_point_source = createKernel("create_point_source") << basic_fermion_code << "spinorfield_point_source.cl";
	else
		create_stochastic_source = createKernel("create_stochastic_source") << basic_fermion_code << "spinorfield_stochastic_source.cl";
	ps_correlator = createKernel("ps_correlator") << basic_fermion_code << "fermionobservables.cl";

	return;
}

void Opencl_Module_Correlator::clear_kernels()
{
	Opencl_Module_Spinors::clear_kernels();
	int clerr = clReleaseKernel(ps_correlator);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);
	clerr = clReleaseKernel(create_point_source);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);
	clerr = clReleaseKernel(create_stochastic_source);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);
	return;
}

void Opencl_Module_Correlator::clear_buffers()
{
	Opencl_Module_Spinors::clear_buffers();
	int clerr = clReleaseMemObject(clmem_source);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clMemObject",__FILE__,__LINE__);
	clerr = clReleaseMemObject(clmem_corr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clMemObject",__FILE__,__LINE__);
	return;
}

void Opencl_Module_Correlator::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
  Opencl_Module_Spinors::get_work_sizes(kernel, dev_type, ls, gs, num_groups);
  
	return;
}

cl_mem Opencl_Module_Correlator::get_clmem_source()
{
	return clmem_source;
}

cl_mem Opencl_Module_Correlator::get_clmem_corr()
{
	return clmem_corr;
}

void Opencl_Module_Correlator::create_point_source_device(cl_mem inout, int i, int spacepos, int timepos)
{
	set_zero_spinorfield_device(inout);

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_point_source, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_point_source, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(create_point_source, 1, sizeof(int), &i);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(create_point_source, 2, sizeof(int), &spacepos);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(create_point_source, 3, sizeof(int), &timepos);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	enqueueKernel( create_point_source, gs2, ls2);
}

void Opencl_Module_Correlator::create_stochastic_source_device(cl_mem inout)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_stochastic_source, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_stochastic_source, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	throw Opencl_Error(clerr,"stochastic source not yet implemented!!",__FILE__,__LINE__);
	enqueueKernel( create_stochastic_source, gs2, ls2);
}

void Opencl_Module_Correlator::ps_correlator_device(cl_mem in){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(ps_correlator, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(ps_correlator,0,sizeof(cl_mem),&in); 
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

  ///@todo improve this
	ls2 = 1;
	gs2 = 1;
	enqueueKernel(ps_correlator , gs2, ls2);
}

#ifdef _PROFILING_
usetimer* Opencl_Module_Correlator::get_timer(char * in){
	usetimer *noop = NULL;
	noop = Opencl::get_timer(in);
	if(noop != NULL) return noop;
	
	if (strcmp(in, "create_point_source") == 0){
    return &this->timer_create_point_source;
	}
	if (strcmp(in, "create_stochastic_source") == 0){
    return &this->timer_create_stochastic_source;
	}
	if (strcmp(in, "ps_correlator") == 0){
    return &this->timer_ps_correlator;
	}
	
	//if the kernelname has not matched, return NULL
	else{
		return NULL;
	}
}

int Opencl_Module_Correlator::get_read_write_size(char * in, inputparameters * parameters){
	Opencl::get_read_write_size(in, parameters);
		//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	int R = (*parameters).get_mat_size();
	int S;
	if((*parameters).get_use_eo() == 1)
	  S = EOPREC_SPINORFIELDSIZE;
	else
	  S = SPINORFIELDSIZE;
	//this is the same as in the function above
	if (strcmp(in, "create_point_source") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "create_stochastic_source") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "ps_correlator") == 0){
    return 1000000000000000000000000;
	}
	return 0;
}

void Opencl_Module_Correlator::print_profiling(std::string filename){
	Opencl::print_profiling(filename);
	char * kernelName;
	kernelName = "create_point_source";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "create_stochastic_source";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "ps_correlator";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	
}
#endif
