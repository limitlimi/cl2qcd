#include "opencl_heatbath.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_heatbath::fill_collect_options(stringstream* collect_options)
{
  Opencl::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters()->get_beta();
  return;
}


void Opencl_heatbath::fill_buffers()
{
  Opencl::fill_buffers();

	return;
}

void Opencl_heatbath::clear_buffers(){

  Opencl::clear_buffers();
  cl_int clerr = clReleaseMemObject(clmem_rndarray);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);

  return;

}

void Opencl_heatbath::fill_kernels()
{
  Opencl::fill_kernels();

  logger.debug() << "Create heatbath kernels...";
  heatbath_even = createKernel("heatbath_even") << basic_opencl_code << "random.cl" << "update_heatbath.cl";
  heatbath_odd = createKernel("heatbath_odd") << basic_opencl_code << "random.cl" << "update_heatbath.cl";
  
  logger.debug() << "Create overrelax kernels...";
  overrelax_even = createKernel("overrelax_even") << basic_opencl_code << "random.cl" << "overrelax.cl";
  overrelax_odd = createKernel("overrelax_odd") << basic_opencl_code << "random.cl" << "overrelax.cl";

#ifdef _PROFILING_
	//init timers
	usetimer noop;
	timer_heatbath_even.reset();
#endif

	return;

}

void Opencl_heatbath::clear_kernels()
{
  Opencl::clear_kernels();

  cl_int clerr = CL_SUCCESS;

  clerr = clReleaseKernel(heatbath_even);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);

  clerr = clReleaseKernel(heatbath_odd);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);

  clerr = clReleaseKernel(overrelax_even);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);

  clerr = clReleaseKernel(overrelax_odd);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);

  return;
}

void Opencl_heatbath::run_heatbath()
{
	cl_int clerr = CL_SUCCESS;

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
	  global_work_size = min(VOLSPACE * NTIME / 2, this->Opencl::get_num_rndstates());
	else
	  global_work_size = min(max_compute_units, (cl_uint) this->Opencl::get_num_rndstates());

	clerr = clSetKernelArg(heatbath_even, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(heatbath_even, 2, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_even, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

		enqueueKernel(heatbath_even, global_work_size);
	}

	clerr = clSetKernelArg(heatbath_odd, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(heatbath_odd, 2, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_odd, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);
		enqueueKernel(heatbath_odd, global_work_size);
	}
	clerr = clFinish(queue);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFinish",__FILE__,__LINE__);
	return;

}

void Opencl_heatbath::run_overrelax()
{
	cl_int clerr = CL_SUCCESS;

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
	  global_work_size = min(VOLSPACE * NTIME / 2, this->Opencl::get_num_rndstates());
	else
	  global_work_size = min(max_compute_units, (cl_uint) this->Opencl::get_num_rndstates());

	clerr = clSetKernelArg(overrelax_even, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(overrelax_even, 2, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_even, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

		enqueueKernel(overrelax_even, global_work_size);
	}

	clerr = clSetKernelArg(overrelax_odd, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(overrelax_odd, 2, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_odd, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

		enqueueKernel(overrelax_odd, global_work_size);
	}
	clerr = clFinish(queue);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFinish",__FILE__,__LINE__);
	return;
}

#ifdef _PROFILING_
usetimer* Opencl_heatbath::get_timer(char * in){
	usetimer *noop = NULL;
	noop = Opencl::get_timer(in);
	if(noop != NULL) return noop;
	if (strcmp(in, "heatbath_even") == 0){
    return &this->timer_heatbath_even;
	}
	if (strcmp(in, "heatbath_odd") == 0){
    return &this->timer_heatbath_odd;
	}
	if (strcmp(in, "overrelax_even") == 0){
    return &this->timer_overrelax_even;
	}
	if (strcmp(in, "overrelax_odd") == 0){
    return &this->timer_overrelax_odd;
	}
	//if the kernelname has not matched, return NULL
	else{
		return NULL;
	}
}

int Opencl_heatbath::get_read_write_size(char * in, inputparameters * parameters){
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
	if (strcmp(in, "heatbath_even") == 0){
    return VOL4D*D*R + 1;
	}
	if (strcmp(in, "heatbath_odd") == 0){
    return VOL4D*D*R + 1;
	}
	if (strcmp(in, "overrelax_even") == 0){
    return 48*VOL4D *D*R + 1;
	}
	if (strcmp(in, "overrelax_odd") == 0){
    return 48*VOL4D *D*R + 1;	
	}
}

void Opencl_heatbath::print_profiling(std::string filename){
	Opencl::print_profiling(filename);
	char * kernelName;
	kernelName = "heatbath_even";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "heatbath_odd";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "overrelax_even";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "overrelax_odd";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
}
#endif
