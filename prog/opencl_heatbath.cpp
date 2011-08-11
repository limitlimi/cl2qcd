#include "opencl_heatbath.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

hmc_error Opencl_heatbath::fill_collect_options(stringstream* collect_options)
{
  Opencl::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters()->get_beta();
  return HMC_SUCCESS;
}


hmc_error Opencl_heatbath::fill_buffers()
{
  Opencl::fill_buffers();

	return HMC_SUCCESS;
}

hmc_error Opencl_heatbath::clear_buffers(){

  Opencl::clear_buffers();
  if(clReleaseMemObject(clmem_rndarray) != CL_SUCCESS) exit(HMC_OCLERROR);

  return HMC_SUCCESS;

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
}

hmc_error Opencl_heatbath::clear_kernels()
{
  Opencl::clear_kernels();

	if(clReleaseKernel(heatbath_even) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(heatbath_odd) != CL_SUCCESS) exit(HMC_OCLERROR);

	if(clReleaseKernel(overrelax_even) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(overrelax_odd) != CL_SUCCESS) exit(HMC_OCLERROR);

	return HMC_SUCCESS;
}

hmc_error Opencl_heatbath::run_heatbath()
{
	cl_int clerr = CL_SUCCESS;

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = min(VOLSPACE * NTIME / 2, NUMRNDSTATES);
	else
		global_work_size = min(max_compute_units, (cl_uint) NUMRNDSTATES);

	clerr = clSetKernelArg(heatbath_even, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 at heatbath_even failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_even, 2, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 at heatbath_even failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_even, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg2 at heatbath_even failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(heatbath_even, global_work_size);
	}

	clerr = clSetKernelArg(heatbath_odd, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 at heatbath_odd failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_odd, 2, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 at heatbath_odd failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_odd, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg2 at heatbath_odd failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(heatbath_odd, global_work_size);
	}
	clFinish(queue);
	return HMC_SUCCESS;

}

hmc_error Opencl_heatbath::run_overrelax()
{
	cl_int clerr = CL_SUCCESS;

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = min(VOLSPACE * NTIME / 2, NUMRNDSTATES);
	else
		global_work_size = min(max_compute_units, (cl_uint) NUMRNDSTATES);

	clerr = clSetKernelArg(overrelax_even, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_even, 2, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_even, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg4 failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(overrelax_even, global_work_size);
	}

	clerr = clSetKernelArg(overrelax_odd, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg5 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_odd, 2, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg7 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_odd, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg8 failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(overrelax_odd, global_work_size);
	}
	clFinish(queue);
	return HMC_SUCCESS;
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
