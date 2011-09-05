#include "opencl_module_kappa.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Kappa::fill_collect_options(stringstream* collect_options){
  Opencl_Module::fill_collect_options(collect_options);
  return;
}


void Opencl_Module_Kappa::fill_buffers()
{

  Opencl_Module::fill_buffers();

	cl_int clerr = CL_SUCCESS;

	cl_uint num_groups;
	size_t local_work_size;
	size_t global_work_size;
	get_work_sizes(&local_work_size, &global_work_size, &num_groups, Opencl_Module::get_device_type() );


	cout << "Create buffer for transport coefficient kappa_clover..." << endl;
	clmem_kappa_clover = clCreateBuffer(get_context(), CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clCreateBuffer",__FILE__,__LINE__);

	cout << "Create scratch buffer..." << endl;
	int global_buf_size_float = sizeof(hmc_float) * num_groups;
	clmem_kappa_clover_buf_glob = clCreateBuffer(get_context(), CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clCreateBuffer",__FILE__,__LINE__);

  return;
}

void Opencl_Module_Kappa::fill_kernels()
{
  Opencl_Module::fill_kernels();

	cout << "Create TK clover kernels..." << endl;
	kappa_clover_gpu = createKernel("kappa_clover_gpu") << basic_opencl_code << "opencl_tk_kappa.cl";

  return;	
}


void Opencl_Module_Kappa::run_kappa_clover(const hmc_float beta)
{
	//variables
	cl_int clerr = CL_SUCCESS;

	size_t local_work_size;
	size_t global_work_size;
	cl_uint num_groups;
	//CP: This has no effect yet!!
	get_work_sizes(&local_work_size, &global_work_size, &num_groups, Opencl_Module::get_device_type());


	clerr = clSetKernelArg(kappa_clover_gpu, 0, sizeof(cl_mem), get_gaugefield());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(kappa_clover_gpu, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(kappa_clover_gpu, 2, sizeof(cl_mem), &clmem_kappa_clover_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clKernelArg",__FILE__,__LINE__);

	enqueueKernel(kappa_clover_gpu, global_work_size, local_work_size);

	/*
	//read out values
	clerr = clEnqueueReadBuffer(get_queue(), clmem_kappa_clover_buf_glob, CL_FALSE, 0, sizeof(hmc_float), kappa_clover_out, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clEnqueueReadBuffer",__FILE__,__LINE__);
	*/
	// wait for results to have been read back
	//don't do that anymore ;-)
	//	clFinish(queue);
	//	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFinish",__FILE__,__LINE__);

	return;

}



void Opencl_Module_Kappa::clear_kernels()
{
  Opencl_Module::clear_kernels();

		cl_int clerr = clReleaseKernel(kappa_clover_gpu);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);

  return;
}

void Opencl_Module_Kappa::clear_buffers()
{
  Opencl_Module::clear_buffers();

  cl_int clerr = CL_SUCCESS;

  clerr = clReleaseMemObject(clmem_kappa_clover);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
  
  clerr = clReleaseMemObject(clmem_kappa_clover_buf_glob);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
  
  return;
}


