/** @file
 * Basic OpenCL functionality
 */
#ifndef _OPENCLMODULEFERMIONSH_
#define _OPENCLMODULEFERMIONSH_

#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "globaldefs.h"
#include "types.h"
#include "types_fermions.h"
#include "host_use_timer.h"
#include "inputparameters.h"
#include "opencl_compiler.hpp"

#include "opencl_module.h"
#include "opencl_module_ran.h"
#include "opencl_module_spinors.h"

#include "exceptions.h"

class Opencl_Module_Fermions;

/**
 * this is a workaround to be able to pass a (fermionmatrix-)function, which are methods in this class,
 * to another function inside this class.
 * This type points to a helper-function, which then calls the wanted function.
 */
typedef void (*matrix_function_call) (Opencl_Module_Fermions* that, cl_mem in, cl_mem out, cl_mem gf);
void M_call(Opencl_Module_Fermions* that, cl_mem in, cl_mem out, cl_mem gf);
void Qplus_call(Opencl_Module_Fermions* that, cl_mem in, cl_mem out, cl_mem gf);
void Qminus_call(Opencl_Module_Fermions* that, cl_mem in, cl_mem out, cl_mem gf);
void QplusQminus_call(Opencl_Module_Fermions* that, cl_mem in, cl_mem out, cl_mem gf);
void Aee_call(Opencl_Module_Fermions* that, cl_mem in, cl_mem out, cl_mem gf);


/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Fermions : public Opencl_Module_Spinors {
public:

	// OpenCL specific methods needed for building/compiling the OpenCL program
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual void fill_collect_options(stringstream* collect_options);
	/**
	 * Collect the buffers to generate for OpenCL.
	 * Virtual method, allows to include more buffers in inherited classes.
	 */
	virtual void fill_buffers();
	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	virtual void fill_kernels();
	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	virtual void clear_kernels();
	/**
	 * Clear out the buffers,
	 * Virtual method, allows to clear additional buffers in inherited classes.
	 */
	virtual void clear_buffers();

	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param dev_type type of device on which the kernel should be executed
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups);


	//    fermionmatrix operations
	//    non-eoprec
	//        compound
	void M(cl_mem in, cl_mem out, cl_mem gf);
	void Qplus(cl_mem in, cl_mem out, cl_mem gf);
	void Qminus(cl_mem in, cl_mem out, cl_mem gf);
	void QplusQminus(cl_mem in, cl_mem out, cl_mem gf);
	//        explicit
	void M_wilson_device(cl_mem in, cl_mem out, cl_mem gf);
	void M_tm_plus_device(cl_mem in, cl_mem out, cl_mem gf);
	void M_tm_minus_device(cl_mem in, cl_mem out, cl_mem gf);
	void gamma5_device(cl_mem inout);
	//    eoprec
	//        compound
	void Qplus_eoprec(cl_mem in, cl_mem out, cl_mem gf);
	void Qminus_eoprec(cl_mem in, cl_mem out, cl_mem gf);
	void Aee(cl_mem in, cl_mem out, cl_mem gf);
	//        explicit
	void gamma5_eoprec_device(cl_mem inout);
	void M_tm_inverse_sitediagonal_device(cl_mem in, cl_mem out);
	void M_tm_sitediagonal_device(cl_mem in, cl_mem out);
	void dslash_eoprec_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd);

	//    solver operations
	//    non-eoprec
	/// this calls the solver according to parameter settings using the fermionmatrix f
	void solver(matrix_function_call f, cl_mem inout, cl_mem source, cl_mem gf, usetimer * solvertimer);
	/// this executes the bicgstab on the device, using the fermionmatrix f
	bool bicgstab(matrix_function_call f, cl_mem inout, cl_mem source, cl_mem gf);
	/// this executes the cg on the device, using the fermionmatrix f
	bool cg(matrix_function_call f, cl_mem inout, cl_mem source, cl_mem gf);
	//    eoprec
	/// this executes the eoprec bicgstab on the device, using the fermionmatrix f
	bool bicgstab_eoprec(matrix_function_call f, cl_mem inout, cl_mem source, cl_mem gf);
	bool cg_eoprec(matrix_function_call f, cl_mem inout, cl_mem source, cl_mem gf);

	/////////////////////////////////////////////////
	//functions to get private variables
	cl_mem get_clmem_inout();
	cl_mem get_clmem_source();
	cl_mem get_clmem_tmp();

	cl_mem get_clmem_inout_eoprec();
	cl_mem get_clmem_tmp_eoprec_1();
	cl_mem get_clmem_source_even();
	cl_mem get_clmem_source_odd();

	cl_mem get_clmem_minusone();

protected:

#ifdef _PROFILING_
	//CP: if PROFILING is activated, one needs a timer for each kernel
	//fermionmatrix
	usetimer timer_M_wilson;
	usetimer timer_gamma5;
	usetimer timer_M_tm_plus;
	usetimer timer_M_tm_minus;
	usetimer timer_gamma5_eoprec;
	usetimer timer_M_tm_sitediagonal;
	usetimer timer_M_tm_inverse_sitediagonal;
	usetimer timer_dslash_eoprec;
	usetimer timer_M_tm_sitediagonal_minus;
	usetimer timer_M_tm_inverse_sitediagonal_minus;

	/**
	 * Return the timer connected to a specific kernel.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual usetimer* get_timer(const char * in);

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual int get_read_write_size(const char * in, inputparameters * parameters);

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(std::string filename);

#endif

protected:
private:
	////////////////////////////////////
	//kernels, sorted roughly by groups
	//fermionmatrix
	cl_kernel M_wilson;
	cl_kernel gamma5;
	cl_kernel M_tm_plus;
	cl_kernel M_tm_minus;
	cl_kernel gamma5_eoprec;
	cl_kernel M_tm_sitediagonal;
	cl_kernel M_tm_inverse_sitediagonal;
	cl_kernel M_tm_sitediagonal_minus;
	cl_kernel M_tm_inverse_sitediagonal_minus;
	cl_kernel dslash_eoprec;
	//CP: variables for normal solver
	cl_mem clmem_inout;
	cl_mem clmem_rn;
	cl_mem clmem_rhat;
	cl_mem clmem_v;
	cl_mem clmem_p;
	cl_mem clmem_s;
	cl_mem clmem_t;
	cl_mem clmem_aux;
	//this is needed in QplusQminus as a temporary field
	cl_mem clmem_tmp;

	//CP: variables for eoprec solver
	cl_mem clmem_inout_eoprec;
	cl_mem clmem_source;
	cl_mem clmem_source_even;
	cl_mem clmem_source_odd;
	cl_mem clmem_rn_eoprec;
	cl_mem clmem_rhat_eoprec;
	cl_mem clmem_v_eoprec;
	cl_mem clmem_p_eoprec;
	cl_mem clmem_s_eoprec;
	cl_mem clmem_t_eoprec;
	cl_mem clmem_aux_eoprec;
	cl_mem clmem_tmp_eoprec_1;
	cl_mem clmem_tmp_eoprec_2;

	cl_mem clmem_rho;
	cl_mem clmem_rho_next;
	cl_mem clmem_alpha;
	cl_mem clmem_omega;
	cl_mem clmem_beta;
	cl_mem clmem_tmp1;
	cl_mem clmem_tmp2;
	cl_mem clmem_one;
	cl_mem clmem_minusone;

	cl_mem clmem_kappa_cmplx;// = {kappa, 0.};
	cl_mem clmem_resid;
	cl_mem clmem_trueresid;


private:

};

#endif //OPENCLMODULEFERMIONSH
