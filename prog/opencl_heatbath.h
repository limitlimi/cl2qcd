/** @file
 * OpenCL device including heatbath.
 */
#ifndef _MYOPENCLHEATBATHH_
#define _MYOPENCLHEATBATHH_

#include "opencl.h"


/**
 * An OpenCL device for heatbath.
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_heatbath : public Opencl {
public:

	//Calculations, calls to kernels
	/**
	 * Perform one heatbath step.
	 */
	hmc_error run_heatbath(const hmc_float beta, usetimer * const timer);

	/**
	 * Perform one overrelaxation step.
	 */
	hmc_error run_overrelax(const hmc_float beta, usetimer * const timer);

	/**
	 * Calculate plaquette and polyakov.
	 *
	 * @param[out] plaq Storage for result of plaquette calculation
	 * @param[out] tplaq Storage for result of plaquette calculation
	 * @param[out] splaq Storage for result of plaquette calculation
	 * @param[out] pol Storage for result of polyakov calculation
	 * @param[in,out] timer1 Timer into which to aggregate plaquette calculation time
	 * @param[in,out] timer2 Timer into which to aggregate polyakov calculation time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	//protected:

	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual hmc_error fill_collect_options(stringstream* collect_options);

	/**
	 * Collect the buffers to generate for OpenCL.
	 * Virtual method, allows to include more buffers in inherited classes.
	 */
	virtual hmc_error fill_buffers();

	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	virtual void fill_kernels();

	//private:


	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	virtual hmc_error clear_kernels();

	/**
	 * Clear out the buffers,
	 * Virtual method, allows to clear additional buffers in inherited classes.
	 */
	virtual hmc_error clear_buffers();

	/**
	 * Copy the RNG state to the appropriate OpenCL buffer.
	 *
	 * @param host_rndarray The RNG state to copy
	 * @param timer The timer to use to measure the copying time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error copy_rndarray_to_device(hmc_rndarray host_rndarray,  usetimer* timer);

	/**
	 * Copy the RNG state from the OpenCL buffer.
	 *
	 * @param[out] rndarray The RNG copy target
	 * @param[in,out] timer The timer to use to measure the copying time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error copy_rndarray_from_device(hmc_rndarray rndarray, usetimer* timer);



	///////////////////////////////////////////////////////////
	//LZ what follows should eventually be private
	//heatbath variables

	cl_mem clmem_rndarray;

	cl_kernel heatbath_odd;
	cl_kernel heatbath_even;
	cl_kernel overrelax_odd;
	cl_kernel overrelax_even;

protected:

private:

};

#endif /* _MYOPENCLHEATBATHH_ */
