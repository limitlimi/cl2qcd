/** @file
 * Heatbath for OpenCL
 */
#ifndef _HARDWARE_CODE_HEATBATH_
#define _HARDWARE_CODE_HEATBATH_

#include "opencl_module.hpp"
#include "../../types.h"
#include "../../meta/inputparameters.hpp"
#include "../buffers/su3.hpp"
#include "../buffers/prng_buffer.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * Adds heatbath to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Heatbath : public hardware::code::Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Heatbath();

	/**
	 * Perform one heatbath step.
	 */
	void run_heatbath(const hardware::buffers::SU3 * gaugefield, const hardware::buffers::PRNGBuffer * prng) const;

	/**
	 * Perform one overrelaxation step.
	 */
	void run_overrelax(const hardware::buffers::SU3 * gaugefield, const hardware::buffers::PRNGBuffer * prng) const;

	/**
	 * Add specific work_size determination for this child class
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;


private:
	/**
	 * @param[in] params points to an instance of inputparameters
	 */
	Heatbath(const meta::Inputparameters& params, hardware::Device * device);

	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	void fill_kernels();

	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	void clear_kernels();

	cl_kernel heatbath_odd;
	cl_kernel heatbath_even;
	cl_kernel overrelax_odd;
	cl_kernel overrelax_even;


	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param parameters inputparameters
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string& in) const override;

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string& in) const override;


};

}

}

#endif // _HARDWARE_CODE_HEATBATH_