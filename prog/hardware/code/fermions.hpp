/** @file
 * Basic OpenCL functionality
 */

#ifndef _HARDWARE_CODE_FERMIONS_
#define _HARDWARE_CODE_FERMIONS_

#include "opencl_module.hpp"

#include "../buffers/plain.hpp"
#include "../buffers/su3.hpp"
#include "../buffers/spinor.hpp"
#include "../../host_use_timer.h"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Fermions : public hardware::code::Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Fermions();


	//    fermionmatrix operations
	//    non-eo
	//        explicit
	void M_wilson_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF) const;
	void M_tm_plus_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	void M_tm_minus_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	void gamma5_device(const hardware::buffers::Plain<spinor> * inout) const;
	//    eo
	//        explicit
	void gamma5_eo_device(const hardware::buffers::Spinor * inout) const;
	void M_tm_inverse_sitediagonal_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	void M_tm_sitediagonal_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	void M_tm_inverse_sitediagonal_minus_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	void M_tm_sitediagonal_minus_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	void dslash_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF) const;
	//        merged
//	void Aee_AND_gamma5_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
//	void Aee_minus_AND_gamma5_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void dslash_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF) const;
	void dslash_AND_M_tm_inverse_sitediagonal_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	void dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	void M_tm_sitediagonal_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	void M_tm_sitediagonal_minus_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
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

	ClSourcePackage get_sources() const noexcept;

protected:
	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;

private:
	/**
	 * Default constructor.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Fermions(const meta::Inputparameters& params, hardware::Device * device);

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

	////////////////////////////////////
	//kernels, sorted roughly by groups
	//fermionmatrix
	cl_kernel M_wilson;
	cl_kernel gamma5;
	cl_kernel M_tm_plus;
	cl_kernel M_tm_minus;
	cl_kernel gamma5_eo;
	cl_kernel M_tm_sitediagonal;
	cl_kernel M_tm_inverse_sitediagonal;
	cl_kernel M_tm_sitediagonal_minus;
	cl_kernel M_tm_inverse_sitediagonal_minus;
	cl_kernel dslash_eo;
	cl_kernel dslash_AND_gamma5_eo;
	cl_kernel dslash_AND_M_tm_inverse_sitediagonal_eo;
	cl_kernel dslash_AND_M_tm_inverse_sitediagonal_minus_eo;
	cl_kernel M_tm_sitediagonal_AND_gamma5_eo;
	cl_kernel M_tm_sitediagonal_minus_AND_gamma5_eo;

	ClSourcePackage sources;
};

}

}

#endif // _HARDWARE_CODE_FERMIONS_
