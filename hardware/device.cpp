/** @file
 * Implementation of the hardware::Device class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @Todo: Refactor:
 *	Here, one only needs: is_ocl_compiler_opt_disabled(),  params.get_nspace(); get_ntime();
 *	from meta::Inputparameters
 *	Unfortunately, also every single hardware::code class takes it as arg!!
 *	This is needed for opencl_module only, where it is used in collect_build_options(..)
 */

#include "device.hpp"
#include "system.hpp"
#include "openClCode.hpp"
#include "../host_functionality/logger.hpp"

static bool retrieve_device_availability(cl_device_id device_id);
static size_4 calculate_local_lattice_size(size_4 grid_size, const unsigned NSPACE, const unsigned NTIME);
static size_4 calculate_mem_lattice_size(size_4 grid_size, size_4 local_lattice_size, unsigned halo_size);

hardware::Device::Device(cl_context context, cl_device_id device_id, size_4 grid_pos, size_4 grid_size, const meta::Inputparameters& params, bool enable_profiling)
	: DeviceInfo(device_id),
	  context(context), params(params), hardwareParameters(nullptr),
	  profiling_enabled(enable_profiling),
	  profiling_data(),
	  gaugefield_code(nullptr),
	  prng_code(nullptr),
	  real_code(nullptr),
	  complex_code(nullptr),
	  spinor_code(nullptr),
	  spinor_staggered_code(nullptr),
	  fermion_code(nullptr),
	  fermion_staggered_code(nullptr),
	  gaugemomentum_code(nullptr),
	  molecular_dynamics_code(nullptr),
	  correlator_code(nullptr),
	  correlator_staggered_code(nullptr),
	  heatbath_code(nullptr),
	  kappa_code(nullptr),
	  buffer_code(nullptr),
	  grid_pos(grid_pos),
	  grid_size(grid_size),
	  local_lattice_size(calculate_local_lattice_size(grid_size, params.get_nspace(), params.get_ntime())),
	  halo_size(2),
	  mem_lattice_size(calculate_mem_lattice_size(grid_size, local_lattice_size, halo_size)),
	  allocated_bytes(0),
	  max_allocated_bytes(0),
	  allocated_hostptr_bytes(0)
{
	hardwareParameters = new hardware::HardwareParameters( &params );
	logger.debug() << "Initializing " << get_name();
	logger.debug() << "Device position: " << grid_pos;
	logger.debug() << "Local lattice size: " << local_lattice_size;
	logger.debug() << "Memory lattice size: " << mem_lattice_size;

	bool available = retrieve_device_availability(device_id);
	if(!available) {
		logger.error() << "Device is not available!";
	}

	cl_int err;
	logger.debug() << context << ' ' << device_id;
	command_queue = clCreateCommandQueue(context, device_id, profiling_enabled ? CL_QUEUE_PROFILING_ENABLE : 0, &err);
	if(err) {
		throw OpenclException(err, "clCreateCommandQueue", __FILE__, __LINE__);
	}

	logger.trace() << "Initial memory usage (" << grid_pos << "): " << allocated_bytes << " bytes - Maximum usage: " << max_allocated_bytes << " - Host backed memory: " << allocated_hostptr_bytes;
}

hardware::Device::~Device()
{
	if (hardwareParameters)
	{
		delete hardwareParameters;
	}
	if(buffer_code) {
		delete buffer_code;
	}
	if(kappa_code) {
		delete kappa_code;
	}
	if(heatbath_code) {
		delete heatbath_code;
	}
	if(correlator_code) {
		delete correlator_code;
	}
	if(correlator_staggered_code) {
		delete correlator_staggered_code;
	}
	if(gaugemomentum_code) {
		delete gaugemomentum_code;
	}
	if(molecular_dynamics_code) {
		delete molecular_dynamics_code;
	}
	if(fermion_code) {
		delete fermion_code;
	}
	if(spinor_code) {
		delete spinor_code;
	}
	if(fermion_staggered_code) {
		delete fermion_staggered_code;
	}
	if(spinor_staggered_code) {
		delete spinor_staggered_code;
	}
	if(real_code) {
		delete real_code;
	}
	if(complex_code) {
		delete complex_code;
	}
	if(prng_code) {
		delete prng_code;
	}
	if(gaugefield_code) {
		delete gaugefield_code;
	}

	clFinish(command_queue);
	clReleaseCommandQueue(command_queue);

	logger.info() << "Maximum memory used (" << grid_pos << "): " << max_allocated_bytes << " bytes";
}

hardware::Device::operator cl_command_queue() const noexcept
{
	return command_queue;
}
cl_command_queue hardware::Device::get_queue() const noexcept
{
	return command_queue;
}

TmpClKernel hardware::Device::create_kernel(const char * const kernel_name, std::string build_opts) const
{
	if(hardwareParameters->disableOpenCLCompilerOptimizations()) {
		build_opts +=  " -cl-opt-disable";
	}
	return TmpClKernel(kernel_name, build_opts, context, get_id());
}

void hardware::Device::enqueue_kernel(cl_kernel kernel)
{
	enqueue_kernel(kernel, get_preferred_global_thread_num());
}

void hardware::Device::enqueue_kernel(cl_kernel kernel, size_t global_threads)
{
	enqueue_kernel(kernel, global_threads, get_preferred_local_thread_num());
}

void hardware::Device::enqueue_kernel(cl_kernel kernel, size_t global_threads, size_t local_threads)
{
	// setup profiling if required
	cl_event profiling_event;
	// we only want to pass the event if we are actually profiling
	// otherwise the API will write back data into a no longer valid object
	cl_event * const profiling_event_p = profiling_enabled ? &profiling_event : 0;

	if(logger.beDebug() ) {
		logger.trace() << "calling clEnqueueNDRangeKernel...";
		logger.trace() << "global_work_size: " << global_threads;
		logger.trace() << "local_work_size:  " << local_threads;

		size_t bytesInKernelName;
		if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName) == CL_SUCCESS) {
			char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
			if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL) == CL_SUCCESS) {
				logger.trace() << "Kernel: " << kernelName;
			} else {
				logger.error() << "Could not retrieve kernel name";
			}
			delete [] kernelName;
		} else {
			logger.error() << "Could not retrieve length of kernel name";
		}
	}

	// queue kernel
	cl_int clerr = clEnqueueNDRangeKernel(command_queue, kernel, 1, 0, &global_threads, &local_threads, 0, 0, profiling_event_p);

	// check for errors
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clEnqueueNDRangeKernel failed, aborting...";
		logger.fatal() << "Some more information:";
		logger.fatal() << "global_work_size: " << global_threads;
		logger.fatal() << "local_work_size:  " << local_threads;

		size_t bytesInKernelName;
		if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName) == CL_SUCCESS) {
			char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
			if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL) == CL_SUCCESS) {
				logger.fatal() << "Failed kernel: " << kernelName;
			} else {
				logger.error() << "Could not retrieve kernel name";
			}
			delete [] kernelName;
		} else {
			logger.error() << "Could not retrieve length of kernel name";
		}

		throw hardware::OpenclException(clerr, "clEnqueueNDRangeKernel", __FILE__, __LINE__);
	}

	// evaluate profiling if required
	if(profiling_enabled) {
		// collect the data of this kernel invocation
		clerr = clWaitForEvents(1, &profiling_event);
		if(clerr) {
			throw hardware::OpenclException(clerr, "clWaitForEvents", __FILE__, __LINE__);
		}

		profiling_data[kernel] += profiling_event;
	}
}

void hardware::Device::enqueue_marker(cl_event * event) const
{
	cl_int err = clEnqueueMarker(command_queue, event);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueMarker()", __FILE__, __LINE__);
	}
}

void hardware::Device::enqueue_barrier(const hardware::SynchronizationEvent& event) const
{
	cl_event const cl_event = event.raw();
#if CL_VERSION_1_2
	// Only supported on OpenCL 1.2 and up
	cl_int err = clEnqueueBarrierWithWaitList(command_queue, 1, &cl_event, 0);
#else
#pragma message "Using legacy function clEnqueueWaitForEvents, upgrade to OpenCL 1.2 to enable the usage of clEnqueueBarrierWithWaitList()."
	cl_int err = clEnqueueWaitForEvents(command_queue, 1, &cl_event);
#endif
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueBarrier()", __FILE__, __LINE__);
	}
}

void hardware::Device::enqueue_barrier(const hardware::SynchronizationEvent& event1, const hardware::SynchronizationEvent& event2) const
{
	cl_event const cl_events[] = {event1.raw(), event2.raw()};
#if CL_VERSION_1_2
	// Only supported on OpenCL 1.2 and up
	cl_int err = clEnqueueBarrierWithWaitList(command_queue, 2, cl_events, 0);
#else
#pragma message "Using legacy function clEnqueueWaitForEvents, upgrade to OpenCL 1.2 to enable the usage of clEnqueueBarrierWithWaitList()."
	cl_int err = clEnqueueWaitForEvents(command_queue, 2, cl_events);
#endif
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueBarrier()", __FILE__, __LINE__);
	}
}

static int get_alignment_badness(size_t bytes)
{
	// this is pretty generic for all GPUs
	return (bytes % 256) ? 1 : 0;
}

static int get_cypress_stride_badness(size_t bytes, size_t lanes)
{
	const size_t CRITICAL_STRIDE = 16 * 1024; // at this stride performance is worst
	const size_t CRITICAL_STRIDE_RANGE = 768; // width of the critical stride

	int badness = get_alignment_badness(bytes);
	for(size_t hops = 1; hops < lanes; ++hops) {
		size_t dist_to_critical = (bytes * hops) % CRITICAL_STRIDE;

		if(dist_to_critical >= (CRITICAL_STRIDE - CRITICAL_STRIDE_RANGE) || dist_to_critical <= CRITICAL_STRIDE_RANGE) {
			++badness;
		}
	}
	return badness;
}

size_t hardware::Device::recommend_stride(size_t elems, size_t type_size, size_t lane_count) const
{
	const size_t MAX_ADD_STRIDE = 8 * 1024; // never add more than 8 KiB per lane
	const auto name = get_name();
	if(name == std::string("Cypress") || name == std::string("Cayman")) {
		logger.debug() << "Using cypress stride";
		// apply advanced stride rules
		for(size_t stride = elems; stride <= elems + MAX_ADD_STRIDE / type_size; ++stride) {
			if(get_cypress_stride_badness(stride * type_size, lane_count) == 0) {
				logger.debug() << "Return stride of " << stride << " elements, which is " << stride * type_size << " bytes.";
				return stride;
			}
		}
		throw OptimizationError();
	} else {
		// simply align to 256 Bytes
		logger.debug() << "Using default stride";
		for(size_t stride = elems; stride <= elems + MAX_ADD_STRIDE / type_size; ++stride) {
			if(get_alignment_badness(stride * type_size) == 0) {
				return stride;
			}
		}
		throw OptimizationError();
	}
}

bool hardware::Device::is_profiling_enabled() const noexcept
{
	return profiling_enabled;
}

void hardware::Device::flush() const
{
	cl_int err = clFlush(command_queue);
	if(err) {
		throw hardware::OpenclException(err, "Failed when waiting for OpenCL device to finish.", __FILE__, __LINE__);
	}
}

void hardware::Device::synchronize() const
{
	cl_int err = clFinish(command_queue);
	if(err) {
		throw hardware::OpenclException(err, "Failed when waiting for OpenCL device to finish.", __FILE__, __LINE__);
	}
}

hardware::ProfilingData hardware::Device::get_profiling_data(const cl_kernel& kernel) noexcept {
	return profiling_data[kernel];
}

const hardware::code::Gaugefield * hardware::Device::get_gaugefield_code()
{
	if(!gaugefield_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		gaugefield_code = codeBuilder.getCode_gaugefield(this).release();
	}
	return gaugefield_code;
}

const hardware::code::PRNG * hardware::Device::get_prng_code()
{
	if(!prng_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		prng_code = codeBuilder.getCode_PRNG(this).release();
	}
	return prng_code;
}

const hardware::code::Real * hardware::Device::get_real_code()
{
	if(!real_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		real_code = codeBuilder.getCode_real(this).release();
	}
	return real_code;
}

const hardware::code::Complex * hardware::Device::get_complex_code()
{
	if(!complex_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		complex_code = codeBuilder.getCode_complex(this).release();
	}
	return complex_code;
}

const hardware::code::Spinors * hardware::Device::get_spinor_code()
{
	if(!spinor_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		spinor_code = codeBuilder.getCode_Spinors(this).release();
	}
	return spinor_code;
}

const hardware::code::Spinors_staggered * hardware::Device::get_spinor_staggered_code()
{
	if(!spinor_staggered_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		spinor_staggered_code = codeBuilder.getCode_Spinors_staggered(this).release();
	}
	return spinor_staggered_code;
}

const hardware::code::Fermions * hardware::Device::get_fermion_code()
{
	if(!fermion_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		fermion_code = codeBuilder.getCode_Fermions(this).release();
	}
	return fermion_code;
}

const hardware::code::Fermions_staggered * hardware::Device::get_fermion_staggered_code()
{
	if(!fermion_staggered_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		fermion_staggered_code = codeBuilder.getCode_Fermions_staggered(this).release();
	}
	return fermion_staggered_code;
}

const hardware::code::Gaugemomentum * hardware::Device::get_gaugemomentum_code()
{
	if(!gaugemomentum_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		gaugemomentum_code = codeBuilder.getCode_Gaugemomentum(this).release();
	}
	return gaugemomentum_code;
}

const hardware::code::Molecular_Dynamics * hardware::Device::get_molecular_dynamics_code()
{
	if(!molecular_dynamics_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		molecular_dynamics_code = codeBuilder.getCode_Molecular_Dynamics(this).release();
	}
	return molecular_dynamics_code;
}

const hardware::code::Correlator * hardware::Device::get_correlator_code()
{
	if(!correlator_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		correlator_code = codeBuilder.getCode_Correlator(this).release();
	}
	return correlator_code;
}

const hardware::code::Correlator_staggered * hardware::Device::get_correlator_staggered_code()
{
	if(!correlator_staggered_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		correlator_staggered_code = codeBuilder.getCode_Correlator_staggered(this).release();
	}
	return correlator_staggered_code;
}

const hardware::code::Heatbath * hardware::Device::get_heatbath_code()
{
	if(!heatbath_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		heatbath_code = codeBuilder.getCode_Heatbath(this).release();
	}
	return heatbath_code;
}

const hardware::code::Kappa * hardware::Device::get_kappa_code()
{
	if(!kappa_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		kappa_code = codeBuilder.getCode_Kappa(this).release();
	}
	return kappa_code;
}

const hardware::code::Buffer * hardware::Device::get_buffer_code()
{
	if(!buffer_code) {
		hardware::OpenClCode_fromMetaInputparameters codeBuilder( params );
		//todo: do not use release here. real_code itself should rather be a smart pointer
		buffer_code = codeBuilder.getCode_Buffer(this).release();
	}
	return buffer_code;
}

void hardware::print_profiling(Device * device, const std::string& filename, int id)
{
	if(device->kappa_code) {
		device->kappa_code->print_profiling(filename, id);
	}
	if(device->heatbath_code) {
		device->heatbath_code->print_profiling(filename, id);
	}
	if(device->correlator_code) {
		device->correlator_code->print_profiling(filename, id);
	}
	if(device->correlator_staggered_code) {
		device->correlator_staggered_code->print_profiling(filename, id);
	}
	if(device->fermion_code) {
		device->fermion_code->print_profiling(filename, id);
	}
	if(device->fermion_staggered_code) {
		device->fermion_staggered_code->print_profiling(filename, id);
	}
	if(device->spinor_code) {
		device->spinor_code->print_profiling(filename, id);
	}
	if(device->spinor_staggered_code) {
		device->spinor_staggered_code->print_profiling(filename, id);
	}
	if(device->complex_code) {
		device->complex_code->print_profiling(filename, id);
	}
	if(device->prng_code) {
		device->prng_code->print_profiling(filename, id);
	}
	if(device->gaugefield_code) {
		device->gaugefield_code->print_profiling(filename, id);
	}
	if(device->buffer_code) {
		device->buffer_code->print_profiling(filename, id);
	}
	if(device->molecular_dynamics_code) {
		device->molecular_dynamics_code->print_profiling(filename, id);
	}
	if(device->gaugemomentum_code) {
		device->gaugemomentum_code->print_profiling(filename, id);
	}
}

static bool retrieve_device_availability(cl_device_id device_id)
{
	using namespace hardware;
	cl_bool available;
	cl_int err = clGetDeviceInfo(device_id, CL_DEVICE_AVAILABLE, sizeof(cl_bool), &available, 0);
	if(err) {
		throw OpenclException(err, "clGetDeviceInfo(CL_DEVICE_AVAILABLE)", __FILE__, __LINE__);
	}
	return available;
}

size_4 hardware::Device::get_grid_pos() const
{
	return grid_pos;
}

size_4 hardware::Device::get_grid_size() const
{
	return grid_size;
}

static size_4 calculate_local_lattice_size(size_4 grid_size, const unsigned NSPACE, const unsigned NTIME)
{
	const size_4 local_size(NSPACE / grid_size.x, NSPACE / grid_size.y, NSPACE / grid_size.z, NTIME / grid_size.t);

	if(local_size.x * grid_size.x != NSPACE
	   || local_size.y * grid_size.y != NSPACE
	   || local_size.z * grid_size.z != NSPACE
	   || local_size.t * grid_size.t != NTIME) {
		throw std::invalid_argument("The lattice cannot be distributed onto the given grid.");
	}

	if(local_size.x % 2 || local_size.y % 2 || local_size.z % 2 || local_size.t % 2) {
		logger.warn() << "Local lattice size is odd. This is known to cause problems!";
	}

	return local_size;
}

size_4 hardware::Device::get_local_lattice_size() const
{
	return local_lattice_size;
}

size_4 hardware::Device::get_mem_lattice_size() const
{
	return mem_lattice_size;
}

unsigned hardware::Device::get_halo_size() const
{
	return halo_size;
}

static size_4 calculate_mem_lattice_size(size_4 grid_size, size_4 local_lattice_size, unsigned halo_size)
{
	if(local_lattice_size.t < halo_size) {
		throw std::invalid_argument("The lattice cannot be distributed onto the given grid.");
	}
	return size_4(
	         local_lattice_size.x + (grid_size.x > 1 ? 2 * halo_size : 0),
	         local_lattice_size.y + (grid_size.y > 1 ? 2 * halo_size : 0),
	         local_lattice_size.z + (grid_size.z > 1 ? 2 * halo_size : 0),
	         local_lattice_size.t + (grid_size.t > 1 ? 2 * halo_size : 0)
	       );
}

void hardware::Device::markMemReleased(const bool host, const size_t size)
{
	logger.trace() << "Released " << size << " bytes on (" << grid_pos << ").";
	if(host) {
		allocated_hostptr_bytes -= size;
	} else {
		allocated_bytes -= size;
	}
	if(allocated_bytes > max_allocated_bytes)
		max_allocated_bytes = allocated_bytes;
	logger.trace() << "Memory usage (" << grid_pos << "): " << allocated_bytes << " bytes - Maximum usage: " << max_allocated_bytes << " - Host backed memory: " << allocated_hostptr_bytes;
}

void hardware::Device::markMemAllocated(const bool host, const size_t size)
{
	logger.trace() << "Allocted " << size << " bytes on (" << grid_pos << ").";
	if(host) {
		allocated_hostptr_bytes += size;
	} else {
		allocated_bytes += size;
	}
	logger.trace() << "Memory usage (" << grid_pos << "): " << allocated_bytes << " bytes - Maximum usage: " << max_allocated_bytes << " - Host backed memory: " << allocated_hostptr_bytes;
}
