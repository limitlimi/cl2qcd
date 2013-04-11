/** @file
 * Implementation of the hardware::buffers::Spinor class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "spinor.hpp"
#include "../device.hpp"
#include "plain.hpp"
#include "../code/spinors.hpp"

#include <stdexcept>

typedef hmc_complex soa_storage_t;
const size_t soa_storage_lanes = 12;

static size_t calculate_spinor_buffer_size(size_t elems, hardware::Device * device);

hardware::buffers::Spinor::Spinor(size_t elems, hardware::Device * device)
	: Buffer(calculate_spinor_buffer_size(elems, device), device),
	  elems(elems),
	  soa(check_Spinor_for_SOA(device))
{
	// nothing to do
}

size_t hardware::buffers::check_Spinor_for_SOA(hardware::Device * device)
{
	return device->get_prefers_soa();
}

static size_t calculate_spinor_buffer_size(size_t elems, hardware::Device * device)
{
	using namespace hardware::buffers;
	if(check_Spinor_for_SOA(device)) {
		size_t stride = get_Spinor_buffer_stride(elems, device);
		return stride * soa_storage_lanes * sizeof(soa_storage_t);
	} else {
		return elems * sizeof(spinor);
	}
}

size_t hardware::buffers::get_Spinor_buffer_stride(size_t elems, Device * device)
{
	return device->recommend_stride(elems, sizeof(soa_storage_t), soa_storage_lanes);
}

size_t hardware::buffers::Spinor::get_elements() const noexcept
{
	return elems;
}

bool hardware::buffers::Spinor::is_soa() const noexcept
{
	return soa;
}

void hardware::buffers::Spinor::load(const spinor * ptr, size_t elems, size_t offset) const
{
	if(is_soa()) {
		auto device = get_device();
		Plain<spinor> plain(get_elements(), device);
		plain.load(ptr, elems * sizeof(spinor), offset * sizeof(spinor));
		device->get_spinor_code()->convertSpinorfieldToSOA_eo_device(this, &plain);
		device->synchronize();
	} else {
		Buffer::load(ptr, elems * sizeof(spinor), offset * sizeof(spinor));
	}
}

void hardware::buffers::Spinor::dump(spinor * ptr, size_t elems, size_t offset) const
{
	if(is_soa()) {
		auto device = get_device();
		Plain<spinor> plain(get_elements(), device);
		device->get_spinor_code()->convertSpinorfieldFromSOA_eo_device(&plain, this);
		plain.dump(ptr, elems * sizeof(spinor), offset * sizeof(spinor));
	} else {
		Buffer::dump(ptr, elems * sizeof(spinor), offset * sizeof(spinor));
	}
}

void hardware::buffers::Spinor::load_raw(const void * ptr, size_t bytes, size_t offset) const
{
	logger.trace() << "Loading raw data into Spinor buffer.";
	Buffer::load(ptr, bytes, offset);
}

void hardware::buffers::Spinor::dump_raw(void * ptr, size_t bytes, size_t offset) const
{
	logger.trace() << "Dumping raw data from Spinor buffer.";
	Buffer::dump(ptr, bytes, offset);
}

size_t hardware::buffers::Spinor::get_storage_type_size() const noexcept
{
	return soa ? sizeof(soa_storage_t) : sizeof(spinor);
}

size_t hardware::buffers::Spinor::get_lane_stride() const noexcept
{
	return soa ? (get_bytes() / sizeof(soa_storage_t) / soa_storage_lanes) : 0;
}

size_t hardware::buffers::Spinor::get_lane_count() const noexcept
{
	return soa ? soa_storage_lanes : 1;
}
