/** @file
* Implementation of the hardware::buffers::6x6 class
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

#include "6x6.hpp"
#include "../device.hpp"
#include "../code/matrix6x6Field.hpp"

#include <stdexcept>

typedef hmc_complex soa_storage_t;
const size_t soa_storage_lanes = 36;

static size_t calculate_6x6_buffer_size(size_t elems, const hardware::Device * device);

size_t hardware::buffers::calculate6x6FieldSize(const LatticeExtents latticeExtentsIn) noexcept
{
    return 	latticeExtentsIn.getLatticeVolume();
}

hardware::buffers::matrix6x6::matrix6x6(const size_t elems, const hardware::Device * device)
: Buffer(calculate_6x6_buffer_size(elems, device), device),
elems(elems),
soa(check_6x6_for_SOA(device))
{
    // nothing to do
}

hardware::buffers::matrix6x6::matrix6x6(const LatticeExtents lE, const hardware::Device * device)
: Buffer(calculate_6x6_buffer_size(calculate6x6FieldSize(lE), device), device),
elems(calculate6x6FieldSize(lE)),
soa(check_6x6_for_SOA(device))
{
    // nothing to do
}

size_t hardware::buffers::check_6x6_for_SOA(const hardware::Device * device)
{
    return device->get_prefers_soa();
}

static size_t calculate_6x6_buffer_size(const size_t elems, const hardware::Device * device)
{
    using namespace hardware::buffers;
    if(check_6x6_for_SOA(device)) {
        size_t stride = get_6x6_buffer_stride(elems, device);
        return stride * soa_storage_lanes * sizeof(soa_storage_t);
    } else {
        return elems * sizeof(::Matrix6x6);
    }
}

size_t hardware::buffers::get_6x6_buffer_stride(const size_t elems, const Device * device)
{
    return device->recommendStride(elems, sizeof(soa_storage_t), soa_storage_lanes);
}

size_t hardware::buffers::matrix6x6::get_elements() const noexcept
{
    return elems;
}

bool hardware::buffers::matrix6x6::is_soa() const noexcept
{
    return soa;
}

void hardware::buffers::matrix6x6::load(const Matrix6x6 * ptr, const size_t elems, const size_t offset) const
{
    if(is_soa()) {
        throw std::logic_error("Data cannot be loaded into SOA buffers.");
    } else {
        Buffer::load(ptr, elems * sizeof(::Matrix6x6), offset * sizeof(::Matrix6x6));
    }
}

void hardware::buffers::matrix6x6::dump(Matrix6x6 * ptr, const size_t elems, const size_t offset) const
{
    if(is_soa()) {
        throw std::logic_error("Data cannot be dumped from SOA buffers.");
    } else {
        Buffer::dump(ptr, elems * sizeof(::Matrix6x6), offset * sizeof(::Matrix6x6));
    }
}

size_t hardware::buffers::matrix6x6::get_storage_type_size() const noexcept
{
    return soa ? sizeof(soa_storage_t) : sizeof(::Matrix6x6);
}

size_t hardware::buffers::matrix6x6::get_lane_stride() const noexcept
{
    return soa ? (get_bytes() / sizeof(soa_storage_t) / soa_storage_lanes) : 0;
}

size_t hardware::buffers::matrix6x6::get_lane_count() const noexcept
{
    return soa ? soa_storage_lanes : 1;
}

