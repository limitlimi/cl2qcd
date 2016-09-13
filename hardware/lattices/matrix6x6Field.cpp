/*
 * Copyright 2016 Francesca Cuteri
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

#include "matrix6x6Field.hpp"
#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "../buffers/halo_update.hpp"
#include "../code/matrix6x6Field.hpp"
#include "../../geometry/parallelization.hpp"

hardware::lattices::Matrix6x6Field::Matrix6x6Field(const hardware::System& system):
system(system), buffers(allocate_buffers())
{}

hardware::lattices::Matrix6x6Field::~Matrix6x6Field()
{
    release_buffers(&buffers);
}

/*void hardware::lattices::Matrix6x6Field::set6x6Field(physics::lattices::Gaugefield * gaugefield, physics::lattices::Matrix6x6Field * matrix6x6Field, bool upperOrLower, double kappa, double csw)
{
	auto gaugefieldBuffers = gaugefield->get_buffers();
	auto matrix6x6FieldBuffers = matrix6x6Field->get_buffers();
	size_t num_devs = matrix6x6FieldBuffers.size();

	if(num_devs ==1){
	auto device = matrix6x6FieldBuffers[0]->get_device();
	upperOrLower ?
	device->getMatrix6x6FieldCode()->clover_eo_inverse_explizit_upper_left_device(matrix6x6FieldBuffers[0], gaugefieldBuffers[0], kappa, csw) :
	device->getMatrix6x6FieldCode()->clover_eo_inverse_explizit_lower_right_device(matrix6x6FieldBuffers[0], gaugefieldBuffers[0], kappa, csw);
	}
	else {
		for(size_t i = 0; i < num_devs; ++i) {
			auto device = matrix6x6FieldBuffers[i]->get_device();
			upperOrLower ?
				device->getMatrix6x6FieldCode()->clover_eo_inverse_explizit_upper_left_device(matrix6x6FieldBuffers[i], gaugefieldBuffers[i], kappa, csw) :
				device->getMatrix6x6FieldCode()->clover_eo_inverse_explizit_lower_right_device(matrix6x6FieldBuffers[i], gaugefieldBuffers[i], kappa, csw);
		}
	}
}*/

const std::vector<const hardware::buffers::matrix6x6 *> hardware::lattices::Matrix6x6Field::get_buffers() const noexcept
{
    return buffers;
}

std::vector<const hardware::buffers::matrix6x6 *> hardware::lattices::Matrix6x6Field::allocate_buffers()
{
    using hardware::buffers::matrix6x6;
    
    std::vector<const matrix6x6 *> buffers;
    
    auto const devices = system.get_devices();
    for(auto device: devices)
    {
        buffers.push_back(new matrix6x6(device->getLocalLatticeMemoryExtents().getLatticeVolume(), device)); //todo: do not calculate here!
    }
    return buffers;
}

void hardware::lattices::Matrix6x6Field::release_buffers(std::vector<const hardware::buffers::matrix6x6 *>* buffers)
{
    for(auto buffer: *buffers)
    {
        delete buffer;
    }
    buffers->clear();
}

void hardware::lattices::Matrix6x6Field::send_matrix6x6_to_buffers(const Matrix6x6 * const gf_host)
{
    logger.trace() << "importing matrix6x6field";
    // 	if(buffers.size() == 1) {
    // 		auto device = buffers[0]->get_device();
    // 		device->getGaugefieldCode()->importGaugefield(buffers[0], gf_host);
    // 		device->synchronize();
    // 	} else {
    for(auto const buffer: buffers) {
        auto device = buffer->get_device();
        TemporalParallelizationHandlerLink tmp2(device->getGridPos(), device->getLocalLatticeExtents(), sizeof(Matrix6x6), device->getHaloExtent());
        
        if(buffers.size() == 1) device->getMatrix6x6FieldCode()->importMatrix6x6Field(buffer, gf_host);
        else{
            Matrix6x6 * mem_host = new Matrix6x6[buffer->get_elements()];
            //				//todo: put these calls into own fct.! With smart pointers?
            memcpy(&mem_host[tmp2.getMainPartIndex_destination()]  , &gf_host[tmp2.getMainPartIndex_source()]  , tmp2.getMainPartSizeInBytes());
            memcpy(&mem_host[tmp2.getFirstHaloIndex_destination()] , &gf_host[tmp2.getFirstHaloPartIndex_source()] , tmp2.getHaloPartSizeInBytes());
            memcpy(&mem_host[tmp2.getSecondHaloIndex_destination()], &gf_host[tmp2.getSecondHaloPartIndex_source()], tmp2.getHaloPartSizeInBytes());
            
            device->getMatrix6x6FieldCode()->importMatrix6x6Field(buffer, mem_host);
            delete[] mem_host;
        }
        device->synchronize();
    }
    // 	}
    logger.trace() << "import complete";
}

void hardware::lattices::Matrix6x6Field::fetch_matrix6x6_from_buffers(Matrix6x6 * const gf_host)
{
    logger.trace() << "fetching matrix6x6field";
    if(buffers.size() == 1) {
        auto device = buffers[0]->get_device();
        device->getMatrix6x6FieldCode()->exportMatrix6x6Field(gf_host, buffers[0]);
        device->synchronize();
    } else {
        for(auto const buffer: buffers) {
            // fetch local part for each device
            auto device = buffer->get_device();
            Matrix6x6 * mem_host = new Matrix6x6[buffer->get_elements()];
            
            device->getMatrix6x6FieldCode()->exportMatrix6x6Field(gf_host, buffers[0]);
            
            TemporalParallelizationHandlerLink tmp2(device->getGridPos(), device->getLocalLatticeExtents(), sizeof(Matrix6x6), device->getHaloExtent());
            memcpy(&gf_host[tmp2.getMainPartIndex_source()]  , &mem_host[tmp2.getMainPartIndex_destination()]  , tmp2.getMainPartSizeInBytes());
            
            delete[] mem_host;
        }
    }
}

void hardware::lattices::Matrix6x6Field::update_halo_aos(const std::vector<const hardware::buffers::matrix6x6 *> buffers, const hardware::System& system) const
{
    // check all buffers are non-soa
    for(auto const buffer: buffers) {
        if(buffer->is_soa()) {
            throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
        }
    }
    
    hardware::buffers::update_halo<Matrix6x6>(buffers, system, NDIM);
}

void hardware::lattices::Matrix6x6Field::update_halo_soa(const std::vector<const hardware::buffers::matrix6x6 *> buffers, const hardware::System& system) const
{
    // check all buffers are non-soa
    for(auto const buffer: buffers) {
        if(!buffer->is_soa()) {
            throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
        }
    }
    
    hardware::buffers::update_halo_soa<Matrix6x6>(buffers, system, .5, 2 * NDIM);
}

void hardware::lattices::Matrix6x6Field::update_halo() const
{
    if(buffers.size() > 1) { // for a single device this will be a noop
        // currently either all or none of the buffers must be SOA
        if(buffers[0]->is_soa()) {
            update_halo_soa(buffers, system);
        } else {
            update_halo_aos(buffers, system);
        }
    }
}
