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

#include "matrix6x6.hpp"
#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "../buffers/halo_update.hpp"
#include "../../geometry/parallelization.hpp"

hardware::lattices::Matrix6x6::Matrix6x6(const hardware::System& system):
system(system), buffers(allocate_buffers())
{}

hardware::lattices::Matrix6x6::~Matrix6x6()
{
    release_buffers(&buffers);
}

const std::vector<const hardware::buffers::6x6 *> hardware::lattices::Matrix6x6::get_buffers() const noexcept
{
    return buffers;
}

std::vector<const hardware::buffers::6x6 *> hardware::lattices::Matrix6x6::allocate_buffers()
{
    using hardware::buffers::6x6;
    
    std::vector<const 6x6 *> buffers;
    
    auto const devices = system.get_devices();
    for(auto device: devices)
    {
        buffers.push_back(new 6x6(device->getLocalLatticeMemoryExtents().getLatticeVolume(), device)); //todo: do not calculate here!
    }
    return buffers;
}

void hardware::lattices::Matrix6x6::release_buffers(std::vector<const hardware::buffers::6x6 *>* buffers)
{
    for(auto buffer: *buffers)
    {
        delete buffer;
    }
    buffers->clear();
}

void hardware::lattices::Matrix6x6::send_matrix6x6_to_buffers(const Matrix6x6 * const gf_host)
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
        
        if(buffers.size() == 1) device->getGaugefieldCode()->importGaugefield(buffer, gf_host);
        else{
            Matrix6x6 * mem_host = new Matrix6x6[buffer->get_elements()];
            //				//todo: put these calls into own fct.! With smart pointers?
            memcpy(&mem_host[tmp2.getMainPartIndex_destination()]  , &gf_host[tmp2.getMainPartIndex_source()]  , tmp2.getMainPartSizeInBytes());
            memcpy(&mem_host[tmp2.getFirstHaloIndex_destination()] , &gf_host[tmp2.getFirstHaloPartIndex_source()] , tmp2.getHaloPartSizeInBytes());
            memcpy(&mem_host[tmp2.getSecondHaloIndex_destination()], &gf_host[tmp2.getSecondHaloPartIndex_source()], tmp2.getHaloPartSizeInBytes());
            
            importMatrix6x6Field(buffer, mem_host);
            delete[] mem_host;
        }
        device->synchronize();
    }
    // 	}
    logger.trace() << "import complete";
}

void hardware::lattices::Matrix6x6::fetch_matrix6x6_from_buffers(Matrix6x6 * const gf_host)
{
    logger.trace() << "fetching matrix6x6field";
    if(buffers.size() == 1) {
        auto device = buffers[0]->get_device();
        exportMatrix6x6Field(gf_host, buffers[0]);
        device->synchronize();
    } else {
        for(auto const buffer: buffers) {
            // fetch local part for each device
            auto device = buffer->get_device();
            Matrix6x6 * mem_host = new Matrix6x6[buffer->get_elements()];
            
            exportMatrix6x6Field(mem_host, buffer);
            
            TemporalParallelizationHandlerLink tmp2(device->getGridPos(), device->getLocalLatticeExtents(), sizeof(Matrix6x6), device->getHaloExtent());
            memcpy(&gf_host[tmp2.getMainPartIndex_source()]  , &mem_host[tmp2.getMainPartIndex_destination()]  , tmp2.getMainPartSizeInBytes());
            
            delete[] mem_host;
        }
    }
}

void hardware::lattices::importMatrix6x6Field(const hardware::buffers::6x6 * matrix6x6, const Matrix6x6 * const data) const
{
    using namespace hardware::buffers;
    
    logger.trace() << "Import matrix6x6field to get_device()";
    if(get_device()->get_prefers_soa()) {
        Plain<Matrix6x6> tmp(6x6->get_elements(), get_device());
        tmp.load(data);
        //convertGaugefieldToSOA_device(matrix6x6, &tmp);
        throw Print_Error_Message("No Soa for Matrix6x6 implemented yet");
    } else {
        6x6->load(data);
    }
}

void hardware::lattices::exportMatrix6x6Field(Matrix6x6 * const dest, const hardware::buffers::6x6 * matrix6x6) const
{
    using namespace hardware::buffers;
    
    logger.trace() << "Exporting matrix6x6field from get_device()";
    if(get_device()->get_prefers_soa()) {
        Plain<Matrix6x6> tmp(6x6->get_elements(), get_device());
        //convertGaugefieldFromSOA_device(&tmp, matrix6x6);
        tmp.dump(dest);
        throw Print_Error_Message("No Soa for Matrix6x6 implemented yet");

    } else {
        6x6->dump(dest);
    }
}

void hardware::lattices::Matrix6x6::update_halo_aos(const std::vector<const hardware::buffers::6x6 *> buffers, const hardware::System& system) const
{
    // check all buffers are non-soa
    for(auto const buffer: buffers) {
        if(buffer->is_soa()) {
            throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
        }
    }
    
    hardware::buffers::update_halo<Matrix6x6>(buffers, system, NDIM);
}

void hardware::lattices::Matrix6x6::update_halo_soa(const std::vector<const hardware::buffers::6x6 *> buffers, const hardware::System& system) const
{
    // check all buffers are non-soa
    for(auto const buffer: buffers) {
        if(!buffer->is_soa()) {
            throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
        }
    }
    
    hardware::buffers::update_halo_soa<Matrix6x6>(buffers, system, .5, 2 * NDIM);
}

void hardware::lattices::Matrix6x6::update_halo() const
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


void hardware::lattices::Matrix6x6::smear(unsigned int smearingSteps)
{
    unsmeared_buffers = allocate_buffers();
    
    for(size_t i = 0; i < buffers.size(); ++i) {
        auto buf = buffers[i];
        auto device = buf->get_device();
        auto gf_code = device->getGaugefieldCode(); // should be like: get_gaugefield_code( HardwareParameters_gaugefield )
        
        hardware::buffers::copyData(unsmeared_buffers[i], buf);
        
        int rho_iter = smearingSteps;
        logger.debug() << "\t\tperform " << rho_iter << " steps of stout-smearing to the gaugefield...";
        
        //one needs a temporary gf to apply the smearing to
        const hardware::buffers::6x6 gf_tmp(buf->get_elements(), device);
        for(int i = 0; i < rho_iter - 1; i += 2) {
            gf_code->stout_smear_device(buf, &gf_tmp);
            gf_code->stout_smear_device(&gf_tmp, buf);
        }
        //if rho_iter is odd one has to copy ones more
        if(rho_iter % 2 == 1) {
            gf_code->stout_smear_device(buf, &gf_tmp);
            hardware::buffers::copyData(buf, &gf_tmp);
        }
    }
}

void hardware::lattices::Matrix6x6::unsmear()
{
    if(unsmeared_buffers.size() == 0) {
        logger.warn() << "Tried to unsmear gaugefield that is not smeared.";
        return;
    }
    
    unsmeared_buffers = allocate_buffers();
    
    for(size_t i = 0; i < buffers.size(); ++i) {
        auto buf = buffers[i];
        hardware::buffers::copyData(buf, unsmeared_buffers[i]);
    }
    
    release_buffers(&unsmeared_buffers);
}