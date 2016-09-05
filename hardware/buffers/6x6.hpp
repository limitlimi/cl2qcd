/** @file
 * Declaration of the hardware::buffers::6x6 class
 *
 * Copyright (c) 2016 Max Theilig
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

#ifndef _HARDWARE_BUFFERS_6X6_
#define _HARDWARE_BUFFERS_6x6_

#include "buffer.hpp"
#include "../../common_header_files/types.h"
#include "../../geometry/latticeExtents.hpp"

namespace hardware {
namespace buffers {
    size_t calculate6x6FieldSize(LatticeExtents latticeExtentsIn) noexcept;
    
    /**
     * Check whether 6x6 should be stored SOA style on this device
     */
    size_t check_6x6_for_SOA(const hardware::Device * device);
    
    /**
     * Get the stride for an 6x6 buffer of the given number of elements on the given device
     */
    size_t get_6x6_buffer_stride(const size_t elems, const Device * device);
    
    /*
     * A buffer storing 6x6s on the device
     */
    class matrix6x6 : public Buffer {
        
    public:
        /**
         * Allocate a buffer with the default number of
         * elemets for this device.
         *
         * \param elems The size of the buffer in elements
         * \param device The device to locate the buffer on
         */
        matrix6x6(const size_t elems, const Device * device);
        matrix6x6(const LatticeExtents lE, const Device * device);
        
        /*
         * Matrix6x6 buffers cannot be copied
         */
        matrix6x6& operator=(const matrix6x6&) = delete;
        matrix6x6(const matrix6x6&) = delete;
        matrix6x6() = delete;
        
        /**
         * Load data from the given pointer into the buffer.
         *
         * This only works for AoS-Buffers. If the buffer is a SoA buffer
         * an std::logic_error will be thrown.
         *
         * \param elems Allows to limit the number of elements loaded from the given pointer
         * \param offset Allows to store the elements at the given offset into the buffer
         */
        void load(const Matrix6x6 *, size_t elems = 0, size_t offset = 0) const;
        
        /**
         * Store data from the buffer into the given pointer.
         *
         * This only works for AoS-Buffers. If the buffer is a SoA buffer
         * an std::logic_error will be thrown.
         *
         * \param elems Allows to limit the number of elements dumped to the given pointer
         * \param offset Allows to read the elements at the given offset into the buffer
         */
        void dump(Matrix6x6 *, const size_t elems = 0, const size_t offset = 0) const;
        
        /**
         * Get the size of the buffer in elements
         */
        size_t get_elements() const noexcept;
        
        /**
         * Check whether this Buffer uses soa layout
         */
        bool is_soa() const noexcept;
        
        /**
         * Get the size of the type used for storage.
         */
        size_t get_storage_type_size() const noexcept;
        
        /**
         * Get the stride between two lanes (in elements).
         *
         * 0 if not a SOA buffer
         */
        size_t get_lane_stride() const noexcept;
        
        /**
         * Get the number of lanes.
         *
         * 1 if not a SOA buffer
         */
        size_t get_lane_count() const noexcept;
        
    private:
        
        /**
         * The size of the buffer in bytes.
         */
        const size_t elems;
        
        /**
         * Whether the data is stored in a soa fashion
         */
        const bool soa;
    };
}
}

#endif /* _HARDWARE_BUFFERS_6X6_ */
