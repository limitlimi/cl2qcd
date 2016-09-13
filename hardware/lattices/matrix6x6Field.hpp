/*
 * Copyright 2016 Max Theilig
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

#ifndef _HARDWARE_LATTICES_MATRIX6X6_
#define _HARDWARE_LATTICES_MATRIX6X6_

#include "../system.hpp"
#include "../buffers/6x6.hpp"
#include "../buffers/plain.hpp"
#include "../code/matrix6x6Field.hpp"


namespace hardware {
    
    namespace lattices {
        
        class Matrix6x6Field
        {
        public:
            
            virtual ~Matrix6x6Field();
            
            Matrix6x6Field(const hardware::System& system);
            
            const std::vector<const hardware::buffers::matrix6x6 *> get_buffers() const noexcept;
            std::vector<const hardware::buffers::matrix6x6 *> allocate_buffers();
            void release_buffers(std::vector<const hardware::buffers::matrix6x6 *>* buffers);
            void send_matrix6x6_to_buffers(const Matrix6x6 * const gf_host);
            void fetch_matrix6x6_from_buffers( Matrix6x6 * const gf_host);
            
            void update_halo() const;

            
        private:
            hardware::System const& system;
            std::vector<const hardware::buffers::matrix6x6 *> buffers;
            
            void update_halo_soa(std::vector<const hardware::buffers::matrix6x6 *> buffers, const hardware::System& system) const;
            void update_halo_aos(std::vector<const hardware::buffers::matrix6x6 *> buffers, const hardware::System& system) const;

        };
        
    }
    
}

#endif /* _HARDWARE_LATTICES_MATRIX6X6_ */
