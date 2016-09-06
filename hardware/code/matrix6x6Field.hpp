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

#ifndef _HARDWARE_CODE_MATRIX6X6FIELD_
#define _HARDWARE_CODE_MATRIX6X6FIELD_

#include "opencl_module.hpp"
#inlcude "../buffers/6x6.cpp"
#include "../buffers/plain.hpp"

namespace hardware {

namespace code {

/**
* An OpenCL device
*
* This class wraps all matrix6x6 operations on a device. Each kernel
* has it's own wrapper function.
*
* @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
*/
    
class matrix6x6Field : public OpenCL_Module {
   
public:
    friend hardware::Device;
    
    virtual ~matrix6x6Field
    
    /**
     * Import the matrix6x6Field data into the OpenCL buffer using the device
     * specific storage format.
     *
     * @param[out] matrix6x6Field The OpenCL buffer to writ the gaugefield data to in the device specific format
     * @param[in]  data       The matrix6x6Field data to import into the OpenCL buffer
     *
     * @todo should not be public
     */
    void importMatrix6x6Field(const hardware::buffers::matrix6x6 * matrix6x6Field, const Matrix6x6 * const data) const;
    
    /**
     * Export the matrix6x6Field from the OpenCL buffer, that uses a device
     * specific storage format, into the given pointer using the generic
     * storage format.
     *
     * @param[out] dest The array to store the matrix6x6Field in
     */
    void exportMatrix6x6Field(Matrix6x6 * const dest, const hardware::buffers::matrix6x6 * matrix6x6Field) const;
    
private:
    //cl_kernel convertMatrix6x6FieldToSOA;
    //cl_kernel convertMatrix6x6FieldFromSOA;
    
    //void convertMatrix6x6FieldToSOA_device(const hardware::buffers::matrix6x6 * out, const hardware::buffers::Plain<Matrix6x6> * in) const;
    //void convertMatrix6x6FieldFromSOA_device(const hardware::buffers::Plain<Matrix6x6> * out, const hardware::buffers::matrix6x6 * in) const;
    
    /**
     * Collect the kernels for OpenCL.
     */
    void fill_kernels();
    /**
     * Clear out the kernels,
     */
    void clear_kernels();
};
    
}
    
}

#endif /* _HARDWARE_CODE_MATRIX6X6FIELD_ */