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

#include "matrix6x6Field.hpp"

#include <cmath>

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
//#include "../buffers/6x6.hpp"
#include "flopUtilities.hpp"
#include "../../geometry/latticeGrid.hpp"

using namespace std;

void hardware::code::matrix6x6Field::fill_kernels()
{
	basic_opencl_code = get_basic_sources() << "operations_geometry.cl" << "operations_complex.cl" << "operations_complex.h" << "types_fermions.h" << "types_hmc.h" << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl" << "operations_su3vec.cl" << "operations_spinor.cl" << "spinorfield.cl" << "operations_gaugemomentum.cl";

	clover_eo_inverse_explizit_upper_left = createKernel("clover_eo_inverse_explizit_upper_left") << basic_opencl_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "operations_matrix6x6.cl" << "fermionmatrix_eo_clover.cl" << "fermionmatrix_eo_clover_explizit.cl" << "fermionmatrix_eo_clover_inverse.cl";
	clover_eo_inverse_explizit_lower_right = createKernel("clover_eo_inverse_explizit_lower_right") << basic_opencl_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "operations_matrix6x6.cl" << "fermionmatrix_eo_clover.cl" << "fermionmatrix_eo_clover_explizit.cl" << "fermionmatrix_eo_clover_inverse.cl";
}

void hardware::code::matrix6x6Field::clear_kernels()
{
	cl_uint clerr = CL_SUCCESS;

	logger.debug() << "Clearing Molecular_Dynamics kernels.." ;

    clerr = clReleaseKernel(clover_eo_inverse_explizit_upper_left);
    clerr = clReleaseKernel(clover_eo_inverse_explizit_lower_right);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}

/*void hardware::code::matrix6x6Field::print_profiling(const std::string& filename, int number) const
{
    Opencl_Module::print_profiling(filename, clover_eo_inverse_explizit_upper_left);
    Opencl_Module::print_profiling(filename, clover_eo_inverse_explizit_lower_right);
    Opencl_Module::print_profiling(filename, convertMatrix6x6FieldToSOA);
    Opencl_Module::print_profiling(filename, convertMatrix6x6FieldFromSOA);
}*/

void hardware::code::matrix6x6Field::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);
}

void hardware::code::matrix6x6Field::clover_eo_inverse_explizit_upper_left_device(const hardware::buffers::matrix6x6 * out, const hardware::buffers::SU3 * gf, hmc_float kappa, hmc_float csw) const
{
	using namespace hardware::buffers;

	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = kernelParameters->getKappa();
	else kappa_tmp = kappa;

	//get csw
	hmc_float csw_tmp;
    if(csw == ARG_DEF) csw_tmp = kernelParameters->getCsw();
    else csw_tmp = csw;

    //query work-sizes for kernel
   	size_t ls2, gs2;
   	cl_uint num_groups;
   	this->get_work_sizes(clover_eo_inverse_explizit_upper_left, &ls2, &gs2, &num_groups);
   	//set arguments
   	int clerr = clSetKernelArg(clover_eo_inverse_explizit_upper_left, 0, sizeof(cl_mem), out->get_cl_buffer());
   	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

   	clerr = clSetKernelArg(clover_eo_inverse_explizit_upper_left, 1, sizeof(cl_mem), gf->get_cl_buffer());
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

   	clerr = clSetKernelArg(clover_eo_inverse_explizit_upper_left, 2, sizeof(hmc_float), &kappa_tmp);
   	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

   	clerr = clSetKernelArg(clover_eo_inverse_explizit_upper_left, 3, sizeof(hmc_float), &csw_tmp);
   	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

   	get_device()->enqueue_kernel( clover_eo_inverse_explizit_upper_left, gs2, ls2);
}

void hardware::code::matrix6x6Field::clover_eo_inverse_explizit_lower_right_device(const hardware::buffers::matrix6x6 * out, const hardware::buffers::SU3 * gf, hmc_float kappa, hmc_float csw) const
{
	using namespace hardware::buffers;

	//get kappa
    hmc_float kappa_tmp;
    if(kappa == ARG_DEF) kappa_tmp = kernelParameters->getKappa();
    else kappa_tmp = kappa;

    //get csw
    hmc_float csw_tmp;
    if(csw == ARG_DEF) csw_tmp = kernelParameters->getCsw();
    else csw_tmp = csw;

    //query work-sizes for kernel
   	size_t ls2, gs2;
   	cl_uint num_groups;
   	this->get_work_sizes(clover_eo_inverse_explizit_lower_right, &ls2, &gs2, &num_groups);
   	//set arguments
   	int clerr = clSetKernelArg(clover_eo_inverse_explizit_lower_right, 0, sizeof(cl_mem), out->get_cl_buffer());
   	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

   	clerr = clSetKernelArg(clover_eo_inverse_explizit_lower_right, 1, sizeof(cl_mem), gf->get_cl_buffer());
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

   	clerr = clSetKernelArg(clover_eo_inverse_explizit_lower_right, 2, sizeof(hmc_float), &kappa_tmp);
   	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

   	clerr = clSetKernelArg(clover_eo_inverse_explizit_lower_right, 3, sizeof(hmc_float), &csw_tmp);
   	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

   	get_device()->enqueue_kernel( clover_eo_inverse_explizit_upper_left, gs2, ls2);
}

uint64_t hardware::code::matrix6x6Field::get_flop_size(const std::string& in) const
{
	return 0;
}

size_t hardware::code::matrix6x6Field::get_read_write_size(const std::string& in) const
{
	return 0;
}

void hardware::code::matrix6x6Field::importMatrix6x6Field(const hardware::buffers::matrix6x6 * matrix6x6Field, const Matrix6x6 * const data) const
{
    using namespace hardware::buffers;

    logger.trace() << "Import matrix6x6Field to get_device()";
    if(get_device()->get_prefers_soa()) {
        Plain<Matrix6x6> tmp(matrix6x6Field->get_elements(), get_device());
        tmp.load(data);
        //convertGaugefieldToSOA_device(gaugefield, &tmp);
    } else {
    	matrix6x6Field->load(data);
    }
}

void hardware::code::matrix6x6Field::exportMatrix6x6Field(Matrix6x6 * const dest, const hardware::buffers::matrix6x6 * matrix6x6Field) const
{
    using namespace hardware::buffers;

    logger.trace() << "Exporting matrix6x6Field from get_device()";
    if(get_device()->get_prefers_soa()) {
        Plain<Matrix6x6> tmp(matrix6x6Field->get_elements(), get_device());
        //convertGaugefieldFromSOA_device(&tmp, gaugefield);
        tmp.dump(dest);
    } else {
    	matrix6x6Field->dump(dest);
    }
}

/*void hardware::code::matrix6x6Field::convertMatrix6x6FieldToSOA_device(const hardware::buffers::matrix6x6 * out, const hardware::buffers::Plain<Matrix6x6> * in) const
{
    if(!out->is_soa()) {
        throw std::invalid_argument("Destination buffer must be a SOA buffer");
    }
    
    size_t ls2, gs2;
    cl_uint num_groups;
    this->get_work_sizes(convertMatrix6x6FieldToSOA, &ls2, &gs2, &num_groups);
    
    //set arguments
    int clerr = clSetKernelArg(convertMatrix6x6FieldToSOA, 0, sizeof(cl_mem), out->get_cl_buffer());
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
    
    clerr = clSetKernelArg(convertMatrix6x6FieldToSOA, 1, sizeof(cl_mem), in->get_cl_buffer());
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
    
    get_device()->enqueue_kernel(convertMatrix6x6FieldToSOA, gs2, ls2);
}

void hardware::code::matrix6x6Field::convertMatrix6x6FieldFromSOA_device(const hardware::buffers::Plain<Matrix6x6> * out, const hardware::buffers::matrix6x6 * in) const
{
    if(!in->is_soa()) {
        throw std::invalid_argument("Source buffer must be a SOA buffer");
    }
    
    size_t ls2, gs2;
    cl_uint num_groups;
    this->get_work_sizes(convertMatrix6x6FieldFromSOA, &ls2, &gs2, &num_groups);
    
    //set arguments
    int clerr = clSetKernelArg(convertMatrix6x6FieldFromSOA, 0, sizeof(cl_mem), out->get_cl_buffer());
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
    
    clerr = clSetKernelArg(convertMatrix6x6FieldFromSOA, 1, sizeof(cl_mem), in->get_cl_buffer());
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
    
    get_device()->enqueue_kernel(convertMatrix6x6FieldFromSOA, gs2, ls2);
}*/

hardware::code::matrix6x6Field::matrix6x6Field(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
: Opencl_Module(kernelParameters, device),
  clover_eo_inverse_explizit_upper_left(0), clover_eo_inverse_explizit_lower_right(0)
{
    fill_kernels();
}

hardware::code::matrix6x6Field::~matrix6x6Field()
{
    clear_kernels();
}

ClSourcePackage hardware::code::matrix6x6Field::get_sources() const noexcept
{
	return basic_opencl_code;
}

