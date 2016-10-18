/** @file
 * Implementation of the physics::lattices::Matrix6x6Field class
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 *
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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
#include "../../hardware/device.hpp"
#include "../../ildg_io/ildgIo.hpp"
#include "../../hardware/code/matrix6x6Field.hpp"
#include "../../hardware/code/Matrix6x6FieldTester.cpp"

physics::lattices::Matrix6x6Field::Matrix6x6Field(const hardware::System& system, const GaugefieldParametersInterface * parameters)
  : system(system), latticeObjectParameters(parameters), matrix6x6Field(system)
{
}

void physics::lattices::Matrix6x6Field::setField(const physics::lattices::Gaugefield * gaugefield, const bool upperOrLower)
{
	if(latticeObjectParameters->getFermact() != common::action::clover)
	{
		throw Invalid_Parameters("The setField method in matrix6x6Field is not to be used with wilson or tm fermion action!", "clover", latticeObjectParameters->getFermact());
	}
	//logger.fatal() << "Matrix6x6Field filled(0=upperleft, 1=lowerright)" << upperOrLower;
	auto gaugefieldBuffers = gaugefield->get_buffers();
	auto matrix6x6FieldBuffers = this->get_buffers();
	size_t num_devs = matrix6x6FieldBuffers.size();

	if(num_devs ==1){
	auto device = matrix6x6FieldBuffers[0]->get_device();
	upperOrLower ?
	device->getMatrix6x6FieldCode()->clover_eo_inverse_explicit_upper_left_device(matrix6x6FieldBuffers[0], gaugefieldBuffers[0], latticeObjectParameters->getKappa(), latticeObjectParameters->getCsw()) :
	device->getMatrix6x6FieldCode()->clover_eo_inverse_explicit_lower_right_device(matrix6x6FieldBuffers[0], gaugefieldBuffers[0], latticeObjectParameters->getKappa(), latticeObjectParameters->getCsw());
	}
	else {
		for(size_t i = 0; i < num_devs; ++i) {
			auto device = matrix6x6FieldBuffers[i]->get_device();
			upperOrLower ?
				device->getMatrix6x6FieldCode()->clover_eo_inverse_explicit_upper_left_device(matrix6x6FieldBuffers[i], gaugefieldBuffers[i], latticeObjectParameters->getKappa(), latticeObjectParameters->getCsw()) :
				device->getMatrix6x6FieldCode()->clover_eo_inverse_explicit_lower_right_device(matrix6x6FieldBuffers[i], gaugefieldBuffers[i], latticeObjectParameters->getKappa(), latticeObjectParameters->getCsw());
		}
	}
}

physics::lattices::Matrix6x6Field::~Matrix6x6Field()
{}

const std::vector<const hardware::buffers::matrix6x6 *> physics::lattices::Matrix6x6Field::get_buffers() const noexcept
{
	return matrix6x6Field.get_buffers();
}

hmc_float physics::lattices::count_Matrix6x6Field(const Matrix6x6Field& field)
{
	hmc_float res = 0;
	auto field_buffers = field.get_buffers();
	size_t num_buffers = field_buffers.size();

	for(size_t i = 0; i < num_buffers; ++i) {
		auto field_buf = field_buffers[i];
		Matrix6x6 * matrix6x6_in = new Matrix6x6[field_buf->get_elements()];
		field_buf->dump(matrix6x6_in);
		res += count_matrix6x6Field(matrix6x6_in, field_buf->get_elements());
	}
	return res;
}

//calculation of clover action S_det = -log(det[(1+T_ee)^2]), see hep-lat/9603008 for details
hmc_float physics::lattices::S_det(const Gaugefield& field, const hmc_float kappa, const hmc_float csw)
{
	const Scalar<hmc_float> res(*field.getSystem());
	S_det(&res, field, kappa, csw);
	return res.get();
}

void physics::lattices::S_det(const Scalar<hmc_float>* res, const Gaugefield& field, const hmc_float kappa, const hmc_float csw)
{
	auto field_buffers = field.get_buffers();
	auto res_buffers = res->get_buffers();
	size_t num_buffers = field_buffers.size();

	// TODO implemente for more than one device
	if(num_buffers != res_buffers.size()) {
		throw std::invalid_argument("The given lattices do not sue the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto field_buf = field_buffers[i];
		auto res_buf = res_buffers[i];
		auto device = field_buf->get_device();
		auto mat6x6_code = device->getMatrix6x6FieldCode();

		mat6x6_code->S_det_device(field_buf, res_buf, kappa, csw);
	}

	res->sum();
}



void physics::lattices::Matrix6x6Field::update_halo() const
{
	matrix6x6Field.update_halo();
}

const hardware::System * physics::lattices::Matrix6x6Field::getSystem() const
{
	return &system;
}

const physics::lattices::GaugefieldParametersInterface * physics::lattices::Matrix6x6Field::getParameters() const
{
	return latticeObjectParameters;
}

