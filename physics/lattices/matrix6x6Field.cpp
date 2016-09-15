/** @file
 * Implementation of the physics::lattices::Gaugefield class
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
//#include "../../host_functionality/host_operations_gaugefield.h"
#include "../../hardware/device.hpp"
#include "../../ildg_io/ildgIo.hpp"
#include "../../hardware/code/matrix6x6Field.hpp"

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

