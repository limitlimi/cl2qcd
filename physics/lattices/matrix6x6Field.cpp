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

physics::lattices::Matrix6x6Field::Matrix6x6Field(const hardware::System& system, const Matrix6x6FieldParametersInterface * parameters, bool upperOrLower)
  : system(system), latticeObjectParameters(parameters), matrix6x6Field(system)
{
}

void physics::lattices::Matrix6x6Field::set6x6Field(physics::lattices::Gaugefield * gaugefield, bool upperOrLower)
{
	//matrix6x6Field.set6x6Field(gaugefield, upperOrLower, latticeObjectParameters->getKappa(), latticeObjectParameters->getCsw());
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

const physics::lattices::Matrix6x6FieldParametersInterface * physics::lattices::Matrix6x6Field::getParameters() const
{
	return latticeObjectParameters;
}

