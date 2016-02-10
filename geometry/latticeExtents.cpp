/**
 * Copyright 2015 Christopher Pinke
 *           2016 Francesca Cuteri
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

#include "latticeExtents.hpp"

unsigned int LatticeExtents::getNs() const
{
	return ns;
}

unsigned int LatticeExtents::getNt() const
{
	return nt;
}

unsigned int LatticeExtents::getLatticeVolume() const
{
	return ns * ns * ns * nt;
}

unsigned int LatticeExtents::getSpatialLatticeVolume() const
{
	return ns * ns * ns;
}
