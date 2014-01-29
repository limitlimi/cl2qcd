# Build specifications for the contents of the physics package
#
# Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
#
# This file is part of CL2QCD.
#
# CL2QCD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CL2QCD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.

#
# Collect all code in this part of the source tree
#

add_library(physics
	prng.cpp
	sources.cpp
	gaugeObservables.cpp
	../host_functionality/host_random.cpp
)
add_modules(physics
	lattices
	fermionmatrix
	algorithms
)
target_link_libraries(physics
	meta
	hardware
)

#
# Definition of tets
#

add_unit_test(physics::PRNG prng_test)
add_unit_test(physics::<sources> sources_test)