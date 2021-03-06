# Build specifications for the contents of the buffers package
#
# Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
# 	        2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
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

add_library(buffers
	buffer.cpp
	prng_buffer.cpp
	su3.cpp
	spinor.cpp
	gaugemomentum.cpp
	3x3.cpp
	halo_update.cpp
	su3vec.cpp
)

target_link_libraries(buffers
	hardware
	code
)

#
# Definition of tests
#

add_unit_test_withLibraries(hardware::buffers::Buffer buffer_test buffers crypto)
add_unit_test_withLibraries(hardware::buffers::PRNGBuffer prng_buffer_test buffers crypto)
add_unit_test_withLibraries(hardware::buffers::SU3 su3_test buffers crypto)
add_unit_test_withLibraries(hardware::buffers::Plain plain_test buffers crypto)
add_unit_test_withLibraries(hardware::buffers::Spinor spinor_test buffers crypto)
add_unit_test_withLibraries(hardware::buffers::Gaugemomentum gaugemomentum_test buffers crypto)
add_unit_test_withLibraries(hardware::buffers::Matrix3x3 3x3_test buffers crypto)
add_unit_test_withLibraries(hardware::buffers::SU3vec su3vec_test buffers crypto)
