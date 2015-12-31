/** @file
 * Testcases for the hardware::buffers::Gaugemomentum class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "gaugemomentum.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::Gaugemomentum
#include <boost/test/unit_test.hpp>

#include "../system.hpp"
#include "../code/mockups.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;
	using namespace hardware::buffers;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	for(Device * device : system.get_devices())
	{
		Gaugemomentum dummy(system.getHardwareParameters()->getLatticeVolume(), device);
		const cl_mem * tmp = dummy;
		BOOST_CHECK(tmp);
		BOOST_CHECK(*tmp);
	}
}

#include "../../meta/type_ops.hpp"
BOOST_AUTO_TEST_CASE(import_export)
{
	using namespace hardware;
	using namespace hardware::buffers;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	const size_t elems = system.getHardwareParameters()->getLatticeVolume() / 2;
	for(Device * device : system.get_devices())
	{
		ae* buf = new ae[elems];
		ae* buf2 = new ae[elems];
		Gaugemomentum dummy(elems, device);
		if(dummy.is_soa()) {
			BOOST_CHECK_THROW(dummy.load(buf), std::logic_error);
			BOOST_CHECK_THROW(dummy.dump(buf), std::logic_error);
		} else {
			fill(buf, elems, 1);
			fill(buf2, elems, 2);
			dummy.load(buf);
			dummy.dump(buf2);
			BOOST_CHECK_EQUAL_COLLECTIONS(buf, buf + elems, buf2, buf2 + elems);
		}
		delete[] buf;
		delete[] buf2;
	}
}


BOOST_AUTO_TEST_CASE(copy)
{
	using namespace hardware;
	using namespace hardware::buffers;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	const size_t elems = system.getHardwareParameters()->getLatticeVolume();
	for(Device * device : system.get_devices())
	{
		if(!check_Gaugemomentum_for_SOA(device))
		{
			ae* buf = new ae[elems];
			ae* buf2 = new ae[elems];
			Gaugemomentum dummy(elems, device);
			Gaugemomentum dummy2(elems, device);

			fill(buf, elems, 1);
			fill(buf2, elems, 2);
			dummy.load(buf);
			copyData(&dummy2, &dummy);
			dummy2.dump(buf2);
			BOOST_CHECK_EQUAL_COLLECTIONS(buf, buf + elems, buf2, buf2 + elems);

			delete[] buf;
			delete[] buf2;
		}
	}
}
