/** @file
 * Testcases for the hardware::buffers::Buffer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "buffer.hpp"
#include "plain.hpp"
#include "../system.hpp"
#include "../device.hpp"
#include "../../meta/type_ops.hpp"
#include <stdexcept>

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::Buffer
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;
	using namespace hardware::buffers;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	System system(params);
	const std::vector<Device*>& devices = system.get_devices();
for(Device * device : devices) {

		Buffer dummy(sizeof(float), device);
		BOOST_REQUIRE_EQUAL(dummy.get_bytes(), sizeof(float));
		const cl_mem * tmp = dummy;
		BOOST_REQUIRE(tmp);
		BOOST_REQUIRE(*tmp);
		BOOST_REQUIRE_EQUAL(device->get_id(), dummy.get_device()->get_id());
	}
}

BOOST_AUTO_TEST_CASE(copy)
{
	// normal buffer doesn't have input/output -> only test exception if different size

	using namespace hardware;
	using namespace hardware::buffers;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	System system(params);
	const std::vector<Device*>& devices = system.get_devices();
for(Device * device : devices) {
		Buffer dummy(sizeof(float), device);
		BOOST_REQUIRE_EQUAL(dummy.get_bytes(), sizeof(float));
		Buffer dummy2(2 * dummy.get_bytes(), device);
		Buffer dummy3(2 * dummy.get_bytes(), device);

		BOOST_CHECK_THROW(copyData(&dummy, &dummy2), std::invalid_argument);
		copyData(&dummy2, &dummy3);
		BOOST_CHECK_THROW(copyData(&dummy3, &dummy), std::invalid_argument);
	}
}

BOOST_AUTO_TEST_CASE(clear)
{
	// normal buffer doesn't have input/output -> only test exception if different size

	using namespace hardware;
	using namespace hardware::buffers;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	System system(params);
	const std::vector<Device*>& devices = system.get_devices();
for(Device * device : devices) {
		{
			const size_t SIZE = 1024;
			float host[SIZE];
			fill(host, SIZE);
			Plain<float> dummy(SIZE, device);
			dummy.load(host);
			dummy.clear();
			dummy.dump(host);
			for(int i = 0; i < SIZE; ++i) {
				BOOST_REQUIRE_EQUAL(host[i], 0.f);
			}
		}

		{
			const size_t SIZE = 71;
			char host[SIZE];
			fill(host, SIZE);
			Plain<char> dummy(SIZE, device);
			dummy.load(host);
			dummy.clear();
			dummy.dump(host);
			for(int i = 0; i < SIZE; ++i) {
				BOOST_REQUIRE_EQUAL(host[i], 0);
			}
		}
	}
}
