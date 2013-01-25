#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"
#include "../host_random.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_SPINORS
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

class TestGaugefield : public Gaugefield_hybrid {

public:
	TestGaugefield(const hardware::System * system) : Gaugefield_hybrid(system), prng(*system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, prng);
		meta::print_info_hmc("test program", inputfile);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	hardware::code::Spinors * get_device();

private:
	physics::PRNG prng;
};

void TestGaugefield::init_tasks()
{
	opencl_modules = new hardware::code::Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = get_device_for_task(0)->get_spinor_code();
}

void TestGaugefield::finalize_opencl()
{
	Gaugefield_hybrid::finalize_opencl();
}

void fill_sf_with_one(spinor * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0 = hmc_complex_one;
		sf_in[i].e0.e1 = hmc_complex_one;
		sf_in[i].e0.e2 = hmc_complex_one;
		sf_in[i].e1.e0 = hmc_complex_one;
		sf_in[i].e1.e1 = hmc_complex_one;
		sf_in[i].e1.e2 = hmc_complex_one;
		sf_in[i].e2.e0 = hmc_complex_one;
		sf_in[i].e2.e1 = hmc_complex_one;
		sf_in[i].e2.e2 = hmc_complex_one;
		sf_in[i].e3.e0 = hmc_complex_one;
		sf_in[i].e3.e1 = hmc_complex_one;
		sf_in[i].e3.e2 = hmc_complex_one;
	}
	return;
}

void fill_sf_with_random(spinor * sf_in, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0.re = prng_double();
		sf_in[i].e0.e1.re = prng_double();
		sf_in[i].e0.e2.re = prng_double();
		sf_in[i].e1.e0.re = prng_double();
		sf_in[i].e1.e1.re = prng_double();
		sf_in[i].e1.e2.re = prng_double();
		sf_in[i].e2.e0.re = prng_double();
		sf_in[i].e2.e1.re = prng_double();
		sf_in[i].e2.e2.re = prng_double();
		sf_in[i].e3.e0.re = prng_double();
		sf_in[i].e3.e1.re = prng_double();
		sf_in[i].e3.e2.re = prng_double();

		sf_in[i].e0.e0.im = prng_double();
		sf_in[i].e0.e1.im = prng_double();
		sf_in[i].e0.e2.im = prng_double();
		sf_in[i].e1.e0.im = prng_double();
		sf_in[i].e1.e1.im = prng_double();
		sf_in[i].e1.e2.im = prng_double();
		sf_in[i].e2.e0.im = prng_double();
		sf_in[i].e2.e1.im = prng_double();
		sf_in[i].e2.e2.im = prng_double();
		sf_in[i].e3.e0.im = prng_double();
		sf_in[i].e3.e1.im = prng_double();
		sf_in[i].e3.e2.im = prng_double();
	}
	return;
}

void fill_sf_with_random(spinor * sf_in, int size)
{
	fill_sf_with_random(sf_in, size, 123456);
}

hardware::code::Spinors* TestGaugefield::get_device()
{
	return static_cast<hardware::code::Spinors*>(opencl_modules[0]);
}

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_spinors";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	logger.info() << "Finalize device";
	cpu.finalize();
	BOOST_MESSAGE("Test done");
}

void test_sf_cold(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "set_spinorfield_cold";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gf_code = device->get_device()->get_gaugefield_code();

	logger.info() << "Run kernel";
	device->set_spinorfield_cold_device(&in);
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_cold_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "set_spinorfield_cold_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gf_code = device->get_device()->get_gaugefield_code();

	logger.info() << "Run kernel";
	device->set_eoprec_spinorfield_cold_device(&in);
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
	test_build("/opencl_module_spinors_build_input_1");
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("/opencl_module_spinors_build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_COLD)

BOOST_AUTO_TEST_CASE( SF_COLD_1 )
{
	test_sf_cold("/sf_cold_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_COLD_EO)

BOOST_AUTO_TEST_CASE( SF_COLD_EO_1 )
{
	test_sf_cold("/sf_cold_eo_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

