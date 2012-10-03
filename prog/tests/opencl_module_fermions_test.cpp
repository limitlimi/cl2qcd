#include "../opencl_module_fermions.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

extern std::string const version;
std::string const version = "0.1";

class TestGaugefield : public Gaugefield_hybrid {

public:
	TestGaugefield(const hardware::System * system) : Gaugefield_hybrid(system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU);
		meta::print_info_hmc("test program", inputfile);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	cl_mem in, out;
  cl_mem in_eo_even, in_eo_odd, out_eo;
	cl_mem sqnorm;
  Opencl_Module_Fermions * get_device();
private:
	void fill_buffers();
	void clear_buffers();
	spinor * sf_in;
	spinor * sf_out;
};

void TestGaugefield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Opencl_Module_Fermions(get_parameters(), get_device_for_task(0));
	opencl_modules[0]->init();

	fill_buffers();
}

void TestGaugefield::finalize_opencl()
{
	clear_buffers();
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

void fill_sf_with_random(spinor * sf_in, int size)
{
	prng_init(123456);
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

void TestGaugefield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	cl_int err;
	cl_context context = opencl_modules[0]->get_context();
	size_t NUM_ELEMENTS_SF =  meta::get_spinorfieldsize(get_parameters());

	Opencl_Module_Fermions * device = this->get_device();

	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_out = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters().get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	size_t sf_buf_size = NUM_ELEMENTS_SF * sizeof(spinor);
	//create buffer for sf on device (and copy sf_in to both for convenience)
	in = clCreateBuffer(context, CL_MEM_READ_WRITE, sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(device->get_queue(), in, CL_TRUE, 0, sf_buf_size, sf_in, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	out = clCreateBuffer(context, CL_MEM_READ_WRITE, sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(device->get_queue(), out, CL_TRUE, 0, sf_buf_size, sf_in, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	if(get_parameters().get_use_eo() ) {
	  size_t eo_buf_size = device->get_eoprec_spinorfield_buffer_size();
	  in_eo_even = clCreateBuffer(context, CL_MEM_READ_WRITE, eo_buf_size, 0, &err );
	  BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	  in_eo_odd = clCreateBuffer(context, CL_MEM_READ_WRITE, eo_buf_size, 0, &err );
	  BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	  out_eo = clCreateBuffer(context, CL_MEM_READ_WRITE, eo_buf_size, 0, &err );
	  BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	  device->convert_to_eoprec_device(in_eo_even, in_eo_odd,in);
	} else {
	  in_eo_even = 0;
	  in_eo_odd = 0;
	  out_eo = 0;
	}

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
}

void TestGaugefield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	clReleaseMemObject(in);
	clReleaseMemObject(out);
	clReleaseMemObject(in_eo_even);
	clReleaseMemObject(in_eo_odd);
	clReleaseMemObject(out_eo);
	clReleaseMemObject(sqnorm);

	delete[] sf_in;
	delete[] sf_out;
}

Opencl_Module_Fermions* TestGaugefield::get_device()
{
	return static_cast<Opencl_Module_Fermions*>(opencl_modules[0]);
}

void test_m_tm_plus(std::string inputfile)
{
  std::string kernelName = "m_tm_plus";
  printKernelInfo(kernelName);
  logger.info() << "Init device";
  meta::Inputparameters params = create_parameters(inputfile);
  hardware::System system(params);
  TestGaugefield cpu(&system);
  cl_int err = CL_SUCCESS;
  Opencl_Module_Fermions * device = cpu.get_device();
  logger.info() << "|phi|^2:";
  hmc_float cpu_back;
  device->set_float_to_global_squarenorm_device(cpu.in, cpu.sqnorm);
  err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back, 0, 0, 0);
  BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
  logger.info() << cpu_back;
  logger.info() << "Run kernel";
  device->M_tm_plus_device(cpu.in, cpu.out,  device->get_gaugefield(), params.get_kappa(), meta::get_mubar(params));
  logger.info() << "result:";
  hmc_float cpu_res;
  device->set_float_to_global_squarenorm_device(cpu.out, cpu.sqnorm);
  err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
  BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
  logger.info() << cpu_res;
  logger.info() << "Finalize device";
  cpu.finalize();
  
  logger.info() << "Choosing reference value and acceptance precision";
  hmc_float ref_val = params.get_test_ref_value();
  logger.info() << "reference value:\t" << ref_val;
  hmc_float prec = params.get_solver_prec();
  logger.info() << "acceptance precision: " << prec;
  
  BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);
  BOOST_MESSAGE("Test done");
}

void test_dslash_eo(std::string inputfile){
  std::string kernelName = "dslash_eo";
  printKernelInfo(kernelName);
  logger.info() << "Init device";
  meta::Inputparameters params = create_parameters(inputfile);
  hardware::System system(params);
  TestGaugefield cpu(&system);
  cl_int err = CL_SUCCESS;
  Opencl_Module_Fermions * device = cpu.get_device();
  logger.info() << "|phi|^2:";
  hmc_float cpu_back;
  device->set_float_to_global_squarenorm_eoprec_device(cpu.in_eo_even, cpu.sqnorm);
  err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back, 0, 0, 0);
  BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
  logger.info() << cpu_back;
  
  //switch according to "use_pointsource"
  hmc_float cpu_res;
  if(params.get_use_pointsource()) {
    device->dslash_eo_device( cpu.in_eo_even, cpu.out_eo, device->get_gaugefield(), EVEN, params.get_kappa() );
  } else {
    device->dslash_eo_device( cpu.in_eo_even, cpu.out_eo, device->get_gaugefield(), ODD, params.get_kappa() );
  }
  device->set_float_to_global_squarenorm_eoprec_device(cpu.out_eo, cpu.sqnorm);
  err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
  BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
  logger.info() << "result:";
  logger.info() << cpu_res;
  logger.info() << "Finalize device";
  cpu.finalize();

  logger.info() << "Choosing reference value and acceptance precision";
  hmc_float ref_val = params.get_test_ref_value();
  logger.info() << "reference value:\t" << ref_val;
  hmc_float prec = params.get_solver_prec();
  logger.info() << "acceptance precision: " << prec;
  
  BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);
  BOOST_MESSAGE("Test done");
}


BOOST_AUTO_TEST_SUITE( M_WILSON ) 

BOOST_AUTO_TEST_CASE( M_WILSON_1){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_MINUS  ) 

BOOST_AUTO_TEST_CASE( M_TM_MINUS_1 ){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_PLUS ) 

BOOST_AUTO_TEST_CASE( M_TM_PLUS_1 ){
  test_m_tm_plus("/m_tm_plus_input_1");
}

BOOST_AUTO_TEST_CASE( M_TM_PLUS_2 ){
  test_m_tm_plus("/m_tm_plus_input_2");
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GAMMA5 ) 

BOOST_AUTO_TEST_CASE( GAMMA5_1){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GAMMA5_EO) 

BOOST_AUTO_TEST_CASE( GAMMA5_EO_1){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL ) 

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL1_){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_INVERSE_SITEDIAGONAL ) 

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_1){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_MINUS ) 

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_1){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_INVERSE_SITEDIAGONAL_MINUS) 

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_1 ){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_EO ) 

BOOST_AUTO_TEST_CASE( DSLASH_EO_1){
  test_dslash_eo("/dslash_eo_input_1");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_2){
  test_dslash_eo("/dslash_eo_input_2");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_3){
  test_dslash_eo("/dslash_eo_input_3");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_4){
  test_dslash_eo("/dslash_eo_input_4");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_AND_GAMMA5_EO ) 

BOOST_AUTO_TEST_CASE( DSLASH_AND_GAMMA5_EO_1){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO ) 

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_1){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO ) 

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_1){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO ) 

BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO_1){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO ) 

BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO_1){
  BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

