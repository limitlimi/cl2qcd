/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HARDWARE_CODE_GAUGEMOMENTUM

#include "kernelTester.hpp"
#include "gaugemomentum.hpp"
#include "../../host_functionality/host_random.h"

class GaugemomentumTester : public KernelTester
{
public:
  GaugemomentumTester(std::string kernelName, std::string inputfile, int numberOfValues = 1) :
    KernelTester(kernelName, getSpecificInputfile(inputfile), numberOfValues)
  {
    code = device->get_gaugemomentum_code();
    doubleBuffer = new hardware::buffers::Plain<double> (1, device);

    NUM_ELEMENTS_AE = meta::get_vol4d(*parameters) * NDIM * meta::get_su3algebrasize();
    numberOfGaugemomentumElements = meta::get_vol4d(*parameters) * NDIM;
    useRandom = (parameters->get_solver() == meta::Inputparameters::cg)  ? false : true;
  }

protected:
  std::string getSpecificInputfile(std::string inputfileIn)
  {
    return "gaugemomentum/" + inputfileIn;
  }

  double * createGaugemomentum(int seed = 123456)
  {
    double * gm_in;
    gm_in = new hmc_float[NUM_ELEMENTS_AE];
    useRandom ? fill_with_random(gm_in, seed) : fill_with_one(gm_in);
    BOOST_REQUIRE(gm_in);
    return gm_in;    
  }

   void fill_with_one(hmc_float * sf_in)
   {
     for(int i = 0; i < (int) NUM_ELEMENTS_AE; ++i) {
       sf_in[i] = 1.;
     }
     return;
   }

   void fill_with_random(hmc_float * sf_in, int seed)
   {
     prng_init(seed);
     for(int i = 0; i < (int) NUM_ELEMENTS_AE; ++i) {
       sf_in[i] = prng_double();
     }
     return;
   }

  void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Gaugemomentum * in)
  {
    code->set_float_to_gaugemomentum_squarenorm_device(in, doubleBuffer);
    doubleBuffer->dump(&kernelResult[0]);
  }

  const hardware::code::Gaugemomentum * code;

  hardware::buffers::Plain<double> * doubleBuffer;

  size_t NUM_ELEMENTS_AE;
  size_t numberOfGaugemomentumElements;
  bool useRandom;
};

BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
  GaugemomentumTester tester("build", "opencl_module_gaugemomentum_build_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SQUARENORM )

class SquarenormTester : public GaugemomentumTester
{
public:
  SquarenormTester(std::string inputfile) :
    GaugemomentumTester("gaugemomenta squarenorm", inputfile, 1)
  {
    hardware::buffers::Gaugemomentum in(numberOfGaugemomentumElements, device);
    code->importGaugemomentumBuffer(&in, reinterpret_cast<ae*>( createGaugemomentum() ));
    calcSquarenormAndStoreAsKernelResult(&in);
  }
};

BOOST_AUTO_TEST_CASE(SQUARENORM_1  )
{
	SquarenormTester tester("squarenorm_input_1");
}

BOOST_AUTO_TEST_CASE(SQUARENORM_2  )
{
	SquarenormTester tester("squarenorm_input_2");
}

BOOST_AUTO_TEST_CASE(SQUARENORM_REDUCTION_1  )
{
	SquarenormTester tester("squarenorm_reduction_input_1");
}

BOOST_AUTO_TEST_CASE(SQUARENORM_REDUCTION_2  )
{
	SquarenormTester tester("squarenorm_reduction_input_2");
}

BOOST_AUTO_TEST_CASE(SQUARENORM_REDUCTION_3  )
{
	SquarenormTester tester("squarenorm_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SET_ZERO )

class SetZeroTester : public GaugemomentumTester
{
public:
  SetZeroTester(std::string inputfile) :
    GaugemomentumTester("set zero", inputfile, 1)
  {
    hardware::buffers::Gaugemomentum in(numberOfGaugemomentumElements, device);
    code->importGaugemomentumBuffer(&in, reinterpret_cast<ae*>( createGaugemomentum() ));
    code->set_zero_gaugemomentum(&in);
    calcSquarenormAndStoreAsKernelResult(&in);
  }
};

BOOST_AUTO_TEST_CASE( SET_ZERO_1 )
{
  SetZeroTester tester("set_zero_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SAXPY )

class SaxpyTester : public GaugemomentumTester
{
public:
  SaxpyTester(std::string inputfile) :
    GaugemomentumTester("saxpy", inputfile, 1)
  {
    hardware::buffers::Gaugemomentum in(numberOfGaugemomentumElements, device);
    hardware::buffers::Gaugemomentum out(numberOfGaugemomentumElements, device);
    code->importGaugemomentumBuffer(&in, reinterpret_cast<ae*>( createGaugemomentum(123456) ));
    code->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>( createGaugemomentum(789101) ));
    double alpha = parameters->get_tau();
    doubleBuffer->load(&alpha);

    code->saxpy_device(&in, &out, doubleBuffer, &out);
    calcSquarenormAndStoreAsKernelResult(&out);
  }
};

BOOST_AUTO_TEST_CASE( SAXPY_1 )
{
	SaxpyTester tester("saxpy_input_1");
}

BOOST_AUTO_TEST_CASE( SAXPY_2 )
{
	SaxpyTester tester("saxpy_input_2");
}

BOOST_AUTO_TEST_CASE( SAXPY_3 )
{
	SaxpyTester tester("saxpy_input_3");
}

BOOST_AUTO_TEST_CASE( SAXPY_4 )
{
	SaxpyTester tester("saxpy_input_4");
}

BOOST_AUTO_TEST_CASE( SAXPY_5 )
{
	SaxpyTester tester("saxpy_input_5");
}

BOOST_AUTO_TEST_SUITE_END()


#include "../../meta/util.hpp"
#include "../../physics/prng.hpp"

#include "../tests/test_util.h"

extern std::string const version;
std::string const version = "0.1";

void fill_with_one(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 1.;
	}
	return;
}

void fill_with_random(hmc_float * sf_in, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in[i] = prng_double();
	}
	return;
}

ae make_ae(hmc_float e1, hmc_float e2, hmc_float e3, hmc_float e4,
           hmc_float e5, hmc_float e6, hmc_float e7, hmc_float e8)
{
	ae tmp = {e1, e2, e3, e4, e5, e6, e7, e8};
	return tmp;
}

void fill_with_zero(ae * ae, int size)
{
	for(int i = 0; i < size; ++i) {
		ae[i] = make_ae(0., 0., 0., 0., 0., 0., 0., 0.);
	}
	return;
}

hmc_float count_gm(ae * ae_in, int size)
{
  hmc_float sum = 0.;
  for (int i = 0; i<size;i++){
    sum +=
       ae_in[i].e0
      + ae_in[i].e1
      + ae_in[i].e2
      + ae_in[i].e3
      + ae_in[i].e4
      + ae_in[i].e5
      + ae_in[i].e6
      + ae_in[i].e7;
  }
  return sum;
}

hmc_float calc_var(hmc_float in, hmc_float mean){
  return (in - mean) * (in - mean);
}

hmc_float calc_var_gm(ae * ae_in, int size, hmc_float sum){
  hmc_float var = 0.;
  for(int k = 0; k<size; k++){
    var +=
      calc_var(   ae_in[k].e0 , sum) 
      + calc_var( ae_in[k].e1 , sum) 
      + calc_var( ae_in[k].e2 , sum)
      + calc_var( ae_in[k].e3 , sum) 
      + calc_var( ae_in[k].e4 , sum) 
      + calc_var( ae_in[k].e5 , sum) 
      + calc_var( ae_in[k].e6 , sum) 
      + calc_var( ae_in[k].e7 , sum) ;
  }
  return var;
}

void test_generate_gaussian_gaugemomenta(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "generate_gaussian_gaugemomentum";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);

	physics::PRNG prng(system);
	cl_int err = CL_SUCCESS;
	auto * device = system.get_devices().at(0)->get_gaugemomentum_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();
	size_t NUM_ELEMENTS_GM = meta::get_vol4d(params) * NDIM;
	hardware::buffers::Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	//CP: run the kernel a couple of times times
	int iterations = params.get_integrationsteps(0);

	ae * gm_out;
	gm_out = new ae[NUM_ELEMENTS_GM * iterations];
	BOOST_REQUIRE(gm_out);

	auto gm_code = device->get_device()->get_gaugemomentum_code();
	auto prng_buf = prng.get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  logger.info() << "Run kernel";
	  device->generate_gaussian_gaugemomenta_device(&out, prng_buf);
	  out.dump(&gm_out[i*NUM_ELEMENTS_GM]);
	  sum += count_gm(&gm_out[i*NUM_ELEMENTS_GM], NUM_ELEMENTS_GM);
	}
	logger.info() << "result: mean";
	hmc_float cpu_res = 0.;
	sum = sum/iterations/NUM_ELEMENTS_GM/8;	
	cpu_res= sum;
	logger.info() << cpu_res;

	if(params.get_read_multiple_configs()  == false){
	  //CP: calc std derivation
	  hmc_float var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_gm(&gm_out[i*NUM_ELEMENTS_GM], NUM_ELEMENTS_GM, sum);
	  }
	  var=var/iterations/NUM_ELEMENTS_GM/8;
	  
	  cpu_res = sqrt(var);
	  logger.info() << "result: variance";
	  logger.info() << cpu_res;
	}

	testFloatSizeAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_gm_convert_to_soa(std::string inputfile)
{

}

void test_gm_convert_from_soa(std::string inputfile)
{

}

BOOST_AUTO_TEST_SUITE(GENERATE_GAUSSIAN_GAUGEMOMENTA  )

BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_1 )
{
	test_generate_gaussian_gaugemomenta("/gm_gaussian_input_1");
}

BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_2 )
{
	test_generate_gaussian_gaugemomenta("/gm_gaussian_input_2");
}

BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_3 )
{
	test_generate_gaussian_gaugemomenta("/gm_gaussian_input_3");
}

BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_4 )
{
	test_generate_gaussian_gaugemomenta("/gm_gaussian_input_4");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GM_CONVERT_TO_SOA )

BOOST_AUTO_TEST_CASE( GM_CONVERT_TO_SOA_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_gm_convert_to_soa("/gm_convert_to_soa_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GM_CONVERT_FROM_SOA )

BOOST_AUTO_TEST_CASE( GM_CONVERT_FROM_SOA_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_gm_convert_from_soa("/gm_convert_from_soa_input_1");
}

BOOST_AUTO_TEST_SUITE_END()
