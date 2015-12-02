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

#include "testUtilities.hpp"

#include "../../meta/util.hpp"
#include "../../host_functionality/host_random.h"
#include "../../physics/prng.hpp"
#include "../device.hpp"
#include "correlator.hpp"
#include "spinors.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_CORRELATORS
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"
#include "testUtilities.hpp"
#include "SpinorTester.hpp"

#include <boost/lexical_cast.hpp>

std::string setArgument_spatialExtent( const int value )
{
	return "--ns=" + boost::lexical_cast<std::string>(value);
}

std::string setArgument_temporalExtent( const int value )
{
	return "--nt=" + boost::lexical_cast<std::string>(value);
}

void fill_sf_with_one_eo(spinor * sf_in, int size, bool eo, meta::Inputparameters &params)
{
  int ns = params.get_nspace();
  int nt = params.get_ntime();
  int x,y,z,t;
  for (x = 0; x<ns;x++){
    for (y = 0; y<ns;y++){
      for (z = 0; z<ns;z++){
	for (t = 0; t<nt;t++){
	  int coord[4];
	  coord[0] = t;
	  coord[1] = x;
	  coord[2] = y;
	  coord[3] = z;
	  int nspace =  get_nspace(coord, params);
	  int global_pos = get_global_pos(nspace, t, params);
	  if (global_pos > size)
	    break;
	  hmc_complex content;
	  if ((x+y+z+t) %2 == 0){
	    if (eo)
	      content = hmc_complex_one;
	    else
	      content = hmc_complex_zero;
	  }
	  else{
	    if (eo)
	      content = hmc_complex_zero;
	    else
	      content = hmc_complex_one;
	  }

	  sf_in[global_pos].e0.e0 = content;
	  sf_in[global_pos].e0.e1 = content;
	  sf_in[global_pos].e0.e2 = content;
	  sf_in[global_pos].e1.e0 = content;
	  sf_in[global_pos].e1.e1 = content;
	  sf_in[global_pos].e1.e2 = content;
	  sf_in[global_pos].e2.e0 = content;
	  sf_in[global_pos].e2.e1 = content;
	  sf_in[global_pos].e2.e2 = content;
	  sf_in[global_pos].e3.e0 = content;
	  sf_in[global_pos].e3.e1 = content;
	  sf_in[global_pos].e3.e2 = content;
	}}}}
  return;
}

hmc_float count_sf_eo(spinor * sf_in, int size, bool eo, meta::Inputparameters &params)
{
  int ns = params.get_nspace();
  int nt = params.get_ntime();
  int x,y,z,t;
  hmc_float sum = 0.;
  for (x = 0; x<ns;x++){
    for (y = 0; y<ns;y++){
      for (z = 0; z<ns;z++){
	for (t = 0; t<nt;t++){
	  int coord[4];
	  coord[0] = t;
	  coord[1] = x;
	  coord[2] = y;
	  coord[3] = z;
	  int nspace =  get_nspace(coord, params);
	  int global_pos = get_global_pos(nspace, t, params);
	  if (global_pos > size)
	    break;
	  if (
	      ( eo ==true && (x+y+z+t) %2 == 0) ||
	      ( eo ==false &&  (x+y+z+t) %2 == 1 )
	      )
	    {
	      int i = global_pos;
	      sum +=
		sf_in[i].e0.e0.re+ sf_in[i].e0.e0.im
		+sf_in[i].e0.e1.re+ sf_in[i].e0.e1.im
		+sf_in[i].e0.e2.re+ sf_in[i].e0.e2.im
		+sf_in[i].e1.e0.re+ sf_in[i].e0.e0.im
		+sf_in[i].e1.e1.re+ sf_in[i].e0.e1.im
		+sf_in[i].e1.e2.re+ sf_in[i].e0.e2.im
		+sf_in[i].e2.e0.re+ sf_in[i].e0.e0.im
		+sf_in[i].e2.e1.re+ sf_in[i].e0.e1.im
		+sf_in[i].e2.e2.re+ sf_in[i].e0.e1.im
		+sf_in[i].e3.e0.re+ sf_in[i].e0.e0.im
		+sf_in[i].e3.e1.re+ sf_in[i].e0.e1.im
		+sf_in[i].e3.e2.re+ sf_in[i].e0.e1.im;
	    }
	  else{
	    continue;
	  }
	}}}}
  return sum;
}

hmc_float count_sf(spinor * sf_in, int size)
{
  hmc_float sum = 0.;
  for (int i = 0; i<size;i++){
    sum +=
       sf_in[i].e0.e0.re+ sf_in[i].e0.e0.im
      +sf_in[i].e0.e1.re+ sf_in[i].e0.e1.im
      +sf_in[i].e0.e2.re+ sf_in[i].e0.e2.im
      +sf_in[i].e1.e0.re+ sf_in[i].e1.e0.im
      +sf_in[i].e1.e1.re+ sf_in[i].e1.e1.im
      +sf_in[i].e1.e2.re+ sf_in[i].e1.e2.im
      +sf_in[i].e2.e0.re+ sf_in[i].e2.e0.im
      +sf_in[i].e2.e1.re+ sf_in[i].e2.e1.im
      +sf_in[i].e2.e2.re+ sf_in[i].e2.e2.im
      +sf_in[i].e3.e0.re+ sf_in[i].e3.e0.im
      +sf_in[i].e3.e1.re+ sf_in[i].e3.e1.im
      +sf_in[i].e3.e2.re+ sf_in[i].e3.e2.im;
  }
  return sum;
}

hmc_float calc_var(hmc_float in, hmc_float mean){ //@todo: this is already defined in spinors test!
  return (in - mean) * (in - mean);
}

hmc_float calc_var_sf(spinor * sf_in, int size, hmc_float sum){
  hmc_float var = 0.;
  for(int k = 0; k<size; k++){
    var +=
      calc_var( sf_in[k].e0.e0.re , sum)
      + calc_var( sf_in[k].e0.e0.im , sum)
      + calc_var( sf_in[k].e0.e1.re , sum)
      + calc_var( sf_in[k].e0.e1.im , sum)
      + calc_var( sf_in[k].e0.e2.re , sum)
      + calc_var( sf_in[k].e0.e2.im , sum)
      + calc_var( sf_in[k].e1.e0.re , sum)
      + calc_var( sf_in[k].e1.e0.im , sum)
      + calc_var( sf_in[k].e1.e1.re , sum)
      + calc_var( sf_in[k].e1.e1.im , sum)
      + calc_var( sf_in[k].e1.e2.re , sum)
      + calc_var( sf_in[k].e1.e2.im , sum)
      + calc_var( sf_in[k].e2.e0.re , sum)
      + calc_var( sf_in[k].e2.e0.im , sum)
      + calc_var( sf_in[k].e2.e1.re , sum)
      + calc_var( sf_in[k].e2.e1.im , sum)
      + calc_var( sf_in[k].e2.e2.re , sum)
      + calc_var( sf_in[k].e2.e2.im , sum)
      + calc_var( sf_in[k].e3.e0.re , sum)
      + calc_var( sf_in[k].e3.e0.im , sum)
      + calc_var( sf_in[k].e3.e1.re , sum)
      + calc_var( sf_in[k].e3.e1.im , sum)
      + calc_var( sf_in[k].e3.e2.re , sum)
      + calc_var( sf_in[k].e3.e2.im , sum);
  }
  return var;
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

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_correlators";
	logger.info() << "Init device";
	auto params = createParameters("correlator" + inputfile);
	hardware::System system(*params);
	for(auto device: system.get_devices()) {
		device->getCorrelatorCode();
	}
	BOOST_MESSAGE("Test done");
}

void test_src_volume(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "create_volume_source";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	auto params = createParameters("correlator" + inputfile);
	hardware::System system(*params);

	physics::ParametersPrng_fromMetaInputparameters prngParameters{&(*params)};
	physics::PRNG prng{system, &prngParameters};
	cl_int err = CL_SUCCESS;
	auto * device = system.get_devices().at(0)->getCorrelatorCode();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = params->get_nspace() * params->get_nspace() * params->get_nspace() * params->get_ntime(); //todo: make proper
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	//CP: run the kernel a couple of times times
	int iterations = params->get_integrationsteps(0);

	spinor * sf_out;
	sf_out = new spinor[NUM_ELEMENTS_SF * iterations];
	BOOST_REQUIRE(sf_out);

	auto prng_buf = prng.get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  logger.info() << "Run kernel";
	  out.clear();
	  device->create_volume_source_device(&out, prng_buf);
	  out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
	  sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
	}
	logger.info() << "result: mean";
	hmc_float cpu_res = 0.;
	sum = sum/iterations/NUM_ELEMENTS_SF/24;
	cpu_res= sum;
	logger.info() << cpu_res;

	if(params->get_read_multiple_configs()  == false){
	  //CP: calc std derivation
	  hmc_float var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
	  }
	  var=var/iterations/NUM_ELEMENTS_SF/24;

	  cpu_res = sqrt(var);
	  logger.info() << "result: variance";
	  logger.info() << cpu_res;
	}

	if(params->get_sourcecontent() == common::one){
	  testFloatAgainstInputparameters(cpu_res, *params);
	} else{
	  testFloatSizeAgainstInputparameters(cpu_res, *params);
	}
	BOOST_MESSAGE("Test done");
}

void test_src_zslice(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "create_zslice_source";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	auto params = createParameters("correlator" + inputfile);
	hardware::System system(*params);

	physics::ParametersPrng_fromMetaInputparameters prngParameters{&(*params)};
	physics::PRNG prng{system, &prngParameters};
	cl_int err = CL_SUCCESS;
	auto device = system.get_devices().at(0)->getCorrelatorCode();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = params->get_nspace() * params->get_nspace() * params->get_nspace() * params->get_ntime(); //todo: make proper
	//CP: this source does have a weight only on one slice
	//todo: must be params->get_ntime() * params->get_nspace() * params->get_nspace();
	size_t NUM_ELEMENTS_SRC = params->get_nspace() * params->get_nspace() * params->get_nspace(); //todo: make proper
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	//CP: run the kernel a couple of times times
	int iterations = params->get_integrationsteps(0);

	spinor * sf_out;
	sf_out = new spinor[NUM_ELEMENTS_SF * iterations];
	BOOST_REQUIRE(sf_out);

	auto prng_buf = prng.get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  logger.info() << "Run kernel";
	  out.clear();
	  device->create_zslice_source_device(&out, prng_buf, params->get_source_z());
	  out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
	  sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
	}
	logger.info() << "result: mean";
	hmc_float cpu_res = 0.;
	sum = sum/iterations/NUM_ELEMENTS_SRC/24;
	cpu_res= sum;
	logger.info() << cpu_res;

	if(params->get_read_multiple_configs()  == false){
	  //CP: calc std derivation
	  hmc_float var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
	  }
	  var=var/iterations/NUM_ELEMENTS_SRC/24;

	  cpu_res = sqrt(var);
	  logger.info() << "result: variance";
	  logger.info() << cpu_res;
	}

	if(params->get_sourcecontent() == common::one){
	  testFloatAgainstInputparameters(cpu_res, *params);
	} else{
	  testFloatSizeAgainstInputparameters(cpu_res, *params);
	}
	BOOST_MESSAGE("Test done");
}


BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
  test_build("/correlator_build_input_1");
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("/correlator_build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SRC_VOLUME)

BOOST_AUTO_TEST_CASE( SRC_VOLUME_1 )
{
  test_src_volume("/src_volume_input_1");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_2 )
{
  test_src_volume("/src_volume_input_2");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_3 )
{
  test_src_volume("/src_volume_input_3");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_4 )
{
  test_src_volume("/src_volume_input_4");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_5 )
{
  test_src_volume("/src_volume_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SRC_ZSLICE)

BOOST_AUTO_TEST_CASE( SRC_ZSLICE_1 )
{
  test_src_zslice("/src_zslice_input_1");
}

BOOST_AUTO_TEST_CASE( SRC_ZSLICE_2 )
{
  test_src_zslice("/src_zslice_input_2");
}

BOOST_AUTO_TEST_CASE( SRC_ZSLICE_3 )
{
  test_src_zslice("/src_zslice_input_3");
}

BOOST_AUTO_TEST_CASE( SRC_ZSLICE_4 )
{
  test_src_zslice("/src_zslice_input_4");
}

BOOST_AUTO_TEST_CASE( SRC_ZSLICE_5 )
{
  test_src_zslice("/src_zslice_input_5");
}

BOOST_AUTO_TEST_SUITE_END()


bool checkVariance( const meta::Inputparameters * parameters)
{
	return parameters->get_read_multiple_configs();
}

int getIterationNumber( const meta::Inputparameters * parameters)
{
	return parameters->get_integrationsteps(0);
}

std::string setCoverArgument_checkVariance( const std::string value )
{
	return "--read_multiple_configs=" + value;
}

std::string setCoverArgument_acceptancePrecision( const std::string value )
{
	return "--solver_prec=" + value;
}

std::string setCoverArgument_iterationSteps( const std::string value )
{
	return "--integrationsteps0=" + value;
}

std::string setCoverArgument_useRandomNumbersAsInput( const bool value = false )
{
	if (value)
	{
		return "--solver=bicgstab";
	}
	else
	{
		return "--solver=cg";
	}

}

BOOST_AUTO_TEST_SUITE(SRC_TSLICE)

	void test_src_tslice( const std::vector<std::string> parameterStrings )
	{
		using namespace hardware::buffers;

		std::string kernelName;
		kernelName = "create_timeslice_source";
		printKernelInfo(kernelName);
		logger.info() << "Init device";
		auto params = createParameters(parameterStrings).release();
		hardware::System system(*params);

		physics::ParametersPrng_fromMetaInputparameters prngParameters{&(*params)};
		physics::PRNG prng{system, &prngParameters};
		cl_int err = CL_SUCCESS;
		auto device = system.get_devices().at(0)->getCorrelatorCode();

		logger.info() << "Fill buffers...";
		size_t NUM_ELEMENTS_SF = params->get_nspace() * params->get_nspace() * params->get_nspace() * params->get_ntime(); //todo: make proper
		//CP: this source does have a weight only on one slice
		size_t NUM_ELEMENTS_SRC = params->get_nspace() * params->get_nspace() * params->get_nspace(); //todo: make proper
		const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
		hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
		BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

		int iterations = getIterationNumber( params );

		spinor * sf_out;
		sf_out = new spinor[NUM_ELEMENTS_SF * iterations];
		BOOST_REQUIRE(sf_out);

		auto prng_buf = prng.get_buffers().at(0);

		hmc_float sum = 0;
		for (int i = 0; i< iterations; i++){
		  logger.info() << "Run kernel";
		  out.clear();
		  device->create_timeslice_source_device(&out, prng_buf, params->get_source_t());
		  out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
		  sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
		}
		logger.info() << "result: mean";
		hmc_float cpu_res = 0.;
		sum = sum/iterations/NUM_ELEMENTS_SRC/24;
		cpu_res= sum;
		logger.info() << cpu_res;

		if( checkVariance( params ) ){
		  //CP: calc std derivation
		  hmc_float var=0.;
		  for (int i=0; i<iterations; i++){
			var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
		  }
		  var=var/iterations/NUM_ELEMENTS_SRC/24;

		  cpu_res = sqrt(var);
		  logger.info() << "result: variance";
		  logger.info() << cpu_res;
		}

		if(params->get_sourcecontent() == common::one){
		  testFloatAgainstInputparameters(cpu_res, *params);
		} else{
		  testFloatSizeAgainstInputparameters(cpu_res, *params);
		}
		BOOST_MESSAGE("Test done");
	}

	BOOST_AUTO_TEST_CASE( SRC_TSLICE_1 )
	{
		std::vector<std::string> parameterStrings {"--nt=6", "--ns=4", setCoverArgument_iterationSteps("1"), "--sourcecontent=one", "--measure_pbp=true",
			setCoverArgument_acceptancePrecision("1e-8"), "--test_ref_val=.5", "--use_eo=false", setCoverArgument_checkVariance("false"), "--sourcetype=timeslice", "--measure_correlators=false"};
		test_src_tslice(parameterStrings);
	}

	BOOST_AUTO_TEST_CASE( SRC_TSLICE_2 )
	{
		std::vector<std::string> parameterStrings {"--nt=4", "--ns=4", setCoverArgument_iterationSteps("1000"), "--sourcecontent=z4",
			setCoverArgument_acceptancePrecision("1e-8"), "--test_ref_val=0.001", "--use_eo=false", setCoverArgument_checkVariance("false"), "--sourcetype=timeslice", "--measure_correlators=false"};
		test_src_tslice(parameterStrings);
	}

	BOOST_AUTO_TEST_CASE( SRC_TSLICE_3 )
	{
		std::vector<std::string> parameterStrings {"--nt=4", "--ns=4", setCoverArgument_iterationSteps("100"), "--sourcecontent=one",
			setCoverArgument_acceptancePrecision("1e-8"), "--test_ref_val=.5", "--use_eo=false", setCoverArgument_checkVariance("false"), "--sourcetype=timeslice", "--measure_correlators=false"};
		test_src_tslice(parameterStrings);
	}

	BOOST_AUTO_TEST_CASE( SRC_TSLICE_4 )
	{
		std::vector<std::string> parameterStrings {"--nt=4", "--ns=4", setCoverArgument_iterationSteps("2000"), "--sourcecontent=gaussian",
			setCoverArgument_acceptancePrecision("1e-8"), "--test_ref_val=1.2", "--use_eo=false", setCoverArgument_checkVariance("true"), "--sourcetype=timeslice", "--measure_correlators=false"};
		test_src_tslice(parameterStrings);
	}

	BOOST_AUTO_TEST_CASE( SRC_TSLICE_5 )
	{
		std::vector<std::string> parameterStrings {"--nt=4", "--ns=4", setCoverArgument_iterationSteps("100"), "--sourcecontent=gaussian",
			setCoverArgument_acceptancePrecision("1e-8"), "--test_ref_val=1.2", "--use_eo=false", setCoverArgument_checkVariance("true"), "--sourcetype=timeslice", "--measure_correlators=false"};
		test_src_tslice(parameterStrings);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SRC_POINT)

	void test_src_point(const std::vector<std::string> parameterStrings )
	{
		using namespace hardware::buffers;

		std::string kernelName;
		kernelName = "create_point_source";
		printKernelInfo(kernelName);
		logger.info() << "Init device";
		auto params = createParameters(parameterStrings).release();
		hardware::System system(*params);

		physics::ParametersPrng_fromMetaInputparameters prngParameters{&(*params)};
		physics::PRNG prng{system, &prngParameters};
		cl_int err = CL_SUCCESS;
		auto device = system.get_devices().at(0)->getCorrelatorCode();

		logger.info() << "Fill buffers...";
		size_t NUM_ELEMENTS_SF = params->get_nspace() * params->get_nspace() * params->get_nspace() * params->get_ntime(); //todo: make proper
		//CP: this source does have a weight only on one site
		size_t NUM_ELEMENTS_SRC = 1;
		const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
		hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
		BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

		int iterations = getIterationNumber( params );

		spinor * sf_out;
		sf_out = new spinor[NUM_ELEMENTS_SF * iterations];
		BOOST_REQUIRE(sf_out);

		hmc_float sum = 0;
		for (int i = 0; i< iterations; i++){
			logger.info() << "Run kernel";
			out.clear();
			device->create_point_source_device(&out,i, get_source_pos_spatial(*params),params->get_source_t());
			out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
			sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
		}
		logger.info() << "result: mean";
		hmc_float cpu_res = 0.;

		sum = sum/iterations/NUM_ELEMENTS_SRC;
		cpu_res= sum;
		logger.info() << cpu_res;

		if( checkVariance( params ) )
		{
			//CP: calc std derivation
			hmc_float var=0.;
			for (int i=0; i<iterations; i++){
				var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
			}
			var=var/iterations/NUM_ELEMENTS_SRC;

			cpu_res = sqrt(var);
			logger.info() << "result: variance";
			logger.info() << cpu_res;
		}

		if(params->get_sourcecontent() == common::one){
			testFloatAgainstInputparameters(cpu_res, *params);
		} else{
			testFloatSizeAgainstInputparameters(cpu_res, *params);
		}
		BOOST_MESSAGE("Test done");
	}

	BOOST_AUTO_TEST_CASE( SRC_POINT_1 )
	{
		std::vector<std::string> parameterStrings {"--nt=4", "--ns=4", setCoverArgument_iterationSteps("12"), "--sourcetype=point",
			setCoverArgument_acceptancePrecision("1e-8"), "--test_ref_val=1.", "--use_eo=false", setCoverArgument_checkVariance("false"), "--measure_correlators=false"};
		test_src_point(parameterStrings);
	}

BOOST_AUTO_TEST_SUITE_END()

enum CorrelatorDirection {temporal=0, spatialX=1, spatialY=2, spatialZ=3}; //@todo: should this be name more general, "direction" or so?
typedef std::string KernelIdentifier;

struct CorrelatorTestParameters : public SpinorTestParameters
{
	CorrelatorTestParameters(LatticeExtents lE, CorrelatorDirection directionIn, SpinorFillTypes sF) :
		TestParameters(lE), SpinorTestParameters(lE, sF), direction(directionIn) {};
	const double kappa=1.;
	CorrelatorDirection direction;
};

struct CorrelatorTester : public NonEvenOddSpinorTester
{
	CorrelatorTester(const KernelIdentifier kI, const ParameterCollection pC, const ReferenceValues rV, const CorrelatorTestParameters tP):
		NonEvenOddSpinorTester(kI, pC, tP, rV), correlatorEntries(rV.size())
	{
		code = device->getCorrelatorCode();
		result = new hardware::buffers::Plain<hmc_float>(correlatorEntries, device);
		result->clear();

		for( unsigned int i=0; i<tP.fillTypes.size(); i++ )
		{
			spinorfields.push_back( new hardware::buffers::Plain<spinor> (elements, device) );
			auto spinorfield = SpinorTester::createSpinorfield(tP.fillTypes.at(i));
			spinorfields.at(i)->load(spinorfield);
			delete[] spinorfield;
		}
	};
	~CorrelatorTester()
	{
		result->dump(&kernelResult.at(0));
	}
protected:
	const int correlatorEntries;
	const hardware::buffers::Plain<hmc_float> * result;
	const hardware::code::Correlator * code;
	std::vector< const hardware::buffers::Plain<spinor>* > spinorfields;
};

struct ComponentwiseCorrelatorTester : public CorrelatorTester
{
	ComponentwiseCorrelatorTester(const KernelIdentifier kI, const ParameterCollection pC, const ReferenceValues rV, const CorrelatorTestParameters tP):
		CorrelatorTester(kI, pC, rV, tP)
	{
		code->correlator(code->get_correlator_kernel(kI), result, spinorfields.at(0) );
	}
};

struct ColorwiseCorrelatorTester : public CorrelatorTester
{
	ColorwiseCorrelatorTester(const KernelIdentifier kI, const ParameterCollection pC, const ReferenceValues rV, const CorrelatorTestParameters tP):
		CorrelatorTester(kI, pC, rV, tP)
	{
		code->correlator(code->get_correlator_kernel(kI), result, spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), spinorfields.at(3));
	}
};

template <class TesterClass>
void callTest(const KernelIdentifier kI, const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	CorrelatorTestParameters parametersForThisTest(lE, cD, sF);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, false);
	hardware::code::OpenClKernelParametersMockupForCorrelators kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.kappa, parametersForThisTest.direction);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(kI, parameterCollection, rV, parametersForThisTest );
}

void testPsCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ComponentwiseCorrelatorTester>("ps", lE, cD, sF, rV);
}

void testAvpsCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ComponentwiseCorrelatorTester>("avps", lE, cD, sF, rV);
}

void testScCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("sc", lE, cD, sF, rV);
}

void testVxCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("vx", lE, cD, sF, rV);
}

void testVyCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("vy", lE, cD, sF, rV);
}

void testVzCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("vz", lE, cD, sF, rV);
}

void testAxCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("ax", lE, cD, sF, rV);
}

void testAyCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("ay", lE, cD, sF, rV);
}

void testAzCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("az", lE, cD, sF, rV);
}

BOOST_AUTO_TEST_SUITE(CORRELATOR_PS_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testPsCorrelator(LatticeExtents{ ns8, nt4}, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero}, ReferenceValues(ns8, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testPsCorrelator(LatticeExtents{ ns8, nt4}, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one}, ReferenceValues(ns8, 48.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_PS_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testPsCorrelator(LatticeExtents {ns8, nt4}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero}, ReferenceValues(nt4, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testPsCorrelator(LatticeExtents {ns8, nt4}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one}, ReferenceValues(nt4, 48.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_SC_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testScCorrelator(LatticeExtents{ ns8, nt4}, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns8, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testScCorrelator(LatticeExtents{ ns8, nt4}, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns8, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testScCorrelator(LatticeExtents{ ns8, nt4}, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns8, 1872));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_SC_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testScCorrelator(LatticeExtents{ ns4, nt8}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt8, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testScCorrelator(LatticeExtents{ ns4, nt8}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt8, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testScCorrelator(LatticeExtents{ ns4, nt8}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt8, 1872));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VX_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns8, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns8, 192.));
	}

	BOOST_AUTO_TEST_CASE( nonZero2 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(ns8, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VX_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 192.));
	}

	BOOST_AUTO_TEST_CASE( nonZero2 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(nt12, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VY_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns8, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns8, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(ns8, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VY_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(nt12, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VZ_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns4, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testVzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(ns4, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VZ_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testVzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(nt12, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AX_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAxCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns4, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAxCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAxCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, 144.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AX_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAxCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAxCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAxCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 144.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AY_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAyCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns4, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAyCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAyCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, -144.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AY_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAyCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAyCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAyCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, -144.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AZ_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns4, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, -432.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AZ_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, -432.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AVPS_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAvpsCorrelator(LatticeExtents {ns8, nt8}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero}, ReferenceValues(nt8, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAvpsCorrelator(LatticeExtents {ns8, nt8}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one}, ReferenceValues(nt8, -48.));
	}

BOOST_AUTO_TEST_SUITE_END()

