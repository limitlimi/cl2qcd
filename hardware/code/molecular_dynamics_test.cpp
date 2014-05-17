/*
 * Copyright 2012, 2013, 2014 Lars Zeidlewicz, Christopher Pinke,
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
#define BOOST_TEST_MODULE HARDWARE_CODE_MOLECULAR_DYNAMICS

#include "GaugemomentumTester.hpp"
#include "SpinorTester.hpp"
#include "molecular_dynamics.hpp"

#include "../physics/lattices/gaugefield.hpp"

class MolecularDynamicsTester : public GaugemomentumTester
{
public:
	MolecularDynamicsTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1) :
		GaugemomentumTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues)
		{
			gaugefieldCode = device->get_gaugefield_code();
			molecularDynamicsCode = device->get_molecular_dynamics_code();
			
			gaugefield = new physics::lattices::Gaugefield(*system, *prng);
		}
		
		~MolecularDynamicsTester()
		{
			molecularDynamicsCode = NULL;
			gaugefieldCode = NULL;
		}
		
protected:
	std::string getSpecificInputfile(std::string inputfileIn)
	{
		//todo: this is ugly, find a better solution.
		// The problem is that the parent class calls a similar fct.
		return "../molecularDynamics/" + inputfileIn;
	}
	
	const hardware::buffers::SU3* getGaugefieldBuffer() {
		return gaugefield->get_buffers()[0];
	}
	
	const hardware::code::Molecular_Dynamics * molecularDynamicsCode;
	const hardware::code::Gaugefield * gaugefieldCode;
	
	physics::PRNG * prng;
	physics::lattices::Gaugefield * gaugefield;
};

BOOST_AUTO_TEST_SUITE(BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
		MolecularDynamicsTester tester("build", "molecular_dynamics_build_input_1");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GF_UPDATE )

	class GaugefieldUpdateTester : public MolecularDynamicsTester
	{
	public:
		GaugefieldUpdateTester(std::string inputfile) : 
			MolecularDynamicsTester("md_update_gaugefield", inputfile)
			{
				code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentum() ));
				double epsilon = parameters->get_tau();
				molecularDynamicsCode->md_update_gaugefield_device(gaugemomentumBuffer, getGaugefieldBuffer(), epsilon);
				
				const hardware::buffers::Plain<hmc_float> plaq(1, device );
				const hardware::buffers::Plain<hmc_float> splaq(1, device);
				const hardware::buffers::Plain<hmc_float> tplaq(1, device);

				gaugefieldCode->plaquette_device(getGaugefieldBuffer(), &plaq, &tplaq, &splaq);
				plaq.dump(&kernelResult[0]);
			}
	};

	BOOST_AUTO_TEST_CASE(GF_UPDATE_1 )
	{
		GaugefieldUpdateTester tester("gf_update_input_1");
	}

	BOOST_AUTO_TEST_CASE( GF_UPDATE_2 )
	{
		GaugefieldUpdateTester tester("gf_update_input_2");
	}

	BOOST_AUTO_TEST_CASE( GF_UPDATE_3 )
	{
		GaugefieldUpdateTester tester("gf_update_input_3");
	}

	BOOST_AUTO_TEST_CASE( GF_UPDATE_4 )
	{
		GaugefieldUpdateTester tester("gf_update_input_4");
	}

	BOOST_AUTO_TEST_CASE(GF_UPDATE_5 )
	{
		BOOST_MESSAGE("THIS TEST HAS TO BE INVESTIGATED!! IT COULD HINT TO AN ERROR IN THE FUNCTION!");
		GaugefieldUpdateTester tester("gf_update_input_5");
	}

	BOOST_AUTO_TEST_CASE(GF_UPDATE_6 )
	{
		GaugefieldUpdateTester tester("gf_update_input_6");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE )

	class FGaugeTester : public MolecularDynamicsTester
	{
	public:
		FGaugeTester(std::string inputfile) : 
			MolecularDynamicsTester("f_gauge", inputfile)
			{
				code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));
				molecularDynamicsCode->gauge_force_device( getGaugefieldBuffer(), gaugemomentumBuffer);
				calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
			}
	};

	BOOST_AUTO_TEST_CASE( F_GAUGE_1 )
	{
		FGaugeTester tester("/f_gauge_input_1");
	}

	BOOST_AUTO_TEST_CASE( F_GAUGE_2 )
	{
		FGaugeTester tester("f_gauge_input_2");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE_TLSYM )

	class FGaugeTlsymTester : public MolecularDynamicsTester
	{
	public:
		FGaugeTlsymTester(std::string inputfile) : 
			MolecularDynamicsTester("f_gauge_tlsym", inputfile)
			{
				code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));
				molecularDynamicsCode->gauge_force_tlsym_device( getGaugefieldBuffer(), gaugemomentumBuffer);
				calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
			}
	};

	BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_1 )
	{
		FGaugeTlsymTester tester("f_gauge_tlsym_input_1");
	}

	BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_2 )
	{
		FGaugeTlsymTester tester("f_gauge_tlsym_input_2");
	}

	BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_3 )
	{
		FGaugeTlsymTester tester("f_gauge_tlsym_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( STOUT_SMEAR_FERMION_FORCE )

	BOOST_AUTO_TEST_CASE(STOUT_SMEAR_FERMION_FORCE_1 )
	{
		BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION )

	class FFermionTester : public MolecularDynamicsTester, public SpinorTester
	{
	public:
		FFermionTester(std::string inputfile) : 
			MolecularDynamicsTester("f_fermion", inputfile), SpinorTester("f_fermion", MolecularDynamicsTester::getSpecificInputfile(inputfile))
			{
				MolecularDynamicsTester::code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));
				const hardware::buffers::Plain<spinor> in1(SpinorTester::spinorfieldElements, MolecularDynamicsTester::device);
				const hardware::buffers::Plain<spinor> in2(SpinorTester::spinorfieldElements, MolecularDynamicsTester::device);

				in1.load(SpinorTester::createSpinorfield(SpinorTester::spinorfieldElements));
				in2.load(SpinorTester::createSpinorfield(SpinorTester::spinorfieldElements));
				
				molecularDynamicsCode->fermion_force_device( &in1, &in2, getGaugefieldBuffer(), gaugemomentumBuffer, MolecularDynamicsTester::parameters->get_kappa());
				MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
				
				SpinorTester::setReferenceValuesToZero();
			}
	};

	BOOST_AUTO_TEST_CASE( F_FERMION_1 )
	{
		FFermionTester tester("f_fermion_input_1");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_2 )
	{
		FFermionTester tester("f_fermion_input_2");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_3 )
	{
		FFermionTester tester("f_fermion_input_3");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_4 )
	{
		FFermionTester tester("f_fermion_input_4");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_5 )
	{
		FFermionTester tester("f_fermion_input_5");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_6 )
	{
		FFermionTester tester("f_fermion_input_6");
	}

BOOST_AUTO_TEST_SUITE_END()

void fill_sf_with_random_eo(spinor * sf_in1, spinor * sf_in2, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in1[i].e0.e0.re = prng_double();
		sf_in1[i].e0.e1.re = prng_double();
		sf_in1[i].e0.e2.re = prng_double();
		sf_in1[i].e1.e0.re = prng_double();
		sf_in1[i].e1.e1.re = prng_double();
		sf_in1[i].e1.e2.re = prng_double();
		sf_in1[i].e2.e0.re = prng_double();
		sf_in1[i].e2.e1.re = prng_double();
		sf_in1[i].e2.e2.re = prng_double();
		sf_in1[i].e3.e0.re = prng_double();
		sf_in1[i].e3.e1.re = prng_double();
		sf_in1[i].e3.e2.re = prng_double();

		sf_in1[i].e0.e0.im = prng_double();
		sf_in1[i].e0.e1.im = prng_double();
		sf_in1[i].e0.e2.im = prng_double();
		sf_in1[i].e1.e0.im = prng_double();
		sf_in1[i].e1.e1.im = prng_double();
		sf_in1[i].e1.e2.im = prng_double();
		sf_in1[i].e2.e0.im = prng_double();
		sf_in1[i].e2.e1.im = prng_double();
		sf_in1[i].e2.e2.im = prng_double();
		sf_in1[i].e3.e0.im = prng_double();
		sf_in1[i].e3.e1.im = prng_double();
		sf_in1[i].e3.e2.im = prng_double();

		sf_in2[i].e0.e0.re = prng_double();
		sf_in2[i].e0.e1.re = prng_double();
		sf_in2[i].e0.e2.re = prng_double();
		sf_in2[i].e1.e0.re = prng_double();
		sf_in2[i].e1.e1.re = prng_double();
		sf_in2[i].e1.e2.re = prng_double();
		sf_in2[i].e2.e0.re = prng_double();
		sf_in2[i].e2.e1.re = prng_double();
		sf_in2[i].e2.e2.re = prng_double();
		sf_in2[i].e3.e0.re = prng_double();
		sf_in2[i].e3.e1.re = prng_double();
		sf_in2[i].e3.e2.re = prng_double();

		sf_in2[i].e0.e0.im = prng_double();
		sf_in2[i].e0.e1.im = prng_double();
		sf_in2[i].e0.e2.im = prng_double();
		sf_in2[i].e1.e0.im = prng_double();
		sf_in2[i].e1.e1.im = prng_double();
		sf_in2[i].e1.e2.im = prng_double();
		sf_in2[i].e2.e0.im = prng_double();
		sf_in2[i].e2.e1.im = prng_double();
		sf_in2[i].e2.e2.im = prng_double();
		sf_in2[i].e3.e0.im = prng_double();
		sf_in2[i].e3.e1.im = prng_double();
		sf_in2[i].e3.e2.im = prng_double();
	}
	return;
}

BOOST_AUTO_TEST_SUITE( F_FERMION_EO )

	class FFermionEvenOddTester : public MolecularDynamicsTester, public SpinorTester
	{
	public:
		FFermionEvenOddTester(std::string inputfile) : 
			MolecularDynamicsTester("f_fermion_eo", inputfile), SpinorTester("f_fermion_eo", MolecularDynamicsTester::getSpecificInputfile(inputfile))
			{
				MolecularDynamicsTester::code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));
				const hardware::buffers::Spinor in1(SpinorTester::spinorfieldEvenOddElements, MolecularDynamicsTester::device);
				const hardware::buffers::Spinor in2(SpinorTester::spinorfieldEvenOddElements, MolecularDynamicsTester::device);

				spinor * sf_in1;
				spinor * sf_in2;
				sf_in1 = new spinor[SpinorTester::spinorfieldEvenOddElements];
				sf_in2 = new spinor[SpinorTester::spinorfieldEvenOddElements];
				//use the variable use_cg to switch between cold and random input sf
				if(MolecularDynamicsTester::parameters->get_solver() == meta::Inputparameters::cg) {
					SpinorTester::fill_with_one(sf_in1, SpinorTester::spinorfieldEvenOddElements);
					SpinorTester::fill_with_one(sf_in2, SpinorTester::spinorfieldEvenOddElements);
				} else {
					fill_sf_with_random_eo(sf_in1, sf_in2, SpinorTester::spinorfieldEvenOddElements, 123456);
				}
				BOOST_REQUIRE(sf_in1);
				BOOST_REQUIRE(sf_in2);
			
				in1.load(sf_in1);
				in2.load(sf_in2);
				delete sf_in1;
				delete sf_in2;
				
				int tmp = ( MolecularDynamicsTester::parameters->get_read_multiple_configs() ) ? EVEN : ODD;
				molecularDynamicsCode->fermion_force_eo_device( &in1, &in2, getGaugefieldBuffer(), gaugemomentumBuffer, tmp,  MolecularDynamicsTester::parameters->get_kappa());
				MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
				
				SpinorTester::setReferenceValuesToZero();
			}
	};
	
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_1 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_2 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_3 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_4 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_4");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_5 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_5");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_6 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_6");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_7 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_7");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_8 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_8");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_9 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_9");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_10 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_10");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_11)
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_11");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_12 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_12");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_13)
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_13");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_14 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_14");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_15 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_15");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_16 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_16");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_17 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_17");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_18 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_18");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_19 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_19");
	}

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_20 )
	{
		FFermionEvenOddTester tester("f_fermion_eo_input_20");
	}

BOOST_AUTO_TEST_SUITE_END()

#include "../../meta/util.hpp"
#include "../../physics/lattices/gaugefield.hpp"
#include "molecular_dynamics.hpp"
#include "spinors_staggered.hpp"

#include "../../tests/test_util.h"
#include "../../tests/test_util_staggered.h"

class TestMolecularDynamics {

public:
	TestMolecularDynamics(const hardware::System * system) : system(system), prng(*system), gf(*system, prng) {
		BOOST_REQUIRE_EQUAL(system->get_devices().size(), 1);
		auto inputfile = system->get_inputparameters();
		meta::print_info_hmc(inputfile);
	};

	const hardware::code::Molecular_Dynamics * get_device();
	const hardware::buffers::SU3 * get_gaugefield();

	void print_gaugeobservables() {
		physics::lattices::print_gaugeobservables(gf, 0);
	}

private:
	const hardware::System * const system;
	physics::PRNG prng;
	const physics::lattices::Gaugefield gf;
};

const hardware::code::Molecular_Dynamics* TestMolecularDynamics::get_device()
{
	return system->get_devices()[0]->get_molecular_dynamics_code();
}

const hardware::buffers::SU3 * TestMolecularDynamics::get_gaugefield()
{
	return gf.get_buffers().at(0);
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

void fill_sf_with_one(su3vec * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0 = hmc_complex_one;
		sf_in[i].e1 = hmc_complex_one;
		sf_in[i].e2 = hmc_complex_one;
	}
	return;
}

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

void fill_sf_with_random(spinor * sf_in, int size, int seed)
{
	//  Random rnd_loc(seed);
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

void fill_sf_with_random(su3vec * sf_in, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.re = prng_double();
		sf_in[i].e1.re = prng_double();
		sf_in[i].e2.re = prng_double();
		
		sf_in[i].e0.im = prng_double();
		sf_in[i].e1.im = prng_double();
		sf_in[i].e2.im = prng_double();
	}
	return;
}

void fill_sf_with_random_eo(su3vec * sf_in1, su3vec * sf_in2, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in1[i].e0.re = prng_double();
		sf_in1[i].e1.re = prng_double();
		sf_in1[i].e2.re = prng_double();

		sf_in1[i].e0.im = prng_double();
		sf_in1[i].e1.im = prng_double();
		sf_in1[i].e2.im = prng_double();
		
		sf_in2[i].e0.re = prng_double();
		sf_in2[i].e1.re = prng_double();
		sf_in2[i].e2.re = prng_double();

		sf_in2[i].e0.im = prng_double();
		sf_in2[i].e1.im = prng_double();
		sf_in2[i].e2.im = prng_double();
	}
	return;
}

/////////////////////////////////////////////////////////////////////////
//    TESTS FOR STAGGERED FERMIONS MOLECULAR DYNAMICS RELATED TOOLS    //
/////////////////////////////////////////////////////////////////////////

void test_f_stagg_fermion_eo(std::string inputfile)
{
	using namespace hardware::buffers;
	
	std::string kernelName = "fermion_staggered_partial_force_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestMolecularDynamics cpu(&system);
	auto * device = cpu.get_device();
	
	su3vec * sf_in1;
	su3vec * sf_in2;
	ae * ae_out;
	
	logger.info() << "fill buffers";
	size_t NUM_ELEMENTS_SF =  hardware::code::get_eoprec_spinorfieldsize(params);
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM;

	sf_in1 = new su3vec[NUM_ELEMENTS_SF];
	sf_in2 = new su3vec[NUM_ELEMENTS_SF];
	ae_out = new ae[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in1, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random(sf_in1, NUM_ELEMENTS_SF, 123); //With these seeds the fields are the same
		fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 456); //as the test_sf_saxpy_staggered_eo
	}
	fill_with_zero(ae_out, NUM_ELEMENTS_AE);
	BOOST_REQUIRE(sf_in1);
	BOOST_REQUIRE(sf_in2);
	BOOST_REQUIRE(ae_out);
	
	auto spinor_code = device->get_device()->get_spinor_staggered_code();
	auto gm_code = device->get_device()->get_gaugemomentum_code();

	const SU3vec in1(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec in2(NUM_ELEMENTS_SF, device->get_device());
	Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	in1.load(sf_in1);
	in2.load(sf_in2);
	gm_code->importGaugemomentumBuffer(&out, ae_out);

	//The following seven lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
	/*
        print_staggeredfield_eo_to_textfile("ref_vec_f_stagg1_eo", sf_in1, params); 
        logger.info() << "Produced the ref_vec_f_stagg1_eo text file with the staggered field for the ref. code."; 
	print_staggeredfield_eo_to_textfile("ref_vec_f_stagg2_eo", sf_in2, params); 
        logger.info() << "Produced the ref_vec_f_stagg2_eo text file with the staggered field for the ref. code. Returning...";   
        return;
	// */

	hmc_float cpu_res, cpu_back, cpu_back2;
	logger.info() << "|in_1|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in1, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "|in_2|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in2, &sqnorm);
	sqnorm.dump(&cpu_back2);
	logger.info() << cpu_back2;
	logger.info() << "Run kernel";

	//switch according to "read_multiple_configs"
	if(params.get_read_multiple_configs()) {
		device->fermion_staggered_partial_force_device(cpu.get_gaugefield(), &in1, &in2, &out, EVEN);
	} else {
		device->fermion_staggered_partial_force_device(cpu.get_gaugefield(), &in1, &in2, &out, ODD);
	}
	
	logger.info() << "|force|^2:";
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	logger.info() << "Clear buffers";
	delete[] sf_in1;
	delete[] sf_in2;
	delete[] ae_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

///////////////////////////////////////////////////////////////////////////////
//    TESTS SUITE FOR STAGGERED FERMIONS MOLECULAR DYNAMICS RELATED TOOLS    //
///////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE( F_STAGG_FERMION_EO )

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_1 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_1");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_2 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_2");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_3 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_3");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_4 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_4");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_5 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_5");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_6 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_6");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_7 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_7");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_8 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_8");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_9 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_9");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_10 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_10");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_11 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_11");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_12 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_12");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_13 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_13");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_14 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_14");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_15 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_15");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_16 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_16");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_17 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_17");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_18 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_18");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_19 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_19");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_20 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_20");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_21 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_21");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_22 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_22");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_23 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_23");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_24 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_24");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_25 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_25");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_26 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_26");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_27 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_27");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_28 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_28");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_29 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_29");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_30 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_30");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_31 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_31");
}

BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_32 )
{
	test_f_stagg_fermion_eo("/f_staggered_fermion_partial_eo_input_32");
}


BOOST_AUTO_TEST_SUITE_END()







////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

void test_f_fermion_compare_noneo_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName = "Compare f_fermion in eo and noneo version";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestMolecularDynamics cpu(&system);
	auto * device = cpu.get_device();

	spinor * sf_in1_noneo;
	spinor * sf_in2_noneo;
	ae * sf_out_noneo;
	ae * sf_out_eo;
	spinor * sf_in1_eo;
	spinor * sf_in2_eo;
	spinor * sf_in3_eo;
	spinor * sf_in4_eo;

	size_t NUM_ELEMENTS_SF_EO =  hardware::code::get_eoprec_spinorfieldsize(params);
	size_t NUM_ELEMENTS_SF =  hardware::code::get_spinorfieldsize(params);
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM;

	sf_in1_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in2_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in3_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in4_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_out_eo = new ae[NUM_ELEMENTS_AE];
	sf_in1_noneo = new spinor[NUM_ELEMENTS_SF];
	sf_in2_noneo = new spinor[NUM_ELEMENTS_SF];
	sf_out_noneo = new ae[NUM_ELEMENTS_AE];

	fill_with_zero(sf_out_eo, NUM_ELEMENTS_AE);
	fill_with_zero(sf_out_noneo, NUM_ELEMENTS_AE);
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in1_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in2_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in3_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in4_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in1_noneo, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in2_noneo, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random_eo(sf_in1_eo, sf_in2_eo, NUM_ELEMENTS_SF_EO, 123456);
		fill_sf_with_random_eo(sf_in3_eo, sf_in4_eo, NUM_ELEMENTS_SF_EO, 789101);
		fill_sf_with_random(sf_in1_noneo, NUM_ELEMENTS_SF, 123456);
		fill_sf_with_random(sf_in2_noneo, NUM_ELEMENTS_SF, 789101);
	}

	BOOST_REQUIRE(sf_in1_eo);
	BOOST_REQUIRE(sf_in2_eo);
	BOOST_REQUIRE(sf_in3_eo);
	BOOST_REQUIRE(sf_in4_eo);
	BOOST_REQUIRE(sf_in1_noneo);
	BOOST_REQUIRE(sf_in2_noneo);
	BOOST_REQUIRE(sf_out_noneo);
	BOOST_REQUIRE(sf_out_eo);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gm_code = device->get_device()->get_gaugemomentum_code();

	const Spinor in1_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in2_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in3_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in4_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	const Plain<spinor> in1_noneo(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in2_noneo(NUM_ELEMENTS_SF, device->get_device());
	const Gaugemomentum out_noneo(NUM_ELEMENTS_AE, device->get_device());
	const Gaugemomentum out_eo(NUM_ELEMENTS_AE, device->get_device());
	const Plain<hmc_float> sqnorm(1, device->get_device());

	gm_code->importGaugemomentumBuffer(&out_eo, sf_out_eo);
	gm_code->importGaugemomentumBuffer(&out_noneo, sf_out_noneo);

	//in case of rnd input, it is nontrivial to supply the same rnd vectors as eo and noneo input.
	//therefore, simply convert the eo input back to noneo
	if(params.get_solver() == meta::Inputparameters::cg) {
		spinor_code->copy_to_eoprec_spinorfield_buffer(&in1_eo, sf_in1_eo);
		spinor_code->copy_to_eoprec_spinorfield_buffer(&in2_eo, sf_in2_eo);
		spinor_code->copy_to_eoprec_spinorfield_buffer(&in3_eo, sf_in3_eo);
		spinor_code->copy_to_eoprec_spinorfield_buffer(&in4_eo, sf_in4_eo);
		in1_noneo.load(sf_in1_noneo);
		in2_noneo.load(sf_in2_noneo);
	} else {
		//one can either convert to or from eoprec, use read_multiple_configs for that
		//NOTE: there is machinery to compare vectors in the old executable
		if(params.get_read_multiple_configs()) {
			spinor_code->copy_to_eoprec_spinorfield_buffer(&in1_eo, sf_in1_eo);
			spinor_code->copy_to_eoprec_spinorfield_buffer(&in2_eo, sf_in2_eo);
			spinor_code->convert_from_eoprec_device(&in1_eo, &in2_eo, &in1_noneo);

			spinor_code->copy_to_eoprec_spinorfield_buffer(&in3_eo, sf_in3_eo);
			spinor_code->copy_to_eoprec_spinorfield_buffer(&in4_eo, sf_in4_eo);
			spinor_code->convert_from_eoprec_device(&in3_eo, &in4_eo, &in2_noneo);
		}  else {
			in1_noneo.load(sf_in1_noneo);
			in2_noneo.load(sf_in2_noneo);
			spinor_code->convert_to_eoprec_device(&in1_eo, &in2_eo, &in1_noneo);
			spinor_code->convert_to_eoprec_device(&in3_eo, &in4_eo, &in2_noneo);
		}
	}

	hmc_float cpu_back_eo, cpu_back_eo2, cpu_back_eo3, cpu_back_eo4;
	logger.info() << "eo input:";
	logger.info() << "|phi_even_1|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in1_eo, &sqnorm);
	sqnorm.dump(&cpu_back_eo);
	logger.info() << cpu_back_eo;
	logger.info() << "|phi_even_2|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in2_eo, &sqnorm);
	sqnorm.dump(&cpu_back_eo2);
	logger.info() << cpu_back_eo2;
	logger.info() << "|phi_odd_1|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in3_eo, &sqnorm);
	sqnorm.dump(&cpu_back_eo3);
	logger.info() << cpu_back_eo3;
	logger.info() << "|phi_odd_2|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in4_eo, &sqnorm);
	sqnorm.dump(&cpu_back_eo4);
	logger.info() << cpu_back_eo4;

	logger.info() << "run eo force on EVEN and ODD sites...";
	device->fermion_force_eo_device(&in1_eo, &in4_eo, cpu.get_gaugefield(), &out_eo, ODD, params.get_kappa() );
	device->fermion_force_eo_device(&in2_eo, &in3_eo, cpu.get_gaugefield(), &out_eo, EVEN, params.get_kappa() );

	logger.info() << "|force_eo (even) + force_eo (odd)|^2:";
	hmc_float cpu_res_eo;
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out_eo, &sqnorm);
	sqnorm.dump(&cpu_res_eo);
	logger.info() << cpu_res_eo;

	logger.info() << "non-eo input:";
	hmc_float cpu_back_noneo, cpu_back2_noneo;
	logger.info() << "|phi_1|^2:";
	spinor_code->set_float_to_global_squarenorm_device(&in1_noneo, &sqnorm);
	sqnorm.dump(&cpu_back_noneo);
	logger.info() << cpu_back_noneo;
	logger.info() << "|phi_2|^2:";
	spinor_code->set_float_to_global_squarenorm_device(&in2_noneo, &sqnorm);
	sqnorm.dump(&cpu_back2_noneo);
	logger.info() << cpu_back2_noneo;
	logger.info() << "run noneo force with noneo input...";
	device->fermion_force_device( &in1_noneo, &in2_noneo, cpu.get_gaugefield(), &out_noneo, params.get_kappa());
	logger.info() << "|force_noneo|^2:";
	hmc_float cpu_res_noneo;
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out_noneo, &sqnorm);
	sqnorm.dump(&cpu_res_noneo);
	logger.info() << cpu_res_noneo;

	logger.info() << "Clear buffers";
	delete[] sf_in1_eo;
	delete[] sf_in2_eo;
	delete[] sf_out_eo;
	delete[] sf_in3_eo;
	delete[] sf_in4_eo;
	delete[] sf_out_noneo;
	delete[] sf_in1_noneo;
	delete[] sf_in2_noneo;

	logger.info() << "Compare eo and non-eo results";
	BOOST_REQUIRE_CLOSE(cpu_res_eo, cpu_res_noneo, params.get_solver_prec() );
	testFloatAgainstInputparameters(cpu_res_eo, params);
	testFloatAgainstInputparameters(cpu_res_noneo, params);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_SUITE( F_FERMION_COMPARE_NONEO_EO )

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_1 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_1");
}

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_2 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_2");
}

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_3 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_3");
}

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_4 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_4");
}

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_5 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_5");
}

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_6 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_6");
}

BOOST_AUTO_TEST_SUITE_END()