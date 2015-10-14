/*
 * Copyright 2012, 2013, 2014 Christopher Pinke, Matthias Bach
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
#define BOOST_TEST_MODULE hardware::code::Gaugefield
#include <boost/test/unit_test.hpp>

#include "testUtilities.hpp"
#include "kernelTester.hpp"
#include "gaugefield.hpp"
#include "../../host_functionality/host_operations_gaugefield.h"
#include "mockups.hpp"

enum FillType {cold = 1, nonTrivial};

//todo: check if this can be moved elsewhere...
void set_cold(Matrixsu3 * field, size_t elems)
{
	for(size_t i = 0; i < elems; ++i) {
		field[i] = unit_matrixsu3();
	}
}

void setGaugefield(Matrixsu3 * field, size_t elems, const FillType fillTypeIn)
{
	for(size_t i = 0; i < elems; ++i)
	{
		switch (fillTypeIn)
		{
		case cold:
			field[i] = unit_matrixsu3();
			break;
		case nonTrivial:
			field[i] = nonTrivialSu3Matrix();
			break;
		default:
			BOOST_ERROR("No valid FillType specified");
		}
	}
}

const Matrixsu3* createGaugefield(const int numberOfElements, const FillType fillTypeIn)
{
	Matrixsu3 * tmp = new Matrixsu3[numberOfElements];
	setGaugefield( tmp, numberOfElements, fillTypeIn);
	return tmp;
}

class GaugefieldTester : public KernelTester {
public:
		GaugefieldTester(std::string kernelName, std::vector<std::string> parameterStrings, const hardware::HardwareParametersMockup & hardwareParameters, const FillType fillType, int numberOfValues = 1):
			KernelTester(kernelName, parameterStrings, numberOfValues), numberOfElements(hardwareParameters.getLatticeVolume() * NDIM) {
		gaugefieldBuffer = new hardware::buffers::SU3( numberOfElements, device);
		const Matrixsu3 * gf_host = createGaugefield(numberOfElements, fillType);
		device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;

		code = device->getGaugefieldCode();
	}

	~GaugefieldTester()
	{
		delete gaugefieldBuffer;
	}

protected:
	const int numberOfElements;
	const hardware::code::Gaugefield * code;
	const hardware::buffers::SU3 * gaugefieldBuffer;

};

BOOST_AUTO_TEST_SUITE ( PLAQUETTE )

	class PlaquetteTester : public GaugefieldTester {
	public:
			PlaquetteTester(std::vector<std::string> parameterStrings, const hardware::HardwareParametersMockup & hardwareParameters, const FillType fillType, double referenceValuePerSite, int typeOfPlaquette = 1):
				GaugefieldTester("plaquette", parameterStrings, hardwareParameters, fillType, 1), typeOfPlaquette(typeOfPlaquette) {
			const hardware::buffers::Plain<hmc_float> plaq(1, device );
			const hardware::buffers::Plain<hmc_float> splaq(1, device);
			const hardware::buffers::Plain<hmc_float> tplaq(1, device);

			code->plaquette_device(gaugefieldBuffer, &plaq, &tplaq, &splaq);

			switch( typeOfPlaquette ) {
				case 1:
					plaq.dump(&kernelResult[0]);
					break;
				case 2:
					tplaq.dump(&kernelResult[0]);
					break;
				case 3:
					splaq.dump(&kernelResult[0]);
					break;
				default:
					throw std::invalid_argument(  "Do not recognize type of plaquette. Should be 1,2 or 3 (normal plaquette, temporal plaquette, spatial plaquette)" );
					break;
			}
		}
	private:
		int typeOfPlaquette;
	};

	BOOST_AUTO_TEST_CASE( PLAQUETTE_1 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		double referenceValuePerSite = 7887186;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS" , "--gaugeact=wilson", "--beta=5.69", "--rho=0", "--solver_prec=1e-8", "--test_ref_val=1536"};
		PlaquetteTester plaquetteTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite, 1);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_2 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::nonTrivial;
		double referenceValuePerSite = 7887186;
		std::vector<std::string> parameterStrings {"--startcondition=hot", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8", "--test_ref_val=1536.002605"};
		PlaquetteTester plaquetteTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite,  1);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_1 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		double referenceValuePerSite = 7887186;
		// todo: add KernelBuilder object
		// add struct which contains: expectedresults (vector) , precision to test, type of comparison, numberOfValues
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8", "--test_ref_val=768."};
		PlaquetteTester plaquetteTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite,  2);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_2 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::nonTrivial;
		double referenceValuePerSite = 7887186;
		std::vector<std::string> parameterStrings {"--startcondition=hot", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8", "--test_ref_val=768.00130250240136"};
		PlaquetteTester plaquetteTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite,  2);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_1 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		double referenceValuePerSite = 7887186;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8", "--test_ref_val=768"};
		PlaquetteTester plaquetteTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite,  3);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_2 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::nonTrivial;
		double referenceValuePerSite = 7887186;
		std::vector<std::string> parameterStrings {"--startcondition=hot", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8", "--test_ref_val=768.00130250240136"};
		PlaquetteTester plaquetteTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite,  3);
	}
	
	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_1 )
	{
		hardware::HardwareParametersMockup hardwareParameters(8,8);
		FillType fillType = FillType::cold;
		double referenceValuePerSite = 7887186;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=8", "--ntime=8", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8", "--test_ref_val=24576"};
		PlaquetteTester plaquetteTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite,  1);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_2 )
	{
		hardware::HardwareParametersMockup hardwareParameters(12,12);
		FillType fillType = FillType::cold;
		double referenceValuePerSite = 7887186;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=12", "--ntime=12", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8", "--test_ref_val=124416"};
		PlaquetteTester plaquetteTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite,  1);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_3 )
	{
		hardware::HardwareParametersMockup hardwareParameters(16,16);
		FillType fillType = FillType::cold;
		double referenceValuePerSite = 7887186;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=16", "--ntime=16", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8", "--test_ref_val=393216"};
		PlaquetteTester plaquetteTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite,  1);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( POLYAKOV )

	class PolyakovloopTester : public GaugefieldTester {
	public:
			PolyakovloopTester(std::vector<std::string> parameterString, const hardware::HardwareParametersMockup & hardwareParameters, const FillType fillType, double referenceValuePerSite, double referenceValuePerSite2):
				GaugefieldTester("polyakov", parameterString, hardwareParameters, fillType, 2) {
			const hardware::buffers::Plain<hmc_complex> pol(1, device);
			code->polyakov_device(gaugefieldBuffer, &pol);

			hmc_complex kernelResult_tmp;
			pol.dump(&kernelResult_tmp);
			kernelResult[0] = kernelResult_tmp.re;
			kernelResult[1] = kernelResult_tmp.im;
		}
	};

	BOOST_AUTO_TEST_CASE( POLYAKOV_1 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		double referenceValuePerSite = 64.;
		double referenceValuePerSite2 = 0.;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=64.", "--test_ref_val2=0."};
		PolyakovloopTester polyakovloopTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite, referenceValuePerSite2);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_2)
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::nonTrivial;
		double referenceValuePerSite = 0.26737144;
		double referenceValuePerSite2 = 0.48554374;
		std::vector<std::string> parameterStrings {"--startcondition=hot", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=-17.1117721375", "--test_ref_val2=-31.0747993518"};
		PolyakovloopTester polyakovloopTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite, referenceValuePerSite2);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_1 )
	{
		hardware::HardwareParametersMockup hardwareParameters(8,8);
		FillType fillType = FillType::cold;
		double referenceValuePerSite = 64.;
		double referenceValuePerSite2 = 0.;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=8", "--ntime=8", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=512.", "--test_ref_val2=0."};
		PolyakovloopTester polyakovloopTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite, referenceValuePerSite2);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_2 )
	{
		hardware::HardwareParametersMockup hardwareParameters(12,12);
		FillType fillType = FillType::cold;
		double referenceValuePerSite = 64.;
		double referenceValuePerSite2 = 0.;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=12", "--ntime=12", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=1728.", "--test_ref_val2=0."};
		PolyakovloopTester polyakovloopTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite, referenceValuePerSite2);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_3 )
	{
		hardware::HardwareParametersMockup hardwareParameters(16,16);
		FillType fillType = FillType::cold;
		double referenceValuePerSite = 64.;
		double referenceValuePerSite2 = 0.;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=16", "--ntime=16", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=4096.", "--test_ref_val2=0."};
		PolyakovloopTester polyakovloopTester(parameterStrings, hardwareParameters, fillType, referenceValuePerSite, referenceValuePerSite2);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( CONVERT_GAUGEFIELD_TO_SOA )

BOOST_AUTO_TEST_CASE( CONVERT_GAUGEFIELD_TO_SOA_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( CONVERT_GAUGEFIELD_FROM_SOA )

BOOST_AUTO_TEST_CASE( CONVERT_GAUGEFIELD_FROM_SOA_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( STOUT_SMEAR )

	class StoutSmearTester : public GaugefieldTester {
	public:
			StoutSmearTester(std::vector<std::string> parameterString, const hardware::HardwareParametersMockup & hardwareParameters, const FillType fillType):
				GaugefieldTester("stout_smear", parameterString, hardwareParameters, fillType) {

			const hardware::buffers::Plain<hmc_float> plaq(1, device );
			const hardware::buffers::Plain<hmc_float> splaq(1, device);
			const hardware::buffers::Plain<hmc_float> tplaq(1, device);
			const hardware::buffers::SU3 out(gaugefieldBuffer->get_elements(), device);

			code->stout_smear_device( gaugefieldBuffer, &out);

			code->plaquette_device( &out, &plaq, &tplaq, &splaq);
			plaq.dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_1 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		std::vector<std::string> parameterStrings {"--startcondition=continue", "--sourcefile=conf.00200", "--use_smearing=true", "--rho_iter=1", "--rho=0.001", "--prec=64", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=882.11113688812929"};
		StoutSmearTester StoutSmearTester(parameterStrings, hardwareParameters, fillType);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_2 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		std::vector<std::string> parameterStrings {"--startcondition=continue", "--sourcefile=conf.00200", "--use_smearing=true", "--rho_iter=1", "--rho=0.", "--prec=64", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=877.17444356279361"};
		StoutSmearTester StoutSmearTester(parameterStrings, hardwareParameters, fillType);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_3 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		std::vector<std::string> parameterStrings {"--startcondition=continue", "--sourcefile=conf.00200", "--use_smearing=true", "--rho_iter=1", "--rho=0.001538", "--prec=64", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=884.76195718059716"};
		StoutSmearTester StoutSmearTester(parameterStrings, hardwareParameters, fillType);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_4 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--use_smearing=true", "--rho_iter=1", "--rho=0.", "--prec=64", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=1536"};
		StoutSmearTester StoutSmearTester(parameterStrings, hardwareParameters, fillType);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_5 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--use_smearing=true", "--rho_iter=1", "--rho=0.001", "--prec=64", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=1536."};
		StoutSmearTester StoutSmearTester(parameterStrings, hardwareParameters, fillType);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( RECTANGLES )

	class RectanglesTester : public GaugefieldTester {
	public:
			RectanglesTester(std::vector<std::string> parameterString, const hardware::HardwareParametersMockup & hardwareParameters, const FillType fillType):
				GaugefieldTester("rectangles", parameterString, hardwareParameters, fillType) {
			const hardware::buffers::Plain<hmc_float> rect(1, device );
			code->rectangles_device(gaugefieldBuffer, &rect);
			rect.dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE( RECTANGLES_1 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--solver=BICGSTAB", "--gaugeact=tlsym", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=3072."};
		RectanglesTester rectanglesTester(parameterStrings, hardwareParameters, fillType);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_2 )
	{
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		FillType fillType = FillType::cold;
		std::vector<std::string> parameterStrings {"--startcondition=continue", "--sourcefile=conf.00200", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--solver=BICGSTAB", "--gaugeact=tlsym", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=1103.2398401620"};
		RectanglesTester rectanglesTester(parameterStrings, hardwareParameters, fillType);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_1 )
	{
		hardware::HardwareParametersMockup hardwareParameters(8,8);
		FillType fillType = FillType::cold;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=8", "--ntime=8", "--fermact=TWISTEDMASS", "--solver=BICGSTAB", "--gaugeact=tlsym", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=49152."};
		RectanglesTester rectanglesTester(parameterStrings, hardwareParameters, fillType);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_2 )
	{
		hardware::HardwareParametersMockup hardwareParameters(12,12);
		FillType fillType = FillType::cold;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=12", "--ntime=12", "--fermact=TWISTEDMASS", "--solver=BICGSTAB", "--gaugeact=tlsym", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=248832."};
		RectanglesTester rectanglesTester(parameterStrings, hardwareParameters, fillType);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_3 )
	{
		hardware::HardwareParametersMockup hardwareParameters(16,16);
		FillType fillType = FillType::cold;
		std::vector<std::string> parameterStrings {"--startcondition=cold", "--nspace=16", "--ntime=16", "--fermact=TWISTEDMASS", "--solver=BICGSTAB", "--gaugeact=tlsym", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=786432."};
		RectanglesTester rectanglesTester(parameterStrings, hardwareParameters, fillType);
	}

BOOST_AUTO_TEST_SUITE_END()

