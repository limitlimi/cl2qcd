/** @file
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 *
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

#include "gaugefield.hpp"
#include "matrix6x6Field.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Gaugefield
#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"
#include "../../meta/type_ops.hpp"
#include <stdexcept>

#include "../observables/gaugeObservables.hpp"
#include "../../interfaceImplementations/latticesParameters.hpp"
#include "../../interfaceImplementations/observablesParameters.hpp"
#include "../../interfaceImplementations/physicsParameters.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

BOOST_AUTO_TEST_CASE(set_field)
{
	using namespace physics::lattices;

	{
		const char * _params[] = {"foo", "--ntime=8", "--fermact=clover", "--csw=0.1"};
		meta::Inputparameters params(4, _params);
		const GaugefieldParametersImplementation parametersTmp{ &params };
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		logger.debug() << "Devices: " << system.get_devices().size();
		physics::PrngParametersImplementation prngParameters(params);
		physics::PRNG prng(system, &prngParameters);
		physics::observables::GaugeObservablesParametersImplementation gaugeobservablesParameters(params);

		Matrix6x6Field mat2(system, &parametersTmp);
		Matrix6x6Field mat3(system, &parametersTmp);
		Matrix6x6Field mat4(system, &parametersTmp);
		Matrix6x6Field mat5(system, &parametersTmp);

		// init gaugefield hot
		Gaugefield gf2(system, &parametersTmp, prng, true);

		mat2.setField(&gf2, true);
		BOOST_CHECK_CLOSE(count_Matrix6x6Field(mat2), 12345., 0.1);

		mat3.setField(&gf2, false);
		BOOST_CHECK_CLOSE(count_Matrix6x6Field(mat3), 12345., 0.1);


		// init gaugefield cold
		Gaugefield gf3(system, &parametersTmp, prng, false);

		mat4.setField(&gf3, true);
		BOOST_CHECK_CLOSE(count_Matrix6x6Field(mat4), 3072., 0.1);

		mat5.setField(&gf3, false);
		BOOST_CHECK_CLOSE(count_Matrix6x6Field(mat5), 3072., 0.1);
	}

	{
		const char * _params[] = {"foo", "--ntime=4", "--fermact=clover", "--csw=0.1"};
		meta::Inputparameters params(4, _params);
		const GaugefieldParametersImplementation parametersTmp{ &params };
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		logger.debug() << "Devices: " << system.get_devices().size();
		physics::PrngParametersImplementation prngParameters(params);
		physics::PRNG prng(system, &prngParameters);
		physics::observables::GaugeObservablesParametersImplementation gaugeobservablesParameters(params);

		Matrix6x6Field mat(system, &parametersTmp);
  	    Matrix6x6Field mat1(system, &parametersTmp);

		// init gaugefield from file
		Gaugefield gf(system, &parametersTmp, prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");

		mat.setField(&gf, true);
		BOOST_CHECK_CLOSE(count_Matrix6x6Field(mat), 12345., 0.1);

		mat1.setField(&gf, false);
		BOOST_CHECK_CLOSE(count_Matrix6x6Field(mat1), 12345., 0.1);
	}
}

BOOST_AUTO_TEST_CASE(S_det_test)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo", "--ntime=4", "--fermact=clover", "--csw=0.1"};
	meta::Inputparameters params(4, _params);
	const GaugefieldParametersImplementation parametersTmp{ &params };
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);
	physics::observables::GaugeObservablesParametersImplementation gaugeobservablesParameters(params);

	// init gaugefield cold
	Gaugefield gf1(system, &parametersTmp, prng, false);

	BOOST_CHECK_EQUAL(S_det(gf1, 0.12345, 0.12345), 0.);

	// init gaugefield hot
	Gaugefield gf2(system, &parametersTmp, prng, true);

	BOOST_CHECK_CLOSE(S_det(gf2, 0.12345, 0.12345), 0.12345, 0.1);

	// init gaugefield from file
	Gaugefield gf3(system, &parametersTmp, prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");

	BOOST_CHECK_CLOSE(S_det(gf3, 0.12345, 0.12345), 0.12345, 0.1);
}
