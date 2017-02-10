/** @file
 * Tests of the fermion force algorithms
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2016 Christopher Czaban <czaban@th.physik.uni-frankfurt.de>
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

#include "fermion_force.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::fermion_force
#include <boost/test/unit_test.hpp>

#include "../lattices/util.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"
#include "find_minmax_eigenvalue.hpp"
#include "../lattices/rooted_spinorfield_eo.hpp"

BOOST_AUTO_TEST_CASE(fermion_force)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 11);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 12);
		gm.zero();

		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, interfacesHandler.getAdditionalParameters<Spinorfield>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 18550.897680606064, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 13);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 14);
		gm.zero();

		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, interfacesHandler.getAdditionalParameters<Spinorfield>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 3561.5546616746424, 0.01);
	}

	{
			using namespace physics::lattices;
			const char * _params[] = {"foo", "--ntime=4"};
			meta::Inputparameters params(2, _params);
	        physics::InterfacesHandlerImplementation interfacesHandler{params};
	        hardware::HardwareParametersImplementation hP(&params);
	        hardware::code::OpenClKernelParametersImplementation kP(params);
	        hardware::System system(hP, kP);
			physics::PrngParametersImplementation prngParameters{params};
			physics::PRNG prng{system, &prngParameters};

			Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
			Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
			Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
			Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

			pseudo_randomize<Spinorfield, spinor>(&sf1, 13);
			gm.zero();

			physics::algorithms::calc_fermion_forces(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<physics::lattices::Spinorfield>());

			BOOST_CHECK_CLOSE(squarenorm(gm), 48291.055206667035, 1e-6);

		}
}

BOOST_AUTO_TEST_CASE(fermion_force_shifted)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;

	const char * _params[] = {"foo", "--ntime=4", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};
	const physics::algorithms::RhmcParametersInterface & parametersInterface = interfacesHandler.getRhmcParametersInterface();

	Rational_Approximation approx(25, 99999999,100000000, 1.e-5,1); // x = 99999999 and y = 100000000 are chosen such that N_f=2 is approximated

	physics::fermionmatrix::QplusQminus Qpm(system, interfacesHandler.getInterface<physics::fermionmatrix::QplusQminus>());
	hmc_float minEigen, maxEigen;

	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

	find_maxmin_eigenvalue(maxEigen, minEigen, Qpm, gf, system, interfacesHandler, parametersInterface.getFindMinMaxPrec(), interfacesHandler.getAdditionalParameters<wilson::Rooted_Spinorfield>());

	wilson::Rooted_Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::wilson::Rooted_Spinorfield>(), approx);

	sf1.Rescale_Coefficients(approx, minEigen, maxEigen);

	pseudo_randomize<Spinorfield, spinor>(&sf1, 13); //it will be A

	gm.zero();
	physics::algorithms::calc_fermion_forces(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<physics::lattices::wilson::Rooted_Spinorfield>());
	//TODO: Result still has to be checked by true analytic test
	BOOST_CHECK_CLOSE(squarenorm(gm), 48291.055042961729, 1.e-6);
}

BOOST_AUTO_TEST_CASE(fermion_force_eo_shifted)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;

	const char * _params[] = {"foo", "--ntime=4", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};
	const physics::algorithms::RhmcParametersInterface & parametersInterface = interfacesHandler.getRhmcParametersInterface();

	Rational_Approximation approx(25, 99999999,100000000, 1.e-5,1); // x = 99999999 and y = 100000000 are chosen such that N_f=2 is approximated

	physics::fermionmatrix::QplusQminus_eo Qpm(system, interfacesHandler.getInterface<physics::fermionmatrix::QplusQminus_eo>());
	hmc_float minEigen, maxEigen;

	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

	find_maxmin_eigenvalue(maxEigen, minEigen, Qpm, gf, system, interfacesHandler, parametersInterface.getFindMinMaxPrec(), interfacesHandler.getAdditionalParameters<wilson::Rooted_Spinorfield_eo>());

	wilson::Rooted_Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::wilson::Rooted_Spinorfield_eo>(), approx);

	sf1.Rescale_Coefficients(approx, minEigen, maxEigen);

	pseudo_randomize<Spinorfield_eo, spinor>(&sf1, 13); //it will be A

	gm.zero();
	physics::algorithms::calc_fermion_forces(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<physics::lattices::wilson::Rooted_Spinorfield_eo>());
	//TODO: Result still has to be checked by true analytic test
	BOOST_CHECK_CLOSE(squarenorm(gm), 48291.055042961729, 1.e-6);
}

BOOST_AUTO_TEST_CASE(fermion_force_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&src, 15);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		BOOST_REQUIRE_SMALL(squarenorm(gm), 0.001);
		physics::algorithms::fermion_force(&gm, sf1, sf2, EVEN, gf, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 6180.0548464319563, 0.01);
		physics::algorithms::fermion_force(&gm, sf1, sf2, ODD, gf, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 18590.240819832132, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&src, 16);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::fermion_force(&gm, sf1, sf2, EVEN, gf, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 1294.880037707632, 0.01);
		physics::algorithms::fermion_force(&gm, sf1, sf2, ODD, gf, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 3659.59932413153, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_fermion_forces)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 21);
		gm.zero();

		physics::algorithms::calc_fermion_forces(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 42199.514415107173, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_fermion_forces_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&src, 22);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::calc_fermion_forces(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 3441.344988280136, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_detratio_forces)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--use_mp=true", "--kappa_mp=.25"};
		meta::Inputparameters params(4, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 21);
		gm.zero();

		physics::algorithms::calc_detratio_forces(&gm, gf, sf1, system, interfacesHandler);
		BOOST_CHECK_CLOSE(squarenorm(gm), 12480.139807156647, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_detratio_forces_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--use_mp=true", "--kappa_mp=.25"};
		meta::Inputparameters params(4, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&src, 22);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::calc_detratio_forces(&gm, gf, sf1, system, interfacesHandler);
		BOOST_CHECK_CLOSE(squarenorm(gm), 33313.511647643441, 0.01);
	}
}











