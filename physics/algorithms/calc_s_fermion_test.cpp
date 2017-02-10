/** @file
 * Tests of the fermion force algorithms
 *
 * Copyright (c) 2016 Christopher Czaban
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

#include "metropolis.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::calc_s_fermion
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(calc_s_fermion_staggered){

	using namespace physics::algorithms;
	using namespace physics::lattices;

	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(4, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};

	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	Rational_Approximation approx(15,1,4,1e-5,1);
	Rooted_Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);

	pseudo_randomize<Staggeredfield_eo, su3vec>(&sf, 123);

	hmc_float s_fermion = calc_s_fermion(gf, sf, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());

	BOOST_CHECK_CLOSE(s_fermion, 260.29470334599131, 1.e-8);
}

BOOST_AUTO_TEST_CASE(calc_s_fermion_rooted_wilson){

	using namespace physics::algorithms;
	using namespace physics::lattices;
	std::cout << "Creating parameter string" << std::endl;
	const char * _params[] = {"foo", "--ntime=4", "--num_dev=1", "--beta=5.69"};
	meta::Inputparameters params(4, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};

	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	//In the following N_f=2 flavours are approximated with a rational approximation
	//Rational_Approximation approx(15,99999999,100000000,1e-5,1,1);
	Rational_Approximation approx("Nf2_approximation_for_fermion_force_shifted_test"); // x = 99999999 and y = 100000000 are chosen such that N_f=2 is approximated
	wilson::Rooted_Spinorfield sf(system, interfacesHandler.getInterface<wilson::Rooted_Spinorfield>(), approx);

	pseudo_randomize<Spinorfield, spinor>(&sf, 123);

	hmc_float s_fermion = calc_s_fermion(gf, sf, system, interfacesHandler, interfacesHandler.getAdditionalParameters<wilson::Rooted_Spinorfield>());

	//TODO: Result still has to be checked by true analytic test
	BOOST_CHECK_CLOSE(s_fermion, 2655.639679253864, 1.e-8);
}

BOOST_AUTO_TEST_CASE(calc_s_fermion_wilson){

	using namespace physics::algorithms;
	using namespace physics::lattices;
	std::cout << "Creating parameter string" << std::endl;
	const char * _params[] = {"foo", "--ntime=4", "--num_dev=1", "--beta=5.69", "--solver=cg"};
	meta::Inputparameters params(5, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};

	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	pseudo_randomize<Spinorfield, spinor>(&sf, 123);

	hmc_float s_fermion = calc_s_fermion(gf, sf, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield>());

	BOOST_CHECK_CLOSE(s_fermion, 2655.6396963076531, 1.e-8);
}

BOOST_AUTO_TEST_CASE(calc_s_fermion_wilson_eo){

	using namespace physics::algorithms;
	using namespace physics::lattices;
	std::cout << "Creating parameter string" << std::endl;
	const char * _params[] = {"foo", "--ntime=4", "--num_dev=1", "--beta=5.69", "--solver=cg"};
	meta::Inputparameters params(5, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};

	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	Spinorfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

	pseudo_randomize<Spinorfield_eo, spinor>(&sf, 123);

	hmc_float s_fermion = calc_s_fermion(gf, sf, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());

	BOOST_CHECK_CLOSE(s_fermion, 1086.431332953646, 1.e-8);
}

BOOST_AUTO_TEST_CASE(calc_s_fermion_rooted_wilson_eo){

	using namespace physics::algorithms;
	using namespace physics::lattices;
	std::cout << "Creating parameter string" << std::endl;
	const char * _params[] = {"foo", "--ntime=4", "--num_dev=1", "--beta=5.69"};
	meta::Inputparameters params(4, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};

	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	//In the following N_f=2 flavours are approximated with a rational approximation
	//Rational_Approximation approx(15,99999999,100000000,1e-5,1,1);
	Rational_Approximation approx("Nf2_approximation_for_fermion_force_shifted_test"); // x = 99999999 and y = 100000000 are chosen such that N_f=2 is approximated
	wilson::Rooted_Spinorfield_eo sf(system, interfacesHandler.getInterface<wilson::Rooted_Spinorfield_eo>(), approx);

	pseudo_randomize<Spinorfield_eo, spinor>(&sf, 123);

	hmc_float s_fermion = calc_s_fermion(gf, sf, system, interfacesHandler, interfacesHandler.getAdditionalParameters<wilson::Rooted_Spinorfield_eo>());

	//TODO: Result still has to be checked by true analytic test
	BOOST_CHECK_CLOSE(s_fermion, 1086.431332953646, 1.e-6);
}
