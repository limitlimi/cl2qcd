/** @file
 * Unit test for the physics::lattices::Spinorfield class
 *
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

#include "spinorfield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Spinorfield
#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"
#include "../../meta/type_ops.hpp"
#include <cmath>
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
}

BOOST_AUTO_TEST_CASE(fields)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);

	auto fields = create_spinorfields(system, 2, interfacesHandler);

	BOOST_CHECK_EQUAL(fields.size(), 2u);

	release_spinorfields(fields);
}

BOOST_AUTO_TEST_CASE(gamma5)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.setZero();
	sf.gamma5();
	BOOST_CHECK_CLOSE(squarenorm(sf), 0., .1);
	sf.cold();
	sf.gamma5();
	BOOST_CHECK_CLOSE(squarenorm(sf), 1., .1);
	physics::lattices::sax(&sf, { -.5, .3}, sf);
	sf.gamma5();
	BOOST_CHECK_CLOSE(squarenorm(sf), .34, .1);
}

BOOST_AUTO_TEST_CASE(zero)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.setGaussian(prng);
	sf.setZero();
	BOOST_CHECK_CLOSE(squarenorm(sf), 0., .1);
}

BOOST_AUTO_TEST_CASE(cold)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.setGaussian(prng);
	sf.cold();
	BOOST_CHECK_CLOSE(squarenorm(sf), 1., .1);
}


BOOST_AUTO_TEST_CASE(gaussian)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.setZero();
	sf.gamma5();
	hmc_float const gamma5 = squarenorm(sf);
	sf.setGaussian(prng);
	BOOST_CHECK_NE(squarenorm(sf), gamma5);
}

BOOST_AUTO_TEST_CASE(squarenorm)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.setZero();
	sf.gamma5();
	hmc_float const gamma5 = physics::lattices::squarenorm(sf);
	BOOST_REQUIRE_CLOSE(gamma5, 0., .1);
	sf.setGaussian(prng);
	BOOST_CHECK_NE(physics::lattices::squarenorm(sf), gamma5);
	sf.setZero();
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0., .1);
}

BOOST_AUTO_TEST_CASE(scalar_product)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gaussian.setGaussian(prng);

	Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	zero.setZero();

	Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	cold.cold();

	Spinorfield gamma(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gamma.setZero();
	gamma.gamma5();
	gamma.gamma5();

	const hmc_complex gaussian_scalar_prod = physics::lattices::scalar_product(gaussian, gaussian);
	const hmc_float gaussian_squarenorm = physics::lattices::squarenorm(gaussian);
	BOOST_CHECK_CLOSE(gaussian_scalar_prod.re, gaussian_squarenorm, .1);
	BOOST_CHECK_CLOSE(gaussian_scalar_prod.im, 0., .1);
	const hmc_complex gaussian_scalar_cold = physics::lattices::scalar_product(gaussian, cold);
	const hmc_complex cold_scalar_gaussian = physics::lattices::scalar_product(cold, gaussian);
	BOOST_CHECK_CLOSE(std::abs(gaussian_scalar_cold.re), std::abs(cold_scalar_gaussian.re), .1);
	BOOST_CHECK_CLOSE(std::abs(gaussian_scalar_cold.im), std::abs(cold_scalar_gaussian.im), .1);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(gamma, zero), hmc_complex_zero);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(zero, gamma), hmc_complex_zero);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(gamma, cold), hmc_complex_zero);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(cold, gamma), hmc_complex_zero);
}

BOOST_AUTO_TEST_CASE(scalar_product_real_part)
{
	using physics::lattices::Spinorfield;

		const char * _params[] = {"foo"};
		meta::Inputparameters params(1, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters(params);
		physics::PRNG prng(system, &prngParameters);

		Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		gaussian.setGaussian(prng);

		Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		zero.setZero();

		Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		cold.cold();

		Spinorfield gamma(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		gamma.setZero();
		gamma.gamma5();
		gamma.gamma5();

		const hmc_float gaussian_scalar_prod = physics::lattices::scalar_product_real_part(gaussian, gaussian);
		const hmc_float gaussian_squarenorm = physics::lattices::squarenorm(gaussian);
		BOOST_CHECK_CLOSE(gaussian_scalar_prod, gaussian_squarenorm, .1);
		const hmc_float gaussian_scalar_cold = physics::lattices::scalar_product_real_part(gaussian, cold);
		const hmc_float cold_scalar_gaussian = physics::lattices::scalar_product_real_part(cold, gaussian);
		BOOST_CHECK_CLOSE(std::abs(gaussian_scalar_cold), std::abs(cold_scalar_gaussian), .1);
		BOOST_CHECK_EQUAL(physics::lattices::scalar_product_real_part(gamma, zero), 0.);
		BOOST_CHECK_EQUAL(physics::lattices::scalar_product_real_part(zero, gamma), 0.);
		BOOST_CHECK_EQUAL(physics::lattices::scalar_product_real_part(gamma, cold), 0.);
		BOOST_CHECK_EQUAL(physics::lattices::scalar_product_real_part(cold, gamma), 0.);
}

BOOST_AUTO_TEST_CASE(sax)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);



	Spinorfield orig_sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	orig_sf.setGaussian(prng);
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	physics::lattices::sax(&sf, {.5, 0}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .25 * physics::lattices::squarenorm(orig_sf), .1);
	physics::lattices::sax(&sf, {2., 0.}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 4 * physics::lattices::squarenorm(orig_sf), .1);
	physics::lattices::sax(&sf, {0., 0.}, orig_sf);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0.);

	orig_sf.cold();
	physics::lattices::sax(&sf, { -.8, .7}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 1.1299999999999968, .1);
	physics::lattices::sax(&sf, {.65, .3}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.51250000000000162, .1);

	physics::lattices::Vector<hmc_float> real_vec(5, system);

	real_vec.store(std::vector<hmc_float>(5, 0.31415));
	for(int i=0; i<5; i++){
		physics::lattices::sax(&sf, real_vec, i, orig_sf);
		BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.09869022250, 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(saxpy)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gaussian.setGaussian(prng);
	Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	cold.cold();
	Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	zero.setZero();
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	physics::lattices::saxpy(&sf, {1., 0.}, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1e-8);
	physics::lattices::saxpy(&sf, {0., 0.}, gaussian, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1e-8);
	physics::lattices::saxpy(&sf, {.3, .1}, cold, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .1, 1e-8);
}

BOOST_AUTO_TEST_CASE(saxpy_real)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);


	Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gaussian.setGaussian(prng);
	Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	cold.cold();
	Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	zero.setZero();
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	const physics::lattices::Scalar<hmc_float> scalar(system);

	scalar.store(1.0);
	physics::lattices::saxpy(&sf, scalar, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1e-8);

	scalar.store(0.0);
	physics::lattices::saxpy(&sf, scalar, gaussian, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1e-8);

	scalar.store(.3);
	physics::lattices::saxpy(&sf, scalar, cold, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .09, 1e-8);
}

BOOST_AUTO_TEST_CASE(saxpy_real_arg)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);


	Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gaussian.setGaussian(prng);
	Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	cold.cold();
	Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	zero.setZero();
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	hmc_float real = 1.;
	physics::lattices::saxpy(&sf, real, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1e-8);

	real = 0.;
	physics::lattices::saxpy(&sf, real, gaussian, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1e-8);

	real = .3;
	physics::lattices::saxpy(&sf, real, cold, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .09, 1e-8);
}

BOOST_AUTO_TEST_CASE(saxpy_real_vec)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);



	Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gaussian.setGaussian(prng);
	Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	cold.cold();
	Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	zero.setZero();
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	physics::lattices::Vector<hmc_float> real_vec(5, system);

	real_vec.store(std::vector<hmc_float>(5,1.0));
	for(int i = 0; i < 5; i++)
	{
		physics::lattices::saxpy(&sf, real_vec, i, gaussian, zero);
		BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1e-8);
	}

	real_vec.store(std::vector<hmc_float>(5,0.0));
	for(int i = 0; i < 5; i++)
	{
		physics::lattices::saxpy(&sf, real_vec, i, gaussian, gaussian);
		BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1e-8);
	}

	real_vec.store(std::vector<hmc_float>(5,0.3));
	for(int i = 0; i < 5; i++)
	{
		physics::lattices::saxpy(&sf, real_vec, i, cold, zero);
		BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .09, 1e-8);
	}
}

BOOST_AUTO_TEST_CASE(saxpby)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	cold.cold();
	Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	zero.setZero();
	Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gaussian.setGaussian(prng);
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	physics::lattices::Vector<hmc_float> real_vec1(5, system);
	physics::lattices::Vector<hmc_float> real_vec2(5, system);
	//Complex
	physics::lattices::saxpby(&sf, {1., 0.}, gaussian, {0., 0.}, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpby(&sf, {0., 0.}, cold, {1., 0.}, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpby(&sf, {0., 0.}, gaussian, {0., 0.}, gaussian);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
	//Vector
	real_vec1.store(std::vector<hmc_float>(5, 0.31415));
	real_vec2.store(std::vector<hmc_float>(5, 0.51413));
	for(int i=0; i<5; i++){
		for(int j=0; j<5; j++){
			physics::lattices::saxpby(&sf, real_vec1, i, cold, real_vec2, j, zero);
			BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.09869022250, 1.e-8);
			physics::lattices::saxpby(&sf, real_vec1, i, cold, real_vec2, j, cold);
			BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.68604775840, 1.e-8);
		}
	}
}

BOOST_AUTO_TEST_CASE(saxsbypz)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gaussian.setGaussian(prng);
	Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	cold.cold();
	Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	zero.setZero();
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	physics::lattices::saxsbypz(&sf, {1., 0.}, gaussian, {0., 0.}, cold, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxsbypz(&sf, {0., 0.}, cold, {1., 0.}, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxsbypz(&sf, {0., 0.}, gaussian, {0., 0.}, gaussian, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxsbypz(&sf, {.3, .7}, cold, {1., 0.}, zero, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 2.18, .1);
	physics::lattices::saxsbypz(&sf, {.1, .3}, zero, {.7, .3}, cold, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 2.98, .1);
}

BOOST_AUTO_TEST_CASE(pseudorandomize)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo", "--num_dev=1"};
	meta::Inputparameters params(2, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.setZero();
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
	physics::lattices::pseudo_randomize<Spinorfield, spinor>(&sf, 123);
	logger.info() << "The squarenorm of the pseudorandomized field is " << physics::lattices::squarenorm(sf);
}
