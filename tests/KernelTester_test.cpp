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

#include <iostream>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE KernelTester_test
#include <boost/test/unit_test.hpp>

#include "kernelTester.hpp"

BOOST_AUTO_TEST_SUITE ( BUILD )

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
	  std::string nameOfKernel = "test";
	  std::string nameOfInputfileThatExists = "plaquette_input_1";
	  BOOST_CHECK_NO_THROW(KernelTester kernelTester(nameOfKernel, nameOfInputfileThatExists) );
	}

	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
	  std::string nameOfKernel = "test";
	  std::string nameOfInputfileThatDoesNotExist = "filethatdoesnotexist";
	  BOOST_REQUIRE_THROW(KernelTester kernelTester(nameOfKernel, nameOfInputfileThatDoesNotExist), meta::Inputparameters::parse_aborted  );
	}
	
BOOST_AUTO_TEST_SUITE_END()
