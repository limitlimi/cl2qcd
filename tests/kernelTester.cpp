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

#include "kernelTester.hpp"
#include <boost/test/unit_test.hpp>

KernelTester::KernelTester(std::string kernelNameIn, std::string inputfileIn)
{
  printKernelInformation(kernelNameIn);
  meta::Inputparameters parameters_tmp = createParameters(inputfileIn);

  parameters = &parameters_tmp;

  testPrecision = parameters->get_solver_prec();

  system = new hardware::System(*parameters);
  device = system->get_devices()[0];

  prng = new physics::PRNG(*system);
  gaugefield = new physics::lattices::Gaugefield(*system, *prng);
}

KernelTesterDouble::KernelTesterDouble(std::string kernelNameIn, std::string inputfileIn):
  KernelTester( kernelNameIn, inputfileIn), kernelResult(0.)
{
  referenceValue = parameters->get_test_ref_value();
}

KernelTesterDouble::~KernelTesterDouble()
{
  BOOST_CHECK_CLOSE(kernelResult, referenceValue, testPrecision);
}

KernelTesterComplex::KernelTesterComplex(std::string kernelNameIn, std::string inputfileIn):
  KernelTester( kernelNameIn, inputfileIn), kernelResult({ 0., 0.})
{
  referenceValue = {parameters->get_test_ref_value(), parameters->get_test_ref_value2() };
}

KernelTesterComplex::~KernelTesterComplex()
{
  BOOST_CHECK_CLOSE(kernelResult.re, referenceValue.re, testPrecision);
  BOOST_CHECK_CLOSE(kernelResult.im, referenceValue.im, testPrecision);
}
