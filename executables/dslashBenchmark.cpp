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

/*
 * @file
 * Declaration of the hmcExecutable class.
 * This class provides features for the generation of gauge configurations
 * according to the Hybrid Monte Carlo (HMC) algorithm.
 */


#include "dslashBenchmark.h"

dslashBenchmark::dslashBenchmark(int argc, const char* argv[]) :
  benchmarkExecutable(argc, argv)
{
  // expect there to be exactly one device
  if(system->get_devices().size() != 1) {
    logger.fatal() << "There must be exactly one device chosen for the dslash benchmark to be performed.";
  }
  // expect that profiling is enabled
  if(! parameters.get_enable_profiling() )
    {
      throw Print_Error_Message( "Profiling is not enabled. Aborting...\n", __FILE__, __LINE__);
    }
  

}

void dslashBenchmark::performBenchmarkForSpecificKernels()
{

}
