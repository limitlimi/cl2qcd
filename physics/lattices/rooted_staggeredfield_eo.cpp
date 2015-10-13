/** @file
 * Implementation of the physics::lattices::Rooted_Staggeredfield_eo class
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#include <algorithm>    // std::max
#include "rooted_staggeredfield_eo.hpp"

physics::lattices::Rooted_Staggeredfield_eo::Rooted_Staggeredfield_eo(const hardware::System& system)
	: Staggeredfield_eo(system), rootedStaggaredfieldEoParametersInterface(new RootedStaggaredfieldEoParametersImplementation(system.get_inputparameters())),
	  physics::algorithms::Rational_Coefficients(std::max(rootedStaggaredfieldEoParametersInterface->getMetropolisRationalApproximationOrder(),
	                                                      rootedStaggaredfieldEoParametersInterface->getMolecularDynamicsRationalApproximationOrder()))
{
}

physics::lattices::Rooted_Staggeredfield_eo::Rooted_Staggeredfield_eo(const physics::algorithms::Rational_Approximation& approx, const hardware::System& system)
	: Staggeredfield_eo(system), rootedStaggaredfieldEoParametersInterface(new RootedStaggaredfieldEoParametersImplementation(system.get_inputparameters())),
	  physics::algorithms::Rational_Coefficients(approx.Get_order(), approx.Get_a0(), approx.Get_a(), approx.Get_b())
{
}

physics::lattices::Rooted_Staggeredfield_eo::~Rooted_Staggeredfield_eo()
{
    //TODO: remove the following delete
    delete rootedStaggaredfieldEoParametersInterface;
}

void physics::lattices::Rooted_Staggeredfield_eo::Rescale_Coefficients(const physics::algorithms::Rational_Approximation& approx, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf, const hardware::System& system, hmc_float prec, bool conservative)
{
	physics::algorithms::Rational_Coefficients aux = approx.Rescale_Coefficients(A, gf, system, prec, conservative);
	
	Set_coeff(aux.Get_a0(), aux.Get_a(), aux.Get_b());
}


