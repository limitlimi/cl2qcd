/** @file
 * Declaration of the physics::lattices::Rooted_Staggeredfield_eo class
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

#ifndef _PHYSICS_LATTICES_ROOTED_STAGGEREDFIELD_EO_
#define _PHYSICS_LATTICES_ROOTED_STAGGEREDFIELD_EO_

#include "../../hardware/system.hpp"
#include "../algorithms/rational_approximation.hpp"
#include "staggeredfield_eo.hpp"
#include "util.hpp" //This is to make the template pseudo_randomize friend of this class

/**
 * This namespace contains the lattices of the various kind,
 * that is storage of the lattice values as a whole.
 */
namespace physics {
	namespace lattices {

		/**
		 * Representation of a rooted staggeredfield (with eo preconditioning).
		 */
		class Rooted_Staggeredfield_eo : public Staggeredfield_eo, public physics::algorithms::Rational_Coefficients{

			public:
				/**
				 * Construct a rooted staggeredfield based on the input-files of the system
				 */
				Rooted_Staggeredfield_eo(const hardware::System&, const RootedStaggeredfieldEoParametersInterface&);
				Rooted_Staggeredfield_eo(const hardware::System&, const RootedStaggeredfieldEoParametersInterface&, const physics::algorithms::Rational_Approximation& approx);

				virtual ~Rooted_Staggeredfield_eo(){}

			/**
			 	 * Rescale coefficients on the basis of a Rational_Approximation objects
			 */
				void Rescale_Coefficients(const physics::algorithms::Rational_Approximation& approx, const hmc_float minEigenvalue, const hmc_float maxEigenvalue);
				/**
				 * Staggeredfield_eo cannot be copied
				 */
				Rooted_Staggeredfield_eo& operator=(const Rooted_Staggeredfield_eo&) = delete;
				Rooted_Staggeredfield_eo(const Rooted_Staggeredfield_eo&) = delete;
				Rooted_Staggeredfield_eo() = delete;

				friend void pseudo_randomize<Rooted_Staggeredfield_eo, su3vec>(const Rooted_Staggeredfield_eo* to, int seed);
		};

	}
}

#endif /*_PHYSICS_LATTICES_ROOTED_STAGGEREDFIELD_EO_ */

