/** @file
 * Declaration of the physics::lattices::wilson:Rooted_Spinorfield_eo class
 *
 * Copyright (c) 2017 Christopher Czaban <czaban@th.physik.uni-frankfurt.de>
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

#ifndef _PHYSICS_LATTICES_ROOTED_SPINORFIELD_EO_
#define _PHYSICS_LATTICES_ROOTED_SPINORFIELD_EO_


#include "../../hardware/system.hpp"
#include "../algorithms/rational_approximation.hpp"
#include "spinorfield_eo.hpp"
#include "util.hpp" //This is to make the template pseudo_randomize friend of this class

namespace physics {
	namespace lattices {
		namespace wilson {
			//Representation of a rooted spinorfield with eo preconditioning
			class Rooted_Spinorfield_eo: public physics::lattices::Spinorfield_eo {

			public:
				Rooted_Spinorfield_eo(const hardware::System&, const RootedSpinorfieldEoParametersInterface&);
				Rooted_Spinorfield_eo(const hardware::System&, const RootedSpinorfieldEoParametersInterface&, const physics::algorithms::Rational_Approximation& approx);

				~Rooted_Spinorfield_eo(){};

				//Rescale coefficients on the basis of a Rational_Approximation objects
				void Rescale_Coefficients(const physics::algorithms::Rational_Approximation& approx, const hmc_float minEigenvalue, const hmc_float maxEigenvalue);

				Rooted_Spinorfield_eo& operator=(const Rooted_Spinorfield_eo&) = delete;
				Rooted_Spinorfield_eo(const Rooted_Spinorfield_eo&) = delete;
				Rooted_Spinorfield_eo() = delete;

				int Get_order() const;
				hmc_float Get_a0() const;
				std::vector<hmc_float> Get_a() const;
				std::vector<hmc_float> Get_b() const;

			private:
				physics::algorithms::Rational_Coefficients rationalCoefficients;

//				friend void pseudo_randomize<Rooted_Spinorfield_eo, spinor>(const Rooted_Spinorfield_eo* to, int seed);
			};
		}
	}
}

#endif /*_PHYSICS_LATTICES_ROOTED_SPINORFIELD_EO_ */
