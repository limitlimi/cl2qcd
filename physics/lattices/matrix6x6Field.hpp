/** @file
 * Declaration of the physics::lattices::Matrix6x6Field class
 *
 * Copyright 2012, 2013, 2015 Lars Zeidlewicz, Christopher Pinke,
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

#ifndef _PHYSICS_LATTICES_MATRIX6X6FIELD_
#define _PHYSICS_LATTICES_MATRIX6X6FIELD_

#include "../../hardware/system.hpp"
#include "../../hardware/buffers/6x6.hpp"
#include "../prng.hpp"
#include "scalar.hpp"
#include "latticesInterfaces.hpp"
#include "../../hardware/lattices/matrix6x6Field.hpp"
#include "gaugefield.hpp"

/**
 * This namespace contains the lattices of the various kind,
 * that is storage of the lattice values as a whole.
 */
namespace physics {
	namespace lattices {

		class Gaugefield;

		/**
		 * Representation of a Matrix6x6Field.
		 */
		class Matrix6x6Field {

		public:
			/**
			 * Construct a Matrix6x6Field based on the gaugefield configuration
			 */
			Matrix6x6Field(const hardware::System&, const GaugefieldParametersInterface * parameters);

			/**
			 * Release resources
			 */
			~Matrix6x6Field();

			/*
			 * Matrix6x6Field cannot be copied
			 */
			Matrix6x6Field& operator=(const Matrix6x6Field&) = delete;
			Matrix6x6Field(const Matrix6x6Field&) = delete;
			Matrix6x6Field() = delete;

			void setField(const physics::lattices::Gaugefield * gaugefield, const bool upperOrLower);

			/**
			 * Get the buffers containing the matrix6x6Field state on the devices.
			 */
			const std::vector<const hardware::buffers::matrix6x6 *> get_buffers() const noexcept;

			/**
			 * Update the halo cells of each buffer from its neighbours.
			 * On a single device this will be a no-op.
			 */
			void update_halo() const;

			const hardware::System * getSystem() const;
			const GaugefieldParametersInterface * getParameters() const;

		private:
			hardware::System const& system;
			const GaugefieldParametersInterface * latticeObjectParameters;
			hardware::lattices::Matrix6x6Field matrix6x6Field;

			friend hmc_float count_Matrix6x6Field(const Matrix6x6Field&);

			friend hmc_float log_det_Matrix6x6Field(const Matrix6x6Field&, const hmc_float kappa, const hmc_float csw);
		};

		hmc_float count_Matrix6x6Field(const Matrix6x6Field& field);

		hmc_float log_det_Matrix6x6Field(const Gaugefield& field, const hmc_float kappa, const hmc_float csw);
		void log_det_Matrix6x6Field(const Scalar<hmc_float>* res, const Gaugefield& field, const hmc_float kappa, const hmc_float csw);
	}
}



#endif /*_PHYSICS_LATTICES_MATRIX6X6FIELD_ */
