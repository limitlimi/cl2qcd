/** @file
 * Declaration of the fermion_force functions
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#ifndef _PHYSICS_ALGORITHMS_FERMION_FORCE_
#define _PHYSICS_ALGORITHMS_FERMION_FORCE_

#include "../lattices/gaugemomenta.hpp"
#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "../lattices/matrix6x6Field.hpp"
#include "../interfacesHandler.hpp"

namespace physics {
    namespace algorithms {

        //These methods really calculate the total fermion force and they add it to the Gaugemomenta field
        void calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
                                 const physics::lattices::Spinorfield& phi, const hardware::System& system,
                                 physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters);
        void calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
                                 const physics::lattices::Spinorfield_eo& phi, const hardware::System& system,
                                 physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters);

        void calc_detratio_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
                                  const physics::lattices::Spinorfield& phi, const hardware::System& system, physics::InterfacesHandler& interfacesHandler);
        void calc_detratio_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
                                  const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfacesHandler);

        //Here, in the following functions, there is the detailed force calculation (these functions
        //are called from those above, that are actually unified by a template in the .cpp file)
        void calc_fermion_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
                                const physics::lattices::Spinorfield& phi, const hardware::System& system,
                                physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters);
        void calc_fermion_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
                                const physics::lattices::Spinorfield_eo& phi, const hardware::System& system,
                                physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters);

        void calc_fermion_force_detratio(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
                                         const physics::lattices::Spinorfield& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfacesHandler);
        void calc_fermion_force_detratio(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
                                         const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfacesHandler);

        //These methods interfaces only the lower level of the code (Molecular_Dynamics class) with the upper one,
        //namely they just call the function that enqueues the kernel
        void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Spinorfield& Y, const physics::lattices::Spinorfield& X,
                           const physics::lattices::Gaugefield& gf, const physics::AdditionalParameters& additionalParameters);
        void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Spinorfield_eo& Y, const physics::lattices::Spinorfield_eo& X,
                           int evenodd, const physics::lattices::Gaugefield& gf, const physics::AdditionalParameters& additionalParameters);
        void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Spinorfield_eo& Y, const physics::lattices::Spinorfield_eo& X,
        				   int evenodd, const physics::lattices::Gaugefield& gf, const physics::lattices::Matrix6x6Field& C, const physics::lattices::Matrix6x6Field& D, const physics::AdditionalParameters& additionalParameters);

        //Here the Kernels for the explizit version of the inverted clover fermionmatrix are called
        //note: the fermionmatrix is blockdiagonal, so one stores 2 6x6 blocks
        void clover_eo_inverse_explizit_upper_left(const physics::lattices::Matrix6x6Field& in, const physics::lattices::Matrix6x6Field * out, const physics::lattices::Gaugefield& gf, const physics::AdditionalParameters& additionalParameters);
        void clover_eo_inverse_explizit_lower_right(const physics::lattices::Matrix6x6Field& in, const physics::lattices::Matrix6x6Field * out, const physics::lattices::Gaugefield& gf, const physics::AdditionalParameters& additionalParameters);
    }
}

#endif /* _PHYSICS_ALGORITHMS_FERMION_FORCE_ */
