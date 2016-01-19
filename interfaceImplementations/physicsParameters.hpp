/** @file
 * physicsParameters.hpp
 *
 * Copyright 2016 Alessandro Sciarra
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

#pragma once

#include "../meta/inputparameters.hpp"
#include "../meta/util.hpp"
#include "../physics/additionalParameters.hpp"
#include "../physics/fermionmatrix/fermionmatrixInterfaces.hpp"
#include "latticesParameters.hpp"
#include "fermionmatrixParameters.hpp"

namespace physics{

    class WilsonAdditionalParameters final : public AdditionalParameters {
        public:
            WilsonAdditionalParameters() = delete;
            WilsonAdditionalParameters(const meta::Inputparameters& paramsIn, const bool withMassPreconditioning)
            : parameters(paramsIn), withMassPreconditioning(withMassPreconditioning)
            {
            }
            virtual ~WilsonAdditionalParameters()
            {
            }
            hmc_float getKappa() const override
            {
                return withMassPreconditioning ? parameters.get_kappa_mp() : parameters.get_kappa();
            }
            hmc_float getMubar() const override
            {
                return withMassPreconditioning ? meta::get_mubar_mp(parameters) : meta::get_mubar(parameters);
            }

        private:
            const meta::Inputparameters& parameters;
            bool withMassPreconditioning;
    };

    class StaggeredAdditionalParameters final : public AdditionalParameters {
        public:
            StaggeredAdditionalParameters() = delete;
            StaggeredAdditionalParameters(const meta::Inputparameters& paramsIn) : parameters(paramsIn)
            {
            }
            virtual ~StaggeredAdditionalParameters()
            {
            }
            hmc_float getMass() const override
            {
                return parameters.get_mass();
            }
            bool getConservative() const override
            {
                return parameters.get_conservative();
            }

        private:
            const meta::Inputparameters& parameters;
    };


    class FermionParametersImplementation final : public FermionParametersInterface,
                                                  private lattices::SpinorfieldParametersImplementation,
                                                  private fermionmatrix::FermionmatrixParametersImplementation {
        public:
            FermionParametersImplementation() = delete;
            FermionParametersImplementation(const meta::Inputparameters& parametersIn)
                    : lattices::SpinorfieldParametersImplementation(parametersIn), fermionmatrix::FermionmatrixParametersImplementation(parametersIn)
            {
            }
            unsigned getNt() const override
            {
                return lattices::SpinorfieldParametersImplementation::getNt();
            }
            unsigned getNs() const override
            {
                return lattices::SpinorfieldParametersImplementation::getNs();
            }
            unsigned getNumberOfElements() const override
            {
                return lattices::SpinorfieldParametersImplementation::getNumberOfElements();
            }
            common::action getFermionicActionType() const override
            {
                return fermionmatrix::FermionmatrixParametersImplementation::getFermionicActionType();
            }
            bool useMergedFermionicKernels() const override
            {
                return fermionmatrix::FermionmatrixParametersImplementation::useMergedFermionicKernels();
            }
    };

    class FermionEoParametersImplementation final : public FermionEoParametersInterface,
                                                    private lattices::SpinorfieldEoParametersImplementation,
                                                    private fermionmatrix::FermionmatrixParametersImplementation {
        public:
            FermionEoParametersImplementation() = delete;
            FermionEoParametersImplementation(const meta::Inputparameters& parametersIn)
                    : lattices::SpinorfieldEoParametersImplementation(parametersIn), fermionmatrix::FermionmatrixParametersImplementation(parametersIn)
            {
            }
            common::action getFermionicActionType() const override
            {
                return fermionmatrix::FermionmatrixParametersImplementation::getFermionicActionType();
            }
            bool useMergedFermionicKernels() const override
            {
                return fermionmatrix::FermionmatrixParametersImplementation::useMergedFermionicKernels();
            }
    };

    class FermionStaggeredEoParametersImplementation: public FermionStaggeredEoParametersInterface,
                                                      private lattices::StaggeredfieldEoParametersImplementation,
                                                      private fermionmatrix::FermionmatrixStaggeredParametersImplementation {
        public:
            FermionStaggeredEoParametersImplementation() = delete;
            FermionStaggeredEoParametersImplementation(const meta::Inputparameters& parametersIn)
                    : lattices::StaggeredfieldEoParametersImplementation(parametersIn), fermionmatrix::FermionmatrixStaggeredParametersImplementation(parametersIn)
            {
            }
            unsigned getNumberOfElements() const override
            {
                return lattices::StaggeredfieldEoParametersImplementation::getNumberOfElements();
            }
    };

}
