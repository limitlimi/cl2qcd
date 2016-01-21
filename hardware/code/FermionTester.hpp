/*
 * Copyright 2015 Christopher Pinke
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

# pragma once

#include "SpinorTester.hpp"
#include "GaugefieldTester.hpp"
#include "fermions.hpp"

struct WilsonMassParameters
{
	WilsonMassParameters(const double kappaIn) : kappa(kappaIn){};
	const double kappa;
};

struct ThetaParameters
{
	ThetaParameters(const double thetaTIn, const double thetaSIn) : thetaT(thetaTIn), thetaS(thetaSIn){};
	const double thetaT, thetaS;
};

typedef hmc_complex ChemicalPotentials;

struct TwistedMassMassParameters : public WilsonMassParameters
{
	TwistedMassMassParameters(const double kappaIn, const double muIn):
		WilsonMassParameters(kappaIn), mu(muIn) {}

	const double mu;

	double getMubar() const { return 2.*kappa*mu; }
};

struct FermionTestParameters : public SpinorTestParameters, GaugefieldTestParameters
{
	FermionTestParameters(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypesIn) :
		TestParameters(lE),
		SpinorTestParameters(lE, SpinorFillTypes{spinorFillTypeIn}),
		GaugefieldTestParameters(lE, gaugefieldFillTypesIn) {};
	FermionTestParameters(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn) :
		TestParameters(lE),
		SpinorTestParameters(lE, SpinorFillTypes{spinorFillTypeIn}),
		GaugefieldTestParameters(lE, GaugefieldFillType::cold) {};
	FermionTestParameters(const LatticeExtents lE, const SpinorFillTypes spinorFillTypesIn) :
		TestParameters(lE),
		SpinorTestParameters(lE, spinorFillTypesIn),
		GaugefieldTestParameters(lE, GaugefieldFillType::cold) {};
};

template< class MassParameters>
 struct FermionMatrixTestParameters : public FermionTestParameters
{
	FermionMatrixTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn,
			const MassParameters massParametersIn) :
		TestParameters(latticeExtentsIn),
		FermionTestParameters(latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn), massParameters(massParametersIn) {};
	FermionMatrixTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const MassParameters massParametersIn ) :
		TestParameters(latticeExtentsIn),
		FermionTestParameters(latticeExtentsIn, spinorFillTypeIn, GaugefieldFillType::cold), massParameters(massParametersIn) {};
	const MassParameters massParameters;
};

typedef FermionMatrixTestParameters<WilsonMassParameters> WilsonTestParameters;
typedef FermionMatrixTestParameters<TwistedMassMassParameters> TwistedMassTestParameters;

template <class BufferType, class CreatorType>
struct FermionmatrixTester : public KernelTester
{
	FermionmatrixTester(std::string kernelName, const ParameterCollection parameterCollection, const FermionTestParameters testParameters, const ReferenceValues rV) :
		KernelTester(kernelName, parameterCollection.hardwareParameters, parameterCollection.kernelParameters, testParameters, rV)
	{
		GaugefieldCreator gf(testParameters.latticeExtents);
		gaugefieldBuffer = new hardware::buffers::SU3(calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), this->device);
		const Matrixsu3 * gf_host = gf.createGaugefield(testParameters.GaugefieldTestParameters::fillType);
		this->device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;

		code = this->device->getFermionCode();

		CreatorType sf(testParameters.latticeExtents);
		elements = sf.numberOfElements;
		in = new const BufferType(elements, this->device);
		out = new const BufferType(elements, this->device);
		out->load(sf.createSpinorfield(SpinorFillType::zero));
		in->load(sf.createSpinorfield(testParameters.SpinorTestParameters::fillTypes.at(0)));
	}
protected:
	const BufferType * in;
	const BufferType * out;
	const hardware::code::Fermions * code;
	const hardware::buffers::SU3 * gaugefieldBuffer;
	size_t elements;
};

struct NonEvenOddFermionmatrixTester : public FermionmatrixTester<hardware::buffers::Plain<spinor>,NonEvenOddSpinorfieldCreator>
{
	NonEvenOddFermionmatrixTester(std::string kernelName, const ParameterCollection parameterCollection, const FermionTestParameters testParameters, const ReferenceValues rV) :
		FermionmatrixTester<hardware::buffers::Plain<spinor>,NonEvenOddSpinorfieldCreator>(kernelName, parameterCollection, testParameters, rV) {}
};

struct EvenOddFermionmatrixTester : public FermionmatrixTester<hardware::buffers::Spinor,EvenOddSpinorfieldCreator>
{
	EvenOddFermionmatrixTester(std::string kernelName, const ParameterCollection parameterCollection, const FermionTestParameters testParameters, const ReferenceValues rV) :
		FermionmatrixTester<hardware::buffers::Spinor,EvenOddSpinorfieldCreator>(kernelName, parameterCollection, testParameters, rV) {}
};

template <typename TesterType >
struct FermionmatrixTesterWithSumAsKernelResult : public TesterType
{
	FermionmatrixTesterWithSumAsKernelResult( const std::string kernelName, const ParameterCollection parameterCollection, const FermionTestParameters testParameters, const ReferenceValues rV) :
		TesterType(kernelName, parameterCollection, testParameters, rV) {}
	~FermionmatrixTesterWithSumAsKernelResult()
	{
		spinor * sf_in;
		sf_in = new spinor[TesterType::elements];
		TesterType::out->dump(sf_in);
		TesterType::kernelResult.at(0) = count_sf(sf_in, TesterType::elements);
		delete sf_in;
	}
};
