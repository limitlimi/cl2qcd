/*
 * Copyright 2012, 2013, 2014 Christopher Pinke, Matthias Bach
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HARDWARE_CODE_FERMIONS

#include "FermionTester.hpp"

/** todo:
 * - the way the reference values are passed is not nice, it should be somehow done automatically.
 * 		Especially that one has to call the lattice volume explicitely is not nice
 */
const ReferenceValues calculateReferenceValues_mWilson(const int latticeVolume, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn)
{
	if( massParametersIn.kappa == 0. and spinorFillTypeIn == SpinorFillType::ascendingComplex)
	{
		return ReferenceValues{latticeVolume * 300.};
	}
	else if (massParametersIn.kappa == nonTrivialParameter)
	{
		if (spinorFillTypeIn == SpinorFillType::ascendingComplex and gaugefieldFillTypeIn == GaugefieldFillType::cold )
		{
			return ReferenceValues{latticeVolume * 3.705600000000023 };
		}
		else if (spinorFillTypeIn == SpinorFillType::one and gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial )
		{
			return ReferenceValues{latticeVolume * 4.590333841920002};
		}
		else if (spinorFillTypeIn == SpinorFillType::ascendingComplex and gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial)
		{
			return ReferenceValues{latticeVolume * 197.113059468288};
		}
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_mTmMinus(const int latticeVolume, const SpinorFillType spinorFillTypeIn,
		const GaugefieldFillType gaugefieldFillTypeIn, const TwistedMassMassParameters massParametersIn)
{
	if( massParametersIn.kappa == nonTrivialParameter and spinorFillTypeIn == SpinorFillType::ascendingComplex )
	{
		if( massParametersIn.mu == 0. )
		{
			return calculateReferenceValues_mWilson(latticeVolume, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn);
		}
		else if (massParametersIn.mu == nonTrivialParameter and gaugefieldFillTypeIn == GaugefieldFillType::cold )
		{
			return ReferenceValues{latticeVolume * 3.705600000000016};
		}
		else if (massParametersIn.mu == nonTrivialParameter and gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial )
		{
			return ReferenceValues{latticeVolume * 197.113059468288};
		}
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_mTmSitediagonal(const int latticeVolume, const TwistedMassMassParameters massParametersIn)
{
	if (massParametersIn.kappa == nonTrivialParameter and massParametersIn.mu == nonTrivialParameter)
	{
		return ReferenceValues{ latticeVolume * 300.};
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_mTmSitediagonalMinus(const int latticeVolume, const TwistedMassMassParameters massParametersIn)
{
	if (massParametersIn.kappa == nonTrivialParameter and massParametersIn.mu == nonTrivialParameter)
	{
		return ReferenceValues{ latticeVolume * 299.9999999999999};
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_mTmInverseSitediagonal(const int latticeVolume, const TwistedMassMassParameters massParametersIn)
{
	if (massParametersIn.kappa == nonTrivialParameter and massParametersIn.mu == nonTrivialParameter)
	{
		return ReferenceValues{ latticeVolume * 299.7214990417087};
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_mTmInverseSitediagonalMinus(const int latticeVolume, const TwistedMassMassParameters massParametersIn)
{
	if (massParametersIn.kappa == nonTrivialParameter and massParametersIn.mu == nonTrivialParameter)
	{
		return ReferenceValues{ latticeVolume * 299.7214990417086};
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValuesDslashEvenOdd(LatticeExtents latticeExtentsIn, const int latticeVolumeIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn, const double thetaT, const double thetaS, const double chemPotImIn)
{
	if (massParametersIn.kappa == nonTrivialParameter and spinorFillTypeIn == SpinorFillType::ascendingComplex)
	{
		if (thetaS == 0. and thetaT == 0.)
		{
			if (chemPotImIn == 0.)
			{
				if (gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial )
				{
					return ReferenceValues{ latticeVolumeIn * -102.886940531712};
				}
				else if (gaugefieldFillTypeIn == GaugefieldFillType::cold )
				{
					return ReferenceValues{ latticeVolumeIn * -296.2943999999999};
				}
			}
			else
			{
				if (gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial )
				{
					return ReferenceValues{ latticeVolumeIn * -97.13406449727512};
				}
				else if (gaugefieldFillTypeIn == GaugefieldFillType::cold )
				{
					return ReferenceValues{ latticeVolumeIn * -295.3657597941546};
				}
			}
		}
		else
		{
			if (gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial )
			{
				return ReferenceValues{ latticeVolumeIn * (-72.61889049538556 * cos(M_PI*thetaS/latticeExtentsIn.ns) + 144.307369760256 * sin(M_PI*thetaS/latticeExtentsIn.ns) - 30.2680500363264 * cos(M_PI*thetaT/latticeExtentsIn.nt) + 44.84641322803199  * sin(M_PI*thetaT/latticeExtentsIn.nt))};
			}
			else if (gaugefieldFillTypeIn == GaugefieldFillType::cold )
			{
				return ReferenceValues{ latticeVolumeIn * (-222.2208 * cos(M_PI*thetaS/latticeExtentsIn.ns) + 35.55532799999998 * sin(M_PI*thetaS/latticeExtentsIn.ns) - 74.07360000000003 * cos(M_PI*thetaT/latticeExtentsIn.nt) + 2.962943999999992  * sin(M_PI*thetaT/latticeExtentsIn.nt))};
			}
		}
	}
	return defaultReferenceValues();
}

ReferenceValues calculateReferenceValues_gamma5(const int latticeVolume, const SpinorFillType spinorFillTypeIn)
{
	return (spinorFillTypeIn == SpinorFillType::ascendingComplex) ? ReferenceValues{latticeVolume*( sumOfIntegers(1, 12, 1) - sumOfIntegers(13, 24, 1) )} : defaultReferenceValues() ;
}

struct DslashEvenOddTestParameters: public WilsonTestParameters
{
	DslashEvenOddTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn) :
		TestParameters(latticeExtentsIn), WilsonTestParameters( latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn),
		massParameters(massParametersIn), thetaT(0.), thetaS(0.), chemPotIm(0.) {};
	DslashEvenOddTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn, const double thetaSIn, const double thetaTIn) :
		TestParameters(latticeExtentsIn), WilsonTestParameters( latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn),
		massParameters(massParametersIn), thetaT(thetaSIn), thetaS(thetaTIn), chemPotIm(0.) {};
	DslashEvenOddTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn, const double chemPotImIn) :
		TestParameters(latticeExtentsIn), WilsonTestParameters( latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn),
		massParameters(massParametersIn), thetaT(0.), thetaS(0.), chemPotIm(chemPotImIn) {};
	const WilsonMassParameters massParameters;
	double thetaT, thetaS;
	double chemPotIm;
};

struct MWilsonTester : public FermionmatrixTesterWithSumAsKernelResult<NonEvenOddFermionmatrixTester>
{
	MWilsonTester(const ParameterCollection parameterCollection, const WilsonTestParameters & tP) :
		FermionmatrixTesterWithSumAsKernelResult<NonEvenOddFermionmatrixTester>("m_wilson", parameterCollection, tP, calculateReferenceValues_mWilson( getSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.fillTypes.at(0), tP.fillType, tP.massParameters) )
	{
		code->M_wilson_device(in, out,  gaugefieldBuffer, tP.massParameters.kappa );
	}
};
struct MTmMinusTester : public FermionmatrixTesterWithSumAsKernelResult<NonEvenOddFermionmatrixTester>
{
	MTmMinusTester(const ParameterCollection parameterCollection, const TwistedMassTestParameters & tP) :
		FermionmatrixTesterWithSumAsKernelResult<NonEvenOddFermionmatrixTester>("m_tm_minus", parameterCollection, tP, calculateReferenceValues_mTmMinus( getSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.fillTypes.at(0), tP.fillType, tP.massParameters))
	{
		code->M_tm_minus_device(in, out,  gaugefieldBuffer, tP.massParameters.kappa, tP.massParameters.getMubar() );
	}
};
struct MTmPlusTester : public FermionmatrixTesterWithSumAsKernelResult<NonEvenOddFermionmatrixTester>
{
	MTmPlusTester(const ParameterCollection parameterCollection, const TwistedMassTestParameters & tP) :
		FermionmatrixTesterWithSumAsKernelResult<NonEvenOddFermionmatrixTester>("m_tm_plus", parameterCollection, tP, calculateReferenceValues_mTmMinus( getSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.fillTypes.at(0), tP.fillType, tP.massParameters))
	{
		code->M_tm_plus_device(in, out,  gaugefieldBuffer, tP.massParameters.kappa, tP.massParameters.getMubar() );
	}
};

template< class bufferType, typename TesterType>
struct Gamma5Tester : public TesterType
{
	Gamma5Tester( const std::string kernelName, const ParameterCollection parameterCollection, const FermionTestParameters & tP, const ReferenceValues rV, const int numberOfElementsIn) :
		TesterType(kernelName, parameterCollection, tP, rV), numberOfElements (numberOfElementsIn)
	{
		sf_in = new spinor[numberOfElements];
		in = new bufferType(numberOfElements, this->device);

		in->load( TesterType::createSpinorfield(tP.SpinorTestParameters::fillTypes.at(0)) );
	}
	~Gamma5Tester()
	{
		in->dump(sf_in);
		TesterType::kernelResult.at(0) = TesterType::count_sf(sf_in, numberOfElements);
		delete sf_in;
	}
protected:
	const int numberOfElements;
	spinor * sf_in;
	const bufferType * in;
};

struct Gamma5NonEvenOddTester : public Gamma5Tester<hardware::buffers::Plain<spinor>, NonEvenOddFermionmatrixTester >
{
	Gamma5NonEvenOddTester(const ParameterCollection parameterCollection, const FermionTestParameters & tP) :
		Gamma5Tester<hardware::buffers::Plain<spinor>, NonEvenOddFermionmatrixTester>("gamma5", parameterCollection, tP,
				calculateReferenceValues_gamma5(calculateSpinorfieldSize( tP.SpinorTestParameters::latticeExtents ), tP.fillTypes.at(0) ), calculateSpinorfieldSize( tP.SpinorTestParameters::latticeExtents ) )
	{
		code->gamma5_device(in);
	}
};

struct Gamma5EvenOddTester : public Gamma5Tester<hardware::buffers::Spinor, EvenOddFermionmatrixTester>
{
	Gamma5EvenOddTester(const ParameterCollection parameterCollection, const FermionTestParameters & tP) :
		Gamma5Tester<hardware::buffers::Spinor, EvenOddFermionmatrixTester>("gamma5_eo", parameterCollection, tP,
				calculateReferenceValues_gamma5(calculateEvenOddSpinorfieldSize( tP.SpinorTestParameters::latticeExtents ), tP.fillTypes.at(0) ), calculateEvenOddSpinorfieldSize( tP.SpinorTestParameters::latticeExtents ) )
		{
			code->gamma5_eo_device(in);
		}
};

//todo: remove ARG_DEF from all the tm diagonal kernel fcts.!
struct MTmSitediagonalTester: public FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>
{
	MTmSitediagonalTester(const ParameterCollection parameterCollection, const TwistedMassTestParameters & tP):
		FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>("m_tm_sitediagonal", parameterCollection, tP,
				calculateReferenceValues_mTmSitediagonal(getEvenOddSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.massParameters))
		{
			code->M_tm_sitediagonal_device( in, out, tP.massParameters.getMubar());
		}
};

struct MTmInverseSitediagonalTester: public FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>
{
	MTmInverseSitediagonalTester(const ParameterCollection parameterCollection, const TwistedMassTestParameters & tP):
		FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>("m_tm_inverse_sitediagonal", parameterCollection, tP,
				calculateReferenceValues_mTmInverseSitediagonal(getEvenOddSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.massParameters))
		{
			code->M_tm_inverse_sitediagonal_device( in, out, tP.massParameters.getMubar());
		}
};
struct MTmSitediagonalMinusTester: public FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>
{
	MTmSitediagonalMinusTester(const ParameterCollection parameterCollection, const TwistedMassTestParameters & tP):
		FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>("m_tm_sitediagonal_minus", parameterCollection, tP,
				calculateReferenceValues_mTmSitediagonalMinus(getEvenOddSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.massParameters))
		{
			code->M_tm_sitediagonal_minus_device( in, out, tP.massParameters.getMubar());
		}
};

struct MTmInverseSitediagonalMinusTester: public FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>
{
	MTmInverseSitediagonalMinusTester(const ParameterCollection parameterCollection, const TwistedMassTestParameters & tP):
		FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>("m_tm_inverse_sitediagonal", parameterCollection, tP,
				calculateReferenceValues_mTmInverseSitediagonalMinus(getEvenOddSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.massParameters))
		{
			code->M_tm_inverse_sitediagonal_minus_device( in, out, tP.massParameters.getMubar());
		}
};

struct DslashEvenOddTester: public FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>
{
	DslashEvenOddTester(const ParameterCollection parameterCollection, const DslashEvenOddTestParameters & tP, const bool evenOrOddIn):
		FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>("dslash_eo", parameterCollection, tP,
				calculateReferenceValuesDslashEvenOdd(tP.latticeExtents, getEvenOddSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.fillTypes.at(0), tP.fillType, tP.massParameters, tP.thetaT, tP.thetaS, tP.chemPotIm))
		{
			evenOrOddIn ?
				code->dslash_eo_device( in, out,  gaugefieldBuffer, EVEN, tP.massParameters.kappa) :
				code->dslash_eo_device( in, out,  gaugefieldBuffer, ODD, tP.massParameters.kappa );
		}
};

template<typename TesterClass, typename MassParameters, typename TestParameters, typename KernelParameterMockup>
void callTest( const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn,
		const MassParameters massParametersIn, const bool needEvenOdd)
{
	TestParameters parametersForThisTest(latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn);
	KernelParameterMockup kernelParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, needEvenOdd); //todo: could also use latticeExtendsIn here!
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, needEvenOdd);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass tester(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void callTest(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const TwistedMassMassParameters massParameters, const bool needEvenOdd)
{
	TwistedMassTestParameters parametersForThisTest(latticeExtentsIn, spinorFillTypeIn, massParameters);
	hardware::code::OpenClKernelParametersMockupForTwistedMass kernelParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, needEvenOdd);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, needEvenOdd);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass tester(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void callTest(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const bool needEvenOdd)
{
	FermionTestParameters parametersForThisTest(latticeExtentsIn, SpinorFillTypes{spinorFillTypeIn} );
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, needEvenOdd);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, needEvenOdd);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass tester(parameterCollection, parametersForThisTest);
}

void testMWilson(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn,
		const WilsonMassParameters massParametersIn)
{
	callTest<MWilsonTester, WilsonMassParameters, WilsonTestParameters, hardware::code::OpenClKernelParametersMockupForSpinorTests>
		(latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn, false);
}

void testMTmMinus(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn,
		const TwistedMassMassParameters massParametersIn)
{
	callTest<MTmMinusTester, TwistedMassMassParameters, TwistedMassTestParameters, hardware::code::OpenClKernelParametersMockupForTwistedMass>
		(latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn, false);
}

//todo: test should not pass using M_tm_minus ref. values!
void testMTmPlus(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn,
		const TwistedMassMassParameters massParametersIn)
{
	callTest<MTmPlusTester, TwistedMassMassParameters, TwistedMassTestParameters, hardware::code::OpenClKernelParametersMockupForTwistedMass>
		(latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn, false);
}

void testNonEvenOddGamma5(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn)
{
	callTest<Gamma5NonEvenOddTester>(latticeExtentsIn, spinorFillTypeIn, false);
}

void testEvenOddGamma5(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn)
{
	callTest<Gamma5EvenOddTester>(latticeExtentsIn, spinorFillTypeIn, true);
}

void testMTmSitediagonal(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const TwistedMassMassParameters massParametersIn)
{
	callTest<MTmSitediagonalTester>(latticeExtentsIn, spinorFillTypeIn, massParametersIn, true);
}

void testMTmInverseSitediagonal(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const TwistedMassMassParameters massParametersIn)
{
	callTest<MTmInverseSitediagonalTester>(latticeExtentsIn, spinorFillTypeIn, massParametersIn, true);
}

void testMTmSitediagonalMinus(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const TwistedMassMassParameters massParametersIn)
{
	callTest<MTmSitediagonalMinusTester>(latticeExtentsIn, spinorFillTypeIn, massParametersIn, true);
}

void testMTmInverseSitediagonalMinus(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const TwistedMassMassParameters massParametersIn)
{
	callTest<MTmInverseSitediagonalMinusTester>(latticeExtentsIn, spinorFillTypeIn, massParametersIn, true);
}

void testDslashEvenOdd(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn,
		const WilsonMassParameters massParametersIn, const bool evenOrOddIn)
{
	DslashEvenOddTestParameters parametersForThisTest(latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, true);
	hardware::code::OpenClKernelParametersMockupForDslashEvenOdd kernelParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	DslashEvenOddTester tester(parameterCollection, parametersForThisTest, evenOrOddIn);
}

void testDslashEvenOddWithSpecificBC(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn,
		const WilsonMassParameters massParametersIn, const double thetaTIn, const double thetaSIn, const bool evenOrOddIn)
{
	DslashEvenOddTestParameters parametersForThisTest(latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn, thetaTIn, thetaSIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, true);
	hardware::code::OpenClKernelParametersMockupForDslashEvenOdd kernelParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, true, parametersForThisTest.thetaT, parametersForThisTest.thetaS);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	DslashEvenOddTester tester(parameterCollection, parametersForThisTest, evenOrOddIn);
}

void testDslashEvenOddWithChemPotIm(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn,
		const WilsonMassParameters massParametersIn, const double chemPotImIn, const bool evenOrOddIn)
{
	DslashEvenOddTestParameters parametersForThisTest(latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn, chemPotImIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, true);
	hardware::code::OpenClKernelParametersMockupForDslashEvenOdd kernelParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, true, true, parametersForThisTest.chemPotIm);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	DslashEvenOddTester tester(parameterCollection, parametersForThisTest, evenOrOddIn);
}

BOOST_AUTO_TEST_SUITE( M_WILSON )

	BOOST_AUTO_TEST_CASE( M_WILSON_1)
	{
		testMWilson(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{0.});
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_2)
	{
		testMWilson(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialParameter});
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_3)
	{
		testMWilson(LatticeExtents{ns4, nt4}, SpinorFillType::one, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter});
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_4)
	{
		testMWilson(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_MINUS  )

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_1 )
	{
		testMTmMinus(LatticeExtents{ns12, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{nonTrivialParameter, 0.});
	}

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_2 )
	{
		testMTmMinus(LatticeExtents{ns4, nt16}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{nonTrivialParameter, nonTrivialParameter});
	}

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_3 )
	{
		testMTmMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{nonTrivialParameter, nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_PLUS )

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_1 )
	{
		testMTmPlus(LatticeExtents{ns12, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{nonTrivialParameter, 0.});
	}

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_2 )
	{
		testMTmPlus(LatticeExtents{ns4, nt16}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{nonTrivialParameter, nonTrivialParameter});
	}

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_3 )
	{
		testMTmPlus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{nonTrivialParameter, nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GAMMA5 )

	BOOST_AUTO_TEST_CASE( GAMMA5_NONEO )
	{
		testNonEvenOddGamma5( LatticeExtents{ns4,nt4}, SpinorFillType::ascendingComplex);
	}

	BOOST_AUTO_TEST_CASE( GAMMA5_EO )
	{
		testEvenOddGamma5( LatticeExtents{ns4,nt4}, SpinorFillType::ascendingComplex);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL )

	BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_1)
	{
		testMTmSitediagonal(LatticeExtents{ns4,nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{nonTrivialParameter, nonTrivialParameter});
	}

	BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_1)
	{
		testMTmInverseSitediagonal(LatticeExtents{ns4,nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{nonTrivialParameter, nonTrivialParameter});
	}

	BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_1)
	{
		testMTmSitediagonalMinus(LatticeExtents{ns4,nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{nonTrivialParameter, nonTrivialParameter} );
	}

	BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_1)
	{
		testMTmInverseSitediagonalMinus(LatticeExtents{ns4,nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{nonTrivialParameter, nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_EO )

	/**
	 * @todo: Missing tests:
	 * - chemical potential (real, imaginary) (need test with dummy activation?)
	 * Think about if one can automatize the evenOrOdd tests, they do not influence the reference values (?)
	 */

	BOOST_AUTO_TEST_CASE( DSLASH_EO_1)
	{
		testDslashEvenOdd(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialParameter}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_2)
	{
		testDslashEvenOdd(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialParameter}, false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_3)
	{
		testDslashEvenOdd(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_4)
	{
		testDslashEvenOdd(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter}, false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_ABC_1)
	{
		testDslashEvenOddWithSpecificBC(LatticeExtents{ns8,nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialParameter}, 1., 0., true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_ABC_2)
	{
		testDslashEvenOddWithSpecificBC(LatticeExtents{ns8,nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialParameter}, 1., 1., false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_ABC_3)
	{
		testDslashEvenOddWithSpecificBC(LatticeExtents{ns8,nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter}, 1., 0., true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_ABC_4)
	{
		testDslashEvenOddWithSpecificBC(LatticeExtents{ns12,nt12}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter}, 1., 1., false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_CHEMPOTIM_1)
	{
		testDslashEvenOddWithChemPotIm(LatticeExtents{ns4,nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialParameter}, nonTrivialParameter, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_CHEMPOTIM_2)
	{
		testDslashEvenOddWithChemPotIm(LatticeExtents{ns8,nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialParameter}, nonTrivialParameter, false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_CHEMPOTIM_3)
	{
		testDslashEvenOddWithChemPotIm(LatticeExtents{ns12,nt12}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter}, nonTrivialParameter, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_CHEMPOTIM_4)
	{
		testDslashEvenOddWithChemPotIm(LatticeExtents{ns4,nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter}, nonTrivialParameter, false);
	}

BOOST_AUTO_TEST_SUITE_END()

