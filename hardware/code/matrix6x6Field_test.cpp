/*
 * Copyright 2015 Max Theilig
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
#define BOOST_TEST_MODULE hardware::code::Gaugefield
#include <boost/test/unit_test.hpp>

#include "Matrix6x6FieldTester.hpp"


struct CloverEvenOddInverseExplizitUpperLeftTester : public Matrix6x6FieldTester
{
	CloverEvenOddInverseExplizitUpperLeftTester(const ParameterCollection & parameterCollection, const Matrix6x6FieldTestParameters testParams, const ReferenceValues rV):
		Matrix6x6FieldTester("clover_eo_inverse_explizit_upper_left", parameterCollection, testParams, rV)
	{
		GaugefieldCreator gf(testParams.latticeExtents);
		matrix6x6FieldBuffer = new hardware::buffers::matrix6x6(testParams.latticeExtents, this->device);
		gaugefieldBuffer = new hardware::buffers::SU3(testParams.latticeExtents, this->device);
		const Matrixsu3 * gf_host = gf.createGaugefield(testParams.fillType);
		device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;

		code->clover_eo_inverse_explizit_upper_left_device(matrix6x6FieldBuffer, gaugefieldBuffer, parameterCollection.kernelParameters.getKappa(), parameterCollection.kernelParameters.getCsw());
	}
protected:
	const hardware::buffers::matrix6x6 * matrix6x6FieldBuffer;
	const hardware::buffers::SU3 * gaugefieldBuffer;
};

struct CloverEvenOddInverseExplizitLowerRightTester : public Matrix6x6FieldTester
{
	CloverEvenOddInverseExplizitLowerRightTester(const ParameterCollection & parameterCollection, const Matrix6x6FieldTestParameters testParams, const ReferenceValues rV):
		Matrix6x6FieldTester("clover_eo_inverse_explizit_lower_right", parameterCollection, testParams, rV)
	{
		GaugefieldCreator gf(testParams.latticeExtents);
		matrix6x6FieldBuffer = new hardware::buffers::matrix6x6(testParams.latticeExtents, this->device);
		gaugefieldBuffer = new hardware::buffers::SU3(testParams.latticeExtents, this->device);
		const Matrixsu3 * gf_host = gf.createGaugefield(testParams.fillType);
		device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;

		code->clover_eo_inverse_explizit_lower_right_device(matrix6x6FieldBuffer, gaugefieldBuffer, parameterCollection.kernelParameters.getKappa(), parameterCollection.kernelParameters.getCsw());
	}
protected:
	const hardware::buffers::matrix6x6 * matrix6x6FieldBuffer;
	const hardware::buffers::SU3 * gaugefieldBuffer;
};

template<typename TesterClass>
void performTest(const ReferenceValues refValuesIn, const LatticeExtents latticeExtentsIn, const GaugefieldFillType fillType)
{
	Matrix6x6FieldTestParameters parametersForThisTest {latticeExtentsIn, fillType};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	TesterClass tester(parameterCollection, parametersForThisTest, refValuesIn);
}

void testCloverEvenOddInverseExplizitUpperLeft(const ReferenceValues rV, const LatticeExtents lE, const GaugefieldFillType fT)
{
	performTest<CloverEvenOddInverseExplizitUpperLeftTester>(rV, lE, fT);
}

void testCloverEvenOddInverseExplizitLowerRight(const ReferenceValues rV, const LatticeExtents lE, const GaugefieldFillType fT)
{
	performTest<CloverEvenOddInverseExplizitLowerRightTester>(rV, lE, fT);
}


BOOST_AUTO_TEST_SUITE ( CLOVER_EO_INVERSE_EXPLIZIT )

	BOOST_AUTO_TEST_CASE( CLOVER_EO_INVERSE_EXPLIZIT_UPPER_LEFT )
	{
		testCloverEvenOddInverseExplizitUpperLeft(ReferenceValues{1234.}, LatticeExtents{ns4, nt4}, GaugefieldFillType::cold );
	}

	BOOST_AUTO_TEST_CASE( CLOVER_EO_INVERSE_EXPLIZIT_LOWER_RIGHT )
	{
		testCloverEvenOddInverseExplizitLowerRight(ReferenceValues{1234.}, LatticeExtents{ns4, nt4}, GaugefieldFillType::cold );
	}

BOOST_AUTO_TEST_SUITE_END()
