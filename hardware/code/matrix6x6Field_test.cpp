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

/*struct CloverEvenOddInverseTestParameters : public Matrix6x6FieldTestParameters
{
	CloverEvenOddInverseTestParameters(const LatticeExtents latticeExtentsIn, Matrix6x6FieldFillType fillTypeIn,) :
		TestParameters(latticeExtentsIn), Matrix6x6FieldTestParameters(latticeExtentsIn, fillTypeIn)
	{}
};*/

struct CloverEvenOddInverseExplizitUpperLeftTester : Matrix6x6FieldTester
{
	CloverEvenOddInverseExplizitUpperLeftTester(const ParameterCollection & parameterCollection, const Matrix6x6FieldTestParameters testParams, const ReferenceValues rV):
		Matrix6x6FieldTester("clover_eo_inverse_explizit_upper_left", parameterCollection, testParams, rv)
	{
		Matrix6x6FieldCreator C(testParams.latticeExtents);
		matrix6x6FieldBuffer = new hardware::buffers::matrix6x6(testParams.latticeExtents, this->device);
		const Matrix6x6 * C_host = C.createMatrix6x6Field(testParams.fillType);
		device->getMatrix6x6FieldCode()->importMatrix6x6Field(matrix6x6FieldBuffer, C_host);
		delete[] C_host;


	}
protected:
	const hardware::buffers::matrix6x6 * matrix6x6FieldBuffer;
};

template<typename TesterClass>
void performTest(const ReferenceValues refValuesIn, const LatticeExtents latticeExtentsIn, const Matrix6x6FieldFillType fillType)
{
	Matrix6x6FieldTestParameters parametersForThisTest {latticeExtentsIn, fillType};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	TesterClass tester(parameterCollection, parametersForThisTest, refValuesIn);
}

void testCloverEvenOddInverseExplizitUpperLeft(const ReferenceValues rV, const LatticeExtents lE, const Matrix6x6FieldFillType fT)
{
	performTest<CloverEvenOddInverseExplizitUpperLeftTester>(rV, lE, fT);
}

BOOST_AUTO_TEST_SUITE ( CLOVER_EO_INVERSE_EXPLIZIT )

	BOOST_AUTO_TEST_CASE( CLOVER_EO_INVERSE_EXPLIZIT_UPPER_LEFT )
	{
		testCloverEvenOddInverseExplizitUpperLeft(ReferenceValues{1234.}, LatticeExtents{ns4, nt4}, Matrix6x6FieldFillType::cold );
	}

BOOST_AUTO_TEST_SUITE_END()
