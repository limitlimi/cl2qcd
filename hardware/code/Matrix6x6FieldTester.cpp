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

#include "Matrix6x6FieldTester.hpp"

//Matrix6x6 unit_matrix6x6()
//{
//	Matrix6x6 out;
//    out.e00.re = 1.; out.e00.im = 0.;
//    out.e01.re = 0.; out.e01.im = 0.;
//    out.e02.re = 0.; out.e02.im = 0.;
//    out.e03.re = 0.; out.e03.im = 0.;
//    out.e04.re = 0.; out.e04.im = 0.;
//    out.e05.re = 0.; out.e05.im = 0.;
//    //1. row
//    out.e10.re = 0.; out.e10.im = 0.;
//    out.e11.re = 1.; out.e11.im = 0.;
//    out.e12.re = 0.; out.e12.im = 0.;
//    out.e13.re = 0.; out.e13.im = 0.;
//    out.e14.re = 0.; out.e14.im = 0.;
//    out.e15.re = 0.; out.e15.im = 0.;
//    //2. row
//    out.e20.re = 0.; out.e20.im = 0.;
//    out.e21.re = 0.; out.e21.im = 0.;
//    out.e22.re = 1.; out.e22.im = 0.;
//    out.e23.re = 0.; out.e23.im = 0.;
//    out.e24.re = 0.; out.e24.im = 0.;
//    out.e25.re = 0.; out.e25.im = 0.;
//    //3. row
//    out.e30.re = 0.; out.e30.im = 0.;
//    out.e31.re = 0.; out.e31.im = 0.;
//    out.e32.re = 0.; out.e32.im = 0.;
//    out.e33.re = 1.; out.e33.im = 0.;
//    out.e34.re = 0.; out.e34.im = 0.;
//    out.e35.re = 0.; out.e35.im = 0.;
//    //4. row
//    out.e40.re = 0.; out.e40.im = 0.;
//    out.e41.re = 0.; out.e41.im = 0.;
//    out.e42.re = 0.; out.e42.im = 0.;
//    out.e43.re = 0.; out.e43.im = 0.;
//    out.e44.re = 1.; out.e44.im = 0.;
//    out.e45.re = 0.; out.e45.im = 0.;
//    //5. row
//    out.e50.re = 0.; out.e50.im = 0.;
//    out.e51.re = 0.; out.e51.im = 0.;
//    out.e52.re = 0.; out.e52.im = 0.;
//    out.e53.re = 0.; out.e53.im = 0.;
//    out.e54.re = 0.; out.e54.im = 0.;
//    out.e55.re = 1.; out.e55.im = 0.;
//
//    return out;
//}

/*Matrix6x6 nonTrivialMatrix6x6()
{
	Matrix6x6 out;
    out.e00.re = 1.; out.e00.im = 0.;
    out.e01.re = 0.; out.e01.im = 0.;
    out.e02.re = 0.; out.e02.im = 0.;
    out.e03.re = 0.; out.e03.im = 0.;
    out.e04.re = 0.; out.e04.im = 0.;
    out.e05.re = 0.; out.e05.im = 0.;
    //1. row
    out.e10.re = 0.; out.e10.im = 0.;
    out.e11.re = 1.; out.e11.im = 0.;
    out.e12.re = 0.; out.e12.im = 0.;
    out.e13.re = 0.; out.e13.im = 0.;
    out.e14.re = 0.; out.e14.im = 0.;
    out.e15.re = 0.; out.e15.im = 0.;
    //2. row
    out.e20.re = 0.; out.e20.im = 0.;
    out.e21.re = 0.; out.e21.im = 0.;
    out.e22.re = 1.; out.e22.im = 0.;
    out.e23.re = 0.; out.e23.im = 0.;
    out.e24.re = 0.; out.e24.im = 0.;
    out.e25.re = 0.; out.e25.im = 0.;
    //3. row
    out.e30.re = 0.; out.e30.im = 0.;
    out.e31.re = 0.; out.e31.im = 0.;
    out.e32.re = 0.; out.e32.im = 0.;
    out.e33.re = 1.; out.e33.im = 0.;
    out.e34.re = 0.; out.e34.im = 0.;
    out.e35.re = 0.; out.e35.im = 0.;
    //4. row
    out.e40.re = 0.; out.e40.im = 0.;
    out.e41.re = 0.; out.e41.im = 0.;
    out.e42.re = 0.; out.e42.im = 0.;
    out.e43.re = 0.; out.e43.im = 0.;
    out.e44.re = 1.; out.e44.im = 0.;
    out.e45.re = 0.; out.e45.im = 0.;
    //5. row
    out.e50.re = 0.; out.e50.im = 0.;
    out.e51.re = 0.; out.e51.im = 0.;
    out.e52.re = 0.; out.e52.im = 0.;
    out.e53.re = 0.; out.e53.im = 0.;
    out.e54.re = 0.; out.e54.im = 0.;
    out.e55.re = 1.; out.e55.im = 0.;

    return out;
}*/


//Matrix6x6FieldTester::Matrix6x6FieldTester(std::string kernelName, const ParameterCollection & parameterCollection, const Matrix6x6FieldTestParameters testParams, const ReferenceValues rV):
//KernelTester(kernelName, parameterCollection.hardwareParameters, parameterCollection.kernelParameters, testParams, rV)
//{
//	code = device->getMatrix6x6FieldCode();
//}

int calculateMatrix6x6FieldSize(const LatticeExtents latticeExtentsIn) noexcept
{
	return 	latticeExtentsIn.getLatticeVolume();
}

//const Matrix6x6* Matrix6x6FieldCreator::createMatrix6x6Field(const Matrix6x6FieldFillType fillTypeIn)
//{
//	Matrix6x6 * tmp = new Matrix6x6[numberOfElements];
//	for(size_t i = 0; i < (size_t)numberOfElements; ++i)
//	{
//		switch (fillTypeIn)
//		{
//		case Matrix6x6FieldFillType::cold:
//			tmp[i] = unit_matrix6x6();
//			break;
//		/*case Matrix6x6FieldFillType::nonTrivial:
//			tmp[i] = nonTrivialMatrix6x6();
//			break;*/
//		default:
//			BOOST_ERROR("No valid Matrix6x6FieldFillType specified");
//		}
//	}
//	return tmp;
//}
