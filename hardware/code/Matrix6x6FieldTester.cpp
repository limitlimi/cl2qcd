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

hmc_float count_matrix6x6Field(Matrix6x6 * in, int size)
{
	hmc_float sum = 0.;
	for(int i = 0; i < size; i++)
	{
		sum += in[i].e00.re + in[i].e00.im +
				in[i].e01.re + in[i].e01.im +
				in[i].e02.re + in[i].e02.im +
				in[i].e03.re + in[i].e03.im +
				in[i].e04.re + in[i].e04.im +
				in[i].e05.re + in[i].e05.im +
				in[i].e10.re + in[i].e10.im +
				in[i].e11.re + in[i].e11.im +
				in[i].e12.re + in[i].e12.im +
				in[i].e13.re + in[i].e13.im +
				in[i].e14.re + in[i].e14.im +
				in[i].e15.re + in[i].e15.im +
				in[i].e20.re + in[i].e20.im +
				in[i].e21.re + in[i].e21.im +
				in[i].e22.re + in[i].e22.im +
				in[i].e23.re + in[i].e23.im +
				in[i].e24.re + in[i].e24.im +
				in[i].e25.re + in[i].e25.im +
				in[i].e30.re + in[i].e30.im +
				in[i].e31.re + in[i].e31.im +
				in[i].e32.re + in[i].e32.im +
				in[i].e33.re + in[i].e33.im +
				in[i].e34.re + in[i].e34.im +
				in[i].e35.re + in[i].e35.im +
				in[i].e40.re + in[i].e40.im +
				in[i].e41.re + in[i].e41.im +
				in[i].e42.re + in[i].e42.im +
				in[i].e43.re + in[i].e43.im +
				in[i].e44.re + in[i].e44.im +
				in[i].e45.re + in[i].e45.im +
				in[i].e50.re + in[i].e50.im +
				in[i].e51.re + in[i].e51.im +
				in[i].e52.re + in[i].e52.im +
				in[i].e53.re + in[i].e53.im +
				in[i].e54.re + in[i].e54.im +
				in[i].e55.re + in[i].e55.im;
	}
	return sum;
}

Matrix6x6 unit_matrix6x6()
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
}

Matrix6x6 zero_matrix6x6()
{
	Matrix6x6 out;
    out.e00.re = 0.; out.e00.im = 0.;
    out.e01.re = 0.; out.e01.im = 0.;
    out.e02.re = 0.; out.e02.im = 0.;
    out.e03.re = 0.; out.e03.im = 0.;
    out.e04.re = 0.; out.e04.im = 0.;
    out.e05.re = 0.; out.e05.im = 0.;
    //1. row
    out.e10.re = 0.; out.e10.im = 0.;
    out.e11.re = 0.; out.e11.im = 0.;
    out.e12.re = 0.; out.e12.im = 0.;
    out.e13.re = 0.; out.e13.im = 0.;
    out.e14.re = 0.; out.e14.im = 0.;
    out.e15.re = 0.; out.e15.im = 0.;
    //2. row
    out.e20.re = 0.; out.e20.im = 0.;
    out.e21.re = 0.; out.e21.im = 0.;
    out.e22.re = 0.; out.e22.im = 0.;
    out.e23.re = 0.; out.e23.im = 0.;
    out.e24.re = 0.; out.e24.im = 0.;
    out.e25.re = 0.; out.e25.im = 0.;
    //3. row
    out.e30.re = 0.; out.e30.im = 0.;
    out.e31.re = 0.; out.e31.im = 0.;
    out.e32.re = 0.; out.e32.im = 0.;
    out.e33.re = 0.; out.e33.im = 0.;
    out.e34.re = 0.; out.e34.im = 0.;
    out.e35.re = 0.; out.e35.im = 0.;
    //4. row
    out.e40.re = 0.; out.e40.im = 0.;
    out.e41.re = 0.; out.e41.im = 0.;
    out.e42.re = 0.; out.e42.im = 0.;
    out.e43.re = 0.; out.e43.im = 0.;
    out.e44.re = 0.; out.e44.im = 0.;
    out.e45.re = 0.; out.e45.im = 0.;
    //5. row
    out.e50.re = 0.; out.e50.im = 0.;
    out.e51.re = 0.; out.e51.im = 0.;
    out.e52.re = 0.; out.e52.im = 0.;
    out.e53.re = 0.; out.e53.im = 0.;
    out.e54.re = 0.; out.e54.im = 0.;
    out.e55.re = 0.; out.e55.im = 0.;

    return out;
}


Matrix6x6 ascendingRealMatrix6x6()
{
	Matrix6x6 out;
    out.e00.re = 1.; out.e00.im = 0.;
    out.e01.re = 2.; out.e01.im = 0.;
    out.e02.re = 3.; out.e02.im = 0.;
    out.e03.re = 4.; out.e03.im = 0.;
    out.e04.re = 5.; out.e04.im = 0.;
    out.e05.re = 6.; out.e05.im = 0.;
    //1. row
    out.e10.re = 7.; out.e10.im = 0.;
    out.e11.re = 8.; out.e11.im = 0.;
    out.e12.re = 9.; out.e12.im = 0.;
    out.e13.re = 10.; out.e13.im = 0.;
    out.e14.re = 11.; out.e14.im = 0.;
    out.e15.re = 12.; out.e15.im = 0.;
    //2. row
    out.e20.re = 13.; out.e20.im = 0.;
    out.e21.re = 14.; out.e21.im = 0.;
    out.e22.re = 15.; out.e22.im = 0.;
    out.e23.re = 16.; out.e23.im = 0.;
    out.e24.re = 17.; out.e24.im = 0.;
    out.e25.re = 18.; out.e25.im = 0.;
    //3. row
    out.e30.re = 19.; out.e30.im = 0.;
    out.e31.re = 20.; out.e31.im = 0.;
    out.e32.re = 21.; out.e32.im = 0.;
    out.e33.re = 22.; out.e33.im = 0.;
    out.e34.re = 23.; out.e34.im = 0.;
    out.e35.re = 24.; out.e35.im = 0.;
    //4. row
    out.e40.re = 25.; out.e40.im = 0.;
    out.e41.re = 26.; out.e41.im = 0.;
    out.e42.re = 27.; out.e42.im = 0.;
    out.e43.re = 28.; out.e43.im = 0.;
    out.e44.re = 29.; out.e44.im = 0.;
    out.e45.re = 30.; out.e45.im = 0.;
    //5. row
    out.e50.re = 31.; out.e50.im = 0.;
    out.e51.re = 32.; out.e51.im = 0.;
    out.e52.re = 33.; out.e52.im = 0.;
    out.e53.re = 34.; out.e53.im = 0.;
    out.e54.re = 35.; out.e54.im = 0.;
    out.e55.re = 36.; out.e55.im = 0.;

    return out;
}

int calculateMatrix6x6FieldSize(const LatticeExtents latticeExtentsIn) noexcept
{
	return 	latticeExtentsIn.getLatticeVolume();
}

const Matrix6x6* Matrix6x6FieldCreator::createMatrix6x6Field(const Matrix6x6FieldFillType fillTypeIn)
{
	Matrix6x6 * tmp = new Matrix6x6[numberOfElements];
	for(size_t i = 0; i < (size_t)numberOfElements; ++i)
	{
		switch (fillTypeIn)
		{
		case Matrix6x6FieldFillType::unity:
			tmp[i] = unit_matrix6x6();
			break;
		case Matrix6x6FieldFillType::ascendingReal6x6:
			tmp[i] = ascendingRealMatrix6x6();
			break;
		case Matrix6x6FieldFillType::zero6x6:
			tmp[i] = zero_matrix6x6();
			break;
		default:
			BOOST_ERROR("No valid Matrix6x6FieldFillType specified");
		}
	}
	return tmp;
}
