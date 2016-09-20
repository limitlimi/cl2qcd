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

#pragma once

#include "kernelTester.hpp"
#include "matrix6x6Field.hpp"
#include "GaugefieldTester.hpp"

int calculateMatrix6x6FieldSize(LatticeExtents latticeExtentsIn) noexcept;

struct Matrix6x6FieldTestParameters : public virtual TestParameters
{
	GaugefieldFillType fillType;

	Matrix6x6FieldTestParameters(const LatticeExtents latticeExtentsIn, GaugefieldFillType fillTypeIn):
		TestParameters(latticeExtentsIn), fillType( fillTypeIn ) {}
};

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

struct Matrix6x6FieldTester : public KernelTester
{
	Matrix6x6FieldTester(std::string kernelName, const ParameterCollection parameterCollection, const Matrix6x6FieldTestParameters testParameters, const ReferenceValues rV):
		KernelTester(kernelName, parameterCollection.hardwareParameters, parameterCollection.kernelParameters, testParameters, rV)
	{
		elements = testParameters.latticeExtents.getLatticeVolume();
		GaugefieldCreator gf(testParameters.latticeExtents);
		out = new hardware::buffers::matrix6x6(testParameters.latticeExtents, this->device);
		gaugefieldBuffer = new hardware::buffers::SU3(testParameters.latticeExtents, this->device);
		const Matrixsu3 * gf_host = gf.createGaugefield(testParameters.fillType);
		device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;

		code = this->device->getMatrix6x6FieldCode();
	}

protected:
	const hardware::buffers::matrix6x6 * out;
	const hardware::code::matrix6x6Field * code;
	const hardware::buffers::SU3 * gaugefieldBuffer;
	size_t elements;
};

struct Matrix6x6FieldTesterWithSumAsKernelResult : public Matrix6x6FieldTester
{
	Matrix6x6FieldTesterWithSumAsKernelResult(const std::string kernelName, const ParameterCollection parameterCollection, const Matrix6x6FieldTestParameters testParameters, const ReferenceValues rV) :
		Matrix6x6FieldTester(kernelName, parameterCollection, testParameters, rV) {}
	~Matrix6x6FieldTesterWithSumAsKernelResult()
	{
		Matrix6x6 * matrix6x6_in;
		matrix6x6_in = new Matrix6x6[elements];
		Matrix6x6FieldTester::out->dump(matrix6x6_in);
		Matrix6x6FieldTester::kernelResult.at(0) = count_matrix6x6Field(matrix6x6_in, elements);
		delete matrix6x6_in;
	}
};

//struct Matrix6x6FieldCreator
//{
//	Matrix6x6FieldCreator(const LatticeExtents lE): numberOfElements(calculateMatrix6x6FieldSize(lE)){};
//	const Matrix6x6* createMatrix6x6Field(const GaugefieldFillType fillTypeIn);
//
//	size_t numberOfElements;
//};
