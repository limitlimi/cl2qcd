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

enum Matrix6x6FieldFillType {unity = 1, toBeAdded};

struct Matrix6x6FieldTestParameters : public virtual TestParameters
{
	const GaugefieldFillType fillType;
	const hmc_float kappa;
	const hmc_float csw;

	Matrix6x6FieldTestParameters(const LatticeExtents latticeExtentsIn, const GaugefieldFillType fillTypeIn, const hmc_float kappaIn, const hmc_float cswIn):
		TestParameters(latticeExtentsIn), fillType( fillTypeIn ), kappa(kappaIn), csw(cswIn){}
};

hmc_float count_matrix6x6Field(Matrix6x6 * in, int size);

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

struct Matrix6x6FieldCreator
{
	Matrix6x6FieldCreator(const LatticeExtents lE): numberOfElements(calculateMatrix6x6FieldSize(lE)){};
	const Matrix6x6* createMatrix6x6Field(const Matrix6x6FieldFillType fillTypeIn);

	size_t numberOfElements;
};
