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

int calculateMatrix6x6FieldSize(LatticeExtents latticeExtentsIn) noexcept;

enum Matrix6x6FieldFillType {cold = 1, nonTrivial};

struct Matrix6x6FieldTestParameters : public virtual TestParameters
{
	Matrix6x6FieldFillType fillType;

	Matrix6x6FieldTestParameters(const LatticeExtents latticeExtentsIn, Matrix6x6FieldFillType fillTypeIn):
		TestParameters(latticeExtentsIn), fillType( fillTypeIn ) {}
};

struct Matrix6x6FieldTester : public KernelTester
{
	Matrix6x6FieldTester(std::string, const ParameterCollection &, const Matrix6x6FieldTestParameters, const ReferenceValues rV);

protected:
	const hardware::code::matrix6x6Field * code;
};

struct Matrix6x6FieldCreator
{
	Matrix6x6FieldCreator(const LatticeExtents lE): numberOfElements(calculateMatrix6x6FieldSize(lE)){};
	const Matrix6x6* createMatrix6x6Field(const Matrix6x6FieldFillType fillTypeIn);

	size_t numberOfElements;
};
