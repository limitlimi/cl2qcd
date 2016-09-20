/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra, 
 * Max Theilig
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


hmc_complex complex_add(const hmc_complex in1, const hmc_complex in2)
{
    hmc_complex out;
    out.re = in1.re + in2.re;
    out.im = in1.im + in2.im;
    return out;
}

hmc_complex complex_sub(const hmc_complex in1, const hmc_complex in2)
{
    hmc_complex out;
    out.re = in1.re - in2.re;
    out.im = in1.im - in2.im;
    return out;
}

hmc_complex complex_mult(const hmc_complex in1, const hmc_complex in2)
{
    hmc_complex out;
    out.re = in1.re * in2.re - in1.im * in2.im;
    out.im = in1.im * in2.re + in1.re * in2.im;
    return out;
}

hmc_complex complex_divid(const hmc_complex numerator, const hmc_complex denominator)
{
    hmc_float norm = denominator.re * denominator.re + denominator.im * denominator.im;
    hmc_complex out;
    out.re = (numerator.re * denominator.re + numerator.im * denominator.im ) / norm;
    out.im = (numerator.im * denominator.re - numerator.re * denominator.im ) / norm;
    return out;
}

hmc_complex complex_conj(hmc_complex in)
{
    in.im = -(in.im);
    return in;
}

hmc_float complex_abs_value(hmc_complex in)
{
    hmc_float out, tmp;
    if(in.re == 0. && in.im == 0.){
	return 0.;}
    else{
    tmp = in.re * in.re + in.im * in.im;
    out = sqrt(tmp);
    return out;}
}

hmc_complex convert_float_to_complex(hmc_float in)
{
    hmc_complex out;
    out.re = in;
    out.im = 0.;
    return out;
}
