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

/** @file
 * Device code implementing 3x3 matrices
 */

//operations_matrix6x6.cl

#ifdef ENABLE_PRINTF
void print_matrix6x6(Matrix6x6 in)
{
	printf("(%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f)\n (%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f) (%f,%f)\n",
	       in.e00.re, in.e00.im, in.e01.re, in.e01.im, in.e02.re, in.e02.im, in.e03.re, in.e03.im, in.e04.re, in.e04.im, in.e05.re, in.e05.im,
	       in.e10.re, in.e10.im, in.e11.re, in.e11.im, in.e12.re, in.e12.im, in.e13.re, in.e13.im, in.e14.re, in.e14.im, in.e15.re, in.e15.im,
	       in.e20.re, in.e20.im, in.e21.re, in.e21.im, in.e22.re, in.e22.im, in.e23.re, in.e23.im, in.e24.re, in.e24.im, in.e25.re, in.e25.im,
	       in.e30.re, in.e30.im, in.e31.re, in.e31.im, in.e32.re, in.e32.im, in.e33.re, in.e33.im, in.e34.re, in.e34.im, in.e35.re, in.e35.im,
	       in.e40.re, in.e40.im, in.e41.re, in.e41.im, in.e42.re, in.e42.im, in.e43.re, in.e43.im, in.e44.re, in.e44.im, in.e45.re, in.e45.im,
	       in.e50.re, in.e50.im, in.e51.re, in.e51.im, in.e52.re, in.e52.im, in.e53.re, in.e53.im, in.e54.re, in.e54.im, in.e55.re, in.e55.im);
	printf("\n");
}
#endif

inline Matrix6x6 get6x6(__global const Matrix6x6StorageType * const restrict in, const uint idx)
{
    return in[idx];
}

inline void put6x6(__global Matrix6x6StorageType * const restrict out, const uint idx, const Matrix6x6 val)
{
    out[idx] = val;
}

Matrix3x3 get_3x3_block_upperleft(Matrix6x6 in)
{
    Matrix3x3 out;
    
    out.e00.re = in.e00.re;
    out.e00.im = in.e00.im;
    out.e01.re = in.e01.re;
    out.e01.im = in.e01.im;
    out.e02.re = in.e02.re;
    out.e02.im = in.e02.im;
    
    out.e10.re = in.e10.re;
    out.e10.im = in.e10.im;
    out.e11.re = in.e11.re;
    out.e11.im = in.e11.im;
    out.e12.re = in.e12.re;
    out.e12.im = in.e12.im;
    
    out.e20.re = in.e20.re;
    out.e20.im = in.e20.im;
    out.e21.re = in.e21.re;
    out.e21.im = in.e21.im;
    out.e22.re = in.e22.re;
    out.e22.im = in.e22.im;
    
    return out;
}


Matrix3x3 get_3x3_block_upperright(Matrix6x6 in)
{
    Matrix3x3 out;
    
    out.e00.re = in.e03.re;
    out.e00.im = in.e03.im;
    out.e01.re = in.e04.re;
    out.e01.im = in.e04.im;
    out.e02.re = in.e05.re;
    out.e02.im = in.e05.im;
    
    out.e10.re = in.e13.re;
    out.e10.im = in.e13.im;
    out.e11.re = in.e14.re;
    out.e11.im = in.e14.im;
    out.e12.re = in.e15.re;
    out.e12.im = in.e15.im;
    
    out.e20.re = in.e23.re;
    out.e20.im = in.e23.im;
    out.e21.re = in.e24.re;
    out.e21.im = in.e24.im;
    out.e22.re = in.e25.re;
    out.e22.im = in.e25.im;
    
    return out;
}

Matrix3x3 get_3x3_block_lowerleft(Matrix6x6 in)
{
    Matrix3x3 out;
    
    out.e00.re = in.e30.re;
    out.e00.im = in.e30.im;
    out.e01.re = in.e31.re;
    out.e01.im = in.e31.im;
    out.e02.re = in.e32.re;
    out.e02.im = in.e32.im;
    
    out.e10.re = in.e40.re;
    out.e10.im = in.e40.im;
    out.e11.re = in.e41.re;
    out.e11.im = in.e41.im;
    out.e12.re = in.e42.re;
    out.e12.im = in.e42.im;
    
    out.e20.re = in.e50.re;
    out.e20.im = in.e50.im;
    out.e21.re = in.e51.re;
    out.e21.im = in.e51.im;
    out.e22.re = in.e52.re;
    out.e22.im = in.e52.im;
    
    return out;
}

Matrix3x3 get_3x3_block_lowerright(Matrix6x6 in)
{
    Matrix3x3 out;
    
    out.e00.re = in.e33.re;
    out.e00.im = in.e33.im;
    out.e01.re = in.e34.re;
    out.e01.im = in.e34.im;
    out.e02.re = in.e35.re;
    out.e02.im = in.e35.im;
    
    out.e10.re = in.e43.re;
    out.e10.im = in.e43.im;
    out.e11.re = in.e44.re;
    out.e11.im = in.e44.im;
    out.e12.re = in.e45.re;
    out.e12.im = in.e45.im;
    
    out.e20.re = in.e53.re;
    out.e20.im = in.e53.im;
    out.e21.re = in.e54.re;
    out.e21.im = in.e54.im;
    out.e22.re = in.e55.re;
    out.e22.im = in.e55.im;
    
    return out;
}


//multiply matrix6x6 times halfspinor using 3x3 blocks and su3vecs, u={{u_11,u_12},{u_21,u_22}}, in={psi_1,psi_2}
//u * in = {u_11 * psi_1 + u_12 * psi_2, u_21 * psi_1 + u_22 * psi_2}
halfspinor matrix6x6_times_halfspinor(Matrix6x6 u, halfspinor in)
{
    halfspinor out;
    su3vec tmp1, tmp2, tmp3, tmp4;
    
    tmp1 = matrix3x3_times_su3vec(get_3x3_block_upperleft(u), in.e0);
    tmp2 = matrix3x3_times_su3vec(get_3x3_block_upperright(u), in.e1);
    tmp3 = matrix3x3_times_su3vec(get_3x3_block_lowerleft(u), in.e0);
    tmp4 = matrix3x3_times_su3vec(get_3x3_block_lowerright(u), in.e1);
    
    out.e0 = su3vec_acc(tmp1, tmp2);
    out.e1 = su3vec_acc(tmp3, tmp4);
    return out;
}


Matrix6x6 put_3x3block_matrix6x6_upperleft(Matrix6x6 in, Matrix3x3 p)
{
	Matrix6x6 out = in;
	
    out.e00.re = p.e00.re;
    out.e00.im = p.e00.im;
    out.e01.re = p.e01.re;
    out.e01.im = p.e01.im;
    out.e02.re = p.e02.re;
    out.e02.im = p.e02.im;
    
    out.e10.re = p.e10.re;
    out.e10.im = p.e10.im;
    out.e11.re = p.e11.re;
    out.e11.im = p.e11.im;
    out.e12.re = p.e12.re;
    out.e12.im = p.e12.im;
    
    out.e20.re = p.e20.re;
    out.e20.im = p.e20.im;
    out.e21.re = p.e21.re;
    out.e21.im = p.e21.im;
    out.e22.re = p.e22.re;
    out.e22.im = p.e22.im;
    
    return out;
}

Matrix6x6 put_3x3block_matrix6x6_upperright(Matrix6x6 in, Matrix3x3 p)
{
	Matrix6x6 out = in;

    out.e03.re = p.e00.re;
    out.e03.im = p.e00.im;
    out.e04.re = p.e01.re;
    out.e04.im = p.e01.im;
    out.e05.re = p.e02.re;
    out.e05.im = p.e02.im;
    
    out.e13.re = p.e10.re;
    out.e13.im = p.e10.im;
    out.e14.re = p.e11.re;
    out.e14.im = p.e11.im;
    out.e15.re = p.e12.re;
    out.e15.im = p.e12.im;
    
    out.e23.re = p.e20.re;
    out.e23.im = p.e20.im;
    out.e24.re = p.e21.re;
    out.e24.im = p.e21.im;
    out.e25.re = p.e22.re;
    out.e25.im = p.e22.im;
    
    return out;
}

Matrix6x6 put_3x3block_matrix6x6_lowerleft(Matrix6x6 in, Matrix3x3 p)
{
	Matrix6x6 out = in; 
	
    out.e30.re = p.e00.re;
    out.e30.im = p.e00.im;
    out.e31.re = p.e01.re;
    out.e31.im = p.e01.im;
    out.e32.re = p.e02.re;
    out.e32.im = p.e02.im;
    
    out.e40.re = p.e10.re;
    out.e40.im = p.e10.im;
    out.e41.re = p.e11.re;
    out.e41.im = p.e11.im;
    out.e42.re = p.e12.re;
    out.e42.im = p.e12.im;
    
    out.e50.re = p.e20.re;
    out.e50.im = p.e20.im;
    out.e51.re = p.e21.re;
    out.e51.im = p.e21.im;
    out.e52.re = p.e22.re;
    out.e52.im = p.e22.im;
    
    return out;
}

Matrix6x6 put_3x3block_matrix6x6_lowerright(Matrix6x6 in, Matrix3x3 p)
{
	Matrix6x6 out = in;
	
    out.e33.re = p.e00.re;
    out.e33.im = p.e00.im;
    out.e34.re = p.e01.re;
    out.e34.im = p.e01.im;
    out.e35.re = p.e02.re;
    out.e35.im = p.e02.im;
    
    out.e43.re = p.e10.re;
    out.e43.im = p.e10.im;
    out.e44.re = p.e11.re;
    out.e44.im = p.e11.im;
    out.e45.re = p.e12.re;
    out.e45.im = p.e12.im;
    
    out.e53.re = p.e20.re;
    out.e53.im = p.e20.im;
    out.e54.re = p.e21.re;
    out.e54.im = p.e21.im;
    out.e55.re = p.e22.re;
    out.e55.im = p.e22.im;
    
    return out;
}

inline Matrix6x6 zero_matrix6x6 ()
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

inline Matrix6x6 identity_matrix6x6 ()
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

inline Matrix6x6 multiply_matrix6x6_by_complex (Matrix6x6 in, hmc_complex factor)
{
    Matrix6x6 out;
    
    out.e00 = complexmult(in.e00, factor);
    out.e01 = complexmult(in.e01, factor);
    out.e02 = complexmult(in.e02, factor);
    out.e03 = complexmult(in.e03, factor);
    out.e04 = complexmult(in.e04, factor);
    out.e05 = complexmult(in.e05, factor);
    out.e10 = complexmult(in.e10, factor);
    out.e11 = complexmult(in.e11, factor);
    out.e12 = complexmult(in.e12, factor);
    out.e13 = complexmult(in.e13, factor);
    out.e14 = complexmult(in.e14, factor);
    out.e15 = complexmult(in.e15, factor);
    out.e20 = complexmult(in.e20, factor);
    out.e21 = complexmult(in.e21, factor);
    out.e22 = complexmult(in.e22, factor);
    out.e23 = complexmult(in.e23, factor);
    out.e24 = complexmult(in.e24, factor);
    out.e25 = complexmult(in.e25, factor);
    out.e30 = complexmult(in.e30, factor);
    out.e31 = complexmult(in.e31, factor);
    out.e32 = complexmult(in.e32, factor);
    out.e33 = complexmult(in.e33, factor);
    out.e34 = complexmult(in.e34, factor);
    out.e35 = complexmult(in.e35, factor);
    out.e40 = complexmult(in.e40, factor);
    out.e41 = complexmult(in.e41, factor);
    out.e42 = complexmult(in.e42, factor);
    out.e43 = complexmult(in.e43, factor);
    out.e44 = complexmult(in.e44, factor);
    out.e45 = complexmult(in.e45, factor);
    out.e50 = complexmult(in.e50, factor);
    out.e51 = complexmult(in.e51, factor);
    out.e52 = complexmult(in.e52, factor);
    out.e53 = complexmult(in.e53, factor);
    out.e54 = complexmult(in.e54, factor);
    out.e55 = complexmult(in.e55, factor);
    
    return out;
}
