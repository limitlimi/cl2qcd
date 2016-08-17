
//multiply matrix6x6 times halfspinor using 3x3 blocks and su3vecs, u={{u_11,u_12},{u_21,u_22}}, in={psi_1,psi_2}
//u * in = {u_11 * psi_1 + u_12 * psi_2, u_21 * psi_1 + u_22 * psi_2}
halfspinor matrix6x6_times_halfspinor(Matrix6x6 u, halfspinor in)
{
    halfspinor out;
    su3vec tmp1, tmp2, tmp3, tmp4;
    
    tmp1 = matrix3x3_times_su3vec(get_3x3_block_upperleft(u), in.e1);
    tmp2 = matrix3x3_times_su3vec(get_3x3_block_upperright(u), in.e2);
    tmp3 = matrix3x3_times_su3vec(get_3x3_block_lowerleft(u), in.e1);
    tmp4 = matrix3x3_times_su3vec(get_3x3_block_lowerright(u), in.e2);
    
    out.e1 = su3vec_acc(tmp1, tmp2);
    out.e2 = su3vec_acc(tmp1, tmp2);
    return out;
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


void put_3x3block_matrix6x6_upperleft(Matrix6x6 inout, Matrix3x3 p)
{
    inout.e00.re = p.e00.re;
    inout.e00.im = p.e00.im;
    inout.e01.re = p.e01.re;
    inout.e01.im = p.e01.im;
    inout.e02.re = p.e02.re;
    inout.e02.im = p.e02.im;
    
    inout.e10.re = p.e10.re;
    inout.e10.im = p.e10.im;
    inout.e11.re = p.e11.re;
    inout.e11.im = p.e11.im;
    inout.e12.re = p.e12.re;
    inout.e12.im = p.e12.im;
    
    inout.e20.re = p.e20.re;
    inout.e20.im = p.e20.im;
    inout.e21.re = p.e21.re;
    inout.e21.im = p.e21.im;
    inout.e22.re = p.e22.re;
    inout.e22.im = p.e22.im;
}

void put_3x3block_matrix6x6_upperright(Matrix6x6 inout, Matrix3x3 p)
{
    inout.e03.re = p.e00.re;
    inout.e03.im = p.e00.im;
    inout.e04.re = p.e01.re;
    inout.e04.im = p.e01.im;
    inout.e05.re = p.e02.re;
    inout.e05.im = p.e02.im;
    
    inout.e13.re = p.e10.re;
    inout.e13.im = p.e10.im;
    inout.e14.re = p.e11.re;
    inout.e14.im = p.e11.im;
    inout.e15.re = p.e12.re;
    inout.e15.im = p.e12.im;
    
    inout.e23.re = p.e20.re;
    inout.e23.im = p.e20.im;
    inout.e24.re = p.e21.re;
    inout.e24.im = p.e21.im;
    inout.e25.re = p.e22.re;
    inout.e25.im = p.e22.im;
}

void put_3x3block_matrix6x6_lowerleft(Matrix6x6 inout, Matrix3x3 p)
{
    inout.e30.re = p.e00.re;
    inout.e30.im = p.e00.im;
    inout.e31.re = p.e01.re;
    inout.e31.im = p.e01.im;
    inout.e32.re = p.e02.re;
    inout.e32.im = p.e02.im;
    
    inout.e40.re = p.e10.re;
    inout.e40.im = p.e10.im;
    inout.e41.re = p.e11.re;
    inout.e41.im = p.e11.im;
    inout.e42.re = p.e12.re;
    inout.e42.im = p.e12.im;
    
    inout.e50.re = p.e20.re;
    inout.e50.im = p.e20.im;
    inout.e51.re = p.e21.re;
    inout.e51.im = p.e21.im;
    inout.e52.re = p.e22.re;
    inout.e52.im = p.e22.im;
}

void put_3x3block_matrix6x6_lowerright(Matrix6x6 inout, Matrix3x3 p)
{
    inout.e33.re = p.e00.re;
    inout.e33.im = p.e00.im;
    inout.e34.re = p.e01.re;
    inout.e34.im = p.e01.im;
    inout.e35.re = p.e02.re;
    inout.e35.im = p.e02.im;
    
    inout.e43.re = p.e10.re;
    inout.e43.im = p.e10.im;
    inout.e44.re = p.e11.re;
    inout.e44.im = p.e11.im;
    inout.e45.re = p.e12.re;
    inout.e45.im = p.e12.im;
    
    inout.e53.re = p.e20.re;
    inout.e53.im = p.e20.im;
    inout.e54.re = p.e21.re;
    inout.e54.im = p.e21.im;
    inout.e55.re = p.e22.re;
    inout.e55.im = p.e22.im;
}

inline Matrix6x6 zero_matrix6x6 ()
{
    return (Matrix3x3) {
        {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.},
        {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.},
        {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.},
        {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.},
        {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.},
        {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.},
    };
}

inline Matrix6x6 identity_matrix6x6 ()
{
    return (Matrix3x3) {
        {1., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.},
        {0., 0.}, {1., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.},
        {0., 0.}, {0., 0.}, {1., 0.}, {0., 0.}, {0., 0.}, {0., 0.},
        {0., 0.}, {0., 0.}, {0., 0.}, {1., 0.}, {0., 0.}, {0., 0.},
        {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {1., 0.}, {0., 0.},
        {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {1., 0.},
    };
}

inline Matrix6x6 multiply_matrix6x6_by_complex (Matrix3x3 in, hmc_complex factor)
{
    Matrix3x3 out;
    
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