
void put_3x3block_matrix6x6_upperleft(Matrix6x6 inout, Matrix3x3 p)
{
    inout.e00.re = p.e00.re;
    inout.e00.im = p.e00.im;
    inout.e01.re = p.e01.re ;
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
    inout.e04.re = p.e01.re ;
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
    inout.e31.re = p.e01.re ;
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
    inout.e34.re = p.e01.re ;
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