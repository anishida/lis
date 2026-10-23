#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include <stdio.h>
#include <math.h>

#include "lis.h"
#include "lis_esolver.h"


static int close_real(
        LIS_REAL a,
        LIS_REAL b,
        LIS_REAL tol)
{
        LIS_REAL scale;

        scale = (LIS_REAL)1.0+fabs(b);

        return fabs(a-b) <= tol*scale;
}


static int test_gauss_legendre(void)
{
        LIS_REAL x[8];
        LIS_REAL w[8];

        LIS_REAL m0;
        LIS_REAL m1;
        LIS_REAL m2;
        LIS_REAL m14;

        LIS_INT i;
        LIS_INT err;

        m0 = 0.0;
        m1 = 0.0;
        m2 = 0.0;
        m14 = 0.0;

        err = lis_efeast_gauss_legendre(
                8,
                x,
                w);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: Gauss-Legendre generation returned %d\n",
                        (int)err);

                return 1;
        }

        for(i=0;i<8;i++)
        {
                LIS_REAL x2;
                LIS_REAL x4;
                LIS_REAL x8;
                LIS_REAL x14;

                if( w[i]<=0.0 )
                {
                        fprintf(
                                stderr,
                                "FAIL: non-positive Gauss weight at %d\n",
                                (int)i);

                        return 1;
                }

                if( i>0 && !(x[i-1]<x[i]) )
                {
                        fprintf(
                                stderr,
                                "FAIL: Gauss nodes are not strictly ordered\n");

                        return 1;
                }

                x2 = x[i]*x[i];
                x4 = x2*x2;
                x8 = x4*x4;
                x14 = x8*x4*x2;

                m0 += w[i];
                m1 += w[i]*x[i];
                m2 += w[i]*x2;
                m14 += w[i]*x14;
        }

        if( !close_real(
                    m0,
                    (LIS_REAL)2.0,
                    (LIS_REAL)2.0e-13) )
        {
                fprintf(
                        stderr,
                        "FAIL: integral 1 = %.17g\n",
                        (double)m0);

                return 1;
        }

        if( !close_real(
                    m1,
                    (LIS_REAL)0.0,
                    (LIS_REAL)2.0e-13) )
        {
                fprintf(
                        stderr,
                        "FAIL: integral x = %.17g\n",
                        (double)m1);

                return 1;
        }

        if( !close_real(
                    m2,
                    (LIS_REAL)(2.0/3.0),
                    (LIS_REAL)2.0e-13) )
        {
                fprintf(
                        stderr,
                        "FAIL: integral x^2 = %.17g\n",
                        (double)m2);

                return 1;
        }

        if( !close_real(
                    m14,
                    (LIS_REAL)(2.0/15.0),
                    (LIS_REAL)3.0e-13) )
        {
                fprintf(
                        stderr,
                        "FAIL: integral x^14 = %.17g\n",
                        (double)m14);

                return 1;
        }

        return 0;
}


static int test_feast_reference_contour(void)
{
        LIS_REAL zr[4];
        LIS_REAL zi[4];
        LIS_REAL wr[4];
        LIS_REAL wi[4];

        const LIS_REAL ref_zr[4] = {
                (LIS_REAL)1.0473908518410222,
                (LIS_REAL)1.9819684290084905,
                (LIS_REAL)4.0180315709915098,
                (LIS_REAL)4.9526091481589782
        };

        const LIS_REAL ref_zi[4] = {
                (LIS_REAL)0.12984058035656038,
                (LIS_REAL)0.51645431050753199,
                (LIS_REAL)0.51645431050753199,
                (LIS_REAL)0.12984058035656015
        };

        const LIS_REAL ref_wr[4] = {
                (LIS_REAL)-0.050941841463511224,
                (LIS_REAL)-0.049792826738941448,
                (LIS_REAL) 0.049792826738941448,
                (LIS_REAL) 0.050941841463511231
        };

        const LIS_REAL ref_wi[4] = {
                (LIS_REAL)0.037638062477073672,
                (LIS_REAL)0.28066931358780345,
                (LIS_REAL)0.28066931358780339,
                (LIS_REAL)0.037638062477073610
        };

        LIS_INT i;
        LIS_INT err;

        err = lis_efeast_contour_gauss(
                (LIS_REAL)1.0,
                (LIS_REAL)5.0,
                4,
                (LIS_REAL)0.30,
                zr,
                zi,
                wr,
                wi);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: FEAST contour generation returned %d\n",
                        (int)err);

                return 1;
        }

        for(i=0;i<4;i++)
        {
                if( !close_real(
                            zr[i],
                            ref_zr[i],
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            zi[i],
                            ref_zi[i],
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            wr[i],
                            ref_wr[i],
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            wi[i],
                            ref_wi[i],
                            (LIS_REAL)5.0e-13) )
                {
                        fprintf(
                                stderr,
                                "FAIL: contour node %d differs from reference\n"
                                "  z=(%.17g, %.17g)\n"
                                "  w=(%.17g, %.17g)\n",
                                (int)i,
                                (double)zr[i],
                                (double)zi[i],
                                (double)wr[i],
                                (double)wi[i]);

                        return 1;
                }
        }

        return 0;
}


static int test_contour_symmetry(void)
{
        LIS_REAL zr[8];
        LIS_REAL zi[8];
        LIS_REAL wr[8];
        LIS_REAL wi[8];

        LIS_INT i;
        LIS_INT j;
        LIS_INT err;

        err = lis_efeast_contour_gauss(
                (LIS_REAL)-3.0,
                (LIS_REAL)7.0,
                8,
                (LIS_REAL)0.30,
                zr,
                zi,
                wr,
                wi);

        if( err )
        {
                return 1;
        }

        for(i=0;i<4;i++)
        {
                j = 7-i;

                if( !close_real(
                            zr[i]+zr[j],
                            (LIS_REAL)4.0,
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            zi[i],
                            zi[j],
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            wr[i],
                            -wr[j],
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            wi[i],
                            wi[j],
                            (LIS_REAL)5.0e-13) )
                {
                        fprintf(
                                stderr,
                                "FAIL: contour symmetry at pair %d/%d\n",
                                (int)i,
                                (int)j);

                        return 1;
                }
        }

        return 0;
}


int main(
        int argc,
        char **argv)
{
        LIS_ESOLVER esolver;
        LIS_INT err;

        esolver = NULL;

        err = lis_initialize(
                &argc,
                &argv);

        if( err )
        {
                return (int)err;
        }

        if( test_gauss_legendre() )
        {
                goto fail;
        }

        if( test_feast_reference_contour() )
        {
                goto fail;
        }

        if( test_contour_symmetry() )
        {
                goto fail;
        }

        err = lis_esolver_create(
                &esolver);

        if( err )
        {
                goto fail;
        }

        err = lis_esolver_set_option(
                "-e feast "
                "-ss 12 "
                "-feast_emin 1.0 "
                "-feast_emax 5.0 "
                "-feast_ratio 0.30 "
                "-feast_points 8",
                esolver);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: valid contour options rejected\n");

                goto fail;
        }

        err = lis_efeast_check_params(
                esolver);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: valid contour configuration rejected\n");

                goto fail;
        }

        lis_esolver_destroy(esolver);
        lis_finalize();

        printf(
                "PASS: native FEAST Gauss contour\n");

        return 0;

fail:

        if( esolver!=NULL )
        {
                lis_esolver_destroy(esolver);
        }

        lis_finalize();

        return 1;
}
