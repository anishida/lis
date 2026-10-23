#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <stdio.h>
#include <math.h>

#ifdef _COMPLEX
#include <complex.h>
#endif

#include "lis.h"
#include "lis_esolver.h"


#ifdef _COMPLEX
#define TEST_ABS(x) cabs(x)
#define TEST_REAL(x) creal(x)
#define TEST_IMAG(x) cimag(x)
#define TEST_CONJ(x) conj(x)
#else
#define TEST_ABS(x) fabs(x)
#define TEST_REAL(x) (x)
#define TEST_IMAG(x) 0.0
#define TEST_CONJ(x) (x)
#endif


static int
check_projected_system(void)
{
    LIS_SCALAR a[9];
    LIS_SCALAR original[9];

    LIS_SCALAR eval[3];
    LIS_SCALAR evec[9];

    LIS_SCALAR dot;
    LIS_SCALAR residual;

    LIS_REAL value[3];
    LIS_REAL expected[3];

    LIS_REAL norm;
    LIS_REAL rnorm;
    LIS_REAL trace;
    LIS_REAL tmp;

    LIS_INT niter;
    LIS_REAL error;

    LIS_INT err;
    LIS_INT i,j,k,l;


    expected[0] =
        3.0-sqrt(3.0);

    expected[1] =
        3.0;

    expected[2] =
        3.0+sqrt(3.0);


#ifdef _COMPLEX

    {
        LIS_REAL s;

        s = 1.0/sqrt(2.0);

        a[0] = 2.0;
        a[1] = s + s*I;
        a[2] = 0.0;

        a[3] = s - s*I;
        a[4] = 3.0;
        a[5] = 1.0*I;

        a[6] = 0.0;
        a[7] = -1.0*I;
        a[8] = 4.0;
    }

#else

    a[0] = 2.0;
    a[1] = 1.0;
    a[2] = 0.0;

    a[3] = 1.0;
    a[4] = 3.0;
    a[5] = 1.0;

    a[6] = 0.0;
    a[7] = 1.0;
    a[8] = 4.0;

#endif


    for(i=0;i<9;i++)
    {
        original[i] = a[i];
    }


    niter = 0;
    error = 0.0;

    err = lis_eigen_projected_solve(
        3,
        a,
        eval,
        evec,
        &niter,
        &error);

    if( err!=LIS_SUCCESS )
    {
        printf(
            "FAIL: projected solve returned %d\n",
            (int)err);

        return 1;
    }


    if( niter<=0 || niter>50 )
    {
        printf(
            "FAIL: invalid iteration count %d\n",
            (int)niter);

        return 1;
    }


    if( !isfinite((double)error) )
    {
        printf(
            "FAIL: non-finite projected error\n");

        return 1;
    }


    trace = 0.0;

    for(k=0;k<3;k++)
    {
        if( !isfinite(
                (double)TEST_REAL(eval[k])) )
        {
            printf(
                "FAIL: non-finite eigenvalue %d\n",
                (int)k);

            return 1;
        }


        if( fabs(
                (double)TEST_IMAG(eval[k])) >
            1.0e-12 )
        {
            printf(
                "FAIL: eigenvalue %d is not real\n",
                (int)k);

            return 1;
        }


        value[k] =
            (LIS_REAL)TEST_REAL(
                eval[k]);

        trace += value[k];


        dot = 0.0;

        for(i=0;i<3;i++)
        {
            dot +=
                TEST_CONJ(
                    evec[i*3+k]) *
                evec[i*3+k];
        }


        norm =
            (LIS_REAL)TEST_ABS(
                dot-1.0);

        if( norm>1.0e-11 )
        {
            printf(
                "FAIL: eigenvector %d norm %.17e\n",
                (int)k,
                (double)norm);

            return 1;
        }


        rnorm = 0.0;

        for(i=0;i<3;i++)
        {
            residual = 0.0;

            for(j=0;j<3;j++)
            {
                residual +=
                    original[i*3+j] *
                    evec[j*3+k];
            }

            residual -=
                eval[k] *
                evec[i*3+k];

            norm =
                (LIS_REAL)TEST_ABS(
                    residual);

            rnorm +=
                norm*norm;
        }


        rnorm =
            sqrt(rnorm);

        if( rnorm>1.0e-10 )
        {
            printf(
                "FAIL: Ritz residual %d %.17e\n",
                (int)k,
                (double)rnorm);

            return 1;
        }
    }


    for(k=0;k<3;k++)
    {
        for(l=k+1;l<3;l++)
        {
            dot = 0.0;

            for(i=0;i<3;i++)
            {
                dot +=
                    TEST_CONJ(
                        evec[i*3+k]) *
                    evec[i*3+l];
            }

            norm =
                (LIS_REAL)TEST_ABS(
                    dot);

            if( norm>1.0e-11 )
            {
                printf(
                    "FAIL: orthogonality %d %d %.17e\n",
                    (int)k,
                    (int)l,
                    (double)norm);

                return 1;
            }
        }
    }


    for(i=0;i<2;i++)
    {
        for(j=i+1;j<3;j++)
        {
            if( value[j]<value[i] )
            {
                tmp = value[i];
                value[i] = value[j];
                value[j] = tmp;
            }
        }
    }


    for(i=0;i<3;i++)
    {
        if( fabs(
                value[i]-expected[i]) >
            1.0e-10 )
        {
            printf(
                "FAIL: eigenvalue %d %.17e expected %.17e\n",
                (int)i,
                (double)value[i],
                (double)expected[i]);

            return 1;
        }
    }


    if( fabs(trace-9.0)>1.0e-11 )
    {
        printf(
            "FAIL: trace %.17e\n",
            (double)trace);

        return 1;
    }


    return 0;
}


int
main(int argc, char **argv)
{
    LIS_INT err;
    int failed;

    err = lis_initialize(
        &argc,
        &argv);

    if( err )
    {
        return (int)err;
    }


    failed =
        check_projected_system();


    lis_finalize();


    if( failed )
    {
        return 1;
    }


#ifdef _COMPLEX

    printf(
        "PASS: complex Hermitian projected eigensystem\n");

#else

    printf(
        "PASS: real symmetric projected eigensystem\n");

#endif


    return 0;
}
