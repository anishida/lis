#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <stdio.h>
#include <math.h>

#include "lis.h"
#include "lis_esolver.h"


#ifndef _COMPLEX

static int
check_real_projected_eigensystem(void)
{
    static const LIS_SCALAR original[9] = {
        4.0, 1.0, 0.0,
        1.0, 3.0, 1.0,
        0.0, 1.0, 2.0
    };

    LIS_SCALAR a[9] = {
        4.0, 1.0, 0.0,
        1.0, 3.0, 1.0,
        0.0, 1.0, 2.0
    };

    LIS_SCALAR eval[3];
    LIS_SCALAR evec[9];

    LIS_INT err;
    LIS_INT niter;

    LIS_REAL error;
    LIS_REAL trace;
    LIS_REAL residual2;
    LIS_REAL residual;
    LIS_REAL norm2;
    LIS_REAL norm;
    LIS_REAL dot;

    LIS_INT i,j,k,p,q;

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
            "FAIL: invalid Jacobi sweep count %d\n",
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
        if( !isfinite((double)eval[k]) )
        {
            printf(
                "FAIL: non-finite eigenvalue %d\n",
                (int)k);

            return 1;
        }

        trace += eval[k];

        norm2 = 0.0;

        for(i=0;i<3;i++)
        {
            norm2 +=
                evec[i*3+k] *
                evec[i*3+k];
        }

        norm = sqrt(norm2);

        if( !isfinite((double)norm) ||
            fabs(norm-1.0)>1.0e-10 )
        {
            printf(
                "FAIL: eigenvector %d norm %.17e\n",
                (int)k,
                (double)norm);

            return 1;
        }

        residual2 = 0.0;

        for(i=0;i<3;i++)
        {
            LIS_SCALAR value;

            value = 0.0;

            for(j=0;j<3;j++)
            {
                value +=
                    original[i*3+j] *
                    evec[j*3+k];
            }

            value -=
                eval[k] *
                evec[i*3+k];

            residual2 += value*value;
        }

        residual = sqrt(residual2);

        if( !isfinite((double)residual) ||
            residual>1.0e-10 )
        {
            printf(
                "FAIL: Ritz residual %d %.17e\n",
                (int)k,
                (double)residual);

            return 1;
        }
    }

    if( fabs(trace-9.0)>1.0e-10 )
    {
        printf(
            "FAIL: eigenvalue trace %.17e\n",
            (double)trace);

        return 1;
    }

    for(p=0;p<3;p++)
    {
        for(q=0;q<3;q++)
        {
            dot = 0.0;

            for(i=0;i<3;i++)
            {
                dot +=
                    evec[i*3+p] *
                    evec[i*3+q];
            }

            if( p==q )
            {
                dot -= 1.0;
            }

            if( fabs(dot)>1.0e-10 )
            {
                printf(
                    "FAIL: eigenvector orthogonality (%d,%d) %.17e\n",
                    (int)p,
                    (int)q,
                    (double)dot);

                return 1;
            }
        }
    }

    return 0;
}

#else

static int
check_complex_projected_eigensystem(void)
{
    LIS_SCALAR a[4];
    LIS_SCALAR eval[2];
    LIS_SCALAR evec[4];

    LIS_INT err;
    LIS_INT niter;

    LIS_REAL error;

    a[0] = 2.0;
    a[1] = 1.0;
    a[2] = 1.0;
    a[3] = 3.0;

    niter = 0;
    error = 0.0;

    err = lis_eigen_projected_solve(
        2,
        a,
        eval,
        evec,
        &niter,
        &error);

    if( err!=LIS_SUCCESS )
    {
        printf(
            "FAIL: complex projected solve returned %d\n",
            (int)err);

        return 1;
    }

    return 0;
}

#endif


int
main(int argc, char **argv)
{
    LIS_INT err;
    int failed;

    failed = 0;

    err = lis_initialize(
        &argc,
        &argv);

    if( err )
    {
        return (int)err;
    }

#ifndef _COMPLEX

    failed =
        check_real_projected_eigensystem();

#else

    failed =
        check_complex_projected_eigensystem();

#endif

    lis_finalize();

    if( failed )
    {
        return 1;
    }

#ifndef _COMPLEX

    printf(
        "PASS: common projected eigensystem solve\n");

#else

    printf(
        "PASS: common projected eigensystem complex path\n");

#endif

    return 0;
}
