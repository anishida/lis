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


#ifdef _COMPLEX

static LIS_REAL
general_complex_test_norm(
    LIS_INT n,
    const LIS_SCALAR *x)
{
    LIS_INT i;
    LIS_REAL sum;

    sum = 0.0;

    for(i=0;i<n;i++)
    {
        LIS_REAL v;

        v = cabs(x[i]);
        sum += v*v;
    }

    return sqrt(sum);
}


static int
check_complex_general_case(
    const char *label,
    LIS_INT n,
    const LIS_SCALAR *input)
{
    LIS_INT i,j,k;
    LIS_INT err;
    LIS_INT niter;

    LIS_REAL solver_error;
    LIS_REAL max_residual;
    LIS_REAL trace_error;

    LIS_SCALAR trace_original;
    LIS_SCALAR trace_eval;

    LIS_SCALAR *a;
    LIS_SCALAR *original;
    LIS_SCALAR *eval;
    LIS_SCALAR *evec;
    LIS_SCALAR *residual;

    a = NULL;
    original = NULL;
    eval = NULL;
    evec = NULL;
    residual = NULL;

    a =
        (LIS_SCALAR *)malloc(
            (size_t)n*n*sizeof(LIS_SCALAR));

    original =
        (LIS_SCALAR *)malloc(
            (size_t)n*n*sizeof(LIS_SCALAR));

    eval =
        (LIS_SCALAR *)malloc(
            (size_t)n*sizeof(LIS_SCALAR));

    evec =
        (LIS_SCALAR *)malloc(
            (size_t)n*n*sizeof(LIS_SCALAR));

    residual =
        (LIS_SCALAR *)malloc(
            (size_t)n*sizeof(LIS_SCALAR));

    if( a==NULL ||
        original==NULL ||
        eval==NULL ||
        evec==NULL ||
        residual==NULL )
    {
        free(residual);
        free(evec);
        free(eval);
        free(original);
        free(a);

        printf(
            "FAIL: %s allocation\n",
            label);

        return 1;
    }

    for(i=0;i<n*n;i++)
    {
        a[i] = input[i];
        original[i] = input[i];
    }

    niter = 0;
    solver_error = 0.0;

    err =
        lis_eigen_general_complex_solve(
            n,
            a,
            eval,
            evec,
            &niter,
            &solver_error);

    if( err!=LIS_SUCCESS )
    {
        printf(
            "FAIL: %s general complex solve returned %d\n",
            label,
            (int)err);

        free(residual);
        free(evec);
        free(eval);
        free(original);
        free(a);

        return 1;
    }

    if( niter<0 ||
        niter>20000 ||
        !isfinite((double)solver_error) ||
        solver_error>1.0e-8 )
    {
        printf(
            "FAIL: %s convergence niter=%d error=%.17e\n",
            label,
            (int)niter,
            (double)solver_error);

        free(residual);
        free(evec);
        free(eval);
        free(original);
        free(a);

        return 1;
    }

    trace_original = 0.0;
    trace_eval = 0.0;
    max_residual = 0.0;

    for(i=0;i<n;i++)
    {
        trace_original +=
            original[i+i*n];

        trace_eval +=
            eval[i];
    }

    trace_error =
        cabs(
            trace_original
            -trace_eval);

    if( trace_error>1.0e-8 )
    {
        printf(
            "FAIL: %s trace error %.17e\n",
            label,
            (double)trace_error);

        free(residual);
        free(evec);
        free(eval);
        free(original);
        free(a);

        return 1;
    }

    for(k=0;k<n;k++)
    {
        LIS_REAL nr;
        LIS_REAL nv;
        LIS_REAL scale;

        for(i=0;i<n;i++)
        {
            residual[i] =
                -eval[k]
                *evec[i+k*n];

            for(j=0;j<n;j++)
            {
                residual[i] +=
                    original[i+j*n]
                    *evec[j+k*n];
            }
        }

        nr =
            general_complex_test_norm(
                n,
                residual);

        nv =
            general_complex_test_norm(
                n,
                &evec[k*n]);

        scale =
            (
                1.0
                +cabs(eval[k])
            )*nv;

        if( scale<1.0 )
            scale = 1.0;

        nr /= scale;

        if( nr>max_residual )
            max_residual = nr;
    }

    if( max_residual>1.0e-8 )
    {
        printf(
            "FAIL: %s eigenpair residual %.17e\n",
            label,
            (double)max_residual);

        free(residual);
        free(evec);
        free(eval);
        free(original);
        free(a);

        return 1;
    }

    printf(
        "PASS: general complex eigensystem %s "
        "n=%d iter=%d error=%.6e "
        "residual=%.6e trace=%.6e\n",
        label,
        (int)n,
        (int)niter,
        (double)solver_error,
        (double)max_residual,
        (double)trace_error);

    free(residual);
    free(evec);
    free(eval);
    free(original);
    free(a);

    return 0;
}


static int
check_complex_general_eigensystem(void)
{
    static const LIS_SCALAR dense4[16] = {
         4.0 + 0.5*I,
         1.0 - 1.0*I,
         0.0 + 2.0*I,
        -1.0 + 0.5*I,

         2.0 + 0.2*I,
        -1.0 + 1.5*I,
         3.0 - 0.7*I,
         0.5 + 0.4*I,

        -0.5 + 1.0*I,
         1.5 + 0.3*I,
         2.0 - 2.0*I,
         4.0 + 0.6*I,

         1.0 - 0.8*I,
        -2.0 + 0.9*I,
         0.7 + 1.1*I,
         3.0 + 0.2*I
    };

    static const LIS_SCALAR clustered[16] = {
         1.0 + 0.2*I,
         0.0,
         0.0,
         0.0,

         0.7 - 0.3*I,
         1.00001 + 0.20002*I,
         0.0,
         0.0,

        -0.4 + 0.8*I,
         0.5 + 0.1*I,
         2.0 - 0.5*I,
         0.0,

         0.2 - 0.6*I,
        -0.3 + 0.7*I,
         0.9 + 0.4*I,
         3.0 + 0.1*I
    };

    int failed;

    failed = 0;

    failed |=
        check_complex_general_case(
            "dense-4x4",
            4,
            dense4);

    failed |=
        check_complex_general_case(
            "clustered-nonnormal-4x4",
            4,
            clustered);

    return failed;
}

#endif

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

    if( !failed )
    {
        failed =
            check_complex_general_eigensystem();
    }

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
