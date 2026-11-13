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


static LIS_REAL
scalar_error(
    LIS_SCALAR a,
    LIS_SCALAR b)
{
#ifdef _COMPLEX
    return cabs(a-b);
#else
    return fabs(a-b);
#endif
}


static int
check_scalar(
    const char *name,
    LIS_SCALAR actual,
    LIS_SCALAR expected,
    LIS_REAL tol)
{
    LIS_REAL error;

    error = scalar_error(
        actual,
        expected);

    if( !isfinite((double)error) ||
        error>tol )
    {
        printf(
            "FAIL: %s error %.17e\n",
            name,
            (double)error);

        return 1;
    }

    return 0;
}


static LIS_INT
create_vector3(
    LIS_VECTOR *v,
    LIS_SCALAR x0,
    LIS_SCALAR x1,
    LIS_SCALAR x2)
{
    LIS_INT err;

    *v = NULL;

    err = lis_vector_create(
        LIS_COMM_WORLD,
        v);

    if( err )
        return err;

    err = lis_vector_set_size(
        *v,
        0,
        3);

    if( err )
        return err;

    err = lis_vector_set_value(
        LIS_INS_VALUE,
        0,
        x0,
        *v);

    if( err )
        return err;

    err = lis_vector_set_value(
        LIS_INS_VALUE,
        1,
        x1,
        *v);

    if( err )
        return err;

    return lis_vector_set_value(
        LIS_INS_VALUE,
        2,
        x2,
        *v);
}


static LIS_INT
create_diagonal_B(
    LIS_MATRIX *B)
{
    LIS_INT err;

    *B = NULL;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        B);

    if( err )
        return err;

    err = lis_matrix_set_size(
        *B,
        0,
        3);

    if( err )
        return err;

    err = lis_matrix_set_value(
        LIS_INS_VALUE,
        0,
        0,
        2.0,
        *B);

    if( err )
        return err;

    err = lis_matrix_set_value(
        LIS_INS_VALUE,
        1,
        1,
        3.0,
        *B);

    if( err )
        return err;

    err = lis_matrix_set_value(
        LIS_INS_VALUE,
        2,
        2,
        5.0,
        *B);

    if( err )
        return err;

    return lis_matrix_assemble(*B);
}


static int
read_vector3(
    LIS_VECTOR v,
    LIS_SCALAR values[3])
{
    LIS_INT err;

    err = lis_vector_get_values(
        v,
        0,
        3,
        values);

    if( err )
    {
        printf(
            "FAIL: lis_vector_get_values returned %d\n",
            (int)err);

        return 1;
    }

    return 0;
}


static int
check_zero_gram(
    LIS_INT nq,
    LIS_VECTOR Q[],
    LIS_INT nw,
    LIS_VECTOR W[],
    LIS_REAL tol)
{
    LIS_SCALAR gram[16];

    LIS_INT err;
    LIS_INT i;

    err = lis_eigen_block_gram(
        nq,
        Q,
        nw,
        W,
        gram);

    if( err )
    {
        printf(
            "FAIL: post-orthogonalization Gram returned %d\n",
            (int)err);

        return 1;
    }

    for(i=0;i<nq*nw;i++)
    {
        if( scalar_error(
                gram[i],
                0.0)>tol )
        {
            printf(
                "FAIL: residual Gram[%d] not zero\n",
                (int)i);

            return 1;
        }
    }

    return 0;
}


static int
check_zero_bgram(
    LIS_MATRIX B,
    LIS_INT nq,
    LIS_VECTOR Q[],
    LIS_INT nw,
    LIS_VECTOR W[],
    LIS_VECTOR BW[],
    LIS_REAL tol)
{
    LIS_SCALAR gram[16];

    LIS_INT err;
    LIS_INT i;

    err = lis_eigen_block_bgram(
        B,
        nq,
        Q,
        nw,
        W,
        BW,
        gram);

    if( err )
    {
        printf(
            "FAIL: post-B-orthogonalization Gram returned %d\n",
            (int)err);

        return 1;
    }

    for(i=0;i<nq*nw;i++)
    {
        if( scalar_error(
                gram[i],
                0.0)>tol )
        {
            printf(
                "FAIL: residual B-Gram[%d] not zero\n",
                (int)i);

            return 1;
        }
    }

    return 0;
}


static int
check_standard_orthogonalization(void)
{
    LIS_VECTOR Q[2];
    LIS_VECTOR W[2];

    LIS_SCALAR coeff[4];
    LIS_SCALAR work[4];

    LIS_SCALAR values[3];

    LIS_SCALAR expected[4];

    LIS_INT err;
    LIS_INT i;

    int failed;

    Q[0] = NULL;
    Q[1] = NULL;

    W[0] = NULL;
    W[1] = NULL;

    failed = 0;

#ifndef _COMPLEX

    err = create_vector3(
        &Q[0],
        1.0,
        0.0,
        0.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Q[1],
        0.0,
        1.0,
        0.0);

    if( err )
        goto cleanup;

    expected[0] = 2.0;
    expected[1] = -4.0;
    expected[2] = -1.0;
    expected[3] = 5.0;

#else

    {
        LIS_REAL s;

        s = 1.0/sqrt(2.0);

        err = create_vector3(
            &Q[0],
            s*(1.0 + 1.0*I),
            0.0,
            0.0);

        if( err )
            goto cleanup;

        err = create_vector3(
            &Q[1],
            0.0,
            1.0*I,
            0.0);

        if( err )
            goto cleanup;

        expected[0] =
            2.0*s*(1.0 - 1.0*I);

        expected[1] =
            -4.0*s*(1.0 - 1.0*I);

        expected[2] =
            1.0*I;

        expected[3] =
            -5.0*I;
    }

#endif

    err = create_vector3(
        &W[0],
        2.0,
        -1.0,
        3.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &W[1],
        -4.0,
        5.0,
        6.0);

    if( err )
        goto cleanup;

    err = lis_eigen_block_orthogonalize(
        2,
        Q,
        2,
        W,
        2,
        coeff,
        work);

    if( err )
    {
        printf(
            "FAIL: standard orthogonalization returned %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    for(i=0;i<4;i++)
    {
        char name[64];

        snprintf(
            name,
            sizeof(name),
            "standard-coeff[%d]",
            (int)i);

        if( check_scalar(
                name,
                coeff[i],
                expected[i],
                1.0e-11) )
        {
            failed = 1;
            goto cleanup;
        }
    }

    if( check_zero_gram(
            2,
            Q,
            2,
            W,
            1.0e-12) )
    {
        failed = 1;
        goto cleanup;
    }

    if( read_vector3(
            W[0],
            values) )
    {
        failed = 1;
        goto cleanup;
    }

    if( check_scalar(
            "W0[0]",
            values[0],
            0.0,
            1.0e-12) ||
        check_scalar(
            "W0[1]",
            values[1],
            0.0,
            1.0e-12) ||
        check_scalar(
            "W0[2]",
            values[2],
            3.0,
            1.0e-12) )
    {
        failed = 1;
        goto cleanup;
    }

    if( read_vector3(
            W[1],
            values) )
    {
        failed = 1;
        goto cleanup;
    }

    if( check_scalar(
            "W1[0]",
            values[0],
            0.0,
            1.0e-12) ||
        check_scalar(
            "W1[1]",
            values[1],
            0.0,
            1.0e-12) ||
        check_scalar(
            "W1[2]",
            values[2],
            6.0,
            1.0e-12) )
    {
        failed = 1;
        goto cleanup;
    }

cleanup:

    for(i=0;i<2;i++)
    {
        if( Q[i] )
            lis_vector_destroy(Q[i]);

        if( W[i] )
            lis_vector_destroy(W[i]);
    }

    return failed;
}


static int
check_b_orthogonalization(void)
{
    LIS_MATRIX B;

    LIS_VECTOR Q[2];
    LIS_VECTOR W[2];
    LIS_VECTOR BW[2];

    LIS_SCALAR coeff[4];
    LIS_SCALAR work[4];

    LIS_SCALAR values[3];

    LIS_SCALAR expected[4];

    LIS_INT err;
    LIS_INT i;

    int failed;

    B = NULL;

    Q[0] = NULL;
    Q[1] = NULL;

    W[0] = NULL;
    W[1] = NULL;

    BW[0] = NULL;
    BW[1] = NULL;

    failed = 0;

    err = create_diagonal_B(&B);

    if( err )
        goto cleanup;

#ifndef _COMPLEX

    err = create_vector3(
        &Q[0],
        1.0/sqrt(2.0),
        0.0,
        0.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Q[1],
        0.0,
        1.0/sqrt(3.0),
        0.0);

    if( err )
        goto cleanup;

    expected[0] = 2.0*sqrt(2.0);
    expected[1] = -4.0*sqrt(2.0);
    expected[2] = -sqrt(3.0);
    expected[3] = 5.0*sqrt(3.0);

#else

    err = create_vector3(
        &Q[0],
        0.5*(1.0 + 1.0*I),
        0.0,
        0.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Q[1],
        0.0,
        (1.0*I)/sqrt(3.0),
        0.0);

    if( err )
        goto cleanup;

    expected[0] =
        2.0*(1.0 - 1.0*I);

    expected[1] =
        -4.0*(1.0 - 1.0*I);

    expected[2] =
        sqrt(3.0)*1.0*I;

    expected[3] =
        -5.0*sqrt(3.0)*1.0*I;

#endif

    err = create_vector3(
        &W[0],
        2.0,
        -1.0,
        3.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &W[1],
        -4.0,
        5.0,
        6.0);

    if( err )
        goto cleanup;

    err = lis_vector_duplicate(
        W[0],
        &BW[0]);

    if( err )
        goto cleanup;

    err = lis_vector_duplicate(
        W[1],
        &BW[1]);

    if( err )
        goto cleanup;

    err = lis_eigen_block_borthogonalize(
        B,
        2,
        Q,
        2,
        W,
        BW,
        2,
        coeff,
        work);

    if( err )
    {
        printf(
            "FAIL: B-orthogonalization returned %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    for(i=0;i<4;i++)
    {
        char name[64];

        snprintf(
            name,
            sizeof(name),
            "B-coeff[%d]",
            (int)i);

        if( check_scalar(
                name,
                coeff[i],
                expected[i],
                1.0e-11) )
        {
            failed = 1;
            goto cleanup;
        }
    }

    if( check_zero_bgram(
            B,
            2,
            Q,
            2,
            W,
            BW,
            1.0e-12) )
    {
        failed = 1;
        goto cleanup;
    }

    if( read_vector3(
            W[0],
            values) )
    {
        failed = 1;
        goto cleanup;
    }

    if( check_scalar(
            "BW-projected-W0[0]",
            values[0],
            0.0,
            1.0e-12) ||
        check_scalar(
            "BW-projected-W0[1]",
            values[1],
            0.0,
            1.0e-12) ||
        check_scalar(
            "BW-projected-W0[2]",
            values[2],
            3.0,
            1.0e-12) )
    {
        failed = 1;
        goto cleanup;
    }

    if( read_vector3(
            BW[0],
            values) )
    {
        failed = 1;
        goto cleanup;
    }

    if( check_scalar(
            "BW0[0]",
            values[0],
            0.0,
            1.0e-11) ||
        check_scalar(
            "BW0[1]",
            values[1],
            0.0,
            1.0e-11) ||
        check_scalar(
            "BW0[2]",
            values[2],
            15.0,
            1.0e-11) )
    {
        failed = 1;
        goto cleanup;
    }

cleanup:

    if( B )
        lis_matrix_destroy(B);

    for(i=0;i<2;i++)
    {
        if( Q[i] )
            lis_vector_destroy(Q[i]);

        if( W[i] )
            lis_vector_destroy(W[i]);

        if( BW[i] )
            lis_vector_destroy(BW[i]);
    }

    return failed;
}


static int
check_identity_fallback(void)
{
    LIS_VECTOR Q[1];
    LIS_VECTOR W[1];

    LIS_SCALAR coeff[1];
    LIS_SCALAR work[1];

    LIS_SCALAR values[3];

    LIS_INT err;

    int failed;

    Q[0] = NULL;
    W[0] = NULL;

    failed = 0;

    err = create_vector3(
        &Q[0],
        1.0,
        0.0,
        0.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &W[0],
        7.0,
        2.0,
        3.0);

    if( err )
        goto cleanup;

    err = lis_eigen_block_borthogonalize(
        NULL,
        1,
        Q,
        1,
        W,
        NULL,
        2,
        coeff,
        work);

    if( err )
    {
        printf(
            "FAIL: identity B fallback returned %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    if( check_scalar(
            "identity coeff",
            coeff[0],
            7.0,
            1.0e-12) )
    {
        failed = 1;
        goto cleanup;
    }

    if( read_vector3(
            W[0],
            values) )
    {
        failed = 1;
        goto cleanup;
    }

    if( check_scalar(
            "identity W[0]",
            values[0],
            0.0,
            1.0e-12) ||
        check_scalar(
            "identity W[1]",
            values[1],
            2.0,
            1.0e-12) ||
        check_scalar(
            "identity W[2]",
            values[2],
            3.0,
            1.0e-12) )
    {
        failed = 1;
        goto cleanup;
    }

cleanup:

    if( Q[0] )
        lis_vector_destroy(Q[0]);

    if( W[0] )
        lis_vector_destroy(W[0]);

    return failed;
}


static int
check_argument_guards(void)
{
    LIS_SCALAR coeff[1];
    LIS_SCALAR work[1];

    LIS_INT err;

    err = lis_eigen_block_orthogonalize(
        0,
        NULL,
        0,
        NULL,
        0,
        coeff,
        work);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "FAIL: standard invalid-argument guard returned %d\n",
            (int)err);

        return 1;
    }

    err = lis_eigen_block_orthogonalize(
        1,
        NULL,
        1,
        NULL,
        1,
        work,
        work);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "FAIL: workspace-alias guard returned %d\n",
            (int)err);

        return 1;
    }

    return 0;
}


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
        return (int)err;

    if( !failed )
        failed =
            check_standard_orthogonalization();

    if( !failed )
        failed =
            check_b_orthogonalization();

    if( !failed )
        failed =
            check_identity_fallback();

    if( !failed )
        failed =
            check_argument_guards();

    lis_finalize();

    if( failed )
        return 1;

    printf(
        "PASS: common block orthogonalization\n");

    return 0;
}
