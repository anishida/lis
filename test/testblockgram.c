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


static LIS_INT
set_vector3(
    LIS_VECTOR v,
    LIS_SCALAR x0,
    LIS_SCALAR x1,
    LIS_SCALAR x2)
{
    LIS_INT err;

    err = lis_vector_set_value(
        LIS_INS_VALUE,
        0,
        x0,
        v);

    if( err )
        return err;

    err = lis_vector_set_value(
        LIS_INS_VALUE,
        1,
        x1,
        v);

    if( err )
        return err;

    return lis_vector_set_value(
        LIS_INS_VALUE,
        2,
        x2,
        v);
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

    return set_vector3(
        *v,
        x0,
        x1,
        x2);
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
check_value(
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
#ifdef _COMPLEX
        printf(
            "FAIL: %s error %.17e\n",
            name,
            (double)error);
#else
        printf(
            "FAIL: %s actual %.17e expected %.17e error %.17e\n",
            name,
            (double)actual,
            (double)expected,
            (double)error);
#endif

        return 1;
    }

    return 0;
}


static int
check_gram(void)
{
    LIS_VECTOR X[2];
    LIS_VECTOR Y[3];

    LIS_SCALAR gram[6];
    LIS_SCALAR expected[6];

    LIS_INT err;
    LIS_INT i;

    int failed;

    X[0] = NULL;
    X[1] = NULL;

    Y[0] = NULL;
    Y[1] = NULL;
    Y[2] = NULL;

    failed = 0;

#ifndef _COMPLEX

    err = create_vector3(
        &X[0],
        1.0,
        2.0,
        -1.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &X[1],
        0.0,
        1.0,
        3.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Y[0],
        2.0,
        0.0,
        1.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Y[1],
        -1.0,
        4.0,
        2.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Y[2],
        3.0,
        -2.0,
        1.0);

    if( err )
        goto cleanup;

    expected[0] = 1.0;
    expected[1] = 5.0;
    expected[2] = -2.0;

    expected[3] = 3.0;
    expected[4] = 10.0;
    expected[5] = 1.0;

#else

    err = create_vector3(
        &X[0],
        1.0 + 1.0*I,
        2.0,
        -1.0*I);

    if( err )
        goto cleanup;

    err = create_vector3(
        &X[1],
        0.0,
        1.0 - 2.0*I,
        3.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Y[0],
        2.0,
        1.0*I,
        1.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Y[1],
        -1.0,
        4.0 + 1.0*I,
        2.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Y[2],
        3.0,
        -2.0,
        1.0 - 1.0*I);

    if( err )
        goto cleanup;

    /*
     * Expected values use X^H Y.
     */
    expected[0] =
        conj(1.0 + 1.0*I)*2.0 +
        conj(2.0)*(1.0*I) +
        conj(-1.0*I)*1.0;

    expected[1] =
        conj(1.0 + 1.0*I)*(-1.0) +
        conj(2.0)*(4.0 + 1.0*I) +
        conj(-1.0*I)*2.0;

    expected[2] =
        conj(1.0 + 1.0*I)*3.0 +
        conj(2.0)*(-2.0) +
        conj(-1.0*I)*(1.0 - 1.0*I);

    expected[3] =
        conj(0.0)*2.0 +
        conj(1.0 - 2.0*I)*(1.0*I) +
        conj(3.0)*1.0;

    expected[4] =
        conj(0.0)*(-1.0) +
        conj(1.0 - 2.0*I)*(4.0 + 1.0*I) +
        conj(3.0)*2.0;

    expected[5] =
        conj(0.0)*3.0 +
        conj(1.0 - 2.0*I)*(-2.0) +
        conj(3.0)*(1.0 - 1.0*I);

#endif

    err = lis_eigen_block_gram(
        2,
        X,
        3,
        Y,
        gram);

    if( err )
    {
        printf(
            "FAIL: block Gram returned %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    for(i=0;i<6;i++)
    {
        char name[64];

        snprintf(
            name,
            sizeof(name),
            "gram[%d]",
            (int)i);

        if( check_value(
                name,
                gram[i],
                expected[i],
                1.0e-12) )
        {
            failed = 1;
            goto cleanup;
        }
    }

cleanup:

    for(i=0;i<2;i++)
    {
        if( X[i] )
            lis_vector_destroy(X[i]);
    }

    for(i=0;i<3;i++)
    {
        if( Y[i] )
            lis_vector_destroy(Y[i]);
    }

    return failed;
}


static int
check_bgram(void)
{
    LIS_MATRIX B;

    LIS_VECTOR X[2];
    LIS_VECTOR Y[2];
    LIS_VECTOR BY[2];

    LIS_SCALAR gram[4];
    LIS_SCALAR identity[4];

    LIS_SCALAR expected[4];
    LIS_SCALAR expected_identity[4];

    LIS_INT err;
    LIS_INT i;

    int failed;

    B = NULL;

    X[0] = NULL;
    X[1] = NULL;

    Y[0] = NULL;
    Y[1] = NULL;

    BY[0] = NULL;
    BY[1] = NULL;

    failed = 0;

    err = create_diagonal_B(&B);

    if( err )
        goto cleanup;

    err = create_vector3(
        &X[0],
        1.0,
        2.0,
        -1.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &X[1],
        0.0,
        1.0,
        3.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Y[0],
        2.0,
        0.0,
        1.0);

    if( err )
        goto cleanup;

    err = create_vector3(
        &Y[1],
        -1.0,
        4.0,
        2.0);

    if( err )
        goto cleanup;

    err = lis_vector_duplicate(
        Y[0],
        &BY[0]);

    if( err )
        goto cleanup;

    err = lis_vector_duplicate(
        Y[1],
        &BY[1]);

    if( err )
        goto cleanup;

    /*
     * B = diag(2,3,5)
     *
     * BY0 = [4,0,5]
     * BY1 = [-2,12,10]
     */
    expected[0] = -1.0;
    expected[1] = 12.0;
    expected[2] = 15.0;
    expected[3] = 42.0;

    expected_identity[0] = 1.0;
    expected_identity[1] = 5.0;
    expected_identity[2] = 3.0;
    expected_identity[3] = 10.0;

    err = lis_eigen_block_bgram(
        B,
        2,
        X,
        2,
        Y,
        BY,
        gram);

    if( err )
    {
        printf(
            "FAIL: block B-Gram returned %d\n",
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
            "bgram[%d]",
            (int)i);

        if( check_value(
                name,
                gram[i],
                expected[i],
                1.0e-12) )
        {
            failed = 1;
            goto cleanup;
        }
    }

    err = lis_eigen_block_bgram(
        NULL,
        2,
        X,
        2,
        Y,
        NULL,
        identity);

    if( err )
    {
        printf(
            "FAIL: identity B-Gram returned %d\n",
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
            "identity-bgram[%d]",
            (int)i);

        if( check_value(
                name,
                identity[i],
                expected_identity[i],
                1.0e-12) )
        {
            failed = 1;
            goto cleanup;
        }
    }

cleanup:

    if( B )
        lis_matrix_destroy(B);

    for(i=0;i<2;i++)
    {
        if( X[i] )
            lis_vector_destroy(X[i]);

        if( Y[i] )
            lis_vector_destroy(Y[i]);

        if( BY[i] )
            lis_vector_destroy(BY[i]);
    }

    return failed;
}


static int
check_argument_guards(void)
{
    LIS_SCALAR gram[1];

    LIS_INT err;

    err = lis_eigen_block_gram(
        0,
        NULL,
        0,
        NULL,
        gram);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "FAIL: Gram invalid-argument guard returned %d\n",
            (int)err);

        return 1;
    }

    err = lis_eigen_block_bgram(
        NULL,
        0,
        NULL,
        0,
        NULL,
        NULL,
        gram);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "FAIL: B-Gram invalid-argument guard returned %d\n",
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
        failed = check_gram();

    if( !failed )
        failed = check_bgram();

    if( !failed )
        failed = check_argument_guards();

    lis_finalize();

    if( failed )
        return 1;

    printf(
        "PASS: common block Gram and B-Gram algebra\n");

    return 0;
}
