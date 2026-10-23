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


#define CAPACITY 4
#define NOMINAL 3
#define RANK_TOL ((LIS_REAL)1.0e-12)
#define TEST_TOL ((LIS_REAL)2.0e-10)


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


static LIS_REAL
scalar_real(
    LIS_SCALAR value)
{
#ifdef _COMPLEX
    return (LIS_REAL)creal(value);
#else
    return (LIS_REAL)value;
#endif
}


static LIS_SCALAR
scalar_conjugate(
    LIS_SCALAR value)
{
#ifdef _COMPLEX
    return conj(value);
#else
    return value;
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
read_vector3(
    LIS_VECTOR v,
    LIS_SCALAR value[3])
{
    LIS_INT err;

    err = lis_vector_get_values(
        v,
        0,
        3,
        value);

    if( err )
    {
        printf(
            "FAIL: vector read returned %d\n",
            (int)err);

        return 1;
    }

    return 0;
}


static int
check_zero_vector(
    const char *name,
    LIS_VECTOR v,
    LIS_REAL tol)
{
    LIS_SCALAR value[3];
    LIS_INT i;

    if( read_vector3(
            v,
            value) )
    {
        return 1;
    }

    for(i=0;i<3;i++)
    {
        if( scalar_error(
                value[i],
                0.0)>tol )
        {
            printf(
                "FAIL: %s[%d] is not zero\n",
                name,
                (int)i);

            return 1;
        }
    }

    return 0;
}


static int
check_identity_gram(
    LIS_INT rank,
    LIS_VECTOR Q[],
    LIS_REAL tol)
{
    LIS_SCALAR gram[CAPACITY*CAPACITY];

    LIS_INT err;
    LIS_INT i,j;

    err = lis_eigen_block_gram(
        rank,
        Q,
        rank,
        Q,
        gram);

    if( err )
    {
        printf(
            "FAIL: identity Gram returned %d\n",
            (int)err);

        return 1;
    }

    for(i=0;i<rank;i++)
    {
        for(j=0;j<rank;j++)
        {
            LIS_SCALAR expected;

            expected =
                i==j ? 1.0 : 0.0;

            if( scalar_error(
                    gram[i*rank+j],
                    expected)>tol )
            {
                printf(
                    "FAIL: Gram[%d,%d] error %.17e\n",
                    (int)i,
                    (int)j,
                    (double)scalar_error(
                        gram[i*rank+j],
                        expected));

                return 1;
            }
        }
    }

    return 0;
}


static int
check_identity_bgram_from_bq(
    LIS_INT rank,
    LIS_VECTOR Q[],
    LIS_VECTOR BQ[],
    LIS_REAL tol)
{
    LIS_SCALAR gram[CAPACITY*CAPACITY];

    LIS_INT err;
    LIS_INT i,j;

    err = lis_eigen_block_gram(
        rank,
        Q,
        rank,
        BQ,
        gram);

    if( err )
    {
        printf(
            "FAIL: Q^H BQ Gram returned %d\n",
            (int)err);

        return 1;
    }

    for(i=0;i<rank;i++)
    {
        for(j=0;j<rank;j++)
        {
            LIS_SCALAR expected;

            expected =
                i==j ? 1.0 : 0.0;

            if( scalar_error(
                    gram[i*rank+j],
                    expected)>tol )
            {
                printf(
                    "FAIL: B-Gram[%d,%d] error %.17e\n",
                    (int)i,
                    (int)j,
                    (double)scalar_error(
                        gram[i*rank+j],
                        expected));

                return 1;
            }
        }
    }

    return 0;
}


static int
check_linear_combination(
    LIS_VECTOR X[],
    LIS_VECTOR Y[])
{
    LIS_SCALAR coeff[4];
    LIS_SCALAR value[3];

    LIS_INT err;

#ifndef _COMPLEX

    coeff[0] = 1.0;
    coeff[1] = 2.0;
    coeff[2] = 3.0;
    coeff[3] = 4.0;

#else

    coeff[0] = 1.0 + 1.0*I;
    coeff[1] = 2.0;
    coeff[2] = 3.0;
    coeff[3] = 4.0 - 1.0*I;

#endif

    err = set_vector3(
        X[0],
        1.0,
        0.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        X[1],
        0.0,
        1.0,
        0.0);

    if( err )
        return 1;

    err = lis_eigen_block_linear_combination(
        2,
        X,
        2,
        Y,
        coeff,
        2);

    if( err )
    {
        printf(
            "FAIL: linear combination returned %d\n",
            (int)err);

        return 1;
    }

    if( read_vector3(
            Y[0],
            value) )
    {
        return 1;
    }

#ifdef _COMPLEX

    if( scalar_error(
            value[0],
            1.0 + 1.0*I)>TEST_TOL ||
        scalar_error(
            value[1],
            3.0)>TEST_TOL ||
        scalar_error(
            value[2],
            0.0)>TEST_TOL )
    {
        printf(
            "FAIL: complex linear combination column 0\n");

        return 1;
    }

#else

    if( scalar_error(
            value[0],
            1.0)>TEST_TOL ||
        scalar_error(
            value[1],
            3.0)>TEST_TOL ||
        scalar_error(
            value[2],
            0.0)>TEST_TOL )
    {
        printf(
            "FAIL: linear combination column 0\n");

        return 1;
    }

#endif

    if( read_vector3(
            Y[1],
            value) )
    {
        return 1;
    }

#ifdef _COMPLEX

    if( scalar_error(
            value[0],
            2.0)>TEST_TOL ||
        scalar_error(
            value[1],
            4.0 - 1.0*I)>TEST_TOL ||
        scalar_error(
            value[2],
            0.0)>TEST_TOL )
    {
        printf(
            "FAIL: complex linear combination column 1\n");

        return 1;
    }

#else

    if( scalar_error(
            value[0],
            2.0)>TEST_TOL ||
        scalar_error(
            value[1],
            4.0)>TEST_TOL ||
        scalar_error(
            value[2],
            0.0)>TEST_TOL )
    {
        printf(
            "FAIL: linear combination column 1\n");

        return 1;
    }

#endif

    return 0;
}


static int
check_standard_rank(
    LIS_VECTOR W[],
    LIS_VECTOR Q[])
{
    LIS_SCALAR gram[NOMINAL*NOMINAL];
    LIS_SCALAR eval[NOMINAL];
    LIS_SCALAR evec[NOMINAL*NOMINAL];
    LIS_SCALAR transform[NOMINAL*NOMINAL];

    LIS_SCALAR value[3];

    LIS_INT active_rank;
    LIS_INT err;
    LIS_INT i;

#ifndef _COMPLEX

    err = set_vector3(
        W[0],
        1.0,
        0.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[1],
        0.0,
        1.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[2],
        1.0,
        2.0,
        0.0);

#else

    err = set_vector3(
        W[0],
        1.0,
        0.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[1],
        0.0,
        1.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[2],
        1.0 + 1.0*I,
        2.0 - 1.0*I,
        0.0);

#endif

    if( err )
        return 1;

    active_rank = -1;

    err = lis_eigen_block_rank_deflate(
        CAPACITY,
        NOMINAL,
        W,
        Q,
        RANK_TOL,
        &active_rank,
        gram,
        eval,
        evec,
        transform);

    if( err )
    {
        printf(
            "FAIL: standard rank deflation returned %d\n",
            (int)err);

        return 1;
    }

    if( active_rank!=2 )
    {
        printf(
            "FAIL: standard active_rank=%d expected 2\n",
            (int)active_rank);

        return 1;
    }

    if( check_identity_gram(
            active_rank,
            Q,
            TEST_TOL) )
    {
        return 1;
    }

    for(i=0;i<active_rank;i++)
    {
        if( read_vector3(
                Q[i],
                value) )
        {
            return 1;
        }

        if( scalar_error(
                value[2],
                0.0)>TEST_TOL )
        {
            printf(
                "FAIL: standard basis left candidate span\n");

            return 1;
        }
    }

    if( check_zero_vector(
            "Q inactive 2",
            Q[2],
            TEST_TOL) ||
        check_zero_vector(
            "Q inactive 3",
            Q[3],
            TEST_TOL) )
    {
        return 1;
    }

    if( scalar_real(eval[0]) <
        scalar_real(eval[1]) )
    {
        printf(
            "FAIL: rank eigenvalues are not descending\n");

        return 1;
    }

    return 0;
}


static int
check_threshold_rank(
    LIS_VECTOR W[],
    LIS_VECTOR Q[])
{
    LIS_SCALAR gram[NOMINAL*NOMINAL];
    LIS_SCALAR eval[NOMINAL];
    LIS_SCALAR evec[NOMINAL*NOMINAL];
    LIS_SCALAR transform[NOMINAL*NOMINAL];

    LIS_INT active_rank;
    LIS_INT err;

    err = set_vector3(
        W[0],
        1.0,
        0.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[1],
        0.0,
        1.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[2],
        0.0,
        0.0,
        1.0e-8);

    if( err )
        return 1;

    active_rank = -1;

    err = lis_eigen_block_rank_deflate(
        CAPACITY,
        NOMINAL,
        W,
        Q,
        RANK_TOL,
        &active_rank,
        gram,
        eval,
        evec,
        transform);

    if( err )
    {
        printf(
            "FAIL: threshold rank deflation returned %d\n",
            (int)err);

        return 1;
    }

    if( active_rank!=2 )
    {
        printf(
            "FAIL: threshold active_rank=%d expected 2\n",
            (int)active_rank);

        return 1;
    }

    if( check_identity_gram(
            active_rank,
            Q,
            TEST_TOL) )
    {
        return 1;
    }

    if( scalar_real(eval[2]) >
        1.0e-12 )
    {
        printf(
            "FAIL: expected small discarded eigenvalue\n");

        return 1;
    }

    return 0;
}


static int
check_zero_rank(
    LIS_VECTOR W[],
    LIS_VECTOR Q[])
{
    LIS_SCALAR gram[NOMINAL*NOMINAL];
    LIS_SCALAR eval[NOMINAL];
    LIS_SCALAR evec[NOMINAL*NOMINAL];
    LIS_SCALAR transform[NOMINAL*NOMINAL];

    LIS_INT active_rank;
    LIS_INT err;
    LIS_INT i;

    for(i=0;i<NOMINAL;i++)
    {
        err = set_vector3(
            W[i],
            0.0,
            0.0,
            0.0);

        if( err )
            return 1;
    }

    active_rank = -1;

    err = lis_eigen_block_rank_deflate(
        CAPACITY,
        NOMINAL,
        W,
        Q,
        RANK_TOL,
        &active_rank,
        gram,
        eval,
        evec,
        transform);

    if( err )
    {
        printf(
            "FAIL: zero-rank deflation returned %d\n",
            (int)err);

        return 1;
    }

    if( active_rank!=0 )
    {
        printf(
            "FAIL: zero block active_rank=%d expected 0\n",
            (int)active_rank);

        return 1;
    }

    for(i=0;i<CAPACITY;i++)
    {
        if( check_zero_vector(
                "zero-rank Q",
                Q[i],
                TEST_TOL) )
        {
            return 1;
        }
    }

    return 0;
}


static int
check_b_rank(
    LIS_MATRIX B,
    LIS_VECTOR W[],
    LIS_VECTOR BW[],
    LIS_VECTOR Q[],
    LIS_VECTOR BQ[])
{
    LIS_SCALAR gram[NOMINAL*NOMINAL];
    LIS_SCALAR eval[NOMINAL];
    LIS_SCALAR evec[NOMINAL*NOMINAL];
    LIS_SCALAR transform[NOMINAL*NOMINAL];

    LIS_SCALAR qvalue[3];
    LIS_SCALAR bqvalue[3];

    LIS_INT active_rank;
    LIS_INT err;
    LIS_INT i;

#ifndef _COMPLEX

    err = set_vector3(
        W[0],
        1.0,
        0.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[1],
        0.0,
        1.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[2],
        1.0,
        2.0,
        0.0);

#else

    err = set_vector3(
        W[0],
        1.0,
        0.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[1],
        0.0,
        1.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[2],
        1.0 + 1.0*I,
        2.0 - 1.0*I,
        0.0);

#endif

    if( err )
        return 1;

    active_rank = -1;

    err = lis_eigen_block_brank_deflate(
        B,
        CAPACITY,
        NOMINAL,
        W,
        BW,
        Q,
        BQ,
        RANK_TOL,
        &active_rank,
        gram,
        eval,
        evec,
        transform);

    if( err )
    {
        printf(
            "FAIL: B-rank deflation returned %d\n",
            (int)err);

        return 1;
    }

    if( active_rank!=2 )
    {
        printf(
            "FAIL: B active_rank=%d expected 2\n",
            (int)active_rank);

        return 1;
    }

    if( check_identity_bgram_from_bq(
            active_rank,
            Q,
            BQ,
            TEST_TOL) )
    {
        return 1;
    }

    for(i=0;i<active_rank;i++)
    {
        if( read_vector3(
                Q[i],
                qvalue) ||
            read_vector3(
                BQ[i],
                bqvalue) )
        {
            return 1;
        }

        if( scalar_error(
                bqvalue[0],
                2.0*qvalue[0])>TEST_TOL ||
            scalar_error(
                bqvalue[1],
                3.0*qvalue[1])>TEST_TOL ||
            scalar_error(
                bqvalue[2],
                5.0*qvalue[2])>TEST_TOL )
        {
            printf(
                "FAIL: BQ is not synchronized with B*Q\n");

            return 1;
        }
    }

    if( check_zero_vector(
            "BQ inactive 2",
            BQ[2],
            TEST_TOL) ||
        check_zero_vector(
            "BQ inactive 3",
            BQ[3],
            TEST_TOL) )
    {
        return 1;
    }

    return 0;
}


static int
check_ritz_reconstruction(
    LIS_VECTOR basis[],
    LIS_VECTOR ritz[])
{
    LIS_SCALAR coeff[6];
    LIS_SCALAR scalar_coeff[1];

    LIS_SCALAR value[3];

    LIS_INT err;


    /*
     * Use a three-vector coordinate basis and reconstruct two
     * Ritz vectors.  This also verifies nritz < basis_rank.
     */
    err = set_vector3(
        basis[0],
        1.0,
        0.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        basis[1],
        0.0,
        1.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        basis[2],
        0.0,
        0.0,
        1.0);

    if( err )
        return 1;


#ifndef _COMPLEX

    coeff[0] = 1.0;
    coeff[1] = 2.0;

    coeff[2] = 3.0;
    coeff[3] = 4.0;

    coeff[4] = 5.0;
    coeff[5] = 6.0;

#else

    coeff[0] = 1.0 + 1.0*I;
    coeff[1] = 2.0;

    coeff[2] = 3.0;
    coeff[3] = 4.0 - 1.0*I;

    coeff[4] = 5.0 + 2.0*I;
    coeff[5] = 6.0;

#endif


    err = lis_eigen_ritz_reconstruct(
        3,
        basis,
        2,
        coeff,
        2,
        ritz);

    if( err )
    {
        printf(
            "FAIL: Ritz reconstruction returned %d\n",
            (int)err);

        return 1;
    }


    if( read_vector3(
            ritz[0],
            value) )
    {
        return 1;
    }


#ifndef _COMPLEX

    if( scalar_error(
            value[0],
            1.0)>TEST_TOL ||
        scalar_error(
            value[1],
            3.0)>TEST_TOL ||
        scalar_error(
            value[2],
            5.0)>TEST_TOL )
    {
        printf(
            "FAIL: REAL Ritz vector 0\n");

        return 1;
    }

#else

    if( scalar_error(
            value[0],
            1.0 + 1.0*I)>TEST_TOL ||
        scalar_error(
            value[1],
            3.0)>TEST_TOL ||
        scalar_error(
            value[2],
            5.0 + 2.0*I)>TEST_TOL )
    {
        printf(
            "FAIL: COMPLEX Ritz vector 0\n");

        return 1;
    }

#endif


    if( read_vector3(
            ritz[1],
            value) )
    {
        return 1;
    }


#ifndef _COMPLEX

    if( scalar_error(
            value[0],
            2.0)>TEST_TOL ||
        scalar_error(
            value[1],
            4.0)>TEST_TOL ||
        scalar_error(
            value[2],
            6.0)>TEST_TOL )
    {
        printf(
            "FAIL: REAL Ritz vector 1\n");

        return 1;
    }

#else

    if( scalar_error(
            value[0],
            2.0)>TEST_TOL ||
        scalar_error(
            value[1],
            4.0 - 1.0*I)>TEST_TOL ||
        scalar_error(
            value[2],
            6.0)>TEST_TOL )
    {
        printf(
            "FAIL: COMPLEX Ritz vector 1\n");

        return 1;
    }

#endif


    /*
     * Scalar basis_rank=1 compatibility.
     */
    err = set_vector3(
        basis[0],
        2.0,
        -1.0,
        0.5);

    if( err )
        return 1;


#ifndef _COMPLEX

    scalar_coeff[0] = 2.0;

#else

    scalar_coeff[0] = 2.0 - 1.0*I;

#endif


    err = lis_eigen_ritz_reconstruct(
        1,
        basis,
        1,
        scalar_coeff,
        1,
        ritz);

    if( err )
    {
        printf(
            "FAIL: scalar Ritz reconstruction returned %d\n",
            (int)err);

        return 1;
    }


    if( read_vector3(
            ritz[0],
            value) )
    {
        return 1;
    }


#ifndef _COMPLEX

    if( scalar_error(
            value[0],
            4.0)>TEST_TOL ||
        scalar_error(
            value[1],
            -2.0)>TEST_TOL ||
        scalar_error(
            value[2],
            1.0)>TEST_TOL )
    {
        printf(
            "FAIL: scalar REAL Ritz reconstruction\n");

        return 1;
    }

#else

    if( scalar_error(
            value[0],
            4.0 - 2.0*I)>TEST_TOL ||
        scalar_error(
            value[1],
            -2.0 + 1.0*I)>TEST_TOL ||
        scalar_error(
            value[2],
            1.0 - 0.5*I)>TEST_TOL )
    {
        printf(
            "FAIL: scalar COMPLEX Ritz reconstruction\n");

        return 1;
    }

#endif


    /*
     * API guard: zero basis rank.
     */
    err = lis_eigen_ritz_reconstruct(
        0,
        basis,
        1,
        scalar_coeff,
        1,
        ritz);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "FAIL: Ritz basis-rank guard returned %d\n",
            (int)err);

        return 1;
    }


    /*
     * API guard: coefficient leading dimension too small.
     */
    err = lis_eigen_ritz_reconstruct(
        3,
        basis,
        2,
        coeff,
        1,
        ritz);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "FAIL: Ritz leading-dimension guard returned %d\n",
            (int)err);

        return 1;
    }


    /*
     * Out-of-place requirement is inherited from the shared
     * linear-combination kernel.
     */
    err = lis_eigen_ritz_reconstruct(
        1,
        basis,
        1,
        scalar_coeff,
        1,
        basis);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "FAIL: Ritz alias guard returned %d\n",
            (int)err);

        return 1;
    }


    return 0;
}


static int
check_rectangular_coupling(
    LIS_VECTOR W[],
    LIS_VECTOR Q[],
    LIS_VECTOR R[])
{
    LIS_SCALAR gram[NOMINAL*NOMINAL];
    LIS_SCALAR gram_reference[NOMINAL*NOMINAL];

    LIS_SCALAR eval[NOMINAL];
    LIS_SCALAR evec[NOMINAL*NOMINAL];
    LIS_SCALAR transform[NOMINAL*NOMINAL];

    LIS_SCALAR coupling[NOMINAL*NOMINAL];

    LIS_SCALAR original_value[3];
    LIS_SCALAR rebuilt_value[3];

    LIS_SCALAR product;

    LIS_INT active_rank;
    LIS_INT err;

    LIS_INT i,j,k;


#ifndef _COMPLEX

    err = set_vector3(
        W[0],
        1.0,
        0.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[1],
        0.0,
        1.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[2],
        1.0,
        2.0,
        0.0);

#else

    err = set_vector3(
        W[0],
        1.0,
        0.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[1],
        0.0,
        1.0,
        0.0);

    if( err )
        return 1;

    err = set_vector3(
        W[2],
        1.0 + 1.0*I,
        2.0 - 1.0*I,
        0.0);

#endif

    if( err )
        return 1;


    err = lis_eigen_block_gram(
        NOMINAL,
        W,
        NOMINAL,
        W,
        gram_reference);

    if( err )
    {
        printf(
            "FAIL: coupling reference Gram returned %d\n",
            (int)err);

        return 1;
    }


    active_rank = -1;

    err = lis_eigen_block_rank_deflate(
        CAPACITY,
        NOMINAL,
        W,
        Q,
        RANK_TOL,
        &active_rank,
        gram,
        eval,
        evec,
        transform);

    if( err )
    {
        printf(
            "FAIL: coupling rank deflation returned %d\n",
            (int)err);

        return 1;
    }


    if( active_rank!=2 )
    {
        printf(
            "FAIL: coupling active_rank=%d expected 2\n",
            (int)active_rank);

        return 1;
    }


    err = lis_eigen_block_rectangular_coupling(
        NOMINAL,
        active_rank,
        eval,
        evec,
        coupling);

    if( err )
    {
        printf(
            "FAIL: rectangular coupling returned %d\n",
            (int)err);

        return 1;
    }


    /*
     * Inactive rows must remain deterministic zero storage.
     */
    for(j=active_rank;j<NOMINAL;j++)
    {
        for(i=0;i<NOMINAL;i++)
        {
            if( scalar_error(
                    coupling[j*NOMINAL+i],
                    0.0)>TEST_TOL )
            {
                printf(
                    "FAIL: inactive coupling row is not zero\n");

                return 1;
            }
        }
    }


    /*
     * C^H C must reproduce the original Gram matrix:
     *
     *     C = Lambda_r^(1/2) U_r^H
     *
     *     C^H C = U_r Lambda_r U_r^H
     */
    for(i=0;i<NOMINAL;i++)
    {
        for(j=0;j<NOMINAL;j++)
        {
            product = 0.0;

            for(k=0;k<active_rank;k++)
            {
                product +=
                    scalar_conjugate(
                        coupling[k*NOMINAL+i]) *
                    coupling[k*NOMINAL+j];
            }

            if( scalar_error(
                    product,
                    gram_reference[i*NOMINAL+j]) >
                TEST_TOL )
            {
                printf(
                    "FAIL: coupling Gram factorization "
                    "[%d,%d] error %.17e\n",
                    (int)i,
                    (int)j,
                    (double)scalar_error(
                        product,
                        gram_reference[i*NOMINAL+j]));

                return 1;
            }
        }
    }


    /*
     * Reconstruct every original candidate:
     *
     *     W = Q C
     */
    err = lis_eigen_block_linear_combination(
        active_rank,
        Q,
        NOMINAL,
        R,
        coupling,
        NOMINAL);

    if( err )
    {
        printf(
            "FAIL: coupling reconstruction returned %d\n",
            (int)err);

        return 1;
    }


    for(j=0;j<NOMINAL;j++)
    {
        if( read_vector3(
                W[j],
                original_value) ||
            read_vector3(
                R[j],
                rebuilt_value) )
        {
            return 1;
        }

        for(i=0;i<3;i++)
        {
            if( scalar_error(
                    original_value[i],
                    rebuilt_value[i]) >
                TEST_TOL )
            {
                printf(
                    "FAIL: candidate reconstruction "
                    "%d component %d error %.17e\n",
                    (int)j,
                    (int)i,
                    (double)scalar_error(
                        original_value[i],
                        rebuilt_value[i]));

                return 1;
            }
        }
    }


    /*
     * Zero-rank blocks must return a zero rectangular coupling.
     */
    for(i=0;i<NOMINAL;i++)
    {
        err = set_vector3(
            W[i],
            0.0,
            0.0,
            0.0);

        if( err )
            return 1;
    }


    active_rank = -1;

    err = lis_eigen_block_rank_deflate(
        CAPACITY,
        NOMINAL,
        W,
        Q,
        RANK_TOL,
        &active_rank,
        gram,
        eval,
        evec,
        transform);

    if( err )
    {
        printf(
            "FAIL: zero-rank coupling deflation returned %d\n",
            (int)err);

        return 1;
    }


    if( active_rank!=0 )
    {
        printf(
            "FAIL: zero coupling active_rank=%d expected 0\n",
            (int)active_rank);

        return 1;
    }


    err = lis_eigen_block_rectangular_coupling(
        NOMINAL,
        active_rank,
        eval,
        evec,
        coupling);

    if( err )
    {
        printf(
            "FAIL: zero-rank rectangular coupling returned %d\n",
            (int)err);

        return 1;
    }


    for(i=0;i<NOMINAL*NOMINAL;i++)
    {
        if( scalar_error(
                coupling[i],
                0.0)>TEST_TOL )
        {
            printf(
                "FAIL: zero-rank coupling storage is not zero\n");

            return 1;
        }
    }


    /*
     * Argument guard.
     */
    err = lis_eigen_block_rectangular_coupling(
        NOMINAL,
        NOMINAL+1,
        eval,
        evec,
        coupling);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "FAIL: rectangular coupling rank guard returned %d\n",
            (int)err);

        return 1;
    }


    return 0;
}


static int
check_argument_guards(
    LIS_VECTOR W[],
    LIS_VECTOR Q[])
{
    LIS_SCALAR gram[NOMINAL*NOMINAL];
    LIS_SCALAR eval[NOMINAL];
    LIS_SCALAR evec[NOMINAL*NOMINAL];
    LIS_SCALAR transform[NOMINAL*NOMINAL];

    LIS_INT active_rank;
    LIS_INT err;

    active_rank = -1;

    err = lis_eigen_block_rank_deflate(
        2,
        3,
        W,
        Q,
        RANK_TOL,
        &active_rank,
        gram,
        eval,
        evec,
        transform);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "FAIL: capacity guard returned %d\n",
            (int)err);

        return 1;
    }

    err = lis_eigen_block_rank_deflate(
        CAPACITY,
        NOMINAL,
        W,
        Q,
        -1.0,
        &active_rank,
        gram,
        eval,
        evec,
        transform);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "FAIL: rank tolerance guard returned %d\n",
            (int)err);

        return 1;
    }

    return 0;
}


int
main(
    int argc,
    char **argv)
{
    LIS_MATRIX B;

    LIS_VECTOR W[NOMINAL];
    LIS_VECTOR BW[NOMINAL];

    LIS_VECTOR Q[CAPACITY];
    LIS_VECTOR BQ[CAPACITY];

    LIS_INT err;
    LIS_INT i;

    int failed;

    B = NULL;

    for(i=0;i<NOMINAL;i++)
    {
        W[i] = NULL;
        BW[i] = NULL;
    }

    for(i=0;i<CAPACITY;i++)
    {
        Q[i] = NULL;
        BQ[i] = NULL;
    }

    failed = 0;

    err = lis_initialize(
        &argc,
        &argv);

    if( err )
        return (int)err;

    err = create_vector3(
        &W[0],
        1.0,
        0.0,
        0.0);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    for(i=1;i<NOMINAL;i++)
    {
        err = lis_vector_duplicate(
            W[0],
            &W[i]);

        if( err )
        {
            failed = 1;
            goto cleanup;
        }
    }

    for(i=0;i<NOMINAL;i++)
    {
        err = lis_vector_duplicate(
            W[0],
            &BW[i]);

        if( err )
        {
            failed = 1;
            goto cleanup;
        }
    }

    for(i=0;i<CAPACITY;i++)
    {
        err = lis_vector_duplicate(
            W[0],
            &Q[i]);

        if( err )
        {
            failed = 1;
            goto cleanup;
        }

        err = lis_vector_duplicate(
            W[0],
            &BQ[i]);

        if( err )
        {
            failed = 1;
            goto cleanup;
        }
    }

    err = create_diagonal_B(&B);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    if( !failed )
    {
        failed =
            check_linear_combination(
                W,
                Q);
    }

    if( !failed )
    {
        failed =
            check_standard_rank(
                W,
                Q);
    }

    if( !failed )
    {
        failed =
            check_threshold_rank(
                W,
                Q);
    }

    if( !failed )
    {
        failed =
            check_zero_rank(
                W,
                Q);
    }

    if( !failed )
    {
        failed =
            check_b_rank(
                B,
                W,
                BW,
                Q,
                BQ);
    }

    if( !failed )
    {
        failed =
            check_ritz_reconstruction(
                W,
                Q);
    }

    if( !failed )
    {
        failed =
            check_rectangular_coupling(
                W,
                Q,
                BQ);
    }

    if( !failed )
    {
        failed =
            check_argument_guards(
                W,
                Q);
    }


cleanup:

    if( B )
        lis_matrix_destroy(B);

    for(i=0;i<NOMINAL;i++)
    {
        if( W[i] )
            lis_vector_destroy(W[i]);

        if( BW[i] )
            lis_vector_destroy(BW[i]);
    }

    for(i=0;i<CAPACITY;i++)
    {
        if( Q[i] )
            lis_vector_destroy(Q[i]);

        if( BQ[i] )
            lis_vector_destroy(BQ[i]);
    }

    lis_finalize();

    if( failed )
        return 1;

    printf(
        "PASS: common block rank-revealing deflation\n");

    return 0;
}
