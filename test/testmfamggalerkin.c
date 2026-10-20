#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <math.h>
#ifdef _COMPLEX
#include <complex.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lislib.h"

#define TEST_N 8
#define MATCH_TOL ((LIS_REAL)1.0e-11)

/*
 * Internal MFAMG regression hook.
 * This is intentionally not part of the public LIS API.
 */
extern LIS_INT
lis_precon_mfamg_build_serial_coarse_for_test(
    LIS_MATRIX A,
    LIS_INT nc,
    LIS_INT aggregate_count,
    LIS_INT *aggregate,
    LIS_INT *aggregate_rank,
    LIS_INT *aggregate_offset,
    LIS_INT *prolong_row_ptr,
    LIS_SCALAR *prolong,
    LIS_INT multivector,
    LIS_MATRIX *coarse_A);


static LIS_SCALAR
galerkin_product(
    LIS_SCALAR pi,
    LIS_SCALAR aij,
    LIS_SCALAR pj)
{
#ifdef _COMPLEX
    return conj(pi) * aij * pj;
#else
    return pi * aij * pj;
#endif
}


static LIS_REAL
scalar_abs(
    LIS_SCALAR value)
{
#ifdef _COMPLEX
    return (LIS_REAL)cabs(value);
#else
    return (LIS_REAL)fabs(value);
#endif
}


static int
scalar_matches(
    LIS_SCALAR a,
    LIS_SCALAR b)
{
    LIS_REAL diff;
    LIS_REAL scale;

    diff =
        scalar_abs(
            a-b);

    scale =
        (LIS_REAL)1.0
        +
        scalar_abs(a)
        +
        scalar_abs(b);

    return
        diff
        <=
        MATCH_TOL
        *
        scale;
}


static LIS_INT
create_test_csr(
    LIS_MATRIX *A)
{
    LIS_INT i;
    LIS_INT err;

    if(A==NULL)
        return LIS_ERR_ILL_ARG;

    *A = NULL;

    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            A);

    if(err)
        return err;

    err =
        lis_matrix_set_size(
            *A,
            0,
            TEST_N);

    if(err)
        goto fail;

    /*
     * Deliberately nonsymmetric sparse matrix.
     * The test is for P^H A P, not only for an SPD special case.
     */
    for(i=0;i<TEST_N;i++)
    {
        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                (LIS_SCALAR)
                (
                    4.0
                    +
                    0.25*(double)i
                ),
                *A);

        if(err)
            goto fail;

        if(i+1<TEST_N)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i+1,
                    (LIS_SCALAR)
                    (
                        -1.0
                        -
                        0.05*(double)i
                    ),
                    *A);

            if(err)
                goto fail;
        }

        if(i>0)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i-1,
                    (LIS_SCALAR)
                    (
                        -0.55
                        -
                        0.03*(double)i
                    ),
                    *A);

            if(err)
                goto fail;
        }

        if(i+2<TEST_N && (i%2)==0)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i+2,
                    (LIS_SCALAR)
                    (
                        0.15
                        +
                        0.02*(double)i
                    ),
                    *A);

            if(err)
                goto fail;
        }
    }

    err =
        lis_matrix_assemble(
            *A);

    if(err)
        goto fail;

    if(
        (*A)->nprocs!=1
        ||
        (*A)->matrix_type!=LIS_MATRIX_CSR
        ||
        (*A)->n!=TEST_N
        ||
        (*A)->gn!=TEST_N)
    {
        err =
            LIS_ERR_ILL_ARG;

        goto fail;
    }

    return LIS_SUCCESS;

fail:

    if(*A)
    {
        lis_matrix_destroy(
            *A);

        *A = NULL;
    }

    return err;
}


static LIS_INT
build_reference(
    LIS_MATRIX A,
    LIS_INT nc,
    LIS_INT aggregate_count,
    const LIS_INT *aggregate,
    const LIS_INT *aggregate_rank,
    const LIS_INT *aggregate_offset,
    const LIS_INT *prolong_row_ptr,
    const LIS_SCALAR *prolong,
    LIS_INT multivector,
    LIS_SCALAR *dense)
{
    LIS_INT i;
    LIS_INT p;
    LIS_INT j;
    LIS_INT ai;
    LIS_INT aj;
    LIS_INT qi;
    LIS_INT qj;
    LIS_INT ri;
    LIS_INT rj;
    LIS_INT ci;
    LIS_INT cj;
    LIS_INT row_i;
    LIS_INT row_j;

    LIS_SCALAR pi;
    LIS_SCALAR pj;

    if(
        A==NULL
        ||
        dense==NULL
        ||
        aggregate==NULL
        ||
        prolong==NULL
        ||
        nc<=0
        ||
        A->nprocs!=1
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->n<=0)
    {
        return LIS_ERR_ILL_ARG;
    }

    if(
        multivector
        &&
        (
            aggregate_count<=0
            ||
            aggregate_rank==NULL
            ||
            aggregate_offset==NULL
            ||
            prolong_row_ptr==NULL
        ))
    {
        return LIS_ERR_ILL_ARG;
    }

    memset(
        dense,
        0,
        (size_t)nc
        *
        (size_t)nc
        *
        sizeof(LIS_SCALAR));

    for(i=0;i<A->n;i++)
    {
        ai =
            aggregate[i];

        if(
            ai<0
            ||
            (
                multivector
                ?
                ai>=aggregate_count
                :
                ai>=nc
            ))
        {
            return LIS_ERR_ILL_ARG;
        }

        if(multivector)
        {
            ri =
                aggregate_rank[ai];

            row_i =
                prolong_row_ptr[i];

            if(
                ri<0
                ||
                row_i<0
                ||
                prolong_row_ptr[i+1]-row_i!=ri)
            {
                return LIS_ERR_ILL_ARG;
            }
        }
        else
        {
            ri = 1;
            row_i = i;
        }

        for(
            p=A->ptr[i];
            p<A->ptr[i+1];
            p++)
        {
            j =
                A->index[p];

            if(j<0 || j>=A->n)
                return LIS_ERR_ILL_ARG;

            aj =
                aggregate[j];

            if(
                aj<0
                ||
                (
                    multivector
                    ?
                    aj>=aggregate_count
                    :
                    aj>=nc
                ))
            {
                return LIS_ERR_ILL_ARG;
            }

            if(multivector)
            {
                rj =
                    aggregate_rank[aj];

                row_j =
                    prolong_row_ptr[j];

                if(
                    rj<0
                    ||
                    row_j<0
                    ||
                    prolong_row_ptr[j+1]-row_j!=rj)
                {
                    return LIS_ERR_ILL_ARG;
                }
            }
            else
            {
                rj = 1;
                row_j = j;
            }

            for(qi=0;qi<ri;qi++)
            {
                if(multivector)
                {
                    ci =
                        aggregate_offset[ai]
                        +
                        qi;

                    pi =
                        prolong[
                            row_i
                            +
                            qi];
                }
                else
                {
                    ci = ai;
                    pi = prolong[i];
                }

                if(ci<0 || ci>=nc)
                    return LIS_ERR_ILL_ARG;

                for(qj=0;qj<rj;qj++)
                {
                    if(multivector)
                    {
                        cj =
                            aggregate_offset[aj]
                            +
                            qj;

                        pj =
                            prolong[
                                row_j
                                +
                                qj];
                    }
                    else
                    {
                        cj = aj;
                        pj = prolong[j];
                    }

                    if(cj<0 || cj>=nc)
                        return LIS_ERR_ILL_ARG;

                    dense[
                        (size_t)ci
                        *
                        (size_t)nc
                        +
                        (size_t)cj]
                        +=
                        galerkin_product(
                            pi,
                            A->value[p],
                            pj);
                }
            }
        }
    }

    return LIS_SUCCESS;
}


static int
compare_coarse_matrix(
    const char *name,
    LIS_MATRIX A,
    LIS_INT nc,
    const LIS_SCALAR *reference)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT p;
    LIS_INT expected_nnz;
    LIS_INT actual_nnz;
    LIS_INT previous_col;

    unsigned char *seen;

    int failed;

    if(
        name==NULL
        ||
        A==NULL
        ||
        reference==NULL
        ||
        nc<=0)
    {
        return 1;
    }

    if(
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs!=1
        ||
        A->n!=nc
        ||
        A->gn!=nc
        ||
        A->ptr==NULL
        ||
        A->index==NULL
        ||
        A->value==NULL)
    {
        fprintf(
            stderr,
            "%s: invalid coarse CSR layout\n",
            name);

        return 1;
    }

    seen =
        (unsigned char *)calloc(
            (size_t)nc*(size_t)nc,
            sizeof(unsigned char));

    if(seen==NULL)
    {
        fprintf(
            stderr,
            "%s: comparison allocation failed\n",
            name);

        return 1;
    }

    failed = 0;
    expected_nnz = 0;

    for(i=0;i<nc;i++)
    {
        for(j=0;j<nc;j++)
        {
            if(
                reference[
                    (size_t)i
                    *
                    (size_t)nc
                    +
                    (size_t)j]
                !=
                (LIS_SCALAR)0.0)
            {
                expected_nnz++;
            }
        }
    }

    actual_nnz =
        A->ptr[nc];

    if(actual_nnz!=expected_nnz)
    {
        fprintf(
            stderr,
            "%s: nnz mismatch: got=%d expected=%d\n",
            name,
            (int)actual_nnz,
            (int)expected_nnz);

        failed = 1;
    }

    for(i=0;i<nc;i++)
    {
        previous_col = -1;

        if(
            A->ptr[i]<0
            ||
            A->ptr[i+1]<A->ptr[i]
            ||
            A->ptr[i+1]>actual_nnz)
        {
            fprintf(
                stderr,
                "%s: invalid CSR row pointer at row %d\n",
                name,
                (int)i);

            failed = 1;
            continue;
        }

        for(
            p=A->ptr[i];
            p<A->ptr[i+1];
            p++)
        {
            j =
                A->index[p];

            if(j<0 || j>=nc)
            {
                fprintf(
                    stderr,
                    "%s: invalid column %d at row %d\n",
                    name,
                    (int)j,
                    (int)i);

                failed = 1;
                continue;
            }

            if(j<=previous_col)
            {
                fprintf(
                    stderr,
                    "%s: coarse CSR columns are not strictly increasing at row %d\n",
                    name,
                    (int)i);

                failed = 1;
            }

            previous_col = j;

            if(
                reference[
                    (size_t)i
                    *
                    (size_t)nc
                    +
                    (size_t)j]
                ==
                (LIS_SCALAR)0.0)
            {
                fprintf(
                    stderr,
                    "%s: unexpected explicit zero/non-reference entry at (%d,%d)\n",
                    name,
                    (int)i,
                    (int)j);

                failed = 1;
            }
            else if(
                !scalar_matches(
                    A->value[p],
                    reference[
                        (size_t)i
                        *
                        (size_t)nc
                        +
                        (size_t)j]))
            {
                fprintf(
                    stderr,
                    "%s: value mismatch at (%d,%d): diff=%e\n",
                    name,
                    (int)i,
                    (int)j,
                    (double)scalar_abs(
                        A->value[p]
                        -
                        reference[
                            (size_t)i
                            *
                            (size_t)nc
                            +
                            (size_t)j]));

                failed = 1;
            }

            seen[
                (size_t)i
                *
                (size_t)nc
                +
                (size_t)j]
                = 1;
        }
    }

    for(i=0;i<nc;i++)
    {
        for(j=0;j<nc;j++)
        {
            if(
                reference[
                    (size_t)i
                    *
                    (size_t)nc
                    +
                    (size_t)j]
                !=
                (LIS_SCALAR)0.0
                &&
                !seen[
                    (size_t)i
                    *
                    (size_t)nc
                    +
                    (size_t)j])
            {
                fprintf(
                    stderr,
                    "%s: missing coarse entry at (%d,%d)\n",
                    name,
                    (int)i,
                    (int)j);

                failed = 1;
            }
        }
    }

    free(
        seen);

    return failed;
}


static int
check_scalar_galerkin(
    LIS_MATRIX A)
{
    LIS_MATRIX Ac;

    LIS_INT aggregate[TEST_N] =
    {
        0,0,
        1,1,
        2,2,
        3,3
    };

    LIS_SCALAR prolong[TEST_N] =
    {
        (LIS_SCALAR)1.00,
        (LIS_SCALAR)0.75,
        (LIS_SCALAR)1.10,
        (LIS_SCALAR)0.90,
        (LIS_SCALAR)0.80,
        (LIS_SCALAR)1.20,
        (LIS_SCALAR)0.95,
        (LIS_SCALAR)1.05
    };

    LIS_SCALAR reference[4*4];

    LIS_INT err;

    int failed;

    Ac = NULL;
    failed = 0;

    err =
        build_reference(
            A,
            4,
            4,
            aggregate,
            NULL,
            NULL,
            NULL,
            prolong,
            LIS_FALSE,
            reference);

    if(err)
    {
        fprintf(
            stderr,
            "scalar reference build failed rc=%d\n",
            (int)err);

        return 1;
    }

    err =
        lis_precon_mfamg_build_serial_coarse_for_test(
            A,
            4,
            4,
            aggregate,
            NULL,
            NULL,
            NULL,
            prolong,
            LIS_FALSE,
            &Ac);

    if(err || Ac==NULL)
    {
        fprintf(
            stderr,
            "scalar sparse Galerkin build failed rc=%d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    failed =
        compare_coarse_matrix(
            "scalar Galerkin",
            Ac,
            4,
            reference);

cleanup:

    if(Ac)
    {
        lis_matrix_destroy(
            Ac);
    }

    if(!failed)
    {
        printf(
            "MFAMG SPARSE GALERKIN SCALAR TEST PASSED\n");
    }

    return failed;
}


static LIS_INT
build_multivector_prolongation(
    const LIS_INT *aggregate,
    const LIS_INT *aggregate_rank,
    LIS_INT *row_ptr,
    LIS_SCALAR *prolong,
    LIS_INT capacity,
    LIS_INT *used)
{
    LIS_INT i;
    LIS_INT q;
    LIS_INT a;
    LIS_INT rank;
    LIS_INT pos;

    if(
        aggregate==NULL
        ||
        aggregate_rank==NULL
        ||
        row_ptr==NULL
        ||
        prolong==NULL
        ||
        used==NULL
        ||
        capacity<=0)
    {
        return LIS_ERR_ILL_ARG;
    }

    pos = 0;
    row_ptr[0] = 0;

    for(i=0;i<TEST_N;i++)
    {
        a =
            aggregate[i];

        if(a<0 || a>=4)
            return LIS_ERR_ILL_ARG;

        rank =
            aggregate_rank[a];

        if(rank<0)
            return LIS_ERR_ILL_ARG;

        if(pos>capacity-rank)
            return LIS_OUT_OF_MEMORY;

        for(q=0;q<rank;q++)
        {
            prolong[pos++] =
                (LIS_SCALAR)
                (
                    0.30
                    +
                    0.07*(double)(i+1)
                    +
                    0.11*(double)(q+1)
                );
        }

        row_ptr[i+1] =
            pos;
    }

    *used =
        pos;

    return LIS_SUCCESS;
}


static int
check_multivector_galerkin(
    LIS_MATRIX A)
{
    LIS_MATRIX Ac;

    LIS_INT aggregate[TEST_N] =
    {
        0,0,
        1,1,
        2,2,
        3,3
    };

    LIS_INT aggregate_rank[4] =
    {
        1,
        2,
        3,
        2
    };

    LIS_INT aggregate_offset[5] =
    {
        0,
        1,
        3,
        6,
        8
    };

    LIS_INT row_ptr[TEST_N+1];

    LIS_SCALAR prolong[32];
    LIS_SCALAR reference[8*8];

    LIS_INT used;
    LIS_INT err;

    int failed;

    Ac = NULL;
    used = 0;
    failed = 0;

    err =
        build_multivector_prolongation(
            aggregate,
            aggregate_rank,
            row_ptr,
            prolong,
            32,
            &used);

    if(err)
    {
        fprintf(
            stderr,
            "multivector prolongation setup failed rc=%d\n",
            (int)err);

        return 1;
    }

    if(used<=0)
    {
        fprintf(
            stderr,
            "multivector prolongation setup produced no entries\n");

        return 1;
    }

    err =
        build_reference(
            A,
            8,
            4,
            aggregate,
            aggregate_rank,
            aggregate_offset,
            row_ptr,
            prolong,
            LIS_TRUE,
            reference);

    if(err)
    {
        fprintf(
            stderr,
            "multivector reference build failed rc=%d\n",
            (int)err);

        return 1;
    }

    err =
        lis_precon_mfamg_build_serial_coarse_for_test(
            A,
            8,
            4,
            aggregate,
            aggregate_rank,
            aggregate_offset,
            row_ptr,
            prolong,
            LIS_TRUE,
            &Ac);

    if(err || Ac==NULL)
    {
        fprintf(
            stderr,
            "multivector sparse Galerkin build failed rc=%d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    failed =
        compare_coarse_matrix(
            "multivector Galerkin",
            Ac,
            8,
            reference);

cleanup:

    if(Ac)
    {
        lis_matrix_destroy(
            Ac);
    }

    if(!failed)
    {
        printf(
            "MFAMG SPARSE GALERKIN MULTIVECTOR TEST PASSED\n");
    }

    return failed;
}


static int
check_zero_rank_multivector_galerkin(
    LIS_MATRIX A)
{
    LIS_MATRIX Ac;

    LIS_INT aggregate[TEST_N] =
    {
        0,0,
        1,1,
        2,2,
        3,3
    };

    LIS_INT aggregate_rank[4] =
    {
        1,
        2,
        0,
        2
    };

    LIS_INT aggregate_offset[5] =
    {
        0,
        1,
        3,
        3,
        5
    };

    LIS_INT row_ptr[TEST_N+1];
    LIS_SCALAR prolong[24];
    LIS_SCALAR reference[5*5];

    LIS_INT used;
    LIS_INT err;
    int failed;

    Ac = NULL;
    used = 0;
    failed = 0;

    err =
        build_multivector_prolongation(
            aggregate,
            aggregate_rank,
            row_ptr,
            prolong,
            24,
            &used);

    if(err)
    {
        fprintf(stderr,"zero-rank prolongation setup failed rc=%d\n",(int)err);
        return 1;
    }

    err =
        build_reference(
            A,5,4,aggregate,aggregate_rank,aggregate_offset,row_ptr,
            prolong,LIS_TRUE,reference);

    if(err)
    {
        fprintf(stderr,"zero-rank reference build failed rc=%d\n",(int)err);
        return 1;
    }

    err =
        lis_precon_mfamg_build_serial_coarse_for_test(
            A,5,4,aggregate,aggregate_rank,aggregate_offset,row_ptr,
            prolong,LIS_TRUE,&Ac);

    if(err || Ac==NULL)
    {
        fprintf(stderr,"zero-rank sparse Galerkin build failed rc=%d\n",(int)err);
        failed = 1;
        goto cleanup;
    }

    failed =
        compare_coarse_matrix(
            "zero-rank multivector Galerkin",Ac,5,reference);

cleanup:

    if(Ac)
        lis_matrix_destroy(Ac);

    if(!failed)
        printf("MFAMG SPARSE GALERKIN ZERO-RANK TEST PASSED\n");

    return failed;
}


static int
check_invalid_multivector_metadata(
    LIS_MATRIX A)
{
    LIS_MATRIX Ac;

    LIS_INT aggregate[TEST_N] =
    {
        0,0,
        1,1,
        2,2,
        3,3
    };

    LIS_INT aggregate_rank[4] =
    {
        1,
        2,
        3,
        2
    };

    /*
     * Deliberately inconsistent:
     * offset[2]-offset[1] == 3, but rank[1] == 2.
     */
    LIS_INT aggregate_offset[5] =
    {
        0,
        1,
        4,
        7,
        9
    };

    LIS_INT row_ptr[TEST_N+1];
    LIS_SCALAR prolong[32];

    LIS_INT used;
    LIS_INT err;

    Ac = NULL;
    used = 0;

    err =
        build_multivector_prolongation(
            aggregate,
            aggregate_rank,
            row_ptr,
            prolong,
            32,
            &used);

    if(err)
    {
        fprintf(
            stderr,
            "invalid metadata test setup failed rc=%d\n",
            (int)err);

        return 1;
    }

    err =
        lis_precon_mfamg_build_serial_coarse_for_test(
            A,
            9,
            4,
            aggregate,
            aggregate_rank,
            aggregate_offset,
            row_ptr,
            prolong,
            LIS_TRUE,
            &Ac);

    if(Ac)
    {
        lis_matrix_destroy(
            Ac);

        Ac = NULL;
    }

    if(err==LIS_SUCCESS)
    {
        fprintf(
            stderr,
            "invalid multivector metadata was accepted\n");

        return 1;
    }

    printf(
        "MFAMG SPARSE GALERKIN INVALID METADATA TEST PASSED\n");

    return 0;
}


static int
check_exact_cancellation(void)
{
    LIS_MATRIX A;
    LIS_MATRIX Ac;
    LIS_INT aggregate[2] = {0,0};
    LIS_SCALAR prolong[2] = {(LIS_SCALAR)1.0,(LIS_SCALAR)-1.0};
    LIS_SCALAR reference[1];
    LIS_INT i;
    LIS_INT j;
    LIS_INT err;
    int failed;

    A = NULL;
    Ac = NULL;
    err = LIS_SUCCESS;
    failed = 0;

    err = lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto cleanup;
    err = lis_matrix_set_size(A,0,2);
    if(err) goto cleanup;

    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,i,j,(LIS_SCALAR)1.0,A);
            if(err) goto cleanup;
        }
    }

    err = lis_matrix_assemble(A);
    if(err) goto cleanup;

    err = build_reference(
        A,1,1,aggregate,NULL,NULL,NULL,prolong,LIS_FALSE,reference);
    if(err) goto cleanup;

    if(reference[0]!=(LIS_SCALAR)0.0)
    {
        fprintf(stderr,"cancellation reference is not exactly zero\n");
        failed = 1;
        goto cleanup;
    }

    err = lis_precon_mfamg_build_serial_coarse_for_test(
        A,1,1,aggregate,NULL,NULL,NULL,prolong,LIS_FALSE,&Ac);

    /*
     * Preserve the previous dense MFAMG behavior: if exact
     * cancellation removes every coarse entry, construction
     * reports LIS_BREAKDOWN rather than creating an empty CSR.
     */
    if(err!=LIS_BREAKDOWN || Ac!=NULL)
    {
        fprintf(
            stderr,
            "exact cancellation expected LIS_BREAKDOWN, got rc=%d Ac=%p\n",
            (int)err,
            (void *)Ac);

        failed = 1;
        goto cleanup;
    }

    /*
     * LIS_BREAKDOWN is the expected result for this regression.
     * Clear it so the common cleanup path reports success.
     */
    err = LIS_SUCCESS;

cleanup:

    if(Ac) lis_matrix_destroy(Ac);
    if(A) lis_matrix_destroy(A);

    if(err && !failed)
    {
        fprintf(stderr,"cancellation test failed rc=%d\n",(int)err);
        failed = 1;
    }

    if(!failed)
        printf("MFAMG SPARSE GALERKIN CANCELLATION TEST PASSED\n");

    return failed;
}


static int
run_serial(void)
{
    LIS_MATRIX A;

    LIS_INT err;

    int failed;

    A = NULL;
    failed = 0;

    err =
        create_test_csr(
            &A);

    if(err)
    {
        fprintf(
            stderr,
            "MFAMG sparse Galerkin matrix setup failed rc=%d\n",
            (int)err);

        return 1;
    }

    failed |=
        check_scalar_galerkin(
            A);

    failed |=
        check_multivector_galerkin(
            A);

    failed |=
        check_zero_rank_multivector_galerkin(
            A);

    failed |=
        check_invalid_multivector_metadata(
            A);

    failed |=
        check_exact_cancellation();

    lis_matrix_destroy(
        A);

    if(failed)
    {
        fprintf(
            stderr,
            "MFAMG SPARSE GALERKIN TEST FAILED\n");
    }
    else
    {
        printf(
            "MFAMG SPARSE GALERKIN TEST PASSED\n");
    }

    return failed;
}


int
main(
    int argc,
    char *argv[])
{
    LIS_INT err;

    int rank;
    int nprocs;
    int failed;

    setvbuf(
        stdout,
        NULL,
        _IONBF,
        0);

    setvbuf(
        stderr,
        NULL,
        _IONBF,
        0);

    err =
        lis_initialize(
            &argc,
            &argv);

    if(err)
        return 1;

#ifdef USE_MPI
    MPI_Comm_rank(
        LIS_COMM_WORLD,
        &rank);

    MPI_Comm_size(
        LIS_COMM_WORLD,
        &nprocs);
#else
    rank = 0;
    nprocs = 1;
#endif

    failed = 0;

    if(nprocs==1)
    {
        failed =
            run_serial();
    }
    else if(rank==0)
    {
        /*
         * The production helper under test is intentionally serial-only.
         * Distributed Galerkin behavior is covered by testmfamg.
         */
        printf(
            "MFAMG SPARSE GALERKIN TEST SKIPPED: serial-only regression\n");
    }

    lis_finalize();

    return
        failed
        ?
        1
        :
        0;
}
