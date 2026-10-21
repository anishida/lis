/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Serial two-level matrix-free aggregation preconditioner.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "lislib.h"


#define LIS_MFAMG_MAX_LEVELS 32
#define LIS_MFAMG_COARSE_LIMIT 8
#define LIS_MFAMG_STRENGTH_THETA ((LIS_REAL)0.05)
#define LIS_MFAMG_DEFAULT_OMEGA ((LIS_REAL)(2.0/3.0))
#define LIS_MFAMG_RANK_TOL ((LIS_REAL)1.0e-12)
#define LIS_MFAMG_COARSE_PIVOT_TOL LIS_MFAMG_RANK_TOL

#define LIS_MFAMG_SMOOTHER_WEIGHTED_JACOBI 0
#define LIS_MFAMG_SMOOTHER_BLOCK_JACOBI 1
#define LIS_MFAMG_BLOCK_MAX_SIZE 16


struct LIS_MFAMG_LEVEL
{
    LIS_MATRIX A;
    LIS_INT owns_A;

    LIS_INT n;
    LIS_INT nc;

    LIS_INT global_nc;
    LIS_INT coarse_is;
    LIS_INT coarse_ie;

    LIS_INT *coarse_ranges;
    LIS_INT *global_aggregate;
    LIS_INT *global_aggregate_rank;
    LIS_SCALAR *global_prolong;
    LIS_INT global_prolong_width;

    LIS_INT aggregate_count;

    LIS_INT *aggregate;
    LIS_INT *aggregate_rank;
    LIS_INT *aggregate_offset;
    LIS_INT *prolong_row_ptr;

    LIS_SCALAR *prolong;

    LIS_INT smoother_type;

    LIS_VECTOR diagonal;
    LIS_VECTOR ax;

    LIS_REAL omega;

    /*
     * Fixed-size local block-Jacobi state.
     *
     * The final local block may be shorter than block_size.
     * Singular numerical blocks fall back to scalar Jacobi.
     */
    LIS_INT block_size;
    LIS_INT block_count;
    LIS_SCALAR *block_lu;
    LIS_INT *block_pivots;
    LIS_INT *block_valid;
    LIS_SCALAR *block_work;
    LIS_REAL *block_row_scale;

    LIS_VECTOR rhs;
    LIS_VECTOR x;
};


struct LIS_MFAMG_STRUCT
{
    LIS_INT nlevels;
    LIS_INT level_capacity;
    LIS_INT smoother_block_size;
    struct LIS_MFAMG_LEVEL *levels;

    LIS_INT coarse_n;

    LIS_VECTOR alias_rhs;

    LIS_INT *refresh_cols;
    LIS_SCALAR *refresh_vals;
    LIS_INT refresh_row_capacity;

    LIS_SCALAR *coarse;
    LIS_SCALAR *lu;
    LIS_INT *pivots;

    LIS_SCALAR *coarse_rhs;
    LIS_SCALAR *coarse_x;
};


static void
lis_mfamg_smoother_destroy(
    struct LIS_MFAMG_LEVEL *level)
{
    if(level==NULL)
        return;


    if(level->diagonal)
    {
        lis_vector_destroy(
            level->diagonal);

        level->diagonal =
            NULL;
    }


    if(level->ax)
    {
        lis_vector_destroy(
            level->ax);

        level->ax =
            NULL;
    }


    if(level->block_lu)
    {
        lis_free(
            level->block_lu);

        level->block_lu =
            NULL;
    }


    if(level->block_pivots)
    {
        lis_free(
            level->block_pivots);

        level->block_pivots =
            NULL;
    }


    if(level->block_valid)
    {
        lis_free(
            level->block_valid);

        level->block_valid =
            NULL;
    }


    if(level->block_work)
    {
        lis_free(
            level->block_work);

        level->block_work =
            NULL;
    }


    if(level->block_row_scale)
    {
        lis_free(
            level->block_row_scale);

        level->block_row_scale =
            NULL;
    }


    level->block_size =
        0;

    level->block_count =
        0;

    level->omega =
        0.0;
}



static void
lis_mfamg_state_free(
    struct LIS_MFAMG_STRUCT *state)
{
    LIS_INT level_id;

    struct LIS_MFAMG_LEVEL *level;


    if(state==NULL)
        return;


    if(state->levels)
    {
        for(
            level_id=0;
            level_id<state->nlevels;
            level_id++)
        {
            level =
                &state->levels[level_id];

            if(level->aggregate)
                lis_free(level->aggregate);

            if(level->aggregate_rank)
                lis_free(level->aggregate_rank);

            if(level->aggregate_offset)
                lis_free(level->aggregate_offset);

            if(level->prolong_row_ptr)
                lis_free(level->prolong_row_ptr);

            if(level->prolong)
                lis_free(level->prolong);

            if(level->coarse_ranges)
                lis_free(level->coarse_ranges);

            if(level->global_aggregate)
                lis_free(level->global_aggregate);

            if(level->global_aggregate_rank)
                lis_free(level->global_aggregate_rank);

            if(level->global_prolong)
                lis_free(level->global_prolong);

            lis_mfamg_smoother_destroy(
                level);

            if(level->rhs)
                lis_vector_destroy(level->rhs);

            if(level->x)
                lis_vector_destroy(level->x);

            if(
                level->owns_A
                &&
                level->A)
            {
                lis_matrix_destroy(level->A);
            }
        }

        lis_free(state->levels);
    }


    if(state->alias_rhs)
        lis_vector_destroy(state->alias_rhs);


    if(state->refresh_cols)
        lis_free(state->refresh_cols);

    if(state->refresh_vals)
        lis_free(state->refresh_vals);


    if(state->coarse)
        lis_free(state->coarse);

    if(state->lu)
        lis_free(state->lu);

    if(state->pivots)
        lis_free(state->pivots);

    if(state->coarse_rhs)
        lis_free(state->coarse_rhs);

    if(state->coarse_x)
        lis_free(state->coarse_x);


    lis_free(state);
}


static LIS_INT
lis_mfamg_row_query(
    LIS_MATRIX A,
    LIS_INT row,
    LIS_INT *nnz)
{
    if(
        A==NULL
        ||
        nnz==NULL
        ||
        row<0
        ||
        row>=A->n)
    {
        return LIS_ERR_ILL_ARG;
    }

    if(A->matrix_type==LIS_MATRIX_CSR)
    {
        *nnz =
            A->ptr[row+1]
            -
            A->ptr[row];

        return LIS_SUCCESS;
    }

    if(A->matrix_type==LIS_MATRIX_USER)
    {
        return
            lis_matrix_get_user_row(
                A,
                row,
                nnz,
                NULL,
                NULL);
    }

    return LIS_ERR_NOT_IMPLEMENTED;
}


static LIS_INT
lis_mfamg_row_fill(
    LIS_MATRIX A,
    LIS_INT row,
    LIS_INT capacity,
    LIS_INT *nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    LIS_INT k;
    LIS_INT p;
    LIS_INT count;
    LIS_INT local_col;


    if(
        A==NULL
        ||
        nnz==NULL
        ||
        cols==NULL
        ||
        vals==NULL
        ||
        capacity<0
        ||
        row<0
        ||
        row>=A->n)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(A->matrix_type==LIS_MATRIX_CSR)
    {
        count =
            A->ptr[row+1]
            -
            A->ptr[row];


        if(count>capacity)
        {
            *nnz = count;

            return LIS_ERR_ILL_ARG;
        }


        k = 0;


        for(
            p=A->ptr[row];
            p<A->ptr[row+1];
            p++)
        {
            local_col =
                A->index[p];


            if(
                local_col<0
                ||
                local_col>=A->np)
            {
                return LIS_ERR_ILL_ARG;
            }


            if(local_col<A->n)
            {
                cols[k] =
                    A->is
                    +
                    local_col;
            }
            else
            {
                if(A->l2g_map==NULL)
                {
                    return LIS_ERR_ILL_ARG;
                }


                cols[k] =
                    A->l2g_map[
                        local_col
                        -
                        A->n];
            }


            if(
                cols[k]<0
                ||
                cols[k]>=A->gn)
            {
                return LIS_ERR_ILL_ARG;
            }


            vals[k] =
                A->value[p];


            k++;
        }


        *nnz = count;


        return LIS_SUCCESS;
    }


    if(A->matrix_type==LIS_MATRIX_USER)
    {
        *nnz = capacity;


        return
            lis_matrix_get_user_row(
                A,
                row,
                nnz,
                cols,
                vals);
    }


    return LIS_ERR_NOT_IMPLEMENTED;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_estimate_smoother_omega"
static LIS_INT
lis_mfamg_estimate_smoother_omega(
    LIS_MATRIX A,
    LIS_VECTOR diagonal,
    LIS_REAL *omega)
{
    LIS_INT i;
    LIS_INT k;
    LIS_INT nnz;
    LIS_INT row_nnz;
    LIS_INT max_row_nnz;
    LIS_INT err;

    LIS_INT *cols;
    LIS_SCALAR *vals;

    LIS_REAL mag;
    LIS_REAL diagonal_mag;
    LIS_REAL row_sum;
    LIS_REAL ratio;
    LIS_REAL local_bound;
    LIS_REAL global_bound;
    LIS_REAL candidate;


    if(
        A==NULL
        ||
        diagonal==NULL
        ||
        omega==NULL
        ||
        A->n<=0
        ||
        diagonal->n!=A->n)
    {
        return LIS_ERR_ILL_ARG;
    }


    cols = NULL;
    vals = NULL;

    max_row_nnz = 0;
    local_bound = 0.0;
    global_bound = 0.0;

    err = LIS_SUCCESS;


    for(i=0;i<A->n;i++)
    {
        row_nnz = 0;

        err =
            lis_mfamg_row_query(
                A,
                i,
                &row_nnz);

        if(err)
            goto cleanup;


        if(row_nnz>max_row_nnz)
        {
            max_row_nnz =
                row_nnz;
        }
    }


    if(max_row_nnz<=0)
    {
        err = LIS_ERR_ILL_ARG;
        goto cleanup;
    }


    cols =
        (LIS_INT *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_INT),
            "lis_mfamg_estimate_smoother_omega::cols");

    vals =
        (LIS_SCALAR *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_SCALAR),
            "lis_mfamg_estimate_smoother_omega::vals");


    if(
        cols==NULL
        ||
        vals==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)max_row_nnz
                *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        err = LIS_OUT_OF_MEMORY;
        goto cleanup;
    }


    for(i=0;i<A->n;i++)
    {
        diagonal_mag =
            fabs(
                diagonal->value[i]);


        if(
            diagonal_mag!=diagonal_mag
            ||
            diagonal_mag>LIS_SCALAR_MAX
            ||
            diagonal_mag==0.0)
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG smoother requires finite nonzero diagonal entries\n");

            err = LIS_ERR_ILL_ARG;
            goto cleanup;
        }


        nnz =
            max_row_nnz;


        err =
            lis_mfamg_row_fill(
                A,
                i,
                max_row_nnz,
                &nnz,
                cols,
                vals);

        if(err)
            goto cleanup;


        row_sum =
            0.0;


        for(k=0;k<nnz;k++)
        {
            mag =
                fabs(
                    vals[k]);


            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "MFAMG smoother encountered a non-finite matrix entry\n");

                err = LIS_ERR_ILL_ARG;
                goto cleanup;
            }


            row_sum +=
                mag;


            if(
                row_sum!=row_sum
                ||
                row_sum>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "MFAMG smoother row bound overflowed\n");

                err = LIS_ERR_ILL_ARG;
                goto cleanup;
            }
        }


        ratio =
            row_sum
            /
            diagonal_mag;


        if(
            ratio!=ratio
            ||
            ratio>LIS_SCALAR_MAX
            ||
            !(ratio>0.0))
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG smoother produced an invalid row bound\n");

            err = LIS_ERR_ILL_ARG;
            goto cleanup;
        }


        if(ratio>local_bound)
        {
            local_bound =
                ratio;
        }
    }


    global_bound =
        local_bound;


#ifdef USE_MPI

    if(A->nprocs>1)
    {
        if(
            MPI_Allreduce(
                &local_bound,
                &global_bound,
                1,
                LIS_MPI_REAL,
                MPI_MAX,
                A->comm)
            !=
            MPI_SUCCESS)
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG smoother damping reduction failed\n");

            err = LIS_ERR_ILL_ARG;
            goto cleanup;
        }
    }

#endif


    if(
        global_bound!=global_bound
        ||
        global_bound>LIS_SCALAR_MAX
        ||
        !(global_bound>0.0))
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG smoother produced an invalid global bound\n");

        err = LIS_ERR_ILL_ARG;
        goto cleanup;
    }


    *omega =
        LIS_MFAMG_DEFAULT_OMEGA;


    candidate =
        (LIS_REAL)(4.0/3.0)
        /
        global_bound;


    if(
        candidate!=candidate
        ||
        candidate>LIS_SCALAR_MAX
        ||
        !(candidate>0.0))
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG smoother produced an invalid damping factor\n");

        err = LIS_ERR_ILL_ARG;
        goto cleanup;
    }


    if(candidate<*omega)
    {
        *omega =
            candidate;
    }


cleanup:

    if(vals)
        lis_free(vals);

    if(cols)
        lis_free(cols);


    return err;
}

static LIS_INT
lis_mfamg_smoother_validate_diagonal(
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_INT i;

    LIS_REAL diagonal_mag;


    if(
        level==NULL
        ||
        level->diagonal==NULL
        ||
        level->n<=0
        ||
        level->diagonal->n!=level->n)
    {
        return LIS_ERR_ILL_ARG;
    }


    for(i=0;i<level->n;i++)
    {
        diagonal_mag =
            fabs(
                level->diagonal->value[i]);

        if(
            diagonal_mag!=diagonal_mag
            ||
            diagonal_mag>LIS_SCALAR_MAX
            ||
            diagonal_mag==0.0)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG requires finite nonzero diagonal entries\n");

            return LIS_BREAKDOWN;
        }
    }


    return LIS_SUCCESS;
}


static LIS_INT
lis_mfamg_block_factor(
    LIS_INT n,
    LIS_INT ld,
    LIS_SCALAR *lu,
    LIS_INT *pivots,
    LIS_REAL *row_scale,
    LIS_INT *usable)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT p;

    LIS_REAL best_ratio;
    LIS_REAL ratio;
    LIS_REAL mag;
    LIS_REAL scale;
    LIS_REAL tmp_scale;

    LIS_SCALAR pivot;
    LIS_SCALAR factor;
    LIS_SCALAR tmp;


    if(
        n<=0
        ||
        ld<n
        ||
        lu==NULL
        ||
        pivots==NULL
        ||
        row_scale==NULL
        ||
        usable==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    *usable =
        LIS_FALSE;


    for(i=0;i<n;i++)
    {
        scale =
            0.0;

        for(j=0;j<n;j++)
        {
            mag =
                fabs(
                    lu[i*ld+j]);

            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG block smoother contains non-finite values\n");

                return LIS_BREAKDOWN;
            }


            if(mag>scale)
            {
                scale =
                    mag;
            }
        }


        /*
         * A completely zero row makes this numerical block
         * unsuitable for dense block inversion.  The caller
         * may still use the scalar fallback for the block.
         */
        if(scale==0.0)
        {
            return LIS_SUCCESS;
        }


        row_scale[i] =
            scale;
    }


    for(k=0;k<n;k++)
    {
        p =
            -1;

        best_ratio =
            -1.0;


        for(i=k;i<n;i++)
        {
            ratio =
                fabs(
                    lu[i*ld+k])
                /
                row_scale[i];

            if(
                ratio!=ratio
                ||
                ratio>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG block pivot search became non-finite\n");

                return LIS_BREAKDOWN;
            }


            if(ratio>best_ratio)
            {
                best_ratio =
                    ratio;

                p =
                    i;
            }
        }


        if(
            p<0
            ||
            best_ratio<=LIS_MFAMG_RANK_TOL)
        {
            return LIS_SUCCESS;
        }


        pivots[k] =
            p;


        if(p!=k)
        {
            for(j=0;j<n;j++)
            {
                tmp =
                    lu[k*ld+j];

                lu[k*ld+j] =
                    lu[p*ld+j];

                lu[p*ld+j] =
                    tmp;
            }


            tmp_scale =
                row_scale[k];

            row_scale[k] =
                row_scale[p];

            row_scale[p] =
                tmp_scale;
        }


        pivot =
            lu[k*ld+k];

        mag =
            fabs(
                pivot);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX
            ||
            mag==0.0)
        {
            return LIS_SUCCESS;
        }


        for(i=k+1;i<n;i++)
        {
            factor =
                lu[i*ld+k]
                /
                pivot;

            mag =
                fabs(
                    factor);

            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG block elimination became non-finite\n");

                return LIS_BREAKDOWN;
            }


            lu[i*ld+k] =
                factor;


            for(j=k+1;j<n;j++)
            {
                lu[i*ld+j] -=
                    factor
                    *
                    lu[k*ld+j];

                mag =
                    fabs(
                        lu[i*ld+j]);

                if(
                    mag!=mag
                    ||
                    mag>LIS_SCALAR_MAX)
                {
                    LIS_SETERR(
                        LIS_BREAKDOWN,
                        "MFAMG block factorization became non-finite\n");

                    return LIS_BREAKDOWN;
                }
            }
        }
    }


    *usable =
        LIS_TRUE;


    return LIS_SUCCESS;
}



static LIS_INT
lis_mfamg_block_solve(
    LIS_INT n,
    LIS_INT ld,
    const LIS_SCALAR *lu,
    const LIS_INT *pivots,
    LIS_SCALAR *x)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT p;

    LIS_REAL mag;

    LIS_SCALAR sum;
    LIS_SCALAR tmp;


    if(
        n<=0
        ||
        ld<n
        ||
        lu==NULL
        ||
        pivots==NULL
        ||
        x==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    for(k=0;k<n;k++)
    {
        p =
            pivots[k];

        if(
            p<0
            ||
            p>=n)
        {
            return LIS_ERR_ILL_ARG;
        }


        if(p!=k)
        {
            tmp =
                x[k];

            x[k] =
                x[p];

            x[p] =
                tmp;
        }
    }


    for(i=1;i<n;i++)
    {
        sum =
            x[i];

        for(j=0;j<i;j++)
        {
            sum -=
                lu[i*ld+j]
                *
                x[j];
        }

        x[i] =
            sum;
    }


    for(i=n-1;i>=0;i--)
    {
        sum =
            x[i];

        for(j=i+1;j<n;j++)
        {
            sum -=
                lu[i*ld+j]
                *
                x[j];
        }


        mag =
            fabs(
                lu[i*ld+i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX
            ||
            mag==0.0)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG block solve encountered an invalid pivot\n");

            return LIS_BREAKDOWN;
        }


        x[i] =
            sum
            /
            lu[i*ld+i];


        mag =
            fabs(
                x[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG block solve became non-finite\n");

            return LIS_BREAKDOWN;
        }
    }


    return LIS_SUCCESS;
}



static LIS_INT
lis_mfamg_block_refresh(
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_INT i;
    LIS_INT k;
    LIS_INT row;
    LIS_INT block;
    LIS_INT start;
    LIS_INT width;
    LIS_INT nnz;
    LIS_INT max_row_nnz;
    LIS_INT usable;
    LIS_INT err;

    LIS_INT *cols;

    LIS_SCALAR *vals;
    LIS_SCALAR *lu;

    LIS_REAL mag;

    LIS_INT global_start;


    if(
        level==NULL
        ||
        level->A==NULL
        ||
        level->diagonal==NULL
        ||
        level->n<=0
        ||
        level->block_size<=1
        ||
        level->block_size>LIS_MFAMG_BLOCK_MAX_SIZE
        ||
        level->block_count<=0
        ||
        level->block_lu==NULL
        ||
        level->block_pivots==NULL
        ||
        level->block_valid==NULL
        ||
        level->block_row_scale==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    err =
        lis_matrix_get_diagonal(
            level->A,
            level->diagonal);

    if(err)
        return err;


    /*
     * Use the historical MFAMG damping as a conservative
     * block-Jacobi damping factor.  Scalar adaptive damping
     * remains unchanged for the default smoother.
     */
    level->omega =
        LIS_MFAMG_DEFAULT_OMEGA;


    max_row_nnz =
        0;


    for(i=0;i<level->n;i++)
    {
        nnz =
            0;

        err =
            lis_mfamg_row_query(
                level->A,
                i,
                &nnz);

        if(err)
            return err;


        if(nnz>max_row_nnz)
        {
            max_row_nnz =
                nnz;
        }
    }


    if(max_row_nnz<=0)
    {
        LIS_SETERR(
            LIS_BREAKDOWN,
            "MFAMG block smoother matrix has no row entries\n");

        return LIS_BREAKDOWN;
    }


    cols =
        (LIS_INT *)lis_malloc(
            (size_t)max_row_nnz
            *
            sizeof(LIS_INT),
            "lis_mfamg_block_refresh::cols");

    vals =
        (LIS_SCALAR *)lis_malloc(
            (size_t)max_row_nnz
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_block_refresh::vals");


    if(
        cols==NULL
        ||
        vals==NULL)
    {
        if(cols)
            lis_free(cols);

        if(vals)
            lis_free(vals);

        LIS_SETERR_MEM(
            (size_t)max_row_nnz
            *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        return LIS_OUT_OF_MEMORY;
    }


    for(
        block=0;
        block<level->block_count;
        block++)
    {
        start =
            block
            *
            level->block_size;

        width =
            level->n
            -
            start;

        if(width>level->block_size)
        {
            width =
                level->block_size;
        }


        lu =
            level->block_lu
            +
            (size_t)block
            *
            (size_t)level->block_size
            *
            (size_t)level->block_size;


        memset(
            lu,
            0,
            (size_t)level->block_size
            *
            (size_t)level->block_size
            *
            sizeof(LIS_SCALAR));


        global_start =
            level->A->is
            +
            start;


        for(i=0;i<width;i++)
        {
            row =
                start
                +
                i;

            nnz =
                max_row_nnz;

            err =
                lis_mfamg_row_fill(
                    level->A,
                    row,
                    max_row_nnz,
                    &nnz,
                    cols,
                    vals);

            if(err)
            {
                lis_free(vals);
                lis_free(cols);

                return err;
            }


            for(k=0;k<nnz;k++)
            {
                if(
                    cols[k]>=global_start
                    &&
                    cols[k]<global_start+width)
                {
                    lu[
                        i*level->block_size
                        +
                        (
                            cols[k]
                            -
                            global_start
                        )]
                        +=
                        vals[k];
                }
            }
        }


        usable =
            LIS_FALSE;

        err =
            lis_mfamg_block_factor(
                width,
                level->block_size,
                lu,
                level->block_pivots
                    +
                    block
                    *
                    level->block_size,
                level->block_row_scale,
                &usable);

        if(err)
        {
            lis_free(vals);
            lis_free(cols);

            return err;
        }


        level->block_valid[block] =
            usable;


        /*
         * A singular or numerically deficient diagonal block
         * is not fatal if all scalar diagonal entries in that
         * block remain valid.  Such a block uses the existing
         * weighted-Jacobi correction.
         */
        if(!usable)
        {
            for(i=0;i<width;i++)
            {
                mag =
                    fabs(
                        level->diagonal->value[
                            start+i]);

                if(
                    mag!=mag
                    ||
                    mag>LIS_SCALAR_MAX
                    ||
                    mag==0.0)
                {
                    lis_free(vals);
                    lis_free(cols);

                    LIS_SETERR(
                        LIS_BREAKDOWN,
                        "MFAMG block fallback requires finite nonzero diagonal entries\n");

                    return LIS_BREAKDOWN;
                }
            }
        }
    }


    lis_free(
        vals);

    lis_free(
        cols);


    return LIS_SUCCESS;
}




static LIS_INT
lis_mfamg_smoother_refresh(
    struct LIS_MFAMG_LEVEL *level);



static LIS_INT
lis_mfamg_smoother_setup(
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_INT err;

    size_t block_entries;
    size_t pivot_entries;


    if(
        level==NULL
        ||
        level->A==NULL
        ||
        level->n<=0
        ||
        level->A->n!=level->n
        ||
        level->diagonal!=NULL
        ||
        level->ax!=NULL
        ||
        level->block_lu!=NULL
        ||
        level->block_pivots!=NULL
        ||
        level->block_valid!=NULL
        ||
        level->block_work!=NULL
        ||
        level->block_row_scale!=NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    switch(level->smoother_type)
    {
    case LIS_MFAMG_SMOOTHER_WEIGHTED_JACOBI:

        err =
            lis_vector_duplicate(
                level->A,
                &level->diagonal);

        if(err)
            goto cleanup;


        err =
            lis_vector_duplicate(
                level->A,
                &level->ax);

        if(err)
            goto cleanup;


        err =
            lis_mfamg_smoother_refresh(
                level);

        if(err)
            goto cleanup;


        return LIS_SUCCESS;


    case LIS_MFAMG_SMOOTHER_BLOCK_JACOBI:

        if(
            level->block_size<=1
            ||
            level->block_size>LIS_MFAMG_BLOCK_MAX_SIZE)
        {
            return LIS_ERR_ILL_ARG;
        }


        level->block_count =
            1
            +
            (
                level->n
                -
                1
            )
            /
            level->block_size;


        block_entries =
            (size_t)level->block_count
            *
            (size_t)level->block_size
            *
            (size_t)level->block_size;

        pivot_entries =
            (size_t)level->block_count
            *
            (size_t)level->block_size;


        if(
            block_entries
            >
            ((size_t)-1)
            /
            sizeof(LIS_SCALAR)
            ||
            pivot_entries
            >
            ((size_t)-1)
            /
            sizeof(LIS_INT))
        {
            LIS_SETERR(
                LIS_OUT_OF_MEMORY,
                "MFAMG block smoother allocation size overflow\n");

            return LIS_OUT_OF_MEMORY;
        }


        err =
            lis_vector_duplicate(
                level->A,
                &level->diagonal);

        if(err)
            goto cleanup;


        err =
            lis_vector_duplicate(
                level->A,
                &level->ax);

        if(err)
            goto cleanup;


        level->block_lu =
            (LIS_SCALAR *)lis_malloc(
                block_entries
                *
                sizeof(LIS_SCALAR),
                "lis_mfamg_smoother_setup::block_lu");

        level->block_pivots =
            (LIS_INT *)lis_malloc(
                pivot_entries
                *
                sizeof(LIS_INT),
                "lis_mfamg_smoother_setup::block_pivots");

        level->block_valid =
            (LIS_INT *)lis_malloc(
                (size_t)level->block_count
                *
                sizeof(LIS_INT),
                "lis_mfamg_smoother_setup::block_valid");

        level->block_work =
            (LIS_SCALAR *)lis_malloc(
                (size_t)level->block_size
                *
                sizeof(LIS_SCALAR),
                "lis_mfamg_smoother_setup::block_work");

        level->block_row_scale =
            (LIS_REAL *)lis_malloc(
                (size_t)level->block_size
                *
                sizeof(LIS_REAL),
                "lis_mfamg_smoother_setup::block_row_scale");


        if(
            level->block_lu==NULL
            ||
            level->block_pivots==NULL
            ||
            level->block_valid==NULL
            ||
            level->block_work==NULL
            ||
            level->block_row_scale==NULL)
        {
            LIS_SETERR_MEM(
                block_entries
                *
                sizeof(LIS_SCALAR)
                +
                pivot_entries
                *
                sizeof(LIS_INT)
                +
                (size_t)level->block_count
                *
                sizeof(LIS_INT)
                +
                (size_t)level->block_size
                *
                (
                    sizeof(LIS_SCALAR)
                    +
                    sizeof(LIS_REAL)
                ));

            err =
                LIS_OUT_OF_MEMORY;

            goto cleanup;
        }


        err =
            lis_mfamg_smoother_refresh(
                level);

        if(err)
            goto cleanup;


        return LIS_SUCCESS;


    default:

        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG smoother type is not implemented\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }


cleanup:

    lis_mfamg_smoother_destroy(
        level);


    return err;
}



static LIS_INT
lis_mfamg_smoother_refresh(
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_INT err;


    if(
        level==NULL
        ||
        level->A==NULL
        ||
        level->n<=0
        ||
        level->A->n!=level->n
        ||
        level->diagonal==NULL
        ||
        level->ax==NULL
        ||
        level->diagonal->n!=level->n
        ||
        level->ax->n!=level->n)
    {
        return LIS_ERR_ILL_ARG;
    }


    switch(level->smoother_type)
    {
    case LIS_MFAMG_SMOOTHER_WEIGHTED_JACOBI:

        err =
            lis_matrix_get_diagonal(
                level->A,
                level->diagonal);

        if(err)
            return err;


        err =
            lis_mfamg_estimate_smoother_omega(
                level->A,
                level->diagonal,
                &level->omega);

        if(err)
            return err;


        return
            lis_mfamg_smoother_validate_diagonal(
                level);


    case LIS_MFAMG_SMOOTHER_BLOCK_JACOBI:

        return
            lis_mfamg_block_refresh(
                level);


    default:

        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG smoother type is not implemented\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }
}



static LIS_INT
lis_mfamg_smoother_apply(
    struct LIS_MFAMG_LEVEL *level,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_INT initial_zero)
{
    LIS_INT i;
    LIS_INT block;
    LIS_INT start;
    LIS_INT width;
    LIS_INT err;

    LIS_SCALAR residual;
    LIS_SCALAR correction;

    LIS_SCALAR *bv;
    LIS_SCALAR *xv;
    LIS_SCALAR *av;
    LIS_SCALAR *diag;
    LIS_SCALAR *lu;
    LIS_SCALAR *work;


    if(
        level==NULL
        ||
        level->A==NULL
        ||
        b==NULL
        ||
        x==NULL
        ||
        level->n<=0
        ||
        b->n!=level->n
        ||
        x->n!=level->n
        ||
        level->diagonal==NULL
        ||
        level->ax==NULL
        ||
        level->diagonal->n!=level->n
        ||
        level->ax->n!=level->n
        ||
        (
            initial_zero!=LIS_FALSE
            &&
            initial_zero!=LIS_TRUE
        ))
    {
        return LIS_ERR_ILL_ARG;
    }


    bv =
        b->value;

    xv =
        x->value;

    diag =
        level->diagonal->value;


    switch(level->smoother_type)
    {
    case LIS_MFAMG_SMOOTHER_WEIGHTED_JACOBI:

        if(initial_zero)
        {
            for(i=0;i<level->n;i++)
            {
                xv[i] =
                    level->omega
                    *
                    bv[i]
                    /
                    diag[i];
            }


            return LIS_SUCCESS;
        }


        err =
            lis_matvec(
                level->A,
                x,
                level->ax);

        if(err)
            return err;


        av =
            level->ax->value;


        for(i=0;i<level->n;i++)
        {
            residual =
                bv[i]
                -
                av[i];

            xv[i] +=
                level->omega
                *
                residual
                /
                diag[i];
        }


        return LIS_SUCCESS;


    case LIS_MFAMG_SMOOTHER_BLOCK_JACOBI:

        if(
            level->block_size<=1
            ||
            level->block_size>LIS_MFAMG_BLOCK_MAX_SIZE
            ||
            level->block_count<=0
            ||
            level->block_lu==NULL
            ||
            level->block_pivots==NULL
            ||
            level->block_valid==NULL
            ||
            level->block_work==NULL)
        {
            return LIS_ERR_ILL_ARG;
        }


        av =
            NULL;


        if(!initial_zero)
        {
            err =
                lis_matvec(
                    level->A,
                    x,
                    level->ax);

            if(err)
                return err;


            av =
                level->ax->value;
        }


        work =
            level->block_work;


        for(
            block=0;
            block<level->block_count;
            block++)
        {
            start =
                block
                *
                level->block_size;

            width =
                level->n
                -
                start;

            if(width>level->block_size)
            {
                width =
                    level->block_size;
            }


            /*
             * Gather the complete block residual before any
             * x entry in this block is modified.  This also
             * preserves the existing alias-safe behavior.
             */
            for(i=0;i<width;i++)
            {
                if(initial_zero)
                {
                    work[i] =
                        bv[start+i];
                }
                else
                {
                    work[i] =
                        bv[start+i]
                        -
                        av[start+i];
                }
            }


            if(level->block_valid[block])
            {
                lu =
                    level->block_lu
                    +
                    (size_t)block
                    *
                    (size_t)level->block_size
                    *
                    (size_t)level->block_size;


                err =
                    lis_mfamg_block_solve(
                        width,
                        level->block_size,
                        lu,
                        level->block_pivots
                            +
                            block
                            *
                            level->block_size,
                        work);

                if(err)
                    return err;


                for(i=0;i<width;i++)
                {
                    correction =
                        level->omega
                        *
                        work[i];

                    if(initial_zero)
                    {
                        xv[start+i] =
                            correction;
                    }
                    else
                    {
                        xv[start+i] +=
                            correction;
                    }
                }
            }
            else
            {
                for(i=0;i<width;i++)
                {
                    correction =
                        level->omega
                        *
                        work[i]
                        /
                        diag[start+i];

                    if(initial_zero)
                    {
                        xv[start+i] =
                            correction;
                    }
                    else
                    {
                        xv[start+i] +=
                            correction;
                    }
                }
            }
        }


        return LIS_SUCCESS;


    default:

        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG smoother type is not implemented\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }
}



LIS_INT
lis_precon_mfamg_apply_block_smoother_for_test(
    LIS_MATRIX A,
    LIS_INT block_size,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_INT initial_zero,
    LIS_INT *active_blocks,
    LIS_INT *fallback_blocks,
    LIS_REAL *omega)
{
    struct LIS_MFAMG_LEVEL level;

    LIS_INT block;
    LIS_INT err;


    if(
        A==NULL
        ||
        b==NULL
        ||
        x==NULL
        ||
        active_blocks==NULL
        ||
        fallback_blocks==NULL
        ||
        omega==NULL
        ||
        block_size<=1
        ||
        block_size>LIS_MFAMG_BLOCK_MAX_SIZE
        ||
        A->n<=0
        ||
        b->n!=A->n
        ||
        x->n!=A->n)
    {
        return LIS_ERR_ILL_ARG;
    }


    memset(
        &level,
        0,
        sizeof(level));


    level.A =
        A;

    level.n =
        A->n;

    level.smoother_type =
        LIS_MFAMG_SMOOTHER_BLOCK_JACOBI;

    level.block_size =
        block_size;


    err =
        lis_mfamg_smoother_setup(
            &level);

    if(err)
        goto cleanup;


    *active_blocks =
        0;

    *fallback_blocks =
        0;


    for(
        block=0;
        block<level.block_count;
        block++)
    {
        if(level.block_valid[block])
        {
            (*active_blocks)++;
        }
        else
        {
            (*fallback_blocks)++;
        }
    }


    *omega =
        level.omega;


    err =
        lis_mfamg_smoother_apply(
            &level,
            b,
            x,
            initial_zero);


cleanup:

    lis_mfamg_smoother_destroy(
        &level);


    return err;
}



/*
 * Internal MFAMG global-row query used by the regression test.
 *
 * local_row is local to the calling process. Returned columns
 * always use zero-based global numbering for both CSR and USER.
 */
LIS_INT
lis_precon_mfamg_get_global_row(
    LIS_MATRIX A,
    LIS_INT local_row,
    LIS_INT *nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    if(nnz==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        (cols==NULL)
        !=
        (vals==NULL))
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        cols==NULL
        &&
        vals==NULL)
    {
        return
            lis_mfamg_row_query(
                A,
                local_row,
                nnz);
    }


    if(*nnz<0)
    {
        return LIS_ERR_ILL_ARG;
    }


    return
        lis_mfamg_row_fill(
            A,
            local_row,
            *nnz,
            nnz,
            cols,
            vals);
}

#undef __FUNC__
#define __FUNC__ "lis_mfamg_dense_factor"
static LIS_INT
lis_mfamg_dense_factor(
    LIS_INT n,
    LIS_SCALAR *lu,
    LIS_INT *pivots)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT p;

    LIS_REAL *row_scale;
    LIS_REAL best_ratio;
    LIS_REAL ratio;
    LIS_REAL mag;
    LIS_REAL scale;
    LIS_REAL tmp_scale;

    LIS_SCALAR tmp;


    row_scale =
        (LIS_REAL *)lis_malloc(
            (size_t)n
                *
            sizeof(LIS_REAL),
            "lis_mfamg_dense_factor::row_scale");

    if(row_scale==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)n
                *
            sizeof(LIS_REAL));

        return LIS_OUT_OF_MEMORY;
    }


    for(i=0;i<n;i++)
    {
        scale = 0.0;

        for(j=0;j<n;j++)
        {
            mag =
                fabs(
                    lu[i*n+j]);

            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                lis_free(
                    row_scale);

                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG coarse matrix contains non-finite values\n");

                return LIS_BREAKDOWN;
            }

            if(mag>scale)
                scale = mag;
        }

        /*
         * A zero row is an exact null direction.  Keep its
         * zero scale and let the relative pivot test below
         * classify the corresponding pivot.
         */
        row_scale[i] =
            scale;
    }


    for(k=0;k<n;k++)
    {
        p = -1;
        best_ratio = -1.0;

        for(i=k;i<n;i++)
        {
            mag =
                fabs(
                    lu[i*n+k]);

            if(row_scale[i]>0.0)
            {
                ratio =
                    mag
                    /
                    row_scale[i];
            }
            else
            {
                ratio =
                    0.0;
            }

            if(ratio>best_ratio)
            {
                best_ratio = ratio;
                p = i;
            }
        }

        if(p<0)
        {
            lis_free(
                row_scale);

            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse pivot search failed\n");

            return LIS_BREAKDOWN;
        }


        /*
         * Scaled partial pivoting already provides a
         * dimensionless measure of pivot quality:
         *
         *     |pivot| / row_scale
         *
         * Treat a sufficiently small ratio as a numerical
         * null direction.  Store an exact zero diagonal as
         * the canonical marker consumed by the solve path.
         *
         * This factorization records the null pivot.  The
         * solve path switches to a rank-revealing row-and-
         * column fallback whenever such a marker is present.
         */
        if(
            best_ratio
            <=
            LIS_MFAMG_COARSE_PIVOT_TOL)
        {
            pivots[k] =
                k;

            lu[k*n+k] =
                0.0;

            for(i=k+1;i<n;i++)
            {
                lu[i*n+k] =
                    0.0;
            }

            continue;
        }


        pivots[k] =
            p;

        if(p!=k)
        {
            for(j=0;j<n;j++)
            {
                tmp =
                    lu[k*n+j];

                lu[k*n+j] =
                    lu[p*n+j];

                lu[p*n+j] =
                    tmp;
            }

            tmp_scale =
                row_scale[k];

            row_scale[k] =
                row_scale[p];

            row_scale[p] =
                tmp_scale;
        }

        for(i=k+1;i<n;i++)
        {
            lu[i*n+k] /=
                lu[k*n+k];

            mag =
                fabs(
                    lu[i*n+k]);

            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                lis_free(
                    row_scale);

                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG coarse factorization became non-finite\n");

                return LIS_BREAKDOWN;
            }

            for(j=k+1;j<n;j++)
            {
                lu[i*n+j] -=
                    lu[i*n+k]
                    *
                    lu[k*n+j];

                mag =
                    fabs(
                        lu[i*n+j]);

                if(
                    mag!=mag
                    ||
                    mag>LIS_SCALAR_MAX)
                {
                    lis_free(
                        row_scale);

                    LIS_SETERR(
                        LIS_BREAKDOWN,
                        "MFAMG coarse factorization became non-finite\n");

                    return LIS_BREAKDOWN;
                }
            }
        }
    }


    lis_free(
        row_scale);

    return LIS_SUCCESS;
}



#undef __FUNC__
#define __FUNC__ "lis_mfamg_dense_solve"
static LIS_INT
lis_mfamg_dense_certify_null_solve(
    LIS_INT n,
    const LIS_SCALAR *matrix,
    const LIS_SCALAR *rhs,
    const LIS_SCALAR *x)
{
    LIS_INT i;
    LIS_INT j;

    LIS_REAL matrix_mag;
    LIS_REAL x_mag;
    LIS_REAL residual_mag;
    LIS_REAL reference;
    LIS_REAL limit;

    LIS_SCALAR ax;
    LIS_SCALAR residual;


    if(
        n<=0
        ||
        matrix==NULL
        ||
        rhs==NULL
        ||
        x==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    for(i=0;i<n;i++)
    {
        ax = 0.0;
        reference =
            fabs(
                rhs[i]);


        if(
            reference!=reference
            ||
            reference>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse RHS contains non-finite values\n");

            return LIS_BREAKDOWN;
        }


        for(j=0;j<n;j++)
        {
            matrix_mag =
                fabs(
                    matrix[i*n+j]);

            x_mag =
                fabs(
                    x[j]);

            if(
                matrix_mag!=matrix_mag
                ||
                matrix_mag>LIS_SCALAR_MAX
                ||
                x_mag!=x_mag
                ||
                x_mag>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG coarse residual check encountered non-finite values\n");

                return LIS_BREAKDOWN;
            }


            ax +=
                matrix[i*n+j]
                *
                x[j];

            reference +=
                matrix_mag
                *
                x_mag;

            if(
                reference!=reference
                ||
                reference>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG coarse residual scale became non-finite\n");

                return LIS_BREAKDOWN;
            }
        }


        residual =
            rhs[i]
            -
            ax;

        residual_mag =
            fabs(
                residual);

        if(
            residual_mag!=residual_mag
            ||
            residual_mag>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse residual became non-finite\n");

            return LIS_BREAKDOWN;
        }


        if(reference==0.0)
        {
            if(residual_mag!=0.0)
            {
                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG singular coarse system has an incompatible RHS\n");

                return LIS_BREAKDOWN;
            }

            continue;
        }


        limit =
            ((LIS_REAL)32.0)
            *
            LIS_MFAMG_COARSE_PIVOT_TOL
            *
            reference;


        if(residual_mag>limit)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG singular coarse system has an incompatible RHS\n");

            return LIS_BREAKDOWN;
        }
    }


    return LIS_SUCCESS;
}



static LIS_INT
lis_mfamg_dense_rank_revealing_solve(
    LIS_INT n,
    const LIS_SCALAR *matrix,
    const LIS_SCALAR *rhs,
    LIS_SCALAR *x,
    LIS_INT *nullity)
{
    LIS_SCALAR *work;
    LIS_SCALAR *work_rhs;
    LIS_SCALAR *y;

    LIS_REAL *row_scale;

    LIS_INT *col_perm;

    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT p;
    LIS_INT q;
    LIS_INT rank;
    LIS_INT err;
    LIS_INT tmp_index;

    LIS_REAL best_ratio;
    LIS_REAL ratio;
    LIS_REAL mag;
    LIS_REAL scale;
    LIS_REAL tmp_scale;

    LIS_SCALAR pivot;
    LIS_SCALAR factor;
    LIS_SCALAR sum;
    LIS_SCALAR tmp;


    if(
        n<=0
        ||
        matrix==NULL
        ||
        rhs==NULL
        ||
        x==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(nullity)
    {
        *nullity =
            0;
    }


    work = NULL;
    work_rhs = NULL;
    y = NULL;
    row_scale = NULL;
    col_perm = NULL;


    work =
        (LIS_SCALAR *)lis_malloc(
            (size_t)n
            *
            (size_t)n
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_dense_rank_revealing_solve::work");

    work_rhs =
        (LIS_SCALAR *)lis_malloc(
            (size_t)n
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_dense_rank_revealing_solve::work_rhs");

    y =
        (LIS_SCALAR *)lis_malloc(
            (size_t)n
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_dense_rank_revealing_solve::y");

    row_scale =
        (LIS_REAL *)lis_malloc(
            (size_t)n
            *
            sizeof(LIS_REAL),
            "lis_mfamg_dense_rank_revealing_solve::row_scale");

    col_perm =
        (LIS_INT *)lis_malloc(
            (size_t)n
            *
            sizeof(LIS_INT),
            "lis_mfamg_dense_rank_revealing_solve::col_perm");


    if(
        work==NULL
        ||
        work_rhs==NULL
        ||
        y==NULL
        ||
        row_scale==NULL
        ||
        col_perm==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)n
            *
            (
                (size_t)n
                *
                sizeof(LIS_SCALAR)
                +
                (size_t)2
                *
                sizeof(LIS_SCALAR)
                +
                sizeof(LIS_REAL)
                +
                sizeof(LIS_INT)
            ));

        err =
            LIS_OUT_OF_MEMORY;

        goto cleanup;
    }


    for(i=0;i<n;i++)
    {
        mag =
            fabs(
                rhs[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse RHS contains non-finite values\n");

            err =
                LIS_BREAKDOWN;

            goto cleanup;
        }


        work_rhs[i] =
            rhs[i];

        y[i] =
            0.0;

        col_perm[i] =
            i;

        scale =
            0.0;


        for(j=0;j<n;j++)
        {
            work[i*n+j] =
                matrix[i*n+j];

            mag =
                fabs(
                    work[i*n+j]);

            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG coarse matrix contains non-finite values\n");

                err =
                    LIS_BREAKDOWN;

                goto cleanup;
            }


            if(mag>scale)
            {
                scale =
                    mag;
            }
        }


        row_scale[i] =
            scale;
    }


    rank =
        0;


    for(k=0;k<n;k++)
    {
        p =
            -1;

        q =
            -1;

        best_ratio =
            -1.0;


        /*
         * Complete pivot search over the remaining submatrix,
         * scaled by the original magnitude of each row.
         *
         * Row and column pivoting are both required here:
         * a free early variable must not consume a row that
         * contains an independent equation in a later column.
         */
        for(i=k;i<n;i++)
        {
            if(row_scale[i]==0.0)
            {
                continue;
            }


            for(j=k;j<n;j++)
            {
                mag =
                    fabs(
                        work[i*n+j]);

                ratio =
                    mag
                    /
                    row_scale[i];

                if(
                    ratio!=ratio
                    ||
                    ratio>LIS_SCALAR_MAX)
                {
                    LIS_SETERR(
                        LIS_BREAKDOWN,
                        "MFAMG rank-revealing pivot search became non-finite\n");

                    err =
                        LIS_BREAKDOWN;

                    goto cleanup;
                }


                if(ratio>best_ratio)
                {
                    best_ratio =
                        ratio;

                    p =
                        i;

                    q =
                        j;
                }
            }
        }


        if(
            p<0
            ||
            q<0
            ||
            best_ratio
            <=
            LIS_MFAMG_COARSE_PIVOT_TOL)
        {
            break;
        }


        if(p!=k)
        {
            for(j=0;j<n;j++)
            {
                tmp =
                    work[k*n+j];

                work[k*n+j] =
                    work[p*n+j];

                work[p*n+j] =
                    tmp;
            }


            tmp =
                work_rhs[k];

            work_rhs[k] =
                work_rhs[p];

            work_rhs[p] =
                tmp;


            tmp_scale =
                row_scale[k];

            row_scale[k] =
                row_scale[p];

            row_scale[p] =
                tmp_scale;
        }


        if(q!=k)
        {
            for(i=0;i<n;i++)
            {
                tmp =
                    work[i*n+k];

                work[i*n+k] =
                    work[i*n+q];

                work[i*n+q] =
                    tmp;
            }


            tmp_index =
                col_perm[k];

            col_perm[k] =
                col_perm[q];

            col_perm[q] =
                tmp_index;
        }


        pivot =
            work[k*n+k];

        mag =
            fabs(
                pivot);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX
            ||
            mag==0.0)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG rank-revealing pivot is invalid\n");

            err =
                LIS_BREAKDOWN;

            goto cleanup;
        }


        for(i=k+1;i<n;i++)
        {
            factor =
                work[i*n+k]
                /
                pivot;

            mag =
                fabs(
                    factor);

            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG rank-revealing elimination became non-finite\n");

                err =
                    LIS_BREAKDOWN;

                goto cleanup;
            }


            work[i*n+k] =
                0.0;


            for(j=k+1;j<n;j++)
            {
                work[i*n+j] -=
                    factor
                    *
                    work[k*n+j];

                mag =
                    fabs(
                        work[i*n+j]);

                if(
                    mag!=mag
                    ||
                    mag>LIS_SCALAR_MAX)
                {
                    LIS_SETERR(
                        LIS_BREAKDOWN,
                        "MFAMG rank-revealing matrix update became non-finite\n");

                    err =
                        LIS_BREAKDOWN;

                    goto cleanup;
                }
            }


            work_rhs[i] -=
                factor
                *
                work_rhs[k];

            mag =
                fabs(
                    work_rhs[i]);

            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "MFAMG rank-revealing RHS update became non-finite\n");

                err =
                    LIS_BREAKDOWN;

                goto cleanup;
            }
        }


        rank++;
    }


    for(i=0;i<n;i++)
    {
        y[i] =
            0.0;
    }


    for(i=rank-1;i>=0;i--)
    {
        sum =
            work_rhs[i];


        for(j=i+1;j<rank;j++)
        {
            sum -=
                work[i*n+j]
                *
                y[j];
        }


        pivot =
            work[i*n+i];

        mag =
            fabs(
                pivot);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX
            ||
            mag==0.0)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG rank-revealing backsolve encountered an invalid pivot\n");

            err =
                LIS_BREAKDOWN;

            goto cleanup;
        }


        y[i] =
            sum
            /
            pivot;

        mag =
            fabs(
                y[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG rank-revealing coarse solve became non-finite\n");

            err =
                LIS_BREAKDOWN;

            goto cleanup;
        }
    }


    for(i=0;i<n;i++)
    {
        x[i] =
            0.0;
    }


    for(i=0;i<n;i++)
    {
        if(
            col_perm[i]<0
            ||
            col_perm[i]>=n)
        {
            err =
                LIS_ERR_ILL_ARG;

            goto cleanup;
        }


        x[col_perm[i]] =
            y[i];
    }


    if(nullity)
    {
        *nullity =
            n
            -
            rank;
    }


    err =
        lis_mfamg_dense_certify_null_solve(
            n,
            matrix,
            rhs,
            x);


cleanup:

    if(work)
        lis_free(work);

    if(work_rhs)
        lis_free(work_rhs);

    if(y)
        lis_free(y);

    if(row_scale)
        lis_free(row_scale);

    if(col_perm)
        lis_free(col_perm);


    return err;
}



static LIS_INT
lis_mfamg_dense_solve_values(
    LIS_INT n,
    const LIS_SCALAR *matrix,
    const LIS_SCALAR *lu,
    const LIS_INT *pivots,
    const LIS_SCALAR *rhs,
    LIS_SCALAR *x)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT p;

    LIS_REAL mag;

    LIS_SCALAR tmp;
    LIS_SCALAR sum;


    if(
        n<=0
        ||
        matrix==NULL
        ||
        lu==NULL
        ||
        pivots==NULL
        ||
        rhs==NULL
        ||
        x==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    /*
     * Preserve the normal LU fast path.  A canonical zero
     * diagonal is the only trigger for the more expensive
     * rank-revealing singular fallback.
     */
    for(i=0;i<n;i++)
    {
        mag =
            fabs(
                lu[i*n+i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse solve contains a non-finite pivot\n");

            return LIS_BREAKDOWN;
        }


        if(mag==0.0)
        {
            return
                lis_mfamg_dense_rank_revealing_solve(
                    n,
                    matrix,
                    rhs,
                    x,
                    NULL);
        }
    }


    memcpy(
        x,
        rhs,
        (size_t)n
        *
        sizeof(LIS_SCALAR));


    for(k=0;k<n;k++)
    {
        p =
            pivots[k];

        if(
            p<0
            ||
            p>=n)
        {
            return LIS_ERR_ILL_ARG;
        }


        if(p!=k)
        {
            tmp =
                x[k];

            x[k] =
                x[p];

            x[p] =
                tmp;
        }
    }


    for(i=1;i<n;i++)
    {
        sum =
            x[i];

        for(j=0;j<i;j++)
        {
            sum -=
                lu[i*n+j]
                *
                x[j];
        }

        x[i] =
            sum;
    }


    for(i=n-1;i>=0;i--)
    {
        sum =
            x[i];

        for(j=i+1;j<n;j++)
        {
            sum -=
                lu[i*n+j]
                *
                x[j];
        }


        mag =
            fabs(
                sum);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse solve residual became non-finite\n");

            return LIS_BREAKDOWN;
        }


        x[i] =
            sum
            /
            lu[i*n+i];


        mag =
            fabs(
                x[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse solve became non-finite\n");

            return LIS_BREAKDOWN;
        }
    }


    return LIS_SUCCESS;
}



static LIS_INT
lis_mfamg_dense_solve(
    struct LIS_MFAMG_STRUCT *state,
    LIS_VECTOR rhs,
    LIS_VECTOR x)
{
    if(
        state==NULL
        ||
        rhs==NULL
        ||
        x==NULL
        ||
        state->coarse_n<=0
        ||
        state->coarse==NULL
        ||
        state->lu==NULL
        ||
        state->pivots==NULL
        ||
        rhs->n!=state->coarse_n
        ||
        x->n!=state->coarse_n)
    {
        return LIS_ERR_ILL_ARG;
    }


    return
        lis_mfamg_dense_solve_values(
            state->coarse_n,
            state->coarse,
            state->lu,
            state->pivots,
            rhs->value,
            x->value);
}


static LIS_INT
lis_mfamg_build_aggregates_sized(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level,
    LIS_INT target_size,
    LIS_INT max_row_nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    LIS_VECTOR diagonal;

    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT nnz;
    LIS_INT best;
    LIS_INT fallback_best;
    LIS_INT aggregate;
    LIS_INT aggregate_size;
    LIS_INT member_id;
    LIS_INT row;
    LIS_INT err;

    LIS_INT *members;

    LIS_REAL best_mag;
    LIS_REAL fallback_mag;
    LIS_REAL mag;
    LIS_REAL diagonal_i;
    LIS_REAL diagonal_j;
    LIS_REAL strength_limit;


    if(
        A==NULL
        ||
        level==NULL
        ||
        level->aggregate==NULL
        ||
        level->n!=A->n
        ||
        level->n<=0
        ||
        target_size<=0
        ||
        max_row_nnz<=0
        ||
        cols==NULL
        ||
        vals==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    diagonal = NULL;
    members = NULL;

    if(target_size>level->n)
        target_size = level->n;


    err =
        lis_vector_duplicate(
            A,
            &diagonal);

    if(err)
        return err;


    err =
        lis_matrix_get_diagonal(
            A,
            diagonal);

    if(err)
    {
        lis_vector_destroy(
            diagonal);

        return err;
    }


    for(i=0;i<level->n;i++)
    {
        diagonal_i =
            fabs(
                diagonal->value[i]);

        if(
            diagonal_i!=diagonal_i
            ||
            diagonal_i>LIS_SCALAR_MAX
            ||
            diagonal_i==0.0)
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG aggregation requires finite nonzero diagonal entries\n");

            lis_vector_destroy(
                diagonal);

            return LIS_ERR_ILL_ARG;
        }
    }


    if(
        (size_t)target_size
        >
        ((size_t)-1)
            /
        sizeof(LIS_INT))
    {
        err = LIS_OUT_OF_MEMORY;
        goto cleanup;
    }


    members =
        (LIS_INT *)lis_malloc(
            (size_t)target_size
                *
            sizeof(LIS_INT),
            "lis_mfamg_build_aggregates_sized::members");

    if(members==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)target_size
                *
            sizeof(LIS_INT));

        err = LIS_OUT_OF_MEMORY;
        goto cleanup;
    }


    for(i=0;i<level->n;i++)
        level->aggregate[i] = -1;


    aggregate = 0;
    err = LIS_SUCCESS;


    for(i=0;i<level->n;i++)
    {
        if(level->aggregate[i]>=0)
            continue;


        level->aggregate[i] =
            aggregate;

        members[0] =
            i;

        aggregate_size =
            1;


        while(aggregate_size<target_size)
        {
            best = -1;
            fallback_best = -1;

            best_mag = 0.0;
            fallback_mag = 0.0;


            for(
                member_id=0;
                member_id<aggregate_size;
                member_id++)
            {
                row =
                    members[member_id];

                nnz =
                    max_row_nnz;

                err =
                    lis_mfamg_row_fill(
                        A,
                        row,
                        max_row_nnz,
                        &nnz,
                        cols,
                        vals);

                if(err)
                    goto cleanup;


                diagonal_i =
                    fabs(
                        diagonal->value[row]);


                for(k=0;k<nnz;k++)
                {
                    j =
                        cols[k]
                        -
                        A->is;


                    /*
                     * Aggregates remain rank-owned. Remote
                     * columns are not grown here; distributed
                     * metadata is exchanged after local
                     * aggregation.
                     */
                    if(
                        j<0
                        ||
                        j>=level->n
                        ||
                        j==row
                        ||
                        level->aggregate[j]>=0)
                    {
                        continue;
                    }


                    mag =
                        fabs(
                            vals[k]);

                    if(
                        mag!=mag
                        ||
                        mag>LIS_SCALAR_MAX)
                    {
                        LIS_SETERR(
                            LIS_ERR_ILL_ARG,
                            "MFAMG aggregation encountered a non-finite matrix entry\n");

                        err =
                            LIS_ERR_ILL_ARG;

                        goto cleanup;
                    }


                    diagonal_j =
                        fabs(
                            diagonal->value[j]);

                    strength_limit =
                        LIS_MFAMG_STRENGTH_THETA
                        *
                        sqrt(
                            diagonal_i)
                        *
                        sqrt(
                            diagonal_j);

                    if(
                        strength_limit!=strength_limit
                        ||
                        strength_limit>LIS_SCALAR_MAX)
                    {
                        LIS_SETERR(
                            LIS_ERR_ILL_ARG,
                            "MFAMG aggregation strength threshold overflowed\n");

                        err =
                            LIS_ERR_ILL_ARG;

                        goto cleanup;
                    }


                    /*
                     * Preserve the original strength-based
                     * aggregation whenever a strong connection
                     * is available.  For larger multivector
                     * aggregates, remember the strongest local
                     * unused graph edge as a fallback so that
                     * aggregate growth does not stop solely
                     * because the coarse operator weakens the
                     * strength test.
                     */
                    if(
                        target_size>2
                        &&
                        mag>fallback_mag)
                    {
                        fallback_mag =
                            mag;

                        fallback_best =
                            j;
                    }


                    if(mag<strength_limit)
                        continue;


                    if(mag>best_mag)
                    {
                        best_mag =
                            mag;

                        best =
                            j;
                    }
                }
            }


            if(
                best<0
                &&
                target_size>2)
            {
                best =
                    fallback_best;
            }


            if(best<0)
                break;


            level->aggregate[best] =
                aggregate;

            members[aggregate_size] =
                best;

            aggregate_size++;
        }


        aggregate++;
    }


    level->nc =
        aggregate;


    if(level->nc<=0)
        err = LIS_BREAKDOWN;


cleanup:

    if(members)
        lis_free(
            members);

    lis_vector_destroy(
        diagonal);


    return err;
}

static LIS_INT
lis_mfamg_build_aggregates(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level,
    LIS_INT max_row_nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    return
        lis_mfamg_build_aggregates_sized(
            A,
            level,
            2,
            max_row_nnz,
            cols,
            vals);
}




/*
 * Exchange scalar metadata using CSR local/halo indexing.
 *
 * This intentionally does not use lis_send_recv().  That
 * routine uses the padded LIS vector layout, while this array
 * is indexed directly by matrix-local CSR column indices
 * 0..A->np-1.
 */
static LIS_INT
lis_mfamg_exchange_scalar_halo(
    LIS_MATRIX A,
    LIS_SCALAR *x)
{
#ifdef USE_MPI

    LIS_COMMTABLE commtable;

    LIS_INT neib;
    LIS_INT i;
    LIS_INT is;
    LIS_INT count;

    LIS_INT send_count;
    LIS_INT recv_count;
    LIS_INT send_alloc;
    LIS_INT recv_alloc;

    LIS_SCALAR *ws;
    LIS_SCALAR *wr;


    if(
        A==NULL
        ||
        x==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        A->nprocs<=1
        ||
        A->np<=A->n)
    {
        return LIS_SUCCESS;
    }


    commtable =
        A->commtable;


    if(
        commtable==NULL
        ||
        commtable->neibpetot<=0
        ||
        commtable->export_ptr==NULL
        ||
        commtable->import_ptr==NULL
        ||
        commtable->export_index==NULL
        ||
        commtable->import_index==NULL
        ||
        commtable->neibpe==NULL
        ||
        commtable->req1==NULL
        ||
        commtable->req2==NULL
        ||
        commtable->sta1==NULL
        ||
        commtable->sta2==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        commtable->export_ptr[0]!=0
        ||
        commtable->import_ptr[0]!=0)
    {
        return LIS_ERR_ILL_ARG;
    }


    for(
        neib=0;
        neib<commtable->neibpetot;
        neib++)
    {
        if(
            commtable->export_ptr[neib]<0
            ||
            commtable->export_ptr[neib+1]
                <
            commtable->export_ptr[neib]
            ||
            commtable->import_ptr[neib]<0
            ||
            commtable->import_ptr[neib+1]
                <
            commtable->import_ptr[neib]
            ||
            commtable->neibpe[neib]<0
            ||
            commtable->neibpe[neib]>=A->nprocs)
        {
            return LIS_ERR_ILL_ARG;
        }
    }


    send_count =
        commtable->export_ptr[
            commtable->neibpetot];

    recv_count =
        commtable->import_ptr[
            commtable->neibpetot];


    if(
        send_count<0
        ||
        recv_count<0)
    {
        return LIS_ERR_ILL_ARG;
    }


    /*
     * Validate the complete communication topology before
     * posting any nonblocking MPI request.
     */
    for(i=0;i<send_count;i++)
    {
        if(
            commtable->export_index[i]<0
            ||
            commtable->export_index[i]>=A->n)
        {
            return LIS_ERR_ILL_ARG;
        }
    }


    for(i=0;i<recv_count;i++)
    {
        if(
            commtable->import_index[i]<A->n
            ||
            commtable->import_index[i]>=A->np)
        {
            return LIS_ERR_ILL_ARG;
        }
    }


    send_alloc =
        send_count>0
        ?
        send_count
        :
        1;

    recv_alloc =
        recv_count>0
        ?
        recv_count
        :
        1;


    ws =
        (LIS_SCALAR *)lis_malloc(
            (size_t)send_alloc
                *
            sizeof(LIS_SCALAR),
            "lis_mfamg_exchange_scalar_halo::ws");

    wr =
        (LIS_SCALAR *)lis_malloc(
            (size_t)recv_alloc
                *
            sizeof(LIS_SCALAR),
            "lis_mfamg_exchange_scalar_halo::wr");


    if(
        ws==NULL
        ||
        wr==NULL)
    {
        if(ws)
            lis_free(ws);

        if(wr)
            lis_free(wr);

        LIS_SETERR_MEM(
            (size_t)(
                send_alloc
                +
                recv_alloc)
                *
            sizeof(LIS_SCALAR));

        return LIS_OUT_OF_MEMORY;
    }


    for(i=0;i<send_count;i++)
    {
        ws[i] =
            x[
                commtable->export_index[i]];
    }


    /*
     * From this point onward there are no locally generated
     * error exits until both request sets have completed.
     */
    for(
        neib=0;
        neib<commtable->neibpetot;
        neib++)
    {
        is =
            commtable->export_ptr[neib];

        count =
            commtable->export_ptr[neib+1]
            -
            is;


        MPI_Isend(
            &ws[is],
            count,
            LIS_MPI_SCALAR,
            commtable->neibpe[neib],
            0,
            commtable->comm,
            &commtable->req1[neib]);
    }


    for(
        neib=0;
        neib<commtable->neibpetot;
        neib++)
    {
        is =
            commtable->import_ptr[neib];

        count =
            commtable->import_ptr[neib+1]
            -
            is;


        MPI_Irecv(
            &wr[is],
            count,
            LIS_MPI_SCALAR,
            commtable->neibpe[neib],
            0,
            commtable->comm,
            &commtable->req2[neib]);
    }


    MPI_Waitall(
        commtable->neibpetot,
        commtable->req2,
        commtable->sta2);


    for(i=0;i<recv_count;i++)
    {
        x[
            commtable->import_index[i]]
            =
            wr[i];
    }


    MPI_Waitall(
        commtable->neibpetot,
        commtable->req1,
        commtable->sta1);


    lis_free(wr);
    lis_free(ws);


    return LIS_SUCCESS;

#else

    (void)A;
    (void)x;

    return LIS_SUCCESS;

#endif
}


/*
 * Exchange integer metadata using the matrix communication
 * topology.  export_index addresses owned fine entries and
 * import_index addresses the corresponding halo entries.
 */
static LIS_INT
lis_mfamg_exchange_int_halo(
    LIS_MATRIX A,
    LIS_INT *x)
{
#ifdef USE_MPI

    LIS_COMMTABLE commtable;

    LIS_INT neib;
    LIS_INT i;
    LIS_INT is;
    LIS_INT count;

    LIS_INT send_count;
    LIS_INT recv_count;
    LIS_INT send_alloc;
    LIS_INT recv_alloc;

    LIS_INT *ws;
    LIS_INT *wr;


    if(
        A==NULL
        ||
        x==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        A->nprocs<=1
        ||
        A->np<=A->n)
    {
        return LIS_SUCCESS;
    }


    commtable =
        A->commtable;


    if(
        commtable==NULL
        ||
        commtable->neibpetot<=0
        ||
        commtable->export_ptr==NULL
        ||
        commtable->import_ptr==NULL
        ||
        commtable->export_index==NULL
        ||
        commtable->import_index==NULL
        ||
        commtable->neibpe==NULL
        ||
        commtable->req1==NULL
        ||
        commtable->req2==NULL
        ||
        commtable->sta1==NULL
        ||
        commtable->sta2==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        commtable->export_ptr[0]!=0
        ||
        commtable->import_ptr[0]!=0)
    {
        return LIS_ERR_ILL_ARG;
    }


    for(
        neib=0;
        neib<commtable->neibpetot;
        neib++)
    {
        if(
            commtable->export_ptr[neib]<0
            ||
            commtable->export_ptr[neib+1]
                <
            commtable->export_ptr[neib]
            ||
            commtable->import_ptr[neib]<0
            ||
            commtable->import_ptr[neib+1]
                <
            commtable->import_ptr[neib]
            ||
            commtable->neibpe[neib]<0
            ||
            commtable->neibpe[neib]>=A->nprocs)
        {
            return LIS_ERR_ILL_ARG;
        }
    }


    send_count =
        commtable->export_ptr[
            commtable->neibpetot];

    recv_count =
        commtable->import_ptr[
            commtable->neibpetot];


    if(
        send_count<0
        ||
        recv_count<0)
    {
        return LIS_ERR_ILL_ARG;
    }


    /*
     * Validate the complete communication topology before
     * posting any nonblocking MPI request.
     */
    for(i=0;i<send_count;i++)
    {
        if(
            commtable->export_index[i]<0
            ||
            commtable->export_index[i]>=A->n)
        {
            return LIS_ERR_ILL_ARG;
        }
    }


    for(i=0;i<recv_count;i++)
    {
        if(
            commtable->import_index[i]<A->n
            ||
            commtable->import_index[i]>=A->np)
        {
            return LIS_ERR_ILL_ARG;
        }
    }


    send_alloc =
        send_count>0
        ?
        send_count
        :
        1;

    recv_alloc =
        recv_count>0
        ?
        recv_count
        :
        1;


    ws =
        (LIS_INT *)lis_malloc(
            (size_t)send_alloc
                *
            sizeof(LIS_INT),
            "lis_mfamg_exchange_int_halo::ws");

    wr =
        (LIS_INT *)lis_malloc(
            (size_t)recv_alloc
                *
            sizeof(LIS_INT),
            "lis_mfamg_exchange_int_halo::wr");


    if(
        ws==NULL
        ||
        wr==NULL)
    {
        if(ws)
            lis_free(ws);

        if(wr)
            lis_free(wr);

        LIS_SETERR_MEM(
            (size_t)(
                send_alloc
                +
                recv_alloc)
                *
            sizeof(LIS_INT));

        return LIS_OUT_OF_MEMORY;
    }


    for(i=0;i<send_count;i++)
    {
        ws[i] =
            x[
                commtable->export_index[i]];
    }


    /*
     * From this point onward there are no locally generated
     * error exits until both request sets have completed.
     */
    for(
        neib=0;
        neib<commtable->neibpetot;
        neib++)
    {
        is =
            commtable->export_ptr[neib];

        count =
            commtable->export_ptr[neib+1]
            -
            is;


        MPI_Isend(
            &ws[is],
            count,
            LIS_MPI_INT,
            commtable->neibpe[neib],
            0,
            commtable->comm,
            &commtable->req1[neib]);
    }


    for(
        neib=0;
        neib<commtable->neibpetot;
        neib++)
    {
        is =
            commtable->import_ptr[neib];

        count =
            commtable->import_ptr[neib+1]
            -
            is;


        MPI_Irecv(
            &wr[is],
            count,
            LIS_MPI_INT,
            commtable->neibpe[neib],
            0,
            commtable->comm,
            &commtable->req2[neib]);
    }


    MPI_Waitall(
        commtable->neibpetot,
        commtable->req2,
        commtable->sta2);


    for(i=0;i<recv_count;i++)
    {
        x[
            commtable->import_index[i]]
            =
            wr[i];
    }


    MPI_Waitall(
        commtable->neibpetot,
        commtable->req1,
        commtable->sta1);


    lis_free(wr);
    lis_free(ws);


    return LIS_SUCCESS;

#else

    (void)A;
    (void)x;

    return LIS_SUCCESS;

#endif
}


static LIS_INT
lis_mfamg_build_distributed_multivector_aggregate_metadata(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_INT local_nc;
    LIS_INT global_nc;
    LIS_INT coarse_is;
    LIS_INT coarse_ie;
    LIS_INT nprocs;
    LIS_INT my_rank;

    LIS_INT *coarse_ranges;
    LIS_INT *global_aggregate;
    LIS_INT *global_aggregate_rank;

    LIS_INT i;
    LIS_INT a;
    LIS_INT rank;
    LIS_INT offset;
    LIS_INT err;


    if(
        A==NULL
        ||
        level==NULL
        ||
        level->aggregate==NULL
        ||
        level->aggregate_rank==NULL
        ||
        level->aggregate_offset==NULL
        ||
        level->prolong_row_ptr==NULL
        ||
        level->n!=A->n
        ||
        level->aggregate_count<=0
        ||
        level->nc<=0
        ||
        level->aggregate_offset[
            level->aggregate_count]
            !=
        level->nc
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs<=1
        ||
        A->np<A->n
        ||
        level->coarse_ranges!=NULL
        ||
        level->global_aggregate!=NULL
        ||
        level->global_aggregate_rank!=NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    local_nc =
        level->nc;

    global_nc = 0;
    coarse_is = 0;
    coarse_ie = 0;
    nprocs = 0;
    my_rank = 0;

    coarse_ranges = NULL;
    global_aggregate = NULL;
    global_aggregate_rank = NULL;


    err =
        lis_ranges_create(
            A->comm,
            &local_nc,
            &global_nc,
            &coarse_ranges,
            &coarse_is,
            &coarse_ie,
            &nprocs,
            &my_rank);

    if(err)
        return err;


    if(
        coarse_ranges==NULL
        ||
        local_nc!=level->nc
        ||
        global_nc<=0
        ||
        coarse_is<0
        ||
        coarse_ie<coarse_is
        ||
        coarse_ie-coarse_is!=level->nc
        ||
        nprocs!=A->nprocs
        ||
        my_rank!=A->my_rank)
    {
        if(coarse_ranges)
            lis_free(coarse_ranges);

        return LIS_ERR_ILL_ARG;
    }


    global_aggregate =
        (LIS_INT *)lis_malloc(
            (size_t)A->np
            *
            sizeof(LIS_INT),
            "lis_mfamg_build_distributed_multivector_aggregate_metadata::global_aggregate");

    global_aggregate_rank =
        (LIS_INT *)lis_malloc(
            (size_t)A->np
            *
            sizeof(LIS_INT),
            "lis_mfamg_build_distributed_multivector_aggregate_metadata::global_aggregate_rank");


    if(
        global_aggregate==NULL
        ||
        global_aggregate_rank==NULL)
    {
        if(global_aggregate)
            lis_free(global_aggregate);

        if(global_aggregate_rank)
            lis_free(global_aggregate_rank);

        lis_free(coarse_ranges);

        LIS_SETERR_MEM(
            (size_t)A->np
            *
            2
            *
            sizeof(LIS_INT));

        return LIS_OUT_OF_MEMORY;
    }


    for(i=0;i<A->np;i++)
    {
        global_aggregate[i] = -1;
        global_aggregate_rank[i] = 0;
    }


    for(i=0;i<A->n;i++)
    {
        a =
            level->aggregate[i];

        if(
            a<0
            ||
            a>=level->aggregate_count)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        rank =
            level->aggregate_rank[a];

        offset =
            level->aggregate_offset[a];


        if(
            rank<0
            ||
            offset<0
            ||
            offset+rank>level->nc)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        global_aggregate[i] =
            coarse_is
            +
            offset;

        global_aggregate_rank[i] =
            rank;
    }


    if(A->np>A->n)
    {
        if(A->commtable==NULL)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        err =
            lis_mfamg_exchange_int_halo(
                A,
                global_aggregate);

        if(err)
            goto fail;


        err =
            lis_mfamg_exchange_int_halo(
                A,
                global_aggregate_rank);

        if(err)
            goto fail;
    }


    for(i=0;i<A->np;i++)
    {
        rank =
            global_aggregate_rank[i];

        if(
            rank<0
            ||
            global_aggregate[i]<0
            ||
            global_aggregate[i]>global_nc
            ||
            global_aggregate[i]+rank>global_nc)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }
    }


    level->global_nc =
        global_nc;

    level->coarse_is =
        coarse_is;

    level->coarse_ie =
        coarse_ie;

    level->coarse_ranges =
        coarse_ranges;

    level->global_aggregate =
        global_aggregate;

    level->global_aggregate_rank =
        global_aggregate_rank;


    return LIS_SUCCESS;


fail:

    if(global_aggregate_rank)
        lis_free(global_aggregate_rank);

    if(global_aggregate)
        lis_free(global_aggregate);

    if(coarse_ranges)
        lis_free(coarse_ranges);


    return err;
}


static LIS_INT
lis_mfamg_build_distributed_aggregate_metadata(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_INT local_nc;
    LIS_INT global_nc;
    LIS_INT coarse_is;
    LIS_INT coarse_ie;
    LIS_INT nprocs;
    LIS_INT my_rank;

    LIS_INT *coarse_ranges;
    LIS_INT *global_aggregate;

    LIS_INT i;
    LIS_INT err;


    if(
        level!=NULL
        &&
        level->prolong_row_ptr!=NULL)
    {
        return
            lis_mfamg_build_distributed_multivector_aggregate_metadata(
                A,
                level);
    }


    if(
        A==NULL
        ||
        level==NULL
        ||
        level->aggregate==NULL
        ||
        level->n!=A->n
        ||
        level->nc<=0
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs<=1
        ||
        level->coarse_ranges!=NULL
        ||
        level->global_aggregate!=NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    local_nc =
        level->nc;

    global_nc = 0;
    coarse_is = 0;
    coarse_ie = 0;
    nprocs = 0;
    my_rank = 0;
    coarse_ranges = NULL;
    global_aggregate = NULL;


    err =
        lis_ranges_create(
            A->comm,
            &local_nc,
            &global_nc,
            &coarse_ranges,
            &coarse_is,
            &coarse_ie,
            &nprocs,
            &my_rank);

    if(err)
        return err;


    if(
        coarse_ranges==NULL
        ||
        local_nc!=level->nc
        ||
        global_nc<=0
        ||
        coarse_is<0
        ||
        coarse_ie<coarse_is
        ||
        coarse_ie-coarse_is!=level->nc
        ||
        nprocs!=A->nprocs
        ||
        my_rank!=A->my_rank)
    {
        if(coarse_ranges)
            lis_free(coarse_ranges);

        return LIS_ERR_ILL_ARG;
    }


    global_aggregate =
        (LIS_INT *)lis_malloc(
            (size_t)A->np
                *
            sizeof(LIS_INT),
            "lis_mfamg_build_distributed_aggregate_metadata::global_aggregate");


    if(global_aggregate==NULL)
    {
        lis_free(
            coarse_ranges);

        LIS_SETERR_MEM(
            (size_t)A->np
                *
            sizeof(LIS_INT));

        return LIS_OUT_OF_MEMORY;
    }


    for(i=0;i<A->np;i++)
    {
        global_aggregate[i] = -1;
    }


    for(i=0;i<A->n;i++)
    {
        if(
            level->aggregate[i]<0
            ||
            level->aggregate[i]>=level->nc)
        {
            lis_free(
                global_aggregate);

            lis_free(
                coarse_ranges);

            return LIS_ERR_ILL_ARG;
        }


        global_aggregate[i] =
            coarse_is
            +
            level->aggregate[i];
    }


    if(A->np>A->n)
    {
        if(A->commtable==NULL)
        {
            lis_free(
                global_aggregate);

            lis_free(
                coarse_ranges);

            return LIS_ERR_ILL_ARG;
        }


        err =
            lis_mfamg_exchange_int_halo(
                A,
                global_aggregate);

        if(err)
        {
            lis_free(
                global_aggregate);

            lis_free(
                coarse_ranges);

            return err;
        }
    }


    for(i=0;i<A->np;i++)
    {
        if(
            global_aggregate[i]<0
            ||
            global_aggregate[i]>=global_nc)
        {
            lis_free(
                global_aggregate);

            lis_free(
                coarse_ranges);

            return LIS_ERR_ILL_ARG;
        }
    }


    level->global_nc =
        global_nc;

    level->coarse_is =
        coarse_is;

    level->coarse_ie =
        coarse_ie;

    level->coarse_ranges =
        coarse_ranges;

    level->global_aggregate =
        global_aggregate;


    return LIS_SUCCESS;
}



/*
 * Build an aggregate-local multivector tentative basis.
 *
 * For each aggregate a, the local near-nullspace block is
 * factorized as
 *
 *     B_a = Q_a R_a
 *
 * using two-pass modified Gram-Schmidt with numerical rank
 * detection.
 *
 * The local rank may differ between aggregates.  Coarse
 * numbering therefore uses prefix offsets:
 *
 *     coarse_n = sum_a rank(B_a)
 *
 * The returned arrays satisfy
 *
 *     B = P * Bc
 *
 * up to LIS_MFAMG_RANK_TOL.
 *
 * This routine provides the serial algebraic foundation for
 * multivector MFAMG.  It does not modify the current scalar
 * production hierarchy path.
 */
#undef __FUNC__
#define __FUNC__ "lis_precon_mfamg_build_multivector_basis"
LIS_INT
lis_precon_mfamg_build_multivector_basis(
    LIS_INT n,
    LIS_INT aggregate_count,
    LIS_INT *aggregate,
    LIS_INT mode_count,
    LIS_VECTOR basis[],
    LIS_INT **aggregate_rank_out,
    LIS_INT **aggregate_offset_out,
    LIS_INT **row_ptr_out,
    LIS_SCALAR **prolong_out,
    LIS_SCALAR **coarse_basis_out,
    LIS_INT *coarse_n_out)
{
    LIS_INT *counts;
    LIS_INT *member_offset;
    LIS_INT *cursor;
    LIS_INT *members;

    LIS_INT *aggregate_rank;
    LIS_INT *aggregate_offset;
    LIS_INT *row_ptr;

    LIS_SCALAR *q_all;
    LIS_SCALAR *r_all;
    LIS_SCALAR *work;
    LIS_SCALAR *prolong;
    LIS_SCALAR *coarse_basis;

    LIS_INT i;
    LIS_INT a;
    LIS_INT m;
    LIS_INT q;
    LIS_INT t;
    LIS_INT pass;
    LIS_INT rank;
    LIS_INT count;
    LIS_INT node;
    LIS_INT coarse_n;
    LIS_INT err;

    LIS_REAL mag;
    LIS_REAL norm2;
    LIS_REAL norm;
    LIS_REAL max_norm;
    LIS_REAL threshold;

    LIS_SCALAR dot;

    size_t q_entries;
    size_t r_entries;
    size_t prolong_entries;
    size_t coarse_entries;
    size_t prefix;
    size_t next_size;


    if(
        n<=0
        ||
        aggregate_count<=0
        ||
        aggregate==NULL
        ||
        mode_count<=0
        ||
        basis==NULL
        ||
        aggregate_rank_out==NULL
        ||
        aggregate_offset_out==NULL
        ||
        row_ptr_out==NULL
        ||
        prolong_out==NULL
        ||
        coarse_basis_out==NULL
        ||
        coarse_n_out==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    *aggregate_rank_out = NULL;
    *aggregate_offset_out = NULL;
    *row_ptr_out = NULL;
    *prolong_out = NULL;
    *coarse_basis_out = NULL;
    *coarse_n_out = 0;


    for(m=0;m<mode_count;m++)
    {
        if(
            basis[m]==NULL
            ||
            basis[m]->n!=n)
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG multivector basis has incompatible vector layout\n");

            return LIS_ERR_ILL_ARG;
        }
    }


    if(
        (size_t)n
        >
        ((size_t)-1)
        /
        (size_t)mode_count)
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG multivector basis size overflow\n");

        return LIS_OUT_OF_MEMORY;
    }


    q_entries =
        (size_t)n
        *
        (size_t)mode_count;


    if(
        (size_t)aggregate_count
        >
        ((size_t)-1)
        /
        (size_t)mode_count)
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG multivector R size overflow\n");

        return LIS_OUT_OF_MEMORY;
    }


    r_entries =
        (size_t)aggregate_count
        *
        (size_t)mode_count;


    if(
        r_entries
        >
        ((size_t)-1)
        /
        (size_t)mode_count)
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG multivector R size overflow\n");

        return LIS_OUT_OF_MEMORY;
    }


    r_entries *=
        (size_t)mode_count;


    counts = NULL;
    member_offset = NULL;
    cursor = NULL;
    members = NULL;

    aggregate_rank = NULL;
    aggregate_offset = NULL;
    row_ptr = NULL;

    q_all = NULL;
    r_all = NULL;
    work = NULL;
    prolong = NULL;
    coarse_basis = NULL;

    err = LIS_SUCCESS;


    counts =
        (LIS_INT *)lis_malloc(
            (size_t)aggregate_count
            *
            sizeof(LIS_INT),
            "lis_precon_mfamg_build_multivector_basis::counts");

    member_offset =
        (LIS_INT *)lis_malloc(
            (size_t)(aggregate_count+1)
            *
            sizeof(LIS_INT),
            "lis_precon_mfamg_build_multivector_basis::member_offset");

    cursor =
        (LIS_INT *)lis_malloc(
            (size_t)aggregate_count
            *
            sizeof(LIS_INT),
            "lis_precon_mfamg_build_multivector_basis::cursor");

    members =
        (LIS_INT *)lis_malloc(
            (size_t)n
            *
            sizeof(LIS_INT),
            "lis_precon_mfamg_build_multivector_basis::members");

    aggregate_rank =
        (LIS_INT *)lis_malloc(
            (size_t)aggregate_count
            *
            sizeof(LIS_INT),
            "lis_precon_mfamg_build_multivector_basis::aggregate_rank");

    aggregate_offset =
        (LIS_INT *)lis_malloc(
            (size_t)(aggregate_count+1)
            *
            sizeof(LIS_INT),
            "lis_precon_mfamg_build_multivector_basis::aggregate_offset");

    row_ptr =
        (LIS_INT *)lis_malloc(
            (size_t)(n+1)
            *
            sizeof(LIS_INT),
            "lis_precon_mfamg_build_multivector_basis::row_ptr");

    q_all =
        (LIS_SCALAR *)lis_malloc(
            q_entries
            *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_build_multivector_basis::q_all");

    r_all =
        (LIS_SCALAR *)lis_malloc(
            r_entries
            *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_build_multivector_basis::r_all");

    work =
        (LIS_SCALAR *)lis_malloc(
            (size_t)n
            *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_build_multivector_basis::work");


    if(
        counts==NULL
        ||
        member_offset==NULL
        ||
        cursor==NULL
        ||
        members==NULL
        ||
        aggregate_rank==NULL
        ||
        aggregate_offset==NULL
        ||
        row_ptr==NULL
        ||
        q_all==NULL
        ||
        r_all==NULL
        ||
        work==NULL)
    {
        LIS_SETERR_MEM(
            q_entries
            *
            sizeof(LIS_SCALAR));

        err = LIS_OUT_OF_MEMORY;
        goto fail;
    }


    memset(
        counts,
        0,
        (size_t)aggregate_count
        *
        sizeof(LIS_INT));

    memset(
        q_all,
        0,
        q_entries
        *
        sizeof(LIS_SCALAR));

    memset(
        r_all,
        0,
        r_entries
        *
        sizeof(LIS_SCALAR));


    for(i=0;i<n;i++)
    {
        a =
            aggregate[i];


        if(
            a<0
            ||
            a>=aggregate_count)
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG multivector basis contains an invalid aggregate index\n");

            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        counts[a]++;
    }


    member_offset[0] =
        0;


    for(a=0;a<aggregate_count;a++)
    {
        next_size =
            (size_t)member_offset[a]
            +
            (size_t)counts[a];


        if(
            (size_t)((LIS_INT)next_size)
            !=
            next_size)
        {
            LIS_SETERR(
                LIS_OUT_OF_MEMORY,
                "MFAMG multivector aggregate membership overflow\n");

            err = LIS_OUT_OF_MEMORY;
            goto fail;
        }


        member_offset[a+1] =
            (LIS_INT)next_size;

        cursor[a] =
            member_offset[a];
    }


    if(
        member_offset[
            aggregate_count]
        !=
        n)
    {
        err = LIS_ERR_ILL_ARG;
        goto fail;
    }


    for(i=0;i<n;i++)
    {
        a =
            aggregate[i];

        members[
            cursor[a]]
            =
            i;

        cursor[a]++;
    }


    for(a=0;a<aggregate_count;a++)
    {
        if(counts[a]<=0)
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG multivector basis contains an empty aggregate\n");

            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        count =
            counts[a];

        max_norm =
            0.0;


        /*
         * Determine the local scale used for relative rank
         * detection.
         */
        for(m=0;m<mode_count;m++)
        {
            norm2 =
                0.0;


            for(t=0;t<count;t++)
            {
                node =
                    members[
                        member_offset[a]
                        +
                        t];

                mag =
                    fabs(
                        basis[m]
                        ->
                        value[node]);


                if(
                    mag!=mag
                    ||
                    mag>LIS_SCALAR_MAX)
                {
                    LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "MFAMG multivector basis contains non-finite values\n");

                    err = LIS_ERR_ILL_ARG;
                    goto fail;
                }


                norm2 +=
                    mag
                    *
                    mag;


                if(
                    norm2!=norm2
                    ||
                    norm2>LIS_SCALAR_MAX)
                {
                    LIS_SETERR(
                        LIS_BREAKDOWN,
                        "MFAMG multivector basis norm overflowed\n");

                    err = LIS_BREAKDOWN;
                    goto fail;
                }
            }


            norm =
                sqrt(
                    norm2);


            if(norm>max_norm)
            {
                max_norm =
                    norm;
            }
        }


        if(max_norm==0.0)
        {
            /*
             * A locally constrained aggregate may see every
             * supplied near-nullspace mode vanish.  It has zero
             * numerical rank and contributes no coarse unknowns.
             *
             * A globally vanishing basis remains invalid because
             * the total coarse dimension is checked below.
             */
            aggregate_rank[a] =
                0;

            continue;
        }


        threshold =
            LIS_MFAMG_RANK_TOL
            *
            max_norm;


        if(
            threshold!=threshold
            ||
            threshold>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG multivector rank threshold overflowed\n");

            err = LIS_BREAKDOWN;
            goto fail;
        }


        rank =
            0;


        for(m=0;m<mode_count;m++)
        {
            for(t=0;t<count;t++)
            {
                node =
                    members[
                        member_offset[a]
                        +
                        t];

                work[t] =
                    basis[m]
                    ->
                    value[node];
            }


            /*
             * Two-pass modified Gram-Schmidt reduces the
             * loss of orthogonality for nearly dependent
             * local modes.
             */
            for(pass=0;pass<2;pass++)
            {
                for(q=0;q<rank;q++)
                {
                    dot =
                        (LIS_SCALAR)0.0;


                    for(t=0;t<count;t++)
                    {
                        node =
                            members[
                                member_offset[a]
                                +
                                t];

                        dot +=
                            conj(
                                q_all[
                                    (size_t)node
                                    *
                                    (size_t)mode_count
                                    +
                                    (size_t)q])
                            *
                            work[t];
                    }


                    r_all[
                        (
                            (size_t)a
                            *
                            (size_t)mode_count
                            +
                            (size_t)q
                        )
                        *
                        (size_t)mode_count
                        +
                        (size_t)m]
                        +=
                        dot;


                    for(t=0;t<count;t++)
                    {
                        node =
                            members[
                                member_offset[a]
                                +
                                t];

                        work[t] -=
                            q_all[
                                (size_t)node
                                *
                                (size_t)mode_count
                                +
                                (size_t)q]
                            *
                            dot;
                    }
                }
            }


            norm2 =
                0.0;


            for(t=0;t<count;t++)
            {
                mag =
                    fabs(
                        work[t]);


                if(
                    mag!=mag
                    ||
                    mag>LIS_SCALAR_MAX)
                {
                    LIS_SETERR(
                        LIS_BREAKDOWN,
                        "MFAMG multivector orthogonalization produced non-finite values\n");

                    err = LIS_BREAKDOWN;
                    goto fail;
                }


                norm2 +=
                    mag
                    *
                    mag;


                if(
                    norm2!=norm2
                    ||
                    norm2>LIS_SCALAR_MAX)
                {
                    LIS_SETERR(
                        LIS_BREAKDOWN,
                        "MFAMG multivector residual norm overflowed\n");

                    err = LIS_BREAKDOWN;
                    goto fail;
                }
            }


            norm =
                sqrt(
                    norm2);


            if(norm<=threshold)
                continue;


            if(
                rank>=count
                ||
                rank>=mode_count)
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "MFAMG multivector local rank is inconsistent\n");

                err = LIS_ERR_ILL_ARG;
                goto fail;
            }


            r_all[
                (
                    (size_t)a
                    *
                    (size_t)mode_count
                    +
                    (size_t)rank
                )
                *
                (size_t)mode_count
                +
                (size_t)m]
                =
                (LIS_SCALAR)norm;


            for(t=0;t<count;t++)
            {
                node =
                    members[
                        member_offset[a]
                        +
                        t];

                q_all[
                    (size_t)node
                    *
                    (size_t)mode_count
                    +
                    (size_t)rank]
                    =
                    work[t]
                    /
                    (LIS_SCALAR)norm;
            }


            rank++;
        }




        aggregate_rank[a] =
            rank;
    }


    prefix =
        0;

    aggregate_offset[0] =
        0;


    for(a=0;a<aggregate_count;a++)
    {
        prefix +=
            (size_t)aggregate_rank[a];


        if(
            (size_t)((LIS_INT)prefix)
            !=
            prefix)
        {
            LIS_SETERR(
                LIS_OUT_OF_MEMORY,
                "MFAMG multivector coarse dimension overflow\n");

            err = LIS_OUT_OF_MEMORY;
            goto fail;
        }


        aggregate_offset[a+1] =
            (LIS_INT)prefix;
    }


    coarse_n =
        aggregate_offset[
            aggregate_count];


    if(coarse_n<=0)
    {
        err = LIS_BREAKDOWN;
        goto fail;
    }


    prefix =
        0;

    row_ptr[0] =
        0;


    for(i=0;i<n;i++)
    {
        a =
            aggregate[i];

        prefix +=
            (size_t)aggregate_rank[a];


        if(
            (size_t)((LIS_INT)prefix)
            !=
            prefix)
        {
            LIS_SETERR(
                LIS_OUT_OF_MEMORY,
                "MFAMG multivector prolongation size overflow\n");

            err = LIS_OUT_OF_MEMORY;
            goto fail;
        }


        row_ptr[i+1] =
            (LIS_INT)prefix;
    }


    prolong_entries =
        prefix;


    if(
        (size_t)coarse_n
        >
        ((size_t)-1)
        /
        (size_t)mode_count)
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG multivector coarse basis size overflow\n");

        err = LIS_OUT_OF_MEMORY;
        goto fail;
    }


    coarse_entries =
        (size_t)coarse_n
        *
        (size_t)mode_count;


    if(
        prolong_entries
        >
        ((size_t)-1)
        /
        sizeof(LIS_SCALAR)
        ||
        coarse_entries
        >
        ((size_t)-1)
        /
        sizeof(LIS_SCALAR))
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG multivector output allocation overflow\n");

        err = LIS_OUT_OF_MEMORY;
        goto fail;
    }


    prolong =
        (LIS_SCALAR *)lis_malloc(
            prolong_entries
            *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_build_multivector_basis::prolong");

    coarse_basis =
        (LIS_SCALAR *)lis_malloc(
            coarse_entries
            *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_build_multivector_basis::coarse_basis");


    if(
        prolong==NULL
        ||
        coarse_basis==NULL)
    {
        LIS_SETERR_MEM(
            coarse_entries
            *
            sizeof(LIS_SCALAR));

        err = LIS_OUT_OF_MEMORY;
        goto fail;
    }


    memset(
        coarse_basis,
        0,
        coarse_entries
        *
        sizeof(LIS_SCALAR));


    for(i=0;i<n;i++)
    {
        a =
            aggregate[i];


        for(
            q=0;
            q<aggregate_rank[a];
            q++)
        {
            prolong[
                (size_t)row_ptr[i]
                +
                (size_t)q]
                =
                q_all[
                    (size_t)i
                    *
                    (size_t)mode_count
                    +
                    (size_t)q];
        }
    }


    for(a=0;a<aggregate_count;a++)
    {
        for(
            q=0;
            q<aggregate_rank[a];
            q++)
        {
            for(m=0;m<mode_count;m++)
            {
                coarse_basis[
                    (
                        (size_t)aggregate_offset[a]
                        +
                        (size_t)q
                    )
                    *
                    (size_t)mode_count
                    +
                    (size_t)m]
                    =
                    r_all[
                        (
                            (size_t)a
                            *
                            (size_t)mode_count
                            +
                            (size_t)q
                        )
                        *
                        (size_t)mode_count
                        +
                        (size_t)m];
            }
        }
    }


    lis_free(
        work);

    lis_free(
        r_all);

    lis_free(
        q_all);

    lis_free(
        members);

    lis_free(
        cursor);

    lis_free(
        member_offset);

    lis_free(
        counts);


    *aggregate_rank_out =
        aggregate_rank;

    *aggregate_offset_out =
        aggregate_offset;

    *row_ptr_out =
        row_ptr;

    *prolong_out =
        prolong;

    *coarse_basis_out =
        coarse_basis;

    *coarse_n_out =
        coarse_n;


    return LIS_SUCCESS;


fail:

    if(coarse_basis)
        lis_free(coarse_basis);

    if(prolong)
        lis_free(prolong);

    if(work)
        lis_free(work);

    if(r_all)
        lis_free(r_all);

    if(q_all)
        lis_free(q_all);

    if(row_ptr)
        lis_free(row_ptr);

    if(aggregate_offset)
        lis_free(aggregate_offset);

    if(aggregate_rank)
        lis_free(aggregate_rank);

    if(members)
        lis_free(members);

    if(cursor)
        lis_free(cursor);

    if(member_offset)
        lis_free(member_offset);

    if(counts)
        lis_free(counts);


    return err;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_build_prolongation"
static LIS_INT
lis_mfamg_build_prolongation(
    LIS_VECTOR z,
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_REAL *norm2;

    LIS_INT i;
    LIS_INT a;

    LIS_REAL mag;
    LIS_REAL denom;


    norm2 =
        (LIS_REAL *)lis_malloc(
            (size_t)level->nc
                *
            sizeof(LIS_REAL),
            "lis_mfamg_build_prolongation::norm2");

    if(norm2==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)level->nc
                *
            sizeof(LIS_REAL));

        return LIS_OUT_OF_MEMORY;
    }


    memset(
        norm2,
        0,
        (size_t)level->nc
            *
        sizeof(LIS_REAL));


    for(i=0;i<level->n;i++)
    {
        a =
            level->aggregate[i];

        mag =
            fabs(
                z->value[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            lis_free(
                norm2);

            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG near-nullspace vector contains non-finite values\n");

            return LIS_ERR_ILL_ARG;
        }

        norm2[a] +=
            mag
            *
            mag;

        if(
            norm2[a]!=norm2[a]
            ||
            norm2[a]>LIS_SCALAR_MAX)
        {
            lis_free(
                norm2);

            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG near-nullspace aggregate norm overflowed\n");

            return LIS_BREAKDOWN;
        }
    }


    for(a=0;a<level->nc;a++)
    {
        if(norm2[a]==0.0)
        {
            lis_free(
                norm2);

            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG near-nullspace vector vanishes on an aggregate\n");

            return LIS_ERR_ILL_ARG;
        }
    }


    for(i=0;i<level->n;i++)
    {
        a =
            level->aggregate[i];

        /*
         * LIS math macros map sqrt/fabs correctly for both
         * real and complex scalar builds.  norm2 is real and
         * non-negative, so fabs(sqrt(...)) yields LIS_REAL.
         */
        denom =
            fabs(
                sqrt(
                    (LIS_SCALAR)
                    norm2[a]));

        level->prolong[i] =
            z->value[i]
            /
            denom;
    }


    lis_free(
        norm2);

    return LIS_SUCCESS;
}


static LIS_INT
lis_mfamg_build_distributed_multivector_prolongation(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_SCALAR *global_prolong;

    LIS_REAL mag;

    LIS_INT local_width;
    LIS_INT width;

    LIS_INT i;
    LIS_INT a;
    LIS_INT q;
    LIS_INT rank;
    LIS_INT row_start;
    LIS_INT err;

    size_t entries;


    if(
        A==NULL
        ||
        level==NULL
        ||
        level->n!=A->n
        ||
        level->aggregate_count<=0
        ||
        level->nc<=0
        ||
        level->global_nc<=0
        ||
        level->aggregate==NULL
        ||
        level->aggregate_rank==NULL
        ||
        level->aggregate_offset==NULL
        ||
        level->prolong_row_ptr==NULL
        ||
        level->prolong==NULL
        ||
        level->global_aggregate==NULL
        ||
        level->global_aggregate_rank==NULL
        ||
        level->coarse_ranges==NULL
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs<=1
        ||
        A->np<A->n
        ||
        level->global_prolong!=NULL
        ||
        level->global_prolong_width!=0)
    {
        return LIS_ERR_ILL_ARG;
    }


    local_width = 0;


    for(a=0;a<level->aggregate_count;a++)
    {
        rank =
            level->aggregate_rank[a];

        if(rank<0)
            return LIS_ERR_ILL_ARG;

        if(rank>local_width)
            local_width = rank;
    }


    width =
        local_width;


#ifdef USE_MPI

    if(
        MPI_Allreduce(
            &local_width,
            &width,
            1,
            LIS_MPI_INT,
            MPI_MAX,
            A->comm)
        !=
        MPI_SUCCESS)
    {
        return LIS_ERR_ILL_ARG;
    }

#endif


    if(width<=0)
        return LIS_ERR_ILL_ARG;


    if(
        (size_t)A->np
        >
        ((size_t)-1)
        /
        (size_t)width)
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG distributed multivector prolongation size overflow\n");

        return LIS_OUT_OF_MEMORY;
    }


    entries =
        (size_t)A->np
        *
        (size_t)width;


    if(
        entries
        >
        ((size_t)-1)
        /
        sizeof(LIS_SCALAR))
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG distributed multivector prolongation allocation overflow\n");

        return LIS_OUT_OF_MEMORY;
    }


    global_prolong =
        (LIS_SCALAR *)lis_malloc(
            entries
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_build_distributed_multivector_prolongation::global_prolong");


    if(global_prolong==NULL)
    {
        LIS_SETERR_MEM(
            entries
            *
            sizeof(LIS_SCALAR));

        return LIS_OUT_OF_MEMORY;
    }


    memset(
        global_prolong,
        0,
        entries
        *
        sizeof(LIS_SCALAR));


    for(i=0;i<A->n;i++)
    {
        a =
            level->aggregate[i];


        if(
            a<0
            ||
            a>=level->aggregate_count)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        rank =
            level->aggregate_rank[a];

        row_start =
            level->prolong_row_ptr[i];


        if(
            rank<0
            ||
            rank>width
            ||
            row_start<0
            ||
            level->prolong_row_ptr[i+1]
                -
            row_start
                !=
            rank)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        if(
            level->global_aggregate_rank[i]
            !=
            rank)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        for(q=0;q<rank;q++)
        {
            global_prolong[
                (size_t)q
                *
                (size_t)A->np
                +
                (size_t)i]
                =
            level->prolong[
                row_start
                +
                q];
        }
    }


    if(A->np>A->n)
    {
        if(A->commtable==NULL)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        for(q=0;q<width;q++)
        {
            err =
                lis_mfamg_exchange_scalar_halo(
                    A,
                    global_prolong
                    +
                    (size_t)q
                    *
                    (size_t)A->np);

            if(err)
                goto fail;
        }
    }


    for(i=0;i<A->np;i++)
    {
        rank =
            level->global_aggregate_rank[i];


        if(
            rank<0
            ||
            rank>width
            ||
            level->global_aggregate[i]<0
            ||
            level->global_aggregate[i]>level->global_nc
            ||
            level->global_aggregate[i]+rank
                >
            level->global_nc)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        for(q=0;q<rank;q++)
        {
            mag =
                fabs(
                    global_prolong[
                        (size_t)q
                        *
                        (size_t)A->np
                        +
                        (size_t)i]);


            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "MFAMG distributed multivector prolongation contains non-finite values\n");

                err = LIS_ERR_ILL_ARG;
                goto fail;
            }
        }


        for(q=rank;q<width;q++)
        {
            if(
                global_prolong[
                    (size_t)q
                    *
                    (size_t)A->np
                    +
                    (size_t)i]
                !=
                (LIS_SCALAR)0.0)
            {
                err = LIS_ERR_ILL_ARG;
                goto fail;
            }
        }
    }


    level->global_prolong =
        global_prolong;

    level->global_prolong_width =
        width;


    return LIS_SUCCESS;


fail:

    lis_free(
        global_prolong);


    return err;
}


static LIS_INT
lis_mfamg_build_distributed_prolongation(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_SCALAR *global_prolong;

    LIS_REAL mag;

    LIS_INT i;
    LIS_INT err;


    if(
        level!=NULL
        &&
        level->prolong_row_ptr!=NULL)
    {
        return
            lis_mfamg_build_distributed_multivector_prolongation(
                A,
                level);
    }


    if(
        A==NULL
        ||
        level==NULL
        ||
        level->n!=A->n
        ||
        level->prolong==NULL
        ||
        level->aggregate==NULL
        ||
        level->global_aggregate==NULL
        ||
        level->coarse_ranges==NULL
        ||
        level->nc<=0
        ||
        level->global_nc<=0
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs<=1
        ||
        A->np<A->n
        ||
        level->global_prolong!=NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    global_prolong =
        (LIS_SCALAR *)lis_malloc(
            (size_t)A->np
                *
            sizeof(LIS_SCALAR),
            "lis_mfamg_build_distributed_prolongation::global_prolong");


    if(global_prolong==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)A->np
                *
            sizeof(LIS_SCALAR));

        return LIS_OUT_OF_MEMORY;
    }


    for(i=0;i<A->np;i++)
    {
        global_prolong[i] =
            (LIS_SCALAR)0.0;
    }


    for(i=0;i<A->n;i++)
    {
        global_prolong[i] =
            level->prolong[i];
    }


    if(A->np>A->n)
    {
        if(A->commtable==NULL)
        {
            lis_free(
                global_prolong);

            return LIS_ERR_ILL_ARG;
        }


        err =
            lis_mfamg_exchange_scalar_halo(
                A,
                global_prolong);

        if(err)
        {
            lis_free(
                global_prolong);

            return err;
        }
    }


    for(i=0;i<A->np;i++)
    {
        if(
            level->global_aggregate[i]<0
            ||
            level->global_aggregate[i]>=level->global_nc)
        {
            lis_free(
                global_prolong);

            return LIS_ERR_ILL_ARG;
        }


        mag =
            fabs(
                global_prolong[i]);


        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            lis_free(
                global_prolong);

            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG distributed prolongation contains non-finite values\n");

            return LIS_ERR_ILL_ARG;
        }
    }


    level->global_prolong =
        global_prolong;


    return LIS_SUCCESS;
}






static int
lis_mfamg_cmp_int(
    const void *a,
    const void *b)
{
    LIS_INT x;
    LIS_INT y;

    x = *(const LIS_INT *)a;
    y = *(const LIS_INT *)b;

    return (x>y) - (x<y);
}


/*
 * Sparse Galerkin product Ac = P^H A P for the serial hierarchy.
 *
 * Replaces the former dense nc x nc accumulation (O(nc^2) memory) by a
 * row-by-row sparse accumulation, one aggregate at a time.  Fine rows are
 * visited in increasing order inside each aggregate, so every coarse entry
 * is summed in exactly the same order as the dense code did and the result
 * is bit-identical.  Exact zeros are dropped, as before.
 *
 * Workspace is O(R*nc + n) with R the largest aggregate rank (1 for a
 * single near-nullspace vector); the output is O(nnz(Ac)).
 */
#undef __FUNC__
#define __FUNC__ "lis_mfamg_create_coarse_csr"
static LIS_INT
lis_mfamg_create_coarse_csr(
    LIS_Comm comm,
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level,
    LIS_INT multivector,
    LIS_INT max_row_nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals,
    struct LIS_MFAMG_LEVEL *out)
{
    LIS_INT nc,nagg,n,R;
    LIS_INT a,i,j,k,q,qi,qj,r,ri,rj,aj,cj,row_i,row_j,base,nnz,err;
    LIS_INT *agg_ptr,*agg_rows,*mark,*tl,*cnt,*ptr,*index,*pp,*pi_,*pv_;
    LIS_SCALAR *acc,*value,*vv,pi,pj,v;
    LIS_INT cap,nzc,*tmp_i;
    LIS_SCALAR *tmp_v;
    size_t ws;
    LIS_MATRIX M;

    if(
        A==NULL || level==NULL || out==NULL || cols==NULL || vals==NULL
        || A->nprocs!=1 || A->is!=0 || A->gn!=A->n
        || out->A!=NULL || level->aggregate==NULL || level->prolong==NULL
        || level->n<=0 || level->nc<=0)
    {
        return LIS_ERR_ILL_ARG;
    }

    n = level->n;
    nc = level->nc;
    nagg = multivector ? level->aggregate_count : nc;
    R = 1;

    if(multivector)
    {
        if(
            level->aggregate_rank==NULL || level->aggregate_offset==NULL
            || level->prolong_row_ptr==NULL || nagg<=0
            || level->aggregate_offset[0]!=0
            || level->aggregate_offset[nagg]!=nc
            || level->prolong_row_ptr[0]!=0)
        {
            return LIS_ERR_ILL_ARG;
        }
        for(a=0;a<nagg;a++)
        {
            if(
                level->aggregate_rank[a]<0
                || level->aggregate_offset[a+1]-level->aggregate_offset[a]
                   !=level->aggregate_rank[a])
            {
                return LIS_ERR_ILL_ARG;
            }
            if(level->aggregate_rank[a]>R) R = level->aggregate_rank[a];
        }
    }

    agg_ptr = agg_rows = mark = tl = cnt = ptr = index = NULL;
    acc = value = NULL;
    err = LIS_SUCCESS;

    if(
        n<=0 || nc<=0 || nagg<=0 || R<=0
        || (size_t)n > ((size_t)-1)/sizeof(LIS_INT)
        || (size_t)R > ((size_t)-1)/sizeof(LIS_INT)
        || (size_t)nagg > (((size_t)-1)/sizeof(LIS_INT))-1
        || (size_t)nc > (((size_t)-1)/sizeof(LIS_INT))-1)
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG sparse Galerkin index workspace overflow\n");
        return LIS_OUT_OF_MEMORY;
    }

    if((size_t)R > ((size_t)-1)/(size_t)nc)
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG sparse Galerkin workspace dimension overflow\n");
        return LIS_OUT_OF_MEMORY;
    }

    ws = (size_t)R*(size_t)nc;

    if(
        ws > ((size_t)-1)/sizeof(LIS_INT)
        || ws > ((size_t)-1)/sizeof(LIS_SCALAR))
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG sparse Galerkin workspace allocation overflow\n");
        return LIS_OUT_OF_MEMORY;
    }

    agg_ptr = (LIS_INT *)lis_malloc(((size_t)nagg+1)*sizeof(LIS_INT),"lis_mfamg_create_coarse_csr::agg_ptr");
    agg_rows = (LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),"lis_mfamg_create_coarse_csr::agg_rows");
    mark = (LIS_INT *)lis_malloc(ws*sizeof(LIS_INT),"lis_mfamg_create_coarse_csr::mark");
    tl = (LIS_INT *)lis_malloc(ws*sizeof(LIS_INT),"lis_mfamg_create_coarse_csr::tl");
    cnt = (LIS_INT *)lis_malloc((size_t)R*sizeof(LIS_INT),"lis_mfamg_create_coarse_csr::cnt");
    ptr = (LIS_INT *)lis_malloc(((size_t)nc+1)*sizeof(LIS_INT),"lis_mfamg_create_coarse_csr::ptr");
    acc = (LIS_SCALAR *)lis_malloc(ws*sizeof(LIS_SCALAR),"lis_mfamg_create_coarse_csr::acc");

    cap = nc;
    if(cap<256) cap = 256;
    if(cap>65536) cap = 65536;
    index = (LIS_INT *)lis_malloc((size_t)cap*sizeof(LIS_INT),"lis_mfamg_create_coarse_csr::index");
    value = (LIS_SCALAR *)lis_malloc((size_t)cap*sizeof(LIS_SCALAR),"lis_mfamg_create_coarse_csr::value");

    if(
        agg_ptr==NULL || agg_rows==NULL || mark==NULL || tl==NULL
        || cnt==NULL || ptr==NULL || acc==NULL || index==NULL || value==NULL)
    {
        LIS_SETERR_MEM(ws*sizeof(LIS_SCALAR));
        err = LIS_OUT_OF_MEMORY;
        goto done;
    }

    for(a=0;a<=nagg;a++) agg_ptr[a] = 0;
    for(i=0;i<n;i++)
    {
        a = level->aggregate[i];
        if(a<0 || a>=nagg) { err = LIS_ERR_ILL_ARG; goto done; }
        agg_ptr[a+1]++;
    }
    for(a=0;a<nagg;a++) agg_ptr[a+1] += agg_ptr[a];
    {
        LIS_INT *fillpos;
        fillpos = (LIS_INT *)lis_malloc((size_t)nagg*sizeof(LIS_INT),"lis_mfamg_create_coarse_csr::fillpos");
        if(fillpos==NULL) { LIS_SETERR_MEM((size_t)nagg*sizeof(LIS_INT)); err = LIS_OUT_OF_MEMORY; goto done; }
        for(a=0;a<nagg;a++) fillpos[a] = agg_ptr[a];
        for(i=0;i<n;i++) agg_rows[fillpos[level->aggregate[i]]++] = i;
        lis_free(fillpos);
    }
    memset(mark,0,ws*sizeof(LIS_INT));

    nzc = 0;
    ptr[0] = 0;

    for(a=0;a<nagg;a++)
    {
        ri = multivector ? level->aggregate_rank[a] : 1;
        for(q=0;q<ri;q++) cnt[q] = 0;

        for(r=agg_ptr[a];r<agg_ptr[a+1];r++)
        {
            i = agg_rows[r];
            row_i = 0;
            if(multivector)
            {
                row_i = level->prolong_row_ptr[i];
                if(row_i<0 || level->prolong_row_ptr[i+1]-row_i!=ri)
                { err = LIS_ERR_ILL_ARG; goto done; }
            }
            nnz = max_row_nnz;
            if(lis_mfamg_row_fill(A,i,max_row_nnz,&nnz,cols,vals)!=LIS_SUCCESS)
            { err = LIS_ERR_ILL_ARG; goto done; }

            for(k=0;k<nnz;k++)
            {
                j = cols[k]-A->is;
                if(j<0 || j>=n)
                {
                    LIS_SETERR(LIS_ERR_ILL_ARG,"MFAMG serial row access returned a nonlocal column\n");
                    err = LIS_ERR_ILL_ARG;
                    goto done;
                }
                aj = level->aggregate[j];
                if(multivector)
                {
                    rj = level->aggregate_rank[aj];
                    row_j = level->prolong_row_ptr[j];
                    if(row_j<0 || level->prolong_row_ptr[j+1]-row_j!=rj)
                    { err = LIS_ERR_ILL_ARG; goto done; }
                    base = level->aggregate_offset[aj];
                }
                else
                {
                    rj = 1;
                    row_j = 0;
                    base = aj;
                }

                for(qi=0;qi<ri;qi++)
                {
                    pi = multivector ? level->prolong[row_i+qi] : level->prolong[i];
                    for(qj=0;qj<rj;qj++)
                    {
                        size_t idx;
                        cj = base+qj;
                        pj = multivector ? level->prolong[row_j+qj] : level->prolong[j];
                        idx = (size_t)qi*(size_t)nc+(size_t)cj;
                        if(mark[idx]!=a+1)
                        {
                            mark[idx] = a+1;
                            acc[idx] = (LIS_SCALAR)0.0;
                            tl[(size_t)qi*(size_t)nc+(size_t)cnt[qi]] = cj;
                            cnt[qi]++;
                        }
                        acc[idx] += conj(pi)*vals[k]*pj;
                    }
                }
            }
        }

        for(qi=0;qi<ri;qi++)
        {
            LIS_INT *list;
            list = tl+(size_t)qi*(size_t)nc;
            qsort(list,(size_t)cnt[qi],sizeof(LIS_INT),lis_mfamg_cmp_int);
            for(q=0;q<cnt[qi];q++)
            {
                v = acc[(size_t)qi*(size_t)nc+(size_t)list[q]];
                if(v==(LIS_SCALAR)0.0) continue;
                if(nzc==cap)
                {
                    size_t new_cap_size;
                    LIS_INT new_cap;

                    if(cap<=0 || (size_t)cap > ((size_t)-1)/2)
                    {
                        LIS_SETERR(
                            LIS_OUT_OF_MEMORY,
                            "MFAMG coarse CSR capacity overflow\n");
                        err = LIS_OUT_OF_MEMORY;
                        goto done;
                    }

                    new_cap_size = (size_t)cap*2;

                    if(
                        new_cap_size > ((size_t)-1)/sizeof(LIS_INT)
                        || new_cap_size > ((size_t)-1)/sizeof(LIS_SCALAR))
                    {
                        LIS_SETERR(
                            LIS_OUT_OF_MEMORY,
                            "MFAMG coarse CSR allocation overflow\n");
                        err = LIS_OUT_OF_MEMORY;
                        goto done;
                    }

                    new_cap = (LIS_INT)new_cap_size;

                    if((size_t)new_cap != new_cap_size)
                    {
                        LIS_SETERR(
                            LIS_OUT_OF_MEMORY,
                            "MFAMG coarse CSR index range overflow\n");
                        err = LIS_OUT_OF_MEMORY;
                        goto done;
                    }

                    tmp_i = (LIS_INT *)lis_malloc(
                        new_cap_size*sizeof(LIS_INT),
                        "lis_mfamg_create_coarse_csr::index");
                    tmp_v = (LIS_SCALAR *)lis_malloc(
                        new_cap_size*sizeof(LIS_SCALAR),
                        "lis_mfamg_create_coarse_csr::value");

                    if(tmp_i==NULL || tmp_v==NULL)
                    {
                        if(tmp_i) lis_free(tmp_i);
                        if(tmp_v) lis_free(tmp_v);
                        LIS_SETERR_MEM(new_cap_size*sizeof(LIS_SCALAR));
                        err = LIS_OUT_OF_MEMORY;
                        goto done;
                    }

                    memcpy(tmp_i,index,(size_t)nzc*sizeof(LIS_INT));
                    memcpy(tmp_v,value,(size_t)nzc*sizeof(LIS_SCALAR));
                    lis_free(index);
                    lis_free(value);
                    index = tmp_i;
                    value = tmp_v;
                    cap = new_cap;
                }
                index[nzc] = list[q];
                value[nzc] = v;
                nzc++;
            }
            ptr[(multivector ? level->aggregate_offset[a] : a)+qi+1] = nzc;
        }
    }

    if(nzc==0)
    {
        LIS_SETERR(LIS_BREAKDOWN,"MFAMG coarse CSR matrix has no entries\n");
        err = LIS_BREAKDOWN;
        goto done;
    }

    M = NULL;
    err = lis_matrix_create(comm,&M);
    if(err) goto done;
    err = lis_matrix_set_size(M,nc,nc);
    if(err) { lis_matrix_destroy(M); goto done; }
    pp = pi_ = NULL;
    vv = NULL;
    pv_ = NULL;
    err = lis_matrix_malloc_csr(nc,nzc,&pp,&pi_,&vv);
    if(err) { lis_matrix_destroy(M); goto done; }
    memcpy(pp,ptr,(size_t)(nc+1)*sizeof(LIS_INT));
    memcpy(pi_,index,(size_t)nzc*sizeof(LIS_INT));
    memcpy(vv,value,(size_t)nzc*sizeof(LIS_SCALAR));
    (void)pv_;
    err = lis_matrix_set_csr(nzc,pp,pi_,vv,M);
    if(err)
    {
        lis_free(pp); lis_free(pi_); lis_free(vv);
        lis_matrix_destroy(M);
        goto done;
    }
    err = lis_matrix_assemble(M);
    if(err) { lis_matrix_destroy(M); goto done; }
    if(M->matrix_type!=LIS_MATRIX_CSR || M->n!=nc || M->gn!=nc)
    {
        lis_matrix_destroy(M);
        LIS_SETERR(LIS_ERR_ILL_ARG,"MFAMG internal coarse CSR matrix is inconsistent\n");
        err = LIS_ERR_ILL_ARG;
        goto done;
    }

    out->A = M;
    out->owns_A = LIS_TRUE;
    out->n = nc;
    out->nc = 0;

done:
    if(agg_ptr) lis_free(agg_ptr);
    if(agg_rows) lis_free(agg_rows);
    if(mark) lis_free(mark);
    if(tl) lis_free(tl);
    if(cnt) lis_free(cnt);
    if(ptr) lis_free(ptr);
    if(acc) lis_free(acc);
    if(index) lis_free(index);
    if(value) lis_free(value);
    return err;
}


static LIS_INT
lis_mfamg_get_max_row_nnz(
    LIS_MATRIX A,
    LIS_INT n,
    LIS_INT *max_row_nnz);

/*
 * Internal serial sparse-Galerkin builder used by testmfamggalerkin.
 *
 * This symbol is intentionally not part of the public LIS API.
 * lis_mfamg_get_max_row_nnz() is forward-declared above because
 * its definition follows this test-only hook.
 */
#undef __FUNC__
#define __FUNC__ "lis_precon_mfamg_build_serial_coarse_for_test"
LIS_INT
lis_precon_mfamg_dense_coarse_solve_for_test(
    LIS_INT n,
    const LIS_SCALAR *matrix,
    const LIS_SCALAR *rhs,
    LIS_SCALAR *x,
    LIS_INT *nullity)
{
    LIS_SCALAR *lu;
    LIS_INT *pivots;

    LIS_INT i;
    LIS_INT has_null;
    LIS_INT err;


    if(
        n<=0
        ||
        matrix==NULL
        ||
        rhs==NULL
        ||
        x==NULL
        ||
        nullity==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    *nullity =
        0;

    lu =
        (LIS_SCALAR *)lis_malloc(
            (size_t)n
            *
            (size_t)n
            *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_dense_coarse_solve_for_test::lu");

    pivots =
        (LIS_INT *)lis_malloc(
            (size_t)n
            *
            sizeof(LIS_INT),
            "lis_precon_mfamg_dense_coarse_solve_for_test::pivots");


    if(
        lu==NULL
        ||
        pivots==NULL)
    {
        if(lu)
            lis_free(lu);

        if(pivots)
            lis_free(pivots);

        LIS_SETERR_MEM(
            (size_t)n
            *
            (
                (size_t)n
                *
                sizeof(LIS_SCALAR)
                +
                sizeof(LIS_INT)
            ));

        return LIS_OUT_OF_MEMORY;
    }


    memcpy(
        lu,
        matrix,
        (size_t)n
        *
        (size_t)n
        *
        sizeof(LIS_SCALAR));


    err =
        lis_mfamg_dense_factor(
            n,
            lu,
            pivots);

    if(err)
        goto cleanup;


    has_null =
        LIS_FALSE;


    for(i=0;i<n;i++)
    {
        if(
            fabs(
                lu[i*n+i])
            ==
            0.0)
        {
            has_null =
                LIS_TRUE;

            break;
        }
    }


    if(has_null)
    {
        err =
            lis_mfamg_dense_rank_revealing_solve(
                n,
                matrix,
                rhs,
                x,
                nullity);
    }
    else
    {
        err =
            lis_mfamg_dense_solve_values(
                n,
                matrix,
                lu,
                pivots,
                rhs,
                x);
    }


cleanup:

    lis_free(
        lu);

    lis_free(
        pivots);


    return err;
}



LIS_INT
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
    LIS_MATRIX *coarse_A)
{
    struct LIS_MFAMG_LEVEL level;
    struct LIS_MFAMG_LEVEL coarse_level;

    LIS_INT max_row_nnz;
    LIS_INT *cols;
    LIS_SCALAR *vals;

    LIS_INT err;

    if(
        A==NULL
        ||
        coarse_A==NULL
        ||
        aggregate==NULL
        ||
        prolong==NULL
        ||
        A->n<=0
        ||
        A->nprocs!=1
        ||
        nc<=0
        ||
        aggregate_count<=0)
    {
        return LIS_ERR_ILL_ARG;
    }

    if(
        multivector
        &&
        (
            aggregate_rank==NULL
            ||
            aggregate_offset==NULL
            ||
            prolong_row_ptr==NULL
        ))
    {
        return LIS_ERR_ILL_ARG;
    }

    *coarse_A = NULL;

    memset(
        &level,
        0,
        sizeof(level));

    memset(
        &coarse_level,
        0,
        sizeof(coarse_level));

    level.A =
        A;

    level.n =
        A->n;

    level.nc =
        nc;

    level.aggregate_count =
        aggregate_count;

    level.aggregate =
        aggregate;

    level.aggregate_rank =
        aggregate_rank;

    level.aggregate_offset =
        aggregate_offset;

    level.prolong_row_ptr =
        prolong_row_ptr;

    level.prolong =
        prolong;

    max_row_nnz = 0;
    cols = NULL;
    vals = NULL;

    err =
        lis_mfamg_get_max_row_nnz(
            A,
            A->n,
            &max_row_nnz);

    if(err)
        goto cleanup;

    cols =
        (LIS_INT *)lis_malloc(
            (size_t)max_row_nnz
            *
            sizeof(LIS_INT),
            "lis_precon_mfamg_build_serial_coarse_for_test::cols");

    vals =
        (LIS_SCALAR *)lis_malloc(
            (size_t)max_row_nnz
            *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_build_serial_coarse_for_test::vals");

    if(
        cols==NULL
        ||
        vals==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)max_row_nnz
            *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        err =
            LIS_OUT_OF_MEMORY;

        goto cleanup;
    }

    err =
        lis_mfamg_create_coarse_csr(
            A->comm,
            A,
            &level,
            multivector,
            max_row_nnz,
            cols,
            vals,
            &coarse_level);

    if(err)
        goto cleanup;

    if(
        coarse_level.A==NULL
        ||
        coarse_level.A->matrix_type!=LIS_MATRIX_CSR
        ||
        coarse_level.A->n!=nc
        ||
        coarse_level.A->gn!=nc)
    {
        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }

    *coarse_A =
        coarse_level.A;

    coarse_level.A =
        NULL;

    coarse_level.owns_A =
        LIS_FALSE;

    err =
        LIS_SUCCESS;

cleanup:

    if(coarse_level.A)
    {
        lis_matrix_destroy(
            coarse_level.A);
    }

    if(vals)
    {
        lis_free(
            vals);
    }

    if(cols)
    {
        lis_free(
            cols);
    }

    return err;
}



#undef __FUNC__
#define __FUNC__ "lis_mfamg_prepare_level_work"
static LIS_INT
lis_mfamg_prepare_level_work(
    struct LIS_MFAMG_LEVEL *level,
    LIS_INT block_size)
{
    if(
        level==NULL
        ||
        block_size<0
        ||
        block_size>LIS_MFAMG_BLOCK_MAX_SIZE)
    {
        return LIS_ERR_ILL_ARG;
    }


    /*
     * Zero is accepted internally for helper states that are
     * initialized with memset.  The public option itself begins
     * at one.
     */
    if(block_size>1)
    {
        level->smoother_type =
            LIS_MFAMG_SMOOTHER_BLOCK_JACOBI;

        level->block_size =
            block_size;
    }
    else
    {
        level->smoother_type =
            LIS_MFAMG_SMOOTHER_WEIGHTED_JACOBI;

        level->block_size =
            0;
    }


    return
        lis_mfamg_smoother_setup(
            level);
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_prepare_cycle_vectors"
static LIS_INT
lis_mfamg_prepare_cycle_vectors(
    struct LIS_MFAMG_STRUCT *state)
{
    struct LIS_MFAMG_LEVEL *level;

    LIS_INT level_id;
    LIS_INT err;


    if(
        state==NULL
        ||
        state->levels==NULL
        ||
        state->nlevels<2
        ||
        state->levels[0].A==NULL
        ||
        state->levels[0].n<=0
        ||
        state->alias_rhs!=NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    err =
        lis_vector_duplicate(
            state->levels[0].A,
            &state->alias_rhs);

    if(err)
        return err;


    for(
        level_id=1;
        level_id<state->nlevels;
        level_id++)
    {
        level =
            &state->levels[level_id];

        if(
            level->A==NULL
            ||
            level->n<=0
            ||
            level->rhs!=NULL
            ||
            level->x!=NULL)
        {
            return LIS_ERR_ILL_ARG;
        }


        err =
            lis_vector_duplicate(
                level->A,
                &level->rhs);

        if(err)
            return err;


        err =
            lis_vector_duplicate(
                level->A,
                &level->x);

        if(err)
            return err;
    }


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_get_max_row_nnz"
static LIS_INT
lis_mfamg_get_max_row_nnz(
    LIS_MATRIX A,
    LIS_INT n,
    LIS_INT *max_row_nnz)
{
    LIS_INT i;
    LIS_INT row_nnz;
    LIS_INT err;


    if(
        A==NULL
        ||
        n<=0
        ||
        max_row_nnz==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    *max_row_nnz = 0;


    for(i=0;i<n;i++)
    {
        row_nnz = 0;

        err =
            lis_mfamg_row_query(
                A,
                i,
                &row_nnz);

        if(err)
            return err;

        if(row_nnz<0)
            return LIS_ERR_ILL_ARG;

        if(row_nnz>*max_row_nnz)
        {
            *max_row_nnz =
                row_nnz;
        }
    }


    if(*max_row_nnz<=0)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG requires matrix rows with entries\n");

        return LIS_ERR_ILL_ARG;
    }


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_build_coarse_nullspace"
static LIS_INT
lis_mfamg_build_coarse_nullspace(
    LIS_VECTOR z,
    struct LIS_MFAMG_LEVEL *level,
    LIS_MATRIX coarse_A,
    LIS_VECTOR *coarse_z)
{
    LIS_VECTOR zc;

    LIS_INT i;
    LIS_INT a;
    LIS_INT err;

    LIS_REAL mag;


    if(
        z==NULL
        ||
        level==NULL
        ||
        coarse_A==NULL
        ||
        coarse_z==NULL
        ||
        level->n<=0
        ||
        level->nc<=0
        ||
        z->n!=level->n
        ||
        coarse_A->n!=level->nc)
    {
        return LIS_ERR_ILL_ARG;
    }


    *coarse_z = NULL;
    zc = NULL;


    err =
        lis_vector_duplicate(
            coarse_A,
            &zc);

    if(err)
        return err;


    memset(
        zc->value,
        0,
        (size_t)level->nc
            *
        sizeof(LIS_SCALAR));


    for(i=0;i<level->n;i++)
    {
        a =
            level->aggregate[i];

        if(
            a<0
            ||
            a>=level->nc)
        {
            lis_vector_destroy(zc);
            return LIS_ERR_ILL_ARG;
        }

        zc->value[a] +=
            conj(
                level->prolong[i])
            *
            z->value[i];
    }


    for(a=0;a<level->nc;a++)
    {
        mag =
            fabs(
                zc->value[a]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX
            ||
            mag==0.0)
        {
            lis_vector_destroy(zc);

            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse near-nullspace vector is invalid\n");

            return LIS_BREAKDOWN;
        }
    }


    *coarse_z =
        zc;

    return LIS_SUCCESS;
}


static LIS_INT
lis_mfamg_build_distributed_coarse(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level,
    LIS_MATRIX *coarse_A);


#undef __FUNC__
#define __FUNC__ "lis_mfamg_install_multivector_transfer"
static LIS_INT
lis_mfamg_install_multivector_transfer(
    LIS_INT mode_count,
    LIS_VECTOR basis[],
    struct LIS_MFAMG_LEVEL *level,
    LIS_SCALAR **coarse_basis_out,
    LIS_INT *coarsened)
{
    LIS_INT aggregate_count;
    LIS_INT coarse_n;
    LIS_INT err;

    LIS_INT *aggregate_rank;
    LIS_INT *aggregate_offset;
    LIS_INT *row_ptr;

    LIS_SCALAR *prolong;
    LIS_SCALAR *coarse_basis;


    if(
        mode_count<=1
        ||
        basis==NULL
        ||
        level==NULL
        ||
        level->aggregate==NULL
        ||
        level->n<=0
        ||
        level->nc<=0
        ||
        coarse_basis_out==NULL
        ||
        coarsened==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        level->aggregate_rank!=NULL
        ||
        level->aggregate_offset!=NULL
        ||
        level->prolong_row_ptr!=NULL
        ||
        level->prolong!=NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    *coarse_basis_out = NULL;
    *coarsened = LIS_FALSE;

    aggregate_count =
        level->nc;

    aggregate_rank = NULL;
    aggregate_offset = NULL;
    row_ptr = NULL;
    prolong = NULL;
    coarse_basis = NULL;
    coarse_n = 0;


    err =
        lis_precon_mfamg_build_multivector_basis(
            level->n,
            aggregate_count,
            level->aggregate,
            mode_count,
            basis,
            &aggregate_rank,
            &aggregate_offset,
            &row_ptr,
            &prolong,
            &coarse_basis,
            &coarse_n);

    if(err)
        return err;


    if(
        coarse_n<=0
        ||
        coarse_n>=level->n)
    {
        lis_free(
            coarse_basis);

        lis_free(
            prolong);

        lis_free(
            row_ptr);

        lis_free(
            aggregate_offset);

        lis_free(
            aggregate_rank);

        return LIS_SUCCESS;
    }


    level->aggregate_count =
        aggregate_count;

    level->aggregate_rank =
        aggregate_rank;

    level->aggregate_offset =
        aggregate_offset;

    level->prolong_row_ptr =
        row_ptr;

    level->prolong =
        prolong;

    level->nc =
        coarse_n;

    *coarse_basis_out =
        coarse_basis;

    *coarsened =
        LIS_TRUE;


    return LIS_SUCCESS;
}




static void
lis_mfamg_destroy_basis(
    LIS_INT mode_count,
    LIS_VECTOR *basis)
{
    LIS_INT m;


    if(basis==NULL)
        return;


    for(m=0;m<mode_count;m++)
    {
        if(basis[m])
        {
            lis_vector_destroy(
                basis[m]);
        }
    }


    lis_free(
        basis);
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_create_basis_from_array"
static LIS_INT
lis_mfamg_create_basis_from_array(
    LIS_MATRIX A,
    LIS_INT mode_count,
    LIS_SCALAR *values,
    LIS_VECTOR **basis_out)
{
    LIS_VECTOR *basis;

    LIS_INT i;
    LIS_INT m;
    LIS_INT err;

    LIS_REAL mag;


    if(
        A==NULL
        ||
        A->n<=0
        ||
        mode_count<=1
        ||
        values==NULL
        ||
        basis_out==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    *basis_out = NULL;


    basis =
        (LIS_VECTOR *)lis_malloc(
            (size_t)mode_count
                *
            sizeof(LIS_VECTOR),
            "lis_mfamg_create_basis_from_array::basis");

    if(basis==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)mode_count
                *
            sizeof(LIS_VECTOR));

        return LIS_OUT_OF_MEMORY;
    }


    memset(
        basis,
        0,
        (size_t)mode_count
            *
        sizeof(LIS_VECTOR));


    for(m=0;m<mode_count;m++)
    {
        err =
            lis_vector_duplicate(
                A,
                &basis[m]);

        if(err)
        {
            lis_mfamg_destroy_basis(
                mode_count,
                basis);

            return err;
        }


        for(i=0;i<A->n;i++)
        {
            basis[m]->value[i] =
                values[
                    (size_t)i
                        *
                    (size_t)mode_count
                        +
                    (size_t)m];

            mag =
                fabs(
                    basis[m]->value[i]);

            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                lis_mfamg_destroy_basis(
                    mode_count,
                    basis);

                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "MFAMG multivector coarse basis contains non-finite values\n");

                return LIS_ERR_ILL_ARG;
            }
        }
    }


    *basis_out =
        basis;


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_extend_hierarchy_multivector"
static LIS_INT
lis_mfamg_extend_hierarchy_multivector(
    struct LIS_MFAMG_STRUCT *state,
    LIS_INT mode_count,
    LIS_INT terminal_limit,
    LIS_SCALAR *fine_coarse_basis)
{
    struct LIS_MFAMG_LEVEL *level;
    struct LIS_MFAMG_LEVEL *next_level;

    LIS_VECTOR *current_basis;
    LIS_VECTOR *next_basis;

    LIS_MATRIX distributed_Ac;

    LIS_INT level_id;
    LIS_INT max_row_nnz;
    LIS_INT coarsened;
    LIS_INT aggregate_target;
    LIS_INT err;

    LIS_INT *cols;
    LIS_SCALAR *vals;
    LIS_SCALAR *next_coarse_basis;



    if(
        state==NULL
        ||
        state->levels==NULL
        ||
        state->nlevels!=2
        ||
        state->level_capacity<2
        ||
        state->levels[0].A==NULL
        ||
        state->levels[1].A==NULL
        ||
        mode_count<=1
        ||
        terminal_limit<=0
        ||
        fine_coarse_basis==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    current_basis = NULL;
    next_basis = NULL;
    distributed_Ac = NULL;

    cols = NULL;
    vals = NULL;
    next_coarse_basis = NULL;


    err =
        lis_mfamg_create_basis_from_array(
            state->levels[1].A,
            mode_count,
            fine_coarse_basis,
            &current_basis);

    if(err)
        return err;


    level_id = 1;


    while(
        level_id<state->level_capacity-1
        &&
        state->levels[level_id].A->gn
            >
        terminal_limit)
    {
        level =
            &state->levels[level_id];

        next_level =
            &state->levels[level_id+1];


        if(
            level->A==NULL
            ||
            level->A->n!=level->n
            ||
            level->A->gn<=0
            ||
            current_basis==NULL)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        max_row_nnz = 0;

        err =
            lis_mfamg_get_max_row_nnz(
                level->A,
                level->n,
                &max_row_nnz);

        if(err)
            goto fail;


        level->aggregate =
            (LIS_INT *)lis_malloc(
                (size_t)level->n
                    *
                sizeof(LIS_INT),
                "lis_mfamg_extend_hierarchy_multivector::aggregate");

        if(level->aggregate==NULL)
        {
            LIS_SETERR_MEM(
                (size_t)level->n
                    *
                sizeof(LIS_INT));

            err = LIS_OUT_OF_MEMORY;
            goto fail;
        }


        cols =
            (LIS_INT *)lis_malloc(
                (size_t)max_row_nnz
                    *
                sizeof(LIS_INT),
                "lis_mfamg_extend_hierarchy_multivector::cols");

        vals =
            (LIS_SCALAR *)lis_malloc(
                (size_t)max_row_nnz
                    *
                sizeof(LIS_SCALAR),
                "lis_mfamg_extend_hierarchy_multivector::vals");

        if(
            cols==NULL
            ||
            vals==NULL)
        {
            LIS_SETERR_MEM(
                (size_t)max_row_nnz
                    *
                (
                    sizeof(LIS_INT)
                    +
                    sizeof(LIS_SCALAR)
                ));

            err = LIS_OUT_OF_MEMORY;
            goto fail;
        }


        aggregate_target =
            2;

        if(mode_count>1)
        {
            if(
                mode_count
                >
                level->n/2)
            {
                aggregate_target =
                    level->n;
            }
            else
            {
                aggregate_target =
                    2
                    *
                    mode_count;
            }
        }


        err =
            lis_mfamg_build_aggregates_sized(
                level->A,
                level,
                aggregate_target,
                max_row_nnz,
                cols,
                vals);

        if(err)
            goto fail;


        if(
            level->nc<=0
            ||
            level->nc>=level->n)
        {
            lis_free(
                level->aggregate);

            level->aggregate = NULL;
            level->nc = 0;

            lis_free(vals);
            lis_free(cols);

            vals = NULL;
            cols = NULL;

            break;
        }


        coarsened =
            LIS_FALSE;

        next_coarse_basis =
            NULL;


        err =
            lis_mfamg_install_multivector_transfer(
                mode_count,
                current_basis,
                level,
                &next_coarse_basis,
                &coarsened);

        if(err)
            goto fail;


        if(!coarsened)
        {
            lis_free(
                level->aggregate);

            level->aggregate = NULL;
            level->aggregate_count = 0;
            level->nc = 0;

            lis_free(vals);
            lis_free(cols);

            vals = NULL;
            cols = NULL;

            break;
        }


        err =
            lis_mfamg_prepare_level_work(
                level,
                state->smoother_block_size);

        if(err)
            goto fail;


        if(level->A->nprocs>1)
        {
            err =
                lis_mfamg_build_distributed_aggregate_metadata(
                    level->A,
                    level);

            if(err)
                goto fail;


            err =
                lis_mfamg_build_distributed_prolongation(
                    level->A,
                    level);

            if(err)
                goto fail;


            distributed_Ac = NULL;

            err =
                lis_mfamg_build_distributed_coarse(
                    level->A,
                    level,
                    &distributed_Ac);

            if(err)
                goto fail;


            if(
                distributed_Ac==NULL
                ||
                distributed_Ac->matrix_type!=LIS_MATRIX_CSR
                ||
                distributed_Ac->n!=level->nc
                ||
                distributed_Ac->gn!=level->global_nc
                ||
                distributed_Ac->is!=level->coarse_is
                ||
                distributed_Ac->ie!=level->coarse_ie)
            {
                err = LIS_ERR_ILL_ARG;
                goto fail;
            }


            next_level->A =
                distributed_Ac;

            next_level->owns_A =
                LIS_TRUE;

            next_level->n =
                distributed_Ac->n;

            distributed_Ac = NULL;
        }
        else
        {
            /* sparse Galerkin product, no dense nc x nc workspace */
            err = lis_mfamg_create_coarse_csr(
                level->A->comm, level->A, level, LIS_TRUE,
                max_row_nnz, cols, vals, next_level);
            if(err)
                goto fail;
        }


        state->nlevels =
            level_id+2;


        err =
            lis_mfamg_create_basis_from_array(
                next_level->A,
                mode_count,
                next_coarse_basis,
                &next_basis);

        if(err)
            goto fail;


        lis_free(
            next_coarse_basis);

        next_coarse_basis = NULL;


        lis_free(
            vals);

        lis_free(
            cols);

        vals = NULL;
        cols = NULL;


        lis_mfamg_destroy_basis(
            mode_count,
            current_basis);

        current_basis =
            next_basis;

        next_basis =
            NULL;


        level_id++;
    }


    lis_mfamg_destroy_basis(
        mode_count,
        current_basis);


    return LIS_SUCCESS;


fail:

    if(distributed_Ac)
        lis_matrix_destroy(distributed_Ac);

    if(next_coarse_basis)
        lis_free(next_coarse_basis);


    if(vals)
        lis_free(vals);

    if(cols)
        lis_free(cols);

    if(next_basis)
    {
        lis_mfamg_destroy_basis(
            mode_count,
            next_basis);
    }

    if(current_basis)
    {
        lis_mfamg_destroy_basis(
            mode_count,
            current_basis);
    }


    return err;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_extend_hierarchy"
static LIS_INT
lis_mfamg_extend_hierarchy(
    struct LIS_MFAMG_STRUCT *state,
    LIS_VECTOR fine_z)
{
    struct LIS_MFAMG_LEVEL *level;
    struct LIS_MFAMG_LEVEL *next_level;

    LIS_VECTOR current_z;
    LIS_VECTOR next_z;

    LIS_MATRIX distributed_Ac;

    LIS_INT level_id;
    LIS_INT max_row_nnz;
    LIS_INT err;

    LIS_INT *cols;
    LIS_SCALAR *vals;



    if(
        state==NULL
        ||
        fine_z==NULL
        ||
        state->levels==NULL
        ||
        state->nlevels!=2
        ||
        state->level_capacity<2
        ||
        state->levels[0].A==NULL
        ||
        state->levels[1].A==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    current_z = NULL;
    next_z = NULL;
    distributed_Ac = NULL;

    cols = NULL;
    vals = NULL;


    err =
        lis_mfamg_build_coarse_nullspace(
            fine_z,
            &state->levels[0],
            state->levels[1].A,
            &current_z);

    if(err)
        return err;


    level_id = 1;


    /*
     * Use the global matrix size as the hierarchy stop
     * criterion.  In serial A->gn == A->n, so the existing
     * serial hierarchy remains unchanged.
     */
    while(
        level_id<state->level_capacity-1
        &&
        state->levels[level_id].A->gn
            >
        LIS_MFAMG_COARSE_LIMIT)
    {
        level =
            &state->levels[level_id];

        next_level =
            &state->levels[level_id+1];


        if(
            level->A==NULL
            ||
            level->A->n!=level->n
            ||
            level->A->gn<=0)
        {
            err = LIS_ERR_ILL_ARG;
            goto fail;
        }


        max_row_nnz = 0;

        err =
            lis_mfamg_get_max_row_nnz(
                level->A,
                level->n,
                &max_row_nnz);

        if(err)
            goto fail;


        level->aggregate =
            (LIS_INT *)lis_malloc(
                (size_t)level->n
                    *
                sizeof(LIS_INT),
                "lis_mfamg_extend_hierarchy::aggregate");

        if(level->aggregate==NULL)
        {
            LIS_SETERR_MEM(
                (size_t)level->n
                    *
                sizeof(LIS_INT));

            err = LIS_OUT_OF_MEMORY;
            goto fail;
        }


        cols =
            (LIS_INT *)lis_malloc(
                (size_t)max_row_nnz
                    *
                sizeof(LIS_INT),
                "lis_mfamg_extend_hierarchy::cols");

        vals =
            (LIS_SCALAR *)lis_malloc(
                (size_t)max_row_nnz
                    *
                sizeof(LIS_SCALAR),
                "lis_mfamg_extend_hierarchy::vals");

        if(
            cols==NULL
            ||
            vals==NULL)
        {
            LIS_SETERR_MEM(
                (size_t)max_row_nnz
                    *
                (
                    sizeof(LIS_INT)
                    +
                    sizeof(LIS_SCALAR)
                ));

            err = LIS_OUT_OF_MEMORY;
            goto fail;
        }


        err =
            lis_mfamg_build_aggregates(
                level->A,
                level,
                max_row_nnz,
                cols,
                vals);

        if(err)
            goto fail;


        if(
            level->nc<=0
            ||
            level->nc>=level->n)
        {
            lis_free(
                level->aggregate);

            level->aggregate = NULL;
            level->nc = 0;

            lis_free(vals);
            lis_free(cols);

            vals = NULL;
            cols = NULL;

            break;
        }


        level->prolong =
            (LIS_SCALAR *)lis_malloc(
                (size_t)level->n
                    *
                sizeof(LIS_SCALAR),
                "lis_mfamg_extend_hierarchy::prolong");

        if(level->prolong==NULL)
        {
            LIS_SETERR_MEM(
                (size_t)level->n
                    *
                sizeof(LIS_SCALAR));

            err = LIS_OUT_OF_MEMORY;
            goto fail;
        }


        err =
            lis_mfamg_build_prolongation(
                current_z,
                level);

        if(err)
            goto fail;


        err =
            lis_mfamg_prepare_level_work(
                level,
                state->smoother_block_size);

        if(err)
            goto fail;


        if(level->A->nprocs>1)
        {
            /*
             * Distributed path:
             *
             *   local aggregates
             *        |
             *        v
             *   global aggregate ownership
             *        |
             *        v
             *   halo tentative prolongation
             *        |
             *        v
             *   distributed sparse P^H A P
             */
            err =
                lis_mfamg_build_distributed_aggregate_metadata(
                    level->A,
                    level);

            if(err)
                goto fail;


            err =
                lis_mfamg_build_distributed_prolongation(
                    level->A,
                    level);

            if(err)
                goto fail;


            distributed_Ac = NULL;

            err =
                lis_mfamg_build_distributed_coarse(
                    level->A,
                    level,
                    &distributed_Ac);

            if(err)
                goto fail;


            if(
                distributed_Ac==NULL
                ||
                distributed_Ac->n!=level->nc
                ||
                distributed_Ac->gn!=level->global_nc)
            {
                err = LIS_ERR_ILL_ARG;
                goto fail;
            }


            next_level->A =
                distributed_Ac;

            next_level->owns_A =
                LIS_TRUE;

            next_level->n =
                distributed_Ac->n;

            distributed_Ac = NULL;
        }
        else
        {
            /* sparse Galerkin product, no dense nc x nc workspace */
            err = lis_mfamg_create_coarse_csr(
                level->A->comm, level->A, level, LIS_FALSE,
                max_row_nnz, cols, vals, next_level);
            if(err)
                goto fail;
        }


        state->nlevels =
            level_id+2;


        err =
            lis_mfamg_build_coarse_nullspace(
                current_z,
                level,
                next_level->A,
                &next_z);

        if(err)
            goto fail;



        lis_free(vals);
        lis_free(cols);

        vals = NULL;
        cols = NULL;


        lis_vector_destroy(
            current_z);

        current_z =
            next_z;

        next_z = NULL;

        level_id++;
    }


    if(current_z)
        lis_vector_destroy(current_z);


    return LIS_SUCCESS;


fail:

    if(distributed_Ac)
        lis_matrix_destroy(distributed_Ac);


    if(vals)
        lis_free(vals);

    if(cols)
        lis_free(cols);

    if(next_z)
        lis_vector_destroy(next_z);

    if(current_z)
        lis_vector_destroy(current_z);


    return err;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_validate_hierarchy"
static LIS_INT
lis_mfamg_validate_hierarchy(
    struct LIS_MFAMG_STRUCT *state)
{
    LIS_INT level_id;

    struct LIS_MFAMG_LEVEL *level;
    struct LIS_MFAMG_LEVEL *next_level;


    if(
        state==NULL
        ||
        state->levels==NULL
        ||
        state->nlevels<2
        ||
        state->nlevels>state->level_capacity)
    {
        return LIS_ERR_ILL_ARG;
    }


    for(
        level_id=0;
        level_id<state->nlevels;
        level_id++)
    {
        level =
            &state->levels[level_id];

        if(
            level->A==NULL
            ||
            level->n<=0
            ||
            level->A->n!=level->n)
        {
            return LIS_ERR_ILL_ARG;
        }

        if(
            level_id==0
            &&
            level->owns_A)
        {
            return LIS_ERR_ILL_ARG;
        }

        if(
            level_id>0
            &&
            !level->owns_A)
        {
            return LIS_ERR_ILL_ARG;
        }


        if(
            level_id>0
            &&
            (
                level->rhs==NULL
                ||
                level->x==NULL
            ))
        {
            return LIS_ERR_ILL_ARG;
        }


        if(level_id==state->nlevels-1)
            continue;


        next_level =
            &state->levels[level_id+1];

        if(
            level->aggregate==NULL
            ||
            level->prolong==NULL
            ||
            level->diagonal==NULL
            ||
            level->ax==NULL
            ||
            level->nc<=0
            ||
            level->nc!=next_level->n)
        {
            return LIS_ERR_ILL_ARG;
        }
    }


    return LIS_SUCCESS;
}



#undef __FUNC__
#define __FUNC__ "lis_mfamg_factor_csr_terminal"
static LIS_INT
lis_mfamg_factor_csr_terminal(
    LIS_MATRIX A,
    struct LIS_MFAMG_STRUCT *state)
{
    LIS_INT i;
    LIS_INT p;

    LIS_INT local_col;
    LIS_INT global_row;
    LIS_INT global_col;

    LIS_INT local_error;
    LIS_INT global_error;

    size_t pos;
    size_t coarse_entries;


    if(
        A==NULL
        ||
        state==NULL
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->n<=0
        ||
        A->gn!=state->coarse_n
        ||
        A->is<0
        ||
        A->ie<A->is
        ||
        A->ie>A->gn
        ||
        A->ie-A->is!=A->n
        ||
        state->coarse_n<=0
        ||
        state->coarse==NULL
        ||
        state->lu==NULL
        ||
        state->pivots==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    coarse_entries =
        (size_t)state->coarse_n
        *
        (size_t)state->coarse_n;


    memset(
        state->coarse,
        0,
        coarse_entries
        *
        sizeof(LIS_SCALAR));


    local_error = 0;
    global_error = 0;


    for(i=0;i<A->n;i++)
    {
        global_row =
            A->is
            +
            i;


        if(
            global_row<0
            ||
            global_row>=state->coarse_n)
        {
            local_error = 1;
            continue;
        }


        for(
            p=A->ptr[i];
            p<A->ptr[i+1];
            p++)
        {
            local_col =
                A->index[p];


            if(
                local_col<0
                ||
                local_col>=A->np)
            {
                local_error = 1;
                continue;
            }


            if(local_col<A->n)
            {
                global_col =
                    A->is
                    +
                    local_col;
            }
            else
            {
                if(A->l2g_map==NULL)
                {
                    local_error = 1;
                    continue;
                }


                global_col =
                    A->l2g_map[
                        local_col
                        -
                        A->n];
            }


            if(
                global_col<0
                ||
                global_col>=state->coarse_n)
            {
                local_error = 1;
                continue;
            }


            pos =
                (size_t)global_row
                *
                (size_t)state->coarse_n
                +
                (size_t)global_col;


            state->coarse[pos] +=
                A->value[p];
        }
    }


#ifdef USE_MPI

    if(A->nprocs>1)
    {
        MPI_Allreduce(
            &local_error,
            &global_error,
            1,
            LIS_MPI_INT,
            MPI_MAX,
            A->comm);
    }
    else
    {
        global_error =
            local_error;
    }

#else

    global_error =
        local_error;

#endif


    if(global_error)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG terminal CSR contains an invalid distributed column\n");

        return LIS_ERR_ILL_ARG;
    }


#ifdef USE_MPI

    if(A->nprocs>1)
    {
        MPI_Allreduce(
            MPI_IN_PLACE,
            state->coarse,
            (int)coarse_entries,
            LIS_MPI_SCALAR,
            MPI_SUM,
            A->comm);
    }

#endif


    memcpy(
        state->lu,
        state->coarse,
        coarse_entries
        *
        sizeof(LIS_SCALAR));


    return
        lis_mfamg_dense_factor(
            state->coarse_n,
            state->lu,
            state->pivots);
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_prepare_terminal"
static LIS_INT
lis_mfamg_prepare_terminal(
    LIS_MATRIX A,
    struct LIS_MFAMG_STRUCT *state,
    LIS_INT terminal_limit)
{
    size_t coarse_entries;
    size_t allocation_size;

    LIS_INT err;


    if(
        A==NULL
        ||
        state==NULL
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->gn<=0
        ||
        terminal_limit<=0
        ||
        A->gn>terminal_limit
        ||
        state->coarse!=NULL
        ||
        state->lu!=NULL
        ||
        state->pivots!=NULL
        ||
        state->coarse_rhs!=NULL
        ||
        state->coarse_x!=NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    state->coarse_n =
        A->gn;


    if(
        (size_t)state->coarse_n
        >
        ((size_t)-1)
        /
        (size_t)state->coarse_n)
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG terminal matrix dimension overflow\n");

        return LIS_OUT_OF_MEMORY;
    }


    coarse_entries =
        (size_t)state->coarse_n
        *
        (size_t)state->coarse_n;


    if(
        coarse_entries
        >
        ((size_t)-1)
        /
        sizeof(LIS_SCALAR))
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG terminal factor allocation overflow\n");

        return LIS_OUT_OF_MEMORY;
    }


    state->coarse =
        (LIS_SCALAR *)lis_malloc(
            coarse_entries
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_prepare_terminal::coarse");

    state->lu =
        (LIS_SCALAR *)lis_malloc(
            coarse_entries
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_prepare_terminal::lu");

    state->pivots =
        (LIS_INT *)lis_malloc(
            (size_t)state->coarse_n
            *
            sizeof(LIS_INT),
            "lis_mfamg_prepare_terminal::pivots");

    state->coarse_rhs =
        (LIS_SCALAR *)lis_malloc(
            (size_t)state->coarse_n
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_prepare_terminal::rhs");

    state->coarse_x =
        (LIS_SCALAR *)lis_malloc(
            (size_t)state->coarse_n
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_prepare_terminal::x");


    if(
        state->coarse==NULL
        ||
        state->lu==NULL
        ||
        state->pivots==NULL
        ||
        state->coarse_rhs==NULL
        ||
        state->coarse_x==NULL)
    {
        allocation_size =
            2
            *
            coarse_entries
            *
            sizeof(LIS_SCALAR)
            +
            2
            *
            (size_t)state->coarse_n
            *
            sizeof(LIS_SCALAR)
            +
            (size_t)state->coarse_n
            *
            sizeof(LIS_INT);


        LIS_SETERR_MEM(
            allocation_size);

        return LIS_OUT_OF_MEMORY;
    }


    err =
        lis_mfamg_factor_csr_terminal(
            A,
            state);

    if(err)
        return err;


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_dense_solve_distributed"
static LIS_INT
lis_mfamg_dense_solve_distributed(
    LIS_MATRIX A,
    struct LIS_MFAMG_STRUCT *state,
    LIS_VECTOR rhs,
    LIS_VECTOR x)
{
    LIS_INT i;
    LIS_INT err;

    LIS_SCALAR *rhs_value;
    LIS_SCALAR *x_value;


    if(
        A==NULL
        ||
        state==NULL
        ||
        rhs==NULL
        ||
        x==NULL
        ||
        A->nprocs<=1
        ||
        A->gn!=state->coarse_n
        ||
        A->n<=0
        ||
        A->is<0
        ||
        A->ie>A->gn
        ||
        A->ie-A->is!=A->n
        ||
        state->coarse_n<=0
        ||
        state->lu==NULL
        ||
        state->pivots==NULL
        ||
        state->coarse_rhs==NULL
        ||
        state->coarse_x==NULL
        ||
        rhs->n!=A->n
        ||
        rhs->gn!=A->gn
        ||
        x->n!=A->n
        ||
        x->gn!=A->gn)
    {
        return LIS_ERR_ILL_ARG;
    }


    rhs_value =
        rhs->value;

    x_value =
        x->value;


    memset(
        state->coarse_rhs,
        0,
        (size_t)state->coarse_n
        *
        sizeof(LIS_SCALAR));


    for(i=0;i<A->n;i++)
    {
        state->coarse_rhs[
            A->is+i] =
            rhs_value[i];
    }


#ifdef USE_MPI

    MPI_Allreduce(
        MPI_IN_PLACE,
        state->coarse_rhs,
        (int)state->coarse_n,
        LIS_MPI_SCALAR,
        MPI_SUM,
        A->comm);

#else

    return LIS_ERR_NOT_IMPLEMENTED;

#endif


    err =
        lis_mfamg_dense_solve_values(
            state->coarse_n,
            state->coarse,
            state->lu,
            state->pivots,
            state->coarse_rhs,
            state->coarse_x);

    if(err)
        return err;


    for(i=0;i<A->n;i++)
    {
        x_value[i] =
            state->coarse_x[
                A->is+i];
    }


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_v_cycle"
static LIS_INT
lis_mfamg_v_cycle(
    struct LIS_MFAMG_STRUCT *state,
    LIS_INT level_id,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    struct LIS_MFAMG_LEVEL *level;
    struct LIS_MFAMG_LEVEL *next_level;

    LIS_INT i;
    LIS_INT a;
    LIS_INT q;
    LIS_INT p;
    LIS_INT rank;
    LIS_INT row_start;
    LIS_INT err;

    LIS_SCALAR residual;

    LIS_SCALAR *bv;
    LIS_SCALAR *xv;
    LIS_SCALAR *av;
    LIS_SCALAR *coarse_b;
    LIS_SCALAR *coarse_x;


    if(
        state==NULL
        ||
        state->levels==NULL
        ||
        b==NULL
        ||
        x==NULL
        ||
        level_id<0
        ||
        level_id>=state->nlevels)
    {
        return LIS_ERR_ILL_ARG;
    }


    level =
        &state->levels[level_id];


    if(
        level->A==NULL
        ||
        level->n<=0
        ||
        b->n!=level->n
        ||
        x->n!=level->n)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(level_id==state->nlevels-1)
    {
        if(level->A->nprocs>1)
        {
            if(level->A->gn!=state->coarse_n)
                return LIS_ERR_ILL_ARG;

            return
                lis_mfamg_dense_solve_distributed(
                    level->A,
                    state,
                    b,
                    x);
        }


        if(level->n!=state->coarse_n)
            return LIS_ERR_ILL_ARG;


        return
            lis_mfamg_dense_solve(
                state,
                b,
                x);
    }


    next_level =
        &state->levels[level_id+1];


    if(
        level->aggregate==NULL
        ||
        level->prolong==NULL
        ||
        level->diagonal==NULL
        ||
        level->ax==NULL
        ||
        level->nc<=0
        ||
        level->nc!=next_level->n
        ||
        next_level->rhs==NULL
        ||
        next_level->x==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    bv =
        b->value;

    xv =
        x->value;

    av =
        level->ax->value;

    coarse_b =
        next_level->rhs->value;

    coarse_x =
        next_level->x->value;


    /*
     * One pre-smoothing step from zero.
     */
    err =
        lis_mfamg_smoother_apply(
            level,
            b,
            x,
            LIS_TRUE);

    if(err)
        return err;


    err =
        lis_matvec(
            level->A,
            x,
            level->ax);

    if(err)
        return err;


    memset(
        coarse_b,
        0,
        (size_t)next_level->n
            *
        sizeof(LIS_SCALAR));


    /*
     * Restriction:
     *
     *     b_(l+1) = P_l^H (b_l - A_l x_l)
     */
    if(level->prolong_row_ptr)
    {
        if(
            level->aggregate_count<=0
            ||
            level->aggregate_rank==NULL
            ||
            level->aggregate_offset==NULL
            ||
            level->aggregate_offset[
                level->aggregate_count]
                !=
                next_level->n
            ||
            level->prolong_row_ptr[0]!=0)
        {
            return LIS_ERR_ILL_ARG;
        }


        for(i=0;i<level->n;i++)
        {
            residual =
                bv[i]
                -
                av[i];

            a =
                level->aggregate[i];

            if(
                a<0
                ||
                a>=level->aggregate_count)
            {
                return LIS_ERR_ILL_ARG;
            }


            rank =
                level->aggregate_rank[a];

            row_start =
                level->prolong_row_ptr[i];


            if(
                rank<0
                ||
                level->aggregate_offset[a]<0
                ||
                level->aggregate_offset[a]+rank
                    >
                next_level->n
                ||
                row_start<0
                ||
                level->prolong_row_ptr[i+1]
                    -
                row_start
                    !=
                rank)
            {
                return LIS_ERR_ILL_ARG;
            }


            for(q=0;q<rank;q++)
            {
                p =
                    level->aggregate_offset[a]
                    +
                    q;

                coarse_b[p] +=
                    conj(
                        level->prolong[
                            row_start
                            +
                            q])
                    *
                    residual;
            }
        }
    }
    else
    {
        for(i=0;i<level->n;i++)
        {
            residual =
                bv[i]
                -
                av[i];

            a =
                level->aggregate[i];

            if(
                a<0
                ||
                a>=next_level->n)
            {
                return LIS_ERR_ILL_ARG;
            }

            coarse_b[a] +=
                conj(
                    level->prolong[i])
                *
                residual;
        }
    }


    /*
     * Recursive coarse-grid correction.
     */
    err =
        lis_mfamg_v_cycle(
            state,
            level_id+1,
            next_level->rhs,
            next_level->x);

    if(err)
        return err;


    /*
     * Prolongation:
     *
     *     x_l <- x_l + P_l x_(l+1)
     */
    if(level->prolong_row_ptr)
    {
        for(i=0;i<level->n;i++)
        {
            a =
                level->aggregate[i];

            if(
                a<0
                ||
                a>=level->aggregate_count)
            {
                return LIS_ERR_ILL_ARG;
            }


            rank =
                level->aggregate_rank[a];

            row_start =
                level->prolong_row_ptr[i];


            if(
                rank<0
                ||
                level->aggregate_offset[a]<0
                ||
                level->aggregate_offset[a]+rank
                    >
                next_level->n
                ||
                row_start<0
                ||
                level->prolong_row_ptr[i+1]
                    -
                row_start
                    !=
                rank)
            {
                return LIS_ERR_ILL_ARG;
            }


            for(q=0;q<rank;q++)
            {
                p =
                    level->aggregate_offset[a]
                    +
                    q;

                xv[i] +=
                    level->prolong[
                        row_start
                        +
                        q]
                    *
                    coarse_x[p];
            }
        }
    }
    else
    {
        for(i=0;i<level->n;i++)
        {
            a =
                level->aggregate[i];

            if(
                a<0
                ||
                a>=next_level->n)
            {
                return LIS_ERR_ILL_ARG;
            }

            xv[i] +=
                level->prolong[i]
                *
                coarse_x[a];
        }
    }


    /*
     * One post-smoothing step.
     */
    err =
        lis_mfamg_smoother_apply(
            level,
            b,
            x,
            LIS_FALSE);

    if(err)
        return err;


    return LIS_SUCCESS;
}


/*
 * Internal smoother query used by the MFAMG regression test.
 */
LIS_INT
lis_precon_mfamg_get_smoother_omega(
    LIS_MATRIX A,
    LIS_REAL *omega)
{
    LIS_VECTOR diagonal;

    LIS_INT err;


    if(
        A==NULL
        ||
        omega==NULL
        ||
        A->n<=0)
    {
        return LIS_ERR_ILL_ARG;
    }


    diagonal = NULL;


    err =
        lis_vector_duplicate(
            A,
            &diagonal);

    if(err)
        return err;


    err =
        lis_matrix_get_diagonal(
            A,
            diagonal);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_estimate_smoother_omega(
            A,
            diagonal,
            omega);


cleanup:

    if(diagonal)
    {
        lis_vector_destroy(
            diagonal);
    }


    return err;
}


/*
 * Internal distributed aggregation query used by the regression
 * test.  It exercises rank-local aggregation, coarse ownership
 * ranges, global coarse IDs, and integer halo exchange without
 * enabling the distributed MFAMG solve path.
 */
LIS_INT
lis_precon_mfamg_get_distributed_aggregate_info(
    LIS_MATRIX A,
    LIS_INT map_capacity,
    LIS_INT *local_nc,
    LIS_INT *global_nc,
    LIS_INT *coarse_is,
    LIS_INT *coarse_ie,
    LIS_INT *map_n,
    LIS_INT *fine_to_coarse)
{
    struct LIS_MFAMG_LEVEL level;

    LIS_INT max_row_nnz;
    LIS_INT *cols;
    LIS_SCALAR *vals;

    LIS_INT i;
    LIS_INT err;


    if(
        A==NULL
        ||
        local_nc==NULL
        ||
        global_nc==NULL
        ||
        coarse_is==NULL
        ||
        coarse_ie==NULL
        ||
        map_n==NULL
        ||
        fine_to_coarse==NULL
        ||
        map_capacity<0
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs<=1)
    {
        return LIS_ERR_ILL_ARG;
    }


    memset(
        &level,
        0,
        sizeof(level));

    cols = NULL;
    vals = NULL;


    level.A =
        A;

    level.n =
        A->n;


    if(level.n<=0)
    {
        return LIS_ERR_ILL_ARG;
    }


    level.aggregate =
        (LIS_INT *)lis_malloc(
            (size_t)level.n
                *
            sizeof(LIS_INT),
            "lis_precon_mfamg_get_distributed_aggregate_info::aggregate");


    if(level.aggregate==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)level.n
                *
            sizeof(LIS_INT));

        return LIS_OUT_OF_MEMORY;
    }


    err =
        lis_mfamg_get_max_row_nnz(
            A,
            level.n,
            &max_row_nnz);

    if(err)
        goto cleanup;


    cols =
        (LIS_INT *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_INT),
            "lis_precon_mfamg_get_distributed_aggregate_info::cols");

    vals =
        (LIS_SCALAR *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_SCALAR),
            "lis_mfamg_get_distributed_aggregate_info::vals");


    if(
        cols==NULL
        ||
        vals==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)max_row_nnz
                *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        err = LIS_OUT_OF_MEMORY;

        goto cleanup;
    }


    err =
        lis_mfamg_build_aggregates(
            A,
            &level,
            max_row_nnz,
            cols,
            vals);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_build_distributed_aggregate_metadata(
            A,
            &level);

    if(err)
        goto cleanup;


    *map_n =
        A->np;


    if(map_capacity<A->np)
    {
        err = LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    *local_nc =
        level.nc;

    *global_nc =
        level.global_nc;

    *coarse_is =
        level.coarse_is;

    *coarse_ie =
        level.coarse_ie;


    for(i=0;i<A->np;i++)
    {
        fine_to_coarse[i] =
            level.global_aggregate[i];
    }


    err = LIS_SUCCESS;


cleanup:

    if(vals)
        lis_free(vals);

    if(cols)
        lis_free(cols);

    if(level.global_aggregate)
        lis_free(level.global_aggregate);

    if(level.coarse_ranges)
        lis_free(level.coarse_ranges);

    if(level.aggregate)
        lis_free(level.aggregate);


    return err;
}


/*
 * Build owned coarse rows from owned fine rows.
 *
 * Aggregates are rank-owned, so every fine row owned by this
 * rank contributes only to a coarse row owned by this rank.
 * Fine columns may be local or halo entries.  Their exchanged
 * global aggregate IDs and tentative prolongation values
 * preserve off-rank coarse couplings.
 */
static LIS_INT
lis_mfamg_build_distributed_coarse_multivector(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level,
    LIS_MATRIX *coarse_A)
{
    LIS_MATRIX Ac;

    LIS_SCALAR pi;
    LIS_SCALAR pj;
    LIS_SCALAR contribution;

    LIS_REAL mag;

    LIS_INT i;
    LIS_INT p;
    LIS_INT j;

    LIS_INT ai;
    LIS_INT qi;
    LIS_INT qj;

    LIS_INT ri;
    LIS_INT rj;

    LIS_INT row_i;

    LIS_INT coarse_row;
    LIS_INT coarse_col;

    LIS_INT err;


    if(
        A==NULL
        ||
        level==NULL
        ||
        coarse_A==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    *coarse_A = NULL;


    if(
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs<=1
        ||
        A->n<=0
        ||
        A->np<A->n
        ||
        A->ptr==NULL
        ||
        A->index==NULL
        ||
        A->value==NULL
        ||
        level->n!=A->n
        ||
        level->aggregate_count<=0
        ||
        level->nc<=0
        ||
        level->global_nc<=0
        ||
        level->coarse_is<0
        ||
        level->coarse_ie<level->coarse_is
        ||
        level->coarse_ie-level->coarse_is!=level->nc
        ||
        level->aggregate==NULL
        ||
        level->aggregate_rank==NULL
        ||
        level->aggregate_offset==NULL
        ||
        level->prolong_row_ptr==NULL
        ||
        level->prolong==NULL
        ||
        level->global_aggregate==NULL
        ||
        level->global_aggregate_rank==NULL
        ||
        level->global_prolong==NULL
        ||
        level->global_prolong_width<=0
        ||
        level->coarse_ranges==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    Ac = NULL;


    err =
        lis_matrix_create(
            A->comm,
            &Ac);

    if(err)
        return err;


    err =
        lis_matrix_set_size(
            Ac,
            level->nc,
            level->global_nc);

    if(err)
        goto cleanup;


    if(
        Ac->n!=level->nc
        ||
        Ac->gn!=level->global_nc
        ||
        Ac->is!=level->coarse_is
        ||
        Ac->ie!=level->coarse_ie
        ||
        Ac->nprocs!=A->nprocs
        ||
        Ac->my_rank!=A->my_rank
        ||
        Ac->ranges==NULL)
    {
        err = LIS_ERR_ILL_ARG;
        goto cleanup;
    }


    for(i=0;i<=A->nprocs;i++)
    {
        if(
            Ac->ranges[i]
            !=
            level->coarse_ranges[i])
        {
            err = LIS_ERR_ILL_ARG;
            goto cleanup;
        }
    }


    for(i=0;i<A->n;i++)
    {
        ai =
            level->aggregate[i];


        if(
            ai<0
            ||
            ai>=level->aggregate_count)
        {
            err = LIS_ERR_ILL_ARG;
            goto cleanup;
        }


        ri =
            level->aggregate_rank[ai];

        row_i =
            level->prolong_row_ptr[i];


        if(
            ri<0
            ||
            ri>level->global_prolong_width
            ||
            level->aggregate_offset[ai]<0
            ||
            level->aggregate_offset[ai]+ri
                >
            level->nc
            ||
            row_i<0
            ||
            level->prolong_row_ptr[i+1]
                -
            row_i
                !=
            ri
            ||
            level->global_aggregate_rank[i]
                !=
            ri
            ||
            level->global_aggregate[i]
                !=
            level->coarse_is
                +
            level->aggregate_offset[ai])
        {
            err = LIS_ERR_ILL_ARG;
            goto cleanup;
        }


        for(
            p=A->ptr[i];
            p<A->ptr[i+1];
            p++)
        {
            j =
                A->index[p];


            if(
                j<0
                ||
                j>=A->np)
            {
                err = LIS_ERR_ILL_ARG;
                goto cleanup;
            }


            rj =
                level->global_aggregate_rank[j];


            if(
                rj<0
                ||
                rj>level->global_prolong_width
                ||
                level->global_aggregate[j]<0
                ||
                level->global_aggregate[j]>level->global_nc
                ||
                level->global_aggregate[j]+rj
                    >
                level->global_nc)
            {
                err = LIS_ERR_ILL_ARG;
                goto cleanup;
            }


            for(qi=0;qi<ri;qi++)
            {
                coarse_row =
                    level->global_aggregate[i]
                    +
                    qi;

                pi =
                    level->prolong[
                        row_i
                        +
                        qi];


                mag =
                    fabs(
                        pi);


                if(
                    mag!=mag
                    ||
                    mag>LIS_SCALAR_MAX)
                {
                    LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "MFAMG distributed multivector row prolongation is non-finite\n");

                    err = LIS_ERR_ILL_ARG;
                    goto cleanup;
                }


                for(qj=0;qj<rj;qj++)
                {
                    coarse_col =
                        level->global_aggregate[j]
                        +
                        qj;

                    pj =
                        level->global_prolong[
                            (size_t)qj
                            *
                            (size_t)A->np
                            +
                            (size_t)j];


                    contribution =
                        conj(pi)
                        *
                        A->value[p]
                        *
                        pj;


                    mag =
                        fabs(
                            contribution);


                    if(
                        mag!=mag
                        ||
                        mag>LIS_SCALAR_MAX)
                    {
                        LIS_SETERR(
                            LIS_ERR_ILL_ARG,
                            "MFAMG distributed multivector Galerkin contribution is non-finite\n");

                        err = LIS_ERR_ILL_ARG;
                        goto cleanup;
                    }


                    if(
                        contribution
                        ==
                        (LIS_SCALAR)0.0)
                    {
                        continue;
                    }


                    err =
                        lis_matrix_set_value(
                            LIS_ADD_VALUE,
                            coarse_row,
                            coarse_col,
                            contribution,
                            Ac);

                    if(err)
                        goto cleanup;
                }
            }
        }
    }


    err =
        lis_matrix_assemble(
            Ac);

    if(err)
        goto cleanup;


    if(
        Ac->matrix_type!=LIS_MATRIX_CSR
        ||
        Ac->n!=level->nc
        ||
        Ac->gn!=level->global_nc
        ||
        Ac->is!=level->coarse_is
        ||
        Ac->ie!=level->coarse_ie
        ||
        Ac->ptr==NULL
        ||
        Ac->index==NULL
        ||
        Ac->value==NULL
        ||
        Ac->nnz<=0
        ||
        Ac->np<Ac->n)
    {
        err = LIS_ERR_ILL_ARG;
        goto cleanup;
    }


    if(
        Ac->np>Ac->n
        &&
        (
            Ac->l2g_map==NULL
            ||
            Ac->commtable==NULL
        ))
    {
        err = LIS_ERR_ILL_ARG;
        goto cleanup;
    }


    *coarse_A =
        Ac;


    return LIS_SUCCESS;


cleanup:

    if(Ac)
    {
        lis_matrix_destroy(
            Ac);
    }


    return err;
}


static LIS_INT
lis_mfamg_build_distributed_coarse(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level,
    LIS_MATRIX *coarse_A)
{
    LIS_MATRIX Ac;

    LIS_SCALAR contribution;

    LIS_REAL mag;

    LIS_INT i;
    LIS_INT p;
    LIS_INT j;

    LIS_INT coarse_row;
    LIS_INT coarse_col;

    LIS_INT err;


    if(
        level!=NULL
        &&
        level->prolong_row_ptr!=NULL)
    {
        return
            lis_mfamg_build_distributed_coarse_multivector(
                A,
                level,
                coarse_A);
    }


    if(
        A==NULL
        ||
        level==NULL
        ||
        coarse_A==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    *coarse_A = NULL;


    if(
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs<=1
        ||
        A->n<=0
        ||
        A->np<A->n
        ||
        A->ptr==NULL
        ||
        A->index==NULL
        ||
        A->value==NULL
        ||
        level->n!=A->n
        ||
        level->nc<=0
        ||
        level->global_nc<=0
        ||
        level->coarse_is<0
        ||
        level->coarse_ie<level->coarse_is
        ||
        level->coarse_ie-level->coarse_is!=level->nc
        ||
        level->aggregate==NULL
        ||
        level->prolong==NULL
        ||
        level->global_aggregate==NULL
        ||
        level->global_prolong==NULL
        ||
        level->coarse_ranges==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    Ac = NULL;


    err =
        lis_matrix_create(
            A->comm,
            &Ac);

    if(err)
        return err;


    err =
        lis_matrix_set_size(
            Ac,
            level->nc,
            level->global_nc);

    if(err)
        goto cleanup;


    /*
     * Do not call lis_matrix_set_type() here.
     *
     * lis_matrix_set_value() starts from LIS_MATRIX_NULL,
     * creates the temporary RCO assembly storage, and
     * lis_matrix_assemble() converts that storage to the
     * matrix's default CSR type.
     */
    if(
        Ac->n!=level->nc
        ||
        Ac->gn!=level->global_nc
        ||
        Ac->is!=level->coarse_is
        ||
        Ac->ie!=level->coarse_ie
        ||
        Ac->nprocs!=A->nprocs
        ||
        Ac->my_rank!=A->my_rank
        ||
        Ac->ranges==NULL)
    {
        err = LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    for(i=0;i<=A->nprocs;i++)
    {
        if(
            Ac->ranges[i]
            !=
            level->coarse_ranges[i])
        {
            err = LIS_ERR_ILL_ARG;

            goto cleanup;
        }
    }


    for(i=0;i<A->n;i++)
    {
        if(
            level->aggregate[i]<0
            ||
            level->aggregate[i]>=level->nc)
        {
            err = LIS_ERR_ILL_ARG;

            goto cleanup;
        }


        coarse_row =
            level->global_aggregate[i];


        if(
            coarse_row
            !=
            level->coarse_is
            +
            level->aggregate[i])
        {
            err = LIS_ERR_ILL_ARG;

            goto cleanup;
        }


        for(
            p=A->ptr[i];
            p<A->ptr[i+1];
            p++)
        {
            j =
                A->index[p];


            if(
                j<0
                ||
                j>=A->np)
            {
                err = LIS_ERR_ILL_ARG;

                goto cleanup;
            }


            coarse_col =
                level->global_aggregate[j];


            if(
                coarse_col<0
                ||
                coarse_col>=level->global_nc)
            {
                err = LIS_ERR_ILL_ARG;

                goto cleanup;
            }


            contribution =
                conj(
                    level->prolong[i])
                *
                A->value[p]
                *
                level->global_prolong[j];


            mag =
                fabs(
                    contribution);


            if(
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "MFAMG distributed Galerkin contribution is non-finite\n");

                err = LIS_ERR_ILL_ARG;

                goto cleanup;
            }


            if(
                contribution
                ==
                (LIS_SCALAR)0.0)
            {
                continue;
            }


            err =
                lis_matrix_set_value(
                    LIS_ADD_VALUE,
                    coarse_row,
                    coarse_col,
                    contribution,
                    Ac);

            if(err)
                goto cleanup;
        }
    }


    err =
        lis_matrix_assemble(
            Ac);

    if(err)
        goto cleanup;


    if(
        Ac->matrix_type!=LIS_MATRIX_CSR
        ||
        Ac->n!=level->nc
        ||
        Ac->gn!=level->global_nc
        ||
        Ac->is!=level->coarse_is
        ||
        Ac->ie!=level->coarse_ie
        ||
        Ac->ptr==NULL
        ||
        Ac->index==NULL
        ||
        Ac->value==NULL
        ||
        Ac->nnz<=0
        ||
        Ac->np<Ac->n)
    {
        err = LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    if(
        Ac->np>Ac->n
        &&
        (
            Ac->l2g_map==NULL
            ||
            Ac->commtable==NULL
        ))
    {
        err = LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    *coarse_A =
        Ac;


    return LIS_SUCCESS;


cleanup:

    if(Ac)
    {
        lis_matrix_destroy(
            Ac);
    }


    return err;
}


/*
 * Internal distributed tentative-prolongation query used by
 * the regression test.  It does not enable the distributed
 * MFAMG solve path.
 */
LIS_INT
lis_precon_mfamg_get_distributed_prolongation_info(
    LIS_MATRIX A,
    LIS_VECTOR z,
    LIS_INT fine_capacity,
    LIS_INT coarse_capacity,
    LIS_INT *map_n,
    LIS_INT *local_nc,
    LIS_INT *global_nc,
    LIS_SCALAR *fine_prolong,
    LIS_SCALAR *coarse_z)
{
    struct LIS_MFAMG_LEVEL level;

    LIS_INT max_row_nnz;
    LIS_INT *cols;
    LIS_SCALAR *vals;

    LIS_INT i;
    LIS_INT a;
    LIS_INT err;

    LIS_REAL mag;


    if(
        A==NULL
        ||
        z==NULL
        ||
        map_n==NULL
        ||
        local_nc==NULL
        ||
        global_nc==NULL
        ||
        fine_prolong==NULL
        ||
        coarse_z==NULL
        ||
        fine_capacity<0
        ||
        coarse_capacity<0
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs<=1
        ||
        z->n!=A->n)
    {
        return LIS_ERR_ILL_ARG;
    }


    memset(
        &level,
        0,
        sizeof(level));

    cols = NULL;
    vals = NULL;


    level.A =
        A;

    level.n =
        A->n;


    if(level.n<=0)
    {
        return LIS_ERR_ILL_ARG;
    }


    level.aggregate =
        (LIS_INT *)lis_malloc(
            (size_t)level.n
                *
            sizeof(LIS_INT),
            "lis_precon_mfamg_get_distributed_prolongation_info::aggregate");

    level.prolong =
        (LIS_SCALAR *)lis_malloc(
            (size_t)level.n
                *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_get_distributed_prolongation_info::prolong");


    if(
        level.aggregate==NULL
        ||
        level.prolong==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)level.n
                *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        err = LIS_OUT_OF_MEMORY;

        goto cleanup;
    }


    err =
        lis_mfamg_get_max_row_nnz(
            A,
            level.n,
            &max_row_nnz);

    if(err)
        goto cleanup;


    cols =
        (LIS_INT *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_INT),
            "lis_precon_mfamg_get_distributed_prolongation_info::cols");

    vals =
        (LIS_SCALAR *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_get_distributed_prolongation_info::vals");


    if(
        cols==NULL
        ||
        vals==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)max_row_nnz
                *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        err = LIS_OUT_OF_MEMORY;

        goto cleanup;
    }


    err =
        lis_mfamg_build_aggregates(
            A,
            &level,
            max_row_nnz,
            cols,
            vals);

    if(err)
        goto cleanup;


    if(
        level.nc<=0
        ||
        level.nc>=level.n)
    {
        err = LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    err =
        lis_mfamg_build_prolongation(
            z,
            &level);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_build_distributed_aggregate_metadata(
            A,
            &level);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_build_distributed_prolongation(
            A,
            &level);

    if(err)
        goto cleanup;


    *map_n =
        A->np;

    *local_nc =
        level.nc;

    *global_nc =
        level.global_nc;


    if(
        fine_capacity<A->np
        ||
        coarse_capacity<level.nc)
    {
        err = LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    for(i=0;i<A->np;i++)
    {
        fine_prolong[i] =
            level.global_prolong[i];
    }


    for(a=0;a<level.nc;a++)
    {
        coarse_z[a] =
            (LIS_SCALAR)0.0;
    }


    for(i=0;i<level.n;i++)
    {
        a =
            level.aggregate[i];


        if(
            a<0
            ||
            a>=level.nc)
        {
            err = LIS_ERR_ILL_ARG;

            goto cleanup;
        }


        coarse_z[a] +=
            conj(
                level.prolong[i])
            *
            z->value[i];
    }


    for(a=0;a<level.nc;a++)
    {
        mag =
            fabs(
                coarse_z[a]);


        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX
            ||
            mag==0.0)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG distributed coarse near-nullspace value is invalid\n");

            err = LIS_BREAKDOWN;

            goto cleanup;
        }
    }


    err = LIS_SUCCESS;


cleanup:

    if(vals)
        lis_free(vals);

    if(cols)
        lis_free(cols);

    if(level.global_prolong)
        lis_free(level.global_prolong);

    if(level.global_aggregate)
        lis_free(level.global_aggregate);

    if(level.coarse_ranges)
        lis_free(level.coarse_ranges);

    if(level.prolong)
        lis_free(level.prolong);

    if(level.aggregate)
        lis_free(level.aggregate);


    return err;
}


/*
 * Internal distributed hierarchy query used by the regression
 * test.  It exercises the production hierarchy extension path
 * while the public MPI solve path remains disabled.
 */
LIS_INT
lis_precon_mfamg_get_distributed_hierarchy_info(
    LIS_MATRIX A,
    LIS_VECTOR z,
    LIS_INT level_capacity,
    LIS_INT terminal_row_capacity,
    LIS_INT terminal_nnz_capacity,
    LIS_INT *nlevels,
    LIS_INT *local_n,
    LIS_INT *global_n,
    LIS_INT *owns_A,
    LIS_INT *terminal_np,
    LIS_INT *terminal_nnz,
    LIS_INT *terminal_row,
    LIS_INT *terminal_global_cols,
    LIS_SCALAR *terminal_values,
    LIS_INT terminal_dense_capacity,
    LIS_INT terminal_vector_capacity,
    LIS_SCALAR *terminal_dense,
    LIS_SCALAR *terminal_solution)
{
    struct LIS_MFAMG_STRUCT *state;
    struct LIS_MFAMG_LEVEL *level;
    struct LIS_MFAMG_LEVEL *terminal;

    LIS_INT max_row_nnz;

    LIS_INT *work_cols;
    LIS_SCALAR *work_vals;

    LIS_MATRIX first_Ac;

    LIS_INT level_id;
    LIS_INT p;
    LIS_INT local_col;

    LIS_INT err;
    LIS_INT global_row;

    LIS_SCALAR terminal_sum;


    if(
        A==NULL
        ||
        z==NULL
        ||
        nlevels==NULL
        ||
        local_n==NULL
        ||
        global_n==NULL
        ||
        owns_A==NULL
        ||
        terminal_np==NULL
        ||
        terminal_nnz==NULL
        ||
        terminal_row==NULL
        ||
        terminal_global_cols==NULL
        ||
        terminal_values==NULL
        ||
        terminal_dense_capacity<=0
        ||
        terminal_vector_capacity<=0
        ||
        terminal_dense==NULL
        ||
        terminal_solution==NULL
        ||
        level_capacity<=0
        ||
        terminal_row_capacity<=0
        ||
        terminal_nnz_capacity<=0
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs<=1
        ||
        A->n<=0
        ||
        A->gn<=0
        ||
        z->n!=A->n
        ||
        z->gn!=A->gn
        ||
        z->is!=A->is
        ||
        z->ie!=A->ie)
    {
        return LIS_ERR_ILL_ARG;
    }


    state = NULL;
    work_cols = NULL;
    work_vals = NULL;
    first_Ac = NULL;


    state =
        (struct LIS_MFAMG_STRUCT *)lis_malloc(
            sizeof(struct LIS_MFAMG_STRUCT),
            "lis_precon_mfamg_get_distributed_hierarchy_info::state");

    if(state==NULL)
    {
        LIS_SETERR_MEM(
            sizeof(struct LIS_MFAMG_STRUCT));

        return LIS_OUT_OF_MEMORY;
    }


    memset(
        state,
        0,
        sizeof(struct LIS_MFAMG_STRUCT));


    state->level_capacity =
        LIS_MFAMG_MAX_LEVELS;

    state->nlevels =
        2;


    state->levels =
        (struct LIS_MFAMG_LEVEL *)lis_malloc(
            (size_t)state->level_capacity
                *
            sizeof(struct LIS_MFAMG_LEVEL),
            "lis_precon_mfamg_get_distributed_hierarchy_info::levels");

    if(state->levels==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)state->level_capacity
                *
            sizeof(struct LIS_MFAMG_LEVEL));

        err = LIS_OUT_OF_MEMORY;

        goto cleanup;
    }


    memset(
        state->levels,
        0,
        (size_t)state->level_capacity
            *
        sizeof(struct LIS_MFAMG_LEVEL));


    level =
        &state->levels[0];

    level->A =
        A;

    level->owns_A =
        LIS_FALSE;

    level->n =
        A->n;


    level->aggregate =
        (LIS_INT *)lis_malloc(
            (size_t)level->n
                *
            sizeof(LIS_INT),
            "lis_precon_mfamg_get_distributed_hierarchy_info::aggregate");

    level->prolong =
        (LIS_SCALAR *)lis_malloc(
            (size_t)level->n
                *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_get_distributed_hierarchy_info::prolong");

    if(
        level->aggregate==NULL
        ||
        level->prolong==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)level->n
                *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        err = LIS_OUT_OF_MEMORY;

        goto cleanup;
    }


    max_row_nnz = 0;

    err =
        lis_mfamg_get_max_row_nnz(
            A,
            level->n,
            &max_row_nnz);

    if(err)
        goto cleanup;


    work_cols =
        (LIS_INT *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_INT),
            "lis_precon_mfamg_get_distributed_hierarchy_info::work_cols");

    work_vals =
        (LIS_SCALAR *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_get_distributed_hierarchy_info::work_vals");

    if(
        work_cols==NULL
        ||
        work_vals==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)max_row_nnz
                *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        err = LIS_OUT_OF_MEMORY;

        goto cleanup;
    }


    err =
        lis_mfamg_build_aggregates(
            A,
            level,
            max_row_nnz,
            work_cols,
            work_vals);

    if(err)
        goto cleanup;


    if(
        level->nc<=0
        ||
        level->nc>=level->n)
    {
        err = LIS_ERR_ILL_ARG;
        goto cleanup;
    }


    err =
        lis_mfamg_build_prolongation(
            z,
            level);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_build_distributed_aggregate_metadata(
            A,
            level);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_build_distributed_prolongation(
            A,
            level);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_build_distributed_coarse(
            A,
            level,
            &first_Ac);

    if(err)
        goto cleanup;


    state->levels[1].A =
        first_Ac;

    state->levels[1].owns_A =
        LIS_TRUE;

    state->levels[1].n =
        first_Ac->n;

    first_Ac = NULL;


    err =
        lis_mfamg_extend_hierarchy(
            state,
            z);

    if(err)
        goto cleanup;


    /*
     * Mirror the production setup sequence before hierarchy
     * validation.  The extension prepares work vectors for
     * newly created nonterminal coarse levels, but level zero
     * and the per-level cycle rhs/x vectors are prepared after
     * hierarchy construction.
     */
    err =
        lis_mfamg_prepare_level_work(
            &state->levels[0],
            state->smoother_block_size);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_prepare_cycle_vectors(
            state);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_validate_hierarchy(
            state);

    if(err)
        goto cleanup;


    if(
        state->nlevels<=0
        ||
        state->nlevels>level_capacity)
    {
        err = LIS_ERR_ILL_ARG;
        goto cleanup;
    }


    *nlevels =
        state->nlevels;


    for(
        level_id=0;
        level_id<state->nlevels;
        level_id++)
    {
        if(
            state->levels[level_id].A==NULL
            ||
            state->levels[level_id].A->n
                !=
            state->levels[level_id].n)
        {
            err = LIS_ERR_ILL_ARG;
            goto cleanup;
        }


        local_n[level_id] =
            state->levels[level_id].A->n;

        global_n[level_id] =
            state->levels[level_id].A->gn;

        owns_A[level_id] =
            state->levels[level_id].owns_A;
    }


    terminal =
        &state->levels[
            state->nlevels-1];


    if(
        terminal->A==NULL
        ||
        terminal->A->matrix_type!=LIS_MATRIX_CSR
        ||
        terminal_row_capacity
            <
        terminal->A->n+1
        ||
        terminal_nnz_capacity
            <
        terminal->A->nnz)
    {
        err = LIS_ERR_ILL_ARG;
        goto cleanup;
    }


    *terminal_np =
        terminal->A->np;

    *terminal_nnz =
        terminal->A->nnz;


    for(
        level_id=0;
        level_id<=terminal->A->n;
        level_id++)
    {
        terminal_row[level_id] =
            terminal->A->ptr[level_id];
    }


    for(
        p=0;
        p<terminal->A->nnz;
        p++)
    {
        local_col =
            terminal->A->index[p];


        if(
            local_col<0
            ||
            local_col>=terminal->A->np)
        {
            err = LIS_ERR_ILL_ARG;
            goto cleanup;
        }


        if(local_col<terminal->A->n)
        {
            terminal_global_cols[p] =
                terminal->A->is
                +
                local_col;
        }
        else
        {
            if(terminal->A->l2g_map==NULL)
            {
                err = LIS_ERR_ILL_ARG;
                goto cleanup;
            }


            terminal_global_cols[p] =
                terminal->A->l2g_map[
                    local_col
                    -
                    terminal->A->n];
        }


        terminal_values[p] =
            terminal->A->value[p];
    }


    if(
        terminal->A->gn<=0
        ||
        terminal->A->gn>LIS_MFAMG_COARSE_LIMIT
        ||
        (size_t)terminal->A->gn
            >
        ((size_t)-1)
            /
        (size_t)terminal->A->gn
        ||
        (size_t)terminal_dense_capacity
            <
        (size_t)terminal->A->gn
            *
        (size_t)terminal->A->gn
        ||
        terminal_vector_capacity
            <
        terminal->A->gn)
    {
        err = LIS_ERR_ILL_ARG;
        goto cleanup;
    }


    err =
        lis_mfamg_prepare_terminal(
            terminal->A,
            state,
            LIS_MFAMG_COARSE_LIMIT);

    if(err)
        goto cleanup;


    for(
        level_id=0;
        level_id<terminal->A->n;
        level_id++)
    {
        global_row =
            terminal->A->is
            +
            level_id;


        if(
            global_row<0
            ||
            global_row>=state->coarse_n)
        {
            err = LIS_ERR_ILL_ARG;
            goto cleanup;
        }


        terminal_sum =
            (LIS_SCALAR)0.0;


        for(
            local_col=0;
            local_col<state->coarse_n;
            local_col++)
        {
            terminal_sum +=
                state->coarse[
                    global_row
                        *
                    state->coarse_n
                        +
                    local_col]
                *
                (LIS_SCALAR)(
                    local_col+1);
        }


        terminal->rhs->value[
            level_id] =
            terminal_sum;

        terminal->x->value[
            level_id] =
            (LIS_SCALAR)0.0;
    }


    err =
        lis_mfamg_dense_solve_distributed(
            terminal->A,
            state,
            terminal->rhs,
            terminal->x);

    if(err)
        goto cleanup;


    memcpy(
        terminal_dense,
        state->coarse,
        (size_t)state->coarse_n
            *
        (size_t)state->coarse_n
            *
        sizeof(LIS_SCALAR));


    memcpy(
        terminal_solution,
        state->coarse_x,
        (size_t)state->coarse_n
            *
        sizeof(LIS_SCALAR));


    err = LIS_SUCCESS;


cleanup:

    if(first_Ac)
        lis_matrix_destroy(first_Ac);

    if(work_vals)
        lis_free(work_vals);

    if(work_cols)
        lis_free(work_cols);

    if(state)
        lis_mfamg_state_free(state);


    return err;
}


/*
 * Internal distributed Galerkin query used by the regression
 * test.  It builds and assembles the first distributed coarse
 * CSR matrix without enabling the distributed MFAMG solve path.
 */
LIS_INT
lis_precon_mfamg_get_distributed_coarse_info(
    LIS_MATRIX A,
    LIS_VECTOR z,
    LIS_INT row_capacity,
    LIS_INT nnz_capacity,
    LIS_INT *local_nc,
    LIS_INT *global_nc,
    LIS_INT *coarse_is,
    LIS_INT *coarse_ie,
    LIS_INT *coarse_np,
    LIS_INT *local_nnz,
    LIS_INT *row,
    LIS_INT *global_cols,
    LIS_SCALAR *values)
{
    struct LIS_MFAMG_LEVEL level;

    LIS_MATRIX Ac;

    LIS_INT max_row_nnz;

    LIS_INT *work_cols;
    LIS_SCALAR *work_vals;

    LIS_INT i;
    LIS_INT p;
    LIS_INT local_col;

    LIS_INT err;


    if(
        A==NULL
        ||
        z==NULL
        ||
        local_nc==NULL
        ||
        global_nc==NULL
        ||
        coarse_is==NULL
        ||
        coarse_ie==NULL
        ||
        coarse_np==NULL
        ||
        local_nnz==NULL
        ||
        row==NULL
        ||
        global_cols==NULL
        ||
        values==NULL
        ||
        row_capacity<0
        ||
        nnz_capacity<0
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs<=1
        ||
        z->n!=A->n
        ||
        z->gn!=A->gn
        ||
        z->is!=A->is
        ||
        z->ie!=A->ie)
    {
        return LIS_ERR_ILL_ARG;
    }


    memset(
        &level,
        0,
        sizeof(level));


    Ac = NULL;
    work_cols = NULL;
    work_vals = NULL;


    level.A =
        A;

    level.n =
        A->n;


    if(level.n<=0)
    {
        return LIS_ERR_ILL_ARG;
    }


    level.aggregate =
        (LIS_INT *)lis_malloc(
            (size_t)level.n
                *
            sizeof(LIS_INT),
            "lis_precon_mfamg_get_distributed_coarse_info::aggregate");

    level.prolong =
        (LIS_SCALAR *)lis_malloc(
            (size_t)level.n
                *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_get_distributed_coarse_info::prolong");


    if(
        level.aggregate==NULL
        ||
        level.prolong==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)level.n
                *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        err = LIS_OUT_OF_MEMORY;

        goto cleanup;
    }


    err =
        lis_mfamg_get_max_row_nnz(
            A,
            level.n,
            &max_row_nnz);

    if(err)
        goto cleanup;


    work_cols =
        (LIS_INT *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_INT),
            "lis_precon_mfamg_get_distributed_coarse_info::work_cols");

    work_vals =
        (LIS_SCALAR *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_SCALAR),
            "lis_precon_mfamg_get_distributed_coarse_info::work_vals");


    if(
        work_cols==NULL
        ||
        work_vals==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)max_row_nnz
                *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        err = LIS_OUT_OF_MEMORY;

        goto cleanup;
    }


    err =
        lis_mfamg_build_aggregates(
            A,
            &level,
            max_row_nnz,
            work_cols,
            work_vals);

    if(err)
        goto cleanup;


    if(
        level.nc<=0
        ||
        level.nc>=level.n)
    {
        err = LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    err =
        lis_mfamg_build_prolongation(
            z,
            &level);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_build_distributed_aggregate_metadata(
            A,
            &level);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_build_distributed_prolongation(
            A,
            &level);

    if(err)
        goto cleanup;


    err =
        lis_mfamg_build_distributed_coarse(
            A,
            &level,
            &Ac);

    if(err)
        goto cleanup;


    if(
        row_capacity<Ac->n+1
        ||
        nnz_capacity<Ac->nnz)
    {
        err = LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    *local_nc =
        Ac->n;

    *global_nc =
        Ac->gn;

    *coarse_is =
        Ac->is;

    *coarse_ie =
        Ac->ie;

    *coarse_np =
        Ac->np;

    *local_nnz =
        Ac->nnz;


    for(i=0;i<=Ac->n;i++)
    {
        row[i] =
            Ac->ptr[i];
    }


    for(p=0;p<Ac->nnz;p++)
    {
        local_col =
            Ac->index[p];


        if(
            local_col<0
            ||
            local_col>=Ac->np)
        {
            err = LIS_ERR_ILL_ARG;

            goto cleanup;
        }


        if(local_col<Ac->n)
        {
            global_cols[p] =
                Ac->is
                +
                local_col;
        }
        else
        {
            if(Ac->l2g_map==NULL)
            {
                err = LIS_ERR_ILL_ARG;

                goto cleanup;
            }


            global_cols[p] =
                Ac->l2g_map[
                    local_col
                    -
                    Ac->n];
        }


        if(
            global_cols[p]<0
            ||
            global_cols[p]>=Ac->gn)
        {
            err = LIS_ERR_ILL_ARG;

            goto cleanup;
        }


        values[p] =
            Ac->value[p];
    }


    err = LIS_SUCCESS;


cleanup:

    if(Ac)
    {
        lis_matrix_destroy(
            Ac);
    }


    if(work_vals)
        lis_free(work_vals);

    if(work_cols)
        lis_free(work_cols);

    if(level.global_prolong)
        lis_free(level.global_prolong);

    if(level.global_aggregate)
        lis_free(level.global_aggregate);

    if(level.coarse_ranges)
        lis_free(level.coarse_ranges);

    if(level.prolong)
        lis_free(level.prolong);

    if(level.aggregate)
        lis_free(level.aggregate);


    return err;
}


LIS_INT
lis_precon_mfamg_get_level_smoother_info(
    LIS_PRECON precon,
    LIS_INT level_id,
    LIS_INT *smoother_type,
    LIS_INT *block_size,
    LIS_INT *block_count,
    LIS_INT *fallback_blocks)
{
    struct LIS_MFAMG_STRUCT *state;
    struct LIS_MFAMG_LEVEL *level;

    LIS_INT block;
    LIS_INT fallback;


    if(
        precon==NULL
        ||
        precon->mfamg==NULL
        ||
        smoother_type==NULL
        ||
        block_size==NULL
        ||
        block_count==NULL
        ||
        fallback_blocks==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    state =
        precon->mfamg;


    if(
        level_id<0
        ||
        level_id>=state->nlevels)
    {
        return LIS_ERR_ILL_ARG;
    }


    level =
        &state->levels[level_id];


    *smoother_type =
        level->smoother_type;

    *block_size =
        level->block_size;

    *block_count =
        level->block_count;


    fallback =
        0;


    if(
        level->smoother_type
            ==
            LIS_MFAMG_SMOOTHER_BLOCK_JACOBI
        &&
        level->block_count>0
        &&
        level->block_valid!=NULL)
    {
        for(
            block=0;
            block<level->block_count;
            block++)
        {
            if(!level->block_valid[block])
            {
                fallback++;
            }
        }
    }


    *fallback_blocks =
        fallback;


    return LIS_SUCCESS;
}



/*
 * Internal MFAMG hierarchy query used by the regression test.
 *
 * This symbol is intentionally not part of the public LIS API.
 */
#undef __FUNC__
#define __FUNC__ "lis_precon_mfamg_get_hierarchy_info"
LIS_INT
lis_precon_mfamg_get_hierarchy_info(
    LIS_PRECON precon,
    LIS_INT max_levels,
    LIS_INT *nlevels,
    LIS_INT *coarse_n,
    LIS_INT *level_n,
    LIS_INT *owns_A)
{
    struct LIS_MFAMG_STRUCT *state;

    LIS_INT level_id;


    if(
        precon==NULL
        ||
        precon->precon_type
            !=
        LIS_PRECON_TYPE_MFAMG
        ||
        precon->mfamg==NULL
        ||
        max_levels<0
        ||
        nlevels==NULL
        ||
        coarse_n==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    state =
        precon->mfamg;


    if(
        state->levels==NULL
        ||
        state->nlevels<=0
        ||
        state->level_capacity<=0
        ||
        state->nlevels>state->level_capacity
        ||
        state->coarse_n<=0)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        state->nlevels>max_levels
        ||
        level_n==NULL
        ||
        owns_A==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    *nlevels =
        state->nlevels;

    *coarse_n =
        state->coarse_n;


    for(
        level_id=0;
        level_id<state->nlevels;
        level_id++)
    {
        if(
            state->levels[level_id].A==NULL
            ||
            state->levels[level_id].n<=0)
        {
            return LIS_ERR_ILL_ARG;
        }


        level_n[level_id] =
            state->levels[level_id].n;

        owns_A[level_id] =
            state->levels[level_id].owns_A;
    }


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_precon_create_mfamg"
LIS_INT
lis_precon_create_mfamg(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{

    struct LIS_MFAMG_STRUCT *state;
    struct LIS_MFAMG_LEVEL *level;
    struct LIS_MFAMG_LEVEL *coarse_level;
    struct LIS_MFAMG_LEVEL *terminal_level;

    LIS_MATRIX A;

#ifdef USE_MPI
    LIS_MATRIX distributed_Ac;
#endif

    LIS_INT capabilities;
    LIS_INT required;

    LIS_INT max_row_nnz;
    LIS_INT row_nnz;

    LIS_INT *cols;
    LIS_SCALAR *vals;
    LIS_SCALAR *coarse_basis;

    LIS_INT i;
    LIS_INT mode_count;
    LIS_INT terminal_limit;
    LIS_INT coarsened;
    LIS_INT err;



    LIS_DEBUG_FUNC_IN;


    coarse_basis = NULL;
    mode_count = 0;
    terminal_limit =
        LIS_MFAMG_COARSE_LIMIT;
    coarsened = LIS_FALSE;


    if(
        solver==NULL
        ||
        precon==NULL
        ||
        solver->A==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    A =
        solver->A;


    if(A->n<=0)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG requires a nonempty matrix\n");

        return LIS_ERR_ILL_ARG;
    }


    if(
        solver->options[LIS_OPTIONS_PRECISION]
        !=
        LIS_PRECISION_DOUBLE)
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG supports double precision\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }


    if(
        solver->options[LIS_OPTIONS_SCALE]
        !=
        LIS_SCALE_NONE)
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG requires scale=none\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }


    mode_count =
        solver->near_nullspace_dim;


    if(
        mode_count<=0
        ||
        solver->near_nullspace==NULL)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG requires at least one near-nullspace vector\n");

        return LIS_ERR_ILL_ARG;
    }


    if(mode_count>1)
    {
        if(
            mode_count
            >
            A->gn
                /
            LIS_MFAMG_COARSE_LIMIT)
        {
            terminal_limit =
                A->gn;
        }
        else
        {
            terminal_limit =
                LIS_MFAMG_COARSE_LIMIT
                *
                mode_count;
        }
    }


    for(i=0;i<mode_count;i++)
    {
        if(
            solver->near_nullspace[i]==NULL
            ||
            solver->near_nullspace[i]->n!=A->n
            ||
            solver->near_nullspace[i]->gn!=A->gn
            ||
            solver->near_nullspace[i]->is!=A->is
            ||
            solver->near_nullspace[i]->ie!=A->ie)
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "MFAMG near-nullspace vectors are incompatible with the matrix\n");

            return LIS_ERR_ILL_ARG;
        }
    }


    if(
        A->matrix_type!=LIS_MATRIX_CSR
        &&
        A->matrix_type!=LIS_MATRIX_USER)
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG supports CSR and LIS_MATRIX_USER\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }


#ifdef USE_MPI

    if(
        A->nprocs>1
        &&
        A->matrix_type!=LIS_MATRIX_CSR)
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG MPI supports CSR matrices\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }

#endif


    if(A->matrix_type==LIS_MATRIX_USER)
    {
        capabilities = 0;

        err =
            lis_matrix_get_capabilities(
                A,
                &capabilities);

        if(err)
            return err;

        required =
            LIS_MATRIX_CAP_MATVEC
            |
            LIS_MATRIX_CAP_DIAGONAL
            |
            LIS_MATRIX_CAP_ROW_ACCESS;

        if(
            (capabilities & required)
            !=
            required)
        {
            LIS_SETERR(
                LIS_ERR_NOT_IMPLEMENTED,
                "MFAMG USER matrices require matvec, diagonal and row access\n");

            return LIS_ERR_NOT_IMPLEMENTED;
        }
    }


    if(
        solver->mfamg_block_size<1
        ||
        solver->mfamg_block_size>LIS_MFAMG_BLOCK_MAX_SIZE)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG block size must satisfy 1 <= value <= 16\n");

        return LIS_ERR_ILL_ARG;
    }


    state =
        (struct LIS_MFAMG_STRUCT *)
        lis_malloc(
            sizeof(struct LIS_MFAMG_STRUCT),
            "lis_precon_create_mfamg::state");

    if(state==NULL)
    {
        LIS_SETERR_MEM(
            sizeof(struct LIS_MFAMG_STRUCT));

        return LIS_OUT_OF_MEMORY;
    }


    memset(
        state,
        0,
        sizeof(struct LIS_MFAMG_STRUCT));


    state->nlevels =
        2;

    state->level_capacity =
        LIS_MFAMG_MAX_LEVELS;

    state->smoother_block_size =
        solver->mfamg_block_size;

    state->levels =
        (struct LIS_MFAMG_LEVEL *)
        lis_malloc(
            (size_t)state->level_capacity
                *
            sizeof(struct LIS_MFAMG_LEVEL),
            "lis_precon_create_mfamg::levels");

    if(state->levels==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)state->level_capacity
                *
            sizeof(struct LIS_MFAMG_LEVEL));

        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);

        return LIS_OUT_OF_MEMORY;
    }


    memset(
        state->levels,
        0,
        (size_t)state->level_capacity
            *
        sizeof(struct LIS_MFAMG_LEVEL));


    level =
        &state->levels[0];

    coarse_level =
        &state->levels[1];


    level->A =
        A;

    level->owns_A =
        LIS_FALSE;

    level->n =
        A->n;



    level->aggregate =
        (LIS_INT *)lis_malloc(
            (size_t)level->n
                *
            sizeof(LIS_INT),
            "lis_precon_create_mfamg::aggregate");


    if(mode_count==1)
    {
        level->prolong =
            (LIS_SCALAR *)lis_malloc(
                (size_t)level->n
                    *
                sizeof(LIS_SCALAR),
                "lis_precon_create_mfamg::prolong");
    }


    if(
        level->aggregate==NULL
        ||
        (
            mode_count==1
            &&
            level->prolong==NULL
        ))
    {
        LIS_SETERR_MEM(
            (size_t)level->n
                *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);

        return LIS_OUT_OF_MEMORY;
    }


    max_row_nnz = 0;

    for(i=0;i<level->n;i++)
    {
        row_nnz = 0;

        err =
            lis_mfamg_row_query(
                A,
                i,
                &row_nnz);

        if(err)
        {
            if(coarse_basis)
            {
                lis_free(coarse_basis);
                coarse_basis = NULL;
            }

            lis_mfamg_state_free(state);
            return err;
        }

        if(row_nnz<0)
        {
            if(coarse_basis)
            {
                lis_free(coarse_basis);
                coarse_basis = NULL;
            }

            lis_mfamg_state_free(state);
            return LIS_ERR_ILL_ARG;
        }

        if(row_nnz>max_row_nnz)
            max_row_nnz = row_nnz;
    }


    if(max_row_nnz<=0)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG requires matrix rows with entries\n");

        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);

        return LIS_ERR_ILL_ARG;
    }


    cols =
        (LIS_INT *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_INT),
            "lis_precon_create_mfamg::cols");

    vals =
        (LIS_SCALAR *)lis_malloc(
            (size_t)max_row_nnz
                *
            sizeof(LIS_SCALAR),
            "lis_precon_create_mfamg::vals");

    if(
        cols==NULL
        ||
        vals==NULL)
    {
        if(cols)
            lis_free(cols);

        if(vals)
            lis_free(vals);

        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);

        LIS_SETERR_MEM(
            (size_t)max_row_nnz
                *
            (
                sizeof(LIS_INT)
                +
                sizeof(LIS_SCALAR)
            ));

        return LIS_OUT_OF_MEMORY;
    }


    err =
        lis_mfamg_build_aggregates(
            A,
            level,
            max_row_nnz,
            cols,
            vals);

    if(err)
    {
        lis_free(vals);
        lis_free(cols);
        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);

        return err;
    }


    if(mode_count==1)
    {
        err =
            lis_mfamg_build_prolongation(
                solver->near_nullspace[0],
                level);
    }
    else
    {
        coarsened =
            LIS_FALSE;

        err =
            lis_mfamg_install_multivector_transfer(
                mode_count,
                solver->near_nullspace,
                level,
                &coarse_basis,
                &coarsened);
    }


    if(err)
    {
        lis_free(vals);
        lis_free(cols);
        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);

        return err;
    }


    if(
        mode_count>1
        &&
        !coarsened)
    {
        lis_free(vals);
        lis_free(cols);

        LIS_SETERR(
            LIS_BREAKDOWN,
            "MFAMG multivector tentative space does not coarsen the matrix\n");

        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);

        return LIS_BREAKDOWN;
    }


#ifdef USE_MPI

    if(A->nprocs>1)
    {
        distributed_Ac = NULL;


        err =
            lis_mfamg_build_distributed_aggregate_metadata(
                A,
                level);

        if(err)
        {
            lis_free(vals);
            lis_free(cols);
            if(coarse_basis)
            {
                lis_free(coarse_basis);
                coarse_basis = NULL;
            }

            lis_mfamg_state_free(state);

            return err;
        }


        err =
            lis_mfamg_build_distributed_prolongation(
                A,
                level);

        if(err)
        {
            lis_free(vals);
            lis_free(cols);
            if(coarse_basis)
            {
                lis_free(coarse_basis);
                coarse_basis = NULL;
            }

            lis_mfamg_state_free(state);

            return err;
        }


        err =
            lis_mfamg_build_distributed_coarse(
                A,
                level,
                &distributed_Ac);


        lis_free(vals);
        lis_free(cols);

        vals = NULL;
        cols = NULL;


        if(err)
        {
            if(coarse_basis)
            {
                lis_free(coarse_basis);
                coarse_basis = NULL;
            }

            lis_mfamg_state_free(state);
            return err;
        }


        if(
            distributed_Ac==NULL
            ||
            distributed_Ac->matrix_type!=LIS_MATRIX_CSR
            ||
            distributed_Ac->n!=level->nc
            ||
            distributed_Ac->gn!=level->global_nc)
        {
            if(distributed_Ac)
            {
                lis_matrix_destroy(
                    distributed_Ac);
            }


            if(coarse_basis)
            {
                lis_free(coarse_basis);
                coarse_basis = NULL;
            }

            lis_mfamg_state_free(state);

            return LIS_ERR_ILL_ARG;
        }


        coarse_level->A =
            distributed_Ac;

        coarse_level->owns_A =
            LIS_TRUE;

        coarse_level->n =
            distributed_Ac->n;


        distributed_Ac = NULL;
    }
    else

#endif
    {
        /* sparse Galerkin product, no dense nc x nc workspace */
        err =
            lis_mfamg_create_coarse_csr(
                A->comm,
                A,
                level,
                (mode_count==1) ? LIS_FALSE : LIS_TRUE,
                max_row_nnz,
                cols,
                vals,
                coarse_level);
        lis_free(vals);
        lis_free(cols);
        if(err)
        {
            if(coarse_basis)
            {
                lis_free(coarse_basis);
                coarse_basis = NULL;
            }
            lis_mfamg_state_free(state);
            return err;
        }
    }


    if(mode_count==1)
    {
        err =
            lis_mfamg_extend_hierarchy(
                state,
                solver->near_nullspace[0]);
    }
    else
    {
        err =
            lis_mfamg_extend_hierarchy_multivector(
                state,
                mode_count,
                terminal_limit,
                coarse_basis);
    }


    if(coarse_basis)
    {
        lis_free(
            coarse_basis);

        coarse_basis =
            NULL;
    }


    if(err)
    {
        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);
        return err;
    }


    err =
        lis_mfamg_prepare_level_work(
                level,
                state->smoother_block_size);

    if(err)
    {
        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);
        return err;
    }


    err =
        lis_mfamg_prepare_cycle_vectors(
            state);

    if(err)
    {
        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);
        return err;
    }


    err =
        lis_mfamg_validate_hierarchy(
            state);

    if(err)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG hierarchy validation failed\n");

        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);

        return err;
    }


    terminal_level =
        &state->levels[
            state->nlevels-1];


    err =
        lis_mfamg_prepare_terminal(
            terminal_level->A,
            state,
            terminal_limit);

    if(err)
    {
        if(coarse_basis)
        {
            lis_free(coarse_basis);
            coarse_basis = NULL;
        }

        lis_mfamg_state_free(state);
        return err;
    }


    precon->A =
        A;

    precon->is_copy =
        LIS_FALSE;

    precon->mfamg =
        state;


    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}


static LIS_INT
lis_mfamg_same_serial_csr_pattern(
    LIS_MATRIX A,
    LIS_MATRIX B,
    LIS_INT *same)
{
    LIS_INT i;


    if(same==NULL)
        return LIS_ERR_ILL_ARG;


    *same =
        LIS_FALSE;


    if(
        A==NULL
        ||
        B==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        B->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs!=1
        ||
        B->nprocs!=1)
    {
        return LIS_SUCCESS;
    }


    if(
        A->n!=B->n
        ||
        A->gn!=B->gn
        ||
        A->np!=B->np
        ||
        A->is!=B->is
        ||
        A->ie!=B->ie
        ||
        A->nnz!=B->nnz)
    {
        return LIS_SUCCESS;
    }


    if(
        A->n<=0
        ||
        A->nnz<=0
        ||
        A->ptr==NULL
        ||
        B->ptr==NULL
        ||
        A->index==NULL
        ||
        B->index==NULL)
    {
        return LIS_SUCCESS;
    }


    for(i=0;i<=A->n;i++)
    {
        if(A->ptr[i]!=B->ptr[i])
            return LIS_SUCCESS;
    }


    for(i=0;i<A->nnz;i++)
    {
        if(A->index[i]!=B->index[i])
            return LIS_SUCCESS;
    }


    *same =
        LIS_TRUE;


    return LIS_SUCCESS;
}


static LIS_INT
lis_mfamg_ensure_refresh_row_workspace(
    struct LIS_MFAMG_STRUCT *state,
    LIS_INT required_capacity)
{
    LIS_INT *new_cols;

    LIS_SCALAR *new_vals;


    if(
        state==NULL
        ||
        required_capacity<=0
        ||
        state->refresh_row_capacity<0)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        state->refresh_row_capacity
        >=
        required_capacity)
    {
        if(
            state->refresh_cols==NULL
            ||
            state->refresh_vals==NULL)
        {
            return LIS_ERR_ILL_ARG;
        }


        return LIS_SUCCESS;
    }


    if(
        state->refresh_row_capacity>0
        &&
        (
            state->refresh_cols==NULL
            ||
            state->refresh_vals==NULL
        ))
    {
        return LIS_ERR_ILL_ARG;
    }


    new_cols =
        (LIS_INT *)lis_malloc(
            (size_t)required_capacity
            *
            sizeof(LIS_INT),
            "lis_mfamg_ensure_refresh_row_workspace::cols");


    if(new_cols==NULL)
    {
        LIS_SETERR_MEM(
            (size_t)required_capacity
            *
            sizeof(LIS_INT));

        return LIS_OUT_OF_MEMORY;
    }


    new_vals =
        (LIS_SCALAR *)lis_malloc(
            (size_t)required_capacity
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_ensure_refresh_row_workspace::vals");


    if(new_vals==NULL)
    {
        lis_free(
            new_cols);

        LIS_SETERR_MEM(
            (size_t)required_capacity
            *
            sizeof(LIS_SCALAR));

        return LIS_OUT_OF_MEMORY;
    }


    if(state->refresh_cols)
    {
        lis_free(
            state->refresh_cols);
    }


    if(state->refresh_vals)
    {
        lis_free(
            state->refresh_vals);
    }


    state->refresh_cols =
        new_cols;

    state->refresh_vals =
        new_vals;

    state->refresh_row_capacity =
        required_capacity;


    return LIS_SUCCESS;
}


static LIS_INT
lis_mfamg_refresh_level_work(
    struct LIS_MFAMG_LEVEL *level)
{
    return
        lis_mfamg_smoother_refresh(
            level);
}


static LIS_INT
lis_mfamg_refresh_serial_scalar_numeric(
    LIS_SOLVER solver,
    struct LIS_MFAMG_STRUCT *state)
{
    struct LIS_MFAMG_LEVEL *level;
    struct LIS_MFAMG_LEVEL *next_level;
    struct LIS_MFAMG_LEVEL replacement;

    LIS_MATRIX old_coarse_A;

    LIS_VECTOR current_z;
    LIS_VECTOR next_z;

    LIS_INT current_z_owned;

    LIS_INT level_id;
    LIS_INT max_row_nnz;
    LIS_INT err;




    if(
        solver==NULL
        ||
        solver->A==NULL
        ||
        solver->near_nullspace_dim!=1
        ||
        solver->near_nullspace==NULL
        ||
        solver->near_nullspace[0]==NULL
        ||
        state==NULL
        ||
        state->levels==NULL
        ||
        state->nlevels<2
        ||
        state->levels[0].A==NULL
        ||
        solver->A->matrix_type!=LIS_MATRIX_CSR
        ||
        solver->A->nprocs!=1)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        state->alias_rhs==NULL
        ||
        state->coarse_n<=0
        ||
        state->coarse==NULL
        ||
        state->lu==NULL
        ||
        state->pivots==NULL
        ||
        state->coarse_rhs==NULL
        ||
        state->coarse_x==NULL
        ||
        state->refresh_row_capacity<0)
    {
        return LIS_ERR_ILL_ARG;
    }


    current_z =
        solver->near_nullspace[0];

    next_z = NULL;

    current_z_owned =
        LIS_FALSE;



    memset(
        &replacement,
        0,
        sizeof(replacement));


    state->levels[0].A =
        solver->A;

    state->levels[0].owns_A =
        LIS_FALSE;

    state->levels[0].n =
        solver->A->n;


    /* The Galerkin products are built sparsely: no dense workspace. */


    for(
        level_id=0;
        level_id<state->nlevels-1;
        level_id++)
    {
        level =
            &state->levels[level_id];

        next_level =
            &state->levels[level_id+1];


        if(
            level->A==NULL
            ||
            level->A->matrix_type!=LIS_MATRIX_CSR
            ||
            level->A->nprocs!=1
            ||
            level->n!=level->A->n
            ||
            level->n<=0
            ||
            level->nc<=0
            ||
            level->nc>=level->n
            ||
            level->aggregate==NULL
            ||
            level->prolong==NULL
            ||
            next_level->A==NULL
            ||
            !next_level->owns_A
            ||
            next_level->n!=level->nc)
        {
            err =
                LIS_ERR_ILL_ARG;

            goto fail;
        }


        err =
            lis_mfamg_build_prolongation(
                current_z,
                level);

        if(err)
            goto fail;


        err =
            lis_mfamg_refresh_level_work(
                level);

        if(err)
            goto fail;


        max_row_nnz = 0;


        err =
            lis_mfamg_get_max_row_nnz(
                level->A,
                level->n,
                &max_row_nnz);

        if(err)
            goto fail;


        err =
            lis_mfamg_ensure_refresh_row_workspace(
                state,
                max_row_nnz);

        if(err)
            goto fail;


        memset(&replacement,0,sizeof(replacement));
        /* sparse Galerkin product, no dense nc x nc workspace */
        err = lis_mfamg_create_coarse_csr(
            level->A->comm, level->A, level, LIS_FALSE,
            state->refresh_row_capacity,
            state->refresh_cols, state->refresh_vals,
            &replacement);
        if(err)
            goto fail;


        if(
            replacement.A==NULL
            ||
            replacement.n!=next_level->n
            ||
            replacement.A->n!=next_level->n
            ||
            replacement.A->gn!=next_level->A->gn)
        {
            err =
                LIS_ERR_ILL_ARG;

            goto fail;
        }


        err =
            lis_mfamg_build_coarse_nullspace(
                current_z,
                level,
                replacement.A,
                &next_z);

        if(err)
            goto fail;


        old_coarse_A =
            next_level->A;


        next_level->A =
            replacement.A;

        next_level->owns_A =
            LIS_TRUE;

        next_level->n =
            replacement.n;


        replacement.A =
            NULL;

        replacement.owns_A =
            LIS_FALSE;

        replacement.n =
            0;


        lis_matrix_destroy(
            old_coarse_A);


        if(current_z_owned)
        {
            lis_vector_destroy(
                current_z);
        }


        current_z =
            next_z;

        current_z_owned =
            LIS_TRUE;

        next_z = NULL;
    }





    if(current_z_owned)
    {
        lis_vector_destroy(
            current_z);

        current_z = NULL;

        current_z_owned =
            LIS_FALSE;
    }


    err =
        lis_mfamg_validate_hierarchy(
            state);

    if(err)
        return err;


    level =
        &state->levels[
            state->nlevels-1];


    if(
        level->A==NULL
        ||
        level->A->matrix_type!=LIS_MATRIX_CSR
        ||
        level->A->nprocs!=1
        ||
        level->A->gn!=state->coarse_n)
    {
        return LIS_ERR_ILL_ARG;
    }


    err =
        lis_mfamg_factor_csr_terminal(
            level->A,
            state);

    if(err)
        return err;


    return LIS_SUCCESS;


fail:

    if(replacement.A)
    {
        lis_matrix_destroy(
            replacement.A);
    }




    if(next_z)
    {
        lis_vector_destroy(
            next_z);
    }


    if(current_z_owned)
    {
        lis_vector_destroy(
            current_z);
    }


    return err;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_psd_update_mfamg"
LIS_INT
lis_precon_psd_update_mfamg(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    struct LIS_PRECON_STRUCT candidate;
    struct LIS_MFAMG_STRUCT *old_state;

    LIS_INT same_pattern;
    LIS_INT err;


    LIS_DEBUG_FUNC_IN;


    if(
        solver==NULL
        ||
        solver->A==NULL
        ||
        precon==NULL
        ||
        precon->mfamg==NULL
        ||
        precon->precon_type!=LIS_PRECON_TYPE_MFAMG)
    {
        return LIS_ERR_ILL_ARG;
    }


    same_pattern =
        LIS_FALSE;


    if(
        solver->near_nullspace_dim==1
        &&
        solver->near_nullspace!=NULL
        &&
        solver->near_nullspace[0]!=NULL)
    {
        err =
            lis_mfamg_same_serial_csr_pattern(
                solver->A,
                precon->mfamg->levels[0].A,
                &same_pattern);

        if(err)
            return err;
    }


    if(same_pattern)
    {
        err =
            lis_mfamg_refresh_serial_scalar_numeric(
                solver,
                precon->mfamg);

        if(err)
            return err;


        precon->A =
            solver->A;

        precon->is_copy =
            LIS_FALSE;


        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
    }


    /*
     * Unsupported reuse cases keep the transactional full
     * rebuild path.
     */
    memset(
        &candidate,
        0,
        sizeof(candidate));


    err =
        lis_precon_init(
            &candidate);

    if(err)
        return err;


    candidate.precon_type =
        precon->precon_type;


    err =
        lis_precon_create_mfamg(
            solver,
            &candidate);

    if(err)
    {
        if(candidate.mfamg)
        {
            lis_mfamg_state_free(
                candidate.mfamg);

            candidate.mfamg =
                NULL;
        }

        return err;
    }


    if(candidate.mfamg==NULL)
        return LIS_FAILS;


    old_state =
        precon->mfamg;


    precon->A =
        candidate.A;

    precon->is_copy =
        candidate.is_copy;

    precon->mfamg =
        candidate.mfamg;


    candidate.mfamg =
        NULL;


    lis_mfamg_state_free(
        old_state);


    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_destroy_mfamg"
LIS_INT
lis_precon_destroy_mfamg(
    LIS_PRECON precon)
{
    LIS_DEBUG_FUNC_IN;

    if(precon==NULL)
        return LIS_SUCCESS;

    lis_mfamg_state_free(
        precon->mfamg);

    precon->mfamg = NULL;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_psolve_mfamg"
LIS_INT
lis_psolve_mfamg(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    struct LIS_MFAMG_STRUCT *state;

    LIS_VECTOR rhs;

    LIS_INT err;


    LIS_DEBUG_FUNC_IN;


    if(
        solver==NULL
        ||
        solver->precon==NULL
        ||
        solver->precon->mfamg==NULL
        ||
        b==NULL
        ||
        x==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    state =
        solver->precon->mfamg;


    if(
        state->nlevels<2
        ||
        state->levels==NULL
        ||
        state->alias_rhs==NULL
        ||
        b->n!=state->levels[0].n
        ||
        x->n!=state->levels[0].n
        ||
        state->alias_rhs->n!=state->levels[0].n)
    {
        return LIS_ERR_ILL_ARG;
    }


    rhs = b;


    if(b==x)
    {
        err =
            lis_vector_copy(
                b,
                state->alias_rhs);

        if(err)
            return err;


        rhs =
            state->alias_rhs;
    }


    err =
        lis_mfamg_v_cycle(
            state,
            0,
            rhs,
            x);

    if(err)
        return err;


    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}
