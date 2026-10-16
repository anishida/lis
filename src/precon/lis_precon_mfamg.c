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

    LIS_VECTOR diagonal;
    LIS_VECTOR ax;

    LIS_REAL omega;

    LIS_VECTOR rhs;
    LIS_VECTOR x;
};


struct LIS_MFAMG_STRUCT
{
    LIS_INT nlevels;
    LIS_INT level_capacity;
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

            if(level->diagonal)
                lis_vector_destroy(level->diagonal);

            if(level->ax)
                lis_vector_destroy(level->ax);

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

        if(scale==0.0)
        {
            lis_free(
                row_scale);

            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse matrix has a zero row\n");

            return LIS_BREAKDOWN;
        }

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

            ratio =
                mag
                /
                row_scale[i];

            if(ratio>best_ratio)
            {
                best_ratio = ratio;
                p = i;
            }
        }

        if(
            p<0
            ||
            fabs(
                lu[p*n+k])
            ==
            0.0)
        {
            lis_free(
                row_scale);

            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse matrix is singular\n");

            return LIS_BREAKDOWN;
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
lis_mfamg_dense_solve(
    struct LIS_MFAMG_STRUCT *state,
    LIS_VECTOR rhs,
    LIS_VECTOR x)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT p;

    LIS_REAL mag;

    LIS_SCALAR tmp;
    LIS_SCALAR sum;

    LIS_SCALAR *rhs_value;
    LIS_SCALAR *x_value;


    if(
        state==NULL
        ||
        rhs==NULL
        ||
        x==NULL
        ||
        state->coarse_n<=0
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


    rhs_value =
        rhs->value;

    x_value =
        x->value;


    memcpy(
        x_value,
        rhs_value,
        (size_t)state->coarse_n
            *
        sizeof(LIS_SCALAR));


    for(k=0;k<state->coarse_n;k++)
    {
        p =
            state->pivots[k];

        if(
            p<0
            ||
            p>=state->coarse_n)
        {
            return LIS_ERR_ILL_ARG;
        }

        if(p!=k)
        {
            tmp =
                x_value[k];

            x_value[k] =
                x_value[p];

            x_value[p] =
                tmp;
        }
    }


    for(i=1;i<state->coarse_n;i++)
    {
        sum =
            x_value[i];

        for(j=0;j<i;j++)
        {
            sum -=
                state->lu[
                    i*state->coarse_n+j]
                *
                x_value[j];
        }

        x_value[i] =
            sum;
    }


    for(i=state->coarse_n-1;i>=0;i--)
    {
        sum =
            x_value[i];

        for(j=i+1;j<state->coarse_n;j++)
        {
            sum -=
                state->lu[
                    i*state->coarse_n+j]
                *
                x_value[j];
        }

        if(
            fabs(
                state->lu[
                    i*state->coarse_n+i])
            ==
            0.0)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse solve encountered a zero pivot\n");

            return LIS_BREAKDOWN;
        }

        x_value[i] =
            sum
            /
            state->lu[
                i*state->coarse_n+i];

        mag =
            fabs(
                x_value[i]);

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
            rank<=0
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
            rank<=0
            ||
            global_aggregate[i]<0
            ||
            global_aggregate[i]>=global_nc
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
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG near-nullspace vanishes on an aggregate\n");

            err = LIS_BREAKDOWN;
            goto fail;
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


        if(rank<=0)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG multivector aggregate has zero numerical rank\n");

            err = LIS_BREAKDOWN;
            goto fail;
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

        if(rank<=0)
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
            rank<=0
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
            rank<=0
            ||
            rank>width
            ||
            level->global_aggregate[i]<0
            ||
            level->global_aggregate[i]>=level->global_nc
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


#undef __FUNC__
#define __FUNC__ "lis_mfamg_build_coarse"
static LIS_INT
lis_mfamg_build_coarse(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level,
    LIS_SCALAR *coarse,
    LIS_INT max_row_nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT nnz;

    LIS_INT ai;
    LIS_INT aj;

    LIS_SCALAR pi;


    memset(
        coarse,
        0,
        (size_t)level->nc
            *
        (size_t)level->nc
            *
        sizeof(LIS_SCALAR));


    for(i=0;i<level->n;i++)
    {
        ai =
            level->aggregate[i];

        pi =
            level->prolong[i];

        nnz =
            max_row_nnz;

        if(
            lis_mfamg_row_fill(
                A,
                i,
                max_row_nnz,
                &nnz,
                cols,
                vals)
            !=
            LIS_SUCCESS)
        {
            return LIS_ERR_ILL_ARG;
        }

        for(k=0;k<nnz;k++)
        {
            j =
                cols[k]
                -
                A->is;

            if(
                j<0
                ||
                j>=level->n)
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "MFAMG serial row access returned a nonlocal column\n");

                return LIS_ERR_ILL_ARG;
            }

            aj =
                level->aggregate[j];

            coarse[
                ai*level->nc
                +
                aj]
                +=
                conj(pi)
                *
                vals[k]
                *
                level->prolong[j];
        }
    }


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_create_csr_from_dense"
static LIS_INT
lis_mfamg_create_csr_from_dense(
    LIS_Comm comm,
    LIS_INT n,
    LIS_SCALAR *dense,
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_MATRIX A;

    LIS_INT *ptr;
    LIS_INT *index;
    LIS_SCALAR *value;

    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT nnz;
    LIS_INT err;

    size_t pos;
    size_t nnz_size;

    LIS_SCALAR aij;


    if(
        n<=0
        ||
        dense==NULL
        ||
        level==NULL
        ||
        level->A!=NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    A = NULL;
    ptr = NULL;
    index = NULL;
    value = NULL;

    nnz_size = 0;


    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            pos =
                (size_t)i
                    *
                (size_t)n
                    +
                (size_t)j;

            aij =
                dense[pos];

            if(aij!=(LIS_SCALAR)0.0)
                nnz_size++;
        }
    }


    if(nnz_size==0)
    {
        LIS_SETERR(
            LIS_BREAKDOWN,
            "MFAMG coarse CSR matrix has no entries\n");

        return LIS_BREAKDOWN;
    }


    nnz =
        (LIS_INT)nnz_size;

    if(
        nnz<0
        ||
        (size_t)nnz!=nnz_size)
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG coarse CSR nonzero count overflow\n");

        return LIS_OUT_OF_MEMORY;
    }


    err =
        lis_matrix_create(
            comm,
            &A);

    if(err)
        return err;


    err =
        lis_matrix_set_size(
            A,
            n,
            n);

    if(err)
    {
        lis_matrix_destroy(A);
        return err;
    }


    err =
        lis_matrix_malloc_csr(
            n,
            nnz,
            &ptr,
            &index,
            &value);

    if(err)
    {
        lis_matrix_destroy(A);
        return err;
    }


    k = 0;
    ptr[0] = 0;

    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            pos =
                (size_t)i
                    *
                (size_t)n
                    +
                (size_t)j;

            aij =
                dense[pos];

            if(aij==(LIS_SCALAR)0.0)
                continue;

            index[k] =
                j;

            value[k] =
                aij;

            k++;
        }

        ptr[i+1] =
            k;
    }


    if(k!=nnz)
    {
        lis_free(ptr);
        lis_free(index);
        lis_free(value);
        lis_matrix_destroy(A);

        return LIS_ERR_ILL_ARG;
    }


    err =
        lis_matrix_set_csr(
            nnz,
            ptr,
            index,
            value,
            A);

    if(err)
    {
        lis_free(ptr);
        lis_free(index);
        lis_free(value);
        lis_matrix_destroy(A);

        return err;
    }


    /*
     * Ownership of ptr/index/value now belongs to A.
     */
    ptr = NULL;
    index = NULL;
    value = NULL;


    err =
        lis_matrix_assemble(A);

    if(err)
    {
        lis_matrix_destroy(A);
        return err;
    }


    if(
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->n!=n
        ||
        A->gn!=n)
    {
        lis_matrix_destroy(A);

        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG internal coarse CSR matrix is inconsistent\n");

        return LIS_ERR_ILL_ARG;
    }


    level->A =
        A;

    level->owns_A =
        LIS_TRUE;

    level->n =
        n;

    level->nc =
        0;


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_mfamg_prepare_level_work"
static LIS_INT
lis_mfamg_prepare_level_work(
    struct LIS_MFAMG_LEVEL *level)
{
    LIS_INT i;
    LIS_INT err;

    LIS_REAL diagonal_mag;


    if(
        level==NULL
        ||
        level->A==NULL
        ||
        level->n<=0)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(
        level->diagonal!=NULL
        ||
        level->ax!=NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    err =
        lis_vector_duplicate(
            level->A,
            &level->diagonal);

    if(err)
        return err;


    err =
        lis_vector_duplicate(
            level->A,
            &level->ax);

    if(err)
        return err;


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


#undef __FUNC__
#define __FUNC__ "lis_mfamg_build_coarse_multivector"
static LIS_INT
lis_mfamg_build_coarse_multivector(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level,
    LIS_SCALAR *coarse,
    LIS_INT max_row_nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT nnz;

    LIS_INT ai;
    LIS_INT aj;

    LIS_INT qi;
    LIS_INT qj;

    LIS_INT ri;
    LIS_INT rj;

    LIS_INT row_i;
    LIS_INT row_j;

    LIS_INT ci;
    LIS_INT cj;

    LIS_SCALAR pi;
    LIS_SCALAR pj;


    if(
        A==NULL
        ||
        level==NULL
        ||
        coarse==NULL
        ||
        cols==NULL
        ||
        vals==NULL
        ||
        A->nprocs!=1
        ||
        level->n<=0
        ||
        level->nc<=0
        ||
        level->aggregate_count<=0
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
        level->aggregate_offset[0]!=0
        ||
        level->aggregate_offset[
            level->aggregate_count]
            !=
            level->nc
        ||
        level->prolong_row_ptr[0]!=0)
    {
        return LIS_ERR_ILL_ARG;
    }


    memset(
        coarse,
        0,
        (size_t)level->nc
            *
        (size_t)level->nc
            *
        sizeof(LIS_SCALAR));


    for(i=0;i<level->n;i++)
    {
        ai =
            level->aggregate[i];

        if(
            ai<0
            ||
            ai>=level->aggregate_count)
        {
            return LIS_ERR_ILL_ARG;
        }


        ri =
            level->aggregate_rank[ai];

        row_i =
            level->prolong_row_ptr[i];


        if(
            ri<=0
            ||
            row_i<0
            ||
            level->prolong_row_ptr[i+1]
                -
            row_i
                !=
            ri
            ||
            level->aggregate_offset[ai]<0
            ||
            level->aggregate_offset[ai]+ri
                >
            level->nc)
        {
            return LIS_ERR_ILL_ARG;
        }


        nnz =
            max_row_nnz;

        if(
            lis_mfamg_row_fill(
                A,
                i,
                max_row_nnz,
                &nnz,
                cols,
                vals)
            !=
            LIS_SUCCESS)
        {
            return LIS_ERR_ILL_ARG;
        }


        for(k=0;k<nnz;k++)
        {
            j =
                cols[k]
                -
                A->is;

            if(
                j<0
                ||
                j>=level->n)
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "MFAMG serial multivector row contains a nonlocal column\n");

                return LIS_ERR_ILL_ARG;
            }


            aj =
                level->aggregate[j];

            if(
                aj<0
                ||
                aj>=level->aggregate_count)
            {
                return LIS_ERR_ILL_ARG;
            }


            rj =
                level->aggregate_rank[aj];

            row_j =
                level->prolong_row_ptr[j];


            if(
                rj<=0
                ||
                row_j<0
                ||
                level->prolong_row_ptr[j+1]
                    -
                row_j
                    !=
                rj
                ||
                level->aggregate_offset[aj]<0
                ||
                level->aggregate_offset[aj]+rj
                    >
                level->nc)
            {
                return LIS_ERR_ILL_ARG;
            }


            for(qi=0;qi<ri;qi++)
            {
                ci =
                    level->aggregate_offset[ai]
                    +
                    qi;

                pi =
                    level->prolong[
                        row_i
                        +
                        qi];


                for(qj=0;qj<rj;qj++)
                {
                    cj =
                        level->aggregate_offset[aj]
                        +
                        qj;

                    pj =
                        level->prolong[
                            row_j
                            +
                            qj];


                    coarse[
                        (size_t)ci
                            *
                        (size_t)level->nc
                            +
                        (size_t)cj]
                        +=
                        conj(pi)
                        *
                        vals[k]
                        *
                        pj;
                }
            }
        }
    }


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
    LIS_SCALAR *dense;
    LIS_SCALAR *next_coarse_basis;

    size_t dense_entries;


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
    dense = NULL;
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
                level);

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
            if(
                (size_t)level->nc
                >
                ((size_t)-1)
                /
                (size_t)level->nc)
            {
                LIS_SETERR(
                    LIS_OUT_OF_MEMORY,
                    "MFAMG multivector hierarchy coarse dimension overflow\n");

                err = LIS_OUT_OF_MEMORY;
                goto fail;
            }


            dense_entries =
                (size_t)level->nc
                    *
                (size_t)level->nc;


            if(
                dense_entries
                >
                ((size_t)-1)
                /
                sizeof(LIS_SCALAR))
            {
                LIS_SETERR(
                    LIS_OUT_OF_MEMORY,
                    "MFAMG multivector hierarchy allocation overflow\n");

                err = LIS_OUT_OF_MEMORY;
                goto fail;
            }


            dense =
                (LIS_SCALAR *)lis_malloc(
                    dense_entries
                        *
                    sizeof(LIS_SCALAR),
                    "lis_mfamg_extend_hierarchy_multivector::dense");

            if(dense==NULL)
            {
                LIS_SETERR_MEM(
                    dense_entries
                        *
                    sizeof(LIS_SCALAR));

                err = LIS_OUT_OF_MEMORY;
                goto fail;
            }


            err =
                lis_mfamg_build_coarse_multivector(
                    level->A,
                    level,
                    dense,
                    max_row_nnz,
                    cols,
                    vals);

            if(err)
                goto fail;


            err =
                lis_mfamg_create_csr_from_dense(
                    level->A->comm,
                    level->nc,
                    dense,
                    next_level);

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


        if(dense)
        {
            lis_free(
                dense);

            dense = NULL;
        }


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

    if(dense)
        lis_free(dense);

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
    LIS_SCALAR *dense;

    size_t dense_entries;


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
    dense = NULL;


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
                level);

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
            /*
             * Preserve the existing serial dense Galerkin
             * construction exactly.
             */
            if(
                (size_t)level->nc
                >
                ((size_t)-1)
                /
                (size_t)level->nc)
            {
                LIS_SETERR(
                    LIS_OUT_OF_MEMORY,
                    "MFAMG hierarchy coarse dimension overflow\n");

                err = LIS_OUT_OF_MEMORY;
                goto fail;
            }


            dense_entries =
                (size_t)level->nc
                    *
                (size_t)level->nc;


            if(
                dense_entries
                >
                ((size_t)-1)
                /
                sizeof(LIS_SCALAR))
            {
                LIS_SETERR(
                    LIS_OUT_OF_MEMORY,
                    "MFAMG hierarchy coarse allocation overflow\n");

                err = LIS_OUT_OF_MEMORY;
                goto fail;
            }


            dense =
                (LIS_SCALAR *)lis_malloc(
                    dense_entries
                        *
                    sizeof(LIS_SCALAR),
                    "lis_mfamg_extend_hierarchy::dense");

            if(dense==NULL)
            {
                LIS_SETERR_MEM(
                    dense_entries
                        *
                    sizeof(LIS_SCALAR));

                err = LIS_OUT_OF_MEMORY;
                goto fail;
            }


            err =
                lis_mfamg_build_coarse(
                    level->A,
                    level,
                    dense,
                    max_row_nnz,
                    cols,
                    vals);

            if(err)
                goto fail;


            err =
                lis_mfamg_create_csr_from_dense(
                    level->A->comm,
                    level->nc,
                    dense,
                    next_level);

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


        if(dense)
        {
            lis_free(dense);
            dense = NULL;
        }

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

    if(dense)
        lis_free(dense);

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
    LIS_INT j;
    LIS_INT k;
    LIS_INT p;

    LIS_SCALAR tmp;
    LIS_SCALAR sum;

    LIS_REAL mag;

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


    memcpy(
        state->coarse_x,
        state->coarse_rhs,
        (size_t)state->coarse_n
        *
        sizeof(LIS_SCALAR));


    for(k=0;k<state->coarse_n;k++)
    {
        p =
            state->pivots[k];


        if(
            p<0
            ||
            p>=state->coarse_n)
        {
            return LIS_ERR_ILL_ARG;
        }


        if(p!=k)
        {
            tmp =
                state->coarse_x[k];

            state->coarse_x[k] =
                state->coarse_x[p];

            state->coarse_x[p] =
                tmp;
        }
    }


    for(i=1;i<state->coarse_n;i++)
    {
        sum =
            state->coarse_x[i];


        for(j=0;j<i;j++)
        {
            sum -=
                state->lu[
                    i*state->coarse_n+j]
                *
                state->coarse_x[j];
        }


        state->coarse_x[i] =
            sum;
    }


    for(
        i=state->coarse_n-1;
        i>=0;
        i--)
    {
        sum =
            state->coarse_x[i];


        for(
            j=i+1;
            j<state->coarse_n;
            j++)
        {
            sum -=
                state->lu[
                    i*state->coarse_n+j]
                *
                state->coarse_x[j];
        }


        if(
            fabs(
                state->lu[
                    i*state->coarse_n+i])
            ==
            0.0)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "MFAMG coarse solve encountered a zero pivot\n");

            return LIS_BREAKDOWN;
        }


        state->coarse_x[i] =
            sum
            /
            state->lu[
                i*state->coarse_n+i];


        mag =
            fabs(
                state->coarse_x[i]);


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
    LIS_SCALAR *diag;
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

    diag =
        level->diagonal->value;

    coarse_b =
        next_level->rhs->value;

    coarse_x =
        next_level->x->value;


    /*
     * One weighted-Jacobi pre-smoothing step from zero.
     */
    for(i=0;i<level->n;i++)
    {
        xv[i] =
            level->omega
            *
            bv[i]
            /
            diag[i];
    }


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
                rank<=0
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
                rank<=0
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


    err =
        lis_matvec(
            level->A,
            x,
            level->ax);

    if(err)
        return err;


    /*
     * One weighted-Jacobi post-smoothing step.
     */
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
            ri<=0
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
                rj<=0
                ||
                rj>level->global_prolong_width
                ||
                level->global_aggregate[j]<0
                ||
                level->global_aggregate[j]>=level->global_nc
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
            &state->levels[0]);

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
    size_t coarse_entries;

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
        if(
            (size_t)level->nc
            >
            ((size_t)-1)
            /
            (size_t)level->nc)
        {
            LIS_SETERR(
                LIS_OUT_OF_MEMORY,
                "MFAMG coarse matrix dimension overflow\n");

            lis_free(vals);
            lis_free(cols);
            if(coarse_basis)
            {
                lis_free(coarse_basis);
                coarse_basis = NULL;
            }

            lis_mfamg_state_free(state);

            return LIS_OUT_OF_MEMORY;
        }


        coarse_entries =
            (size_t)level->nc
                *
            (size_t)level->nc;


        if(
            coarse_entries
            >
            ((size_t)-1)
            /
            sizeof(LIS_SCALAR))
        {
            LIS_SETERR(
                LIS_OUT_OF_MEMORY,
                "MFAMG coarse matrix allocation size overflow\n");

            lis_free(vals);
            lis_free(cols);
            if(coarse_basis)
            {
                lis_free(coarse_basis);
                coarse_basis = NULL;
            }

            lis_mfamg_state_free(state);

            return LIS_OUT_OF_MEMORY;
        }


        state->coarse =
            (LIS_SCALAR *)lis_malloc(
                coarse_entries
                    *
                sizeof(LIS_SCALAR),
                "lis_precon_create_mfamg::coarse");

        if(state->coarse==NULL)
        {
            lis_free(vals);
            lis_free(cols);

            LIS_SETERR_MEM(
                coarse_entries
                    *
                sizeof(LIS_SCALAR));

            if(coarse_basis)
            {
                lis_free(coarse_basis);
                coarse_basis = NULL;
            }

            lis_mfamg_state_free(state);

            return LIS_OUT_OF_MEMORY;
        }


        if(mode_count==1)
        {
            err =
                lis_mfamg_build_coarse(
                    A,
                    level,
                    state->coarse,
                    max_row_nnz,
                    cols,
                    vals);
        }
        else
        {
            err =
                lis_mfamg_build_coarse_multivector(
                    A,
                    level,
                    state->coarse,
                    max_row_nnz,
                    cols,
                    vals);
        }

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


        err =
            lis_mfamg_create_csr_from_dense(
                A->comm,
                level->nc,
                state->coarse,
                coarse_level);

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


        lis_free(
            state->coarse);

        state->coarse = NULL;
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
            level);

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
    LIS_INT i;
    LIS_INT err;

    LIS_REAL diagonal_mag;


    if(
        level==NULL
        ||
        level->A==NULL
        ||
        level->n<=0
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

    LIS_SCALAR *dense;

    size_t dense_entries;
    size_t level_dense_entries;


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

    dense = NULL;

    dense_entries = 0;

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


    /*
     * The hierarchy topology is fixed on this reuse path, so
     * the maximum dense Galerkin workspace size is known before
     * numeric refresh begins.  Keep this buffer transient to
     * avoid retaining O(nc^2) storage between updates.
     */
    for(
        level_id=0;
        level_id<state->nlevels-1;
        level_id++)
    {
        level =
            &state->levels[level_id];


        if(
            level->A==NULL
            ||
            level->n<=0
            ||
            level->nc<=0
            ||
            level->nc>=level->n)
        {
            err =
                LIS_ERR_ILL_ARG;

            goto fail;
        }


        if(
            (size_t)level->nc
            >
            ((size_t)-1)
            /
            (size_t)level->nc)
        {
            LIS_SETERR(
                LIS_OUT_OF_MEMORY,
                "MFAMG numeric refresh coarse dimension overflow\n");

            err =
                LIS_OUT_OF_MEMORY;

            goto fail;
        }


        level_dense_entries =
            (size_t)level->nc
            *
            (size_t)level->nc;


        if(
            level_dense_entries
            >
            ((size_t)-1)
            /
            sizeof(LIS_SCALAR))
        {
            LIS_SETERR(
                LIS_OUT_OF_MEMORY,
                "MFAMG numeric refresh coarse allocation overflow\n");

            err =
                LIS_OUT_OF_MEMORY;

            goto fail;
        }


        if(
            level_dense_entries
            >
            dense_entries)
        {
            dense_entries =
                level_dense_entries;
        }
    }


    if(dense_entries==0)
    {
        err =
            LIS_ERR_ILL_ARG;

        goto fail;
    }


    dense =
        (LIS_SCALAR *)lis_malloc(
            dense_entries
            *
            sizeof(LIS_SCALAR),
            "lis_mfamg_refresh_serial_scalar_numeric::dense");


    if(dense==NULL)
    {
        LIS_SETERR_MEM(
            dense_entries
            *
            sizeof(LIS_SCALAR));

        err =
            LIS_OUT_OF_MEMORY;

        goto fail;
    }


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


        err =
            lis_mfamg_build_coarse(
                level->A,
                level,
                dense,
                state->refresh_row_capacity,
                state->refresh_cols,
                state->refresh_vals);

        if(err)
            goto fail;


        memset(
            &replacement,
            0,
            sizeof(replacement));


        err =
            lis_mfamg_create_csr_from_dense(
                level->A->comm,
                level->nc,
                dense,
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


    lis_free(
        dense);

    dense = NULL;


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


    if(dense)
    {
        lis_free(
            dense);
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
