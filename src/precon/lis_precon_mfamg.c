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


struct LIS_MFAMG_LEVEL
{
    LIS_MATRIX A;
    LIS_INT owns_A;

    LIS_INT n;
    LIS_INT nc;

    LIS_INT *aggregate;
    LIS_SCALAR *prolong;

    LIS_VECTOR diagonal;
    LIS_VECTOR ax;
};


struct LIS_MFAMG_STRUCT
{
    LIS_INT nlevels;
    struct LIS_MFAMG_LEVEL *levels;

    LIS_INT coarse_n;

    LIS_SCALAR *coarse;
    LIS_SCALAR *lu;
    LIS_INT *pivots;

    LIS_SCALAR *coarse_rhs;
    LIS_SCALAR *coarse_x;

    LIS_REAL omega;
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

            if(level->prolong)
                lis_free(level->prolong);

            if(level->diagonal)
                lis_vector_destroy(level->diagonal);

            if(level->ax)
                lis_vector_destroy(level->ax);

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
            /*
             * Stage 3 is serial, so the internal CSR column index
             * maps directly to the global range starting at A->is.
             */
            cols[k] =
                A->index[p]
                +
                A->is;

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
    struct LIS_MFAMG_STRUCT *state)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT p;

    LIS_REAL mag;

    LIS_SCALAR tmp;
    LIS_SCALAR sum;


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


    for(i=state->coarse_n-1;i>=0;i--)
    {
        sum =
            state->coarse_x[i];

        for(j=i+1;j<state->coarse_n;j++)
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


    return LIS_SUCCESS;
}


static LIS_INT
lis_mfamg_build_aggregates(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level,
    LIS_INT max_row_nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT nnz;
    LIS_INT best;
    LIS_INT aggregate;

    LIS_REAL best_mag;
    LIS_REAL mag;


    for(i=0;i<level->n;i++)
        level->aggregate[i] = -1;


    aggregate = 0;

    for(i=0;i<level->n;i++)
    {
        if(level->aggregate[i]>=0)
            continue;

        level->aggregate[i] =
            aggregate;

        best = -1;
        best_mag = 0.0;

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
                j>=level->n
                ||
                j==i
                ||
                level->aggregate[j]>=0)
            {
                continue;
            }

            mag =
                fabs(
                    vals[k]);

            if(mag>best_mag)
            {
                best_mag = mag;
                best = j;
            }
        }

        if(best>=0)
        {
            level->aggregate[best] =
                aggregate;
        }

        aggregate++;
    }

    level->nc =
        aggregate;

    if(level->nc<=0)
        return LIS_BREAKDOWN;

    return LIS_SUCCESS;
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


#undef __FUNC__
#define __FUNC__ "lis_mfamg_build_coarse"
static LIS_INT
lis_mfamg_build_coarse(
    LIS_MATRIX A,
    struct LIS_MFAMG_LEVEL *level,
    struct LIS_MFAMG_STRUCT *state,
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
        state->coarse,
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

            state->coarse[
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


    memcpy(
        state->lu,
        state->coarse,
        (size_t)level->nc
            *
        (size_t)level->nc
            *
        sizeof(LIS_SCALAR));


    return
        lis_mfamg_dense_factor(
            level->nc,
            state->lu,
            state->pivots);
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

    LIS_MATRIX A;

    LIS_INT capabilities;
    LIS_INT required;

    LIS_INT max_row_nnz;
    LIS_INT row_nnz;

    LIS_INT *cols;
    LIS_SCALAR *vals;

    LIS_INT i;
    LIS_INT err;

    LIS_REAL diagonal_mag;

    size_t coarse_entries;


    LIS_DEBUG_FUNC_IN;


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


#ifdef USE_MPI
    if(A->nprocs!=1)
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG Stage 3 supports one MPI process\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }
#endif


    if(
        solver->options[LIS_OPTIONS_PRECISION]
        !=
        LIS_PRECISION_DOUBLE)
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG Stage 3 supports double precision\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }


    if(
        solver->options[LIS_OPTIONS_SCALE]
        !=
        LIS_SCALE_NONE)
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG Stage 3 requires scale=none\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }


    if(
        solver->near_nullspace_dim!=1
        ||
        solver->near_nullspace==NULL
        ||
        solver->near_nullspace[0]==NULL)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG Stage 3 requires exactly one near-nullspace vector\n");

        return LIS_ERR_ILL_ARG;
    }


    if(
        solver->near_nullspace[0]->n!=A->n
        ||
        solver->near_nullspace[0]->gn!=A->gn
        ||
        solver->near_nullspace[0]->is!=A->is
        ||
        solver->near_nullspace[0]->ie!=A->ie)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG near-nullspace vector is incompatible with the matrix\n");

        return LIS_ERR_ILL_ARG;
    }


    if(
        A->matrix_type!=LIS_MATRIX_CSR
        &&
        A->matrix_type!=LIS_MATRIX_USER)
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "MFAMG Stage 3 supports CSR and LIS_MATRIX_USER\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }


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
        1;

    state->levels =
        (struct LIS_MFAMG_LEVEL *)
        lis_malloc(
            sizeof(struct LIS_MFAMG_LEVEL),
            "lis_precon_create_mfamg::levels");

    if(state->levels==NULL)
    {
        LIS_SETERR_MEM(
            sizeof(struct LIS_MFAMG_LEVEL));

        lis_mfamg_state_free(state);

        return LIS_OUT_OF_MEMORY;
    }


    memset(
        state->levels,
        0,
        sizeof(struct LIS_MFAMG_LEVEL));


    level =
        &state->levels[0];

    level->A =
        A;

    level->owns_A =
        LIS_FALSE;

    level->n =
        A->n;

    state->omega =
        (LIS_REAL)(2.0/3.0);


    level->aggregate =
        (LIS_INT *)lis_malloc(
            (size_t)level->n
                *
            sizeof(LIS_INT),
            "lis_precon_create_mfamg::aggregate");

    level->prolong =
        (LIS_SCALAR *)lis_malloc(
            (size_t)level->n
                *
            sizeof(LIS_SCALAR),
            "lis_precon_create_mfamg::prolong");

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
            lis_mfamg_state_free(state);
            return err;
        }

        if(row_nnz<0)
        {
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
        lis_mfamg_state_free(state);

        return err;
    }


    err =
        lis_mfamg_build_prolongation(
            solver->near_nullspace[0],
            level);

    if(err)
    {
        lis_free(vals);
        lis_free(cols);
        lis_mfamg_state_free(state);

        return err;
    }


    state->coarse_n =
        level->nc;


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
        lis_mfamg_state_free(state);

        return LIS_OUT_OF_MEMORY;
    }


    state->coarse =
        (LIS_SCALAR *)lis_malloc(
            coarse_entries
                *
            sizeof(LIS_SCALAR),
            "lis_precon_create_mfamg::coarse");

    state->lu =
        (LIS_SCALAR *)lis_malloc(
            coarse_entries
                *
            sizeof(LIS_SCALAR),
            "lis_precon_create_mfamg::lu");

    state->pivots =
        (LIS_INT *)lis_malloc(
            (size_t)level->nc
                *
            sizeof(LIS_INT),
            "lis_precon_create_mfamg::pivots");

    state->coarse_rhs =
        (LIS_SCALAR *)lis_malloc(
            (size_t)level->nc
                *
            sizeof(LIS_SCALAR),
            "lis_precon_create_mfamg::coarse_rhs");

    state->coarse_x =
        (LIS_SCALAR *)lis_malloc(
            (size_t)level->nc
                *
            sizeof(LIS_SCALAR),
            "lis_precon_create_mfamg::coarse_x");


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
        lis_free(vals);
        lis_free(cols);
        lis_mfamg_state_free(state);

        return LIS_OUT_OF_MEMORY;
    }


    err =
        lis_mfamg_build_coarse(
            A,
            level,
            state,
            max_row_nnz,
            cols,
            vals);

    lis_free(vals);
    lis_free(cols);

    if(err)
    {
        lis_mfamg_state_free(state);
        return err;
    }


    err =
        lis_vector_duplicate(
            A,
            &level->diagonal);

    if(err)
    {
        lis_mfamg_state_free(state);
        return err;
    }


    err =
        lis_vector_duplicate(
            A,
            &level->ax);

    if(err)
    {
        lis_mfamg_state_free(state);
        return err;
    }


    err =
        lis_matrix_get_diagonal(
            A,
            level->diagonal);

    if(err)
    {
        lis_mfamg_state_free(state);
        return err;
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

            lis_mfamg_state_free(state);

            return LIS_BREAKDOWN;
        }
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
    struct LIS_MFAMG_LEVEL *level;

    LIS_INT i;
    LIS_INT a;
    LIS_INT err;

    LIS_SCALAR residual;

    LIS_SCALAR *bv;
    LIS_SCALAR *xv;
    LIS_SCALAR *av;
    LIS_SCALAR *diag;


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
        state->nlevels<=0
        ||
        state->levels==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    level =
        &state->levels[0];


    if(
        b->n!=level->n
        ||
        x->n!=level->n)
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


    /*
     * One weighted-Jacobi pre-smoothing step from zero.
     */
    for(i=0;i<level->n;i++)
    {
        xv[i] =
            state->omega
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
        state->coarse_rhs,
        0,
        (size_t)level->nc
            *
        sizeof(LIS_SCALAR));


    /*
     * Restriction r_c = P^H (b - A x).
     */
    for(i=0;i<level->n;i++)
    {
        residual =
            bv[i]
            -
            av[i];

        a =
            level->aggregate[i];

        state->coarse_rhs[a] +=
            conj(
                level->prolong[i])
            *
            residual;
    }


    err =
        lis_mfamg_dense_solve(
            state);

    if(err)
        return err;


    /*
     * Coarse correction x <- x + P e_c.
     */
    for(i=0;i<level->n;i++)
    {
        a =
            level->aggregate[i];

        xv[i] +=
            level->prolong[i]
            *
            state->coarse_x[a];
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
            state->omega
            *
            residual
            /
            diag[i];
    }


    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}
