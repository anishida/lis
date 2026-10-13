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

    LIS_VECTOR rhs;
    LIS_VECTOR x;
};


struct LIS_MFAMG_STRUCT
{
    LIS_INT nlevels;
    LIS_INT level_capacity;
    struct LIS_MFAMG_LEVEL *levels;

    LIS_INT coarse_n;

    LIS_SCALAR *coarse;
    LIS_SCALAR *lu;
    LIS_INT *pivots;

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


    if(state->coarse)
        lis_free(state->coarse);

    if(state->lu)
        lis_free(state->lu);

    if(state->pivots)
        lis_free(state->pivots);


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
        state->nlevels<2)
    {
        return LIS_ERR_ILL_ARG;
    }


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


    while(
        level_id<state->level_capacity-1
        &&
        state->levels[level_id].n
            >
        LIS_MFAMG_COARSE_LIMIT)
    {
        level =
            &state->levels[level_id];

        next_level =
            &state->levels[level_id+1];


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


        if(level->nc>=level->n)
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


        lis_free(dense);
        lis_free(vals);
        lis_free(cols);

        dense = NULL;
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
    LIS_INT j;
    LIS_INT p;

    size_t pos;


    if(
        A==NULL
        ||
        state==NULL
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->n!=state->coarse_n
        ||
        state->coarse_n<=0
        ||
        state->lu==NULL
        ||
        state->pivots==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    memset(
        state->lu,
        0,
        (size_t)state->coarse_n
            *
        (size_t)state->coarse_n
            *
        sizeof(LIS_SCALAR));


    for(i=0;i<A->n;i++)
    {
        for(
            p=A->ptr[i];
            p<A->ptr[i+1];
            p++)
        {
            j =
                A->index[p]
                    -
                A->is;

            if(
                j<0
                ||
                j>=state->coarse_n)
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "MFAMG internal coarse CSR contains a nonlocal column\n");

                return LIS_ERR_ILL_ARG;
            }

            pos =
                (size_t)i
                    *
                (size_t)state->coarse_n
                    +
                (size_t)j;

            state->lu[pos] +=
                A->value[p];
        }
    }


    return
        lis_mfamg_dense_factor(
            state->coarse_n,
            state->lu,
            state->pivots);
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

    LIS_INT capabilities;
    LIS_INT required;

    LIS_INT max_row_nnz;
    LIS_INT row_nnz;

    LIS_INT *cols;
    LIS_SCALAR *vals;

    LIS_INT i;
    LIS_INT err;

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
            "MFAMG supports one MPI process\n");

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


    if(
        solver->near_nullspace_dim!=1
        ||
        solver->near_nullspace==NULL
        ||
        solver->near_nullspace[0]==NULL)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "MFAMG requires exactly one near-nullspace vector\n");

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
            "MFAMG supports CSR and LIS_MATRIX_USER\n");

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

    if(state->coarse==NULL)
    {
        lis_free(vals);
        lis_free(cols);

        LIS_SETERR_MEM(
            coarse_entries
                *
            sizeof(LIS_SCALAR));

        lis_mfamg_state_free(state);

        return LIS_OUT_OF_MEMORY;
    }


    err =
        lis_mfamg_build_coarse(
            A,
            level,
            state->coarse,
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
        lis_mfamg_create_csr_from_dense(
            A->comm,
            level->nc,
            state->coarse,
            coarse_level);

    if(err)
    {
        lis_mfamg_state_free(state);
        return err;
    }


    lis_free(
        state->coarse);

    state->coarse = NULL;


    err =
        lis_mfamg_extend_hierarchy(
            state,
            solver->near_nullspace[0]);

    if(err)
    {
        lis_mfamg_state_free(state);
        return err;
    }


    err =
        lis_mfamg_prepare_level_work(
            level);

    if(err)
    {
        lis_mfamg_state_free(state);
        return err;
    }


    err =
        lis_mfamg_prepare_cycle_vectors(
            state);

    if(err)
    {
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

        lis_mfamg_state_free(state);

        return err;
    }


    terminal_level =
        &state->levels[
            state->nlevels-1];


    state->coarse_n =
        terminal_level->n;


    if(
        state->coarse_n<=0
        ||
        (size_t)state->coarse_n
        >
        ((size_t)-1)
        /
        (size_t)state->coarse_n)
    {
        LIS_SETERR(
            LIS_OUT_OF_MEMORY,
            "MFAMG terminal matrix dimension overflow\n");

        lis_mfamg_state_free(state);

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

        lis_mfamg_state_free(state);

        return LIS_OUT_OF_MEMORY;
    }


    state->lu =
        (LIS_SCALAR *)lis_malloc(
            coarse_entries
                *
            sizeof(LIS_SCALAR),
            "lis_precon_create_mfamg::terminal_lu");

    state->pivots =
        (LIS_INT *)lis_malloc(
            (size_t)state->coarse_n
                *
            sizeof(LIS_INT),
            "lis_precon_create_mfamg::terminal_pivots");


    if(
        state->lu==NULL
        ||
        state->pivots==NULL)
    {
        LIS_SETERR_MEM(
            coarse_entries
                *
            sizeof(LIS_SCALAR)
            +
            (size_t)state->coarse_n
                *
            sizeof(LIS_INT));

        lis_mfamg_state_free(state);

        return LIS_OUT_OF_MEMORY;
    }


    err =
        lis_mfamg_factor_csr_terminal(
            terminal_level->A,
            state);

    if(err)
    {
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
        b->n!=state->levels[0].n
        ||
        x->n!=state->levels[0].n)
    {
        return LIS_ERR_ILL_ARG;
    }


    err =
        lis_mfamg_v_cycle(
            state,
            0,
            b,
            x);

    if(err)
        return err;


    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}
