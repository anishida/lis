/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Near-nullspace coarse-space support.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>

#include "lislib.h"


static LIS_INT
lis_solver_near_nullspace_matrix_compatible(
    LIS_VECTOR v,
    LIS_MATRIX A)
{
    LIS_INT i;

#ifdef USE_MPI
    LIS_INT mpi_err;
    LIS_INT comm_result;
#endif


    if(
        v==NULL
        ||
        A==NULL)
    {
        return LIS_FALSE;
    }


    if(
        v->gn!=A->gn
        ||
        v->n!=A->n
        ||
        v->np!=A->np
        ||
        v->pad!=A->pad
        ||
        v->origin!=A->origin
        ||
        v->my_rank!=A->my_rank
        ||
        v->nprocs!=A->nprocs
        ||
        v->is!=A->is
        ||
        v->ie!=A->ie)
    {
        return LIS_FALSE;
    }


#ifdef USE_MPI

    mpi_err = MPI_Comm_compare(
        v->comm,
        A->comm,
        &comm_result);


    if(
        mpi_err!=MPI_SUCCESS
        ||
        (
            comm_result!=MPI_IDENT
            &&
            comm_result!=MPI_CONGRUENT
        ))
    {
        return LIS_FALSE;
    }


    if(
        v->ranges==NULL
        ||
        A->ranges==NULL)
    {
        if(v->ranges!=A->ranges)
            return LIS_FALSE;
    }
    else
    {
        for(i=0;i<=A->nprocs;i++)
        {
            if(v->ranges[i]!=A->ranges[i])
                return LIS_FALSE;
        }
    }

#else

    (void)i;

    if(v->comm!=A->comm)
        return LIS_FALSE;

#endif


    return LIS_TRUE;
}


static void
lis_solver_near_nullspace_coarse_free(
    LIS_NEAR_NULLSPACE_COARSE coarse)
{
    LIS_INT i;


    if(coarse==NULL)
        return;


    if(coarse->Z)
    {
        for(i=0;i<coarse->dim;i++)
        {
            if(coarse->Z[i])
                lis_vector_destroy(coarse->Z[i]);
        }

        lis_free(coarse->Z);
    }


    if(coarse->AZ)
    {
        for(i=0;i<coarse->dim;i++)
        {
            if(coarse->AZ[i])
                lis_vector_destroy(coarse->AZ[i]);
        }

        lis_free(coarse->AZ);
    }


    if(coarse->work)
        lis_vector_destroy(coarse->work);


    if(coarse->E)
        lis_free(coarse->E);

    if(coarse->LU)
        lis_free(coarse->LU);

    if(coarse->pivots)
        lis_free(coarse->pivots);

    if(coarse->rhs)
        lis_free(coarse->rhs);

    if(coarse->coeff)
        lis_free(coarse->coeff);


    lis_free(coarse);
}


LIS_INT
lis_solver_near_nullspace_coarse_destroy(
    LIS_SOLVER solver)
{
    if(solver==NULL)
        return LIS_ERR_ILL_ARG;


    lis_solver_near_nullspace_coarse_free(
        solver->near_nullspace_coarse);


    solver->near_nullspace_coarse = NULL;


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_near_nullspace_factor"
static LIS_INT
lis_solver_near_nullspace_factor(
    LIS_INT n,
    LIS_SCALAR *lu,
    LIS_INT *pivots)
{
    LIS_INT i,j,k,p;

    LIS_SCALAR temp;

    LIS_REAL mag;
    LIS_REAL best;


    for(i=0;i<n*n;i++)
    {
        mag = fabs(lu[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "near-nullspace coarse matrix contains non-finite values\n");

            return LIS_BREAKDOWN;
        }
    }


    for(k=0;k<n;k++)
    {
        p = k;

        best =
            fabs(lu[k*n+k]);


        for(i=k+1;i<n;i++)
        {
            mag =
                fabs(lu[i*n+k]);

            if(mag>best)
            {
                best = mag;
                p = i;
            }
        }


        /*
         * No absolute near-zero threshold is used here.
         *
         * A near-nullspace coarse operator may legitimately
         * have a very small absolute scale. Only an exact
         * zero or non-finite pivot is rejected.
         */
        if(
            best==0.0
            ||
            best!=best
            ||
            best>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "near-nullspace coarse matrix is singular\n");

            return LIS_BREAKDOWN;
        }


        pivots[k] = p;


        if(p!=k)
        {
            for(j=0;j<n;j++)
            {
                temp =
                    lu[k*n+j];

                lu[k*n+j] =
                    lu[p*n+j];

                lu[p*n+j] =
                    temp;
            }
        }


        for(i=k+1;i<n;i++)
        {
            lu[i*n+k] /=
                lu[k*n+k];


            for(j=k+1;j<n;j++)
            {
                lu[i*n+j] -=
                    lu[i*n+k]
                    *
                    lu[k*n+j];
            }
        }
    }


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_near_nullspace_coarse_solve"
LIS_INT
lis_solver_near_nullspace_coarse_solve(
    LIS_SOLVER solver,
    const LIS_SCALAR *rhs,
    LIS_SCALAR *x)
{
    LIS_NEAR_NULLSPACE_COARSE coarse;

    LIS_INT i,j,k,p,n;

    LIS_SCALAR temp;

    LIS_REAL mag;


    if(
        solver==NULL
        ||
        rhs==NULL
        ||
        x==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    coarse =
        solver->near_nullspace_coarse;


    if(
        coarse==NULL
        ||
        !coarse->ready)
    {
        return LIS_ERR_ILL_ARG;
    }


    n =
        coarse->dim;


    for(i=0;i<n;i++)
        x[i] = rhs[i];


    /*
     * Apply the row permutations generated by
     * partial pivoting during factorization.
     */
    for(k=0;k<n;k++)
    {
        p =
            coarse->pivots[k];

        if(p!=k)
        {
            temp = x[k];
            x[k] = x[p];
            x[p] = temp;
        }
    }


    /*
     * Forward substitution with unit diagonal L.
     */
    for(i=0;i<n;i++)
    {
        for(j=0;j<i;j++)
        {
            x[i] -=
                coarse->LU[i*n+j]
                *
                x[j];
        }
    }


    /*
     * Back substitution with U.
     */
    for(i=n-1;i>=0;i--)
    {
        for(j=i+1;j<n;j++)
        {
            x[i] -=
                coarse->LU[i*n+j]
                *
                x[j];
        }


        x[i] /=
            coarse->LU[i*n+i];


        mag =
            fabs(x[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "near-nullspace coarse solve became non-finite\n");

            return LIS_BREAKDOWN;
        }
    }


    return LIS_SUCCESS;
}



#undef __FUNC__
#define __FUNC__ "lis_solver_near_nullspace_apply"
LIS_INT
lis_solver_near_nullspace_apply(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_PSOLVE_XXX psolve)
{
    LIS_NEAR_NULLSPACE_COARSE coarse;

    LIS_INT i,l,k,n,err;

    LIS_SCALAR sum;


    if(
        solver==NULL
        ||
        b==NULL
        ||
        x==NULL
        ||
        psolve==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    coarse =
        solver->near_nullspace_coarse;


    if(
        coarse==NULL
        ||
        !coarse->ready
        ||
        coarse->work==NULL
        ||
        coarse->rhs==NULL
        ||
        coarse->coeff==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(solver->A==NULL)
        return LIS_ERR_ILL_ARG;


    k =
        coarse->dim;

    n =
        solver->A->n;


#ifdef _OPENMP
#pragma omp parallel for private(l,sum)
#endif
    for(i=0;i<k;i++)
    {
        sum = 0.0;

        for(l=0;l<n;l++)
        {
            sum +=
                conj(coarse->Z[i]->value[l])
                *
                b->value[l];
        }

        coarse->rhs[i] =
            sum;
    }


#ifdef USE_MPI

    MPI_Allreduce(
        coarse->rhs,
        coarse->coeff,
        (int)k,
        LIS_MPI_SCALAR,
        MPI_SUM,
        solver->A->comm);

#else

    for(i=0;i<k;i++)
    {
        coarse->coeff[i] =
            coarse->rhs[i];
    }

#endif


    err =
        lis_solver_near_nullspace_coarse_solve(
            solver,
            coarse->coeff,
            coarse->rhs);

    if(err)
        return err;


    err =
        lis_vector_copy(
            b,
            coarse->work);

    if(err)
        return err;


    for(i=0;i<k;i++)
    {
        err =
            lis_vector_axpy(
                -coarse->rhs[i],
                coarse->AZ[i],
                coarse->work);

        if(err)
            return err;
    }


    /*
     * Apply the already-resolved raw preconditioner.
     * Do not call lis_psolve() here.
     */
    err =
        psolve(
            solver,
            coarse->work,
            x);

    if(err)
        return err;


    for(i=0;i<k;i++)
    {
        err =
            lis_vector_axpy(
                coarse->rhs[i],
                coarse->Z[i],
                x);

        if(err)
            return err;
    }


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_near_nullspace_coarse_setup"
LIS_INT
lis_solver_near_nullspace_coarse_setup(
    LIS_SOLVER solver,
    LIS_MATRIX A,
    LIS_INT scale)
{
    LIS_NEAR_NULLSPACE_COARSE coarse;

    LIS_INT i,j,l,k,n,err;

    LIS_SCALAR sum;


    if(
        solver==NULL
        ||
        A==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    /*
     * A changed system always invalidates AZ and E.
     * Never cache these by matrix pointer: callers may
     * modify an existing matrix in place between solves.
     */
    lis_solver_near_nullspace_coarse_destroy(
        solver);


    k =
        solver->near_nullspace_dim;


    if(k==0)
        return LIS_SUCCESS;


    if(
        solver->near_nullspace==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    /*
     * Stage 7B supports ordinary double-precision
     * coarse state. Quad/switch integration is handled
     * separately rather than silently truncating data.
     */
    if(
        solver->precision!=LIS_PRECISION_DOUBLE)
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "near-nullspace coarse setup currently requires double precision\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }


    if(
        scale!=LIS_SCALE_NONE
        &&
        scale!=LIS_SCALE_JACOBI
        &&
        scale!=LIS_SCALE_SYMM_DIAG)
    {
        return LIS_ERR_ILL_ARG;
    }


    for(i=0;i<k;i++)
    {
        if(
            !lis_solver_near_nullspace_matrix_compatible(
                solver->near_nullspace[i],
                A))
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "near-nullspace basis is incompatible with the effective operator\n");

            return LIS_ERR_ILL_ARG;
        }
    }


    if(
        scale==LIS_SCALE_SYMM_DIAG
        &&
        solver->d==NULL)
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "symmetric near-nullspace transformation requires scaling vector\n");

        return LIS_ERR_ILL_ARG;
    }


    coarse =
        (LIS_NEAR_NULLSPACE_COARSE)
        lis_malloc(
            sizeof(struct LIS_NEAR_NULLSPACE_COARSE_STRUCT),
            "lis_solver_near_nullspace_coarse_setup::coarse");


    if(coarse==NULL)
    {
        LIS_SETERR_MEM(
            sizeof(struct LIS_NEAR_NULLSPACE_COARSE_STRUCT));

        return LIS_OUT_OF_MEMORY;
    }


    coarse->dim = k;
    coarse->Z = NULL;
    coarse->AZ = NULL;
    coarse->E = NULL;
    coarse->LU = NULL;
    coarse->pivots = NULL;
    coarse->work = NULL;
    coarse->rhs = NULL;
    coarse->coeff = NULL;
    coarse->ready = LIS_FALSE;


    coarse->Z =
        (LIS_VECTOR *)
        lis_malloc(
            k*sizeof(LIS_VECTOR),
            "lis_solver_near_nullspace_coarse_setup::Z");


    if(coarse->Z)
    {
        for(i=0;i<k;i++)
            coarse->Z[i] = NULL;
    }


    coarse->AZ =
        (LIS_VECTOR *)
        lis_malloc(
            k*sizeof(LIS_VECTOR),
            "lis_solver_near_nullspace_coarse_setup::AZ");


    if(coarse->AZ)
    {
        for(i=0;i<k;i++)
            coarse->AZ[i] = NULL;
    }


    coarse->E =
        (LIS_SCALAR *)
        lis_malloc(
            k*k*sizeof(LIS_SCALAR),
            "lis_solver_near_nullspace_coarse_setup::E");


    coarse->LU =
        (LIS_SCALAR *)
        lis_malloc(
            k*k*sizeof(LIS_SCALAR),
            "lis_solver_near_nullspace_coarse_setup::LU");


    coarse->pivots =
        (LIS_INT *)
        lis_malloc(
            k*sizeof(LIS_INT),
            "lis_solver_near_nullspace_coarse_setup::pivots");


    coarse->rhs =
        (LIS_SCALAR *)
        lis_malloc(
            k*sizeof(LIS_SCALAR),
            "lis_solver_near_nullspace_coarse_setup::rhs");


    coarse->coeff =
        (LIS_SCALAR *)
        lis_malloc(
            k*sizeof(LIS_SCALAR),
            "lis_solver_near_nullspace_coarse_setup::coeff");


    if(
        coarse->Z==NULL
        ||
        coarse->AZ==NULL
        ||
        coarse->E==NULL
        ||
        coarse->LU==NULL
        ||
        coarse->pivots==NULL
        ||
        coarse->rhs==NULL
        ||
        coarse->coeff==NULL)
    {
        lis_solver_near_nullspace_coarse_free(
            coarse);

        return LIS_OUT_OF_MEMORY;
    }


    err =
        lis_vector_duplicate(
            A,
            &coarse->work);

    if(err)
        goto fail;


    n =
        A->n;


    for(i=0;i<k;i++)
    {
        err = lis_vector_duplicate(
            A,
            &coarse->Z[i]);

        if(err)
            goto fail;


        err = lis_vector_duplicate(
            A,
            &coarse->AZ[i]);

        if(err)
            goto fail;


        if(scale==LIS_SCALE_SYMM_DIAG)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(l=0;l<n;l++)
            {
                /*
                 * For D*A*D:
                 *
                 *     x = D*y
                 *
                 * so a physical near-null vector z is
                 * represented in Krylov coordinates as
                 *
                 *     z_eff = D^-1 z.
                 */
                coarse->Z[i]->value[l] =
                    solver->near_nullspace[i]->value[l]
                    /
                    solver->d->value[l];
            }
        }
        else
        {
            err = lis_vector_copy(
                solver->near_nullspace[i],
                coarse->Z[i]);

            if(err)
                goto fail;
        }


        err = lis_matvec(
            A,
            coarse->Z[i],
            coarse->AZ[i]);

        if(err)
            goto fail;
    }


    /*
     * Form the entire local block
     *
     *     E_local = Z^H * A * Z
     *
     * before any MPI communication.
     */
#ifdef _OPENMP
#pragma omp parallel for private(j,l,sum)
#endif
    for(i=0;i<k;i++)
    {
        for(j=0;j<k;j++)
        {
            sum = 0.0;


            for(l=0;l<n;l++)
            {
                sum +=
                    conj(coarse->Z[i]->value[l])
                    *
                    coarse->AZ[j]->value[l];
            }


            coarse->E[i*k+j] =
                sum;
        }
    }


#ifdef USE_MPI

    /*
     * One global synchronization for the whole
     * k-by-k coarse matrix, instead of k*k separate
     * lis_vector_dot() reductions.
     */
    MPI_Allreduce(
        coarse->E,
        coarse->LU,
        (int)(k*k),
        LIS_MPI_SCALAR,
        MPI_SUM,
        A->comm);


    for(i=0;i<k*k;i++)
    {
        coarse->E[i] =
            coarse->LU[i];
    }

#else

    for(i=0;i<k*k;i++)
    {
        coarse->LU[i] =
            coarse->E[i];
    }

#endif


    /*
     * In MPI builds LU currently contains the global E.
     * In serial builds it was copied directly above.
     */
#ifdef USE_MPI
    for(i=0;i<k*k;i++)
    {
        coarse->LU[i] =
            coarse->E[i];
    }
#endif


    err =
        lis_solver_near_nullspace_factor(
            k,
            coarse->LU,
            coarse->pivots);


    if(err)
        goto fail;


    coarse->ready =
        LIS_TRUE;


    solver->near_nullspace_coarse =
        coarse;


    return LIS_SUCCESS;


fail:

    lis_solver_near_nullspace_coarse_free(
        coarse);


    solver->near_nullspace_coarse =
        NULL;


    return err;
}
