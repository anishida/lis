/* LIS_MATRIX / LIS_VECTOR adapter for the native symmetric eigensolver.
 * Copyright (C) 2026, contributors to LIS.
 * Redistribution and use are governed by the LIS project BSD license.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include <limits.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "lislib.h"
#include "lis_krylov.h"

typedef struct {
    LIS_MATRIX A, B, shifted;
    LIS_SOLVER linear_solver;
    LIS_PRECON linear_precon;
    LIS_VECTOR input, output, rhs, solution;
    LIS_INT last_error;
} LIS_KRYLOV_ADAPTER;

static int
lis_ekrylov_apply_matrix(
        LIS_KRYLOV_ADAPTER *adapter,
        LIS_MATRIX matrix,
        int n,
        const double *x,
        double *y)
{
    LIS_INT err;
    if( matrix==NULL || n!=adapter->A->n ) return -1;
    memcpy(adapter->input->value,x,(size_t)n*sizeof(double));
    err = lis_matvec(matrix,adapter->input,adapter->output);
    if( err )
    {
        adapter->last_error = err;
        return -1;
    }
    memcpy(y,adapter->output->value,(size_t)n*sizeof(double));
    return 0;
}

static int
lis_ekrylov_apply_a(void *state,int n,const double *x,double *y)
{
    LIS_KRYLOV_ADAPTER *ctx = (LIS_KRYLOV_ADAPTER *)state;
    return lis_ekrylov_apply_matrix(ctx,ctx->A,n,x,y);
}

static int
lis_ekrylov_apply_b(void *state,int n,const double *x,double *y)
{
    LIS_KRYLOV_ADAPTER *ctx = (LIS_KRYLOV_ADAPTER *)state;
    return lis_ekrylov_apply_matrix(ctx,ctx->B,n,x,y);
}

static int
lis_ekrylov_solve(void *state,int n,const double *x,double *y)
{
    LIS_KRYLOV_ADAPTER *ctx = (LIS_KRYLOV_ADAPTER *)state;
    LIS_INT status;
    if( ctx->shifted==NULL || n!=ctx->A->n ) return -1;
    memcpy(ctx->rhs->value,x,(size_t)n*sizeof(double));
    lis_vector_set_all(0.0,ctx->solution);
    status = lis_esolver_inner_solve_attempt(
            ctx->shifted,ctx->rhs,ctx->solution,
            ctx->linear_solver,ctx->linear_precon);
    if( status!=LIS_SUCCESS )
    {
        ctx->last_error = status;
        return -1;
    }
    memcpy(y,ctx->solution->value,(size_t)n*sizeof(double));
    return 0;
}

#undef __FUNC__
#define __FUNC__ "lis_ekrylov_check_params"
LIS_INT lis_ekrylov_check_params(LIS_ESOLVER esolver)
{
    LIS_INT ss, nev;
    if( esolver==NULL ) return LIS_ERR_ILL_ARG;
    ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
    nev = esolver->nev>0 ? esolver->nev : ss;
    if( ss<1 || nev<1 || nev>ss || esolver->krylov_ncv<0 ||
        (esolver->krylov_ncv>0 && esolver->krylov_ncv<=nev) )
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
                   "invalid native Krylov eigenpair/subspace configuration\n");
        return LIS_ERR_ILL_ARG;
    }
    if( esolver->options[LIS_EOPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE )
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                   "native Krylov is currently real double precision only\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
#ifdef _COMPLEX
    LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
               "native Krylov requires real arithmetic\n");
    return LIS_ERR_NOT_IMPLEMENTED;
#else
    if( sizeof(LIS_SCALAR)!=sizeof(double) )
        return LIS_ERR_NOT_IMPLEMENTED;
#endif
    return LIS_SUCCESS;
}

LIS_INT lis_ekrylov_malloc_work(LIS_ESOLVER esolver)
{
    esolver->work = NULL;
    esolver->worklen = 0;
    return LIS_SUCCESS;
}

static LIS_INT
lis_ekrylov_status(int status,LIS_KRYLOV_ADAPTER *ctx)
{
    if( status==LIS_KRYLOV_SUCCESS ) return LIS_SUCCESS;
    if( status==LIS_KRYLOV_NO_CONVERGENCE ) return LIS_MAXITER;
    if( status==LIS_KRYLOV_BAD_ARGUMENT ) return LIS_ERR_ILL_ARG;
    if( status==LIS_KRYLOV_NO_MEMORY ) return LIS_ERR_OUT_OF_MEMORY;
    if( status==LIS_KRYLOV_CALLBACK_ERROR && ctx->last_error!=LIS_SUCCESS )
        return ctx->last_error;
    return LIS_BREAKDOWN;
}

#undef __FUNC__
#define __FUNC__ "lis_ekrylov"
LIS_INT lis_ekrylov(LIS_ESOLVER esolver)
{
    LIS_KRYLOV_ADAPTER ctx;
    lis_krylov_problem p;
    lis_krylov_options options;
    lis_krylov_stats stats;
    LIS_INT err,ss,nev,i,return_code;
    LIS_SCALAR shift;
    double *eigenvalues=NULL,*eigenvectors=NULL,*residuals=NULL;
    size_t vector_count;
    int status;

    memset(&ctx,0,sizeof(ctx));
    ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
    nev = esolver->nev>0 ? esolver->nev : ss;
    ctx.A = esolver->A;
    ctx.B = esolver->B;
    err = lis_ekrylov_check_params(esolver);
    if( err ) return err;
    if( esolver->krylov_allow_semidefinite &&
        (ctx.B==NULL ||
         esolver->transform==LIS_EIGEN_TRANSFORM_BUCKLING) )
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
                   "semidefinite metric mode requires the generalized B metric\n");
        return LIS_ERR_ILL_ARG;
    }
    if( ctx.A->nprocs!=1 || ctx.A->n!=ctx.A->gn || ctx.A->gn>INT_MAX ||
        nev>=ctx.A->gn ||
        (esolver->krylov_ncv>0 && esolver->krylov_ncv>ctx.A->gn) ||
        (ctx.B!=NULL &&
         (ctx.B->nprocs!=1 || ctx.B->n!=ctx.B->gn ||
          ctx.B->n!=ctx.A->n)) )
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                   "native Krylov requires serial, matching matrices\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }


    shift = esolver->ishift!=0.0 ? esolver->ishift :
        esolver->params[LIS_EPARAMS_SHIFT-LIS_EOPTIONS_LEN];
    memset(&p,0,sizeof(p));
    p.n = (int)ctx.A->n;
    p.sigma = (double)shift;
    p.ctx = &ctx;
    p.apply_a = lis_ekrylov_apply_a;
    p.apply_b = lis_ekrylov_apply_b;
    p.solve = lis_ekrylov_solve;
    p.transform = ctx.B==NULL ? LIS_KRYLOV_STANDARD :
        esolver->transform==LIS_EIGEN_TRANSFORM_BUCKLING ?
        LIS_KRYLOV_SHIFT_INVERT_A : LIS_KRYLOV_SHIFT_INVERT_B;
    if( ctx.B==NULL && shift!=0.0 ) return LIS_ERR_ILL_ARG;

    lis_krylov_options_default(&options);
    options.nev = (int)nev;
    options.ncv = (int)esolver->krylov_ncv;
    options.maxiter = (int)esolver->options[LIS_EOPTIONS_MAXITER];
    options.tol = (double)esolver->params[
        LIS_EPARAMS_RESID-LIS_EOPTIONS_LEN];
    options.residual_scale = (double)esolver->krylov_residual_scale;
    options.inner_tol = (double)esolver->krylov_inner_tol;
    options.which = esolver->krylov_which;
    options.initial_vector = (const double *)esolver->x->value;
    options.allow_semidefinite_metric = (int)esolver->krylov_allow_semidefinite;

    /* Allocate the result arrays once; callbacks do not allocate or
     * reassemble matrices in the Krylov iteration. */
    vector_count = (size_t)p.n*(size_t)nev;
    if( vector_count>SIZE_MAX/sizeof(double) ) return LIS_ERR_ILL_ARG;
    eigenvalues = (double *)lis_malloc((size_t)nev*sizeof(double),
                                      "lis_ekrylov::eigenvalues");
    eigenvectors = (double *)lis_malloc(vector_count*sizeof(double),
                                       "lis_ekrylov::eigenvectors");
    residuals = (double *)lis_malloc((size_t)nev*sizeof(double),
                                    "lis_ekrylov::residuals");
    if( !eigenvalues || !eigenvectors || !residuals )
    {
        err = LIS_ERR_OUT_OF_MEMORY;
        goto cleanup;
    }
    err = lis_vector_duplicate(ctx.A,&ctx.input);
    if( err ) goto cleanup;
    err = lis_vector_duplicate(ctx.A,&ctx.output);
    if( err ) goto cleanup;
    err = lis_vector_duplicate(ctx.A,&ctx.rhs);
    if( err ) goto cleanup;
    err = lis_vector_duplicate(ctx.A,&ctx.solution);
    if( err ) goto cleanup;

    if( ctx.B!=NULL )
    {
        err = lis_matrix_create_operator(1.0,ctx.A,-shift,ctx.B,
                                         &ctx.shifted);
        if( err ) goto cleanup;
        err = lis_solver_create(&ctx.linear_solver);
        if( err ) goto cleanup;
        err = lis_solver_set_option("-i bicg -p none",ctx.linear_solver);
        if( err ) goto cleanup;
        if( esolver->linear_options && esolver->linear_options[0] )
        {
            err = lis_solver_set_option(esolver->linear_options,
                                        ctx.linear_solver);
            if( err ) goto cleanup;
        }
        err = lis_solver_set_optionC(ctx.linear_solver);
        if( err ) goto cleanup;
        err = lis_esolver_rebuild_precon(
                ctx.shifted,ctx.linear_solver,&ctx.linear_precon);
        if( err ) goto cleanup;
    }

    status = lis_krylov_symmetric(&p,&options,
                                  eigenvalues,eigenvectors,residuals,&stats);
    return_code = lis_ekrylov_status(status,&ctx);
    for( i=0;i<stats.nreturned && i<nev;i++ )
    {
        err = lis_vector_duplicate(ctx.A,&esolver->evector[i]);
        if( err ) goto cleanup;
        memcpy(esolver->evector[i]->value,
               eigenvectors+(size_t)i*p.n,
               (size_t)p.n*sizeof(double));
        esolver->evalue[i] = (LIS_SCALAR)eigenvalues[i];
        esolver->resid[i] = (LIS_REAL)residuals[i];
        esolver->iter[i] = (LIS_INT)stats.iterations;
    }
    esolver->nev = (LIS_INT)stats.nreturned;
    if( stats.nreturned>0 )
        memcpy(esolver->x->value,eigenvectors,(size_t)p.n*sizeof(double));
    esolver->retcode = return_code;
    err = return_code;

cleanup:
    if( ctx.linear_precon ) lis_precon_destroy(ctx.linear_precon);
    if( ctx.linear_solver ) lis_solver_destroy(ctx.linear_solver);
    if( ctx.shifted ) lis_matrix_destroy(ctx.shifted);
    if( ctx.input ) lis_vector_destroy(ctx.input);
    if( ctx.output ) lis_vector_destroy(ctx.output);
    if( ctx.rhs ) lis_vector_destroy(ctx.rhs);
    if( ctx.solution ) lis_vector_destroy(ctx.solution);
    if( eigenvalues ) lis_free(eigenvalues);
    if( eigenvectors ) lis_free(eigenvectors);
    if( residuals ) lis_free(residuals);
    return err;
}
