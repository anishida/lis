/* Copyright (C) 2005 The Scalable Software Infrastructure Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
	#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
	#include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_esolver_init
 * lis_esolver_create
 * lis_esolver_destroy
 * lis_iesolver_destroy
 * lis_esolver_work_destroy
 * lis_esolver_set_option
 * lis_esolver_get_option
 * lis_esolve
 ************************************************/

LIS_ESOLVER_EXECUTE lis_esolver_execute[] = {
  NULL,
  lis_epi, lis_eii, lis_erqi, lis_ecg, lis_ecr,
  lis_esi, lis_eli, lis_eai,
  lis_egpi, lis_egii, lis_egrqi, lis_egcg, lis_egcr,
  lis_egsi, lis_egli, lis_egai,
  lis_efeast
};

#ifdef USE_QUAD_PRECISION
LIS_ESOLVER_EXECUTE lis_esolver_execute_quad[] = {
  NULL,
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
LIS_ESOLVER_EXECUTE lis_esolver_execute_switch[] = {
  NULL,
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
#endif

LIS_ESOLVER_CHECK_PARAMS lis_esolver_check_params[] = {
  NULL,
  lis_epi_check_params,  lis_eii_check_params, lis_erqi_check_params, 
  lis_ecg_check_params,  lis_ecr_check_params,
  lis_esi_check_params,  lis_eli_check_params, lis_eai_check_params,
  lis_egpi_check_params, lis_egii_check_params, lis_egrqi_check_params,
  lis_egcg_check_params, lis_egcr_check_params, lis_egsi_check_params,
  lis_egli_check_params, lis_egai_check_params,
  lis_efeast_check_params
};

LIS_ESOLVER_MALLOC_WORK lis_esolver_malloc_work[] = {
  NULL,
  lis_epi_malloc_work,  lis_eii_malloc_work, lis_erqi_malloc_work, 
  lis_ecg_malloc_work,  lis_ecr_malloc_work, 
  lis_esi_malloc_work,  lis_eli_malloc_work, lis_eai_malloc_work,
  lis_egpi_malloc_work, lis_egii_malloc_work, lis_egrqi_malloc_work,
  lis_egcg_malloc_work, lis_egcr_malloc_work, lis_egsi_malloc_work,
  lis_egli_malloc_work, lis_egai_malloc_work,
  lis_efeast_malloc_work
};

#define LIS_ESOLVER_OPTION_LEN		15
#define LIS_EPRINT_LEN			 4
#define LIS_TRUEFALSE_LEN		 2
#define LIS_ESTORAGE_LEN		11
#define LIS_PRECISION_LEN		 3

char *LIS_ESOLVER_OPTNAME[] = {
  "-emaxiter", "-etol", "-e", "-ss", "-m",
  "-shift", "-shift_im", "-eprint", "-initx_ones", "-ie",
  "-ige", "-estorage", "-estorage_block", "-ef", "-rval"
};

LIS_INT LIS_ESOLVER_OPTACT[] = {
LIS_EOPTIONS_MAXITER, LIS_EPARAMS_RESID, LIS_EOPTIONS_ESOLVER, 
LIS_EOPTIONS_SUBSPACE, LIS_EOPTIONS_MODE, LIS_EPARAMS_SHIFT,
LIS_EPARAMS_SHIFT_IM, LIS_EOPTIONS_OUTPUT, LIS_EOPTIONS_INITGUESS_ONES, 
LIS_EOPTIONS_INNER_ESOLVER, LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER,
LIS_EOPTIONS_STORAGE, LIS_EOPTIONS_STORAGE_BLOCK, LIS_EOPTIONS_PRECISION,
LIS_EOPTIONS_RVAL
};

char *lis_esolver_atoi[] = {"pi", "ii", "rqi", "cg", "cr", "si", "li", "ai", "gpi", "gii", "grqi", "gcg", "gcr", "gsi", "gli", "gai", "feast"};
char *lis_eprint_atoi[] = {"none", "mem", "out", "all"};
char *lis_etruefalse_atoi[] = {"false", "true"};
char *lis_estorage_atoi[] = {"csr", "csc", "msr", "dia", "ell", "jad", "bsr", "bsc", "vbr", "coo", "dns"};
char *lis_eprecision_atoi[] = {"double", "quad", "switch"};

char *lis_esolvername[] = {"", "Power", "Inverse", "Rayleigh Quotient", "CG", "CR", "Subspace", "Lanczos", "Arnoldi", "Generalized Power", "Generalized Inverse", "Generalized Rayleigh Quotient", "Generalized CG", "Generalized CR", "Generalized Subspace", "Generalized Lanczos", "Generalized Arnoldi", "FEAST"};

char *lis_estoragename[]   = {"CSR", "CSC", "MSR", "DIA", "ELL", "JAD", "BSR", "BSC", "VBR", "COO", "DNS"};

char *lis_ereturncode[] = {"LIS_SUCCESS", "LIS_ILL_OPTION", "LIS_BREAKDOWN", "LIS_OUT_OF_MEMORY", "LIS_MAXITER", "LIS_NOT_IMPLEMENTED", "LIS_ERR_FILE_IO"};

char *lis_eprecisionname[] = {"double", "quad", "switch"};

LIS_VECTOR lis_esolver_evalue = NULL;
LIS_INT lis_esolver_evalue_count = 0;

LIS_VECTOR lis_esolver_evector = NULL;
LIS_INT lis_esolver_evector_count = 0;

LIS_VECTOR lis_esolver_residual_history = NULL;
LIS_INT lis_esolver_residual_history_count = 0;

LIS_INT lis_esolver_is_finite_scalar(
        LIS_SCALAR value)
{
        LIS_REAL magnitude;

#ifdef _COMPLEX
        magnitude = cabs(value);
#else
        magnitude = fabs(value);
#endif

        return isfinite((double)magnitude);
}

LIS_INT lis_esolver_is_valid_denominator(
        LIS_SCALAR value)
{
        LIS_REAL magnitude;

#ifdef _COMPLEX
        magnitude = cabs(value);
#else
        magnitude = fabs(value);
#endif

        return isfinite((double)magnitude) &&
               magnitude>0.0;
}

LIS_INT lis_esolver_is_valid_norm(
        LIS_REAL value)
{
        return isfinite((double)value) &&
               value>0.0;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_inner_solve_attempt"
LIS_INT lis_esolver_inner_solve_attempt(
        LIS_MATRIX A,
        LIS_VECTOR b,
        LIS_VECTOR x,
        LIS_SOLVER solver,
        LIS_PRECON precon)
{
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        if( A==NULL || b==NULL || x==NULL ||
            solver==NULL || precon==NULL )
        {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                        "invalid eigensolver inner solve arguments\n");
                return LIS_ERR_ILL_ARG;
        }

        solver->retcode = LIS_SUCCESS;

        err = lis_solve_kernel(
                A,b,x,solver,precon);
        if( err )
        {
                solver->retcode = err;
                return err;
        }

        if( solver->retcode!=LIS_SUCCESS )
        {
                return solver->retcode;
        }

        if( !isfinite((double)solver->resid) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite eigensolver inner solve residual\n");
                solver->retcode = LIS_BREAKDOWN;
                return LIS_BREAKDOWN;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}

LIS_INT lis_esolver_inner_solve_retryable(
        LIS_INT status)
{
        return status==LIS_BREAKDOWN ||
               status==LIS_MAXITER;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_inner_solve_certify"
LIS_INT lis_esolver_inner_solve_certify(
        LIS_MATRIX A,
        LIS_VECTOR b,
        LIS_VECTOR x,
        LIS_SOLVER solver,
        LIS_INT status,
        LIS_REAL tolerance)
{
        LIS_INT err;
        LIS_REAL residual_norm;
        LIS_REAL rhs_norm;
        LIS_REAL relative_residual;

        LIS_DEBUG_FUNC_IN;

        if( status==LIS_SUCCESS )
        {
                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        if( !lis_esolver_inner_solve_retryable(status) )
        {
                LIS_DEBUG_FUNC_OUT;
                return status;
        }

        if( A==NULL || b==NULL || x==NULL ||
            solver==NULL )
        {
                LIS_DEBUG_FUNC_OUT;
                return status;
        }

        if( !isfinite((double)tolerance) ||
            tolerance<=0.0 )
        {
                LIS_DEBUG_FUNC_OUT;
                return status;
        }

        residual_norm = 0.0;
        rhs_norm = 0.0;

        err = lis_solver_compute_true_residual(
                A,b,x,
                &residual_norm,
                &rhs_norm);
        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return status;
        }

        if( !isfinite((double)residual_norm) ||
            !isfinite((double)rhs_norm) ||
            residual_norm<0.0 ||
            rhs_norm<0.0 )
        {
                LIS_DEBUG_FUNC_OUT;
                return status;
        }

        relative_residual =
                rhs_norm>0.0 ?
                residual_norm/rhs_norm :
                residual_norm;

        if( !isfinite((double)relative_residual) ||
            relative_residual>tolerance )
        {
                LIS_DEBUG_FUNC_OUT;
                return status;
        }

        solver->resid = relative_residual;
        solver->retcode = LIS_SUCCESS;

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}

LIS_SCALAR lis_esolver_shift_perturbation(
        LIS_SCALAR shift)
{
        LIS_REAL magnitude;
        LIS_REAL scale;
        LIS_REAL step;

#ifdef _COMPLEX
        magnitude = cabs(shift);
#else
        magnitude = fabs(shift);
#endif

        if( !isfinite((double)magnitude) )
        {
                return 0.0;
        }

        scale = magnitude>1.0 ? magnitude : 1.0;

#if defined(_LONG__DOUBLE)
#ifdef HAVE_QUADMATH_H
        step =
                (LIS_REAL)(sqrtq)(FLT128_EPSILON) *
                scale;
#else
        step =
                (LIS_REAL)(sqrtl)(LDBL_EPSILON) *
                scale;
#endif
#else
        /*
         * Use the real sqrt function explicitly.  In complex builds
         * lis.h maps the sqrt macro to csqrt, while the perturbation
         * scale is always a real quantity.
         */
        step =
                (LIS_REAL)(sqrt)(DBL_EPSILON) *
                scale;
#endif

#ifdef _COMPLEX
        if( magnitude>0.0 )
        {
                return shift/magnitude*step;
        }
#endif

        return (LIS_SCALAR)step;
}

/*
 * Materialize alpha*A + beta*B only when the operator has serial CSR
 * operands with identical sparsity. This preserves the established fast
 * setup path for factorizing preconditioners. Other row-capable operators
 * fall back to the generic direct-row path.
 */
static LIS_INT
lis_esolver_operator_pc_matrix_csr(
        LIS_MATRIX op,
        LIS_MATRIX *P)
{
        LIS_MATRIX A,B,C;
        LIS_INT i,j,err;

        if( op==NULL || P==NULL )
                return LIS_ERR_ILL_ARG;

        *P = NULL;

        if( op->matrix_type!=LIS_MATRIX_OPERATOR ||
            op->operator_scale!=NULL ||
            op->operator_A==NULL ||
            op->operator_B==NULL )
                return LIS_ERR_NOT_IMPLEMENTED;

        A = op->operator_A;
        B = op->operator_B;

        if( A->matrix_type!=LIS_MATRIX_CSR ||
            B->matrix_type!=LIS_MATRIX_CSR ||
            A->nprocs!=1 || B->nprocs!=1 ||
            A->n!=A->gn || B->n!=B->gn ||
            A->n!=B->n || A->gn!=B->gn ||
            A->nnz!=B->nnz ||
            A->ptr==NULL || B->ptr==NULL ||
            A->index==NULL || B->index==NULL ||
            A->value==NULL || B->value==NULL )
                return LIS_ERR_NOT_IMPLEMENTED;

        for(i=0;i<=A->n;i++)
                if( A->ptr[i]!=B->ptr[i] )
                        return LIS_ERR_NOT_IMPLEMENTED;

        for(j=0;j<A->nnz;j++)
                if( A->index[j]!=B->index[j] )
                        return LIS_ERR_NOT_IMPLEMENTED;

        C = NULL;
        err = lis_matrix_duplicate(A,&C);
        if( err ) return err;

        lis_matrix_set_type(C,LIS_MATRIX_CSR);
        err = lis_matrix_convert(A,C);
        if( err )
        {
                lis_matrix_destroy(C);
                return err;
        }

        if( C->value==NULL || C->nnz!=A->nnz )
        {
                lis_matrix_destroy(C);
                return LIS_ERR_ILL_ARG;
        }

        for(j=0;j<C->nnz;j++)
                C->value[j] =
                        op->operator_alpha * A->value[j] +
                        op->operator_beta  * B->value[j];

        *P = C;
        return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_rebuild_precon"
LIS_INT lis_esolver_rebuild_precon(
        LIS_MATRIX A,
        LIS_SOLVER solver,
        LIS_PRECON *precon)
{
        LIS_INT err;
        LIS_INT precon_type;
        LIS_MATRIX pc_matrix;
        LIS_MATRIX build_matrix;

        LIS_DEBUG_FUNC_IN;

        if( A==NULL || solver==NULL || precon==NULL )
        {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                        "invalid eigensolver preconditioner rebuild arguments\n");
                return LIS_ERR_ILL_ARG;
        }

        if( *precon )
        {
                if( solver->precon==*precon ) solver->precon = NULL;
                lis_precon_destroy(*precon);
                *precon = NULL;
        }

        lis_solver_work_destroy(solver);
        solver->A = A;
        solver->retcode = LIS_SUCCESS;

        precon_type = solver->options[LIS_OPTIONS_PRECON];
        pc_matrix = NULL;
        build_matrix = A;

        if( A->matrix_type==LIS_MATRIX_OPERATOR &&
            (precon_type==LIS_PRECON_TYPE_ILU ||
             precon_type==LIS_PRECON_TYPE_ILUT ||
             precon_type==LIS_PRECON_TYPE_ILUC) )
        {
                err = lis_esolver_operator_pc_matrix_csr(A,&pc_matrix);
                if( err==LIS_SUCCESS )
                {
                        build_matrix = pc_matrix;
                }
                else if( err!=LIS_ERR_NOT_IMPLEMENTED )
                {
                        solver->retcode = err;
                        return err;
                }
                /* NOT_IMPLEMENTED means use the generic row-backed path. */
        }

        solver->A = build_matrix;
        err = lis_precon_create(solver,precon);
        solver->A = A;

        if( pc_matrix )
        {
                lis_matrix_destroy(pc_matrix);
                pc_matrix = NULL;
        }

        if( err )
        {
                *precon = NULL;
                solver->precon = NULL;
                solver->retcode = err;
                return err;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_standard_residual"
LIS_INT lis_esolver_standard_residual(
        LIS_MATRIX A,
        LIS_VECTOR x,
        LIS_SCALAR lambda,
        LIS_VECTOR Ax,
        LIS_VECTOR r,
        LIS_REAL *resid)
{
        LIS_INT err;
        LIS_REAL anrm,xnrm,rnrm,denom;
        LIS_REAL lmag;

        LIS_DEBUG_FUNC_IN;

        err = lis_matvec(A,x,Ax);
        if( err ) return err;

        err = lis_vector_nrm2(Ax,&anrm);
        if( err ) return err;

        err = lis_vector_nrm2(x,&xnrm);
        if( err ) return err;

        lis_vector_axpyz(-lambda,x,Ax,r);

        err = lis_vector_nrm2(r,&rnrm);
        if( err ) return err;

#ifdef _COMPLEX
        lmag = cabs(lambda);
#else
        lmag = fabs(lambda);
#endif

        if( !isfinite((double)anrm) ||
            !isfinite((double)xnrm) ||
            !isfinite((double)rnrm) ||
            !isfinite((double)lmag) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite standard eigenpair certification data\n");
                return LIS_BREAKDOWN;
        }

        denom = anrm+lmag*xnrm;
        if( !isfinite((double)denom) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite standard eigenpair certification denominator\n");
                return LIS_BREAKDOWN;
        }

        if( denom>0.0 )
        {
                *resid = rnrm/denom;
        }
        else
        {
                *resid = rnrm;
        }

        if( !isfinite((double)*resid) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite standard eigenpair certification residual\n");
                return LIS_BREAKDOWN;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_generalized_residual"
LIS_INT lis_esolver_generalized_residual(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_VECTOR x,
        LIS_SCALAR lambda,
        LIS_VECTOR Ax,
        LIS_VECTOR Bx,
        LIS_VECTOR r,
        LIS_REAL *resid)
{
        LIS_INT err;
        LIS_REAL anrm,bnrm,rnrm,denom;
        LIS_REAL lmag;

        LIS_DEBUG_FUNC_IN;

        err = lis_matvec(A,x,Ax);
        if( err ) return err;

        err = lis_matvec(B,x,Bx);
        if( err ) return err;

        err = lis_vector_nrm2(Ax,&anrm);
        if( err ) return err;

        err = lis_vector_nrm2(Bx,&bnrm);
        if( err ) return err;

        lis_vector_axpyz(-lambda,Bx,Ax,r);

        err = lis_vector_nrm2(r,&rnrm);
        if( err ) return err;

#ifdef _COMPLEX
        lmag = cabs(lambda);
#else
        lmag = fabs(lambda);
#endif

        if( !isfinite((double)anrm) ||
            !isfinite((double)bnrm) ||
            !isfinite((double)rnrm) ||
            !isfinite((double)lmag) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite generalized eigenpair certification data\n");
                return LIS_BREAKDOWN;
        }

        denom = anrm+lmag*bnrm;
        if( !isfinite((double)denom) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite generalized eigenpair certification denominator\n");
                return LIS_BREAKDOWN;
        }

        if( denom>0.0 )
        {
                *resid = rnrm/denom;
        }
        else
        {
                *resid = rnrm;
        }

        if( !isfinite((double)*resid) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite generalized eigenpair certification residual\n");
                return LIS_BREAKDOWN;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_init"
LIS_INT lis_esolver_init(LIS_ESOLVER esolver)
{

	LIS_DEBUG_FUNC_IN;

	esolver->A          = NULL;
	esolver->x          = NULL;
	esolver->xx         = NULL;
	esolver->evalue     = NULL;
	esolver->evector    = NULL;
	esolver->evector_owned_count = 0;
	esolver->resid      = NULL;
	esolver->work       = NULL;
	esolver->rhistory   = NULL;
	esolver->iter       = NULL;
	esolver->iter2      = NULL;

	esolver->worklen    = 0;
	esolver->time       = 0;
	esolver->itime      = 0;
	esolver->ptime      = 0;
	esolver->p_i_time   = 0;
	esolver->p_c_time   = 0;

	esolver->ishift     = 0;
	esolver->transform  = LIS_EIGEN_TRANSFORM_AUTO;
	esolver->nev        = 0;
	esolver->requested_nev = 0;
	esolver->requested_nev_active = LIS_FALSE;
	esolver->block_size = 1;
	esolver->linear_options = NULL;
	esolver->refine_linear_options = NULL;
	esolver->refine_gap_fraction = 0.01;
	esolver->subspace_retry_max = 0;
	esolver->subspace_retry_growth = 2.0;
	esolver->subspace_retry_min = 0;
	esolver->subspace_retry_attempts = 0;
	esolver->subspace_retry_final_size = 0;
	esolver->tol        = 0;
	esolver->eprecision = LIS_PRECISION_DOUBLE;

        esolver->feast_emin          = 0.0;
        esolver->feast_emax          = 0.0;
        esolver->feast_ratio         = 0.30;
        esolver->feast_points        = 8;
        esolver->feast_solver        = LIS_SOLVER_GMRES;
        esolver->feast_interval_mask = 0;

	esolver->options[LIS_EOPTIONS_ESOLVER]                = LIS_ESOLVER_CR;
	esolver->options[LIS_EOPTIONS_MAXITER]                = 1000;
	esolver->options[LIS_EOPTIONS_SUBSPACE]               = 1;
	esolver->options[LIS_EOPTIONS_MODE]                   = 0;
	esolver->options[LIS_EOPTIONS_OUTPUT]                 = LIS_FALSE;
	esolver->options[LIS_EOPTIONS_INITGUESS_ONES]         = LIS_TRUE;
	esolver->options[LIS_EOPTIONS_INNER_ESOLVER]          = LIS_ESOLVER_II;
	esolver->options[LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER] = LIS_ESOLVER_GII;	
	esolver->options[LIS_EOPTIONS_STORAGE]                = 0;
	esolver->options[LIS_EOPTIONS_STORAGE_BLOCK]          = 2;
	esolver->options[LIS_EOPTIONS_PRECISION]              = LIS_PRECISION_DOUBLE;
	esolver->options[LIS_EOPTIONS_RVAL]                   = LIS_FALSE;
	esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN] = 1.0e-12;
	esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN] = 0.0;
	esolver->params[LIS_EPARAMS_SHIFT_IM - LIS_EOPTIONS_LEN] = 0.0;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


static LIS_INT lis_esolver_result_count(
        LIS_ESOLVER esolver)
{
        LIS_INT ss;

        ss = esolver->options[
                LIS_EOPTIONS_SUBSPACE];

        /*
         * Block Lanczos may terminate with fewer active basis
         * vectors than the nominal subspace dimension after
         * numerical rank deflation or happy breakdown.
         * Expose only the active reconstructed Ritz results.
         */
        if( esolver->options[
                    LIS_EOPTIONS_ESOLVER]
                ==LIS_ESOLVER_LI &&
            esolver->block_size>1 &&
            esolver->nev>0 &&
            esolver->nev<=ss )
        {
                return esolver->nev;
        }

        /*
         * Native FEAST returns only eigenpairs certified inside its
         * requested spectral interval.  nev may legitimately be
         * zero for an empty converged interval.
         */
        if( esolver->options[
                    LIS_EOPTIONS_ESOLVER]
                ==LIS_ESOLVER_FEAST )
        {
                if( esolver->nev>=0 &&
                    esolver->nev<=ss )
                {
                        return esolver->nev;
                }

                return ss;
        }

        if( esolver->options[
                    LIS_EOPTIONS_ESOLVER]
                ==LIS_ESOLVER_GLI &&
            esolver->transform
                !=LIS_EIGEN_TRANSFORM_AUTO &&
            esolver->nev>0 &&
            esolver->nev<ss )
        {
                return esolver->nev;
        }

        return ss;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_create"
LIS_INT lis_esolver_create(LIS_ESOLVER *esolver)
{
	LIS_DEBUG_FUNC_IN;

	*esolver = NULL;

	*esolver = (LIS_ESOLVER)lis_malloc( sizeof(struct LIS_ESOLVER_STRUCT),"lis_esolver_create::esolver" );
	if( NULL==*esolver )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_ESOLVER_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_esolver_init(*esolver);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_work_destroy"
LIS_INT lis_esolver_work_destroy(LIS_ESOLVER esolver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( esolver && esolver->work )
	{
		for(i=0;i<esolver->worklen;i++) lis_vector_destroy(esolver->work[i]);
		lis_free(esolver->work);
		esolver->work    = NULL;
		esolver->worklen = 0;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_destroy"
LIS_INT lis_esolver_destroy(LIS_ESOLVER esolver)
{
	LIS_INT i;
	LIS_DEBUG_FUNC_IN;

	if( esolver )
	{
		lis_esolver_work_destroy(esolver);
   	        if( esolver->rhistory ) lis_free(esolver->rhistory);
	        if( esolver->evalue ) lis_free(esolver->evalue);
 	        if( esolver->resid ) lis_free(esolver->resid);
 	        if( esolver->iter) lis_free(esolver->iter);
 	        if( esolver->iter2) lis_free(esolver->iter2);
                if( esolver->evector )
                {
                        /*
                         * Ownership follows the saved result-array metadata,
                         * not the eigensolver options currently selected.
                         */
                        for(i=0;i<esolver->evector_owned_count;i++)
                        {
                                if( esolver->evector[i] )
                                {
                                        lis_vector_destroy(
                                                esolver->evector[i]);
                                }
                        }

                        lis_free(esolver->evector);
                        esolver->evector = NULL;
                        esolver->evector_owned_count = 0;
                }

                if( esolver->linear_options )
                {
                        lis_free(esolver->linear_options);
                        esolver->linear_options = NULL;
                }

                if( esolver->refine_linear_options )
                {
                        lis_free(esolver->refine_linear_options);
                        esolver->refine_linear_options = NULL;
                }

		lis_free(esolver);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

static LIS_INT
lis_esolver_replace_option_string(
        char **slot,
        const char *text)
{
        char *copy;
        size_t length;

        if( slot==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "invalid eigensolver option-string target\n");
                return LIS_ERR_ILL_ARG;
        }

        if( *slot )
        {
                lis_free(*slot);
                *slot = NULL;
        }

        if( text==NULL || text[0]=='\0' )
        {
                return LIS_SUCCESS;
        }

        length = strlen(text);
        copy = (char *)lis_malloc(
                (length+1)*sizeof(char),
                "lis_esolver_replace_option_string::copy");

        if( copy==NULL )
        {
                LIS_SETERR_MEM((length+1)*sizeof(char));
                return LIS_OUT_OF_MEMORY;
        }

        memcpy(copy,text,length+1);
        *slot = copy;
        return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_esolver_set_linear_options"
LIS_INT lis_esolver_set_linear_options(
        const char *text,
        LIS_ESOLVER esolver)
{
        LIS_DEBUG_FUNC_IN;

        if( esolver==NULL )
        {
                LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver is undefined\n");
                return LIS_ERR_ILL_ARG;
        }

        LIS_DEBUG_FUNC_OUT;
        return lis_esolver_replace_option_string(
                &esolver->linear_options,text);
}


#undef __FUNC__
#define __FUNC__ "lis_esolver_set_refine_linear_options"
LIS_INT lis_esolver_set_refine_linear_options(
        const char *text,
        LIS_ESOLVER esolver)
{
        LIS_DEBUG_FUNC_IN;

        if( esolver==NULL )
        {
                LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver is undefined\n");
                return LIS_ERR_ILL_ARG;
        }

        LIS_DEBUG_FUNC_OUT;
        return lis_esolver_replace_option_string(
                &esolver->refine_linear_options,text);
}


#undef __FUNC__
#define __FUNC__ "lis_esolver_set_refine_gap_fraction"
LIS_INT lis_esolver_set_refine_gap_fraction(
        LIS_REAL fraction,
        LIS_ESOLVER esolver)
{
        LIS_DEBUG_FUNC_IN;

        if( esolver==NULL )
        {
                LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver is undefined\n");
                return LIS_ERR_ILL_ARG;
        }

        if( !isfinite((double)fraction) ||
            fraction<=0.0 || fraction>=0.5 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "eigensolver refine gap fraction must be in (0,0.5)\n");
                return LIS_ERR_ILL_ARG;
        }

        esolver->refine_gap_fraction = fraction;
        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}



#undef __FUNC__
#define __FUNC__ "lis_esolver_set_subspace_retry"
LIS_INT lis_esolver_set_subspace_retry(
        LIS_ESOLVER esolver,
        LIS_INT max_retries,
        LIS_REAL growth_factor,
        LIS_INT min_subspace)
{
        LIS_DEBUG_FUNC_IN;

        if( esolver==NULL || max_retries<0 || min_subspace<0 ||
            !isfinite((double)growth_factor) || growth_factor<=1.0 )
        {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                        "invalid eigensolver subspace retry policy\n");
                return LIS_ERR_ILL_ARG;
        }

        esolver->subspace_retry_max = max_retries;
        esolver->subspace_retry_growth = growth_factor;
        esolver->subspace_retry_min = min_subspace;
        esolver->subspace_retry_attempts = 0;
        esolver->subspace_retry_final_size =
                esolver->options[LIS_EOPTIONS_SUBSPACE];

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_subspace_retry"
LIS_INT lis_esolver_get_subspace_retry(
        LIS_ESOLVER esolver,
        LIS_INT *max_retries,
        LIS_REAL *growth_factor,
        LIS_INT *min_subspace,
        LIS_INT *attempts,
        LIS_INT *final_subspace)
{
        LIS_DEBUG_FUNC_IN;

        if( esolver==NULL )
        {
                LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver is undefined\n");
                return LIS_ERR_ILL_ARG;
        }

        if( max_retries ) *max_retries = esolver->subspace_retry_max;
        if( growth_factor ) *growth_factor = esolver->subspace_retry_growth;
        if( min_subspace ) *min_subspace = esolver->subspace_retry_min;
        if( attempts ) *attempts = esolver->subspace_retry_attempts;
        if( final_subspace ) *final_subspace = esolver->subspace_retry_final_size;

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}

static LIS_INT
lis_esolver_subspace_retry_supported(LIS_INT nesolver)
{
        return nesolver==LIS_ESOLVER_SI ||
               nesolver==LIS_ESOLVER_LI ||
               nesolver==LIS_ESOLVER_AI ||
               nesolver==LIS_ESOLVER_GSI ||
               nesolver==LIS_ESOLVER_GLI ||
               nesolver==LIS_ESOLVER_GAI;
}

#undef __FUNC__
#define __FUNC__ "lis_esolve"
LIS_INT lis_esolve(LIS_MATRIX A, LIS_VECTOR x, LIS_SCALAR *evalue0, LIS_ESOLVER esolver)
{
	LIS_SCALAR *evalue;
	LIS_INT	err;

	LIS_DEBUG_FUNC_IN;

	/* call generalized eigensolver */
	err = lis_gesolve(A,NULL,x,evalue0,esolver);

	if( err )
	{
		esolver->retcode = err;
		return err;
	}
	
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_gesolve_once"
static LIS_INT lis_gesolve_once(LIS_MATRIX A, LIS_MATRIX B, LIS_VECTOR x, LIS_SCALAR *evalue0, LIS_ESOLVER esolver)
{
	LIS_Comm comm;   
	LIS_INT	nesolver,niesolver,nigesolver,emaxiter; 
	LIS_SCALAR *evalue;
	LIS_VECTOR *evector;
	LIS_REAL *resid;
	LIS_REAL *rhistory;
	LIS_INT	*iter,*iter2;
	LIS_INT	err;
	LIS_INT output;
	LIS_INT i;
	LIS_INT ss, mode;
	double time;
	LIS_SCALAR oshift,ishift;
	LIS_INT	estorage,eblock;
	LIS_MATRIX A0;
	LIS_INT eprecision;
	LIS_INT rval;
	LIS_VECTOR xx;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	/* begin parameter check */
	err = lis_matrix_check(A,LIS_MATRIX_CHECK_ALL);
	if( B!=NULL ) err = lis_matrix_check(B,LIS_MATRIX_CHECK_ALL);	

	if( err ) return err;
	if( x==NULL )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"vector x is undefined\n");
		return LIS_ERR_ILL_ARG;
	}
	if( A->n!=x->n )
	{
		return LIS_ERR_ILL_ARG;
	}
	if( B!=NULL && B->n!=x->n )
	{
		return LIS_ERR_ILL_ARG;
	}
	if( A->gn<=0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Size n(=%D) of matrix A is less than 0\n",A->gn);
		return LIS_ERR_ILL_ARG;
	}
	if( B!=NULL && B->gn<=0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Size n(=%D) of matrix B is less than 0\n",B->gn);
		return LIS_ERR_ILL_ARG;
	}

	nesolver = esolver->options[LIS_EOPTIONS_ESOLVER];
	niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];
	nigesolver = esolver->options[LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER];	
	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
	mode = esolver->options[LIS_EOPTIONS_MODE];
	emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
#ifdef _COMPLEX
	oshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN] + esolver->params[LIS_EPARAMS_SHIFT_IM - LIS_EOPTIONS_LEN] * _Complex_I;
#else
	oshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN];
#endif	
	ishift = esolver->ishift;
	output = esolver->options[LIS_EOPTIONS_OUTPUT];
	estorage = esolver->options[LIS_EOPTIONS_STORAGE];
	eblock = esolver->options[LIS_EOPTIONS_STORAGE_BLOCK];
	eprecision = esolver->options[LIS_EOPTIONS_PRECISION];
	esolver->eprecision = eprecision;
	rval = esolver->options[LIS_EOPTIONS_RVAL];

	if( esolver->transform < LIS_EIGEN_TRANSFORM_AUTO ||
	    esolver->transform > LIS_EIGEN_TRANSFORM_BUCKLING )
	{
		LIS_SETERR1(
			LIS_ERR_ILL_ARG,
			"Unknown spectral transformation %D\n",
			esolver->transform);
		return LIS_ERR_ILL_ARG;
	}

	if( esolver->transform != LIS_EIGEN_TRANSFORM_AUTO )
	{
		if( B==NULL )
		{
			LIS_SETERR(
				LIS_ERR_ILL_ARG,
				"spectral transformations require "
" a generalized eigenproblem\n");
			return LIS_ERR_ILL_ARG;
		}

		/*
		 * GII is the reference implementation for explicit
		 * spectral transformations. Transformed Krylov
		 * recurrences are enabled separately.
		 */
		if( nesolver != LIS_ESOLVER_GII &&
		    nesolver != LIS_ESOLVER_GLI )
		{
			LIS_SETERR(
				LIS_ERR_NOT_IMPLEMENTED,
				"spectral transformation is not implemented "
"for this eigensolver\n");
			return LIS_ERR_NOT_IMPLEMENTED;
		}

#ifdef _COMPLEX
            /*
             * The transformed generalized Lanczos path currently uses
             * a real symmetric projected eigensystem. Reject complex
             * arithmetic before allocating Lanczos work or starting
             * shifted solves.
             */
            if( nesolver == LIS_ESOLVER_GLI )
            {
                    LIS_SETERR(
                            LIS_ERR_NOT_IMPLEMENTED,
                            "transformed generalized Lanczos is not "
"implemented for complex arithmetic\n");
                    return LIS_ERR_NOT_IMPLEMENTED;
            }
#endif

		if( esolver->transform ==
		        LIS_EIGEN_TRANSFORM_BUCKLING &&
		    oshift==0.0 && ishift==0.0 )
		{
			LIS_SETERR(
				LIS_ERR_ILL_ARG,
				"buckling transformation requires "
" a nonzero shift\n");
			return LIS_ERR_ILL_ARG;
		}

	}

	/*
	 * LIS_MATRIX_USER is an operator, not materialized LIS storage.
	 * Explicit storage conversion and matrix shifts are separate features.
	 */
	if( A->matrix_type==LIS_MATRIX_USER && estorage>0 )
	{
		LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
		           "storage conversion is unavailable for LIS_MATRIX_USER\n");
		return LIS_ERR_NOT_IMPLEMENTED;
	}
	if( (A->matrix_type==LIS_MATRIX_USER ||
	     (B!=NULL && B->matrix_type==LIS_MATRIX_USER)) &&
	    (oshift!=0.0 || ishift!=0.0) &&
	    !(B!=NULL &&
	  (nesolver==LIS_ESOLVER_GPI ||
	   nesolver==LIS_ESOLVER_GII ||
	   nesolver==LIS_ESOLVER_GCR ||
	   (nesolver==LIS_ESOLVER_GLI &&
	    esolver->transform!=LIS_EIGEN_TRANSFORM_AUTO))) )
	{
		LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
		           "matrix shifts are unavailable for LIS_MATRIX_USER\n");
		return LIS_ERR_NOT_IMPLEMENTED;
	}

	if( nesolver < 1 || nesolver > LIS_ESOLVER_LEN )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %D (Set between 1 to %D)\n",nesolver, LIS_ESOLVER_LEN);
		return LIS_ERR_ILL_ARG;
	}
	/*
	 * Block eigensolver mode is introduced incrementally.
	 *
	 * Standard Lanczos has a complete block path.  Generalized
	 * Lanczos owns the next block-dispatch boundary and rejects
	 * incomplete numerical stages inside lis_egli_block().
	 * Other eigensolvers retain their scalar behavior.
	 */
	if( esolver->block_size>1 &&
	    nesolver!=LIS_ESOLVER_LI &&
	    nesolver!=LIS_ESOLVER_GLI )
	{
		LIS_SETERR(
			LIS_ERR_NOT_IMPLEMENTED,
			"block eigensolver mode is currently available "
			"only for standard or generalized Lanczos\n");
		return LIS_ERR_NOT_IMPLEMENTED;
	}

	if( esolver->block_size>1 &&
	    esolver->block_size>ss )
	{
		LIS_SETERR2(
			LIS_ERR_ILL_ARG,
			"Parameter -eblock is %D "
			"(Set less than or equal to subspace size %D)\n",
			esolver->block_size,
			ss);
		return LIS_ERR_ILL_ARG;
	}

	if( B!=NULL &&
	    (A->matrix_type==LIS_MATRIX_USER || B->matrix_type==LIS_MATRIX_USER) &&
	    nesolver==LIS_ESOLVER_GRQI )
	{
		LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
		           "generalized Rayleigh quotient iteration requires matrix shifts for LIS_MATRIX_USER\n");
		return LIS_ERR_NOT_IMPLEMENTED;
	}

	/*
	if( niesolver < 1 || niesolver > 6 ) 
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_INNER_ESOLVER is %D (Set between 1 to 6)\n", niesolver);
		return LIS_ERR_ILL_ARG;
	}

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_SI && niesolver > 6 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_INNER_ESOLVER is %D (Set Power or Inverse for Subspace)\n", niesolver);
		return LIS_ERR_ILL_ARG;
	}

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_LI && niesolver == LIS_ESOLVER_PI )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_INNER_ESOLVER is %D (Power cannot used for Lanczos)\n", niesolver);
		return LIS_ERR_ILL_ARG;
	}

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_AI && niesolver == LIS_ESOLVER_PI )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_INNER_ESOLVER is %D (Power cannot used for Arnoldi)\n", niesolver);
		return LIS_ERR_ILL_ARG;
	}
	*/

	if (( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_SI || esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GSI ) && ss > A->gn )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_SUBSPACE is %D (Set less than or equal to matrix size %D for Subspace)\n", ss, A->gn);
		return LIS_ERR_ILL_ARG;
	}

	if (( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_LI || esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_AI || esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GLI || esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GAI) && ss > A->gn )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_SUBSPACE is %D (Set less than or equal to matrix size %D for Lanczos and Arnoldi)\n", ss, A->gn);
		return LIS_ERR_ILL_ARG;
	}

	if (( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_SI || esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GSI ) && mode >= ss )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_MODE is %D (Set less than subspace size %D for Subspace)\n", mode, ss);
		return LIS_ERR_ILL_ARG;
	}

	if (( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_LI || esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_AI || esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GLI || esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GAI) && mode >= ss )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_MODE is %D (Set less than subspace size %D for Lanczos and Arnoldi)\n", mode, ss);
		return LIS_ERR_ILL_ARG;
	}

	#ifdef USE_QUAD_PRECISION
		if( eprecision==LIS_PRECISION_QUAD && lis_esolver_execute_quad[nesolver]==NULL )
		{
			LIS_SETERR1(LIS_ERR_NOT_IMPLEMENTED,"Quad precision eigensolver %s is not implemented\n",lis_esolvername[nesolver]);
			return LIS_ERR_NOT_IMPLEMENTED;
		}
		else if( eprecision==LIS_PRECISION_SWITCH && lis_esolver_execute_switch[nesolver]==NULL )
		{
			LIS_SETERR1(LIS_ERR_NOT_IMPLEMENTED,"Switch esolver %s is not implemented\n",lis_esolvername[nesolver]);
			return LIS_ERR_NOT_IMPLEMENTED;
		}
		if( esolver->options[LIS_EOPTIONS_SWITCH_MAXITER]==-1 )
		{
			esolver->options[LIS_EOPTIONS_SWITCH_MAXITER] = emaxiter;
		}
	#endif

	err = lis_esolver_check_params[nesolver](esolver);
	if( err )
	{
		esolver->retcode = err;
		return err;
	}
	/* end parameter check */

	/* create eigenvalue array */
	if( esolver->evalue ) lis_free(esolver->evalue);
	evalue = (LIS_SCALAR *)lis_malloc((ss+2)*sizeof(LIS_SCALAR),"lis_esolve::evalue");
	if( evalue==NULL )
	{
		LIS_SETERR_MEM((ss+2)*sizeof(LIS_SCALAR));
		esolver->retcode = err;
		return err;
	}
	evalue[0] = 1.0;
	evalue[ss-1] = 1.0;

	/* create residual norm array */
	if( esolver->resid ) lis_free(esolver->resid);
	resid = (LIS_REAL *)lis_malloc((ss+2)*sizeof(LIS_REAL),"lis_esolve::resid");
	if( resid==NULL )
	{
		LIS_SETERR_MEM((ss+2)*sizeof(LIS_SCALAR));
		esolver->retcode = err;
		return err;
	}

	/* create number of iterations array */
	if( esolver->iter ) lis_free(esolver->iter);
	iter = (LIS_INT *)lis_malloc((ss+2)*sizeof(LIS_SCALAR),"lis_esolve::iter");
	if( iter==NULL )
	{
		LIS_SETERR_MEM((ss+2)*sizeof(LIS_SCALAR));
		esolver->retcode = err;
		return err;
	}

	/* create quad precision number of iterations array */
	if( esolver->iter2 ) lis_free(esolver->iter2);
	iter2 = (LIS_INT *)lis_malloc((ss+2)*sizeof(LIS_SCALAR),"lis_esolve::iter2");
	if( iter2==NULL )
	{
		LIS_SETERR_MEM((ss+2)*sizeof(LIS_SCALAR));
		esolver->retcode = err;
		return err;
	}

	/* create initial vector */
	#ifndef USE_QUAD_PRECISION
		err = lis_vector_duplicate(A,&xx);
	#else
		if( eprecision==LIS_PRECISION_DOUBLE )
		{
			err = lis_vector_duplicate(A,&xx);
		}
		else
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,A,&xx);
		}
	#endif
	if( err )
	{
		esolver->retcode = err;
		return err;
	}
	if( esolver->options[LIS_EOPTIONS_INITGUESS_ONES] )
	{
	  if( output ) lis_printf(comm,"initial vector x      : all components set to 1\n");
		#ifndef USE_QUAD_PRECISION
			lis_vector_set_all(1.0,xx);
		#else
			if( eprecision==LIS_PRECISION_DOUBLE )
			{
				lis_vector_set_all(1.0,xx);
			}
			else
			{
				lis_vector_set_allex_nm(1.0,xx);
			}
		#endif
	}
	else
	{
	  if( output ) lis_printf(comm,"initial vector x      : user defined\n"); 
		#ifndef USE_QUAD_PRECISION
			lis_vector_copy(x,xx);
		#else
			if( eprecision==LIS_PRECISION_DOUBLE )
			{
				lis_vector_copy(x,xx);
			}
			else
			{
				lis_vector_copyex_nm(x,xx);
			}
		#endif
	}

	/* create eigenvector array */
	if( esolver->evector )
	{
	        /*
	         * Release vectors owned by the previous result array
	         * before replacing its pointer storage.
	         */
	        for(i=0;i<esolver->evector_owned_count;i++)
	        {
	                if( esolver->evector[i] )
	                {
	                        lis_vector_destroy(
	                                esolver->evector[i]);
	                }
	        }

	        lis_free(esolver->evector);
	        esolver->evector = NULL;
	        esolver->evector_owned_count = 0;
	}
	evector = (LIS_VECTOR *)lis_malloc((ss+2)*sizeof(LIS_VECTOR),"lis_esolve::evector");
	if( evector==NULL )
	{
		LIS_SETERR_MEM((ss+2)*sizeof(LIS_VECTOR));
		esolver->retcode = err;
		return err;
	}

	    memset(evector,0,(ss+2)*sizeof(LIS_VECTOR));

/* create residual history array */
	if( esolver->rhistory ) lis_free(esolver->rhistory);
	rhistory = (LIS_REAL *)lis_malloc((emaxiter+2)*sizeof(LIS_REAL),"lis_esolve::rhistory");
	if( rhistory==NULL )
	{
		LIS_SETERR_MEM((emaxiter+2)*sizeof(LIS_SCALAR));
		lis_vector_destroy(xx);
		esolver->retcode = err;
		return err;
	}

	/* convert matrix */
	if( estorage>0 && A->matrix_type!=estorage )
	{
		err = lis_matrix_duplicate(A,&A0);
		if( err ) return err;
		lis_matrix_set_blocksize(A0,eblock,eblock,NULL,NULL);
		lis_matrix_set_type(A0,estorage);
		err = lis_matrix_convert(A,A0);
		if( err ) return err;
		lis_matrix_storage_destroy(A);
		lis_matrix_DLU_destroy(A);
		lis_matrix_diag_destroy(A->WD);
		if( A->l2g_map ) lis_free( A->l2g_map );
		if( A->commtable ) lis_commtable_destroy( A->commtable );
		if( A->ranges ) lis_free( A->ranges );
		err = lis_matrix_copy_struct(A0,A);
		if( err ) return err;
		lis_free(A0);
	}

	esolver->A        = A;
	esolver->B        = B;
	esolver->evalue   = evalue;
	esolver->x        = x;
	esolver->evector  = evector;

	/*
	 * These eigensolvers own every non-NULL vector they place
	 * in the zero-initialized evector[0..ss+1] array.
	 */
	if( nesolver==LIS_ESOLVER_LI ||
	    nesolver==LIS_ESOLVER_AI ||
	    nesolver==LIS_ESOLVER_SI ||
	    nesolver==LIS_ESOLVER_GLI ||
	    nesolver==LIS_ESOLVER_GSI ||
	    nesolver==LIS_ESOLVER_GAI ||
	    nesolver==LIS_ESOLVER_FEAST )
	{
	        esolver->evector_owned_count = ss+2;
	}
	else
	{
	        esolver->evector_owned_count = 0;
	}
	rhistory[0]       = 1.0;
	esolver->rhistory = rhistory;
	esolver->resid    = resid;
	esolver->iter     = iter;
	esolver->iter2    = iter2;

#ifdef _LONG__DOUBLE
	if ( output ) lis_printf(comm,"precision             : long double\n");
#else
	if ( output ) lis_printf(comm,"precision             : %s\n", lis_eprecisionname[eprecision]);
#endif
	if ( output ) lis_printf(comm,"eigensolver           : %s\n", lis_esolvername[nesolver]);

	if ( output ) lis_printf(comm,"convergence condition : ||lx-(B^-1)Ax||_2 <= %6.1e * ||lx||_2\n", (double)esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]); 

	if( A->matrix_type==LIS_MATRIX_BSR || A->matrix_type==LIS_MATRIX_BSC )
	  {
	    if ( output ) lis_printf(comm,"matrix storage format : %s(%D x %D)\n", lis_estoragename[A->matrix_type-1],eblock,eblock); 
	  }
	else
	  {
	    if( A->matrix_type==LIS_MATRIX_USER )
	      {
	        if ( output ) lis_printf(comm,"matrix storage format : user/shell\n");
	      }
	    else
	      {
	        if ( output ) lis_printf(comm,"matrix storage format : %s\n", lis_estoragename[A->matrix_type-1]);
	      }
	  }
	
	time = lis_wtime();

	esolver->ptime = 0;
	esolver->itime = 0;
	esolver->p_c_time = 0;
	esolver->p_i_time = 0;

	/*
         * A subsequent solve on the same eigensolver replaces its scratch
         * vectors. Release the previous work array before allocating a new
         * one; otherwise repeated solves leak the old work vectors.
         */
        if(esolver->work!=NULL)
        {
                err=lis_esolver_work_destroy(esolver);
                if(err)
                {
                        esolver->x=NULL;
                        esolver->xx=NULL;
                        lis_vector_destroy(xx);
                        esolver->retcode=err;
                        return err;
                }
        }

	/* create work vector */
	err = lis_esolver_malloc_work[nesolver](esolver);
	if( err )
	{
	  esolver->x        = NULL;
	  esolver->xx       = NULL;
	  lis_vector_destroy(xx);
	  esolver->retcode = err;
	  return err;
	}

	esolver->x        = xx;
	esolver->xx       = x;

	/* execute esolver */
	#ifndef USE_QUAD_PRECISION
		err = lis_esolver_execute[nesolver](esolver);
	#else
		if( eprecision==LIS_PRECISION_DOUBLE )
		{
			err = lis_esolver_execute[nesolver](esolver);
		}
		else if( eprecision==LIS_PRECISION_QUAD )
		{
			err = lis_esolver_execute_quad[nesolver](esolver);
		}
		else if( eprecision==LIS_PRECISION_SWITCH )
		{
			err = lis_esolver_execute_switch[nesolver](esolver);
		}
	#endif
	esolver->retcode = err;

	*evalue0 = esolver->evalue[0];
	lis_vector_copy(esolver->x, x);

	esolver->time = lis_wtime() - time; 

	if( err )
	  {
	    if ( output ) lis_printf(comm,"eigensolver status    : %s(code=%D)\n\n",lis_ereturncode[err],err);
	  }
	else
	  {
	    if ( output ) lis_printf(comm,"eigensolver status    : normal end\n\n");
	  }

	if( eprecision==LIS_PRECISION_DOUBLE )
	  {
	    esolver->iter2[mode] = esolver->iter[mode];
	  }
	else if( eprecision==LIS_PRECISION_QUAD )
	  {
	    esolver->iter2[mode] = 0;
	  }

	/* generic solve scratch aliases end here */
	esolver->x        = NULL;
	esolver->xx       = NULL;
	lis_vector_destroy(xx);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_gesolve"
LIS_INT lis_gesolve(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_VECTOR x,
        LIS_SCALAR *evalue0,
        LIS_ESOLVER esolver)
{
        LIS_INT err,current_ss,next_ss,nesolver,attempt;
        LIS_INT original_ss,result;
        LIS_REAL grown;
        double total_time;
#ifdef USE_MPI
        LIS_INT local_status,status_min,status_max;
#endif

        LIS_DEBUG_FUNC_IN;

        if( esolver==NULL )
        {
                LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver is undefined\n");
                return LIS_ERR_ILL_ARG;
        }

        original_ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
        result = LIS_SUCCESS;
        total_time = 0.0;
        esolver->subspace_retry_attempts = 0;
        esolver->subspace_retry_final_size = original_ss;

        for(attempt=0;;attempt++)
        {
                /*
                 * Some solvers write their result count to nev. Restore the
                 * stable typed/text request before EVERY attempt, including
                 * retries and later lis_gesolve() invocations.
                 * Legacy direct nev writes remain valid when no explicit
                 * request was set through the public option APIs.
                 */
                if(esolver->requested_nev_active)
                {
                        if(esolver->requested_nev<0 ||
                           esolver->requested_nev>
                               esolver->options[LIS_EOPTIONS_SUBSPACE])
                        {
                                LIS_SETERR(LIS_ERR_ILL_ARG,
                                    "requested eigenpair count exceeds subspace\n");
                                result=LIS_ERR_ILL_ARG;
                                break;
                        }
                        esolver->nev=esolver->requested_nev;
                }
                err = lis_gesolve_once(A,B,x,evalue0,esolver);
                if( err )
                {
                        result = err;
                        break;
                }

                total_time += esolver->time;
                current_ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
                esolver->subspace_retry_final_size = current_ss;

#ifdef USE_MPI
                /*
                 * Every rank must make the same retry decision.  A divergent
                 * eigensolver status would otherwise split the communicator
                 * between retry and return paths.
                 */
                local_status = esolver->retcode;
                err = MPI_Allreduce(&local_status,&status_min,1,LIS_MPI_INT,MPI_MIN,A->comm);
                if( err==MPI_SUCCESS )
                        err = MPI_Allreduce(&local_status,&status_max,1,LIS_MPI_INT,MPI_MAX,A->comm);
                if( err!=MPI_SUCCESS || status_min!=status_max )
                {
                        LIS_SETERR(LIS_ERR_ILL_ARG,
                                "inconsistent eigensolver retry status across MPI ranks\n");
                        result = LIS_ERR_ILL_ARG;
                        break;
                }
#endif

                if( esolver->retcode!=LIS_MAXITER ||
                    attempt>=esolver->subspace_retry_max ||
                    A==NULL || current_ss>=A->gn )
                {
                        break;
                }

                nesolver = esolver->options[LIS_EOPTIONS_ESOLVER];
                if( !lis_esolver_subspace_retry_supported(nesolver) )
                {
                        break;
                }

                grown = (LIS_REAL)current_ss * esolver->subspace_retry_growth;
                if( !isfinite((double)grown) )
                {
                        /* An infinite growth request simply saturates at n. */
                        next_ss = A->gn;
                }
                else if( grown>=(LIS_REAL)A->gn )
                {
                        /* Bound in floating point before converting to LIS_INT. */
                        next_ss = A->gn;
                }
                else
                {
                        next_ss = (LIS_INT)ceil((double)grown);
                }

                if( esolver->subspace_retry_min>A->gn )
                        next_ss = A->gn;
                else if( next_ss<esolver->subspace_retry_min )
                        next_ss = esolver->subspace_retry_min;

                if( next_ss<=current_ss )
                {
                        if( current_ss>=A->gn ) break;
                        next_ss = current_ss+1;
                }

                if( next_ss>A->gn ) next_ss = A->gn;
                if( next_ss<=current_ss ) break;

                lis_esolver_work_destroy(esolver);
                esolver->options[LIS_EOPTIONS_SUBSPACE] = next_ss;
                esolver->subspace_retry_attempts++;

                /*
                 * x intentionally keeps the result of the previous attempt.
                 * The enlarged-subspace retry therefore uses a warm start.
                 */
        }

        /*
         * The retry policy is per-call.  Preserve the user's configured
         * subspace for the next lis_gesolve() while keeping the actual final
         * size in subspace_retry_final_size.  iter[] and resid[] describe the
         * final attempt; time covers the complete retry sequence.
         */
        esolver->options[LIS_EOPTIONS_SUBSPACE] = original_ss;
        esolver->time = total_time;

        LIS_DEBUG_FUNC_OUT;
        return result;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_optionC"
LIS_INT lis_esolver_set_optionC(LIS_ESOLVER esolver)
{
	LIS_INT err;    
	LIS_ARGS p;

	LIS_DEBUG_FUNC_IN;

	p = cmd_args->next;
	while( p!=cmd_args )
	{
		err = lis_esolver_set_option2(p->arg1,p->arg2,esolver);
		if( err )
		  {
		    lis_esolver_work_destroy(esolver);
		    esolver->retcode = err;
		    LIS_DEBUG_FUNC_OUT;
		    return err;
		  }
		p = p->next;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option"
LIS_INT lis_esolver_set_option(char *text, LIS_ESOLVER esolver)
{
	LIS_INT err;   
	LIS_ARGS args,p;

	LIS_DEBUG_FUNC_IN;

	lis_text2args(text,&args);
	p = args->next;
	while( p!=args )
	{
		err = lis_esolver_set_option2(p->arg1,p->arg2,esolver);
		if( err )
		  {
		    lis_esolver_work_destroy(esolver);
		    esolver->retcode = err;
		    LIS_DEBUG_FUNC_OUT;
		    return err;
		  }		
		p = p->next;
	}
	lis_args_free(args);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option2"
LIS_INT lis_esolver_set_option2(char* arg1, char *arg2, LIS_ESOLVER esolver)
{
	LIS_INT err;  
	LIS_INT i;
	double double_value;

	LIS_DEBUG_FUNC_IN;

	err = 0;

        /*
         * Native FEAST options are kept outside LIS_EOPTIONS_*
         * and LIS_EPARAMS_* to preserve their public identifiers.
         */
        if( strncmp(arg1,"-feast_",7)==0 )
        {
                err = lis_efeast_set_option(
                        arg1,
                        arg2,
                        esolver);

                if( err )
                {
                        lis_esolver_work_destroy(esolver);
                        esolver->retcode = err;
                }

                LIS_DEBUG_FUNC_OUT;
                return err;
        }

	/*
	 * Adaptive subspace retry options are kept outside LIS_EOPTIONS_* so
	 * existing public numeric option identifiers remain stable.
	 */
	if( strcmp(arg1,"-ess_retry")==0 )
	{
		LIS_INT value = (LIS_INT)atoi(arg2);
		err = lis_esolver_set_subspace_retry(
			esolver,value,esolver->subspace_retry_growth,
			esolver->subspace_retry_min);
		LIS_DEBUG_FUNC_OUT;
		return err;
	}
	if( strcmp(arg1,"-ess_retry_growth")==0 )
	{
		LIS_REAL value = (LIS_REAL)atof(arg2);
		err = lis_esolver_set_subspace_retry(
			esolver,esolver->subspace_retry_max,value,
			esolver->subspace_retry_min);
		LIS_DEBUG_FUNC_OUT;
		return err;
	}
	if( strcmp(arg1,"-ess_retry_min")==0 )
	{
		LIS_INT value = (LIS_INT)atoi(arg2);
		err = lis_esolver_set_subspace_retry(
			esolver,esolver->subspace_retry_max,
			esolver->subspace_retry_growth,value);
		LIS_DEBUG_FUNC_OUT;
		return err;
	}

	/*
	 * Block size is kept outside LIS_EOPTIONS_* so existing public
	 * eigensolver option and parameter indices remain unchanged.
	 *
	 * A value of one preserves the historical scalar eigensolver
	 * behavior.  Values greater than one are currently reserved for
	 * the standard Block Lanczos path.
	 */
	if( strcmp(arg1,"-eblock")==0 )
	{
		esolver->block_size = atoi(arg2);

		if( esolver->block_size<1 )
		{
			LIS_SETERR(
				LIS_ERR_ILL_ARG,
				"Parameter -eblock must be at least 1\n");
			esolver->retcode = LIS_ERR_ILL_ARG;
			LIS_DEBUG_FUNC_OUT;
			return LIS_ERR_ILL_ARG;
		}

		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}

	/*
	 * Spectral transformations are kept outside LIS_EOPTIONS_* so
	 * existing public option and parameter indices remain unchanged.
	 */
	/*
	 * Requested eigenpair count is independent from the Lanczos
	 * subspace dimension. A zero internal value preserves the
	 * historical behavior nev == ss.
	 */
	if( strcmp(arg1,"-nev")==0 )
	{
		LIS_INT parsed_nev = (LIS_INT)atoi(arg2);

		if( parsed_nev<0 )
		{
			LIS_SETERR(
				LIS_ERR_ILL_ARG,
				"Parameter -nev must not be negative\n");
			esolver->retcode = LIS_ERR_ILL_ARG;
			LIS_DEBUG_FUNC_OUT;
			return LIS_ERR_ILL_ARG;
		}

                /*
                 * Historical text options may specify -nev before -ss.
                 * Preserve that ordering and validate the final relation
                 * at solve time; the typed setter validates immediately.
                 */
                esolver->requested_nev=parsed_nev;
                esolver->requested_nev_active=LIS_TRUE;
                esolver->nev=parsed_nev;

		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}

	if( strcmp(arg1,"-etransform")==0 )
	{
		if( strcmp(arg2,"auto")==0 || strcmp(arg2,"0")==0 )
		{
			esolver->transform = LIS_EIGEN_TRANSFORM_AUTO;
		}
		else if( strcmp(arg2,"shift-invert")==0 ||
		         strcmp(arg2,"shiftinvert")==0 ||
		         strcmp(arg2,"shift_invert")==0 ||
		         strcmp(arg2,"1")==0 )
		{
			esolver->transform =
				LIS_EIGEN_TRANSFORM_SHIFT_INVERT;
		}
		else if( strcmp(arg2,"buckling")==0 ||
		         strcmp(arg2,"2")==0 )
		{
			esolver->transform =
				LIS_EIGEN_TRANSFORM_BUCKLING;
		}
		else
		{
			LIS_SETERR(
				LIS_ERR_ILL_ARG,
				"Parameter -etransform is not correct\n");
			esolver->retcode = LIS_ERR_ILL_ARG;
			LIS_DEBUG_FUNC_OUT;
			return LIS_ERR_ILL_ARG;
		}

		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	
	for(i=0;i<LIS_ESOLVER_OPTION_LEN;i++)
	{
		if( strcmp(arg1, LIS_ESOLVER_OPTNAME[i])==0 )
		{
			switch( LIS_ESOLVER_OPTACT[i] )
			{
			case LIS_EOPTIONS_ESOLVER:
			  err = lis_esolver_set_option_esolver(arg2,esolver);
			  break;
                        case LIS_EOPTIONS_OUTPUT:
			  err = lis_esolver_set_option_print(arg2,esolver);
			  break;
                        case LIS_EOPTIONS_INITGUESS_ONES:
			  err = lis_esolver_set_option_truefalse(arg2,LIS_EOPTIONS_INITGUESS_ONES,esolver);
			  break;
			case LIS_EOPTIONS_INNER_ESOLVER:
			  err = lis_esolver_set_option_iesolver(arg2,esolver);
			  break;
			case LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER:
			  err = lis_esolver_set_option_igesolver(arg2,esolver);
			  break;
			case LIS_EOPTIONS_STORAGE:
			  err = lis_esolver_set_option_storage(arg2,esolver);
			  break;
                        case LIS_EOPTIONS_PRECISION:
			  err = lis_esolver_set_option_eprecision(arg2,LIS_EOPTIONS_PRECISION,esolver);
			  break;
                        case LIS_EOPTIONS_RVAL:
			  err = lis_esolver_set_option_truefalse(arg2,LIS_EOPTIONS_RVAL,esolver);
			  break;
			default:
			  if( LIS_ESOLVER_OPTACT[i] < LIS_EOPTIONS_LEN )
			    {
#ifdef _LONG__LONG
			      sscanf(arg2, "%lld", &esolver->options[LIS_ESOLVER_OPTACT[i]]);
#else
			      sscanf(arg2, "%d", &esolver->options[LIS_ESOLVER_OPTACT[i]]);
#endif
			    }
			  else
			    {
			      sscanf(arg2, "%lg", &double_value);
			      esolver->params[LIS_ESOLVER_OPTACT[i]-LIS_EOPTIONS_LEN] = double_value;
			    }
			  break;
			}
		}
		if( err )
		  {
		    lis_esolver_work_destroy(esolver);
		    esolver->retcode = err;
		    LIS_DEBUG_FUNC_OUT;
		    return err;
		  }
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_esolver"
LIS_INT lis_esolver_set_option_esolver(char *argv, LIS_ESOLVER esolver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &esolver->options[LIS_EOPTIONS_ESOLVER]);
#else
		sscanf(argv, "%d", &esolver->options[LIS_EOPTIONS_ESOLVER]);
#endif
	}
	else
	{
		for(i=0;i<LIS_ESOLVER_LEN;i++)
		{
			if( strcmp(argv,lis_esolver_atoi[i])==0 )
			{
				esolver->options[LIS_EOPTIONS_ESOLVER] = i+1;
				break;
			}
			else if( i==LIS_ESOLVER_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_iesolver"
LIS_INT lis_esolver_set_option_iesolver(char *argv, LIS_ESOLVER esolver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &esolver->options[LIS_EOPTIONS_INNER_ESOLVER]);
#else
		sscanf(argv, "%d", &esolver->options[LIS_EOPTIONS_INNER_ESOLVER]);
#endif
	}
	else
	{
		for(i=0;i<LIS_ESOLVER_LEN;i++)
		{
			if( strcmp(argv,lis_esolver_atoi[i])==0 )
			{
				esolver->options[LIS_EOPTIONS_INNER_ESOLVER] = i+1;
				break;
			}
			else if( i==LIS_ESOLVER_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_INNER_ESOLVER is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_igesolver"
LIS_INT lis_esolver_set_option_igesolver(char *argv, LIS_ESOLVER esolver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &esolver->options[LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER]);
#else
		sscanf(argv, "%d", &esolver->options[LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER]);
#endif
	}
	else
	{
		for(i=0;i<LIS_ESOLVER_LEN;i++)
		{
			if( strcmp(argv,lis_esolver_atoi[i])==0 )
			{
				esolver->options[LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER] = i+1;
				break;
			}
			else if( i==LIS_ESOLVER_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_print"
LIS_INT lis_esolver_set_option_print(char *argv, LIS_ESOLVER esolver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='3' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &esolver->options[LIS_EOPTIONS_OUTPUT]);
#else
		sscanf(argv, "%d", &esolver->options[LIS_EOPTIONS_OUTPUT]);
#endif
	}
	else
	{
		for(i=0;i<LIS_EPRINT_LEN;i++)
		{
			if( strcmp(argv,lis_eprint_atoi[i])==0 )
			{
				esolver->options[LIS_EOPTIONS_OUTPUT] = i;
				break;
			}
			else if( i==LIS_EPRINT_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_OUTPUT is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_truefalse"
LIS_INT lis_esolver_set_option_truefalse(char *argv, LIS_INT opt, LIS_ESOLVER esolver)
{
        LIS_INT i;

        LIS_DEBUG_FUNC_IN;

        if( argv[0]>='0' && argv[0]<='1' )
        {
#ifdef _LONG__LONG
                sscanf(argv, "%lld", &esolver->options[opt]);
#else
                sscanf(argv, "%d", &esolver->options[opt]);
#endif
        }
        else
        {
                for(i=0;i<LIS_TRUEFALSE_LEN;i++)
                {
                        if( strcmp(argv,lis_etruefalse_atoi[i])==0 )
                        {
                                esolver->options[opt] = i;
                                break;
                        }
			else if( i==LIS_TRUEFALSE_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_TRUEFALSE is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
                }
        }
        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_eprecision"
LIS_INT lis_esolver_set_option_eprecision(char *argv, LIS_INT opt, LIS_ESOLVER esolver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='1' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &esolver->options[opt]);
#else
		sscanf(argv, "%d", &esolver->options[opt]);
#endif
	}
	else
	{
		for(i=0;i<LIS_PRECISION_LEN;i++)
		{
			if( strcmp(argv,lis_eprecision_atoi[i])==0 )
			{
				esolver->options[opt] = i;
				break;
			}
			else if( i==LIS_PRECISION_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_PRECISION is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_storage"
LIS_INT lis_esolver_set_option_storage(char *argv, LIS_ESOLVER esolver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &esolver->options[LIS_EOPTIONS_STORAGE]);
#else
		sscanf(argv, "%d", &esolver->options[LIS_EOPTIONS_STORAGE]);
#endif
	}
	else
	{
		for(i=0;i<LIS_ESTORAGE_LEN;i++)
		{
			if( strcmp(argv,lis_estorage_atoi[i])==0 )
			{
				esolver->options[LIS_EOPTIONS_STORAGE] = i+1;
				break;
			}
			else if( i==LIS_ESTORAGE_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_STORAGE is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_iter"
LIS_INT lis_esolver_get_iter(LIS_ESOLVER esolver, LIS_INT *iter)
{
	LIS_DEBUG_FUNC_IN;

	*iter = esolver->iter[0];

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_iterex"
LIS_INT lis_esolver_get_iterex(LIS_ESOLVER esolver, LIS_INT *iter, LIS_INT *iter_double, LIS_INT *iter_quad)
{
	LIS_DEBUG_FUNC_IN;

	*iter = esolver->iter[0];
	*iter_double = esolver->iter2[0];
	*iter_quad = esolver->iter[0] - esolver->iter2[0];

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_time"
LIS_INT lis_esolver_get_time(LIS_ESOLVER esolver, double *time)
{
	LIS_DEBUG_FUNC_IN;

	*time  = esolver->time;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_timeex"
LIS_INT lis_esolver_get_timeex(LIS_ESOLVER esolver, double *time, double *itime, double *ptime, double *p_c_time, double *p_i_time)
{
	LIS_DEBUG_FUNC_IN;

	*time  = esolver->time;
	if( itime ) *itime = esolver->itime;
	if( ptime ) *ptime = esolver->ptime;
	if( p_c_time ) *p_c_time = esolver->p_c_time;
	if( p_i_time ) *p_i_time = esolver->p_i_time;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_residualnorm"
LIS_INT lis_esolver_get_residualnorm(LIS_ESOLVER esolver, LIS_REAL *residual)
{
	LIS_DEBUG_FUNC_IN;

	*residual  = esolver->resid[0];

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_status"
LIS_INT lis_esolver_get_status(LIS_ESOLVER esolver, LIS_INT *status)
{
	LIS_DEBUG_FUNC_IN;

	*status = esolver->retcode;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_rhistory"
LIS_INT lis_esolver_get_rhistory(LIS_ESOLVER esolver, LIS_VECTOR v)
{
	LIS_INT	i,n,maxiter;

	maxiter = esolver->iter[0]+1;
        if( esolver->retcode!=LIS_SUCCESS )
	  {
	    maxiter--;
	  }
	n = _min(v->n,maxiter);
	for(i=0;i<n;i++)
	{
		v->value[i] = esolver->rhistory[i];
	}
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_evalues"
LIS_INT lis_esolver_get_evalues(LIS_ESOLVER esolver, LIS_VECTOR v)
{
	LIS_INT i,ii,ss,n,gn,is,ie;

	LIS_DEBUG_FUNC_IN;

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_SI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_LI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_AI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GSI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GLI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GAI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_FEAST)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %D (Set Subspace, Lanczos, Arnoldi, Generalized Subspace, Generalized Lanczos, or Generalized Arnoldi)\n", esolver->options[LIS_EOPTIONS_ESOLVER]);
		return LIS_ERR_ILL_ARG;
	}

	ss = lis_esolver_result_count(esolver);
	lis_vector_set_size(v,0,ss);
	lis_vector_get_size(v,&n,&gn);
	lis_vector_get_range(v,&is,&ie);
	
	for(i=0;i<n;i++)
	{
	  ii=i;
	  if( v->origin ) ii++;
	  lis_vector_set_value(LIS_INS_VALUE,ii+is,esolver->evalue[i+is],v);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_specific_evalue"
LIS_INT lis_esolver_get_specific_evalue(LIS_ESOLVER esolver, LIS_INT mode, LIS_SCALAR *evalue)
{
   LIS_INT count;
	LIS_DEBUG_FUNC_IN;

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_SI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_LI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_AI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GSI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GLI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GAI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_FEAST)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %D (Set Subspace, Lanczos, Arnoldi, Generalized Subspace, Generalized Lanczos, or Generalized Arnoldi)\n", esolver->options[LIS_EOPTIONS_ESOLVER]);
		return LIS_ERR_ILL_ARG;
	}

	   count = lis_esolver_result_count(esolver);

   if( mode<0 || mode>=count )
   {
           LIS_SETERR2(
                   LIS_ERR_ILL_ARG,
                   "Parameter mode is %D "
                   "(Set between 0 and %D)\n",
                   mode,
                   count-1);
           return LIS_ERR_ILL_ARG;
   }

*evalue = esolver->evalue[mode];

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_evectors"
LIS_INT lis_esolver_get_evectors(LIS_ESOLVER esolver, LIS_MATRIX M)
{
	LIS_INT i,ii,j,jj,n,gn,is,ie,js;
	LIS_INT ss,lis_esolver_evector_size;

	LIS_DEBUG_FUNC_IN;

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_SI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_LI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_AI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GSI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GLI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GAI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_FEAST)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %D (Set Subspace, Lanczos, Arnoldi, Generalized Subspace, Generalized Lanczos, or Generalized Arnoldi)\n", esolver->options[LIS_EOPTIONS_ESOLVER]);
		return LIS_ERR_ILL_ARG;
	}

	ss = lis_esolver_result_count(esolver);

   if( ss<1 ||
       esolver->evector==NULL ||
       esolver->evector[0]==NULL )
   {
           LIS_SETERR(
                   LIS_ERR_ILL_ARG,
                   "Eigenvectors are not available\n");
           return LIS_ERR_ILL_ARG;
   }

	lis_esolver_evector_size = esolver->evector[0]->gn; 
	lis_matrix_set_size(M,0,lis_esolver_evector_size);
	lis_matrix_get_size(M,&n,&gn);
	lis_matrix_get_range(M,&is,&ie);
	js=0;
	if( esolver->evector[0]->origin ) 
	  {
	    is++;
	    js++;
	  }
	for(j=0;j<ss;j++)
	  {
	    for(i=0;i<n;i++)
	      {
		ii=i+is;
		jj=j+js;
		lis_matrix_set_value(LIS_INS_VALUE,ii,jj,esolver->evector[j]->value[i],M);
	      }
	  }
	lis_matrix_set_type(M,LIS_MATRIX_COO);
	lis_matrix_assemble(M);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_specific_evector"
LIS_INT lis_esolver_get_specific_evector(LIS_ESOLVER esolver, LIS_INT mode, LIS_VECTOR x)
{
   LIS_INT count;
	LIS_DEBUG_FUNC_IN;

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_SI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_LI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_AI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GSI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GLI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GAI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_FEAST)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %D (Set Subspace, Lanczos, Arnoldi, Generalized Subspace, Generalized Lanczos, or Generalized Arnoldi)\n", esolver->options[LIS_EOPTIONS_ESOLVER]);
		return LIS_ERR_ILL_ARG;
	}

	   count = lis_esolver_result_count(esolver);

   if( mode<0 || mode>=count )
   {
           LIS_SETERR2(
                   LIS_ERR_ILL_ARG,
                   "Parameter mode is %D "
                   "(Set between 0 and %D)\n",
                   mode,
                   count-1);
           return LIS_ERR_ILL_ARG;
   }

   if( esolver->evector==NULL ||
       esolver->evector[mode]==NULL )
   {
           LIS_SETERR1(
                   LIS_ERR_ILL_ARG,
                   "Eigenvector mode %D is not available\n",
                   mode);
           return LIS_ERR_ILL_ARG;
   }

lis_vector_copy(esolver->evector[mode],x);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_residualnorms"
LIS_INT lis_esolver_get_residualnorms(LIS_ESOLVER esolver, LIS_VECTOR v)
{
	LIS_INT i,ii,ss,n,gn,is,ie;

	LIS_DEBUG_FUNC_IN;

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_SI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_LI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_AI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GSI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GLI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GAI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_FEAST)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %D (Set Subspace, Lanczos, Arnoldi, Generalized Subspace, Generalized Lanczos, or Genralized Arnoldi)\n", esolver->options[LIS_EOPTIONS_ESOLVER]);
		return LIS_ERR_ILL_ARG;
	}

	ss = lis_esolver_result_count(esolver);
	lis_vector_set_size(v,0,ss);
	lis_vector_get_size(v,&n,&gn);
	lis_vector_get_range(v,&is,&ie);
	
	for(i=0;i<n;i++)
	{
	  ii=i;
	  if( v->origin ) ii++;
	  lis_vector_set_value(LIS_INS_VALUE,ii+is,esolver->resid[i+is],v);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_specific_residualnorm"
LIS_INT lis_esolver_get_specific_residualnorm(LIS_ESOLVER esolver, LIS_INT mode, LIS_REAL *residual)
{
   LIS_INT count;
	LIS_DEBUG_FUNC_IN;

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_SI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_LI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_AI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GSI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GLI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GAI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_FEAST)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %D (Set Subspace, Lanczos, Arnoldi, Generalized Subspace, Generalized Lanczos, or Genralized Arnoldi)\n", esolver->options[LIS_EOPTIONS_ESOLVER]);
		return LIS_ERR_ILL_ARG;
	}

	   count = lis_esolver_result_count(esolver);

   if( mode<0 || mode>=count )
   {
           LIS_SETERR2(
                   LIS_ERR_ILL_ARG,
                   "Parameter mode is %D "
                   "(Set between 0 and %D)\n",
                   mode,
                   count-1);
           return LIS_ERR_ILL_ARG;
   }

*residual  = esolver->resid[mode];

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_iters"
LIS_INT lis_esolver_get_iters(LIS_ESOLVER esolver, LIS_VECTOR v)
{
	LIS_INT i,ii,ss,n,gn,is,ie;

	LIS_DEBUG_FUNC_IN;

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_SI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_LI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_AI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GSI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GLI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GAI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_FEAST)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %D (Set Subspace, Lanczos, Arnoldi, Generalized Subspace, Generalized Lanczos, or Genralized Arnoldi)\n", esolver->options[LIS_EOPTIONS_ESOLVER]);
		return LIS_ERR_ILL_ARG;
	}

	ss = lis_esolver_result_count(esolver);
	lis_vector_set_size(v,0,ss);
	lis_vector_get_size(v,&n,&gn);
	lis_vector_get_range(v,&is,&ie);
	
	for(i=0;i<n;i++)
	{
	  ii=i;
	  if( v->origin ) ii++;
	  lis_vector_set_value(LIS_INS_VALUE,ii+is,esolver->iter[i+is],v);
	}

	v->intvalue = 1;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_specific_iter"
LIS_INT lis_esolver_get_specific_iter(LIS_ESOLVER esolver, LIS_INT mode, LIS_INT *iter)
{
   LIS_INT count;
	LIS_DEBUG_FUNC_IN;

	if ( esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_SI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_LI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_AI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GSI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GLI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_GAI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_FEAST)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %D (Set Subspace, Lanczos, Arnoldi, Generalized Subspace, Generalized Lanczos, or Genralized Arnoldi)\n", esolver->options[LIS_EOPTIONS_ESOLVER]);
		return LIS_ERR_ILL_ARG;
	}

	   count = lis_esolver_result_count(esolver);

   if( mode<0 || mode>=count )
   {
           LIS_SETERR2(
                   LIS_ERR_ILL_ARG,
                   "Parameter mode is %D "
                   "(Set between 0 and %D)\n",
                   mode,
                   count-1);
           return LIS_ERR_ILL_ARG;
   }

*iter = esolver->iter[mode];

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_esolver"
LIS_INT lis_esolver_get_esolver(LIS_ESOLVER esolver, LIS_INT *nesol)
{
	LIS_DEBUG_FUNC_IN;

	*nesol = esolver->options[LIS_EOPTIONS_ESOLVER];

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_esolvername"
LIS_INT lis_esolver_get_esolvername(LIS_INT esolver, char *esolvername)
{
	LIS_DEBUG_FUNC_IN;

	if( esolver < 1 || esolver > LIS_ESOLVER_LEN )
	{
		esolvername = NULL;
		return LIS_FAILS;
	}
	strcpy(esolvername,lis_esolvername[esolver]);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

