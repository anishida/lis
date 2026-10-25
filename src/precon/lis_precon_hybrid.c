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
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_precon_create
 * lis_psolve
 * lis_psolveh
 ************************************************/

extern LIS_SOLVER_MALLOC_WORK lis_solver_malloc_work[];
extern LIS_SOLVER_EXECUTE lis_solver_execute[];

/*
 * Hybrid is LIS's native nested solver/preconditioner composition path:
 *
 *     outer solver -> Hybrid PC -> inner solver -> inner PC
 *
 * Keep the option-driven implementation one level deep.  The current option
 * namespace contains only one set of -hybrid_* parameters, so recursive
 * Hybrid composition has no independent configuration and is rejected.
 */
#undef __FUNC__
#define __FUNC__ "lis_precon_hybrid_validate_config"
static LIS_INT lis_precon_hybrid_validate_config(LIS_SOLVER solver)
{
	LIS_INT nsolver;
	LIS_INT maxiter;
	LIS_INT inner_precon;

	if( solver==NULL || solver->A==NULL )
	{
		return LIS_ERR_ILL_ARG;
	}

	nsolver = solver->options[LIS_OPTIONS_PSOLVER];
	maxiter = solver->options[LIS_OPTIONS_PMAXITER];
	inner_precon = solver->options[LIS_OPTIONS_PPRECON];

	if( nsolver<1 || nsolver>LIS_SOLVER_LEN )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"invalid Hybrid inner solver\n");
		return LIS_ERR_ILL_ARG;
	}

	if( maxiter<1 )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"invalid Hybrid inner maximum iteration count\n");
		return LIS_ERR_ILL_ARG;
	}

	if( inner_precon==LIS_PRECON_TYPE_HYBRID )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"recursive Hybrid preconditioner is not supported\n");
		return LIS_ERR_ILL_ARG;
	}

	return LIS_SUCCESS;
}

/*
 * lis_solver_destroy() owns solver work and residual-history storage, while
 * Hybrid owns the inner solution vector and inner preconditioner separately.
 */
static void lis_precon_hybrid_destroy_partial(LIS_SOLVER psolver)
{
	if( psolver==NULL ) return;

	if( psolver->x!=NULL )
	{
		lis_vector_destroy(psolver->x);
		psolver->x = NULL;
	}

	if( psolver->precon!=NULL )
	{
		lis_precon_destroy(psolver->precon);
		psolver->precon = NULL;
	}

	lis_solver_destroy(psolver);
}

#undef __FUNC__
#define __FUNC__ "lis_precon_create_hybrid"
LIS_INT lis_precon_create_hybrid(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT nsolver,maxiter,precision;
	LIS_INT err;
	LIS_REAL *rhistory;
	LIS_VECTOR xx;
	LIS_SOLVER psolver;
	LIS_MATRIX A;
	LIS_PRECON pprecon;

	LIS_DEBUG_FUNC_IN;

	if( solver==NULL || precon==NULL ) return LIS_ERR_ILL_ARG;

	err = lis_precon_hybrid_validate_config(solver);
	if( err )
	{
		solver->retcode = err;
		return err;
	}

	A = solver->A;
	psolver = NULL;
	pprecon = NULL;
	xx = NULL;
	rhistory = NULL;

	err = lis_solver_create(&psolver);
	if( err )
	{
		solver->retcode = err;
		return err;
	}

	psolver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN] =
		solver->params[LIS_PARAMS_PRESID-LIS_OPTIONS_LEN];
	psolver->params[LIS_PARAMS_SSOR_OMEGA-LIS_OPTIONS_LEN] =
		solver->params[LIS_PARAMS_POMEGA-LIS_OPTIONS_LEN];
	psolver->options[LIS_OPTIONS_MAXITER] = solver->options[LIS_OPTIONS_PMAXITER];
	psolver->options[LIS_OPTIONS_ELL] = solver->options[LIS_OPTIONS_PELL];
	psolver->options[LIS_OPTIONS_RESTART] = solver->options[LIS_OPTIONS_PRESTART];
	psolver->options[LIS_OPTIONS_OUTPUT] = 0;
	psolver->options[LIS_OPTIONS_SOLVER] = solver->options[LIS_OPTIONS_PSOLVER];
	psolver->options[LIS_OPTIONS_PRECON] = solver->options[LIS_OPTIONS_PPRECON];
	psolver->options[LIS_OPTIONS_INITGUESS_ZEROS] =
		solver->options[LIS_OPTIONS_INITGUESS_ZEROS];
	psolver->options[LIS_OPTIONS_PRECISION] = solver->options[LIS_OPTIONS_PRECISION];
	psolver->A = solver->A;
	psolver->Ah = NULL;
	psolver->precision = solver->precision;

	/*
	 * Propagate solver configuration used by GCRO-DR when it is selected
	 * as the Hybrid inner Krylov solver.  Only policy/configuration is
	 * copied here; recycle vectors and lifecycle state remain owned by the
	 * inner solver.
	 *
	 * This keeps the inner solver independently owned while allowing
	 * -gcrodr_keep and the recycle-quality thresholds to configure
	 * GCRO-DR selected through -hybrid_i gcrodr.
	 */
	psolver->recycle_target_dim =
		solver->recycle_target_dim;
	psolver->recycle_keep_threshold =
		solver->recycle_keep_threshold;
	psolver->recycle_reject_threshold =
		solver->recycle_reject_threshold;

	psolver->options[LIS_OPTIONS_MAXITER_NO_IMP] =
		solver->options[LIS_OPTIONS_MAXITER_NO_IMP];

	/*
	 * Propagate solver-owned context required by inner preconditioners.
	 * lis_solver_set_near_nullspace() creates inner-solver-owned copies, so
	 * both solver lifetimes remain independent.  This is required for MFAMG
	 * when it is selected through -hybrid_p mfamg.
	 */
	psolver->mfamg_block_size = solver->mfamg_block_size;

	if( solver->near_nullspace_dim>0 && solver->near_nullspace!=NULL )
	{
		err = lis_solver_set_near_nullspace(
			psolver,
			solver->near_nullspace_dim,
			solver->near_nullspace);
		if( err ) goto fail;
	}

	nsolver = psolver->options[LIS_OPTIONS_SOLVER];
	maxiter = psolver->options[LIS_OPTIONS_MAXITER];
	precision = psolver->options[LIS_OPTIONS_PRECISION];
	A = psolver->A;

#ifndef USE_QUAD_PRECISION
	err = lis_vector_duplicate(A,&xx);
#else
	if( precision==LIS_PRECISION_DOUBLE )
		err = lis_vector_duplicate(A,&xx);
	else
		err = lis_vector_duplicateex(LIS_PRECISION_QUAD,A,&xx);
#endif
	if( err ) goto fail;

	psolver->x = xx;
	xx = NULL;

	rhistory = (LIS_REAL *)lis_malloc(
		(maxiter+2)*sizeof(LIS_REAL),
		"lis_precon_create_hybrid::rhistory");
	if( rhistory==NULL )
	{
		LIS_SETERR_MEM((maxiter+2)*sizeof(LIS_REAL));
		err = LIS_OUT_OF_MEMORY;
		goto fail;
	}
	psolver->rhistory = rhistory;
	rhistory = NULL;

	err = lis_precon_create(psolver,&pprecon);
	if( err )
	{
		/* lis_precon_create() cleans a partially-created preconditioner. */
		pprecon = NULL;
		goto fail;
	}
	psolver->precon = pprecon;
	pprecon = NULL;

	err = lis_solver_malloc_work[nsolver](psolver);
	if( err ) goto fail;

	precon->solver = psolver;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;

fail:
	if( xx!=NULL ) lis_vector_destroy(xx);
	if( pprecon!=NULL ) lis_precon_destroy(pprecon);
	if( rhistory!=NULL ) lis_free(rhistory);
	lis_precon_hybrid_destroy_partial(psolver);
	solver->retcode = err;

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_hybrid"
LIS_INT lis_psolve_hybrid(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_VECTOR xx;
	LIS_SOLVER solver2;
	LIS_INT nsolver;
	LIS_INT err;
	LIS_PRECON precon;

	LIS_DEBUG_FUNC_IN;

	if( solver==NULL || solver->precon==NULL || B==NULL || X==NULL )
		return LIS_ERR_ILL_ARG;

	precon = solver->precon;
	solver2 = precon->solver;
	if( solver2==NULL || solver2->x==NULL ) return LIS_ERR_ILL_ARG;

	xx = solver2->x;
	nsolver = solver2->options[LIS_OPTIONS_SOLVER];
	solver2->b = B;

	if( nsolver<1 || nsolver>LIS_SOLVER_LEN ) return LIS_ERR_ILL_ARG;

	if( solver2->options[LIS_OPTIONS_INITGUESS_ZEROS] )
	{
#ifdef USE_QUAD_PRECISION
		if( solver->precision==LIS_PRECISION_DEFAULT )
		{
#endif
			err = lis_vector_set_all(0,xx);
#ifdef USE_QUAD_PRECISION
		}
		else
		{
			err = lis_vector_set_allex_nm(0,xx);
		}
#endif
	}
	else
	{
#ifdef USE_QUAD_PRECISION
		if( solver->precision==LIS_PRECISION_DEFAULT )
		{
#endif
			err = lis_vector_copy(B,xx);
#ifdef USE_QUAD_PRECISION
		}
		else
		{
			err = lis_vector_copyex_mm(B,xx);
		}
#endif
	}

	if( err )
	{
		LIS_DEBUG_FUNC_OUT;
		return err;
	}

	err = lis_solver_execute[nsolver](solver2);
	/*
	 * Hybrid is an iterative preconditioner.  Reaching the configured
	 * inner iteration cap is therefore a usable truncated application,
	 * not a fatal preconditioner failure: the inner Krylov solver has
	 * already updated solver2->x with its best available iterate.
	 *
	 * Preserve fail-closed behavior for real failures such as BREAKDOWN,
	 * invalid configuration, allocation errors and NOT_IMPLEMENTED.
	 */
	if( err && err!=LIS_MAXITER )
	{
		LIS_DEBUG_FUNC_OUT;
		return err;
	}

#ifdef USE_QUAD_PRECISION
	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
#endif
		err = lis_vector_copy(solver2->x,X);
#ifdef USE_QUAD_PRECISION
	}
	else
	{
		err = lis_vector_copyex_mm(solver2->x,X);
	}
#endif

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_hybrid"
LIS_INT lis_psolveh_hybrid(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_DEBUG_FUNC_IN;

	if( solver==NULL || solver->precon==NULL || B==NULL || X==NULL )
		return LIS_ERR_ILL_ARG;

	/*
	 * The legacy implementation attempted to obtain A^H by swapping the
	 * global LIS_MATVEC/LIS_MATVECH function pointers.  Krylov solvers such
	 * as GMRES call lis_matvec() directly, so that swap does not change the
	 * operator used by the inner solve and silently produces a wrong result
	 * for nonsymmetric matrices.
	 *
	 * Do not return a mathematically incorrect preconditioner application.
	 * Proper Hybrid adjoint support requires a solver-local adjoint operator
	 * (or an explicitly owned transposed matrix) and is intentionally left
	 * for a separate change.
	 */
	LIS_SETERR(
		LIS_ERR_NOT_IMPLEMENTED,
		"Hybrid adjoint preconditioning is not implemented for nested solves\n");

	LIS_DEBUG_FUNC_OUT;
	return LIS_ERR_NOT_IMPLEMENTED;
}
