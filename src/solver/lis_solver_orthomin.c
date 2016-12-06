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

#define NWORK 3
/************************************************
 * lis_orthomin_check_params
 * lis_orthomin_malloc_work
 * lis_orthomin
 ************************************************/
#undef __FUNC__
#define __FUNC__ "lis_orthomin_check_params"
LIS_INT lis_orthomin_check_params(LIS_SOLVER solver)
{
	LIS_INT restart;

	LIS_DEBUG_FUNC_IN;

	restart = solver->options[LIS_OPTIONS_RESTART];
	if( restart<0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_RESTART(=%D) is less than 0\n",restart);
		return LIS_ERR_ILL_ARG;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_orthomin_malloc_work"
LIS_INT lis_orthomin_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,restart,worklen,err;

	LIS_DEBUG_FUNC_IN;

	restart = solver->options[LIS_OPTIONS_RESTART];
	worklen = NWORK + 3*(restart+1);
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_orthomin_malloc_work::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicate(solver->A,&work[i]);
			if( err ) break;
		}
	}
	else
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,solver->A,&work[i]);
			if( err ) break;
			memset(work[i]->value_lo,0,solver->A->np*sizeof(LIS_SCALAR));
		}
	}
	if( i<worklen )
	{
		for(j=0;j<i;j++) lis_vector_destroy(work[j]);
		lis_free(work);
		return err;
	}
	solver->worklen = worklen;
	solver->work    = work;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_orthomin"
LIS_INT lis_orthomin(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_PRECON M;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld,*p,*ap,*aptld;
	LIS_SCALAR *dotsave;
	LIS_SCALAR alpha, beta;
	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output,conv;
	double time,ptime;

	LIS_INT m,l,lmax,ip,ip0;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	M       = solver->precon;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	m       = solver->options[LIS_OPTIONS_RESTART];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;

	r       = solver->work[0];
	rtld    = solver->work[1];
	p       = &solver->work[2];
	ap      = &solver->work[  (m+1)+2];
	aptld   = &solver->work[2*(m+1)+2];

	dotsave = (LIS_SCALAR *)lis_malloc( sizeof(LIS_SCALAR) * (m+1),"lis_orthomin::dotsave" );

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,M,r,rtld,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	
	iter=1;
	while( iter<=maxiter )
	{
		ip = (iter-1) % (m+1);

		/* p[ip] = rtld */
		lis_vector_copy(rtld,p[ip]);

		/* ap[ip]    = A*p[ip] */
		/* aptld[ip] = M^-1 ap[ip] */
		lis_matvec(A,p[ip],ap[ip]);
		time = lis_wtime();
		lis_psolve(solver, ap[ip], aptld[ip]);
		ptime += lis_wtime()-time;

		lmax = _min(m,iter-1);
		for(l=1;l<=lmax;l++)
		{
			ip0 = (ip+m+1-l) % (m+1);
			/* beta = -<Ar[ip],Ap[ip0]> / <Ap[ip0],Ap[ip0]> */
			lis_vector_dot(aptld[ip],aptld[ip0],&beta);
			beta = -beta * dotsave[l-1];

			lis_vector_axpy(beta,p[ip0]    ,p[ip]);
			lis_vector_axpy(beta,ap[ip0]   ,ap[ip]);
			lis_vector_axpy(beta,aptld[ip0],aptld[ip]);
		}
		for(l=m-1;l>0;l--)
		{
			dotsave[l] = dotsave[l-1];
		}

		lis_vector_dot(aptld[ip],aptld[ip],&dotsave[0]);
		/* test breakdown */
		if( dotsave[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			lis_free(dotsave);
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		dotsave[0] = 1.0/dotsave[0];

		/* alpha = <rtld,Aptld[ip]> */
		lis_vector_dot(rtld,aptld[ip],&alpha);
		alpha = alpha * dotsave[0];

		lis_vector_axpy( alpha,p[ip],x);
		lis_vector_axpy(-alpha,ap[ip],r);
		lis_vector_axpy(-alpha,aptld[ip],rtld);

		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);
		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
			if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
		}

		if( tol >= nrm2 )
		{
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			lis_free(dotsave);
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		iter++;
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	lis_free(dotsave);
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_orthomin_quad"
LIS_INT lis_orthomin_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_PRECON M;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, *p, *ap, *aptld;
	LIS_QUAD *dotsave;
	LIS_QUAD_PTR alpha, beta, tmp, one;

	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output,conv;
	double time,ptime;

	LIS_INT m,l,lmax,ip,ip0;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	M       = solver->precon;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	m       = solver->options[LIS_OPTIONS_RESTART];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(tmp,3,1);
	LIS_QUAD_SCALAR_MALLOC(one,4,1);

	r       = solver->work[0];
	rtld    = solver->work[1];
	p       = &solver->work[2];
	ap      = &solver->work[  (m+1)+2];
	aptld   = &solver->work[2*(m+1)+2];

	one.hi[0] = 1.0;
	one.lo[0] = 0.0;

	dotsave = (LIS_QUAD *)lis_malloc( sizeof(LIS_QUAD) * (m+1),"lis_orthomin_quad::dotsave" );

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,M,r,rtld,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	
	iter=1;
	while( iter<=maxiter )
	{
		ip = (iter-1) % (m+1);

		/* p[ip] = rtld */
		lis_vector_copyex_mm(rtld,p[ip]);

		/* ap[ip]    = A*p[ip] */
		/* aptld[ip] = M^-1 ap[ip] */
		lis_matvec(A,p[ip],ap[ip]);
		time = lis_wtime();
		lis_psolve(solver, ap[ip], aptld[ip]);
		ptime += lis_wtime()-time;

		lmax = _min(m,iter-1);
		for(l=1;l<=lmax;l++)
		{
			ip0 = (ip+m+1-l) % (m+1);
			/* beta = -<Ar[ip],Ap[ip0]> / <Ap[ip0],Ap[ip0]> */
			lis_vector_dotex_mmm(aptld[ip],aptld[ip0],&beta);
			lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,&dotsave[l-1]);
			lis_quad_minus((LIS_QUAD *)beta.hi);

			lis_vector_axpyex_mmm(beta,p[ip0]    ,p[ip]);
			lis_vector_axpyex_mmm(beta,ap[ip0]   ,ap[ip]);
			lis_vector_axpyex_mmm(beta,aptld[ip0],aptld[ip]);
		}
		for(l=m-1;l>0;l--)
		{
			dotsave[l] = dotsave[l-1];
		}

		lis_vector_dotex_mmm(aptld[ip],aptld[ip],&tmp);
		dotsave[0].hi = tmp.hi[0];
		dotsave[0].lo = tmp.lo[0];
		/* test breakdown */
		if( tmp.hi[0]==0.0 && tmp.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			lis_free(dotsave);
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		lis_quad_div(&dotsave[0],(LIS_QUAD *)one.hi,&dotsave[0]);

		/* alpha = <rtld,Aptld[ip]> */
		lis_vector_dotex_mmm(rtld,aptld[ip],&alpha);
		lis_quad_mul((LIS_QUAD *)alpha.hi,(LIS_QUAD *)alpha.hi,&dotsave[0]);

		lis_vector_axpyex_mmm( alpha,p[ip],x);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,ap[ip],r);
		lis_vector_axpyex_mmm(alpha,aptld[ip],rtld);
		lis_quad_minus((LIS_QUAD *)alpha.hi);

		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);
		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
			if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
		}

		if( tol > nrm2 )
		{
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			lis_free(dotsave);
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		iter++;
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	lis_free(dotsave);
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}
#endif

