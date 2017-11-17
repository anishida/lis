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

/**********************************************
 * Preconditioned Conjugate Gradient Squared  *
 **********************************************
 r(0)    = b - Ax(0)
 rtld(0) = conj(r(0)) or random
 rho(-1) = 1
 p(0)    = (0,...,0)^T
 q(0)    = (0,...,0)^T
 **********************************************
 for k=1,2,...
   rho(k-1)  = <rtld,r(k-1)> 
   beta      = rho(k-1) / rho(k-2) 
   u(k)      = r(k-1) + beta*q(k-1) 
   p(k)      = u(k) + beta*(q(k-1) + beta*p(k-1)) 
   phat(k)   = M^-1 * p(k) 
   vhat(k)   = A * phat(k) 
   tmpdot1   = <rtld,vhat(k-1)> 
   alpha     = rho(k-1) / tmpdot1 
   q(k)      = u(k) - alpha*vhat(k) 
   phat(k)   = u(k) + q(k)
   uhat(k)   = M^-1 * (u(k) + q(k))
   x(k)      = x(k-1) + alpha*uhat(k)
   qhat(k)   = A * uhat(k)
   r(k)      = r(k-1) - alpha*qhat(k)
 **********************************************/

#define NWORK 7
#undef __FUNC__
#define __FUNC__ "lis_cgs_check_params"
LIS_INT lis_cgs_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_cgs_malloc_work"
LIS_INT lis_cgs_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_cgs_malloc_work::work" );
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
#define __FUNC__ "lis_cgs"
LIS_INT lis_cgs(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, p,phat, q, qhat, u, uhat, vhat;
	LIS_SCALAR alpha, beta, rho, rho_old, tmpdot1;
	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output,conv;
	double time,ptime;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;

	r       = solver->work[0];
	rtld    = solver->work[1];
	p       = solver->work[2];
	phat    = solver->work[3];
	q       = solver->work[4];
	qhat    = solver->work[5];
	u       = solver->work[5];
	uhat    = solver->work[6];
	vhat    = solver->work[6];
	alpha   = (LIS_SCALAR)1.0;
	rho_old = (LIS_SCALAR)1.0;


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_all(0,q);
	lis_vector_set_all(0,p);

	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/* rho = <rtld,r> */
		lis_vector_dot(rtld,r,&rho);

		/* test breakdown */
		if( rho==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta = (rho / rho_old) */
		beta = rho / rho_old;

		/* u = r + beta*q */
		lis_vector_axpyz(beta,q,r,u);

		/* p = u + beta*(q + beta*p) */
		lis_vector_xpay(q,beta,p);
		lis_vector_xpay(u,beta,p);
		
		/* phat = M^-1 * p */
		time = lis_wtime();
		lis_psolve(solver, p, phat);
		ptime += lis_wtime()-time;

		/* v = A * phat */
		lis_matvec(A,phat,vhat);
		
		/* tmpdot1 = <rtld,vhat> */
		lis_vector_dot(rtld,vhat,&tmpdot1);
		/* test breakdown */
		if( tmpdot1==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		
		/* alpha = rho / tmpdot1 */
		alpha = rho / tmpdot1;
		
		/* q = u - alpha*vhat */
		lis_vector_axpyz(-alpha,vhat,u,q);

		/* phat = u + q          */
		/* uhat = M^-1 * (u + q) */
		lis_vector_axpyz(1,u,q,phat);
		time = lis_wtime();
		lis_psolve(solver, phat, uhat);
		ptime += lis_wtime()-time;

		/* x = x + alpha*uhat */
		lis_vector_axpy(alpha,uhat,x);

		/* qhat = A * uhat */
		lis_matvec(A,uhat,qhat);

		/* r = r - alpha*qhat */
		lis_vector_axpy(-alpha,qhat,r);

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
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}
		
		rho_old = rho;
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_cgs_quad"
LIS_INT lis_cgs_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, p,phat, q, qhat, u, uhat, vhat;
	LIS_QUAD_PTR alpha, beta, rho, rho_old, tmpdot1, one;
	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output,conv;
	double time,ptime;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;

	r       = solver->work[0];
	rtld    = solver->work[1];
	p       = solver->work[2];
	phat    = solver->work[3];
	q       = solver->work[4];
	qhat    = solver->work[5];
	u       = solver->work[5];
	uhat    = solver->work[6];
	vhat    = solver->work[6];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot1,4,1);
	LIS_QUAD_SCALAR_MALLOC(one,6,1);
	rho_old.hi[0] = 1.0;
	rho_old.lo[0] = 0.0;
	alpha.hi[0]   = 1.0;
	alpha.lo[0]   = 0.0;
	one.hi[0]   = 1.0;
	one.lo[0]   = 0.0;


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_allex_nm(0.0, q);
	lis_vector_set_allex_nm(0.0, p);


	for( iter=1; iter<=maxiter; iter++ )
	{
		/* rho = <rtld,r> */
		lis_vector_dotex_mmm(rtld,r,&rho);

		/* test breakdown */
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta = (rho / rho_old) */
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);

		/* u = r + beta*q */
		lis_vector_axpyzex_mmmm(beta,q,r,u);

		/* p = u + beta*(q + beta*p) */
		lis_vector_xpayex_mmm(q,beta,p);
		lis_vector_xpayex_mmm(u,beta,p);
		
		/* phat = M^-1 * p */
		time = lis_wtime();
		lis_psolve(solver, p, phat);
		ptime += lis_wtime()-time;

		/* v = A * phat */
		lis_matvec(A,phat,vhat);
		
		/* tmpdot1 = <rtld,vhat> */
		lis_vector_dotex_mmm(rtld,vhat,&tmpdot1);
		/* test breakdown */
		if( tmpdot1.hi[0]==0.0 && tmpdot1.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		
		/* alpha = rho / tmpdot1 */
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)tmpdot1.hi);
		
		/* q = u - alpha*vhat */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyzex_mmmm(alpha,vhat,u,q);

		/* phat = u + q          */
		/* uhat = M^-1 * (u + q) */
		lis_vector_axpyzex_mmmm(one,u,q,phat);
		time = lis_wtime();
		lis_psolve(solver, phat, uhat);
		ptime += lis_wtime()-time;

		/* x = x + alpha*uhat */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,uhat,x);

		/* qhat = A * uhat */
		lis_matvec(A,uhat,qhat);

		/* r = r - alpha*qhat */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,qhat,r);

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
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}
		
		rho_old.hi[0] = rho.hi[0];
		rho_old.lo[0] = rho.lo[0];
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}


#undef __FUNC__
#define __FUNC__ "lis_cgs_switch"
LIS_INT lis_cgs_switch(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, p,phat, q, qhat, u, uhat, vhat;
	LIS_QUAD_PTR alpha, beta, rho, rho_old, tmpdot1, one;
	LIS_REAL bnrm2, nrm2, tol, tol2;
	LIS_INT iter,maxiter,output,conv;
	LIS_INT iter2,maxiter2;
	double time,ptime;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	x       = solver->x;
	maxiter  = solver->options[LIS_OPTIONS_MAXITER];
	maxiter2 = solver->options[LIS_OPTIONS_SWITCH_MAXITER];
	output   = solver->options[LIS_OPTIONS_OUTPUT];
	conv     = solver->options[LIS_OPTIONS_CONV_COND];
	tol      = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	tol2     = solver->params[LIS_PARAMS_SWITCH_RESID-LIS_OPTIONS_LEN];
	ptime    = 0.0;

	r       = solver->work[0];
	rtld    = solver->work[1];
	p       = solver->work[2];
	phat    = solver->work[3];
	q       = solver->work[4];
	qhat    = solver->work[5];
	u       = solver->work[5];
	uhat    = solver->work[6];
	vhat    = solver->work[6];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot1,4,1);
	LIS_QUAD_SCALAR_MALLOC(one,6,1);
	rho_old.hi[0] = 1.0;
	rho_old.lo[0] = 0.0;
	alpha.hi[0]   = 1.0;
	alpha.lo[0]   = 0.0;
	one.hi[0]   = 1.0;
	one.lo[0]   = 0.0;

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol2     = solver->tol_switch;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_allex_nm(0.0, q);
	lis_vector_set_allex_nm(0.0, p);

	uhat->precision = LIS_PRECISION_DEFAULT;
	p->precision = LIS_PRECISION_DEFAULT;
	phat->precision = LIS_PRECISION_DEFAULT;

	for( iter=1; iter<=maxiter2; iter++ )
	{
			/* rho = <rtld,r> */
			lis_vector_dot(rtld,r,&rho.hi[0]);

			/* test breakdown */
			if( rho.hi[0]==0.0 )
			{
				solver->retcode   = LIS_BREAKDOWN;
				solver->iter      = iter;
				solver->iter2     = iter;
				solver->resid     = nrm2;
				LIS_DEBUG_FUNC_OUT;
				return LIS_BREAKDOWN;
			}

			/* beta = (rho / rho_old) */
			beta.hi[0] = (rho.hi[0] / rho_old.hi[0]);

			/* u = r + beta*q */
			lis_vector_axpyz(beta.hi[0],q,r,u);

			/* p = u + beta*(q + beta*p) */
			lis_vector_xpay(q,beta.hi[0],p);
			lis_vector_xpay(u,beta.hi[0],p);
			
			/* phat = M^-1 * p */
			time = lis_wtime();
			lis_psolve(solver, p, phat);
			ptime += lis_wtime()-time;

			/* v = A * phat */
			lis_matvec(A,phat,vhat);
			
			/* tmpdot1 = <rtld,vhat> */
			lis_vector_dot(rtld,vhat,&tmpdot1.hi[0]);
			/* test breakdown */
			if( tmpdot1.hi[0]==0.0 )
			{
				solver->retcode   = LIS_BREAKDOWN;
				solver->iter      = iter;
				solver->iter2     = iter;
				solver->resid     = nrm2;
				LIS_DEBUG_FUNC_OUT;
				return LIS_BREAKDOWN;
			}
			
			/* alpha = rho / tmpdot1 */
			alpha.hi[0] = rho.hi[0] / tmpdot1.hi[0];
			
			/* q = u - alpha*vhat */
			lis_vector_axpyz(-alpha.hi[0],vhat,u,q);

			/* phat = u + q          */
			/* uhat = M^-1 * (u + q) */
			lis_vector_axpyz(1.0,u,q,phat);
			time = lis_wtime();
			lis_psolve(solver, phat, uhat);
			ptime += lis_wtime()-time;

			/* x = x + alpha*uhat */
			lis_vector_axpy(alpha.hi[0],uhat,x);

			/* qhat = A * uhat */
			lis_matvec(A,uhat,qhat);

			/* r = r - alpha*qhat */
			lis_vector_axpy(-alpha.hi[0],qhat,r);

			/* convergence check */
			lis_solver_get_residual[conv](r,solver,&nrm2);
			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			if( nrm2 <= tol2 )
			{
				solver->iter       = iter;
				solver->iter2      = iter;
				solver->ptime      = ptime;
				break;
			}
			
			rho_old.hi[0] = rho.hi[0];
	}

	uhat->precision = LIS_PRECISION_QUAD;
	p->precision = LIS_PRECISION_QUAD;
	phat->precision = LIS_PRECISION_QUAD;

	solver->options[LIS_OPTIONS_INITGUESS_ZEROS] = LIS_FALSE;
	lis_vector_copyex_mn(x,solver->xx);
	rho_old.hi[0] = 1.0;

	lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2);
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_allex_nm(0.0, q);
	lis_vector_set_allex_nm(0.0, p);


	for( iter2=iter+1; iter2<=maxiter; iter2++ )
	{
			/* rho = <rtld,r> */
			lis_vector_dotex_mmm(rtld,r,&rho);

			/* test breakdown */
			if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
			{
				solver->retcode   = LIS_BREAKDOWN;
				solver->iter       = iter2;
				solver->iter2      = iter;
				solver->resid     = nrm2;
				LIS_DEBUG_FUNC_OUT;
				return LIS_BREAKDOWN;
			}

			/* beta = (rho / rho_old) */
			lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);

			/* u = r + beta*q */
			lis_vector_axpyzex_mmmm(beta,q,r,u);

			/* p = u + beta*(q + beta*p) */
			lis_vector_xpayex_mmm(q,beta,p);
			lis_vector_xpayex_mmm(u,beta,p);
			
			/* phat = M^-1 * p */
			time = lis_wtime();
			lis_psolve(solver, p, phat);
			ptime += lis_wtime()-time;

			/* v = A * phat */
			lis_matvec(A,phat,vhat);
			
			/* tmpdot1 = <rtld,vhat> */
			lis_vector_dotex_mmm(rtld,vhat,&tmpdot1);
			/* test breakdown */
			if( tmpdot1.hi[0]==0.0 && tmpdot1.lo[0]==0.0 )
			{
				solver->retcode   = LIS_BREAKDOWN;
				solver->iter       = iter2;
				solver->iter2      = iter;
				solver->resid     = nrm2;
				LIS_DEBUG_FUNC_OUT;
				return LIS_BREAKDOWN;
			}
			
			/* alpha = rho / tmpdot1 */
			lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)tmpdot1.hi);
			
			/* q = u - alpha*vhat */
			lis_quad_minus((LIS_QUAD *)alpha.hi);
			lis_vector_axpyzex_mmmm(alpha,vhat,u,q);

			/* phat = u + q          */
			/* uhat = M^-1 * (u + q) */
			lis_vector_axpyzex_mmmm(one,u,q,phat);
			time = lis_wtime();
			lis_psolve(solver, phat, uhat);
			ptime += lis_wtime()-time;

			/* x = x + alpha*uhat */
			lis_quad_minus((LIS_QUAD *)alpha.hi);
			lis_vector_axpyex_mmm(alpha,uhat,x);

			/* qhat = A * uhat */
			lis_matvec(A,uhat,qhat);

			/* r = r - alpha*qhat */
			lis_quad_minus((LIS_QUAD *)alpha.hi);
			lis_vector_axpyex_mmm(alpha,qhat,r);

			/* convergence check */
			lis_solver_get_residual[conv](r,solver,&nrm2);
			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter2] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}
			
			if( tol > nrm2 )
			{
				solver->retcode    = LIS_SUCCESS;
				solver->iter       = iter2;
				solver->iter2      = iter;
				solver->resid      = nrm2;
				solver->ptime      = ptime;
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}
			
			rho_old.hi[0] = rho.hi[0];
			rho_old.lo[0] = rho.lo[0];
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter2;
	solver->iter2     = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}
#endif

/**********************************************
 * Preconditioned Conjugate Residual Squared  *
 **********************************************
 r(0)    = b - Ax(0)
 rtld(0) = conj(r(0)) or random
 rtld(0) = A^H * rtld(0)
 rho(0)  = 1
 p(0)    = (0,...,0)^T
 q(0)    = (0,...,0)^T
 **********************************************
 for k=1,2,...
   z(k)      = M^-1 * r(k)
   rho(k)    = <rtld,z(k)> 
   beta      = rho(k) / rho(k-1)
   u(k)      = z(k) + beta*q(k-1) 
   p(k)      = u(k) + beta*(q(k-1) + beta*p(k-1)) 
   ap(k)     = A * p(k)
   map(k)    = M^-1 * ap(k)
   tmpdot1   = <rtld,map(k)>
   alpha     = rho(k-1) / tmpdot1 
   q(k)      = u(k) - alpha*map(k)
   uq(k)     = u(k) + q(k)
   auq(k)    = A * uq(k)
   x(k)      = x(k-1) + alpha*uq(k)
   r(k)      = r(k-1) - alpha*auq(k)
 **********************************************
   z = u = uq, ap = q, map = auq
 **********************************************/
#undef NWORK
#define NWORK 6
#undef __FUNC__
#define __FUNC__ "lis_crs_check_params"
LIS_INT lis_crs_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_crs_malloc_work"
LIS_INT lis_crs_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_cgs_malloc_work::work" );
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
#define __FUNC__ "lis_crs"
LIS_INT lis_crs(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, p, q, u, z, ap, map, uq, auq;
	LIS_SCALAR alpha, beta, rho, rho_old, tmpdot1;
	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output,conv;
	double time,ptime;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;

	r       = solver->work[0];
	rtld    = solver->work[1];
	p       = solver->work[2];
	z       = solver->work[3];
	u       = solver->work[3];
	uq      = solver->work[3];
	q       = solver->work[4];
	ap      = solver->work[4];
	map     = solver->work[5];
	auq     = solver->work[5];


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,p);

	lis_matvech(A,p,rtld);
	rho_old = 1.0;
	lis_vector_set_all(0,q);
	lis_vector_set_all(0,p);

	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/* z   = M^-1 * r  */
		/* rho = <rtld,z>  */
		time = lis_wtime();
		lis_psolve(solver, r, z);
		ptime += lis_wtime()-time;
		lis_vector_dot(rtld,z,&rho);

		/* test breakdown */
		if( rho==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta    = rho / rho_old         */
		/* u       = z + beta*q            */
		/* p       = u + beta*(q + beta*p) */
		/* ap      = A * p                 */
		/* map     = M^-1 * ap             */
		/* tmpdot1 = <rtld,map>            */
		beta = rho / rho_old;
		lis_vector_axpyz(beta,q,z,u);
		lis_vector_xpay(q,beta,p);
		lis_vector_xpay(u,beta,p);
		lis_matvec(A,p,ap);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;
		lis_vector_dot(rtld,map,&tmpdot1);
		/* test breakdown */
		if( tmpdot1==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		
		/* alpha = rho / tmpdot1 */
		/* q     = u - alpha*map */
		/* uq    = u + q         */
		/* auq   = A * uq        */
		/* x     = x + alpha*uq  */
		/* r     = r - alpha*auq */
		alpha = rho / tmpdot1;
		lis_vector_axpyz(-alpha,map,u,q);
		lis_vector_axpyz(1,u,q,uq);
		lis_matvec(A,uq,auq);
		lis_vector_axpy(alpha,uq,x);
		lis_vector_axpy(-alpha,auq,r);

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
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}
		
		rho_old = rho;
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_crs_quad"
LIS_INT lis_crs_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, p, q, u, z, ap, map, uq, auq;
	LIS_QUAD_PTR alpha, beta, rho, rho_old, tmpdot1, one;
	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output,conv;
	double time,ptime;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;

	r       = solver->work[0];
	rtld    = solver->work[1];
	p       = solver->work[2];
	z       = solver->work[3];
	u       = solver->work[3];
	uq      = solver->work[3];
	q       = solver->work[4];
	ap      = solver->work[4];
	map     = solver->work[5];
	auq     = solver->work[5];
	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot1,4,1);
	LIS_QUAD_SCALAR_MALLOC(one,6,1);

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,p);

	lis_matvech(A,p,rtld);
	lis_vector_set_allex_nm(0.0,q);
	lis_vector_set_allex_nm(0.0,p);
	rho_old.hi[0] = 1.0;
	rho_old.lo[0] = 0.0;
	one.hi[0]   = 1.0;
	one.lo[0]   = 0.0;

	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/* z   = M^-1 * r  */
		/* rho = <rtld,z>  */
		time = lis_wtime();
		lis_psolve(solver, r, z);
		ptime += lis_wtime()-time;
		lis_vector_dotex_mmm(rtld,z,&rho);

		/* test breakdown */
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta    = rho / rho_old         */
		/* u       = z + beta*q            */
		/* p       = u + beta*(q + beta*p) */
		/* ap      = A * p                 */
		/* map     = M^-1 * ap             */
		/* tmpdot1 = <rtld,map>            */
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
		lis_vector_axpyzex_mmmm(beta,q,z,u);
		lis_vector_xpayex_mmm(q,beta,p);
		lis_vector_xpayex_mmm(u,beta,p);
		lis_matvec(A,p,ap);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;
		lis_vector_dotex_mmm(rtld,map,&tmpdot1);
		/* test breakdown */
		if( tmpdot1.hi[0]==0.0 && tmpdot1.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		
		/* alpha = rho / tmpdot1 */
		/* q     = u - alpha*map */
		/* uq    = u + q         */
		/* auq   = A * uq        */
		/* x     = x + alpha*uq  */
		/* r     = r - alpha*auq */
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)tmpdot1.hi);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyzex_mmmm(alpha,map,u,q);
		lis_vector_axpyzex_mmmm(one,u,q,uq);
		lis_matvec(A,uq,auq);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,uq,x);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,auq,r);

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
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}
		
		rho_old.hi[0] = rho.hi[0];
		rho_old.lo[0] = rho.lo[0];
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}
#endif
