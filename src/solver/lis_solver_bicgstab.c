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

/***************************************************
 * Preconditioned BiConjugate Gradient STABilized  *
 ***************************************************
 r(0)    = b - Ax(0)
 rtld(0) = conj(r(0)) or random
 rho(-1) = 1
 p(0)    = (0,...,0)^T
 s(0)    = (0,...,0)^T
 phat(0) = (0,...,0)^T
 shat(0) = (0,...,0)^T
 ***************************************************
 for k=1,2,...
   rho = <rtld,r>
   beta = (rho / rho_old) * (alpha / omega)
   p = r + beta*(p - omega*v)
   phat = M^-1 * p
   v = A * phat
   tmpdot1 = <rtld,v>
   alpha = rho / tmpdot1
   s = r - alpha*v
   shat = M^-1 * s
   t = A * shat
   tmpdot1 = <t,s>
   tmpdot2 = <t,t>
   omega   = tmpdot1 / tmpdot2
   x = x + alpha*phat + omega*shat
   r = s - omega*t
   ***************************************************/

#define NWORK 7
#undef __FUNC__
#define __FUNC__ "lis_bicgstab_check_params"
LIS_INT lis_bicgstab_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_bicgstab_malloc_work"
LIS_INT lis_bicgstab_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_bicgstab_malloc_work::work" );
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
#define __FUNC__ "lis_bicgstab"
LIS_INT lis_bicgstab(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, t,p,v, s, phat, shat;
	LIS_SCALAR alpha, beta, omega, rho, rho_old, tmpdot1, tmpdot2;
	LIS_REAL   bnrm2, nrm2, tol;
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

	rtld    = solver->work[0];
	r       = solver->work[1];
	s       = solver->work[1];
	t       = solver->work[2];
	p       = solver->work[3];
	v       = solver->work[4];
	phat    = solver->work[5];
	shat    = solver->work[6];
	alpha   = (LIS_SCALAR)1.0;
	omega   = (LIS_SCALAR)1.0;
	rho_old = (LIS_SCALAR)1.0;

	lis_vector_set_all(0.0,p);
	lis_vector_set_all(0.0,phat);
	lis_vector_set_all(0.0,s);
	lis_vector_set_all(0.0,shat);

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	
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

		if( iter==1 )
		{
			lis_vector_copy(r,p);
		}
		else
		{
			/* beta = (rho / rho_old) * (alpha / omega) */
			beta = (rho / rho_old) * (alpha / omega);
	
			/* p = r + beta*(p - omega*v) */
			lis_vector_axpy(-omega,v,p);
			lis_vector_xpay(r,beta,p);
		}
		
		/* phat = M^-1 * p */
		time = lis_wtime();
		lis_psolve(solver, p, phat);
		ptime += lis_wtime()-time;

		/* v = A * phat */
		lis_matvec(A,phat,v);

		/* tmpdot1 = <rtld,v> */
		lis_vector_dot(rtld,v,&tmpdot1);
		/* test breakdown */
		/* */
		
		/* alpha = rho / tmpdot1 */
		alpha = rho / tmpdot1;
		
		/* s = r - alpha*v */
		lis_vector_axpy(-alpha,v,r);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](s,solver,&nrm2);
/*		lis_vector_nrm2(s,&nrm2);
		nrm2 = nrm2 * bnrm2;*/
		if( nrm2 <= tol )
		{
			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			lis_vector_axpy(alpha,phat,x);
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		/* shat = M^-1 * s */
		time = lis_wtime();
		lis_psolve(solver, s, shat);
		ptime += lis_wtime()-time;

		/* t = A * shat */
		lis_matvec(A,shat,t);

		/* tmpdot1 = <t,s> */
		/* tmpdot2 = <t,t> */
		/* omega   = tmpdot1 / tmpdot2 */
		lis_vector_dot(t,s,&tmpdot1);
		lis_vector_dot(t,t,&tmpdot2);
		omega   = tmpdot1 / tmpdot2;

		/* x = x + alpha*phat + omega*shat */
		lis_vector_axpy(alpha,phat,x);
		lis_vector_axpy(omega,shat,x);
		
		/* r = s - omega*t */
		lis_vector_axpy(-omega,t,r);
		
		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);
/*		lis_vector_nrm2(r,&nrm2);
		nrm2 = nrm2 * bnrm2;*/

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
		
		if( omega==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
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
#define __FUNC__ "lis_bicgstab_quad"
LIS_INT lis_bicgstab_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, t,p,v, s, phat, shat;
	LIS_QUAD_PTR alpha, beta, omega, rho, rho_old, tmpdot1, tmpdot2;
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

	rtld    = solver->work[0];
	r       = solver->work[1];
	s       = solver->work[1];
	t       = solver->work[2];
	p       = solver->work[3];
	v       = solver->work[4];
	phat    = solver->work[5];
	shat    = solver->work[6];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot1,4,1);
	LIS_QUAD_SCALAR_MALLOC(omega,6,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot2,7,1);
	rho_old.hi[0] = 1.0;
	rho_old.lo[0] = 0.0;
	alpha.hi[0] = 1.0;
	alpha.lo[0] = 0.0;
	omega.hi[0] = 1.0;
	omega.lo[0] = 0.0;

	lis_vector_set_allex_nm(0.0, p);
	lis_vector_set_allex_nm(0.0, phat);
	lis_vector_set_allex_nm(0.0, s);
	lis_vector_set_allex_nm(0.0, shat);

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	
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

		if( iter==1 )
		{
			lis_vector_copyex_mm(r,p);
		}
		else
		{
			/* beta = (rho / rho_old) * (alpha / omega) */
			lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
			lis_quad_div((LIS_QUAD *)tmpdot1.hi,(LIS_QUAD *)alpha.hi,(LIS_QUAD *)omega.hi);
			lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,(LIS_QUAD *)tmpdot1.hi);
	
			/* p = r + beta*(p - omega*v) */
			lis_quad_minus((LIS_QUAD *)omega.hi);
			lis_vector_axpyex_mmm(omega,v,p);
			lis_vector_xpayex_mmm(r,beta,p);
		}
		
		/* phat = M^-1 * p */
		time = lis_wtime();
		lis_psolve(solver, p, phat);
		ptime += lis_wtime()-time;

		/* v = A * phat */
		lis_matvec(A,phat,v);

		/* tmpdot1 = <rtld,v> */
		lis_vector_dotex_mmm(rtld,v,&tmpdot1);
		/* test breakdown */
		/* */
		
		/* alpha = rho / tmpdot1 */
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)tmpdot1.hi);
		
		/* s = r - alpha*v */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,v,r);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](s,solver,&nrm2);
		if( tol > nrm2 )
		{
			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			lis_quad_minus((LIS_QUAD *)alpha.hi);
			lis_vector_axpyex_mmm(alpha,phat,x);
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		/* shat = M^-1 * s */
		time = lis_wtime();
		lis_psolve(solver, s, shat);
		ptime += lis_wtime()-time;

		/* t = A * shat */
		lis_matvec(A,shat,t);

		/* tmpdot1 = <t,s> */
		/* tmpdot2 = <t,t> */
		/* omega   = tmpdot1 / tmpdot2 */
		lis_vector_dotex_mmm(t,s,&tmpdot1);
		lis_vector_dotex_mmm(t,t,&tmpdot2);
		lis_quad_div((LIS_QUAD *)omega.hi,(LIS_QUAD *)tmpdot1.hi,(LIS_QUAD *)tmpdot2.hi);

		/* x = x + alpha*phat + omega*shat */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,phat,x);
		lis_vector_axpyex_mmm(omega,shat,x);
		
		/* r = s - omega*t */
		lis_quad_minus((LIS_QUAD *)omega.hi);
		lis_vector_axpyex_mmm(omega,t,r);
		lis_quad_minus((LIS_QUAD *)omega.hi);
		
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
		
		if( omega.hi[0]==0.0 && omega.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
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
#define __FUNC__ "lis_bicgstab_switch"
LIS_INT lis_bicgstab_switch(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, t,p,v, s, phat, shat;
	LIS_QUAD_PTR alpha, beta, omega, rho, rho_old, tmpdot1, tmpdot2;
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

	rtld    = solver->work[0];
	r       = solver->work[1];
	s       = solver->work[1];
	t       = solver->work[2];
	p       = solver->work[3];
	v       = solver->work[4];
	phat    = solver->work[5];
	shat    = solver->work[6];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot1,4,1);
	LIS_QUAD_SCALAR_MALLOC(omega,6,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot2,7,1);
	rho_old.hi[0] = 1.0;
	rho_old.lo[0] = 0.0;
	alpha.hi[0] = 1.0;
	alpha.lo[0] = 0.0;
	omega.hi[0] = 1.0;
	omega.lo[0] = 0.0;

	lis_vector_set_allex_nm(0.0, p);
	lis_vector_set_allex_nm(0.0, phat);
	lis_vector_set_allex_nm(0.0, s);
	lis_vector_set_allex_nm(0.0, shat);

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol2     = solver->tol_switch;

	lis_solver_set_shadowresidual(solver,r,rtld);

	s->precision = LIS_PRECISION_DEFAULT;
	shat->precision = LIS_PRECISION_DEFAULT;
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

			if( iter==1 )
			{
				lis_vector_copy(r,p);
			}
			else
			{
				/* beta = (rho / rho_old) * (alpha / omega) */
				beta.hi[0] = (rho.hi[0] / rho_old.hi[0]) * (alpha.hi[0] / omega.hi[0]);
		
				/* p = r + beta*(p - omega*v) */
				lis_vector_axpy(-omega.hi[0],v,p);
				lis_vector_xpay(r,beta.hi[0],p);
			}
			
			/* phat = M^-1 * p */
			time = lis_wtime();
			lis_psolve(solver, p, phat);
			ptime += lis_wtime()-time;

			/* v = A * phat */
			lis_matvec(A,phat,v);

			/* tmpdot1 = <rtld,v> */
			lis_vector_dot(rtld,v,&tmpdot1.hi[0]);
			/* test breakdown */
			/* */
			
			/* alpha = rho / tmpdot1 */
			alpha.hi[0] = rho.hi[0] / tmpdot1.hi[0];
			
			/* s = r - alpha*v */
			lis_vector_axpy(-alpha.hi[0],v,r);

			/* Early check for tolerance */
			lis_solver_get_residual[conv](s,solver,&nrm2);
			if( nrm2 <= tol2 )
			{
				if( output )
				{
					if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
					if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
				}

				lis_vector_axpy(alpha.hi[0],phat,x);
				solver->iter       = iter;
				solver->iter2      = iter;
				solver->ptime      = ptime;
				break;
			}

			/* shat = M^-1 * s */
			time = lis_wtime();
			lis_psolve(solver, s, shat);
			ptime += lis_wtime()-time;

			/* t = A * shat */
			lis_matvec(A,shat,t);

			/* tmpdot1 = <t,s> */
			/* tmpdot2 = <t,t> */
			/* omega   = tmpdot1 / tmpdot2 */
			lis_vector_dot(t,s,&tmpdot1.hi[0]);
			lis_vector_dot(t,t,&tmpdot2.hi[0]);
			omega.hi[0]   = tmpdot1.hi[0] / tmpdot2.hi[0];

			/* x = x + alpha*phat + omega*shat */
			lis_vector_axpy(alpha.hi[0],phat,x);
			lis_vector_axpy(omega.hi[0],shat,x);
			
			/* r = s - omega*t */
			lis_vector_axpy(-omega.hi[0],t,r);
			
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
			
			if( omega.hi[0]==0.0 )
			{
				solver->retcode   = LIS_BREAKDOWN;
				solver->iter      = iter;
				solver->iter2     = iter;
				solver->resid     = nrm2;
				LIS_DEBUG_FUNC_OUT;
				return LIS_BREAKDOWN;
			}
			rho_old.hi[0] = rho.hi[0];
	}

	s->precision = LIS_PRECISION_QUAD;
	shat->precision = LIS_PRECISION_QUAD;
	p->precision = LIS_PRECISION_QUAD;
	phat->precision = LIS_PRECISION_QUAD;

	solver->options[LIS_OPTIONS_INITGUESS_ZEROS] = LIS_FALSE;
	lis_vector_copyex_mn(x,solver->xx);
	rho_old.hi[0] = 1.0;
	alpha.hi[0] = 1.0;
	omega.hi[0] = 1.0;

	lis_vector_set_allex_nm(0.0, p);
	lis_vector_set_allex_nm(0.0, phat);
	lis_vector_set_allex_nm(0.0, s);
	lis_vector_set_allex_nm(0.0, shat);

	/* Initial Residual */
	lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2);
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

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

			if( iter2==1 )
			{
				lis_vector_copyex_mm(r,p);
			}
			else
			{
				/* beta = (rho / rho_old) * (alpha / omega) */
				lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
				lis_quad_div((LIS_QUAD *)tmpdot1.hi,(LIS_QUAD *)alpha.hi,(LIS_QUAD *)omega.hi);
				lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,(LIS_QUAD *)tmpdot1.hi);
		
				/* p = r + beta*(p - omega*v) */
				lis_quad_minus((LIS_QUAD *)omega.hi);
				lis_vector_axpyex_mmm(omega,v,p);
				lis_vector_xpayex_mmm(r,beta,p);
			}
			
			/* phat = M^-1 * p */
			time = lis_wtime();
			lis_psolve(solver, p, phat);
			ptime += lis_wtime()-time;

			/* v = A * phat */
			lis_matvec(A,phat,v);

			/* tmpdot1 = <rtld,v> */
			lis_vector_dotex_mmm(rtld,v,&tmpdot1);
			/* test breakdown */
			/* */
			
			/* alpha = rho / tmpdot1 */
			lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)tmpdot1.hi);
			
			/* s = r - alpha*v */
			lis_quad_minus((LIS_QUAD *)alpha.hi);
			lis_vector_axpyex_mmm(alpha,v,r);

			/* Early check for tolerance */
			lis_solver_get_residual[conv](s,solver,&nrm2);
			if( tol > nrm2 )
			{
				if( output )
				{
					if( output & LIS_PRINT_MEM ) solver->rhistory[iter2] = nrm2;
					if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
				}

				lis_quad_minus((LIS_QUAD *)alpha.hi);
				lis_vector_axpyex_mmm(alpha,phat,x);
				solver->retcode    = LIS_SUCCESS;
				solver->iter       = iter2;
				solver->iter2      = iter;
				solver->resid      = nrm2;
				solver->ptime      = ptime;
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}

			/* shat = M^-1 * s */
			time = lis_wtime();
			lis_psolve(solver, s, shat);
			ptime += lis_wtime()-time;

			/* t = A * shat */
			lis_matvec(A,shat,t);

			/* tmpdot1 = <t,s> */
			/* tmpdot2 = <t,t> */
			/* omega   = tmpdot1 / tmpdot2 */
			lis_vector_dotex_mmm(t,s,&tmpdot1);
			lis_vector_dotex_mmm(t,t,&tmpdot2);
			lis_quad_div((LIS_QUAD *)omega.hi,(LIS_QUAD *)tmpdot1.hi,(LIS_QUAD *)tmpdot2.hi);

			/* x = x + alpha*phat + omega*shat */
			lis_quad_minus((LIS_QUAD *)alpha.hi);
			lis_vector_axpyex_mmm(alpha,phat,x);
			lis_vector_axpyex_mmm(omega,shat,x);
			
			/* r = s - omega*t */
			lis_quad_minus((LIS_QUAD *)omega.hi);
			lis_vector_axpyex_mmm(omega,t,r);
			lis_quad_minus((LIS_QUAD *)omega.hi);
			
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
			
			if( omega.hi[0]==0.0 && omega.lo[0]==0.0 )
			{
				solver->retcode   = LIS_BREAKDOWN;
				solver->iter       = iter2;
				solver->iter2      = iter;
				solver->resid     = nrm2;
				LIS_DEBUG_FUNC_OUT;
				return LIS_BREAKDOWN;
			}
			rho_old.hi[0] = rho.hi[0];
			rho_old.lo[0] = rho.lo[0];
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter       = iter2;
	solver->iter2      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}
#endif

/***************************************************
 * Preconditioned BiConjugate Residual STABilized  *
 ***************************************************
 r(0)    = b - Ax(0)
 rtld(0) = conj(r(0)) or random
 rtld(0) = A^H * rtld(0)
 z(0)    = M^-1 * r(0)
 p(0)    = z(0)
 rho(0)  = <rtld(0),z(0)>
 ***************************************************
 for k=1,2,...
   ap(k-1)  = A * p(k-1)
   map(k-1) = M^-1 * ap(k-1)
   tmpdot1  = <rtld(0),map(k-1)>
   alpha    = rho(k-1) / tmpdot1
   s(k-1)   = r(k-1) - alpha*ap(k-1)
   ms(k-1)  = z(k-1) - alpha*map(k-1)
   ams(k-1) = A * ms(k-1)
   tmpdot1  = <ams(k-1),s(k-1)>
   tmpdot2  = <ams(k-1),ams(k-1)>
   omega    = tmpdot1 / tmpdot2
   x(k)     = x(k-1) + alpha*p(k-1) + omega*ms(k-1)
   r(k)     = s(k-1) - omega*ams(k-1)
   z(k)     = M^-1 * r(k)
   rho(k)   = <rtld(0),z(k)>
   beta     = (rho(k) / rho(k-1)) * (alpha / omega)
   p(k)     = z(k) + beta*(p(k-1) - omega*map(k-1))
   ***************************************************/
#undef NWORK
#define NWORK 9
#undef __FUNC__
#define __FUNC__ "lis_bicrstab_check_params"
LIS_INT lis_bicrstab_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_bicrstab_malloc_work"
LIS_INT lis_bicrstab_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_bicgstab_malloc_work::work" );
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
#define __FUNC__ "lis_bicrstab"
LIS_INT lis_bicrstab(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, p, s, ap, ms, map, ams, z;
	LIS_SCALAR alpha, beta, omega, rho, rho_old, tmpdot1, tmpdot2;
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

	rtld    = solver->work[0];
	r       = solver->work[1];
	s       = solver->work[2];
	ms      = solver->work[3];
	ams     = solver->work[4];
	p       = solver->work[5];
	ap      = solver->work[6];
	map     = solver->work[7];
	z       = solver->work[8];

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,p);

	lis_matvech(A,p,rtld);
	time = lis_wtime();
	lis_psolve(solver, r, z);
	ptime += lis_wtime()-time;
	lis_vector_copy(z,p);
	lis_vector_dot(rtld,z,&rho_old);
	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/* ap      = A * p             */
		/* map     = M^-1 * ap         */
		/* tmpdot1 = <rtld,map>        */
		/* alpha   = rho_old / tmpdot1 */
		/* s       = r - alpha*ap      */
		lis_matvec(A,p,ap);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;
		lis_vector_dot(rtld,map,&tmpdot1);
		alpha = rho_old / tmpdot1;
		lis_vector_axpyz(-alpha,ap,r,s);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](s,solver,&nrm2);
		if( nrm2 <= tol )
		{
			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			lis_vector_axpy(alpha,p,x);
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		/* ms      = z - alpha*map     */
		/* ams     = A * ms            */
		/* tmpdot1 = <ams,s>           */
		/* tmpdot2 = <ams,ams>         */
		/* omega   = tmpdot1 / tmpdot2 */
		lis_vector_axpyz(-alpha,map,z,ms);
		lis_matvec(A,ms,ams);
		lis_vector_dot(ams,s,&tmpdot1);
		lis_vector_dot(ams,ams,&tmpdot2);
		omega   = tmpdot1 / tmpdot2;

		/* x = x + alpha*p  + omega*ms  */
		/* r = s - omega*ams            */
		lis_vector_axpy(alpha,p,x);
		lis_vector_axpy(omega,ms,x);
		lis_vector_axpyz(-omega,ams,s,r);
		
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
		
		/* z   = M^-1 * r */
		/* rho = <rtld,z> */
		time = lis_wtime();
		lis_psolve(solver, r, z);
		ptime += lis_wtime()-time;
		lis_vector_dot(rtld,z,&rho);
		if( rho==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta = (rho / rho_old) * (alpha / omega) */
		/* p    = z + beta*(p - omega*map)          */
		beta = (rho / rho_old) * (alpha / omega);
		lis_vector_axpy(-omega,map,p);
		lis_vector_xpay(z,beta,p);

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
#define __FUNC__ "lis_bicrstab_quad"
LIS_INT lis_bicrstab_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, p, s, ap, ms, map, ams, z;
	LIS_QUAD_PTR alpha, beta, omega, rho, rho_old, tmpdot1, tmpdot2;
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

	rtld    = solver->work[0];
	r       = solver->work[1];
	s       = solver->work[2];
	ms      = solver->work[3];
	ams     = solver->work[4];
	p       = solver->work[5];
	ap      = solver->work[6];
	map     = solver->work[7];
	z       = solver->work[8];
	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot1,4,1);
	LIS_QUAD_SCALAR_MALLOC(omega,6,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot2,7,1);

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,p);

	lis_matvech(A,p,rtld);
	time = lis_wtime();
	lis_psolve(solver, r, z);
	ptime += lis_wtime()-time;
	lis_vector_copyex_mm(z,p);
	lis_vector_dotex_mmm(rtld,z,&rho_old);
	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/* ap      = A * p             */
		/* map     = M^-1 * ap         */
		/* tmpdot1 = <rtld,map>        */
		/* alpha   = rho_old / tmpdot1 */
		/* s       = r - alpha*ap      */
		lis_matvec(A,p,ap);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;
		lis_vector_dotex_mmm(rtld,map,&tmpdot1);
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho_old.hi,(LIS_QUAD *)tmpdot1.hi);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyzex_mmmm(alpha,ap,r,s);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](s,solver,&nrm2);
		if( nrm2 <= tol )
		{
			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			lis_quad_minus((LIS_QUAD *)alpha.hi);
			lis_vector_axpyex_mmm(alpha,p,x);
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		/* ms      = z - alpha*map     */
		/* ams     = A * ms            */
		/* tmpdot1 = <ams,s>           */
		/* tmpdot2 = <ams,ams>         */
		/* omega   = tmpdot1 / tmpdot2 */
		lis_vector_axpyzex_mmmm(alpha,map,z,ms);
		lis_matvec(A,ms,ams);
		lis_vector_dotex_mmm(ams,s,&tmpdot1);
		lis_vector_dotex_mmm(ams,ams,&tmpdot2);
		lis_quad_div((LIS_QUAD *)omega.hi,(LIS_QUAD *)tmpdot1.hi,(LIS_QUAD *)tmpdot2.hi);

		/* x = x + alpha*p  + omega*ms  */
		/* r = s - omega*ams            */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,p,x);
		lis_vector_axpyex_mmm(omega,ms,x);
		lis_quad_minus((LIS_QUAD *)omega.hi);
		lis_vector_axpyzex_mmmm(omega,ams,s,r);
		
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
		
		/* z   = M^-1 * r */
		/* rho = <rtld,z> */
		time = lis_wtime();
		lis_psolve(solver, r, z);
		ptime += lis_wtime()-time;
		lis_vector_dotex_mmm(rtld,z,&rho);
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta = (rho / rho_old) * (alpha / omega) */
		/* p    = z + beta*(p - omega*map)          */
		lis_quad_minus((LIS_QUAD *)omega.hi);
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
		lis_quad_div((LIS_QUAD *)tmpdot1.hi,(LIS_QUAD *)alpha.hi,(LIS_QUAD *)omega.hi);
		lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,(LIS_QUAD *)tmpdot1.hi);
		lis_quad_minus((LIS_QUAD *)omega.hi);
		lis_vector_axpyex_mmm(omega,map,p);
		lis_vector_xpayex_mmm(z,beta,p);

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
