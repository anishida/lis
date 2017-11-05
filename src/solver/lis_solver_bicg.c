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

/*****************************************
 * Preconditioned BiConjugate Gradient   *
 *****************************************
 r(0)    = b - Ax(0)
 rtld(0) = conj(r(0)) or random
 rho(-1) = 1
 p(0)    = (0,...,0)^T
 ptld(0) = (0,...,0)^T
 *****************************************
 for k=1,2,...
   z(k-1)    = M^-1 * r(k-1)
   ztld(k-1) = M^-H * rtld(k-1)
   rho(k-1)  = <rtld(k-1),z(k-1)>
   beta      = rho(k-1) / rho(k-2)
   p(k)      = z(k-1) + beta*p(k-1)
   ptld(k)   = ztld(k-1) + conj(beta)*ptld(k-1)
   q(k)      = A * p(k)
   qtld(k)   = A^H * ptld(k)
   tmpdot1   = <ptld(k),q(k)>
   alpha     = rho(k-1) / tmpdot1
   x(k)      = x(k-1) + alpha*p(k)
   r(k)      = r(k-1) - alpha*q(k)
   rtld(k)   = rtld(k-1) - conj(alpha)*qtld(k)
 *****************************************/

#define NWORK 6
#undef __FUNC__
#define __FUNC__ "lis_bicg_check_params"
LIS_INT lis_bicg_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_bicg_malloc_work"
LIS_INT lis_bicg_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	/*
	err = lis_matrix_convert(solver->A,&solver->Ah,LIS_MATRIX_CCS);
	if( err ) return err;
	*/

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_bicg_malloc_work::work" );
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
#define __FUNC__ "lis_bicg"
LIS_INT lis_bicg(LIS_SOLVER solver)
{
	LIS_Comm comm;
	LIS_MATRIX A,Ah;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, z,ztld,p, ptld, q, qtld;
	LIS_SCALAR alpha, beta, rho, rho_old, tmpdot1;
	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output,conv;
	double time,ptime;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	Ah      = solver->A;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;

	r       = solver->work[0];
	rtld    = solver->work[1];
	z       = solver->work[2];
	ztld    = solver->work[3];
	p       = solver->work[4];
	ptld    = solver->work[5];
	q       = solver->work[2];
	qtld    = solver->work[3];
	rho_old = (LIS_SCALAR)1.0;



	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);
	
	lis_vector_set_all(0, p);
	lis_vector_set_all(0, ptld);

	for( iter=1; iter<=maxiter; iter++ )
	{
		/* z    = M^-1 * r */
		/* ztld = M^-H * rtld */
		time = lis_wtime();
		lis_psolve(solver, r, z);
		lis_psolveh(solver, rtld, ztld);
		ptime += lis_wtime()-time;

		/* rho = <rtld,z> */
		lis_vector_dot(rtld,z,&rho);
/*		printf("rho = %e\n",rho);*/

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

		/* p    = z    + beta*p    */
		/* ptld = ztld + conj(beta)*ptld */
		
		/* q    = A   * p    */
		/* qtld = A^H * ptld */
		lis_vector_xpay(z,beta,p);
		lis_matvec(A,p,q);

		lis_vector_xpay(ztld,conj(beta),ptld);
		lis_matvech(Ah,ptld,qtld);

		
		/* tmpdot1 = <ptld,q> */
		lis_vector_dot(ptld,q,&tmpdot1);
/*		printf("tmpdot1 = %e\n",tmpdot1);*/

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
		
		/* x = x + alpha*p */
		lis_vector_axpy(alpha,p,x);
		
		/* r    = r    - alpha*q    */
		lis_vector_axpy(-alpha,q,r);
		
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
		
		/* rtld = rtld - conj(alpha)*qtld */
		lis_vector_axpy(-conj(alpha),qtld,rtld);

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
#define __FUNC__ "lis_bicg_quad"
LIS_INT lis_bicg_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A,Ah;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, z,ztld,p, ptld, q, qtld;
	LIS_QUAD_PTR alpha, beta, rho, rho_old, tmpdot1;
	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output,conv;
	double time,ptime;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	Ah      = solver->A;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;

	r       = solver->work[0];
	rtld    = solver->work[1];
	z       = solver->work[2];
	ztld    = solver->work[3];
	p       = solver->work[4];
	ptld    = solver->work[5];
	q       = solver->work[2];
	qtld    = solver->work[3];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot1,4,1);
	rho_old.hi[0] = 1.0;
	rho_old.lo[0] = 0.0;

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_allex_nm(0.0, p);
	lis_vector_set_allex_nm(0.0, ptld);

	for( iter=1; iter<=maxiter; iter++ )
	{
		/* z    = M^-1 * r */
		/* ztld = M^-H * rtld */
		time = lis_wtime();
		lis_psolve(solver, r, z);
		lis_psolveh(solver, rtld, ztld);
		ptime += lis_wtime()-time;

		/* rho = <z,rtld> */
		lis_vector_dotex_mmm(z,rtld,&rho);
/*		printf("rho = %e %e\n",rho.hi[0],rho.lo[0]);*/

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
/*		printf("beta = %e %e\n",beta.hi[0],beta.lo[0]);*/

		/* p    = z    + beta*p    */
		/* ptld = ztld + beta*ptld */
		
		/* q    = A   * p    */
		/* qtld = A^H * ptld */
		lis_vector_xpayex_mmm(z,beta,p);
		lis_matvec(A,p,q);

		lis_vector_xpayex_mmm(ztld,beta,ptld);
		lis_matvech(Ah,ptld,qtld);

		
		/* tmpdot1 = <ptld,q> */
		lis_vector_dotex_mmm(ptld,q,&tmpdot1);
/*		printf("tmpdot1 = %e %e\n",tmpdot1.hi[0],tmpdot1.lo[0]);*/

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
/*		printf("alpha = %e %e\n",alpha.hi[0],alpha.lo[0]);*/
		
		/* x = x + alpha*p */
		lis_vector_axpyex_mmm(alpha,p,x);
		
		/* r    = r    - alpha*q    */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,q,r);

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
		
		/* rtld = rtld - alpha*qtld */
		lis_vector_axpyex_mmm(alpha,qtld,rtld);

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
#define __FUNC__ "lis_bicg_switch"
LIS_INT lis_bicg_switch(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A,Ah;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, z,ztld,p, ptld, q, qtld;
	LIS_QUAD_PTR alpha, beta, rho, rho_old, tmpdot1;
	LIS_REAL bnrm2, nrm2, tol, tol2;
	LIS_INT iter,maxiter,output,conv;
	LIS_INT iter2,maxiter2;
	double time,ptime;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;
	
	A       = solver->A;
	Ah      = solver->A;
	x        = solver->x;
	maxiter  = solver->options[LIS_OPTIONS_MAXITER];
	maxiter2 = solver->options[LIS_OPTIONS_SWITCH_MAXITER];
	output   = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	tol      = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	tol2     = solver->params[LIS_PARAMS_SWITCH_RESID-LIS_OPTIONS_LEN];
	ptime    = 0.0;

	r        = solver->work[0];
	rtld     = solver->work[1];
	z        = solver->work[2];
	ztld     = solver->work[3];
	p        = solver->work[4];
	ptld     = solver->work[5];
	q        = solver->work[2];
	qtld     = solver->work[3];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot1,4,1);
	rho_old.hi[0] = 1.0;
	rho_old.lo[0] = 0.0;


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol2     = solver->tol_switch;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_allex_nm(0.0, p);
	lis_vector_set_allex_nm(0.0, ptld);

	r->precision = LIS_PRECISION_DEFAULT;
	rtld->precision = LIS_PRECISION_DEFAULT;
	p->precision = LIS_PRECISION_DEFAULT;
	ptld->precision = LIS_PRECISION_DEFAULT;

	for( iter=1; iter<=maxiter2; iter++ )
	{
		/* z    = M^-1 * r */
		/* ztld = M^-H * rtld */
		time = lis_wtime();
		lis_psolve(solver, r, z);
		lis_psolveh(solver, rtld, ztld);
		ptime += lis_wtime()-time;

		/* rho = <z,rtld> */
		lis_vector_dot(z,rtld,&rho.hi[0]);

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
		beta.hi[0] = rho.hi[0] / rho_old.hi[0];

		/* p    = z    + beta*p    */
		/* ptld = ztld + beta*ptld */
		
		/* q    = A   * p    */
		/* qtld = A^H * ptld */
		lis_vector_xpay(z,beta.hi[0],p);
		lis_matvec(A,p,q);

		lis_vector_xpay(ztld,beta.hi[0],ptld);
		lis_matvech(Ah,ptld,qtld);

		
		/* tmpdot1 = <ptld,q> */
		lis_vector_dot(ptld,q,&tmpdot1.hi[0]);

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
		
		/* x = x + alpha*p */
		lis_vector_axpy(alpha.hi[0],p,x);
		
		/* r    = r    - alpha*q    */
		lis_vector_axpy(-alpha.hi[0],q,r);
		
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
		
		/* rtld = rtld - alpha*qtld */
		lis_vector_axpy(-alpha.hi[0],qtld,rtld);

		rho_old.hi[0] = rho.hi[0];
	}

	r->precision = LIS_PRECISION_QUAD;
	rtld->precision = LIS_PRECISION_QUAD;
	p->precision = LIS_PRECISION_QUAD;
	ptld->precision = LIS_PRECISION_QUAD;

/*	solver->precon->precon_type = 0;*/
	solver->options[LIS_OPTIONS_INITGUESS_ZEROS] = LIS_FALSE;
	lis_vector_copyex_mn(x,solver->xx);
	rho_old.hi[0] = 1.0;
	lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2);
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_allex_nm(0.0, p);
	lis_vector_set_allex_nm(0.0, ptld);

	for( iter2=iter+1; iter2<=maxiter; iter2++ )
	{
		/* z    = M^-1 * r */
		/* ztld = M^-H * rtld */
		time = lis_wtime();
		lis_psolve(solver, r, z);
		lis_psolveh(solver, rtld, ztld);
/*		memset(z->value_lo,0,n*sizeof(LIS_SCALAR));
		memset(ztld->value_lo,0,n*sizeof(LIS_SCALAR));*/
		ptime += lis_wtime()-time;

		/* rho = <z,rtld> */
		lis_vector_dotex_mmm(z,rtld,&rho);
/*		printf("rho = %e %e\n",rho.hi[0],rho.lo[0]);*/

		/* test breakdown */
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter2;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta = (rho / rho_old) */
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);

		/* p    = z    + beta*p    */
		/* ptld = ztld + beta*ptld */
		
		/* q    = A   * p    */
		/* qtld = A^H * ptld */
		lis_vector_xpayex_mmm(z,beta,p);
		lis_matvec(A,p,q);

		lis_vector_xpayex_mmm(ztld,beta,ptld);
		lis_matvech(Ah,ptld,qtld);

		
		/* tmpdot1 = <ptld,q> */
		lis_vector_dotex_mmm(ptld,q,&tmpdot1);

		/* test breakdown */
		if( tmpdot1.hi[0]==0.0 && tmpdot1.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter2;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		
		/* alpha = rho / tmpdot1 */
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)tmpdot1.hi);
		
		/* x = x + alpha*p */
		lis_vector_axpyex_mmm(alpha,p,x);
		
		/* r    = r    - alpha*q    */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,q,r);
		
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
		
		/* rtld = rtld - alpha*qtld */
		lis_vector_axpyex_mmm(alpha,qtld,rtld);

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

/*****************************************
 * Preconditioned BiConjugate Residual   *
 *****************************************
 r(0)    = b - Ax(0)
 rtld(0) = conj(r(0)) or random
 z(0)    = M^-1 * r(0)
 ztld(0) = M^-1 * rtld(0)
 p(0)    = z(0)
 ptld(0) = ztld(0)
 ap(0)   = A * z(0)
 rho(0)  = <ztld(0),ap(0)>
 *****************************************
 for k=1,2,...
   aptld(k-1) = A^H * ptld(k-1)
   map(k-1)   = M^-1 * ap(k-1)
   tmpdot1   = <aptld(k-1),map(k-1)>
   alpha     = rho(k-1) / tmpdot1
   x(k)      = x(k-1) + alpha*ap(k-1)
   r(k)      = r(k-1) - alpha*ap(k-1)
   rtld(k)   = rtld(k-1) - conj(alpha)*aptld(k-1)
   z(k)      = z(k-1) - alpha * map(k-1)
   ztld(k)   = M^-H * rtld(k-1)
   az(k)     = A * z(k)
   rho(k)    = <ztld(k),az(k)>
   beta      = rho(k) / rho(k-1)
   p(k)      = z(k) + beta*p(k-1)
   ptld(k)   = ztld(k) + conj(beta)*ptld(k-1)
   ap(k)     = az(k) + beta*ap(k-1)
 *****************************************/
#undef NWORK
#define NWORK 10
#undef __FUNC__
#define __FUNC__ "lis_bicr_check_params"
LIS_INT lis_bicr_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_bicr_malloc_work"
LIS_INT lis_bicr_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	/*
	err = lis_matrix_convert(solver->A,&solver->Ah,LIS_MATRIX_CCS);
	if( err ) return err;
	*/

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_bicg_malloc_work::work" );
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
#define __FUNC__ "lis_bicr"
LIS_INT lis_bicr(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, z,ztld,p, ptld, ap, map, az, aptld;
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
	z       = solver->work[2];
	ztld    = solver->work[3];
	p       = solver->work[4];
	ptld    = solver->work[5];
	ap      = solver->work[6];
	az      = solver->work[7];
	map     = solver->work[8];
	aptld   = solver->work[9];



	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_psolve(solver, r, z);
	lis_psolveh(solver, rtld, ztld);
	lis_vector_copy(z,p);
	lis_vector_copy(ztld,ptld);
	lis_matvec(A,z,ap);
	lis_vector_dot(ztld,ap,&rho_old);

	for( iter=1; iter<=maxiter; iter++ )
	{
		/* aptld = A^H * ptld */
		/* map   = M^-1 * ap  */
		lis_matvech(A,ptld,aptld);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;

		/* tmpdot1 = <aptld,map> */
		lis_vector_dot(aptld,map,&tmpdot1);
		/* test breakdown */
		if( tmpdot1==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* alpha = rho_old / tmpdot1 */
		/* x     = x + alpha*p   */
		/* r     = r - alpha*ap  */
		alpha = rho_old / tmpdot1;
		lis_vector_axpy(alpha,p,x);
		lis_vector_axpy(-alpha,ap,r);
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
		
		/* rtld = rtld - conj(alpha)*aptld */
		/* z    = z - alpha*map            */
		/* ztld = M^-H * rtld              */
		/* az   = A * z                    */
		/* rho = <ztld,az>                 */
		lis_vector_axpy(-conj(alpha),aptld,rtld);
		lis_vector_axpy(-alpha,map,z);
		time = lis_wtime();
		lis_psolveh(solver, rtld, ztld);
		ptime += lis_wtime()-time;
		lis_matvec(A,z,az);
		lis_vector_dot(ztld,az,&rho);

		/* test breakdown */
		if( rho==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta = rho / rho_old          */
		/* p    = z    + beta*p          */
		/* ptld = ztld + conj(beta)*ptld */
		/* ap   = az   + beta*ap         */
		beta = rho / rho_old;
		lis_vector_xpay(z,beta,p);
		lis_vector_xpay(ztld,conj(beta),ptld);
		lis_vector_xpay(az,beta,ap);

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
#define __FUNC__ "lis_bicr_quad"
LIS_INT lis_bicr_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,rtld, z,ztld,p, ptld, ap, map, az, aptld;
	LIS_QUAD_PTR alpha, beta, rho, rho_old, tmpdot1;
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
	z       = solver->work[2];
	ztld    = solver->work[3];
	p       = solver->work[4];
	ptld    = solver->work[5];
	ap      = solver->work[6];
	az      = solver->work[7];
	map     = solver->work[8];
	aptld   = solver->work[9];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot1,4,1);

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_psolve(solver, r, z);
	lis_psolveh(solver, rtld, ztld);
	lis_vector_copyex_mm(z,p);
	lis_vector_copyex_mm(ztld,ptld);
	lis_matvec(A,z,ap);
	lis_vector_dotex_mmm(ap,ztld,&rho_old);

	for( iter=1; iter<=maxiter; iter++ )
	{
		/* aptld = A^H * ptld */
		/* map   = M^-1 * ap  */
		lis_matvech(A,ptld,aptld);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;

		/* tmpdot1 = <map,aptld> */
		lis_vector_dotex_mmm(map,aptld,&tmpdot1);
		/* test breakdown */
		if( tmpdot1.hi[0]==0.0 && tmpdot1.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* alpha = rho_old / tmpdot1 */
		/* x     = x + alpha*p   */
		/* r     = r - alpha*ap  */
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho_old.hi,(LIS_QUAD *)tmpdot1.hi);
		lis_vector_axpyex_mmm(alpha,p,x);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,ap,r);
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
		
		/* rtld = rtld - alpha*aptld */
		/* z    = z - alpha*map      */
		/* ztld = M^-H * rtld        */
		/* az   = A * z              */
		/* rho = <az,ztld>           */
		lis_vector_axpyex_mmm(alpha,aptld,rtld);
		lis_vector_axpyex_mmm(alpha,map,z);
		time = lis_wtime();
		lis_psolveh(solver, rtld, ztld);
		ptime += lis_wtime()-time;
		lis_matvec(A,z,az);
		lis_vector_dotex_mmm(az,ztld,&rho);

		/* test breakdown */
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta = rho / rho_old    */
		/* p    = z    + beta*p    */
		/* ptld = ztld + beta*ptld */
		/* ap   = az   + beta*ap   */
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
		lis_vector_xpayex_mmm(z,beta,p);
		lis_vector_xpayex_mmm(ztld,beta,ptld);
		lis_vector_xpayex_mmm(az,beta,ap);

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
