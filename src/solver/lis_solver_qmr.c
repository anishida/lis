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
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

#define NWORK 9
/************************************************
 * lis_tfqmr_check_params
 * lis_tfqmr_malloc_work
 * lis_tfqmr
 ************************************************/
#undef __FUNC__
#define __FUNC__ "lis_tfqmr_check_params"
LIS_INT lis_tfqmr_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_tfqmr_malloc_work"
LIS_INT lis_tfqmr_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR), "lis_tfqmr_malloc_work::work" );
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
#define __FUNC__ "lis_tfqmr"
LIS_INT lis_tfqmr(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, u, p, d, t, t1, q, v;
	LIS_REAL tau,w;
	LIS_SCALAR rho, rhoold, theta,eta,beta,alpha,ww,wold,s,c;
	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output;
	double time,ptime;
	LIS_INT m;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	ptime   = 0.0;

	r       = solver->work[0];
	rtld    = solver->work[1];
	u       = solver->work[2];
	p       = solver->work[3];
	d       = solver->work[4];
	t       = solver->work[5];
	t1      = solver->work[6];
	q       = solver->work[7];
	v       = solver->work[8];


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	/* rtld = r
	   p    = r
	   u    = r
	   d    = 0
	   v    = Ap
	*/
	lis_vector_copy(r,p);
	lis_vector_copy(r,u);
	lis_vector_set_all(0,d);
	time = lis_wtime();
	lis_psolve(solver, p, t);
	ptime += lis_wtime()-time;
	lis_matvec(A,t,v);

	/* rhoold = (r,rtld) */
	lis_vector_dot(r,rtld,&rhoold);
	lis_vector_nrm2(r,&tau);
	wold     = tau;
	theta    = 0.0;
	eta      = 0.0;

	/*	tol = tol * bnrm2; */

	iter=1;
	while( iter<=maxiter )
	{
		/* alpha = rho / (v,rtld) */
		lis_vector_dot(v,rtld,&s);
		/* test breakdown */
		if( s==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		alpha = rhoold / s;

		/* q = u - alpha*v         */
		/* r = r - alpha*A(u + q)  */
		lis_vector_axpyz(-alpha,v,u,q);
		lis_vector_axpyz(1.0,u,q,t);

		time = lis_wtime();
		lis_psolve(solver, t, t1);
		ptime += lis_wtime()-time;
		lis_matvec(A,t1,v);
		lis_vector_axpy(-alpha,v,r);

		lis_vector_nrm2(r,&w);
		for(m=0;m<2;m++)
/*		for(m=2*iter-1;m<=2*iter;m++)*/
		{
			if( m==0 )
			{
				ww = sqrt(w*wold);
				lis_vector_xpay(u,theta*theta*eta/alpha,d);
			}
			else
			{
				ww = w;
				lis_vector_xpay(q,theta*theta*eta/alpha,d);
			}
			/*
			if( m%2!=0 )
			{
				ww = w;
				lis_vector_xpay(u,theta*theta*eta/alpha,d);
			}
			else
			{
				ww = sqrt(w*wold);
				lis_vector_xpay(q,theta*theta*eta/alpha,d);
			}
			*/

			theta = ww / tau; 
			c     = 1.0 / sqrt(1.0 + theta*theta);
			eta   = c * c * alpha;
			tau   = tau * theta * c;

			time = lis_wtime();
			lis_psolve(solver, d, t1);
			ptime += lis_wtime()-time;
			lis_vector_axpy(eta,t1,x);

			nrm2 = tau * sqrt(1.0+m) * bnrm2;
			/* convergence check */

			if( m==0 && output )
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
		}

		/* rho  = (r,rtld)              */
		/* beta = rho / rhoold          */
		/* u    = r + beta*q            */
		/* p    = u + beta*(q + beta*p) */
		lis_vector_dot(r,rtld,&rho);
		/* test breakdown */
		if( rho==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		beta = rho / rhoold;
		lis_vector_axpyz(beta,q,r,u);
/*		lis_vector_axpy(beta,p,q);
		lis_vector_axpyz(beta,q,u,p);*/
		lis_vector_xpay(q,beta,p);
		lis_vector_xpay(u,beta,p);
		time = lis_wtime();
		lis_psolve(solver, p, t1);
		ptime += lis_wtime()-time;
		lis_matvec(A,t1,v);

		rhoold = rho;
		wold   = w;
		iter++;
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_tfqmr_quad"
LIS_INT lis_tfqmr_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_PRECON M;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, u, p, d, t, t1, q, v;
	LIS_QUAD_PTR tau,rho, rhoold, theta,eta,beta,alpha,w,ww,wold,s,c,etaold,thetaold,one;
	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output;
	double time,ptime;

	LIS_INT m;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	M       = solver->precon;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	ptime   = 0.0;

	r       = solver->work[0];
	rtld    = solver->work[1];
	u       = solver->work[2];
	p       = solver->work[3];
	d       = solver->work[4];
	t       = solver->work[5];
	t1      = solver->work[6];
	q       = solver->work[7];
	v       = solver->work[8];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rhoold,3,1);
	LIS_QUAD_SCALAR_MALLOC(tau,4,1);
	LIS_QUAD_SCALAR_MALLOC(theta,5,1);
	LIS_QUAD_SCALAR_MALLOC(eta,6,1);
	LIS_QUAD_SCALAR_MALLOC(w,7,1);
	LIS_QUAD_SCALAR_MALLOC(ww,8,1);
	LIS_QUAD_SCALAR_MALLOC(wold,9,1);
	LIS_QUAD_SCALAR_MALLOC(s,10,1);
	LIS_QUAD_SCALAR_MALLOC(c,11,1);
	LIS_QUAD_SCALAR_MALLOC(etaold,12,1);
	LIS_QUAD_SCALAR_MALLOC(thetaold,13,1);
	LIS_QUAD_SCALAR_MALLOC(one,15,1);


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,M,r,rtld,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,rtld,r);

	/* rtld = r
	   p    = r
	   u    = r
	   d    = 0
	   v    = Ap
	*/
	lis_vector_copyex_mm(r,p);
	lis_vector_copyex_mm(r,u);
	lis_vector_set_allex_nm(0.0, d);

	lis_matvec(A,p,t);
	time = lis_wtime();
	lis_psolve(solver, t, v);
	ptime += lis_wtime()-time;

	/* rhoold = (r,rtld) */
	lis_vector_dotex_mmm(r,rtld,&rhoold);
	lis_vector_nrm2ex_mm(r,&tau);
	wold.hi[0]  = tau.hi[0];
	wold.lo[0]  = tau.lo[0];
	theta.hi[0]    = 0.0;
	theta.lo[0]    = 0.0;
	eta.hi[0]    = 0.0;
	eta.lo[0]    = 0.0;
	thetaold.hi[0]    = 0.0;
	thetaold.lo[0]    = 0.0;
	etaold.hi[0]    = 0.0;
	etaold.lo[0]    = 0.0;
	one.hi[0]    = 1.0;
	one.lo[0]    = 0.0;

	/*	tol = tol * bnrm2; */

	iter=1;
	while( iter<=maxiter )
	{
		/* alpha = rho / (v,rtld) */
		lis_vector_dotex_mmm(v,rtld,&s);
		/* test breakdown */
		if( s.hi[0]==0.0 && s.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rhoold.hi,(LIS_QUAD *)s.hi);

		/* q = u - alpha*v         */
		/* r = r - alpha*A(u + q)  */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyzex_mmmm(alpha,v,u,q);
		lis_vector_axpyzex_mmmm(one,u,q,t);

		lis_matvec(A,t,t1);
		time = lis_wtime();
		lis_psolve(solver, t1, v);
		ptime += lis_wtime()-time;
		lis_vector_axpyex_mmm(alpha,v,r);

		lis_vector_nrm2ex_mm(r,&w);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		for(m=0;m<2;m++)
/*		for(m=2*iter-1;m<=2*iter;m++)*/
		{
			if( m==0 )
			{
				lis_quad_mul((LIS_QUAD *)ww.hi,(LIS_QUAD *)w.hi,(LIS_QUAD *)wold.hi);
				lis_quad_sqrt((LIS_QUAD *)ww.hi,(LIS_QUAD *)ww.hi);
				lis_quad_sqr((LIS_QUAD *)theta.hi,(LIS_QUAD *)theta.hi);
				lis_quad_mul((LIS_QUAD *)theta.hi,(LIS_QUAD *)theta.hi,(LIS_QUAD *)eta.hi);
				lis_quad_div((LIS_QUAD *)theta.hi,(LIS_QUAD *)theta.hi,(LIS_QUAD *)alpha.hi);
				lis_vector_xpayex_mmm(q,theta,d);
			}
			else
			{
				ww.hi[0] = w.hi[0];
				ww.lo[0] = w.lo[0];
				lis_quad_sqr((LIS_QUAD *)theta.hi,(LIS_QUAD *)theta.hi);
				lis_quad_mul((LIS_QUAD *)theta.hi,(LIS_QUAD *)theta.hi,(LIS_QUAD *)eta.hi);
				lis_quad_div((LIS_QUAD *)theta.hi,(LIS_QUAD *)theta.hi,(LIS_QUAD *)alpha.hi);
				lis_vector_xpayex_mmm(u,theta,d);
			}

			lis_quad_div((LIS_QUAD *)theta.hi,(LIS_QUAD *)ww.hi,(LIS_QUAD *)tau.hi);

			lis_quad_sqr((LIS_QUAD *)c.hi,(LIS_QUAD *)theta.hi);
			lis_quad_add((LIS_QUAD *)c.hi,(LIS_QUAD *)c.hi,(LIS_QUAD *)one.hi);
			lis_quad_sqrt((LIS_QUAD *)c.hi,(LIS_QUAD *)c.hi);
			lis_quad_div((LIS_QUAD *)c.hi,(LIS_QUAD *)one.hi,(LIS_QUAD *)c.hi);

			lis_quad_sqr((LIS_QUAD *)eta.hi,(LIS_QUAD *)c.hi);
			lis_quad_mul((LIS_QUAD *)eta.hi,(LIS_QUAD *)eta.hi,(LIS_QUAD *)alpha.hi);

			lis_quad_mul((LIS_QUAD *)tau.hi,(LIS_QUAD *)tau.hi,(LIS_QUAD *)theta.hi);
			lis_quad_mul((LIS_QUAD *)tau.hi,(LIS_QUAD *)tau.hi,(LIS_QUAD *)c.hi);
			lis_vector_axpyex_mmm(eta,d,x);


			nrm2 = tau.hi[0] * sqrt(1.0+m) * bnrm2;
			/* convergence check */

			if( m==0 && output )
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
		}

		/* rho  = (r,rtld)              */
		/* beta = rho / rhoold          */
		/* u    = r + beta*q            */
		/* p    = u + beta*(q + beta*p) */
		lis_vector_dotex_mmm(r,rtld,&rho);
		/* test breakdown */
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rhoold.hi);
		lis_vector_axpyzex_mmmm(beta,q,r,u);
		lis_vector_xpayex_mmm(q,beta,p);
		lis_vector_xpayex_mmm(u,beta,p);
		lis_matvec(A,p,t1);
		time = lis_wtime();
		lis_psolve(solver, t1, v);
		ptime += lis_wtime()-time;

		rhoold.hi[0] = rho.hi[0];
		rhoold.lo[0] = rho.lo[0];
		wold.hi[0] = w.hi[0];
		wold.lo[0] = w.lo[0];
		iter++;
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}
#endif
