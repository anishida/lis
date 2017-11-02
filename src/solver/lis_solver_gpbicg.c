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

/********************************************************************
 * Preconditioned Generalized Product type of BiConjugate Gradient  *
 ********************************************************************
 r(0)    = b - Ax(0)
 rtld(0) = conj(r(0)) or random
 p(0)    = M^-1 * r(0)
 rho(0)  = <rtld,r(0)>
 t(-1)   = (0,...,0)^T
 w(-1)   = (0,...,0)^T
 ********************************************************************
 for k=1,2,...
   ap(k-1) = A * p(k-1)
   map(k-1)= M^-1 * ap(k-1)
   tmpdot0 = <rtld(0),ap(k-1)>
   alpha   = rho(k-1) / tmpdot0
   y(k-1)  = t(k-2) - r(k-1) - alpha*w(k-2) + alpha*ap(k-1)
   t(k-1)  = r(k-1) - alpha*ap(k-1)
   mt(k-1) = mr(k-1) - alpha*map(k-1)
   amt(k-1)= A * mt(k-1)
   tmpdot0 = <y(k-1),y(k-1)>
   tmpdot1 = <amt(k-1),t(k-1)>
   tmpdot2 = <y(k-1),t(k-1)>
   tmpdot3 = <amt(k-1),y(k-1)>
   tmpdot4 = <amt(k-1),amt(k-1)>
   tmp     = tmpdot4*tmpdot0-tmpdot3*tmpdot3
   qsi     = (tmpdot0*tmpdot1-tmpdot2*tmpdot3) / tmp
   eta     = (tmpdot4*tmpdot2-tmpdot3*tmpdot1) / tmp
   u(k-1)  = qsi*map(k-1) + eta*(mt(k-2) - mr(k-1) + beta*u(k-2))
   z(k-1)  = qsi*mr(k-1) + eta*z(k-2) - alpha*u(k-1)
   x(k)    = x(k-1) + alpha*p(k-1) + z(k-1)
   r(k)    = t(k-1) - eta*y(k-1) - qsi*amt(k-1)
   mr(k)   = M^-1 * r(k)
   rho(k)  = <rtld,r(k)>
   beta    = (rho(k) / rho(k-1)) * (alpha / qsi)
   w(k)    = amt(k) + beta*ap(k)
   p(k)    = mr(k) + beta*(p(k-1) - u(k-1))
 ********************************************************************/

#define NWORK 14
#undef __FUNC__
#define __FUNC__ "lis_gpbicg_check_params"
LIS_INT lis_gpbicg_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_gpbicg_malloc_work"
LIS_INT lis_gpbicg_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_gpbicg_malloc_work::work" );
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
#define __FUNC__ "lis_gpbicg"
LIS_INT lis_gpbicg(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, mr, p, ap, map, t, mt_old, mt, amt, y, u, w, z;
	LIS_SCALAR alpha, beta, rho, rho_old;
	LIS_SCALAR qsi, eta;
	LIS_SCALAR tmp, tmpdot[5];
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
	mr      = solver->work[2];
	p       = solver->work[3];
	ap      = solver->work[4];
	map     = solver->work[5];
	t       = solver->work[6];
	mt      = solver->work[7];
	amt     = solver->work[8];
	u       = solver->work[9];
	y       = solver->work[10];
	w       = solver->work[11];
	z       = solver->work[12];
	mt_old  = solver->work[13];



	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	time = lis_wtime();
	lis_psolve(solver, r, p);
	ptime += lis_wtime()-time;
	lis_vector_dot(rtld,r,&rho_old);
	lis_vector_set_all(0,t);
	lis_vector_set_all(0,w);
	beta = 0.0;
	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/*   ap(k-1) = A * p(k-1)         */
		/*   map(k-1)= M^-1 * ap(k-1)     */
		/*   tmpdot0 = <rtld(0),map(k-1)> */
		lis_matvec(A,p,ap);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;
		lis_vector_dot(rtld,ap,&tmpdot[0]);

		if( tmpdot[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   alpha   = rho(k-1) / tmpdot0                             */
		/*   y(k-1)  = t(k-2) - r(k-1) - alpha*w(k-2) + alpha*ap(k-1) */
		/*   t(k-1)  = r(k-1) - alpha*ap(k-1)                         */
		alpha = rho_old / tmpdot[0];
		lis_vector_axpyz(-1,w,ap,y);
		lis_vector_xpay(t,alpha,y);
		lis_vector_axpy(-1,r,y);
		lis_vector_axpyz(-alpha,ap,r,t);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](t,solver,&nrm2);
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

		/*   mt(k-1) = mr(k-1) - alpha*map(k-1)                       */
		/*   amt(k-1)= A * mt(k-1)                                    */
		/*   tmpdot0 = <y(k-1),y(k-1)>                                */
		/*   tmpdot1 = <amt(k-1),t(k-1)>                              */
		/*   tmpdot2 = <y(k-1),t(k-1)>                                */
		/*   tmpdot3 = <amt(k-1),y(k-1)>                              */
		/*   tmpdot4 = <amt(k-1),amt(k-1)>                            */
		/*   tmp     = tmpdot4*tmpdot0-tmpdot3*tmpdot3                */
		/*   qsi     = (tmpdot0*tmpdot1-tmpdot2*tmpdot3) / tmp        */
		/*   eta     = (tmpdot4*tmpdot2-tmpdot3*tmpdot1) / tmp        */
		lis_vector_axpyz(-alpha,map,mr,mt);
		lis_matvec(A,mt,amt);
		lis_vector_dot(y,y,&tmpdot[0]);
		lis_vector_dot(amt,t,&tmpdot[1]);
		lis_vector_dot(y,t,&tmpdot[2]);
		lis_vector_dot(amt,y,&tmpdot[3]);
		lis_vector_dot(amt,amt,&tmpdot[4]);
		if(iter==1)
		{
			qsi = tmpdot[1] / tmpdot[4];
			eta = 0.0;
		}
		else
		{
			tmp = tmpdot[4]*tmpdot[0] - tmpdot[3]*tmpdot[3];
			qsi = (tmpdot[0]*tmpdot[1] - tmpdot[2]*tmpdot[3]) / tmp;
			eta = (tmpdot[4]*tmpdot[2] - tmpdot[3]*tmpdot[1]) / tmp;
		}

		/*   u(k-1)  = qsi*map(k-1) + eta*(mt(k-2) - mr(k-1) + beta*u(k-2)) */
		lis_vector_xpay(mt_old,beta,u);
		lis_vector_axpy(-1,mr,u);
		lis_vector_scale(eta,u);
		lis_vector_axpy(qsi,map,u);

		/*   z(k-1)  = qsi*mr(k-1) + eta*z(k-2) - alpha*u(k-1)              */
		lis_vector_scale(eta,z);
		lis_vector_axpy(qsi,mr,z);
		lis_vector_axpy(-alpha,u,z);

		/*   x(k)    = x(k-1) + alpha*p(k-1) + z(k-1)                       */
		lis_vector_axpy(alpha,p,x);
		lis_vector_axpy(1,z,x);
		
		/*   r(k)    = t(k-1) - eta*y(k-1) - qsi*amt(k-1)                   */
		lis_vector_axpyz(-qsi,amt,t,r);
		lis_vector_axpy(-eta,y,r);
		
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

		/*   mr(k)   = M^-1 * r(k)                                          */
		/*   rho(k)  = <rtld,mr(k)>                                         */
		time = lis_wtime();
		lis_psolve(solver, r, mr);
		ptime += lis_wtime()-time;
		lis_vector_dot(rtld,r,&rho);
		if( rho==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   beta    = (rho(k) / rho(k-1)) * (alpha / qsi) */
		beta = (rho / rho_old) * (alpha / qsi);

		/*   w(k)    = amt(k) + beta*ap(k)            */
		/*   p(k)    = mr(k) + beta*(p(k-1) - u(k-1)) */
		lis_vector_axpyz(beta,ap,amt,w);
		lis_vector_axpy(-1,u,p);
		lis_vector_xpay(mr,beta,p);
		
		lis_vector_copy(mt,mt_old);
		rho_old = rho;
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

#if 0
#undef __FUNC__
#define __FUNC__ "lis_gpbicg"
LIS_INT lis_gpbicg(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_PRECON M;
	LIS_VECTOR b,x;
	LIS_VECTOR r, rtld, rhat, p, ptld, phat;
	LIS_VECTOR t, ttld, that, t0, t0hat;
	LIS_VECTOR y, w, u, z;
	LIS_SCALAR alpha, beta, rho, rho_old;
	LIS_SCALAR qsi, eta;
	LIS_SCALAR tmp, tmpdot[5];
	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,n,output,conv;
	double time,ptime;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	M       = solver->precon;
	b       = solver->b;
	x       = solver->x;
	n       = A->n;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;

	rtld    = solver->work[0];
	r       = solver->work[1];
	rhat    = solver->work[2];
	p       = solver->work[3];
	ptld    = solver->work[4];
	phat    = solver->work[5];
	t       = solver->work[6];
	ttld    = solver->work[7];
	that    = solver->work[8];
	t0      = solver->work[9];
	t0hat   = solver->work[10];
	y       = solver->work[11];
	w       = solver->work[12];
	u       = solver->work[13];
	z       = solver->work[14];

	alpha   = (LIS_SCALAR)1.0;
	qsi     = (LIS_SCALAR)1.0;
	rho_old = (LIS_SCALAR)1.0;


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_all(0,ttld);
	lis_vector_set_all(0,ptld);
	lis_vector_set_all(0,p);
	lis_vector_set_all(0,u);
	lis_vector_set_all(0,t);
	lis_vector_set_all(0,t0);
	
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

		/* beta = (rho / rho_old) * (alpha / qsi) */
		beta = (rho / rho_old) * (alpha / qsi);

		/* w = ttld + beta*ptld */
		lis_vector_axpyz(beta,ptld,ttld,w);

		/* rhat = M^-1 * r */
		time = lis_wtime();
		lis_psolve(solver, r, rhat);
		ptime += lis_wtime()-time;

		/* p = rhat + beta*(p - u) */
		lis_vector_axpy(-1,u,p);
		lis_vector_xpay(rhat,beta,p);
		
		/* ptld = A * p */
		lis_matvec(A,p,ptld);

		/* tmpdot[0] = <rtld,ptld> */
		lis_vector_dot(rtld,ptld,&tmpdot[0]);
		/* test breakdown */
		/* */
		
		/* alpha = rho / tmpdot[0] */
		alpha = rho / tmpdot[0];

		/* y = t - r + alpha*(-w + ptld) */
		lis_vector_axpyz(-1,w,ptld,y);
		lis_vector_xpay(t,alpha,y);
		lis_vector_axpy(-1,r,y);

		/* t = r - alpha*ptld */
		lis_vector_axpyz(-alpha,ptld,r,t);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](t,solver,&nrm2);
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

		/* that  = M^-1 * t */
		/* phat  = M^-1 * ptld */
		/* t0hat = M^-1 * t0 */
		time = lis_wtime();
		lis_psolve(solver, t, that);
		lis_psolve(solver, ptld, phat);
		lis_psolve(solver, t0, t0hat);
		ptime += lis_wtime()-time;


		/* ttld = A * that */
		lis_matvec(A,that,ttld);

		/* tmpdot[0] = <y,y>       */
		/* tmpdot[1] = <ttld,t>    */
		/* tmpdot[2] = <y,t>       */
		/* tmpdot[3] = <ttld,y>    */
		/* tmpdot[4] = <ttld,ttld> */
		lis_vector_dot(y,y,&tmpdot[0]);
		lis_vector_dot(ttld,t,&tmpdot[1]);
		lis_vector_dot(y,t,&tmpdot[2]);
		lis_vector_dot(ttld,y,&tmpdot[3]);
		lis_vector_dot(ttld,ttld,&tmpdot[4]);
		if(iter==1)
		{
			qsi = tmpdot[1] / tmpdot[4];
			eta = 0.0;
		}
		else
		{
			tmp = tmpdot[4]*tmpdot[0] - tmpdot[3]*tmpdot[3];
			qsi = (tmpdot[0]*tmpdot[1] - tmpdot[2]*tmpdot[3]) / tmp;
			eta = (tmpdot[4]*tmpdot[2] - tmpdot[3]*tmpdot[1]) / tmp;
		}

		/* u = qsi*phat + eta*(t0hat - rhat + beta*u) */
		lis_vector_xpay(t0hat,beta,u);
		lis_vector_axpy(-1,rhat,u);
		lis_vector_scale(eta,u);
		lis_vector_axpy(qsi,phat,u);

		/* z = qsi*rhat + eta*z - alpha*u */
		lis_vector_scale(eta,z);
		lis_vector_axpy(qsi,rhat,z);
		lis_vector_axpy(-alpha,u,z);

		/* x = x + alpha*p + z */
		lis_vector_axpy(alpha,p,x);
		lis_vector_axpy(1,z,x);
		
		/* r = t - eta*y - qsi*ttld */
		lis_vector_axpyz(-eta,y,t,r);
		lis_vector_axpy(-qsi,ttld,r);
		
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

		lis_vector_copy(t,t0);
		rho_old = rho;
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}
#endif

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_gpbicg_quad"
LIS_INT lis_gpbicg_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, mr, p, ap, map, t, mt_old, mt, amt, y, u, w, z;
	LIS_QUAD_PTR alpha, beta, rho, rho_old;
	LIS_QUAD_PTR qsi, eta, one;
	LIS_QUAD_PTR tmp, tmpdot[5];
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
	mr      = solver->work[2];
	p       = solver->work[3];
	ap      = solver->work[4];
	map     = solver->work[5];
	t       = solver->work[6];
	mt      = solver->work[7];
	amt     = solver->work[8];
	u       = solver->work[9];
	y       = solver->work[10];
	w       = solver->work[11];
	z       = solver->work[12];
	mt_old  = solver->work[13];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(qsi,4,1);
	LIS_QUAD_SCALAR_MALLOC(eta,5,1);
	LIS_QUAD_SCALAR_MALLOC(tmp,6,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[0],7,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[1],8,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[2],9,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[3],10,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[4],11,1);
	LIS_QUAD_SCALAR_MALLOC(one,13,1);


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	time = lis_wtime();
	lis_psolve(solver, r, p);
	ptime += lis_wtime()-time;
	lis_vector_dotex_mmm(rtld,r,&rho_old);
	lis_vector_set_allex_nm(0.0,t);
	lis_vector_set_allex_nm(0.0,w);
	one.hi[0] = -1.0;
	one.lo[0] = 0.0;
	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/*   ap(k-1) = A * p(k-1)         */
		/*   map(k-1)= M^-1 * ap(k-1)     */
		/*   tmpdot0 = <map(k-1),rtld(0)> */
		lis_matvec(A,p,ap);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;
		lis_vector_dotex_mmm(rtld,ap,&tmpdot[0]);

		if( tmpdot[0].hi[0]==0.0 && tmpdot[0].lo[0] )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   alpha   = rho(k-1) / tmpdot0                             */
		/*   y(k-1)  = t(k-2) - r(k-1) - alpha*w(k-2) + alpha*ap(k-1) */
		/*   t(k-1)  = r(k-1) - alpha*ap(k-1)                         */
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho_old.hi,(LIS_QUAD *)tmpdot[0].hi);
		lis_vector_axpyzex_mmmm(one,w,ap,y);
		lis_vector_xpayex_mmm(t,alpha,y);
		lis_vector_axpyex_mmm(one,r,y);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyzex_mmmm(alpha,ap,r,t);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](t,solver,&nrm2);
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

		/*   mt(k-1) = mr(k-1) - alpha*map(k-1)                       */
		/*   amt(k-1)= A * mt(k-1)                                    */
		/*   tmpdot0 = <y(k-1),y(k-1)>                                */
		/*   tmpdot1 = <amt(k-1),t(k-1)>                              */
		/*   tmpdot2 = <y(k-1),t(k-1)>                                */
		/*   tmpdot3 = <amt(k-1),y(k-1)>                              */
		/*   tmpdot4 = <amt(k-1),amt(k-1)>                            */
		/*   tmp     = tmpdot4*tmpdot0-tmpdot3*tmpdot3                */
		/*   qsi     = (tmpdot0*tmpdot1-tmpdot2*tmpdot3) / tmp        */
		/*   eta     = (tmpdot4*tmpdot2-tmpdot3*tmpdot1) / tmp        */
		lis_vector_axpyzex_mmmm(alpha,map,mr,mt);
		lis_matvec(A,mt,amt);
		lis_vector_dotex_mmm(y,y,&tmpdot[0]);
		lis_vector_dotex_mmm(amt,t,&tmpdot[1]);
		lis_vector_dotex_mmm(y,t,&tmpdot[2]);
		lis_vector_dotex_mmm(amt,y,&tmpdot[3]);
		lis_vector_dotex_mmm(amt,amt,&tmpdot[4]);
		if(iter==1)
		{
			lis_quad_div((LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmpdot[1].hi,(LIS_QUAD *)tmpdot[4].hi);
			eta.hi[0] = 0.0;
			eta.lo[0] = 0.0;
		}
		else
		{
			lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmpdot[4].hi,(LIS_QUAD *)tmpdot[0].hi);
			lis_quad_sqr((LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmpdot[3].hi);
			lis_quad_sub((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)qsi.hi);

			lis_quad_mul((LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmpdot[0].hi,(LIS_QUAD *)tmpdot[1].hi);
			lis_quad_mul((LIS_QUAD *)eta.hi,(LIS_QUAD *)tmpdot[2].hi,(LIS_QUAD *)tmpdot[3].hi);
			lis_quad_sub((LIS_QUAD *)qsi.hi,(LIS_QUAD *)qsi.hi,(LIS_QUAD *)eta.hi);
			lis_quad_div((LIS_QUAD *)qsi.hi,(LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmp.hi);

			lis_quad_mul((LIS_QUAD *)eta.hi,(LIS_QUAD *)tmpdot[4].hi,(LIS_QUAD *)tmpdot[2].hi);
			lis_quad_mul((LIS_QUAD *)tmpdot[0].hi,(LIS_QUAD *)tmpdot[3].hi,(LIS_QUAD *)tmpdot[1].hi);
			lis_quad_sub((LIS_QUAD *)eta.hi,(LIS_QUAD *)eta.hi,(LIS_QUAD *)tmpdot[0].hi);
			lis_quad_div((LIS_QUAD *)eta.hi,(LIS_QUAD *)eta.hi,(LIS_QUAD *)tmp.hi);
		}

		/*   u(k-1)  = qsi*map(k-1) + eta*(mt(k-2) - mr(k-1) + beta*u(k-2)) */
		lis_vector_xpayex_mmm(mt_old,beta,u);
		lis_vector_axpyex_mmm(one,mr,u);
		lis_vector_scaleex_mm(eta,u);
		lis_vector_axpyex_mmm(qsi,map,u);

		/*   z(k-1)  = qsi*mr(k-1) + eta*z(k-2) - alpha*u(k-1)              */
		lis_vector_scaleex_mm(eta,z);
		lis_vector_axpyex_mmm(qsi,mr,z);
		lis_vector_axpyex_mmm(alpha,u,z);

		/*   x(k)    = x(k-1) + alpha*p(k-1) + z(k-1)                       */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_quad_minus((LIS_QUAD *)one.hi);
		lis_vector_axpyex_mmm(alpha,p,x);
		lis_vector_axpyex_mmm(one,z,x);
		lis_quad_minus((LIS_QUAD *)one.hi);
		
		/*   r(k)    = t(k-1) - eta*y(k-1) - qsi*amt(k-1)                   */
		lis_quad_minus((LIS_QUAD *)eta.hi);
		lis_quad_minus((LIS_QUAD *)qsi.hi);
		lis_vector_axpyzex_mmmm(qsi,amt,t,r);
		lis_vector_axpyex_mmm(eta,y,r);
		
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

		/*   mr(k)   = M^-1 * r(k)                                          */
		/*   rho(k)  = <rtld,mr(k)>                                         */
		time = lis_wtime();
		lis_psolve(solver, r, mr);
		ptime += lis_wtime()-time;
		lis_vector_dotex_mmm(rtld,r,&rho);
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   beta    = (rho(k) / rho(k-1)) * (alpha / qsi) */
		lis_quad_minus((LIS_QUAD *)eta.hi);
		lis_quad_minus((LIS_QUAD *)qsi.hi);
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)alpha.hi,(LIS_QUAD *)qsi.hi);
		lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,(LIS_QUAD *)tmp.hi);

		/*   w(k)    = amt(k) + beta*ap(k)            */
		/*   p(k)    = mr(k) + beta*(p(k-1) - u(k-1)) */
		lis_vector_axpyzex_mmmm(beta,ap,amt,w);
		lis_vector_axpyex_mmm(one,u,p);
		lis_vector_xpayex_mmm(mr,beta,p);
		
		lis_vector_copyex_mm(mt,mt_old);
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
#define __FUNC__ "lis_gpbicg_switch"
LIS_INT lis_gpbicg_switch(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, mr, p, ap, map, t, mt_old, mt, amt, y, u, w, z;
	LIS_QUAD_PTR alpha, beta, rho, rho_old;
	LIS_QUAD_PTR qsi, eta, one;
	LIS_QUAD_PTR tmp, tmpdot[5];
	LIS_REAL bnrm2, nrm2, tol;
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
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	tol      = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	ptime   = 0.0;

	rtld    = solver->work[0];
	r       = solver->work[1];
	mr      = solver->work[2];
	p       = solver->work[3];
	ap      = solver->work[4];
	map     = solver->work[5];
	t       = solver->work[6];
	mt      = solver->work[7];
	amt     = solver->work[8];
	u       = solver->work[9];
	y       = solver->work[10];
	w       = solver->work[11];
	z       = solver->work[12];
	mt_old  = solver->work[13];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(qsi,4,1);
	LIS_QUAD_SCALAR_MALLOC(eta,5,1);
	LIS_QUAD_SCALAR_MALLOC(tmp,6,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[0],7,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[1],8,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[2],9,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[3],10,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[4],11,1);
	LIS_QUAD_SCALAR_MALLOC(one,13,1);
	one.hi[0] = -1.0;
	one.lo[0] = 0.0;

	r->precision = LIS_PRECISION_DEFAULT;
	p->precision = LIS_PRECISION_DEFAULT;
	mt->precision = LIS_PRECISION_DEFAULT;
	ap->precision = LIS_PRECISION_DEFAULT;


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	time = lis_wtime();
	lis_psolve(solver, r, p);
	ptime += lis_wtime()-time;
	lis_vector_dot(rtld,r,&rho_old.hi[0]);
	lis_vector_set_all(0,t);
	lis_vector_set_all(0,w);
	
	for( iter=1; iter<=maxiter2; iter++ )
	{
		/*   ap(k-1) = A * p(k-1)         */
		/*   map(k-1)= M^-1 * ap(k-1)     */
		/*   tmpdot0 = <map(k-1),rtld(0)> */
		lis_matvec(A,p,ap);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;
		lis_vector_dot(rtld,ap,&tmpdot[0].hi[0]);

		if( tmpdot[0].hi[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->iter2     = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   alpha   = rho(k-1) / tmpdot0                             */
		/*   y(k-1)  = t(k-2) - r(k-1) - alpha*w(k-2) + alpha*ap(k-1) */
		/*   t(k-1)  = r(k-1) - alpha*ap(k-1)                         */
		alpha.hi[0] = rho_old.hi[0] / tmpdot[0].hi[0];
		lis_vector_axpyz(-1,w,ap,y);
		lis_vector_xpay(t,alpha.hi[0],y);
		lis_vector_axpy(-1,r,y);
		lis_vector_axpyz(-alpha.hi[0],ap,r,t);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](t,solver,&nrm2);
		if( nrm2 <= tol )
		{
			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			lis_vector_axpy(alpha.hi[0],p,x);
			solver->iter       = iter;
			solver->iter2      = iter;
			solver->ptime      = ptime;
			break;
		}

		/*   mt(k-1) = mr(k-1) - alpha*map(k-1)                       */
		/*   amt(k-1)= A * mt(k-1)                                    */
		/*   tmpdot0 = <y(k-1),y(k-1)>                                */
		/*   tmpdot1 = <amt(k-1),t(k-1)>                              */
		/*   tmpdot2 = <y(k-1),t(k-1)>                                */
		/*   tmpdot3 = <amt(k-1),y(k-1)>                              */
		/*   tmpdot4 = <amt(k-1),amt(k-1)>                            */
		/*   tmp     = tmpdot4*tmpdot0-tmpdot3*tmpdot3                */
		/*   qsi     = (tmpdot0*tmpdot1-tmpdot2*tmpdot3) / tmp        */
		/*   eta     = (tmpdot4*tmpdot2-tmpdot3*tmpdot1) / tmp        */
		lis_vector_axpyz(-alpha.hi[0],map,mr,mt);
		lis_matvec(A,mt,amt);
		lis_vector_dot(y,y,&tmpdot[0].hi[0]);
		lis_vector_dot(amt,t,&tmpdot[1].hi[0]);
		lis_vector_dot(y,t,&tmpdot[2].hi[0]);
		lis_vector_dot(amt,y,&tmpdot[3].hi[0]);
		lis_vector_dot(amt,amt,&tmpdot[4].hi[0]);
		if(iter==1)
		{
			qsi.hi[0] = tmpdot[1].hi[0] / tmpdot[4].hi[0];
			eta.hi[0] = 0.0;
		}
		else
		{
			tmp.hi[0] = tmpdot[4].hi[0]*tmpdot[0].hi[0]  - tmpdot[3].hi[0]*tmpdot[3].hi[0];
			qsi.hi[0] = (tmpdot[0].hi[0]*tmpdot[1].hi[0] - tmpdot[2].hi[0]*tmpdot[3].hi[0]) / tmp.hi[0];
			eta.hi[0] = (tmpdot[4].hi[0]*tmpdot[2].hi[0] - tmpdot[3].hi[0]*tmpdot[1].hi[0]) / tmp.hi[0];
		}

		/*   u(k-1)  = qsi*map(k-1) + eta*(mt(k-2) - mr(k-1) + beta*u(k-2)) */
		lis_vector_xpay(mt_old,beta.hi[0],u);
		lis_vector_axpy(-1,mr,u);
		lis_vector_scale(eta.hi[0],u);
		lis_vector_axpy(qsi.hi[0],map,u);

		/*   z(k-1)  = qsi*mr(k-1) + eta*z(k-2) - alpha*u(k-1)              */
		lis_vector_scale(eta.hi[0],z);
		lis_vector_axpy(qsi.hi[0],mr,z);
		lis_vector_axpy(-alpha.hi[0],u,z);

		/*   x(k)    = x(k-1) + alpha*p(k-1) + z(k-1)                       */
		lis_vector_axpy(alpha.hi[0],p,x);
		lis_vector_axpy(1,z,x);
		
		/*   r(k)    = t(k-1) - eta*y(k-1) - qsi*amt(k-1)                   */
		lis_vector_axpyz(-qsi.hi[0],amt,t,r);
		lis_vector_axpy(-eta.hi[0],y,r);
		
		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);
		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
			if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
		}

		if( tol >= nrm2 )
		{
			solver->iter       = iter;
			solver->iter2      = iter;
			solver->ptime      = ptime;
			break;
		}

		/*   mr(k)   = M^-1 * r(k)                                          */
		/*   rho(k)  = <rtld,mr(k)>                                         */
		time = lis_wtime();
		lis_psolve(solver, r, mr);
		ptime += lis_wtime()-time;
		lis_vector_dot(rtld,r,&rho.hi[0]);
		if( rho.hi[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->iter2     = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   beta    = (rho(k) / rho(k-1)) * (alpha / qsi) */
		beta.hi[0] = (rho.hi[0] / rho_old.hi[0]) * (alpha.hi[0] / qsi.hi[0]);

		/*   w(k)    = amt(k) + beta*ap(k)            */
		/*   p(k)    = mr(k) + beta*(p(k-1) - u(k-1)) */
		lis_vector_axpyz(beta.hi[0],ap,amt,w);
		lis_vector_axpy(-1,u,p);
		lis_vector_xpay(mr,beta.hi[0],p);
		
		lis_vector_copy(mt,mt_old);
		rho_old.hi[0] = rho.hi[0];
	}


	r->precision = LIS_PRECISION_QUAD;
	p->precision = LIS_PRECISION_QUAD;
	mt->precision = LIS_PRECISION_QUAD;
	ap->precision = LIS_PRECISION_QUAD;

	solver->options[LIS_OPTIONS_INITGUESS_ZEROS] = LIS_FALSE;
	lis_vector_copyex_mn(x,solver->xx);


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	time = lis_wtime();
	lis_psolve(solver, r, p);
	ptime += lis_wtime()-time;
	lis_vector_dotex_mmm(rtld,r,&rho_old);
	lis_vector_set_allex_nm(0.0,t);
	lis_vector_set_allex_nm(0.0,w);

	for( iter2=iter+1; iter2<=maxiter; iter2++ )
	{
		/*   ap(k-1) = A * p(k-1)         */
		/*   map(k-1)= M^-1 * ap(k-1)     */
		/*   tmpdot0 = <map(k-1),rtld(0)> */
		lis_matvec(A,p,ap);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;
		lis_vector_dotex_mmm(rtld,ap,&tmpdot[0]);

		if( tmpdot[0].hi[0]==0.0 && tmpdot[0].lo[0] )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter2;
			solver->iter2     = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   alpha   = rho(k-1) / tmpdot0                             */
		/*   y(k-1)  = t(k-2) - r(k-1) - alpha*w(k-2) + alpha*ap(k-1) */
		/*   t(k-1)  = r(k-1) - alpha*ap(k-1)                         */
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho_old.hi,(LIS_QUAD *)tmpdot[0].hi);
		lis_vector_axpyzex_mmmm(one,w,ap,y);
		lis_vector_xpayex_mmm(t,alpha,y);
		lis_vector_axpyex_mmm(one,r,y);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyzex_mmmm(alpha,ap,r,t);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](t,solver,&nrm2);
		if( nrm2 <= tol )
		{
			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter2] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			lis_quad_minus((LIS_QUAD *)alpha.hi);
			lis_vector_axpyex_mmm(alpha,p,x);
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter2;
			solver->iter2      = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		/*   mt(k-1) = mr(k-1) - alpha*map(k-1)                       */
		/*   amt(k-1)= A * mt(k-1)                                    */
		/*   tmpdot0 = <y(k-1),y(k-1)>                                */
		/*   tmpdot1 = <amt(k-1),t(k-1)>                              */
		/*   tmpdot2 = <y(k-1),t(k-1)>                                */
		/*   tmpdot3 = <amt(k-1),y(k-1)>                              */
		/*   tmpdot4 = <amt(k-1),amt(k-1)>                            */
		/*   tmp     = tmpdot4*tmpdot0-tmpdot3*tmpdot3                */
		/*   qsi     = (tmpdot0*tmpdot1-tmpdot2*tmpdot3) / tmp        */
		/*   eta     = (tmpdot4*tmpdot2-tmpdot3*tmpdot1) / tmp        */
		lis_vector_axpyzex_mmmm(alpha,map,mr,mt);
		lis_matvec(A,mt,amt);
		lis_vector_dotex_mmm(y,y,&tmpdot[0]);
		lis_vector_dotex_mmm(amt,t,&tmpdot[1]);
		lis_vector_dotex_mmm(y,t,&tmpdot[2]);
		lis_vector_dotex_mmm(amt,y,&tmpdot[3]);
		lis_vector_dotex_mmm(amt,amt,&tmpdot[4]);
		if(iter==1)
		{
			lis_quad_div((LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmpdot[1].hi,(LIS_QUAD *)tmpdot[4].hi);
			eta.hi[0] = 0.0;
			eta.lo[0] = 0.0;
		}
		else
		{
			lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmpdot[4].hi,(LIS_QUAD *)tmpdot[0].hi);
			lis_quad_sqr((LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmpdot[3].hi);
			lis_quad_sub((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)qsi.hi);

			lis_quad_mul((LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmpdot[0].hi,(LIS_QUAD *)tmpdot[1].hi);
			lis_quad_mul((LIS_QUAD *)eta.hi,(LIS_QUAD *)tmpdot[2].hi,(LIS_QUAD *)tmpdot[3].hi);
			lis_quad_sub((LIS_QUAD *)qsi.hi,(LIS_QUAD *)qsi.hi,(LIS_QUAD *)eta.hi);
			lis_quad_div((LIS_QUAD *)qsi.hi,(LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmp.hi);

			lis_quad_mul((LIS_QUAD *)eta.hi,(LIS_QUAD *)tmpdot[4].hi,(LIS_QUAD *)tmpdot[2].hi);
			lis_quad_mul((LIS_QUAD *)tmpdot[0].hi,(LIS_QUAD *)tmpdot[3].hi,(LIS_QUAD *)tmpdot[1].hi);
			lis_quad_sub((LIS_QUAD *)eta.hi,(LIS_QUAD *)eta.hi,(LIS_QUAD *)tmpdot[0].hi);
			lis_quad_div((LIS_QUAD *)eta.hi,(LIS_QUAD *)eta.hi,(LIS_QUAD *)tmp.hi);
		}

		/*   u(k-1)  = qsi*map(k-1) + eta*(mt(k-2) - mr(k-1) + beta*u(k-2)) */
		lis_vector_xpayex_mmm(mt_old,beta,u);
		lis_vector_axpyex_mmm(one,mr,u);
		lis_vector_scaleex_mm(eta,u);
		lis_vector_axpyex_mmm(qsi,map,u);

		/*   z(k-1)  = qsi*mr(k-1) + eta*z(k-2) - alpha*u(k-1)              */
		lis_vector_scaleex_mm(eta,z);
		lis_vector_axpyex_mmm(qsi,mr,z);
		lis_vector_axpyex_mmm(alpha,u,z);

		/*   x(k)    = x(k-1) + alpha*p(k-1) + z(k-1)                       */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_quad_minus((LIS_QUAD *)one.hi);
		lis_vector_axpyex_mmm(alpha,p,x);
		lis_vector_axpyex_mmm(one,z,x);
		lis_quad_minus((LIS_QUAD *)one.hi);
		
		/*   r(k)    = t(k-1) - eta*y(k-1) - qsi*amt(k-1)                   */
		lis_quad_minus((LIS_QUAD *)eta.hi);
		lis_quad_minus((LIS_QUAD *)qsi.hi);
		lis_vector_axpyzex_mmmm(qsi,amt,t,r);
		lis_vector_axpyex_mmm(eta,y,r);
		
		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);
		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[iter2] = nrm2;
			if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
		}

		if( tol >= nrm2 )
		{
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter2;
			solver->iter2      = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		/*   mr(k)   = M^-1 * r(k)                                          */
		/*   rho(k)  = <rtld,mr(k)>                                         */
		time = lis_wtime();
		lis_psolve(solver, r, mr);
		ptime += lis_wtime()-time;
		lis_vector_dotex_mmm(rtld,r,&rho);
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter2;
			solver->iter2     = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   beta    = (rho(k) / rho(k-1)) * (alpha / qsi) */
		lis_quad_minus((LIS_QUAD *)eta.hi);
		lis_quad_minus((LIS_QUAD *)qsi.hi);
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)alpha.hi,(LIS_QUAD *)qsi.hi);
		lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,(LIS_QUAD *)tmp.hi);

		/*   w(k)    = amt(k) + beta*ap(k)            */
		/*   p(k)    = mr(k) + beta*(p(k-1) - u(k-1)) */
		lis_vector_axpyzex_mmmm(beta,ap,amt,w);
		lis_vector_axpyex_mmm(one,u,p);
		lis_vector_xpayex_mmm(mr,beta,p);
		
		lis_vector_copyex_mm(mt,mt_old);
		rho_old.hi[0] = rho.hi[0];
		rho_old.lo[0] = rho.lo[0];
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->iter2     = iter2;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}
#endif

/********************************************************************
 * Preconditioned Generalized Product type of BiConjugate Residual  *
 ********************************************************************
 r(0)    = b - Ax(0)
 rtld(0) = conj(r(0)) or random
 rtld(0) = A^H * rtld(0)
 p(0)    = M^-1 * r(0)
 rho(0)  = <rtld,p(0)>
 t(-1)   = (0,...,0)^T
 w(-1)   = (0,...,0)^T
 ********************************************************************
 for k=1,2,...
   ap(k-1) = A * p(k-1)
   map(k-1)= M^-1 * ap(k-1)
   tmpdot0 = <map(k-1),rtld(0)>
   alpha   = rho(k-1) / tmpdot0
   y(k-1)  = t(k-2) - r(k-1) - alpha*w(k-2) + alpha*ap(k-1)
   t(k-1)  = r(k-1) - alpha*ap(k-1)
   mt(k-1) = mr(k-1) - alpha*map(k-1)
   amt(k-1)= A * mt(k-1)
   tmpdot0 = <y(k-1),y(k-1)>
   tmpdot1 = <amt(k-1),t(k-1)>
   tmpdot2 = <y(k-1),t(k-1)>
   tmpdot3 = <amt(k-1),y(k-1)>
   tmpdot4 = <amt(k-1),amt(k-1)>
   tmp     = tmpdot4*tmpdot0-tmpdot3*tmpdot3
   qsi     = (tmpdot0*tmpdot1-tmpdot2*tmpdot3) / tmp
   eta     = (tmpdot4*tmpdot2-tmpdot3*tmpdot1) / tmp
   u(k-1)  = qsi*map(k-1) + eta*(mt(k-2) - mr(k-1) + beta*u(k-2))
   z(k-1)  = qsi*mr(k-1) + eta*z(k-2) - alpha*u(k-1)
   x(k)    = x(k-1) + alpha*p(k-1) + z(k-1)
   r(k)    = t(k-1) - eta*y(k-1) - qsi*amt(k-1)
   mr(k)   = M^-1 * r(k)
   rho(k)  = <rtld,mr(k)>
   beta    = (rho(k) / rho(k-1)) * (alpha / qsi)
   w(k)    = amt(k) + beta*ap(k)
   p(k)    = mr(k) + beta*(p(k-1) - u(k-1))
 ********************************************************************/
#undef NWORK
#define NWORK 14
#undef __FUNC__
#define __FUNC__ "lis_gpbicr_check_params"
LIS_INT lis_gpbicr_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_gpbicg_malloc_work"
LIS_INT lis_gpbicr_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_gpbicg_malloc_work::work" );
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
#define __FUNC__ "lis_gpbicg"
LIS_INT lis_gpbicr(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, mr, p, ap, map, t, mt_old, mt, amt, y, u, w, z;
	LIS_SCALAR alpha, beta, rho, rho_old;
	LIS_SCALAR qsi, eta;
	LIS_SCALAR tmp, tmpdot[5];
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
	mr      = solver->work[2];
	p       = solver->work[3];
	ap      = solver->work[4];
	map     = solver->work[5];
	t       = solver->work[6];
	mt      = solver->work[7];
	amt     = solver->work[8];
	u       = solver->work[9];
	y       = solver->work[10];
	w       = solver->work[11];
	z       = solver->work[12];
	mt_old  = solver->work[13];



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
	lis_psolve(solver, r, p);
	ptime += lis_wtime()-time;
	lis_vector_dot(rtld,p,&rho_old);
	lis_vector_set_all(0,t);
	lis_vector_set_all(0,w);
	beta = 0.0;
	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/*   ap(k-1) = A * p(k-1)         */
		/*   map(k-1)= M^-1 * ap(k-1)     */
		/*   tmpdot0 = <map(k-1),rtld(0)> */
		lis_matvec(A,p,ap);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;
		lis_vector_dot(rtld,map,&tmpdot[0]);

		if( tmpdot[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   alpha   = rho(k-1) / tmpdot0                             */
		/*   y(k-1)  = t(k-2) - r(k-1) - alpha*w(k-2) + alpha*ap(k-1) */
		/*   t(k-1)  = r(k-1) - alpha*ap(k-1)                         */
		alpha = rho_old / tmpdot[0];
		lis_vector_axpyz(-1,w,ap,y);
		lis_vector_xpay(t,alpha,y);
		lis_vector_axpy(-1,r,y);
		lis_vector_axpyz(-alpha,ap,r,t);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](t,solver,&nrm2);
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

		/*   mt(k-1) = mr(k-1) - alpha*map(k-1)                       */
		/*   amt(k-1)= A * mt(k-1)                                    */
		/*   tmpdot0 = <y(k-1),y(k-1)>                                */
		/*   tmpdot1 = <amt(k-1),t(k-1)>                              */
		/*   tmpdot2 = <y(k-1),t(k-1)>                                */
		/*   tmpdot3 = <amt(k-1),y(k-1)>                              */
		/*   tmpdot4 = <amt(k-1),amt(k-1)>                            */
		/*   tmp     = tmpdot4*tmpdot0-tmpdot3*tmpdot3                */
		/*   qsi     = (tmpdot0*tmpdot1-tmpdot2*tmpdot3) / tmp        */
		/*   eta     = (tmpdot4*tmpdot2-tmpdot3*tmpdot1) / tmp        */
		lis_vector_axpyz(-alpha,map,mr,mt);
		lis_matvec(A,mt,amt);
		lis_vector_dot(y,y,&tmpdot[0]);
		lis_vector_dot(amt,t,&tmpdot[1]);
		lis_vector_dot(y,t,&tmpdot[2]);
		lis_vector_dot(amt,y,&tmpdot[3]);
		lis_vector_dot(amt,amt,&tmpdot[4]);
		if(iter==1)
		{
			qsi = tmpdot[1] / tmpdot[4];
			eta = 0.0;
		}
		else
		{
			tmp = tmpdot[4]*tmpdot[0] - tmpdot[3]*tmpdot[3];
			qsi = (tmpdot[0]*tmpdot[1] - tmpdot[2]*tmpdot[3]) / tmp;
			eta = (tmpdot[4]*tmpdot[2] - tmpdot[3]*tmpdot[1]) / tmp;
		}

		/*   u(k-1)  = qsi*map(k-1) + eta*(mt(k-2) - mr(k-1) + beta*u(k-2)) */
		lis_vector_xpay(mt_old,beta,u);
		lis_vector_axpy(-1,mr,u);
		lis_vector_scale(eta,u);
		lis_vector_axpy(qsi,map,u);

		/*   z(k-1)  = qsi*mr(k-1) + eta*z(k-2) - alpha*u(k-1)              */
		lis_vector_scale(eta,z);
		lis_vector_axpy(qsi,mr,z);
		lis_vector_axpy(-alpha,u,z);

		/*   x(k)    = x(k-1) + alpha*p(k-1) + z(k-1)                       */
		lis_vector_axpy(alpha,p,x);
		lis_vector_axpy(1,z,x);
		
		/*   r(k)    = t(k-1) - eta*y(k-1) - qsi*amt(k-1)                   */
		lis_vector_axpyz(-qsi,amt,t,r);
		lis_vector_axpy(-eta,y,r);
		
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

		/*   mr(k)   = M^-1 * r(k)                                          */
		/*   rho(k)  = <rtld,mr(k)>                                         */
		time = lis_wtime();
		lis_psolve(solver, r, mr);
		ptime += lis_wtime()-time;
		lis_vector_dot(rtld,mr,&rho);
		if( rho==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   beta    = (rho(k) / rho(k-1)) * (alpha / qsi) */
		beta = (rho / rho_old) * (alpha / qsi);

		/*   w(k)    = amt(k) + beta*ap(k)            */
		/*   p(k)    = mr(k) + beta*(p(k-1) - u(k-1)) */
		lis_vector_axpyz(beta,ap,amt,w);
		lis_vector_axpy(-1,u,p);
		lis_vector_xpay(mr,beta,p);
		
		lis_vector_copy(mt,mt_old);
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
#define __FUNC__ "lis_gpbicg_quad"
LIS_INT lis_gpbicr_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, mr, p, ap, map, t, mt_old, mt, amt, y, u, w, z;
	LIS_QUAD_PTR alpha, beta, rho, rho_old;
	LIS_QUAD_PTR qsi, eta, one;
	LIS_QUAD_PTR tmp, tmpdot[5];
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
	mr      = solver->work[2];
	p       = solver->work[3];
	ap      = solver->work[4];
	map     = solver->work[5];
	t       = solver->work[6];
	mt      = solver->work[7];
	amt     = solver->work[8];
	u       = solver->work[9];
	y       = solver->work[10];
	w       = solver->work[11];
	z       = solver->work[12];
	mt_old  = solver->work[13];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho_old,3,1);
	LIS_QUAD_SCALAR_MALLOC(qsi,4,1);
	LIS_QUAD_SCALAR_MALLOC(eta,5,1);
	LIS_QUAD_SCALAR_MALLOC(tmp,6,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[0],7,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[1],8,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[2],9,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[3],10,1);
	LIS_QUAD_SCALAR_MALLOC(tmpdot[4],11,1);
	LIS_QUAD_SCALAR_MALLOC(one,13,1);


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
	lis_psolve(solver, r, p);
	ptime += lis_wtime()-time;
	lis_vector_dotex_mmm(rtld,p,&rho_old);
	lis_vector_set_allex_nm(0.0,t);
	lis_vector_set_allex_nm(0.0,w);
	one.hi[0] = -1.0;
	one.lo[0] = 0.0;
	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/*   ap(k-1) = A * p(k-1)         */
		/*   map(k-1)= M^-1 * ap(k-1)     */
		/*   tmpdot0 = <map(k-1),rtld(0)> */
		lis_matvec(A,p,ap);
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;
		lis_vector_dotex_mmm(rtld,map,&tmpdot[0]);

		if( tmpdot[0].hi[0]==0.0 && tmpdot[0].lo[0] )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   alpha   = rho(k-1) / tmpdot0                             */
		/*   y(k-1)  = t(k-2) - r(k-1) - alpha*w(k-2) + alpha*ap(k-1) */
		/*   t(k-1)  = r(k-1) - alpha*ap(k-1)                         */
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho_old.hi,(LIS_QUAD *)tmpdot[0].hi);
		lis_vector_axpyzex_mmmm(one,w,ap,y);
		lis_vector_xpayex_mmm(t,alpha,y);
		lis_vector_axpyex_mmm(one,r,y);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyzex_mmmm(alpha,ap,r,t);

		/* Early check for tolerance */
		lis_solver_get_residual[conv](t,solver,&nrm2);
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

		/*   mt(k-1) = mr(k-1) - alpha*map(k-1)                       */
		/*   amt(k-1)= A * mt(k-1)                                    */
		/*   tmpdot0 = <y(k-1),y(k-1)>                                */
		/*   tmpdot1 = <amt(k-1),t(k-1)>                              */
		/*   tmpdot2 = <y(k-1),t(k-1)>                                */
		/*   tmpdot3 = <amt(k-1),y(k-1)>                              */
		/*   tmpdot4 = <amt(k-1),amt(k-1)>                            */
		/*   tmp     = tmpdot4*tmpdot0-tmpdot3*tmpdot3                */
		/*   qsi     = (tmpdot0*tmpdot1-tmpdot2*tmpdot3) / tmp        */
		/*   eta     = (tmpdot4*tmpdot2-tmpdot3*tmpdot1) / tmp        */
		lis_vector_axpyzex_mmmm(alpha,map,mr,mt);
		lis_matvec(A,mt,amt);
		lis_vector_dotex_mmm(y,y,&tmpdot[0]);
		lis_vector_dotex_mmm(amt,t,&tmpdot[1]);
		lis_vector_dotex_mmm(y,t,&tmpdot[2]);
		lis_vector_dotex_mmm(amt,y,&tmpdot[3]);
		lis_vector_dotex_mmm(amt,amt,&tmpdot[4]);
		if(iter==1)
		{
			lis_quad_div((LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmpdot[1].hi,(LIS_QUAD *)tmpdot[4].hi);
			eta.hi[0] = 0.0;
			eta.lo[0] = 0.0;
		}
		else
		{
			lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmpdot[4].hi,(LIS_QUAD *)tmpdot[0].hi);
			lis_quad_sqr((LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmpdot[3].hi);
			lis_quad_sub((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)qsi.hi);

			lis_quad_mul((LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmpdot[0].hi,(LIS_QUAD *)tmpdot[1].hi);
			lis_quad_mul((LIS_QUAD *)eta.hi,(LIS_QUAD *)tmpdot[2].hi,(LIS_QUAD *)tmpdot[3].hi);
			lis_quad_sub((LIS_QUAD *)qsi.hi,(LIS_QUAD *)qsi.hi,(LIS_QUAD *)eta.hi);
			lis_quad_div((LIS_QUAD *)qsi.hi,(LIS_QUAD *)qsi.hi,(LIS_QUAD *)tmp.hi);

			lis_quad_mul((LIS_QUAD *)eta.hi,(LIS_QUAD *)tmpdot[4].hi,(LIS_QUAD *)tmpdot[2].hi);
			lis_quad_mul((LIS_QUAD *)tmpdot[0].hi,(LIS_QUAD *)tmpdot[3].hi,(LIS_QUAD *)tmpdot[1].hi);
			lis_quad_sub((LIS_QUAD *)eta.hi,(LIS_QUAD *)eta.hi,(LIS_QUAD *)tmpdot[0].hi);
			lis_quad_div((LIS_QUAD *)eta.hi,(LIS_QUAD *)eta.hi,(LIS_QUAD *)tmp.hi);
		}

		/*   u(k-1)  = qsi*map(k-1) + eta*(mt(k-2) - mr(k-1) + beta*u(k-2)) */
		lis_vector_xpayex_mmm(mt_old,beta,u);
		lis_vector_axpyex_mmm(one,mr,u);
		lis_vector_scaleex_mm(eta,u);
		lis_vector_axpyex_mmm(qsi,map,u);

		/*   z(k-1)  = qsi*mr(k-1) + eta*z(k-2) - alpha*u(k-1)              */
		lis_vector_scaleex_mm(eta,z);
		lis_vector_axpyex_mmm(qsi,mr,z);
		lis_vector_axpyex_mmm(alpha,u,z);

		/*   x(k)    = x(k-1) + alpha*p(k-1) + z(k-1)                       */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_quad_minus((LIS_QUAD *)one.hi);
		lis_vector_axpyex_mmm(alpha,p,x);
		lis_vector_axpyex_mmm(one,z,x);
		lis_quad_minus((LIS_QUAD *)one.hi);
		
		/*   r(k)    = t(k-1) - eta*y(k-1) - qsi*amt(k-1)                   */
		lis_quad_minus((LIS_QUAD *)eta.hi);
		lis_quad_minus((LIS_QUAD *)qsi.hi);
		lis_vector_axpyzex_mmmm(qsi,amt,t,r);
		lis_vector_axpyex_mmm(eta,y,r);
		
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

		/*   mr(k)   = M^-1 * r(k)                                          */
		/*   rho(k)  = <rtld,mr(k)>                                         */
		time = lis_wtime();
		lis_psolve(solver, r, mr);
		ptime += lis_wtime()-time;
		lis_vector_dotex_mmm(rtld,mr,&rho);
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/*   beta    = (rho(k) / rho(k-1)) * (alpha / qsi) */
		lis_quad_minus((LIS_QUAD *)eta.hi);
		lis_quad_minus((LIS_QUAD *)qsi.hi);
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)alpha.hi,(LIS_QUAD *)qsi.hi);
		lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,(LIS_QUAD *)tmp.hi);

		/*   w(k)    = amt(k) + beta*ap(k)            */
		/*   p(k)    = mr(k) + beta*(p(k-1) - u(k-1)) */
		lis_vector_axpyzex_mmmm(beta,ap,amt,w);
		lis_vector_axpyex_mmm(one,u,p);
		lis_vector_xpayex_mmm(mr,beta,p);
		
		lis_vector_copyex_mm(mt,mt_old);
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
