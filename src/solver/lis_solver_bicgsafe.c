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

/********************************************
 * Preconditioned BiConjugate Gradient Safe *
 ********************************************
 r(0)    = b - Ax(0)
 rtld(0) = conj(r(0)) or random
 mr(0)   = M^-1 * r(0)
 amr(0)  = A * mr(0)
 p(0)    = mr(0)
 ap(0)   = amr(0)
 rho(0)  = <rtld,r(0)>
 ********************************************
 for k=1,2,...
   tmpdot0 = <ap(k-1),rtld(0)>
   alpha   = rho(k-1) / tmpdot0
   tmpdot0 = <y(k-1),y(k-1)>
   tmpdot1 = <amr(k-1),r(k-1)>
   tmpdot2 = <y(k-1),r(k-1)>
   tmpdot3 = <amr(k-1),y(k-1)>
   tmpdot4 = <amr(k-1),amr(k-1)>
   tmp     = tmpdot4*tmpdot0-tmpdot3*tmpdot3
   qsi     = (tmpdot0*tmpdot1-tmpdot2*tmpdot3) / tmp
   eta     = (tmpdot4*tmpdot2-tmpdot3*tmpdot1) / tmp
   t(k-1)  = qsi*ap(k-1) + eta*y(k-1)
   mt(k-1) = M^-1 * t(k-1)
   u(k-1)  = mt(k-1) + eta*beta*u(k-2)
   au(k-1) = A * u(k-1)
   z(k-1)  = qsi*mr(k-1) + eta*z(k-2) - alpha*u(k-1)
   y(k)    = qsi*amr(k-1) + eta*y(k-1) - alpha*au(k-1)
   x(k)    = x(k-1) + alpha*p(k-1) + z(k-1)
   r(k)    = r(k-1) - alpha*ap(k-1) - y(k)
   mr(k)   = M^-1 * r(k)
   amr(k)  = A * mr(k)
   rho(k)  = <rtld,r(k)>
   beta    = (rho(k) / rho(k-1)) * (alpha / qsi)
   p(k)    = mr(k) + beta*(p(k-1) - u(k-1))
   ap(k-1) = amr(k) + beta*(ap(k-1) - au(k-1))
 ********************************************/

#define NWORK 12
#undef __FUNC__
#define __FUNC__ "lis_bicgsafe_check_params"
LIS_INT lis_bicgsafe_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_bicgsafe_malloc_work"
LIS_INT lis_bicgsafe_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_bicgsafe_malloc_work::work" );
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
#define __FUNC__ "lis_bicgsafe"
LIS_INT lis_bicgsafe(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, mr, amr, t, mt, p, ap;
	LIS_VECTOR y, u, au, z;
	LIS_SCALAR alpha, beta;
	LIS_SCALAR rho, rho_old;
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
	amr     = solver->work[3];
	p       = solver->work[4];
	ap      = solver->work[5];
	t       = solver->work[6];
	mt      = solver->work[7];
	y       = solver->work[8];
	u       = solver->work[9];
	z       = solver->work[10];
	au      = solver->work[11];


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	time = lis_wtime();
	lis_psolve(solver, r, mr);
	ptime += lis_wtime()-time;
	lis_matvec(A,mr,amr);
	lis_vector_dot(rtld,r,&rho_old);
	lis_vector_copy(amr,ap);
	lis_vector_copy(mr,p);
	beta = 0.0;

	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/* tmpdot[0] = <rtld,ap> */
		/* alpha = rho_old / tmpdot[0] */
		lis_vector_dot(rtld,ap,&tmpdot[0]);
		alpha = rho_old / tmpdot[0];


		/* tmpdot[0] = <y,y>           */
		/* tmpdot[1] = <amr,r>         */
		/* tmpdot[2] = <y,r>           */
		/* tmpdot[3] = <amr,y>         */
		/* tmpdot[4] = <amr,amr>       */
		lis_vector_dot(y,y,&tmpdot[0]);
		lis_vector_dot(amr,r,&tmpdot[1]);
		lis_vector_dot(y,r,&tmpdot[2]);
		lis_vector_dot(amr,y,&tmpdot[3]);
		lis_vector_dot(amr,amr,&tmpdot[4]);
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

		/* t = qsi*ap + eta*y */
		lis_vector_copy(y,t);
		lis_vector_scale(eta,t);
		lis_vector_axpy(qsi,ap,t);

		/* mt  = M^-1 * t */
		time = lis_wtime();
		lis_psolve(solver, t, mt);
		ptime += lis_wtime()-time;

		/* u    = mt + eta*beta*u */
		/* au = A * u             */
		lis_vector_xpay(mt,eta*beta,u);
		lis_matvec(A,u,au);

		/* z = qsi*mr + eta*z - alpha*u */
		lis_vector_scale(eta,z);
		lis_vector_axpy(qsi,mr,z);
		lis_vector_axpy(-alpha,u,z);

		/* y = qsi*amr + eta*y - alpha*au */
		lis_vector_scale(eta,y);
		lis_vector_axpy(qsi,amr,y);
		lis_vector_axpy(-alpha,au,y);

		/* x = x + alpha*p + z */
		lis_vector_axpy(alpha,p,x);
		lis_vector_axpy(1.0,z,x);
		
		/* r = r - alpha*ap - y */
		lis_vector_axpy(-alpha,ap,r);
		lis_vector_axpy(-1.0,y,r);
		
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

		/* rho = <rtld,r> */
		lis_vector_dot(rtld,r,&rho);
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

		/* mr  = M^-1 * r */
		/* amr = A * mr   */
		time = lis_wtime();
		lis_psolve(solver, r, mr);
		ptime += lis_wtime()-time;
		lis_matvec(A,mr,amr);

		/* p  = mr + beta*(p - u)    */
		/* ap = amr + beta*(ap - au) */
		lis_vector_axpy(-1.0,u,p);
		lis_vector_xpay(mr,beta,p);
		lis_vector_axpy(-1.0,au,ap);
		lis_vector_xpay(amr,beta,ap);

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
#define __FUNC__ "lis_bicgsafe_quad"
LIS_INT lis_bicgsafe_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, rhat, p, ptld;
	LIS_VECTOR t, ttld;
	LIS_VECTOR y, v, u, utld, z;
	LIS_QUAD_PTR alpha, beta, rho, rho_old;
	LIS_QUAD_PTR qsi, eta;
	LIS_QUAD_PTR tmp, tmpdot[5],one;
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
	rhat    = solver->work[2];
	p       = solver->work[3];
	ptld    = solver->work[4];
	t       = solver->work[5];
	ttld    = solver->work[6];
	y       = solver->work[7];
	v       = solver->work[8];
	u       = solver->work[9];
	z       = solver->work[10];
	utld    = solver->work[11];

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

	rho_old.hi[0] = 1.0;
	rho_old.lo[0] = 0.0;
	alpha.hi[0] = 1.0;
	alpha.lo[0] = 0.0;
	qsi.hi[0] = 1.0;
	qsi.lo[0] = 0.0;
	one.hi[0] = -1.0;
	one.lo[0] = 0.0;


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_allex_nm(0.0,p);
	lis_vector_set_allex_nm(0.0,u);
	lis_vector_set_allex_nm(0.0,ptld);
	lis_vector_set_allex_nm(0.0,utld);
	
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

		/* beta = (rho / rho_old) * (alpha / qsi) */
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)alpha.hi,(LIS_QUAD *)qsi.hi);
		lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,(LIS_QUAD *)tmp.hi);

		/* rhat = M^-1 * r */
		/* v    = A * rhat */
		time = lis_wtime();
		lis_psolve(solver, r, rhat);
		ptime += lis_wtime()-time;
		lis_matvec(A,rhat,v);

		/* p = rhat + beta*(p - u) */
		lis_vector_axpyex_mmm(one,u,p);
		lis_vector_xpayex_mmm(rhat,beta,p);
		
		/* ptld = v + beta*(ptld - utld) */
		lis_vector_axpyex_mmm(one,utld,ptld);
		lis_vector_xpayex_mmm(v,beta,ptld);

		/* tmpdot[0] = <rtld,ptld> */
		lis_vector_dotex_mmm(rtld,ptld,&tmpdot[0]);
		/* test breakdown */
		/* */
		
		/* alpha = rho / tmpdot[0] */
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)tmpdot[0].hi);


		/* tmpdot[0] = <y,y>       */
		/* tmpdot[1] = <v,r>       */
		/* tmpdot[2] = <y,r>       */
		/* tmpdot[3] = <v,y>       */
		/* tmpdot[4] = <v,v>       */
		lis_vector_dotex_mmm(y,y,&tmpdot[0]);
		lis_vector_dotex_mmm(v,r,&tmpdot[1]);
		lis_vector_dotex_mmm(y,r,&tmpdot[2]);
		lis_vector_dotex_mmm(v,y,&tmpdot[3]);
		lis_vector_dotex_mmm(v,v,&tmpdot[4]);
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

		/* t = qsi*ptld + eta*y */
		lis_vector_copyex_mm(y,t);
		lis_vector_scaleex_mm(eta,t);
		lis_vector_axpyex_mmm(qsi,ptld,t);

		/* ttld  = M^-1 * t */
		time = lis_wtime();
		lis_psolve(solver, t, ttld);
		ptime += lis_wtime()-time;

		/* u    = ttld + eta*beta*u */
		/* utld = A * u             */
		lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)eta.hi,(LIS_QUAD *)beta.hi);
		lis_vector_xpayex_mmm(ttld,tmp,u);
		lis_matvec(A,u,utld);

		/* z = qsi*rhat + eta*z - alpha*u */
		lis_vector_scaleex_mm(eta,z);
		lis_vector_axpyex_mmm(qsi,rhat,z);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,u,z);

		/* y = qsi*v + eta*y - alpha*utld */
		lis_vector_scaleex_mm(eta,y);
		lis_vector_axpyex_mmm(qsi,v,y);
		lis_vector_axpyex_mmm(alpha,utld,y);
		lis_quad_minus((LIS_QUAD *)alpha.hi);

		/* x = x + alpha*p + z */
		lis_vector_axpyex_mmm(alpha,p,x);
		lis_quad_minus((LIS_QUAD *)one.hi);
		lis_vector_axpyex_mmm(one,z,x);
		lis_quad_minus((LIS_QUAD *)one.hi);
		
		/* r = r - alpha*ptld - y */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,ptld,r);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(one,y,r);
		
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
#define __FUNC__ "lis_bicgsafe_switch"
LIS_INT lis_bicgsafe_switch(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, rhat, p, ptld, phat;
	LIS_VECTOR t, ttld, that, t0, t0hat;
	LIS_VECTOR y, w, u, z;
	LIS_QUAD_PTR alpha, beta, rho, rho_old;
	LIS_QUAD_PTR qsi, eta, one;
	LIS_QUAD_PTR tmp, tmpdot[5];
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
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	tol      = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	tol2     = solver->params[LIS_PARAMS_SWITCH_RESID-LIS_OPTIONS_LEN];
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

	rho_old.hi[0] = 1.0;
	rho_old.lo[0] = 0.0;
	alpha.hi[0] = 1.0;
	alpha.lo[0] = 0.0;
	qsi.hi[0] = 1.0;
	qsi.lo[0] = 0.0;
	one.hi[0] = -1.0;
	one.lo[0] = 0.0;


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol2     = solver->tol_switch;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_allex_nm(0.0, ttld);
	lis_vector_set_allex_nm(0.0, ptld);
	lis_vector_set_allex_nm(0.0, p);
	lis_vector_set_allex_nm(0.0, u);
	lis_vector_set_allex_nm(0.0, t);
	lis_vector_set_allex_nm(0.0, t0);

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

		/* beta = (rho / rho_old) * (alpha / qsi) */
		beta.hi[0] = (rho.hi[0] / rho_old.hi[0]) * (alpha.hi[0] / qsi.hi[0]);

		/* w = ttld + beta*ptld */
		lis_vector_axpyz(beta.hi[0],ptld,ttld,w);

		/* rhat = M^-1 * r */
		time = lis_wtime();
		lis_psolve(solver, r, rhat);
		ptime += lis_wtime()-time;

		/* p = rhat + beta*(p - u) */
		lis_vector_axpy(-1,u,p);
		lis_vector_xpay(rhat,beta.hi[0],p);
		
		/* ptld = A * p */
		lis_matvec(A,p,ptld);

		/* tmpdot[0] = <rtld,ptld> */
		lis_vector_dot(rtld,ptld,&tmpdot[0].hi[0]);
		/* test breakdown */
		/* */
		
		/* alpha = rho / tmpdot[0] */
		alpha.hi[0] = rho.hi[0] / tmpdot[0].hi[0];

		/* y = t - r + alpha*(-w + ptld) */
		lis_vector_axpyz(-1,w,ptld,y);
		lis_vector_xpay(t,alpha.hi[0],y);
		lis_vector_axpy(-1,r,y);

		/* t = r - alpha*ptld */
		lis_vector_axpyz(-alpha.hi[0],ptld,r,t);

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
		lis_vector_dot(y,y,&tmpdot[0].hi[0]);
		lis_vector_dot(ttld,t,&tmpdot[1].hi[0]);
		lis_vector_dot(y,t,&tmpdot[2].hi[0]);
		lis_vector_dot(ttld,y,&tmpdot[3].hi[0]);
		lis_vector_dot(ttld,ttld,&tmpdot[4].hi[0]);
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

		/* u = qsi*phat + eta*(t0hat - rhat + beta*u) */
		lis_vector_xpay(t0hat,beta.hi[0],u);
		lis_vector_axpy(-1,rhat,u);
		lis_vector_scale(eta.hi[0],u);
		lis_vector_axpy(qsi.hi[0],phat,u);

		/* z = qsi*rhat + eta*z - alpha*u */
		lis_vector_scale(eta.hi[0],z);
		lis_vector_axpy(qsi.hi[0],rhat,z);
		lis_vector_axpy(-alpha.hi[0],u,z);

		/* x = x + alpha*p + z */
		lis_vector_axpy(alpha.hi[0],p,x);
		lis_vector_axpy(1,z,x);
		
		/* r = t - eta*y - qsi*ttld */
		lis_vector_axpyz(-eta.hi[0],y,t,r);
		lis_vector_axpy(-qsi.hi[0],ttld,r);
		
		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);
		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
			if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
		}

		if( tol2 >= nrm2 )
		{
			solver->iter       = iter;
			solver->iter2      = iter;
			solver->ptime      = ptime;
			break;
		}

		lis_vector_copy(t,t0);
		rho_old.hi[0] = rho.hi[0];
	}

	r->precision = LIS_PRECISION_QUAD;
	p->precision = LIS_PRECISION_QUAD;
	t->precision = LIS_PRECISION_QUAD;
	t0->precision = LIS_PRECISION_QUAD;
	ptld->precision = LIS_PRECISION_QUAD;
	that->precision = LIS_PRECISION_QUAD;

	solver->options[LIS_OPTIONS_INITGUESS_ZEROS] = LIS_FALSE;
	lis_vector_copyex_mn(x,solver->xx);

	rho_old.hi[0] = 1.0;
	alpha.hi[0] = 1.0;
	qsi.hi[0] = 1.0;
	one.hi[0] = -1.0;

	/* Initial Residual */
	lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2);
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_vector_set_allex_nm(0.0, ttld);
	lis_vector_set_allex_nm(0.0, ptld);
	lis_vector_set_allex_nm(0.0, p);
	lis_vector_set_allex_nm(0.0, u);
	lis_vector_set_allex_nm(0.0, t);
	lis_vector_set_allex_nm(0.0, t0);

	for( iter2=iter+1; iter2<=maxiter; iter2++ )
	{
		/* rho = <rtld,r> */
		lis_vector_dotex_mmm(rtld,r,&rho);

		/* test breakdown */
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter2;
			solver->iter2     = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta = (rho / rho_old) * (alpha / qsi) */
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)alpha.hi,(LIS_QUAD *)qsi.hi);
		lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,(LIS_QUAD *)tmp.hi);

		/* w = ttld + beta*ptld */
		lis_vector_axpyzex_mmmm(beta,ptld,ttld,w);

		/* rhat = M^-1 * r */
		time = lis_wtime();
		lis_psolve(solver, r, rhat);
		ptime += lis_wtime()-time;

		/* p = rhat + beta*(p - u) */
		lis_vector_axpyex_mmm(one,u,p);
		lis_vector_xpayex_mmm(rhat,beta,p);
		
		/* ptld = A * p */
		lis_matvec(A,p,ptld);

		/* tmpdot[0] = <rtld,ptld> */
		lis_vector_dotex_mmm(rtld,ptld,&tmpdot[0]);
		/* test breakdown */
		/* */
		
		/* alpha = rho / tmpdot[0] */
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)tmpdot[0].hi);

		/* y = t - r + alpha*(-w + ptld) */
		lis_vector_axpyzex_mmmm(one,w,ptld,y);
		lis_vector_xpayex_mmm(t,alpha,y);
		lis_vector_axpyex_mmm(one,r,y);

		/* t = r - alpha*ptld */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyzex_mmmm(alpha,ptld,r,t);

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
		lis_vector_dotex_mmm(y,y,&tmpdot[0]);
		lis_vector_dotex_mmm(ttld,t,&tmpdot[1]);
		lis_vector_dotex_mmm(y,t,&tmpdot[2]);
		lis_vector_dotex_mmm(ttld,y,&tmpdot[3]);
		lis_vector_dotex_mmm(ttld,ttld,&tmpdot[4]);
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

		/* u = qsi*phat + eta*(t0hat - rhat + beta*u) */
		lis_vector_xpayex_mmm(t0hat,beta,u);
		lis_vector_axpyex_mmm(one,rhat,u);
		lis_vector_scaleex_mm(eta,u);
		lis_vector_axpyex_mmm(qsi,phat,u);

		/* z = qsi*rhat + eta*z - alpha*u */
		lis_vector_scaleex_mm(eta,z);
		lis_vector_axpyex_mmm(qsi,rhat,z);
		lis_vector_axpyex_mmm(alpha,u,z);

		/* x = x + alpha*p + z */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_quad_minus((LIS_QUAD *)one.hi);
		lis_vector_axpyex_mmm(alpha,p,x);
		lis_vector_axpyex_mmm(one,z,x);
		lis_quad_minus((LIS_QUAD *)one.hi);
		
		/* r = t - eta*y - qsi*ttld */
		lis_quad_minus((LIS_QUAD *)eta.hi);
		lis_quad_minus((LIS_QUAD *)qsi.hi);
		lis_vector_axpyzex_mmmm(eta,y,t,r);
		lis_vector_axpyex_mmm(qsi,ttld,r);
		lis_quad_minus((LIS_QUAD *)eta.hi);
		lis_quad_minus((LIS_QUAD *)qsi.hi);
		
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

		lis_vector_copyex_mm(t,t0);
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

/********************************************
 * Preconditioned BiConjugate Residual Safe *
 ********************************************
 r(0)    = b - Ax(0)
 rtld(0) = conj(r(0)) or random
 artld(0)= A^H * rtld(0)
 mr(0)   = M^-1 * r(0)
 amr(0)  = A * mr(0)
 p(0)    = mr(0)
 ap(0)   = amr(0)
 rho(0)  = <rtld,amr(0)>
 ********************************************
 for k=1,2,...
   map(k-1)= M^-1 * ap(k-1)
   tmpdot0 = <map(k-1),artld(0)>
   alpha   = rho(k-1) / tmpdot0
   tmpdot0 = <y(k-1),y(k-1)>
   tmpdot1 = <amr(k-1),r(k-1)>
   tmpdot2 = <y(k-1),r(k-1)>
   tmpdot3 = <amr(k-1),y(k-1)>
   tmpdot4 = <amr(k-1),amr(k-1)>
   tmp     = tmpdot4*tmpdot0-tmpdot3*tmpdot3
   qsi     = (tmpdot0*tmpdot1-tmpdot2*tmpdot3) / tmp
   eta     = (tmpdot4*tmpdot2-tmpdot3*tmpdot1) / tmp
   u(k-1)  = qsi*map(k-1) + eta*my(k) + eta*beta*u(k-2)
   au(k-1) = A * u(k-1)
   z(k-1)  = qsi*mr(k-1) + eta*z(k-2) - alpha*u(k-1)
   y(k)    = qsi*amr(k-1) + eta*y(k-1) - alpha*au(k-1)
   my(k)   = M^-1 * y(k)
   x(k)    = x(k-1) + alpha*p(k-1) + z(k-1)
   r(k)    = r(k-1) - alpha*ap(k-1) - y(k)
   mr(k)   = mr(k-1) - alpha*map(k-1) - my(k)
   amr(k)  = A * mr(k)
   rho(k)  = <rtld,amr(k)>
   beta    = (rho(k) / rho(k-1)) * (alpha / qsi)
   p(k)    = mr(k) + beta*(p(k-1) - u(k-1))
   ap(k-1) = amr(k) + beta*(ap(k-1) - au(k-1))
 ********************************************/
#undef NWORK
#define NWORK 13
#undef __FUNC__
#define __FUNC__ "lis_bicrsafe_check_params"
LIS_INT lis_bicrsafe_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_bicrsafe_malloc_work"
LIS_INT lis_bicrsafe_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_bicgsafe_malloc_work::work" );
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
#define __FUNC__ "lis_bicrsafe"
LIS_INT lis_bicrsafe(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, artld, mr, amr, p, ap, map;
	LIS_VECTOR y, my, u, au, z;
	LIS_SCALAR alpha, beta;
	LIS_SCALAR rho, rho_old;
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
	amr     = solver->work[3];
	p       = solver->work[4];
	ap      = solver->work[5];
	map     = solver->work[6];
	my      = solver->work[7];
	y       = solver->work[8];
	u       = solver->work[9];
	z       = solver->work[10];
	au      = solver->work[11];
	artld   = solver->work[12];


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r,rtld);

	lis_matvech(A,rtld,artld);
	time = lis_wtime();
	lis_psolve(solver, r, mr);
	ptime += lis_wtime()-time;
	lis_matvec(A,mr,amr);
	lis_vector_dot(rtld,amr,&rho_old);
	lis_vector_copy(amr,ap);
	lis_vector_copy(mr,p);
	beta = 0.0;

	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/* map  = M^-1 * ap */
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;

		/* tmpdot[0] = <artld,map> */
		/* alpha = rho_old / tmpdot[0] */
		lis_vector_dot(artld,map,&tmpdot[0]);
		alpha = rho_old / tmpdot[0];


		/* tmpdot[0] = <y,y>           */
		/* tmpdot[1] = <amr,r>         */
		/* tmpdot[2] = <y,r>           */
		/* tmpdot[3] = <amr,y>         */
		/* tmpdot[4] = <amr,amr>       */
		lis_vector_dot(y,y,&tmpdot[0]);
		lis_vector_dot(amr,r,&tmpdot[1]);
		lis_vector_dot(y,r,&tmpdot[2]);
		lis_vector_dot(amr,y,&tmpdot[3]);
		lis_vector_dot(amr,amr,&tmpdot[4]);
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

		/* u    = qsi*map + eta*my + eta*beta*u */
		/* au   = A * u                         */
		lis_vector_scale(eta*beta,u);
		lis_vector_axpy(qsi,map,u);
		lis_vector_axpy(eta,my,u);
		lis_matvec(A,u,au);

		/* z = qsi*mr + eta*z - alpha*u */
		lis_vector_scale(eta,z);
		lis_vector_axpy(qsi,mr,z);
		lis_vector_axpy(-alpha,u,z);

		/* y  = qsi*amr + eta*y - alpha*au */
		/* my = M^-1 * y */
		lis_vector_scale(eta,y);
		lis_vector_axpy(qsi,amr,y);
		lis_vector_axpy(-alpha,au,y);
		time = lis_wtime();
		lis_psolve(solver, y, my);
		ptime += lis_wtime()-time;

		/* x = x + alpha*p + z */
		lis_vector_axpy(alpha,p,x);
		lis_vector_axpy(1.0,z,x);
		
		/* r = r - alpha*ap - y */
		lis_vector_axpy(-alpha,ap,r);
		lis_vector_axpy(-1.0,y,r);
		
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

		/* mr  = mr - alpha*map - my */
		/* amr = A * mr              */
		/* rho = <rtld,amr> */
		lis_vector_axpy(-alpha,map,mr);
		lis_vector_axpy(-1.0,my,mr);
		lis_matvec(A,mr,amr);
		lis_vector_dot(rtld,amr,&rho);
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


		/* p  = mr + beta*(p - u)    */
		/* ap = amr + beta*(ap - au) */
		lis_vector_axpy(-1.0,u,p);
		lis_vector_xpay(mr,beta,p);
		lis_vector_axpy(-1.0,au,ap);
		lis_vector_xpay(amr,beta,ap);

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
#define __FUNC__ "lis_bicrsafe_quad"
LIS_INT lis_bicrsafe_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r, rtld, artld, mr, amr, p, ap, map;
	LIS_VECTOR y, my, u, au, z;
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
	amr     = solver->work[3];
	p       = solver->work[4];
	ap      = solver->work[5];
	map     = solver->work[6];
	my      = solver->work[7];
	y       = solver->work[8];
	u       = solver->work[9];
	z       = solver->work[10];
	au      = solver->work[11];
	artld   = solver->work[12];

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

	lis_matvech(A,rtld,artld);
	time = lis_wtime();
	lis_psolve(solver, r, mr);
	ptime += lis_wtime()-time;
	lis_matvec(A,mr,amr);
	lis_vector_dotex_mmm(rtld,amr,&rho_old);
	lis_vector_copyex_mm(amr,ap);
	lis_vector_copyex_mm(mr,p);
	one.hi[0] = -1.0;
	one.lo[0] = 0.0;

	
	for( iter=1; iter<=maxiter; iter++ )
	{
		/* map  = M^-1 * ap */
		time = lis_wtime();
		lis_psolve(solver, ap, map);
		ptime += lis_wtime()-time;

		/* tmpdot[0] = <artld,map> */
		/* alpha = rho_old / tmpdot[0] */
		lis_vector_dotex_mmm(artld,map,&tmpdot[0]);
		lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho_old.hi,(LIS_QUAD *)tmpdot[0].hi);


		/* tmpdot[0] = <y,y>           */
		/* tmpdot[1] = <amr,r>         */
		/* tmpdot[2] = <y,r>           */
		/* tmpdot[3] = <amr,y>         */
		/* tmpdot[4] = <amr,amr>       */
		lis_vector_dotex_mmm(y,y,&tmpdot[0]);
		lis_vector_dotex_mmm(amr,r,&tmpdot[1]);
		lis_vector_dotex_mmm(y,r,&tmpdot[2]);
		lis_vector_dotex_mmm(amr,y,&tmpdot[3]);
		lis_vector_dotex_mmm(amr,amr,&tmpdot[4]);
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

		/* u    = qsi*map + eta*my + eta*beta*u */
		/* au   = A * u                         */
		lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)eta.hi,(LIS_QUAD *)beta.hi);
		lis_vector_scaleex_mm(tmp,u);
		lis_vector_axpyex_mmm(qsi,map,u);
		lis_vector_axpyex_mmm(eta,my,u);
		lis_matvec(A,u,au);

		/* z = qsi*mr + eta*z - alpha*u */
		lis_vector_scaleex_mm(eta,z);
		lis_vector_axpyex_mmm(qsi,mr,z);
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,u,z);

		/* y  = qsi*amr + eta*y - alpha*au */
		/* my = M^-1 * y */
		lis_vector_scaleex_mm(eta,y);
		lis_vector_axpyex_mmm(qsi,amr,y);
		lis_vector_axpyex_mmm(alpha,au,y);
		time = lis_wtime();
		lis_psolve(solver, y, my);
		ptime += lis_wtime()-time;

		/* x = x + alpha*p + z */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_vector_axpyex_mmm(alpha,p,x);
		lis_quad_minus((LIS_QUAD *)one.hi);
		lis_vector_axpyex_mmm(one,z,x);
		
		/* r = r - alpha*ap - y */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_quad_minus((LIS_QUAD *)one.hi);
		lis_vector_axpyex_mmm(alpha,ap,r);
		lis_vector_axpyex_mmm(one,y,r);
		
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

		/* mr  = mr - alpha*map - my */
		/* amr = A * mr              */
		/* rho = <rtld,amr> */
		lis_vector_axpyex_mmm(alpha,map,mr);
		lis_vector_axpyex_mmm(one,my,mr);
		lis_matvec(A,mr,amr);
		lis_vector_dotex_mmm(rtld,amr,&rho);
		if( rho.hi[0]==0.0 && rho.lo[0]==0.0 )
		{
			solver->retcode   = LIS_BREAKDOWN;
			solver->iter      = iter;
			solver->resid     = nrm2;
			LIS_DEBUG_FUNC_OUT;
			return LIS_BREAKDOWN;
		}

		/* beta = (rho / rho_old) * (alpha / qsi) */
		lis_quad_minus((LIS_QUAD *)alpha.hi);
		lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho.hi,(LIS_QUAD *)rho_old.hi);
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)alpha.hi,(LIS_QUAD *)qsi.hi);
		lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,(LIS_QUAD *)tmp.hi);


		/* p  = mr + beta*(p - u)    */
		/* ap = amr + beta*(ap - au) */
		lis_vector_axpyex_mmm(one,u,p);
		lis_vector_xpayex_mmm(mr,beta,p);
		lis_vector_axpyex_mmm(one,au,ap);
		lis_vector_xpayex_mmm(amr,beta,ap);

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
