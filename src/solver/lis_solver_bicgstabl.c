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

#define NWORK 4
/************************************************
 * lis_bicgstabl_check_params
 * lis_bicgstabl_malloc_work
 * lis_bicgstabl
 ************************************************/
#undef __FUNC__
#define __FUNC__ "lis_bicgstabl_check_params"
LIS_INT lis_bicgstabl_check_params(LIS_SOLVER solver)
{
	LIS_INT	ell;

	LIS_DEBUG_FUNC_IN;

	ell = solver->options[LIS_OPTIONS_ELL];
	if( ell<1 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_ELL(=%D) is less than 1\n",ell);
		return LIS_ERR_ILL_ARG;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_bicgstabl_malloc_work"
LIS_INT lis_bicgstabl_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,l,worklen,err;

	LIS_DEBUG_FUNC_IN;

	l       = solver->options[LIS_OPTIONS_ELL];
	worklen = NWORK + 2*(l+1);
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_bicgstabl_malloc_work::work" );
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
#define __FUNC__ "lis_bicgstabl"
LIS_INT lis_bicgstabl(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR rtld, bp,t,xp, *r,*u;
	LIS_SCALAR *tau, *gamma, *gamma1, *gamma2;
	LIS_SCALAR *sigma;
	LIS_SCALAR alpha, beta, omega, rho0, rho1;
	LIS_REAL rnorm0, rnorm;
	LIS_REAL normx, normr;
	LIS_SCALAR nu;

	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,output,conv;
	double time,ptime;

	LIS_INT l,i,j;
	LIS_INT z_dim;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	l       = solver->options[LIS_OPTIONS_ELL];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;
	z_dim   = l+1;

	rtld    = solver->work[0];
	xp      = solver->work[1];
	bp      = solver->work[2];
	t       = solver->work[3];
	r       = &solver->work[4];
	u       = &solver->work[l+1+4];

	tau     = (LIS_SCALAR *)lis_malloc( sizeof(LIS_SCALAR) * z_dim * (4+l+1),"lis_bicgstabl::tau" );
	if( tau==NULL )
	{
		LIS_SETERR_MEM(sizeof(LIS_SCALAR) * z_dim * (4+l+1));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	gamma   = &tau[z_dim*z_dim];
	gamma1  = &gamma[z_dim];
	gamma2  = &gamma1[z_dim];
	sigma   = &gamma2[z_dim];

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r[0],&bnrm2) )
	{
		lis_free(tau);
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	lis_solver_set_shadowresidual(solver,r[0],rtld);

	lis_vector_copy(r[0],bp);
	lis_vector_copy(x,xp);
	lis_vector_set_all(0.0,u[0]);
	lis_vector_nrm2(r[0],&rnorm0);
	rnorm  = rnorm0;
	normx  = rnorm0;
	normr  = rnorm0;

	alpha = (LIS_SCALAR)0.0;
	omega = (LIS_SCALAR)1.0;
	rho0  = (LIS_SCALAR)1.0;

	iter = 0;
	while( iter<=maxiter )
	{
		/*           */
		/* BiCG PART */
		/*           */

		/* rho0 = -w*rho0 */
		rho0 = -omega*rho0;

		for( j=0; j<l; j++)
		{
			iter++;
			/* rho1 = <rtld,r[j]> */
			lis_vector_dot(rtld,r[j],&rho1);

			/* test breakdown */
			if( rho1==0.0 )
			{
				time = lis_wtime();
				lis_psolve(solver, x, t);
				lis_vector_copy(t, x);
				ptime += lis_wtime()-time;
				lis_vector_axpy(1.0,xp,x);
				solver->retcode   = LIS_BREAKDOWN;
				solver->iter      = iter;
				solver->resid     = nrm2;
				lis_free(tau);
				LIS_DEBUG_FUNC_OUT;
				return LIS_BREAKDOWN;
			}

			/* beta = alpha * (rho1/rho0) */
			/* rho0 = rho1                */
			beta = alpha*(rho1/rho0);
			rho0 = rho1;

			/* u[i] = r[i] - beta*u[i] (i=0,j) */
			for( i=0; i<=j; i++)
			{
				lis_vector_xpay(r[i],-beta,u[i]);
			}

			/* u[j+1] = A    * u[j]   */
			/* u[j+1] = M^-1 * u[j+1] */
			time = lis_wtime();
			lis_psolve(solver, u[j], t);
			ptime += lis_wtime()-time;
			lis_matvec(A,t,u[j+1]);

			/* nu = <rtld, u[j+1]> */
			lis_vector_dot(rtld,u[j+1],&nu);

			/* test breakdown */
			if( nu==0.0 )
			{
				time = lis_wtime();
				lis_psolve(solver, x, t);
				lis_vector_copy(t, x);
				ptime += lis_wtime()-time;
				lis_vector_axpy(1.0,xp,x);
				solver->retcode   = LIS_BREAKDOWN;
				solver->iter      = iter;
				solver->resid     = nrm2;
				lis_free(tau);
				LIS_DEBUG_FUNC_OUT;
				return LIS_BREAKDOWN;
			}

			/* alpha = rho1 / nu */
			alpha = rho1 / nu;

			/* x = x + alpha*u[0] */
			lis_vector_axpy(alpha,u[0],x);

			/* r[i] = r[i] - alpha*u[i+1] (i=0,j) */
			for( i=0; i<=j; i++)
			{
				lis_vector_axpy(-alpha,u[i+1],r[i]);
			}

			/* convergence check */
			lis_solver_get_residual[conv](r[0],solver,&nrm2);
			if( iter%l!=0 )
			{
				if( output )
				{
					if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
					if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
				}
			}

			if( tol >= nrm2 )
			{
				if( output )
				{
					if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
					if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
				}	
				time = lis_wtime();
				lis_psolve(solver, x, t);
				lis_vector_copy(t, x);
				ptime += lis_wtime()-time;
				lis_vector_axpy(1.0,xp,x);
				solver->retcode    = LIS_SUCCESS;
				solver->iter       = iter;
				solver->resid      = nrm2;
				solver->ptime     = ptime;
				lis_free(tau);
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}

			/* r[j+1] = A    * r[j]   */
			/* r[j+1] = M^-1 * r[j+1] */
			time = lis_wtime();
			lis_psolve(solver, r[j], t);
			ptime += lis_wtime()-time;
			lis_matvec(A,t,r[j+1]);

			lis_vector_nrm2(r[0],&rnorm);
			normx = _max(normx,rnorm);
			normr = _max(normr,rnorm);
		}



		/*         */
		/* MR PART */
		/*         */

		for(j=1;j<=l;j++)
		{
			for(i=1;i<=j-1;i++)
			{
				lis_vector_dot(r[j],r[i],&nu);
				nu               = nu / sigma[i];
				tau[i*z_dim + j] = nu;
				lis_vector_axpy(-nu,r[i],r[j]);
			}
			lis_vector_dot(r[j],r[j],&sigma[j]);
			lis_vector_dot(r[0],r[j],&nu);
			gamma1[j] = nu / sigma[j]; 
		}

		gamma[l] = gamma1[l];
		omega    = gamma[l];
		for(j=l-1;j>=1;j--)
		{
			nu = 0.0;
			for(i=j+1;i<=l;i++)
			{
				nu += tau[j*z_dim + i]*gamma[i];
			}
			gamma[j] = gamma1[j] - nu;
		}
		for(j=1;j<=l-1;j++)
		{
			nu = 0.0;
			for(i=j+1;i<=l-1;i++)
			{
				nu += tau[j*z_dim + i]*gamma[i+1];
			}
			gamma2[j] = gamma[j+1] + nu;
		}

		/* UPDATE */
		lis_vector_axpy( gamma[1] ,r[0],x);
		lis_vector_axpy(-gamma1[l],r[l],r[0]);
		lis_vector_axpy(-gamma[l] ,u[l],u[0]);
		for(j=1;j<=l-1;j++)
		{
			lis_vector_axpy(-gamma[j] ,u[j],u[0]);
			lis_vector_axpy( gamma2[j],r[j],x);
			lis_vector_axpy(-gamma1[j],r[j],r[0]);
		}
		lis_solver_get_residual[conv](r[0],solver,&nrm2);
		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
			if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
		}

		if( tol >= nrm2 )
		{
			time = lis_wtime();
			lis_psolve(solver, x, t);
			lis_vector_copy(t, x);
			ptime += lis_wtime()-time;
			lis_vector_axpy(1.0,xp,x);
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			lis_free(tau);
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	lis_free(tau);
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_bicgstabl_quad"
LIS_INT lis_bicgstabl_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR rtld, bp,t,xp, *r,*u;
	LIS_QUAD *tau, *gamma, *gamma1, *gamma2, *sigma;
	LIS_QUAD_PTR alpha, beta, omega, rho0, rho1,sigtmp;
	LIS_QUAD_PTR rnorm0, rnorm;
	LIS_QUAD_PTR normx, normr;
	LIS_QUAD_PTR delta, one, zero, onem;
	LIS_QUAD_PTR nu;

	LIS_REAL bnrm2, nrm2, tol;
	LIS_INT iter,maxiter,n,output,conv;
	double time,ptime;

	LIS_INT l,i,j;
	LIS_INT z_dim;


	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	x       = solver->x;
	n       = A->n;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	l       = solver->options[LIS_OPTIONS_ELL];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	ptime   = 0.0;
	z_dim   = l+1;

	rtld    = solver->work[0];
	xp      = solver->work[1];
	bp      = solver->work[2];
	t       = solver->work[3];
	r       = &solver->work[4];
	u       = &solver->work[l+1+4];

	tau     = (LIS_QUAD *)lis_malloc( sizeof(LIS_QUAD) * z_dim * (4+l+1),"lis_bicgstabl_quad::tau" );
	if( tau==NULL )
	{
		LIS_SETERR_MEM(sizeof(LIS_QUAD) * z_dim * (4+l+1));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	gamma   = &tau[z_dim*z_dim];
	gamma1  = &gamma[z_dim];
	gamma2  = &gamma1[z_dim];
	sigma   = &gamma2[z_dim];

	LIS_QUAD_SCALAR_MALLOC(alpha,0,1);
	LIS_QUAD_SCALAR_MALLOC(beta,1,1);
	LIS_QUAD_SCALAR_MALLOC(rho0,2,1);
	LIS_QUAD_SCALAR_MALLOC(rho1,3,1);
	LIS_QUAD_SCALAR_MALLOC(rnorm0,4,1);
	LIS_QUAD_SCALAR_MALLOC(rnorm,5,1);
	LIS_QUAD_SCALAR_MALLOC(normx,6,1);
	LIS_QUAD_SCALAR_MALLOC(normr,7,1);
	LIS_QUAD_SCALAR_MALLOC(nu,9,1);
	LIS_QUAD_SCALAR_MALLOC(delta,10,1);
	LIS_QUAD_SCALAR_MALLOC(one,11,1);
	LIS_QUAD_SCALAR_MALLOC(zero,12,1);
	LIS_QUAD_SCALAR_MALLOC(onem,13,1);
	LIS_QUAD_SCALAR_MALLOC(sigtmp,14,1);
	LIS_QUAD_SCALAR_MALLOC(omega,15,1);

	one.hi[0]     = 1.0;
	one.lo[0]     = 0.0;
	onem.hi[0]    = -1.0;
	onem.lo[0]    = 0.0;
	zero.hi[0]    = 0.0;
	zero.lo[0]    = 0.0;


	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,t,&bnrm2) )
	{
		lis_free(tau);
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	time = lis_wtime();
	lis_psolve(solver,t,r[0]);
	ptime += lis_wtime()-time;
	lis_solver_set_shadowresidual(solver,r[0],rtld);

	lis_vector_copyex_mm(r[0],bp);
	lis_vector_copyex_mm(x,xp);
	lis_vector_set_allex_nm(0.0,u[0]);
	lis_vector_nrm2ex_mm(r[0],&rnorm0);
	rnorm.hi[0]  = rnorm0.hi[0];
	rnorm.lo[0]  = rnorm0.lo[0];
	normx.hi[0]  = rnorm0.hi[0];
	normx.lo[0]  = rnorm0.lo[0];
	normr.hi[0]  = rnorm0.hi[0];
	normr.lo[0]  = rnorm0.lo[0];

	alpha.hi[0] = 0.0;
	alpha.lo[0] = 0.0;
	omega.hi[0] = 1.0;
	omega.lo[0] = 0.0;
	rho0.hi[0]  = 1.0;
	rho0.lo[0]  = 0.0;
	delta.hi[0] = 1.0e-2;
	delta.lo[0] = 0.0;


	iter = 0;
	while( iter<=maxiter )
	{
		/*           */
		/* BiCG PART */
		/*           */

		/* rho0 = -w*rho0 */
		lis_quad_minus((LIS_QUAD *)omega.hi);
		lis_quad_mul((LIS_QUAD *)rho0.hi,(LIS_QUAD *)omega.hi,(LIS_QUAD *)rho0.hi);

		for( j=0; j<l; j++)
		{
			iter++;
			/* rho1 = <rtld,r[j]> */
			lis_vector_dotex_mmm(rtld,r[j],&rho1);

			/* test breakdown */
			if( rho1.hi[0]==0.0 && rho1.lo[0]==0.0 )
			{
				solver->retcode   = LIS_BREAKDOWN;
				solver->iter      = iter;
				solver->resid     = nrm2;
				lis_free(tau);
				LIS_DEBUG_FUNC_OUT;
				return LIS_BREAKDOWN;
			}

			/* beta = alpha * (rho1/rho0) */
			/* rho0 = rho1                */
			lis_quad_div((LIS_QUAD *)beta.hi,(LIS_QUAD *)rho1.hi,(LIS_QUAD *)rho0.hi);
			lis_quad_mul((LIS_QUAD *)beta.hi,(LIS_QUAD *)beta.hi,(LIS_QUAD *)alpha.hi);
			rho0.hi[0] = rho1.hi[0];
			rho0.lo[0] = rho1.lo[0];

			/* u[i] = r[i] - beta*u[i] (i=0,j) */
			lis_quad_minus((LIS_QUAD *)beta.hi);
			for( i=0; i<=j; i++)
			{
				lis_vector_xpayex_mmm(r[i],beta,u[i]);
			}

			/* u[j+1] = A    * u[j]   */
			/* u[j+1] = M^-1 * u[j+1] */
			lis_matvec(A,u[j],t);
			time = lis_wtime();
			lis_psolve(solver, t, u[j+1]);
			ptime += lis_wtime()-time;

			/* nu = <rtld, u[j+1]> */
			lis_vector_dotex_mmm(rtld,u[j+1],&nu);

			/* test breakdown */
			if( nu.hi[0]==0.0 && nu.lo[0]==0.0 )
			{
				solver->retcode   = LIS_BREAKDOWN;
				solver->iter      = iter;
				solver->resid     = nrm2;
				lis_free(tau);
				LIS_DEBUG_FUNC_OUT;
				return LIS_BREAKDOWN;
			}

			/* alpha = rho1 / nu */
			lis_quad_div((LIS_QUAD *)alpha.hi,(LIS_QUAD *)rho1.hi,(LIS_QUAD *)nu.hi);

			/* x = x + alpha*u[0] */
			lis_vector_axpyex_mmm(alpha,u[0],x);

			/* r[i] = r[i] - alpha*u[i+1] (i=0,j) */
			lis_quad_minus((LIS_QUAD *)alpha.hi);
			for( i=0; i<=j; i++)
			{
				lis_vector_axpyex_mmm(alpha,u[i+1],r[i]);
			}
			lis_quad_minus((LIS_QUAD *)alpha.hi);

			/* convergence check */
			lis_solver_get_residual[conv](r[0],solver,&nrm2);
			if( iter%l!=0 )
			{
				if( output )
				{
					if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
					if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
				}
			}

			if( tol > nrm2 )
			{
				if( output )
				{
					if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
					if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
				}	
				lis_vector_axpyex_mmm(one,xp,x);
				solver->retcode    = LIS_SUCCESS;
				solver->iter       = iter;
				solver->resid      = nrm2;
				solver->ptime      = ptime;
				lis_free(tau);
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}

			/* r[j+1] = A    * r[j]   */
			/* r[j+1] = M^-1 * r[j+1] */
			lis_matvec(A,r[j],t);
			time = lis_wtime();
			lis_psolve(solver, t, r[j+1]);
			ptime += lis_wtime()-time;

			lis_vector_nrm2ex_mm(r[0],&rnorm);
			lis_quad_max((LIS_QUAD *)normx.hi,(LIS_QUAD *)rnorm.hi,(LIS_QUAD *)normx.hi);
			lis_quad_max((LIS_QUAD *)normr.hi,(LIS_QUAD *)rnorm.hi,(LIS_QUAD *)normr.hi);
		}



		/*         */
		/* MR PART */
		/*         */

		for(j=1;j<=l;j++)
		{
			for(i=1;i<=j-1;i++)
			{
				lis_vector_dotex_mmm(r[j],r[i],&nu);
				lis_quad_div((LIS_QUAD *)nu.hi,(LIS_QUAD *)nu.hi,&sigma[i]);
				tau[i*z_dim + j].hi = nu.hi[0];
				tau[i*z_dim + j].lo = nu.lo[0];
/*				nu               = nu / sigma[i];
				tau[i*z_dim + j] = nu;*/
				lis_quad_minus((LIS_QUAD *)nu.hi);
				lis_vector_axpyex_mmm(nu,r[i],r[j]);
			}
			lis_vector_dotex_mmm(r[j],r[j],&sigtmp);
			sigma[j].hi = sigtmp.hi[0];
			sigma[j].lo = sigtmp.lo[0];
			lis_vector_dotex_mmm(r[0],r[j],&nu);
			lis_quad_div(&gamma1[j],(LIS_QUAD *)nu.hi,&sigma[j]);
		}

		gamma[l].hi = gamma1[l].hi;
		gamma[l].lo = gamma1[l].lo;
		omega.hi[0] = gamma[l].hi;
		omega.lo[0] = gamma[l].lo;
		for(j=l-1;j>=1;j--)
		{
			nu.hi[0] = 0.0;
			nu.lo[0] = 0.0;
			for(i=j+1;i<=l;i++)
			{
				lis_quad_mul((LIS_QUAD *)sigtmp.hi,&tau[j*z_dim + i],&gamma[i]);
				lis_quad_add((LIS_QUAD *)nu.hi,(LIS_QUAD *)nu.hi,(LIS_QUAD *)sigtmp.hi);
/*				nu += tau[j*z_dim + i]*gamma[i];*/
			}
			lis_quad_sub(&gamma[j],&gamma1[j],(LIS_QUAD *)nu.hi);
/*			gamma[j] = gamma1[j] - nu;*/
		}
		for(j=1;j<=l-1;j++)
		{
			nu.hi[0] = 0.0;
			nu.lo[0] = 0.0;
			for(i=j+1;i<=l-1;i++)
			{
				lis_quad_mul((LIS_QUAD *)sigtmp.hi,&tau[j*z_dim + i],&gamma[i+1]);
				lis_quad_add((LIS_QUAD *)nu.hi,(LIS_QUAD *)nu.hi,(LIS_QUAD *)sigtmp.hi);
/*				nu += tau[j*z_dim + i]*gamma[i+1];*/
			}
			lis_quad_add(&gamma2[j],&gamma[j+1],(LIS_QUAD *)nu.hi);
/*			gamma2[j] = gamma[j+1] + nu;*/
		}

		/* UPDATE */
		sigtmp.hi[0] = gamma[1].hi;
		sigtmp.lo[0] = gamma[1].lo;
		lis_vector_axpyex_mmm(sigtmp,r[0],x);
		sigtmp.hi[0] = -gamma1[l].hi;
		sigtmp.lo[0] = -gamma1[l].lo;
		lis_vector_axpyex_mmm(sigtmp,r[l],r[0]);
		sigtmp.hi[0] = -gamma[l].hi;
		sigtmp.lo[0] = -gamma[l].lo;
		lis_vector_axpyex_mmm(sigtmp,u[l],u[0]);
		for(j=1;j<=l-1;j++)
		{
			sigtmp.hi[0] = -gamma[j].hi;
			sigtmp.lo[0] = -gamma[j].lo;
			lis_vector_axpyex_mmm(sigtmp,u[j],u[0]);
			sigtmp.hi[0] = gamma2[j].hi;
			sigtmp.lo[0] = gamma2[j].lo;
			lis_vector_axpyex_mmm(sigtmp,r[j],x);
			sigtmp.hi[0] = -gamma1[j].hi;
			sigtmp.lo[0] = -gamma1[j].lo;
			lis_vector_axpyex_mmm(sigtmp,r[j],r[0]);
		}

		lis_solver_get_residual[conv](r[0],solver,&nrm2);
		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
			if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
		}

		if( tol > nrm2 )
		{
			lis_vector_axpyex_mmm(one,xp,x);
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			lis_free(tau);
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	lis_free(tau);
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}
#endif
