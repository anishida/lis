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

#define NWORK 4
/************************************************
 * lis_gmres_check_params
 * lis_gmres_malloc_work
 * lis_gmres
 ************************************************/
#undef __FUNC__
#define __FUNC__ "lis_gmres_check_params"
LIS_INT lis_gmres_check_params(LIS_SOLVER solver)
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
#define __FUNC__ "lis_gmres_malloc_work"
LIS_INT lis_gmres_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,restart,worklen,err;

	LIS_DEBUG_FUNC_IN;

	restart = solver->options[LIS_OPTIONS_RESTART];
	worklen = NWORK + (restart+1);
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_gmres_malloc_work::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
		for(i=1;i<worklen;i++)
		{
			err = lis_vector_duplicate(solver->A,&work[i]);
			if( err ) break;
		}
	}
	else
	{
		for(i=1;i<worklen;i++)
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,solver->A,&work[i]);
			if( err ) break;
			memset(work[i]->value_lo,0,solver->A->np*sizeof(LIS_SCALAR));
		}
	}
	if( i<worklen )
	{
		for(j=1;j<i;j++) lis_vector_destroy(work[j]);
		lis_free(work);
		return err;
	}
	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
		lis_vector_create(solver->A->comm,&work[0]);
	}
	else
	{
		lis_vector_createex(LIS_PRECISION_QUAD,solver->A->comm,&work[0]);
	}
	lis_vector_set_size(work[0],restart+1,0);
	solver->worklen = worklen;
	solver->work    = work;


	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_gmres"
LIS_INT lis_gmres(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR b,x;
	LIS_VECTOR r,s,z,*v;
	LIS_SCALAR *h;
	LIS_SCALAR aa,bb,rr,a2,b2;
	LIS_SCALAR t;

	LIS_REAL bnrm2,nrm2,tol;
	LIS_INT iter,maxiter,n,output;
	double time,ptime;

	LIS_REAL rnorm;
	LIS_INT i,j,k,m;
	LIS_INT ii,i1,iiv,i1v,iih,jj;
	LIS_INT h_dim;
	LIS_INT cs,sn;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	b       = solver->b;
	x       = solver->x;
	n       = A->n;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	m       = solver->options[LIS_OPTIONS_RESTART];
	h_dim   = m+1;
	ptime   = 0.0;

	s       = solver->work[0];
	r       = solver->work[1];
	z       = solver->work[2];
	v       = &solver->work[3];

	h       = (LIS_SCALAR *)lis_malloc( sizeof(LIS_SCALAR)*(h_dim+1)*(h_dim+2),"lis_gmres::h" );
	cs      = (m+1)*h_dim;
	sn      = (m+2)*h_dim;

	/* r = M^-1 * (b - A * x) */
	lis_matvec(A,x,z);
	lis_vector_xpay(b,-1.0,z);
	lis_psolve(solver,z,v[0]);
	
	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,v[0],&bnrm2) )
	{
		lis_free(h);
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;


	iter=0;
	while( iter<maxiter )
	{
		/* first column of V */
		/* v = r / ||r||_2 */
		lis_vector_nrm2(v[0],&rnorm);
		lis_vector_scale(1.0/rnorm,v[0]);

		/* s = ||r||_2 e_1 */
		lis_vector_set_all(0,s);
		s->value[0] = rnorm;

		i = 0;
		do
		{
			iter++;
			i++;
			ii  = i-1;
			i1  = i;
			iiv = i-1;
			i1v = i;
			iih = (i-1)*h_dim;


			/* z = M^-1 * v */
			time = lis_wtime();
			lis_psolve(solver,v[iiv],z);
			ptime += lis_wtime()-time;

			/* w = A * z */
			lis_matvec(A,z,v[i1v]);

			for(k=0;k<i;k++)
			{
				/* h[k,i]   = <w,v[k]>          */
				/* w        = w - h[k,i] * v[k] */
				lis_vector_dot(v[i1v],v[k],&t);
				h[k+iih] = t;
				lis_vector_axpy(-t,v[k],v[i1v]);
			}
			/* h[i+1,i] = ||w||          */
			/* v[i+1]   = w / h[i+1,i]   */
			lis_vector_nrm2(v[i1v],(LIS_REAL *)&t);
			h[i1+iih] = t;
			lis_vector_scale(1.0/t,v[i1v]);

			for(k=1;k<=ii;k++)
			{
				jj  =  k-1;
				t   =  h[jj+iih];
				aa  =  h[jj+cs]*t;
				aa +=  h[jj+sn]*h[k+iih];
				bb  = -h[jj+sn]*t;
				bb +=  h[jj+cs]*h[k+iih];
				h[jj+iih] = aa;
				h[k+iih] = bb;
			}
			aa = h[ii+iih];
			bb = h[i1+iih];
			a2 = aa*aa;
			b2 = bb*bb;
			rr = sqrt(a2+b2);
			if( rr==0.0 ) rr=1.0e-17;
			h[ii+cs] = aa/rr;
			h[ii+sn] = bb/rr;
			s->value[i1] = -h[ii+sn]*s->value[ii];
			s->value[ii] =  h[ii+cs]*s->value[ii];

			aa  =  h[ii+cs]*h[ii+iih];
			aa +=  h[ii+sn]*h[i1+iih];
			h[ii+iih] = aa;

			/* convergence check */
			nrm2 = fabs(s->value[i1])*bnrm2;

			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			if( tol >= nrm2 ) break;
		} while( i<m && iter <maxiter );

		/* Solve H * Y = S for upper Hessenberg matrix H */
		s->value[ii] = s->value[ii]/h[ii+iih];
		for(k=1;k<=ii;k++)
		{
			jj = ii-k;
			t  = s->value[jj];
			for(j=jj+1;j<=ii;j++)
			{
				t -= h[jj+j*h_dim]*s->value[j];
			}
			s->value[jj] = t/h[jj+jj*h_dim];
		}
		/* z = z + y * v */
		#ifdef _OPENMP
		#pragma omp parallel for private(k)
		#endif
		for(k=0;k<n;k++)
		{
			z->value[k] = s->value[0]*v[0]->value[k];
		}
		for(j=1;j<=ii;j++)
		{
			lis_vector_axpy(s->value[j],v[j],z);
		}

		/* r = M^-1 * z */
		time = lis_wtime();
		lis_psolve(solver,z,r);
		ptime += lis_wtime()-time;

		/* x = x + r */
		lis_vector_axpy(1,r,x);

		if( tol >= nrm2 )
		{
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			lis_free(h);
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		for(j=1;j<=i;j++)
		{
			jj = i1-j+1;
			s->value[jj-1] = -h[jj-1+sn]*s->value[jj];
			s->value[jj]   =  h[jj-1+cs]*s->value[jj];
		}

		for(j=0;j<=i1;j++)
		{
			t = s->value[j];
			if( j==0 ) t = t-1.0;
			lis_vector_axpy(t,v[j],v[0]);
		}
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter+1;
	solver->resid     = nrm2;
	lis_free(h);
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_gmres_quad"
LIS_INT lis_gmres_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR b,x;
	LIS_VECTOR r,s,z,*v;
	LIS_QUAD *h;
	LIS_QUAD_PTR aa,bb,rr,a2,b2,t,one,tmp;
	LIS_QUAD_PTR rnorm;

	LIS_REAL bnrm2,nrm2,tol;
	LIS_INT iter,maxiter,n,output;
	double time,ptime;

	LIS_INT i,j,k,m;
	LIS_INT ii,i1,iiv,i1v,iih,jj;
	LIS_INT h_dim;
	LIS_INT cs,sn;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	b       = solver->b;
	x       = solver->x;
	n       = A->n;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	m       = solver->options[LIS_OPTIONS_RESTART];
	h_dim   = m+1;
	ptime   = 0.0;

	s       = solver->work[0];
	r       = solver->work[1];
	z       = solver->work[2];
	v       = &solver->work[3];

	LIS_QUAD_SCALAR_MALLOC(aa,0,1);
	LIS_QUAD_SCALAR_MALLOC(bb,1,1);
	LIS_QUAD_SCALAR_MALLOC(rr,2,1);
	LIS_QUAD_SCALAR_MALLOC(a2,3,1);
	LIS_QUAD_SCALAR_MALLOC(b2,4,1);
	LIS_QUAD_SCALAR_MALLOC(t,5,1);
	LIS_QUAD_SCALAR_MALLOC(tmp,6,1);
	LIS_QUAD_SCALAR_MALLOC(one,7,1);
	LIS_QUAD_SCALAR_MALLOC(rnorm,8,1);

	h       = (LIS_QUAD *)lis_malloc( sizeof(LIS_QUAD)*(h_dim+1)*(h_dim+2),"lis_gmres_quad::h" );
	cs      = (m+1)*h_dim;
	sn      = (m+2)*h_dim;
	one.hi[0]   = 1.0;
	one.lo[0]   = 0.0;

	/* r = M^-1 * (b - A * x) */
	lis_matvec(A,x,z);
	lis_vector_xpay(b,-1.0,z);
	lis_psolve(solver,z,v[0]);

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,v[0],&bnrm2) )
	{
		lis_free(h);
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;


	iter=0;
	while( iter<maxiter )
	{
		/* first column of V */
		/* v = r / ||r||_2 */
		lis_vector_nrm2ex_mm(v[0],&rnorm);
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)one.hi,(LIS_QUAD *)rnorm.hi);
		lis_vector_scaleex_mm(tmp,v[0]);

		/* s = ||r||_2 e_1 */
		lis_vector_set_allex_nm(0.0,s);
		s->value[0]    = rnorm.hi[0];
		s->value_lo[0] = rnorm.lo[0];

		i = 0;
		do
		{
			iter++;
			i++;
			ii  = i-1;
			i1  = i;
			iiv = i-1;
			i1v = i;
			iih = (i-1)*h_dim;


			/* z = M^-1 * v */
			time = lis_wtime();
			lis_psolve(solver,v[iiv],z);
			ptime += lis_wtime()-time;

			/* w = A * z */
			lis_matvec(A,z, v[i1v]);

			for(k=0;k<i;k++)
			{
				/* h[k,i]   = <w,v[k]>          */
				/* w        = w - h[k,i] * v[k] */
				lis_vector_dotex_mmm(v[i1v],v[k],&t);
				h[k+iih].hi = t.hi[0];
				h[k+iih].lo = t.lo[0];
				lis_quad_minus((LIS_QUAD *)t.hi);
				lis_vector_axpyex_mmm(t,v[k],v[i1v]);
			}
			/* h[i+1,i] = ||w||          */
			/* v[i+1]   = w / h[i+1,i]   */
			lis_vector_nrm2ex_mm(v[i1v],&t);
			h[i1+iih].hi = t.hi[0];
			h[i1+iih].lo = t.lo[0];
			lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)one.hi,(LIS_QUAD *)t.hi);
			lis_vector_scaleex_mm(tmp,v[i1v]);

			for(k=1;k<=ii;k++)
			{
				jj  = k-1;
				t.hi[0]   =  h[jj+iih].hi;
				t.lo[0]   =  h[jj+iih].lo;
				lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)&h[jj+cs],(LIS_QUAD *)t.hi);
				lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj+sn],(LIS_QUAD *)&h[k+iih]);
				lis_quad_add((LIS_QUAD *)aa.hi,(LIS_QUAD *)aa.hi,(LIS_QUAD *)tmp.hi);
				lis_quad_mul((LIS_QUAD *)bb.hi,(LIS_QUAD *)&h[jj+sn],(LIS_QUAD *)t.hi);
				lis_quad_minus((LIS_QUAD *)bb.hi);
				lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj+cs],(LIS_QUAD *)&h[k+iih]);
				lis_quad_add((LIS_QUAD *)bb.hi,(LIS_QUAD *)bb.hi,(LIS_QUAD *)tmp.hi);
				h[jj+iih].hi = aa.hi[0];
				h[jj+iih].lo = aa.lo[0];
				h[k+iih].hi = bb.hi[0];
				h[k+iih].lo = bb.lo[0];
			}
			aa.hi[0] = h[ii+iih].hi;
			aa.lo[0] = h[ii+iih].lo;
			bb.hi[0] = h[i1+iih].hi;
			bb.lo[0] = h[i1+iih].lo;
			lis_quad_sqr((LIS_QUAD *)a2.hi,(LIS_QUAD *)aa.hi);
			lis_quad_sqr((LIS_QUAD *)b2.hi,(LIS_QUAD *)bb.hi);
			lis_quad_add((LIS_QUAD *)rr.hi,(LIS_QUAD *)a2.hi,(LIS_QUAD *)b2.hi);
			lis_quad_sqrt((LIS_QUAD *)rr.hi,(LIS_QUAD *)rr.hi);
			if( rr.hi[0]==0.0 )
			{
				rr.hi[0]=1.0e-17;
				rr.lo[0]=0.0;
			}
			lis_quad_div((LIS_QUAD *)&h[ii+cs],(LIS_QUAD *)aa.hi,(LIS_QUAD *)rr.hi);
			lis_quad_div((LIS_QUAD *)&h[ii+sn],(LIS_QUAD *)bb.hi,(LIS_QUAD *)rr.hi);
			tmp.hi[0] = s->value[ii];
			tmp.lo[0] = s->value_lo[ii];
			lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)&h[ii+sn],(LIS_QUAD *)tmp.hi);
			lis_quad_mul((LIS_QUAD *)bb.hi,(LIS_QUAD *)&h[ii+cs],(LIS_QUAD *)tmp.hi);
			lis_quad_minus((LIS_QUAD *)aa.hi);
			s->value[i1] = aa.hi[0];
			s->value_lo[i1] = aa.lo[0];
			s->value[ii] = bb.hi[0];
			s->value_lo[ii] = bb.lo[0];

			lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)&h[ii+cs],(LIS_QUAD *)&h[ii+iih]);
			lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[ii+sn],(LIS_QUAD *)&h[i1+iih]);
			lis_quad_add((LIS_QUAD *)aa.hi,(LIS_QUAD *)aa.hi,(LIS_QUAD *)tmp.hi);
			h[ii+iih].hi = aa.hi[0];
			h[ii+iih].lo = aa.lo[0];

			/* convergence check */
			nrm2 = fabs(s->value[i1])*bnrm2;

			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			if( tol >= nrm2 ) break;
		} while( i<m && iter <maxiter );

		/* Solve H * Y = S for upper Hessenberg matrix H */
		tmp.hi[0] = s->value[ii];
		tmp.lo[0] = s->value_lo[ii];
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[ii+iih]);
		s->value[ii] = tmp.hi[0];
		s->value_lo[ii] = tmp.lo[0];
		for(k=1;k<=ii;k++)
		{
			jj = ii-k;
			t.hi[0]  = s->value[jj];
			t.lo[0]  = s->value_lo[jj];
			for(j=jj+1;j<=ii;j++)
			{
				tmp.hi[0] = s->value[j];
				tmp.lo[0] = s->value_lo[j];
				lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj+j*h_dim]);
				lis_quad_sub((LIS_QUAD *)t.hi,(LIS_QUAD *)t.hi,(LIS_QUAD *)tmp.hi);
			}
			lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)t.hi,(LIS_QUAD *)&h[jj+jj*h_dim]);
			s->value[jj] = tmp.hi[0];
			s->value_lo[jj] = tmp.lo[0];
		}
		/* z = z + y * v */
		for(k=0;k<n;k++)
		{
			aa.hi[0] = s->value[0];
			aa.lo[0] = s->value_lo[0];
			bb.hi[0] = v[0]->value[k];
			bb.lo[0] = v[0]->value_lo[k];
			lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)aa.hi,(LIS_QUAD *)bb.hi);
			z->value[k] = tmp.hi[0];
			z->value_lo[k] = tmp.lo[0];
		}
		for(j=1;j<=ii;j++)
		{
			aa.hi[0] = s->value[j];
			aa.lo[0] = s->value_lo[j];
			lis_vector_axpyex_mmm(aa,v[j],z);
		}
		/* r = M^-1 * z */
		time = lis_wtime();
		lis_psolve(solver,z,r);
		ptime += lis_wtime()-time;

		/* x = x + r */
		lis_vector_axpyex_mmm(one,r,x);

		if( tol >= nrm2 )
		{
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->iter2      = 0;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			lis_free(h);
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		for(j=1;j<=i;j++)
		{
			jj = i1-j+1;
			tmp.hi[0] = s->value[jj];
			tmp.lo[0] = s->value_lo[jj];
			lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj-1+sn]);
			lis_quad_mul((LIS_QUAD *)bb.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj-1+cs]);
			lis_quad_minus((LIS_QUAD *)aa.hi);
			s->value[jj-1] = aa.hi[0];
			s->value_lo[jj-1] = aa.lo[0];
			s->value[jj] = bb.hi[0];
			s->value_lo[jj] = bb.lo[0];
		}

		for(j=0;j<=i1;j++)
		{
			t.hi[0] = s->value[j];
			t.lo[0] = s->value_lo[j];
			if( j==0 )
			{
				lis_quad_sub((LIS_QUAD *)t.hi,(LIS_QUAD *)t.hi,(LIS_QUAD *)one.hi);
			}
			lis_vector_axpyex_mmm(t,v[j],v[0]);
		}
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter+1;
	solver->iter2     = 0;
	solver->resid     = nrm2;
	lis_free(h);
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

#undef __FUNC__
#define __FUNC__ "lis_gmres_switch"
LIS_INT lis_gmres_switch(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR b,x;
	LIS_VECTOR r,s,z,*v;
	LIS_QUAD *h;
	LIS_SCALAR *hd;
	LIS_QUAD_PTR aa,bb,rr,a2,b2,t,one,tmp;
	LIS_QUAD_PTR rnorm;
	LIS_REAL bnrm2,nrm2,tol,tol2;
	LIS_INT iter,maxiter,n,output;
	LIS_INT iter2,maxiter2;
	double time,ptime;

	LIS_INT i,j,k,m;
	LIS_INT ii,i1,iiv,i1v,iih,jj;
	LIS_INT h_dim;
	LIS_INT cs,sn;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	b       = solver->b;
	x       = solver->x;
	n       = A->n;
	maxiter  = solver->options[LIS_OPTIONS_MAXITER];
	maxiter2 = solver->options[LIS_OPTIONS_SWITCH_MAXITER];
	output   = solver->options[LIS_OPTIONS_OUTPUT];
	tol      = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	tol2     = solver->params[LIS_PARAMS_SWITCH_RESID-LIS_OPTIONS_LEN];
	m        = solver->options[LIS_OPTIONS_RESTART];
	h_dim    = m+1;
	ptime    = 0.0;

	s       = solver->work[0];
	r       = solver->work[1];
	z       = solver->work[2];
	v       = &solver->work[3];

	LIS_QUAD_SCALAR_MALLOC(aa,0,1);
	LIS_QUAD_SCALAR_MALLOC(bb,1,1);
	LIS_QUAD_SCALAR_MALLOC(rr,2,1);
	LIS_QUAD_SCALAR_MALLOC(a2,3,1);
	LIS_QUAD_SCALAR_MALLOC(b2,4,1);
	LIS_QUAD_SCALAR_MALLOC(t,5,1);
	LIS_QUAD_SCALAR_MALLOC(tmp,6,1);
	LIS_QUAD_SCALAR_MALLOC(one,7,1);
	LIS_QUAD_SCALAR_MALLOC(rnorm,8,1);

	h       = (LIS_QUAD *)lis_malloc( sizeof(LIS_QUAD)*(h_dim+1)*(h_dim+2),"lis_gmres_switch::h" );
	hd      = (LIS_SCALAR *)h;
	cs      = (m+1)*h_dim;
	sn      = (m+2)*h_dim;
	one.hi[0]   = 1.0;
	one.lo[0]   = 0.0;

	z->precision = LIS_PRECISION_DEFAULT;

	/* r = M^-1 * (b - A * x) */
	lis_matvec(A,x,z);
	lis_vector_xpay(b,-1.0,z);
	lis_psolve(solver,z,v[0]);

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,v[0],&bnrm2) )
	{
		lis_free(h);
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol2     = solver->tol_switch;


	iter=0;
	while( iter<maxiter2 )
	{
		/* first column of V */
		/* v = r / ||r||_2 */
		lis_vector_nrm2(v[0],&rnorm.hi[0]);
		lis_vector_scale(1.0/rnorm.hi[0],v[0]);

		/* s = ||r||_2 e_1 */
		lis_vector_set_all(0,s);
		s->value[0] = rnorm.hi[0];

		i = 0;
		do
		{
			iter++;
			i++;
			ii  = i-1;
			i1  = i;
			iiv = i-1;
			i1v = i;
			iih = (i-1)*h_dim;


			/* z = M^-1 * v */
			time = lis_wtime();
			lis_psolve(solver,v[iiv],z);
			ptime += lis_wtime()-time;

			/* w = A * z */
			lis_matvec(A,z, v[i1v]);

			for(k=0;k<i;k++)
			{
				/* h[k,i]   = <w,v[k]>          */
				/* w        = w - h[k,i] * v[k] */
				lis_vector_dot(v[i1v],v[k],&t.hi[0]);
				hd[k+iih] = t.hi[0];
				lis_vector_axpy(-t.hi[0],v[k],v[i1v]);
			}
			/* h[i+1,i] = ||w||          */
			/* v[i+1]   = w / h[i+1,i]   */
			lis_vector_nrm2(v[i1v],&t.hi[0]);
			hd[i1+iih] = t.hi[0];
			lis_vector_scale(1.0/t.hi[0],v[i1v]);

			for(k=1;k<=ii;k++)
			{
				jj        = k-1;
				t.hi[0]   =  hd[jj+iih];
				aa.hi[0]  =  hd[jj+cs]*t.hi[0];
				aa.hi[0] +=  hd[jj+sn]*hd[k+iih];
				bb.hi[0]  = -hd[jj+sn]*t.hi[0];
				bb.hi[0] +=  hd[jj+cs]*hd[k+iih];
				hd[jj+iih] = aa.hi[0];
				hd[k+iih] = bb.hi[0];
			}
			aa.hi[0] = hd[ii+iih];
			bb.hi[0] = hd[i1+iih];
			a2.hi[0] = aa.hi[0]*aa.hi[0];
			b2.hi[0] = bb.hi[0]*bb.hi[0];
			rr.hi[0] = sqrt(a2.hi[0]+b2.hi[0]);
			if( rr.hi[0]==0.0 ) rr.hi[0]=1.0e-17;
			hd[ii+cs] = aa.hi[0]/rr.hi[0];
			hd[ii+sn] = bb.hi[0]/rr.hi[0];
			s->value[i1] = -hd[ii+sn]*s->value[ii];
			s->value[ii] =  hd[ii+cs]*s->value[ii];

			aa.hi[0]  =  hd[ii+cs]*hd[ii+iih];
			aa.hi[0] +=  hd[ii+sn]*hd[i1+iih];
			hd[ii+iih] = aa.hi[0];

			/* convergence check */
			nrm2 = fabs(s->value[i1])*bnrm2;

			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			if( tol2 >= nrm2 ) break;
		} while( i<m && iter <maxiter2 );

		/* Solve H * Y = S for upper Hessenberg matrix H */
		s->value[ii] = s->value[ii]/hd[ii+iih];
		for(k=1;k<=ii;k++)
		{
			jj = ii-k;
			t.hi[0]  = s->value[jj];
			for(j=jj+1;j<=ii;j++)
			{
				t.hi[0] -= hd[jj+j*h_dim]*s->value[j];
			}
			s->value[jj] = t.hi[0]/hd[jj+jj*h_dim];
		}
		/* z = z + y * v */
		for(k=0;k<n;k++)
		{
			z->value[k] = s->value[0]*v[0]->value[k];
		}
		for(j=1;j<=ii;j++)
		{
			lis_vector_axpy(s->value[j],v[j],z);
		}
		/* r = M^-1 * z */
		time = lis_wtime();
		lis_psolve(solver,z,r);
		ptime += lis_wtime()-time;

		/* x = x + r */
		lis_vector_axpy(1,r,x);

		if( tol2 >= nrm2 )
		{
			solver->iter       = iter;
			solver->iter2      = iter;
			solver->ptime      = ptime;
			break;
		}

		for(j=1;j<=i;j++)
		{
			jj = i1-j+1;
			s->value[jj-1] = -hd[jj-1+sn]*s->value[jj];
			s->value[jj]   =  hd[jj-1+cs]*s->value[jj];
		}

		for(j=0;j<=i1;j++)
		{
			t.hi[0] = s->value[j];
			if( j==0 ) t.hi[0] = t.hi[0]-1.0;
			lis_vector_axpy(t.hi[0],v[j],v[0]);
		}
	}

	/* Initial Residual */
	z->precision = LIS_PRECISION_QUAD;

	solver->options[LIS_OPTIONS_INITGUESS_ZEROS] = LIS_FALSE;
	lis_vector_copyex_mn(x,solver->xx);

	lis_solver_get_initial_residual(solver,NULL,NULL,v[0],&bnrm2);
	tol     = solver->tol;


	iter2=iter;
	while( iter2<maxiter )
	{
		/* first column of V */
		/* v = r / ||r||_2 */
		lis_vector_nrm2ex_mm(v[0],&rnorm);
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)one.hi,(LIS_QUAD *)rnorm.hi);
		lis_vector_scaleex_mm(tmp,v[0]);

		/* s = ||r||_2 e_1 */
		lis_vector_set_allex_nm(0.0,s);
		s->value[0]    = rnorm.hi[0];
		s->value_lo[0] = rnorm.lo[0];

		i = 0;
		do
		{
			iter2++;
			i++;
			ii  = i-1;
			i1  = i;
			iiv = i-1;
			i1v = i;
			iih = (i-1)*h_dim;


			/* z = M^-1 * v */
			time = lis_wtime();
			lis_psolve(solver,v[iiv],z);
			ptime += lis_wtime()-time;

			/* w = A * z */
			lis_matvec(A,z, v[i1v]);

			for(k=0;k<i;k++)
			{
				/* h[k,i]   = <w,v[k]>          */
				/* w        = w - h[k,i] * v[k] */
				lis_vector_dotex_mmm(v[i1v],v[k],&t);
				h[k+iih].hi = t.hi[0];
				h[k+iih].lo = t.lo[0];
				lis_quad_minus((LIS_QUAD *)t.hi);
				lis_vector_axpyex_mmm(t,v[k],v[i1v]);
			}
			/* h[i+1,i] = ||w||          */
			/* v[i+1]   = w / h[i+1,i]   */
			lis_vector_nrm2ex_mm(v[i1v],&t);
			h[i1+iih].hi = t.hi[0];
			h[i1+iih].lo = t.lo[0];
			lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)one.hi,(LIS_QUAD *)t.hi);
			lis_vector_scaleex_mm(tmp,v[i1v]);

			for(k=1;k<=ii;k++)
			{
				jj  = k-1;
				t.hi[0]   =  h[jj+iih].hi;
				t.lo[0]   =  h[jj+iih].lo;
				lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)&h[jj+cs],(LIS_QUAD *)t.hi);
				lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj+sn],(LIS_QUAD *)&h[k+iih]);
				lis_quad_add((LIS_QUAD *)aa.hi,(LIS_QUAD *)aa.hi,(LIS_QUAD *)tmp.hi);
				lis_quad_mul((LIS_QUAD *)bb.hi,(LIS_QUAD *)&h[jj+sn],(LIS_QUAD *)t.hi);
				lis_quad_minus((LIS_QUAD *)bb.hi);
				lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj+cs],(LIS_QUAD *)&h[k+iih]);
				lis_quad_add((LIS_QUAD *)bb.hi,(LIS_QUAD *)bb.hi,(LIS_QUAD *)tmp.hi);
				h[jj+iih].hi = aa.hi[0];
				h[jj+iih].lo = aa.lo[0];
				h[k+iih].hi = bb.hi[0];
				h[k+iih].lo = bb.lo[0];
			}
			aa.hi[0] = h[ii+iih].hi;
			aa.lo[0] = h[ii+iih].lo;
			bb.hi[0] = h[i1+iih].hi;
			bb.lo[0] = h[i1+iih].lo;
			lis_quad_sqr((LIS_QUAD *)a2.hi,(LIS_QUAD *)aa.hi);
			lis_quad_sqr((LIS_QUAD *)b2.hi,(LIS_QUAD *)bb.hi);
			lis_quad_add((LIS_QUAD *)rr.hi,(LIS_QUAD *)a2.hi,(LIS_QUAD *)b2.hi);
			lis_quad_sqrt((LIS_QUAD *)rr.hi,(LIS_QUAD *)rr.hi);
			lis_quad_div((LIS_QUAD *)&h[ii+cs],(LIS_QUAD *)aa.hi,(LIS_QUAD *)rr.hi);
			lis_quad_div((LIS_QUAD *)&h[ii+sn],(LIS_QUAD *)bb.hi,(LIS_QUAD *)rr.hi);
			tmp.hi[0] = s->value[ii];
			tmp.lo[0] = s->value_lo[ii];
			lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)&h[ii+sn],(LIS_QUAD *)tmp.hi);
			lis_quad_mul((LIS_QUAD *)bb.hi,(LIS_QUAD *)&h[ii+cs],(LIS_QUAD *)tmp.hi);
			lis_quad_minus((LIS_QUAD *)aa.hi);
			s->value[i1] = aa.hi[0];
			s->value_lo[i1] = aa.lo[0];
			s->value[ii] = bb.hi[0];
			s->value_lo[ii] = bb.lo[0];

			lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)&h[ii+cs],(LIS_QUAD *)&h[ii+iih]);
			lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[ii+sn],(LIS_QUAD *)&h[i1+iih]);
			lis_quad_add((LIS_QUAD *)aa.hi,(LIS_QUAD *)aa.hi,(LIS_QUAD *)tmp.hi);
			h[ii+iih].hi = aa.hi[0];
			h[ii+iih].lo = aa.lo[0];

			/* convergence check */
			nrm2 = fabs(s->value[i1])*bnrm2;

			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter2] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			if( tol >= nrm2 ) break;
		} while( i<m && iter2 <maxiter );

		/* Solve H * Y = S for upper Hessenberg matrix H */
		tmp.hi[0] = s->value[ii];
		tmp.lo[0] = s->value_lo[ii];
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[ii+iih]);
		s->value[ii] = tmp.hi[0];
		s->value_lo[ii] = tmp.lo[0];
		for(k=1;k<=ii;k++)
		{
			jj = ii-k;
			t.hi[0]  = s->value[jj];
			t.lo[0]  = s->value_lo[jj];
			for(j=jj+1;j<=ii;j++)
			{
				tmp.hi[0] = s->value[j];
				tmp.lo[0] = s->value_lo[j];
				lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj+j*h_dim]);
				lis_quad_sub((LIS_QUAD *)t.hi,(LIS_QUAD *)t.hi,(LIS_QUAD *)tmp.hi);
			}
			lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)t.hi,(LIS_QUAD *)&h[jj+jj*h_dim]);
			s->value[jj] = tmp.hi[0];
			s->value_lo[jj] = tmp.lo[0];
		}
		/* z = z + y * v */
		for(k=0;k<n;k++)
		{
			aa.hi[0] = s->value[0];
			aa.lo[0] = s->value_lo[0];
			bb.hi[0] = v[0]->value[k];
			bb.lo[0] = v[0]->value_lo[k];
			lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)aa.hi,(LIS_QUAD *)bb.hi);
			z->value[k] = tmp.hi[0];
			z->value_lo[k] = tmp.lo[0];
		}
		for(j=1;j<=ii;j++)
		{
			aa.hi[0] = s->value[j];
			aa.lo[0] = s->value_lo[j];
			lis_vector_axpyex_mmm(aa,v[j],z);
		}
		/* r = M^-1 * z */
		time = lis_wtime();
		lis_psolve(solver,z,r);
		ptime += lis_wtime()-time;

		/* x = x + r */
		lis_vector_axpyex_mmm(one,r,x);

		if( tol >= nrm2 )
		{
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter2;
			solver->iter2      = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			lis_free(h);
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		for(j=1;j<=i;j++)
		{
			jj = i1-j+1;
			tmp.hi[0] = s->value[jj];
			tmp.lo[0] = s->value_lo[jj];
			lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj-1+sn]);
			lis_quad_mul((LIS_QUAD *)bb.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj-1+cs]);
			lis_quad_minus((LIS_QUAD *)aa.hi);
			s->value[jj-1] = aa.hi[0];
			s->value_lo[jj-1] = aa.lo[0];
			s->value[jj] = bb.hi[0];
			s->value_lo[jj] = bb.lo[0];
		}

		for(j=0;j<=i1;j++)
		{
			t.hi[0] = s->value[j];
			t.lo[0] = s->value_lo[j];
			if( j==0 )
			{
				lis_quad_sub((LIS_QUAD *)t.hi,(LIS_QUAD *)t.hi,(LIS_QUAD *)one.hi);
			}
			lis_vector_axpyex_mmm(t,v[j],v[0]);
		}
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter       = iter2+1;
	solver->iter2      = iter;
	solver->resid     = nrm2;
	lis_free(h);
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}
#endif

#undef NWORK
#define NWORK 4
#undef __FUNC__
#define __FUNC__ "lis_fgmres_check_params"
LIS_INT lis_fgmres_check_params(LIS_SOLVER solver)
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
#define __FUNC__ "lis_fgmres_malloc_work"
LIS_INT lis_fgmres_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,restart,worklen,err;

	LIS_DEBUG_FUNC_IN;

	restart = solver->options[LIS_OPTIONS_RESTART];
	worklen = NWORK+(2*restart+1);
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_gmres_malloc_work::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
		for(i=1;i<worklen;i++)
		{
			err = lis_vector_duplicate(solver->A,&work[i]);
			if( err ) break;
		}
	}
	else
	{
		for(i=1;i<worklen;i++)
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
	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
		lis_vector_create(solver->A->comm,&work[0]);
	}
	else
	{
		lis_vector_createex(LIS_PRECISION_QUAD,solver->A->comm,&work[0]);
	}
	lis_vector_set_size(work[0],restart+1,0);
	solver->worklen = worklen;
	solver->work    = work;


	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_fgmres"
LIS_INT lis_fgmres(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR b,x;
	LIS_VECTOR s,*z,*v;
	LIS_SCALAR *h;
	LIS_SCALAR aa,bb,rr,a2,b2;
	LIS_SCALAR t;

	LIS_REAL bnrm2,nrm2,tol;
	LIS_INT iter,maxiter,output;
	double time,ptime;

	LIS_REAL rnorm;
	LIS_INT i,j,k,m;
	LIS_INT ii,i1,iiv,i1v,iih,jj;
	LIS_INT h_dim;
	LIS_INT cs,sn;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	b       = solver->b;
	x       = solver->x;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	m       = solver->options[LIS_OPTIONS_RESTART];
	h_dim   = m+1;
	ptime   = 0.0;

	s       = solver->work[0];
	z       = &solver->work[2];
	v       = &solver->work[m+2];

	h       = (LIS_SCALAR *)lis_malloc( sizeof(LIS_SCALAR)*(h_dim+1)*(h_dim+2),"lis_gmres::h" );
	cs      = (m+1)*h_dim;
	sn      = (m+2)*h_dim;

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,v[0],&bnrm2) )
	{
		lis_free(h);
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;
	rnorm   = 1.0/bnrm2;


	iter=0;
	while( iter<maxiter )
	{
		/* first column of V */
		/* v = r / ||r||_2 */
		lis_vector_scale(bnrm2,v[0]);

		/* s = ||r||_2 e_1 */
		lis_vector_set_all(0,s);
		s->value[0] = rnorm;

		i = 0;
		do
		{
			iter++;
			i++;
			ii  = i-1;
			i1  = i;
			iiv = i-1;
			i1v = i;
			iih = (i-1)*h_dim;

			/* z = M^-1 * v */
			time = lis_wtime();
			lis_psolve(solver,v[iiv],z[iiv]);
			ptime += lis_wtime()-time;

			/* w = A * z */
			lis_matvec(A,z[iiv], v[i1v]);

			for(k=0;k<i;k++)
			{
				/* h[k,i]   = <w,v[k]>          */
				/* w        = w - h[k,i] * v[k] */
				lis_vector_dot(v[i1v],v[k],&t);
				h[k+iih] = t;
				lis_vector_axpy(-t,v[k],v[i1v]);
			}
			/* h[i+1,i] = ||w||          */
			/* v[i+1]   = w / h[i+1,i]   */
			lis_vector_nrm2(v[i1v],(LIS_REAL *)&t);
			h[i1+iih] = t;
			lis_vector_scale(1.0/t,v[i1v]);

			for(k=1;k<=ii;k++)
			{
				jj  = k-1;
				t   =  h[jj+iih];
				aa  =  h[jj+cs]*t;
				aa +=  h[jj+sn]*h[k+iih];
				bb  = -h[jj+sn]*t;
				bb +=  h[jj+cs]*h[k+iih];
				h[jj+iih] = aa;
				h[k+iih] = bb;
			}
			aa = h[ii+iih];
			bb = h[i1+iih];
			a2 = aa*aa;
			b2 = bb*bb;
			rr = sqrt(a2+b2);
			if( rr==0.0 ) rr=1.0e-17;
			h[ii+cs] = aa/rr;
			h[ii+sn] = bb/rr;
			s->value[i1] = -h[ii+sn]*s->value[ii];
			s->value[ii] =  h[ii+cs]*s->value[ii];

			aa  =  h[ii+cs]*h[ii+iih];
			aa +=  h[ii+sn]*h[i1+iih];
			h[ii+iih] = aa;

			/* convergence check */
			nrm2 = fabs(s->value[i1]);

			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			if( tol >= nrm2 ) break;
		} while( i<m && iter <maxiter );

		/* Solve H * Y = S for upper Hessenberg matrix H */
		s->value[ii] = s->value[ii]/h[ii+iih];
		for(k=1;k<=ii;k++)
		{
			jj = ii-k;
			t  = s->value[jj];
			for(j=jj+1;j<=ii;j++)
			{
				t -= h[jj+j*h_dim]*s->value[j];
			}
			s->value[jj] = t/h[jj+jj*h_dim];
		}
		/* x = x + z * y */
		for(j=0;j<=ii;j++)
		{
			lis_vector_axpy(s->value[j],z[j],x);
		}

		if( tol >= nrm2 )
		{
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			lis_free(h);
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		lis_matvec(A,x,v[0]);
		lis_vector_xpay(b,-1.0,v[0]);
		lis_vector_nrm2(v[0],&rnorm);
		bnrm2 = 1.0/rnorm;
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter+1;
	solver->resid     = nrm2;
	lis_free(h);
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_fgmres_quad"
LIS_INT lis_fgmres_quad(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR b,x;
	LIS_VECTOR r,s,*z,*v;
	LIS_QUAD *h;
	LIS_QUAD_PTR aa,bb,rr,a2,b2,t,one,tmp;

	LIS_REAL bnrm2,nrm2,tol;
	LIS_INT iter,maxiter,n,output;
	double time,ptime;

	LIS_REAL rnorm;
	LIS_INT i,j,k,m;
	LIS_INT ii,i1,iiv,i1v,iih,jj;
	LIS_INT h_dim;
	LIS_INT cs,sn;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	b       = solver->b;
	x       = solver->x;
	n       = A->n;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	m       = solver->options[LIS_OPTIONS_RESTART];
	h_dim   = m+1;
	ptime   = 0.0;

	s       = solver->work[0];
	r       = solver->work[1];
	z       = &solver->work[2];
	v       = &solver->work[m+2];

	h       = (LIS_QUAD *)lis_malloc( sizeof(LIS_QUAD)*(h_dim+1)*(h_dim+2),"lis_fgmres_quad::h" );
	cs      = (m+1)*h_dim;
	sn      = (m+2)*h_dim;

	LIS_QUAD_SCALAR_MALLOC(aa,0,1);
	LIS_QUAD_SCALAR_MALLOC(bb,1,1);
	LIS_QUAD_SCALAR_MALLOC(rr,2,1);
	LIS_QUAD_SCALAR_MALLOC(a2,3,1);
	LIS_QUAD_SCALAR_MALLOC(b2,4,1);
	LIS_QUAD_SCALAR_MALLOC(t,5,1);
	LIS_QUAD_SCALAR_MALLOC(tmp,6,1);
	LIS_QUAD_SCALAR_MALLOC(one,7,1);

	one.hi[0]   = 1.0;
	one.lo[0]   = 0.0;

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,v[0],&bnrm2) )
	{
		lis_free(h);
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;
	rnorm   = 1.0/bnrm2;


	iter=0;
	while( iter<maxiter )
	{
		/* first column of V */
		/* v = r / ||r||_2 */
		lis_vector_scaleex_nm(bnrm2,v[0]);

		/* s = ||r||_2 e_1 */
		lis_vector_set_allex_nm(0.0,s);
		s->value[0]    = rnorm;
		s->value_lo[0] = 0.0;

		i = 0;
		do
		{
			iter++;
			i++;
			ii  = i-1;
			i1  = i;
			iiv = i-1;
			i1v = i;
			iih = (i-1)*h_dim;


			/* z = M^-1 * v */
			time = lis_wtime();
			lis_psolve(solver,v[iiv],z[iiv]);
			ptime += lis_wtime()-time;

			/* w = A * z */
			lis_matvec(A,z[iiv], v[i1v]);

			for(k=0;k<i;k++)
			{
				/* h[k,i]   = <w,v[k]>          */
				/* w        = w - h[k,i] * v[k] */
				lis_vector_dotex_mmm(v[i1v],v[k],&t);
				h[k+iih].hi = t.hi[0];
				h[k+iih].lo = t.lo[0];
				lis_quad_minus((LIS_QUAD *)t.hi);
				lis_vector_axpyex_mmm(t,v[k],v[i1v]);
			}
			/* h[i+1,i] = ||w||          */
			/* v[i+1]   = w / h[i+1,i]   */
			lis_vector_nrm2ex_mm(v[i1v],&t);
			h[i1+iih].hi = t.hi[0];
			h[i1+iih].lo = t.lo[0];
			lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)one.hi,(LIS_QUAD *)t.hi);
			lis_vector_scaleex_mm(tmp,v[i1v]);

			for(k=1;k<=ii;k++)
			{
				jj  = k-1;
				t.hi[0]   =  h[jj+iih].hi;
				t.lo[0]   =  h[jj+iih].lo;
				lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)&h[jj+cs],(LIS_QUAD *)t.hi);
				lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj+sn],(LIS_QUAD *)&h[k+iih]);
				lis_quad_add((LIS_QUAD *)aa.hi,(LIS_QUAD *)aa.hi,(LIS_QUAD *)tmp.hi);
				lis_quad_mul((LIS_QUAD *)bb.hi,(LIS_QUAD *)&h[jj+sn],(LIS_QUAD *)t.hi);
				lis_quad_minus((LIS_QUAD *)bb.hi);
				lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj+cs],(LIS_QUAD *)&h[k+iih]);
				lis_quad_add((LIS_QUAD *)bb.hi,(LIS_QUAD *)bb.hi,(LIS_QUAD *)tmp.hi);
				h[jj+iih].hi = aa.hi[0];
				h[jj+iih].lo = aa.lo[0];
				h[k+iih].hi = bb.hi[0];
				h[k+iih].lo = bb.lo[0];
			}
			aa.hi[0] = h[ii+iih].hi;
			aa.lo[0] = h[ii+iih].lo;
			bb.hi[0] = h[i1+iih].hi;
			bb.lo[0] = h[i1+iih].lo;
			lis_quad_sqr((LIS_QUAD *)a2.hi,(LIS_QUAD *)aa.hi);
			lis_quad_sqr((LIS_QUAD *)b2.hi,(LIS_QUAD *)bb.hi);
			lis_quad_add((LIS_QUAD *)rr.hi,(LIS_QUAD *)a2.hi,(LIS_QUAD *)b2.hi);
			lis_quad_sqrt((LIS_QUAD *)rr.hi,(LIS_QUAD *)rr.hi);
			if( rr.hi[0]==0.0 )
			{
				rr.hi[0]=1.0e-17;
				rr.lo[0]=0.0;
			}
			lis_quad_div((LIS_QUAD *)&h[ii+cs],(LIS_QUAD *)aa.hi,(LIS_QUAD *)rr.hi);
			lis_quad_div((LIS_QUAD *)&h[ii+sn],(LIS_QUAD *)bb.hi,(LIS_QUAD *)rr.hi);
			tmp.hi[0] = s->value[ii];
			tmp.lo[0] = s->value_lo[ii];
			lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)&h[ii+sn],(LIS_QUAD *)tmp.hi);
			lis_quad_mul((LIS_QUAD *)bb.hi,(LIS_QUAD *)&h[ii+cs],(LIS_QUAD *)tmp.hi);
			lis_quad_minus((LIS_QUAD *)aa.hi);
			s->value[i1] = aa.hi[0];
			s->value_lo[i1] = aa.lo[0];
			s->value[ii] = bb.hi[0];
			s->value_lo[ii] = bb.lo[0];

			lis_quad_mul((LIS_QUAD *)aa.hi,(LIS_QUAD *)&h[ii+cs],(LIS_QUAD *)&h[ii+iih]);
			lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[ii+sn],(LIS_QUAD *)&h[i1+iih]);
			lis_quad_add((LIS_QUAD *)aa.hi,(LIS_QUAD *)aa.hi,(LIS_QUAD *)tmp.hi);
			h[ii+iih].hi = aa.hi[0];
			h[ii+iih].lo = aa.lo[0];

			/* convergence check */
			nrm2 = fabs(s->value[i1]);

			if( output )
			{
				if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
				if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
			}

			if( tol >= nrm2 ) break;
		} while( i<m && iter <maxiter );

		/* Solve H * Y = S for upper Hessenberg matrix H */
		tmp.hi[0] = s->value[ii];
		tmp.lo[0] = s->value_lo[ii];
		lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[ii+iih]);
		s->value[ii] = tmp.hi[0];
		s->value_lo[ii] = tmp.lo[0];
		for(k=1;k<=ii;k++)
		{
			jj = ii-k;
			t.hi[0]  = s->value[jj];
			t.lo[0]  = s->value_lo[jj];
			for(j=jj+1;j<=ii;j++)
			{
				tmp.hi[0] = s->value[j];
				tmp.lo[0] = s->value_lo[j];
				lis_quad_mul((LIS_QUAD *)tmp.hi,(LIS_QUAD *)tmp.hi,(LIS_QUAD *)&h[jj+j*h_dim]);
				lis_quad_sub((LIS_QUAD *)t.hi,(LIS_QUAD *)t.hi,(LIS_QUAD *)tmp.hi);
			}
			lis_quad_div((LIS_QUAD *)tmp.hi,(LIS_QUAD *)t.hi,(LIS_QUAD *)&h[jj+jj*h_dim]);
			s->value[jj] = tmp.hi[0];
			s->value_lo[jj] = tmp.lo[0];
		}
		/* x = x + y * z */
		for(j=0;j<=ii;j++)
		{
			aa.hi[0] = s->value[j];
			aa.lo[0] = s->value_lo[j];
			lis_vector_axpyex_mmm(aa,z[j],x);
		}

		if( tol >= nrm2 )
		{
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			lis_free(h);
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		lis_matvec(A,x,v[0]);
		lis_vector_xpay(b,-1.0,v[0]);
		memset(v[0]->value_lo,0,n*sizeof(LIS_SCALAR));
		lis_vector_nrm2(v[0],&rnorm);
		bnrm2 = 1.0/rnorm;
	}

	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter+1;
	solver->resid     = nrm2;
	lis_free(h);
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}
#endif
