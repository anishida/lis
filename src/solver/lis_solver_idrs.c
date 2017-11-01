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
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

#define PRE_RIGHT

#define NWORK 4
/************************************************
 * lis_idrs_check_params
 * lis_idrs_malloc_work
 * lis_idrs
 ************************************************/
#undef __FUNC__
#define __FUNC__ "lis_idrs_check_params"
LIS_INT lis_idrs_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_idrs_malloc_work"
LIS_INT lis_idrs_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,s,worklen,err;

	LIS_DEBUG_FUNC_IN;

	/*
	err = lis_matrix_convert(solver->A,&solver->Ah,LIS_MATRIX_CCS);
	if( err ) return err;
	*/

	s       = solver->options[LIS_OPTIONS_IDRS_RESTART];
	worklen = NWORK + 3*s;
	work    = (LIS_VECTOR *)lis_malloc(
worklen*sizeof(LIS_VECTOR),"lis_idrs_malloc_work::work" );
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
			err =
lis_vector_duplicateex(LIS_PRECISION_QUAD,solver->A,&work[i]);
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

#define NWORK 4
#undef __FUNC__
#define __FUNC__ "lis_idr1_check_params"
LIS_INT lis_idr1_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_idr1_malloc_work"
LIS_INT lis_idr1_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,s,worklen,err;

	LIS_DEBUG_FUNC_IN;

	/*
	err = lis_matrix_convert(solver->A,&solver->Ah,LIS_MATRIX_CCS);
	if( err ) return err;
	*/

	s       = solver->options[LIS_OPTIONS_IDRS_RESTART];
	worklen = NWORK + 3*s;
	work    = (LIS_VECTOR *)lis_malloc(
worklen*sizeof(LIS_VECTOR),"lis_idrs_malloc_work::work" );
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
			err =
lis_vector_duplicateex(LIS_PRECISION_QUAD,solver->A,&work[i]);
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
#define __FUNC__ "lis_idrs_omega"
void lis_idrs_omega(LIS_VECTOR t, LIS_VECTOR s, LIS_SCALAR angle, LIS_SCALAR
*om)
{
	LIS_SCALAR nt;
	LIS_SCALAR ts;

	lis_vector_dot(t,t,&nt);
	lis_vector_dot(t,s,&ts);
	*om = ts / nt;
}

#undef __FUNC__
#define __FUNC__ "lis_idrs_orth"
void lis_idrs_orth(LIS_INT s, LIS_VECTOR *P)
{
	LIS_INT i,j;
	LIS_REAL r;
	LIS_SCALAR d;

	for(j=0;j<s;j++)
	{
		lis_vector_nrm2(P[j],&r);
		r = 1.0/r;
		lis_vector_scale(r,P[j]);
		for(i=j+1;i<s;i++)
		{
			lis_vector_dot(P[j],P[i],&d);
			lis_vector_axpy(-d,P[j],P[i]);
		}
	}
}

#undef __FUNC__
#define __FUNC__ "lis_idr1"
LIS_INT lis_idr1(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,t,v,av,*dX,*dR,*P;
	LIS_SCALAR om, h;
	LIS_SCALAR M,m,c;
	LIS_REAL bnrm2,nrm2,tol;
	LIS_INT i,s;
	LIS_INT iter,maxiter,n,output,conv;
	double time,ptime,tim;
        unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	x       = solver->x;
	n       = A->n;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	s       = 1;
	ptime   = 0.0;

	r       = solver->work[0];
	t       = solver->work[1];
	v       = solver->work[2];
	av      = solver->work[3];
	P       = &solver->work[4];
	dX      = &solver->work[4+s];
	dR      = &solver->work[4+2*s];

	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	init_by_array(init, length);
		for(i=0;i<n;i++)
		{
			P[0]->value[i] = genrand_real1();
		}
		/*
	lis_vector_copy(r,P[0]);
		*/
	lis_idrs_orth(s,P);

		#ifdef PRE_RIGHT
			time = lis_wtime();
			lis_psolve(solver, r, dX[0]);
			ptime += lis_wtime()-time;
			lis_matvec(A,dX[0],dR[0]);
		#else
		#ifdef PRE_BOTH
			time = lis_wtime();
			lis_psolve_right(solver, r, t);
			ptime += lis_wtime()-time;
			lis_matvec(A,t,av);
			lis_vector_print(av);
			time = lis_wtime();
			lis_psolve_left(solver, av, v);
			ptime += lis_wtime()-time;
		#endif
		#endif

			/*
		lis_idrs_omega(dR[k],r,angle,&om);
			*/
		lis_vector_dot(dR[0],dR[0],&h);
		lis_vector_dot(dR[0],r,&om);
		om = om / h;
		lis_vector_scale(om,dX[0]);
		lis_vector_scale(-om,dR[0]);

		lis_vector_axpy(1.0,dX[0],x);
		lis_vector_axpy(1.0,dR[0],r);


		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);

		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[1] = nrm2;
			if( output & LIS_PRINT_OUT ) lis_printf(comm,"iteration: %5d  relative residual = %e\n", 1, (double)nrm2);
		}

		if( tol >= nrm2 )
		{

			solver->retcode    = LIS_SUCCESS;
			solver->iter       = 1;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		lis_vector_dot(P[0],dR[0],&M);

	iter = s;
	lis_vector_dot(P[0],r,&m);

	while( iter<=maxiter )
	{
		tim = lis_wtime();

		/* solve Mc=m */
		c = m/M;

		for(i=0;i<n;i++)
		{
			v->value[i] = r->value[i] + -c*dR[0]->value[i];
		}
		/*
		lis_vector_copy(r,v);
		lis_vector_axpy(-c,dR[0],v);
		*/

			#ifdef PRE_RIGHT
				time = lis_wtime();
				lis_psolve(solver, v, av);
				ptime += lis_wtime()-time;
				lis_matvec(A,av,t);
			#else
			#ifdef PRE_BOTH
				time = lis_wtime();
				lis_psolve_right(solver, v, t);
				ptime += lis_wtime()-time;
				lis_matvec(A,t,av);
				time = lis_wtime();
				lis_psolve_left(solver, av, t);
				ptime += lis_wtime()-time;
			#endif
			#endif

				/*
			lis_idrs_omega(t,v,angle,&om);
			lis_vector_dot(t,t,&h);
			lis_vector_dot(t,v,&om);
				*/
			h  = t->value[0]*t->value[0];
			om = t->value[0]*v->value[0];
			for(i=1;i<n;i++)
			{
				h  += t->value[i]*t->value[i];
				om += t->value[i]*v->value[i];
			}
			om = om / h;
			/*
			lis_printf(comm,"i=%D om = %lf\n",iter,om);
			*/
			#if 0
				lis_vector_scale(-om,t);
				for(j=0;j<s;j++)
				{
					lis_vector_axpy(-c[j],dR[j],t);
				}
				lis_vector_copy(t,dR[oldest]);
				lis_vector_scale(om,av);
				for(j=0;j<s;j++)
				{
					lis_vector_axpy(-c[j],dX[j],av);
				}
				lis_vector_copy(av,dX[oldest]);
			#else
				for(i=0;i<n;i++)
				{
					h = om*av->value[i];
					h -= dX[0]->value[i] * c;
					dX[0]->value[i] = h;
					h = -om*t->value[i];
					h -= dR[0]->value[i] * c;
					dR[0]->value[i] = h;
				}
			#endif

		lis_vector_axpy(1.0,dR[0],r);
		lis_vector_axpy(1.0,dX[0],x);

		iter++;

		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);

		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
#ifdef _LONG__LONG
			if( output & LIS_PRINT_OUT ) lis_printf(comm,"iteration: %5lld  relative residual = %e\n", iter, (double)nrm2);
#else
			if( output & LIS_PRINT_OUT ) lis_printf(comm,"iteration: %5d  relative residual = %e\n", iter, (double)nrm2);
#endif
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

		lis_vector_dot(P[0],dR[0],&h);
		m += h;
		M = h;




		/* solve Mc=m */
		c = m/M;

		for(i=0;i<n;i++)
		{
			v->value[i] = r->value[i] + -c*dR[0]->value[i];
		}
		/*
		lis_vector_copy(r,v);
		lis_vector_axpy(-c,dR[0],v);
		*/

			#ifdef PRE_RIGHT
				time = lis_wtime();
				lis_psolve(solver, v, av);
				ptime += lis_wtime()-time;
			#endif

			#if 0
				lis_vector_scale(om,av);
				for(j=0;j<s;j++)
				{
					lis_vector_axpy(-c[j],dX[j],av);
				}
				lis_vector_copy(av,dX[oldest]);
			#else
				for(i=0;i<n;i++)
				{
					h = om*av->value[i];
					h -= dX[0]->value[i] * c;
					dX[0]->value[i] = h;
				}
			#endif

			lis_matvec(A,dX[0],dR[0]);
			lis_vector_scale(-1.0,dR[0]);

		lis_vector_axpy(1.0,dR[0],r);
		lis_vector_axpy(1.0,dX[0],x);

		iter++;

		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);

		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
#ifdef _LONG__LONG
			if( output & LIS_PRINT_OUT ) lis_printf(comm,"iteration: %5lld  relative residual = %e\n", iter, (double)nrm2);
#else
			if( output & LIS_PRINT_OUT ) lis_printf(comm,"iteration: %5d  relative residual = %e\n", iter, (double)nrm2);
#endif
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

		lis_vector_dot(P[0],dR[0],&h);
		m += h;
		M = h;

		tim = lis_wtime() - tim;
		/*
		lis_printf(comm,"update m,M: %e\n",(double)tim);
		*/
	}
	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

#undef __FUNC__
#define __FUNC__ "lis_idrs"
LIS_INT lis_idrs(LIS_SOLVER solver)
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x;
	LIS_VECTOR r,t,v,av,*dX,*dR,*P;
	LIS_SCALAR om, h;
	LIS_SCALAR *M,*m,*c,*MM;
	LIS_REAL bnrm2,nrm2,tol;
	LIS_INT i,j,k,s,oldest;
	LIS_INT iter,maxiter,n,output,conv;
	double time,ptime,tim;
        unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A       = solver->A;
	x       = solver->x;
	n       = A->n;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	s       = solver->options[LIS_OPTIONS_IDRS_RESTART];
	ptime   = 0.0;

	r       = solver->work[0];
	t       = solver->work[1];
	v       = solver->work[2];
	av      = solver->work[3];
	dX      = &solver->work[4];
	P       = &solver->work[4+s];
	dR      = &solver->work[4+2*s];

	m = (LIS_SCALAR *)lis_malloc(s*sizeof(LIS_SCALAR), "lis_idrs::m");
	c = (LIS_SCALAR *)lis_malloc(s*sizeof(LIS_SCALAR), "lis_idrs::c");
	M = (LIS_SCALAR *)lis_malloc(s*s*sizeof(LIS_SCALAR), "lis_idrs::M");
	MM = (LIS_SCALAR *)lis_malloc(s*s*sizeof(LIS_SCALAR),
"lis_idrs::M");



	/* Initial Residual */
	if( lis_solver_get_initial_residual(solver,NULL,NULL,r,&bnrm2) )
	{
		lis_free2(4,m,c,M,MM);
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	tol     = solver->tol;

	init_by_array(init, length);
	for(k=0;k<s;k++)
	{
		for(i=0;i<n;i++)
		{
			P[k]->value[i] = genrand_real1();
		}
	}
	lis_idrs_orth(s,P);

	for( k=0; k<s; k++ )
	{
		#ifdef PRE_RIGHT
			time = lis_wtime();
			lis_psolve(solver, r, dX[k]);
			ptime += lis_wtime()-time;
			lis_matvec(A,dX[k],dR[k]);
		#endif

		lis_vector_dot(dR[k],dR[k],&h);
		lis_vector_dot(dR[k],r,&om);
		om = om / h;
		lis_vector_scale(om,dX[k]);
		lis_vector_scale(-om,dR[k]);

		lis_vector_axpy(1.0,dX[k],x);
		lis_vector_axpy(1.0,dR[k],r);


		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);

		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[k+1] = nrm2;
#ifdef _LONG__LONG
			if( output & LIS_PRINT_OUT ) lis_printf(comm,"iteration: %5lld  relative residual = %e\n", k+1, (double)nrm2);
#else
			if( output & LIS_PRINT_OUT ) lis_printf(comm,"iteration: %5d  relative residual = %e\n", k+1, (double)nrm2);
#endif
		}

		if( tol >= nrm2 )
		{
			lis_free2(4,m,c,M,MM);

			solver->retcode    = LIS_SUCCESS;
			solver->iter       = k+1;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		for(i=0;i<s;i++)
		{
			lis_vector_dot(P[i],dR[k],&M[k*s+i]);
		}
	}

	iter = s;
	oldest = 0;
	for(i=0;i<s;i++)
	{
		lis_vector_dot(P[i],r,&m[i]);
	}

	while( iter<=maxiter )
	{
		tim = lis_wtime();
		lis_array_solve(s,M,m,c,MM); /* solve Mc=m */

		lis_vector_copy(r,v);
		for(j=0;j<s;j++)
		{
			lis_vector_axpy(-c[j],dR[j],v);
		}

		if( (iter%(s+1))==s )
		{
			#ifdef PRE_RIGHT
				time = lis_wtime();
				lis_psolve(solver, v, av);
				ptime += lis_wtime()-time;
				lis_matvec(A,av,t);
			#endif

			lis_vector_dot(t,t,&h);
			lis_vector_dot(t,v,&om);
			om = om / h;
			#if 0
				lis_vector_scale(-om,t);
				for(j=0;j<s;j++)
				{
					lis_vector_axpy(-c[j],dR[j],t);
				}
				lis_vector_copy(t,dR[oldest]);
				lis_vector_scale(om,av);
				for(j=0;j<s;j++)
				{
					lis_vector_axpy(-c[j],dX[j],av);
				}
				lis_vector_copy(av,dX[oldest]);
			#else
				for(i=0;i<n;i++)
				{
					h = om*av->value[i];
					for(j=0;j<s;j++)
					{
						h -= dX[j]->value[i] * c[j];
					}
					dX[oldest]->value[i] = h;
				}
				for(i=0;i<n;i++)
				{
					h = -om*t->value[i];
					for(j=0;j<s;j++)
					{
						h -= dR[j]->value[i] * c[j];
					}
					dR[oldest]->value[i] = h;
				}
			#endif
		}
		else
		{
			#ifdef PRE_RIGHT
				time = lis_wtime();
				lis_psolve(solver, v, av);
				ptime += lis_wtime()-time;
			#endif

			#if 0
				lis_vector_scale(om,av);
				for(j=0;j<s;j++)
				{
					lis_vector_axpy(-c[j],dX[j],av);
				}
				lis_vector_copy(av,dX[oldest]);
			#else
				for(i=0;i<n;i++)
				{
					h = om*av->value[i];
					for(j=0;j<s;j++)
					{
						h -= dX[j]->value[i] * c[j];
					}
					dX[oldest]->value[i] = h;
				}
			#endif

			lis_matvec(A,dX[oldest],dR[oldest]);
			lis_vector_scale(-1.0,dR[oldest]);
		}

		lis_vector_axpy(1.0,dR[oldest],r);
		lis_vector_axpy(1.0,dX[oldest],x);

		iter++;

		/* convergence check */
		lis_solver_get_residual[conv](r,solver,&nrm2);

		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
#ifdef _LONG__LONG
			if( output & LIS_PRINT_OUT ) lis_printf(comm,"iteration: %5lld  relative residual = %e\n", iter, (double)nrm2);
#else
			if( output & LIS_PRINT_OUT ) lis_printf(comm,"iteration: %5d  relative residual = %e\n", iter, (double)nrm2);
#endif
		}

		if( tol >= nrm2 )
		{
			lis_free2(4,m,c,M,MM);

			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptime      = ptime;
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}

		for(i=0;i<s;i++)
		{
			lis_vector_dot(P[i],dR[oldest],&h);
			m[i] += h;
			M[oldest*s+i] = h;
		}

		oldest++;
		if( oldest==s ) oldest = 0;
		tim = lis_wtime() - tim;
		/*
		lis_printf(comm,"update m,M: %e\n",(double)tim);
		*/
	}
	lis_free2(4,m,c,M,MM);
	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

