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
#include <math.h>
#include <string.h>
#include <stdarg.h>
#ifdef USE_SSE2
	#include <emmintrin.h>
#endif
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/***************************************
 * Conjugate Gradient Method           *
 ***************************************
 x(0)=(1,...,1)^T 
 x(0)=x(0)/||x(0)||_2
 ***************************************
 for k=0,1,...
   mu=<x(k),x(k)>/<x(k),A*x(k)>
   r=x(k)-mu*A*x(k)
   w=M^-1*r
   use Rayleigh-Ritz method for I-mu*A on span {w,x(k),x(k-1)}
   x(k+1)=Ritz vector correspoinding to the smallest Ritz value
 ***************************************/

#define NWORK 6
#undef __FUNC__
#define __FUNC__ "lis_ecg_check_params"
LIS_INT lis_ecg_check_params(LIS_ESOLVER esolver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_ecg_malloc_work"
LIS_INT lis_ecg_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;

	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_ecg_malloc_work::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	if( esolver->eprecision==LIS_PRECISION_DEFAULT )
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicate(esolver->A,&work[i]);
			if( err ) break;
		}
	}
	else
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,esolver->A,&work[i]);
			if( err ) break;
		}
	}
	if( i<worklen )
	{
		for(j=0;j<i;j++) lis_vector_destroy(work[j]);
		lis_free(work);
		return err;
	}
	esolver->worklen = worklen;
	esolver->work    = work;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_ecg"
LIS_INT lis_ecg(LIS_ESOLVER esolver)
{
  LIS_MATRIX A;
  LIS_VECTOR x;
  LIS_SCALAR evalue;
  LIS_INT emaxiter;
  LIS_REAL tol;
  LIS_INT iter,iter3,output;
  LIS_REAL nrm2,resid,resid3;
  LIS_SCALAR lshift;
  LIS_VECTOR r,w,p,Aw,Ax,Ap;
  LIS_SCALAR *A3,*B3,*W3,*v3,*A3v3,*B3v3,*z3,*q3,*B3z3,ievalue3;
  LIS_SOLVER solver;
  LIS_PRECON precon;
  double time,itime,ptime,p_c_time,p_i_time;
  LIS_INT nsol,precon_type;
  char solvername[128],preconname[128];

  A = esolver->A;
  x = esolver->x;
  if (esolver->options[LIS_EOPTIONS_INITGUESS_ONES] ) 
    {
      lis_vector_set_all(1.0,x);
    }
  
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  lshift = esolver->lshift;

#ifdef _COMPLEX
#ifdef _LONG__DOUBLE
  if( A->my_rank==0 ) printf("local shift           : (%Le, %Le)\n", creall(lshift), cimagl(lshift));
#else
  if( A->my_rank==0 ) printf("local shift           : (%e, %e)\n", creal(lshift), cimag(lshift));
#endif
#else  
#ifdef _LONG__DOUBLE
  if( A->my_rank==0 ) printf("local shift           : %Le\n", lshift);
#else
  if( A->my_rank==0 ) printf("local shift           : %e\n", lshift);
#endif
#endif  
  if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);

  A3 = (LIS_SCALAR *)lis_malloc(3*3*sizeof(LIS_SCALAR), "lis_ecg::A3");
  B3 = (LIS_SCALAR *)lis_malloc(3*3*sizeof(LIS_SCALAR), "lis_ecg::B3");
  W3 = (LIS_SCALAR *)lis_malloc(3*3*sizeof(LIS_SCALAR), "lis_ecg::W3");
  v3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::v3");
  z3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::z3");
  q3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::q3");
  A3v3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::A3v3");
  B3v3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::B3v3");
  B3z3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::B3z3");

  r = esolver->work[0];
  w = esolver->work[1];
  p = esolver->work[2];
  Ax = esolver->work[3];
  Aw = esolver->work[4];
  Ap = esolver->work[5];

  ptime = 0;

  lis_vector_nrm2(x, &nrm2);
  lis_vector_scale(1/nrm2, x);
  lis_matvec(A,x,Ax);
  lis_solver_create(&solver);
  lis_solver_set_option("-i cg -p none",solver);
  lis_solver_set_optionC(solver);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);

  if( A->my_rank==0 ) printf("linear solver         : %s\n", solvername);
  if( A->my_rank==0 ) printf("preconditioner        : %s\n", preconname);

  /* lis_solve must be called before lis_precon_create */
  /* p=A^-1*x */

  lis_solve(A,x,p,solver);
  lis_vector_copy(x,Ap);

  lis_precon_create(solver,&precon);
  solver->precon = precon;

  iter=0;

  while (iter<emaxiter)
    {
      iter=iter+1;

      /* mu=<x,x>/<x,A*x>, where mu=1/lambda */
      lis_vector_dot(x,Ax,&evalue);

      /* r=x-mu*A*x */
      lis_vector_axpyz(-1/evalue,Ax,x,r); 
      lis_vector_nrm2(r,&nrm2);
 
      /* convergence check */
      resid = nrm2;
      if( output )
	{
	  if( output & LIS_EPRINT_MEM ) esolver->rhistory[iter] = resid;
	  if( output & LIS_EPRINT_OUT && A->my_rank==0 ) lis_print_rhistory(iter,resid);
	}
      if (resid<tol) break;  

      /* w=M^-1*r */
      time = lis_wtime();
      lis_psolve(solver,r,w);
      ptime += lis_wtime() - time;
      lis_matvec(A,w,Aw);

      /* Rayleigh-Ritz method for I-mu*A on span {w,x(k),x(k-1)} */

      /* solve eigenproblem A_3*x=lambda*B_3*x, 
         where A_3 and B_3 are matrices of size 3*3 */
      lis_vector_dot(w,Aw,&A3[0]);
      lis_vector_dot(x,Aw,&A3[3]);
      lis_vector_dot(p,Aw,&A3[6]);
      /* (w,A*x)=(x,A*w) */
      A3[1]=A3[3];
      lis_vector_dot(x,Ax,&A3[4]);
      lis_vector_dot(p,Ax,&A3[7]);
      A3[2]=A3[6];
      A3[5]=A3[7];
      lis_vector_dot(p,Ap,&A3[8]);

      lis_vector_dot(w,w,&B3[0]);
      lis_vector_dot(x,w,&B3[3]);
      lis_vector_dot(p,w,&B3[6]);
      B3[1]=B3[3];
      lis_vector_dot(x,x,&B3[4]);
      lis_vector_dot(p,x,&B3[7]);
      B3[2]=B3[6];
      B3[5]=B3[7];
      lis_vector_dot(p,p,&B3[8]);
      
      /* compute eigenvector v_3 of size 3 using inverse iteration */
      lis_array_set_all(3,1.0,v3);
      iter3=0;
      while (iter3<emaxiter)
	{
	  iter3=iter3+1;
	  lis_array_nrm2(3,v3,&nrm2); 
	  lis_array_scale(3,1/nrm2,v3);
	  lis_array_matvec(3,B3,v3,B3v3,LIS_INS_VALUE);
	  lis_array_solve(3,A3,B3v3,z3,W3);
	  lis_array_dot(3,B3v3,z3,&ievalue3);
	  lis_array_axpyz(3,-ievalue3,B3v3,z3,q3);
	  lis_array_nrm2(3,q3,&resid3); 
	  if (resid3<tol) break;   
	  lis_array_copy(3,z3,v3);
	}

      /* update w, x, and p */
      lis_vector_scale(v3[0],w);  
      lis_vector_axpy(v3[2],p,w);
      lis_vector_xpay(w,v3[1],x);
      lis_vector_copy(w,p);

      /* update A*w, A*x, and A*p */      
      lis_vector_scale(v3[0],Aw);  
      lis_vector_axpy(v3[2],Ap,Aw);
      lis_vector_xpay(Aw,v3[1],Ax);
      lis_vector_copy(Aw,Ap);
      
      /* compute Ritz vector x corresponding to the smallest Ritz value */
      lis_vector_nrm2(x,&nrm2);
      lis_vector_scale(1/nrm2,x);
      lis_vector_scale(1/nrm2,Ax);
      
      lis_vector_nrm2(p,&nrm2);
      lis_vector_scale(1/nrm2,p);
      lis_vector_scale(1/nrm2,Ap);
      
    }

  lis_precon_destroy(precon);
  lis_solver_destroy(solver);

  esolver->iter[0]    = iter;
  esolver->resid[0]   = resid;
  esolver->evalue[0]  = evalue;

  lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
  esolver->ptime = ptime;
  esolver->itime = solver->itime;
  esolver->p_c_time = solver->p_c_time;
  esolver->p_i_time = solver->p_i_time;

  if (lshift != 0) lis_matrix_shift_diagonal(A,-lshift);

  lis_free(A3);
  lis_free(B3);
  lis_free(W3);
  lis_free(v3);
  lis_free(z3);
  lis_free(q3);
  lis_free(A3v3);
  lis_free(B3v3);
  lis_free(B3z3);

  if (resid<tol) 
    {
      esolver->retcode = LIS_SUCCESS;
      return LIS_SUCCESS;
    }
  else
    {
      esolver->retcode = LIS_MAXITER;
      return LIS_MAXITER;
    }
}

/***************************************
 * Conjugate Residual Method           *
 ***************************************
 x(0)=(1,...,1)^T 
 x(0)=x(0)/||x(0)||_2
 ***************************************
 lambda=<A*x(0),x(0)>/<x(0),x(0)> 
 r=lambda*x(0)-A*x(0)
 p=r
 for k=0,1,...
   alpha=(<r,A*p>-lambda*<r,p>) 
        /(<A*p,A*p>-2*lambda*<A*p,p>+lambda^2*<p,p>)
   x(k+1)=x(k)+alpha*p 
   lambda=<A*x(k+1),x(k+1)>/<x(k+1),x(k+1)>
   r=lambda*x(k+1)-A*x(k+1) 
   w=M^-1*r
   beta=-[<A*w,A*p>-lambda*(<A*w,p>+<A*p,w>)+lambda^2*<w,p>] 
        /[<A*p,A*p>-2*lambda*<A*p,p>+lambda^2*<p,p>] 
   p=w+beta*p
   A*p=A*w+beta*A*p
   convergence check
 ***************************************/

#undef NWORK
#define NWORK 6
#undef __FUNC__
#define __FUNC__ "lis_ecr_check_params"
LIS_INT lis_ecr_check_params(LIS_ESOLVER esolver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
 
#undef __FUNC__
#define __FUNC__ "lis_ecr_malloc_work"
LIS_INT lis_ecr_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;

	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_ecr_malloc_work::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	if( esolver->eprecision==LIS_PRECISION_DEFAULT )
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicate(esolver->A,&work[i]);
			if( err ) break;
		}
	}
	else
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,esolver->A,&work[i]);
			if( err ) break;
		}
	}
	if( i<worklen )
	{
		for(j=0;j<i;j++) lis_vector_destroy(work[j]);
		lis_free(work);
		return err;
	}
	esolver->worklen = worklen;
	esolver->work    = work;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_ecr"
LIS_INT lis_ecr(LIS_ESOLVER esolver)
{
  LIS_MATRIX A;
  LIS_VECTOR x;
  LIS_SCALAR evalue;
  LIS_INT emaxiter;
  LIS_REAL tol;
  LIS_INT iter,output;
  LIS_REAL nrm2,resid;
  LIS_SCALAR lshift;
  LIS_VECTOR r,p,w,Ax,Ap,Aw;
  LIS_SCALAR alpha,beta;
  LIS_SCALAR rAp,rp,ApAp,pAp,pp,AwAp,pAw,wAp,wp;
  LIS_SOLVER solver;
  LIS_PRECON precon;
  double time,itime,ptime,p_c_time,p_i_time;
  LIS_INT nsol,precon_type;
  char solvername[128],preconname[128];

  A = esolver->A;
  x = esolver->x;
  if (esolver->options[LIS_EOPTIONS_INITGUESS_ONES] ) 
    {
      lis_vector_set_all(1.0,x);
    }

  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  lshift = esolver->lshift;

#ifdef _COMPLEX
#ifdef _LONG__DOUBLE
  if( A->my_rank==0 ) printf("local shift           : (%Le, %Le)\n", creall(lshift), cimagl(lshift));
#else
  if( A->my_rank==0 ) printf("local shift           : (%e, %e)\n", creal(lshift), cimag(lshift));
#endif
#else  
#ifdef _LONG__DOUBLE
  if( A->my_rank==0 ) printf("local shift           : %Le\n", lshift);
#else
  if( A->my_rank==0 ) printf("local shift           : %e\n", lshift);
#endif
#endif  
  if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);

  r = esolver->work[0];
  p = esolver->work[1];
  w = esolver->work[2];  
  Ax = esolver->work[3];
  Ap = esolver->work[4];
  Aw = esolver->work[5];

  ptime = 0;

  lis_vector_nrm2(x, &nrm2);
  lis_vector_scale(1/nrm2, x);
  lis_matvec(A,x,Ax);
  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p none",solver);
  lis_solver_set_optionC(solver);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);

  if( A->my_rank==0 ) printf("linear solver         : %s\n", solvername);
  if( A->my_rank==0 ) printf("preconditioner        : %s\n", preconname);

  /* lis_solve must be called before lis_precon_create */
  /* p=A^-1*x */
  lis_solve(A,x,p,solver);
  lis_vector_copy(x,Ap);

  lis_precon_create(solver,&precon);
  solver->precon = precon;

  /* lambda=<A*x,x>/<x,x> */
  lis_vector_dot(x,Ax,&evalue); 
  
  /* r=lambda*x-A*x */
  lis_vector_axpyz(-evalue,x,Ax,r);
  lis_vector_scale(-1.0,r);

  iter=0;

  while (iter<emaxiter)
    {
      iter = iter + 1;

      /* alpha=(<r,A*p>-lambda*<r,p>) 
              /(<A*p,A*p>-2*lambda*<A*p,p>+lambda^2*<p,p>) */
      lis_vector_dot(r,Ap,&rAp);
      lis_vector_dot(r,p,&rp);
      lis_vector_dot(Ap,Ap,&ApAp);
      lis_vector_dot(p,Ap,&pAp); 
      lis_vector_dot(p,p,&pp);
      alpha = (rAp - evalue * rp)
	    / (ApAp - evalue * (2.0 * pAp - evalue * pp));

      /* x=x+alpha*p */
      lis_vector_axpy(alpha,p,x);

      /* lambda=<A*x,x>/<x,x> */
      lis_matvec(A,x,Ax);
      lis_vector_dot(x,Ax,&evalue);
      lis_vector_nrm2(x, &nrm2);
      evalue = evalue / (nrm2 * nrm2);

      /* r=lambda*x-A*x */
      lis_vector_axpyz(-evalue,x,Ax,r); 
      lis_vector_scale(-1.0,r);

      /* w=M^-1*r */
      time = lis_wtime();
      lis_psolve(solver,r,w);
      ptime += lis_wtime() - time;
      lis_matvec(A,w,Aw);

      /* beta=-[<A*w,A*p>-lambda*(<A*w,p>+<A*p.w>)+lambda^2*<w,p>] 
              /[<A*p,A*p>-2*lambda*<A*p,p>+lambda^2*<p,p>] */
      lis_vector_dot(Aw,Ap,&AwAp);
      lis_vector_dot(p,Aw,&pAw);
      lis_vector_dot(w,Ap,&wAp);      
      lis_vector_dot(w,p,&wp);
      beta = - (AwAp - evalue * ((pAw + wAp) - evalue * wp))
	     / (ApAp - evalue * (2.0 * pAp - evalue * pp));

      /* p=w+beta*p */
      lis_vector_xpay(w,beta,p);

      /* A*p=A*w+beta*A*p */
      lis_vector_xpay(Aw,beta,Ap);

      /* convergence check */
      lis_vector_nrm2(r,&nrm2);
      resid = fabs(nrm2 / (evalue));

      if( output )
	{
	  if( output & LIS_EPRINT_MEM ) esolver->rhistory[iter] = resid;
	  if( output & LIS_EPRINT_OUT && A->my_rank==0 ) lis_print_rhistory(iter,resid);
	}

      if (resid<tol) break;  
      
    }

  lis_precon_destroy(precon);
  lis_solver_destroy(solver);

  esolver->iter[0]    = iter;
  esolver->resid[0]   = resid;
  esolver->evalue[0]  = evalue;
  lis_vector_nrm2(x, &nrm2);
  lis_vector_scale(1/nrm2, x);

  lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
  esolver->ptime = ptime;
  esolver->itime = solver->itime;
  esolver->p_c_time = solver->p_c_time;
  esolver->p_i_time = solver->p_i_time;

  if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);

  if (resid<tol) 
    {
      esolver->retcode = LIS_SUCCESS;
      return LIS_SUCCESS;
    }
  else
    {
      esolver->retcode = LIS_MAXITER;
      return LIS_MAXITER;
    }
}
