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

/*
 * This subroutine is made based on the contribution by Hirofumi Notsu.
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

#define NWORK 7
/************************************************
 * lis_minres_check_params
 * lis_minres_malloc_work
 * lis_minres
 ************************************************/
#undef __FUNC__
#define __FUNC__ "lis_minres_check_params"
LIS_INT lis_minres_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_minres_malloc_work"
LIS_INT lis_minres_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR *work;
	LIS_INT i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_minres_malloc_work::work" );
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
#define __FUNC__ "lis_minres"
LIS_INT lis_minres(LIS_SOLVER solver)
{
  LIS_Comm comm;  
  LIS_MATRIX A;
  LIS_VECTOR b,x;
  LIS_VECTOR v1,v2,v3,v4,w0,w1,w2;
  LIS_REAL nrm2,tol;
  LIS_SCALAR alpha;
  LIS_REAL beta2,beta3;
  LIS_SCALAR gamma1,gamma2,gamma3;
  LIS_SCALAR delta,eta;
  LIS_SCALAR sigma1,sigma2,sigma3;
  LIS_SCALAR rho1,rho2,rho3;
  LIS_REAL r0_euc,r_euc; 
  LIS_INT iter,maxiter,output;
  double time,ptime;

  LIS_DEBUG_FUNC_IN;

  comm = LIS_COMM_WORLD;
  
  A       = solver->A;
  b       = solver->b;
  x       = solver->x;
  tol     = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
  maxiter = solver->options[LIS_OPTIONS_MAXITER];
  output  = solver->options[LIS_OPTIONS_OUTPUT];
  ptime   = 0.0;

  v1       = solver->work[0];
  v2       = solver->work[1];
  v3       = solver->work[2];
  v4       = solver->work[3];
  w0       = solver->work[4];
  w1       = solver->work[5];
  w2       = solver->work[6];

  /* Lanczos algorithm */
  lis_matvec(A,x,v2); 
  lis_vector_xpay(b,-1.0,v2);

  time = lis_wtime();
  lis_psolve(solver,v2,v3);
  ptime += lis_wtime()-time;
  lis_vector_copy(v3,v2);

  /* Compute elements of Hermitian tridiagonal matrix */
  lis_vector_nrm2(v2,&r_euc); 
  eta = beta2 = r0_euc = r_euc; 
  gamma2 = gamma1 = 1.0; 
  sigma2 = sigma1 = 0.0;

  lis_vector_set_all(0.0,v1); 
  lis_vector_set_all(0.0,w0); 
  lis_vector_set_all(0.0,w1);

  nrm2 = r_euc / r0_euc; 

  for(iter=1;iter<=maxiter;iter++)
    {

      /* Lanczos algorithm */
      lis_vector_scale(1.0 / beta2,v2); 

      lis_matvec(A,v2,v3); 
      time = lis_wtime();

      lis_psolve(solver,v3,v4);
      ptime += lis_wtime()-time;

      lis_vector_dot(v2,v4,&alpha);
      lis_vector_axpy(-alpha,v2,v4);
      lis_vector_axpy(-beta2,v1,v4);
      lis_vector_nrm2(v4,&beta3);

      /* Compute elements of Hermitian tridiagonal matrix */
      delta = gamma2 * alpha - gamma1 * sigma2 * beta2;
      rho1 = sqrt(delta * delta + beta3 * beta3); 
      rho2 = sigma2 * alpha + gamma1 * gamma2 * beta2; 
      rho3 = sigma1 * beta2;
      gamma3 = delta / rho1; 
      sigma3 = beta3 / rho1;

      lis_vector_axpyz(-rho3,w0,v2,w2); 
      lis_vector_axpy(-rho2,w1,w2); 
      lis_vector_scale(1.0 / rho1,w2);

      lis_vector_axpy(gamma3 * eta,w2,x);

      /* convergence check */
      r_euc *= fabs(sigma3);
      nrm2 = r_euc / r0_euc;
      
      if( output )
	{
	  if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
	  if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
	}
      
      if( nrm2 <= tol )
	{ 
	  solver->retcode    = LIS_SUCCESS;
	  solver->iter       = iter;
	  solver->resid      = nrm2;
	  solver->ptime      = ptime;
	  LIS_DEBUG_FUNC_OUT;
	  return LIS_SUCCESS;
	}

      eta *= -sigma3;

      lis_vector_copy(v2,v1); 
      lis_vector_copy(v4,v2);
      lis_vector_copy(w1,w0); 
      lis_vector_copy(w2,w1);

      beta2 = beta3;
      gamma1 = gamma2; 
      gamma2 = gamma3; 
      sigma1 = sigma2; 
      sigma2 = sigma3;

    }

  lis_vector_destroy(v1);
  lis_vector_destroy(v2); 
  lis_vector_destroy(v3);
  lis_vector_destroy(v4);
  lis_vector_destroy(w0); 
  lis_vector_destroy(w1); 
  lis_vector_destroy(w2);

  solver->retcode   = LIS_MAXITER;
  solver->iter      = iter;
  solver->resid     = nrm2;
  LIS_DEBUG_FUNC_OUT;
  return LIS_MAXITER;
}
