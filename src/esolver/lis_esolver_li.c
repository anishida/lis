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
 * Lanczos Iteration                   *
 ***************************************
 v(0)    = (0,...,0)^T
 ***************************************
 r = v(0)
 beta(0) = ||r||_2 
 for j=1,2,...
   v(j)      = r / beta(j-1)
   r         = A * v(j)
   r         = r - beta(j-1) * v(j-1)
   alpha(j)  = <v(j), r>
   r         = r - alpha(j) * v(j)
   reorthogonalization 
   beta(j)   = ||r||_2
 end for
 compute eigenvalues of a symmetric tridiagonal matrix 
 T(j) = ST'(j)S^*, where
     (alpha(1) beta(1)                      )
     (beta(1)  alpha(2)                     )
 T = (               ...                    )
     (                  alpha(j-1) beta(j-1))                       
     (                  beta(j-1)  alpha(j) )
 compute refined eigenpairs
 ***************************************/

#define NWORK 2
#undef __FUNC__
#define __FUNC__ "lis_eli_check_params"
LIS_INT lis_eli_check_params(LIS_ESOLVER esolver)
{
        LIS_INT ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
	if( ss<0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SUBSPACE(=%D) is less than 0\n",ss);
		return LIS_ERR_ILL_ARG;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_eli_malloc_work"
LIS_INT lis_eli_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err,ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

	worklen = NWORK + ss;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_eli_malloc_work::work" );
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
#define __FUNC__ "lis_eli"
LIS_INT lis_eli(LIS_ESOLVER esolver)
{
  LIS_Comm comm; 
  LIS_MATRIX A;
  LIS_INT err;
  LIS_INT ss,ic;
  LIS_INT emaxiter,iter0,qriter;
  LIS_REAL tol,qrerr;
  double time,time0;    
  LIS_INT i,j,k;
  LIS_INT output, niesolver;
  LIS_REAL nrm2,resid0;
  LIS_SCALAR dot;
  LIS_VECTOR *v,r;
  LIS_SCALAR *t,*tq,*tr,evalue,evalue0;
  LIS_SOLVER solver;
  LIS_ESOLVER esolver2;
  char esolvername[128],solvername[128],preconname[128];
  LIS_INT nsol,precon_type;
  LIS_INT rval;

  LIS_DEBUG_FUNC_IN;

  comm = LIS_COMM_WORLD;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];
  rval = esolver->options[LIS_EOPTIONS_RVAL];  

  t = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::t");
  tq = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::tq");
  tr = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::tr");
  
  A = esolver->A;
  r = esolver->work[0];
  v = &esolver->work[1];
  lis_vector_set_all(0.0,v[0]);
  lis_vector_set_all(1.0,r);
  lis_vector_nrm2(r, &nrm2);

  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p none",solver);
  err = lis_solver_set_optionC(solver);
  CHKERR(err);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);
  lis_esolver_get_esolvername(niesolver, esolvername);
  if( output )
    {
      lis_printf(comm,"inner eigensolver     : %s\n", esolvername);
      lis_printf(comm,"linear solver         : %s\n", solvername);
      lis_printf(comm,"preconditioner        : %s\n", preconname);
    }

  for (i=0;i<ss*ss;i++) t[i] = 0.0;

  j=0;
  while (j<ss-1)
    {
      j = j+1;
      lis_vector_copy(r, v[j]);

      /* v(j) = r / beta(j-1) */
      /* r = A * v(j) */
      /* r = r - beta(j-1) * v(j-1) */
      if (j==1) {
	lis_vector_scale(1.0/nrm2, v[j]);
	lis_matvec(A, v[j], r);
      }
      else {
	lis_vector_scale(1.0/t[(j-2)*ss+j-1], v[j]);
	lis_matvec(A, v[j], r);
	lis_vector_axpy(-t[(j-2)*ss+j-1], v[j-1], r); 
      }

      /* alpha(j)  = <v(j), r> */
      lis_vector_dot(v[j], r, &t[(j-1)*ss+j-1]); 

      /* r = r - alpha(j) * v(j) */
      lis_vector_axpy(-t[(j-1)*ss+j-1], v[j], r); 

      /* reorthogonalization */
      for (k=1;k<j;k++)
	{ 
	  lis_vector_dot(v[j], v[k], &dot); 
	  lis_vector_axpy(-dot, v[k], v[j]);
	}

      /* beta(j) = ||r||_2 */
      lis_vector_nrm2(r, (LIS_REAL *)&t[(j-1)*ss+j]);

      /* convergence check */
      if (fabs(t[(j-1)*ss+j])<tol) break;  
      t[j*ss+j-1] = t[(j-1)*ss+j];
    }

  /* compute eigenvalues of a symmetric tridiagonal matrix 
     T(j) = ST'(j)S^* */
  time0 = lis_wtime();
  lis_array_qr(ss,t,tq,tr,&qriter,&qrerr);
  time = lis_wtime() - time0;

  for (i=0;i<ss;i++)
    {
      esolver->evalue[i] = t[i*ss+i];
    }

  if( output ) 
    {
      lis_printf(comm,"size of subspace      : %D\n\n", ss);
      lis_printf(comm,"Ritz values:\n\n");
      for (i=0;i<ss;i++)
	{
	  lis_printf(comm,"Lanczos: mode number          = %D\n", i);
#ifdef _COMPLEX	  
	  lis_printf(comm,"Lanczos: Ritz value           = (%e, %e)\n", (double)creal(esolver->evalue[i]), (double)cimag(esolver->evalue[i]));
#else
	  lis_printf(comm,"Lanczos: Ritz value           = %e\n", (double)esolver->evalue[i]);
#endif	  
	}
      lis_printf(comm,"Lanczos: elapsed time         = %e sec.\n\n", time);
    }
  if( rval )
    {
      lis_free(t); 
      lis_free(tq);
      lis_free(tr);
      lis_solver_destroy(solver);
      return LIS_SUCCESS;
    }
  
  if( output ) 
    {
      lis_printf(comm,"computing refined eigenpairs using inner eigensolver:\n\n");      
    }

  lis_esolver_create(&esolver2);
  esolver2->options[LIS_EOPTIONS_ESOLVER] = niesolver;
  esolver2->options[LIS_EOPTIONS_SUBSPACE] = 1;
  esolver2->options[LIS_EOPTIONS_MAXITER] = emaxiter;
  esolver2->options[LIS_EOPTIONS_OUTPUT] = esolver->options[LIS_EOPTIONS_OUTPUT];
  esolver2->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN] = tol; 

  /* compute refined eigenpairs */
  for (i=0;i<ss;i++)
    {
      lis_vector_duplicate(A, &esolver->evector[i]); 
      esolver2->ishift = esolver->evalue[i];
      lis_esolve(A, esolver->evector[i], &evalue, esolver2);
      lis_esolver_work_destroy(esolver2); 
      esolver->evalue[i] = evalue;
      esolver->iter[i] = esolver2->iter[0];      
      esolver->resid[i] = esolver2->resid[0];

      if (i==0) 
	{
	  evalue0 = esolver->evalue[0];
	  iter0 = esolver2->iter[0];
	  resid0 = esolver2->resid[0];
	  if( output & LIS_EPRINT_MEM ) 
	    {
	      for (ic=0;ic<iter0+1;ic++)
		{
		  esolver->rhistory[ic] = esolver2->rhistory[ic]; 
		}
	    }
	  esolver->time = esolver2->time;	  
	  esolver->ptime += esolver2->ptime;
	  esolver->itime += esolver2->itime;
	  esolver->p_c_time += esolver2->p_c_time;
	  esolver->p_i_time += esolver2->p_i_time;
	}

      if (output) 
	{
	  lis_printf(comm,"Lanczos: mode number          = %D\n", i);
#ifdef _COMPLEX
	  lis_printf(comm,"Lanczos: eigenvalue           = (%e, %e)\n", (double)creal(esolver->evalue[i]), (double)cimag(esolver->evalue[i]));
#else	  
	  lis_printf(comm,"Lanczos: eigenvalue           = %e\n", (double)esolver->evalue[i]);
#endif
	  lis_printf(comm,"Lanczos: elapsed time         = %e sec.\n", esolver2->time);	  	  
	  lis_printf(comm,"Lanczos: number of iterations = %D\n",esolver2->iter[0]);
	  lis_printf(comm,"Lanczos: relative residual    = %e\n\n",(double)esolver2->resid[0]);
	}
    }
  esolver->evalue[0] = evalue0; 
  esolver->iter[0] = iter0;
  esolver->resid[0] = resid0;

  lis_vector_copy(esolver->evector[0], esolver->x);

  lis_esolver_destroy(esolver2); 

  lis_free(t); 
  lis_free(tq);
  lis_free(tr);

  lis_solver_destroy(solver);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

/***************************************
 * Generalized Lanczos Iteration       *
 ***************************************
 q       = (1,...,1)^T
 ***************************************
 r = B * q
 beta(0) = |<q,r>|^1/2
 for j=1,2,...
   w(j)      = r / beta(j-1)
   v(j)      = q / beta(j-1)
   r         = A * v(j)
   r         = r - beta(j-1) * w(j-1)
   alpha(j)  = <v(j), r>
   r         = r - alpha(j) * w(j)
   reorthogonalization 
   solve B * q = r for q
   beta(j)   = |<q,r>|^1/2
 end for
 compute eigenvalues of a symmetric tridiagonal matrix 
 T(j) = ST'(j)S^*, where
     (alpha(1) beta(1)                      )
     (beta(1)  alpha(2)                     )
 T = (               ...                    )
     (                  alpha(j-1) beta(j-1))                       
     (                  beta(j-1)  alpha(j) )
 compute refined eigenpairs
 ***************************************/

#undef NWORK
#define NWORK 4
#undef __FUNC__
#define __FUNC__ "lis_egli_check_params"
LIS_INT lis_egli_check_params(LIS_ESOLVER esolver)
{
        LIS_INT ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
	if( ss<0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SUBSPACE(=%D) is less than 0\n",ss);
		return LIS_ERR_ILL_ARG;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_egli_malloc_work"
LIS_INT lis_egli_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err,ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

	worklen = NWORK + ss;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_egli_malloc_work::work" );
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
#define __FUNC__ "lis_egli"
LIS_INT lis_egli(LIS_ESOLVER esolver)
{
  LIS_Comm comm;  
  LIS_INT err;
  LIS_MATRIX A,B;
  LIS_INT ss,ic;
  LIS_INT emaxiter,iter0,qriter;
  LIS_REAL tol,qrerr;
  double time,time0;    
  LIS_INT i,j,k;
  LIS_INT output, nigesolver;
  LIS_REAL resid0;
  LIS_SCALAR beta,dot;
  LIS_VECTOR q,r,*w,*v;
  LIS_SCALAR *t,*tq,*tr,evalue,evalue0;
  LIS_SOLVER solver;
  LIS_PRECON precon;
  LIS_ESOLVER esolver2;
  char esolvername[128],solvername[128],preconname[128];
  LIS_INT nsol,precon_type;
  LIS_INT rval;  

  LIS_DEBUG_FUNC_IN;

  comm = LIS_COMM_WORLD;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  nigesolver = esolver->options[LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER];
  rval = esolver->options[LIS_EOPTIONS_RVAL];    

  t = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_egli::t");
  tq = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_egli::tq");
  tr = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_egli::tr");
  
  A = esolver->A;
  B = esolver->B;
  q = esolver->work[0];  
  r = esolver->work[1];
  w = &esolver->work[2];  
  v = &esolver->work[3];
  lis_vector_set_all(0.0,w[0]);  
  lis_vector_set_all(0.0,v[0]);
  lis_vector_set_all(1.0,q);

  /* create linear solver */
  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p none",solver);
  err = lis_solver_set_optionC(solver);
  CHKERR(err);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);
  lis_esolver_get_esolvername(nigesolver, esolvername);
  if( output )
    {
      lis_printf(comm,"inner eigensolver     : %s\n", esolvername);
      lis_printf(comm,"linear solver         : %s\n", solvername);
      lis_printf(comm,"preconditioner        : %s\n", preconname);
    }

  /* create preconditioner */
  solver->A = B;
  err = lis_precon_create(solver, &precon);
  if( err )
    {
      lis_solver_work_destroy(solver);
      solver->retcode = err;
      return err;
    }

  for (i=0;i<ss*ss;i++) t[i] = 0.0;

  j=0;
  while (j<ss-1)
    {
      j = j+1;

      /* w(j) = r / beta(j-1) */      
      /* v(j) = r / beta(j-1) */
      /* r = A * v(j) */
      /* r = r - beta(j-1) * w(j-1) */
      
      if (j==1) {
 	lis_matvec(B, q, r);
	lis_vector_dot(q, r, &beta);
	beta = sqrt(fabs(beta));
	lis_vector_copy(r, w[j]);
	lis_vector_copy(q, v[j]);	
	lis_vector_scale(1.0/beta, w[j]);
	lis_vector_scale(1.0/beta, v[j]);	
	lis_matvec(A, v[j], r);
	lis_vector_axpy(-beta, w[j-1], r); 	
      }
      else {
	lis_vector_copy(r, w[j]);
	lis_vector_copy(q, v[j]);	
	lis_vector_scale(1.0/t[(j-2)*ss+j-1], w[j]);	
	lis_vector_scale(1.0/t[(j-2)*ss+j-1], v[j]);
	lis_matvec(A, v[j], r);
	lis_vector_axpy(-t[(j-2)*ss+j-1], w[j-1], r); 
      }

      /* alpha(j)  = <v(j), r> */
      lis_vector_dot(v[j], r, &t[(j-1)*ss+j-1]); 

      /* r = r - alpha(j) * w(j) */
      lis_vector_axpy(-t[(j-1)*ss+j-1], w[j], r); 

      /* reorthogonalization */
      for (k=1;k<j;k++)
	{ 
	  lis_vector_dot(v[j], v[k], &dot); 
	  lis_vector_axpy(-dot, v[k], v[j]);
	}

      /* solve B * q = r */
      err = lis_solve_kernel(B, r, q, solver, precon);
      if( err )
	{
	  lis_solver_work_destroy(solver);	  
	  solver->retcode = err;
	  return err;
	}

      /* beta(j) = |<q,r>|^1/2 */
      lis_vector_dot(q, r, &beta);
      beta = sqrt(fabs(beta));
      t[(j-1)*ss+j] = (LIS_REAL)beta;
      
      /* convergence check */
      if (fabs(t[(j-1)*ss+j])<tol) break;  
      t[j*ss+j-1] = t[(j-1)*ss+j];
    }

  /* compute eigenvalues of a symmetric tridiagonal matrix 
     T(j) = ST'(j)S^* */
  time0 = lis_wtime();
  lis_array_qr(ss,t,tq,tr,&qriter,&qrerr);
  time = lis_wtime() - time0;

  for (i=0;i<ss;i++)
    {
      esolver->evalue[i] = t[i*ss+i];
    }

  if( output ) 
    {
      lis_printf(comm,"size of subspace      : %D\n\n", ss);
      lis_printf(comm,"Ritz values :\n\n");
      for (i=0;i<ss;i++)
	{
	  lis_printf(comm,"Generalized Lanczos: mode number          = %D\n", i);
#ifdef _COMPLEX	  
	  lis_printf(comm,"Generalized Lanczos: Ritz value           = (%e, %e)\n", (double)creal(esolver->evalue[i]), (double)cimag(esolver->evalue[i]));
#else
	  lis_printf(comm,"Generalized Lanczos: Ritz value           = %e\n", (double)esolver->evalue[i]);
#endif	  
	}
      lis_printf(comm,"Lanczos: elapsed time         = %e sec.\n\n", time);  
    }
  if( rval ) return LIS_SUCCESS;
      
  if( output ) 
    {
      lis_printf(comm,"computing refined eigenpairs using inner eigensolver:\n\n");      
    }

  lis_esolver_create(&esolver2);
  esolver2->options[LIS_EOPTIONS_ESOLVER] = nigesolver;
  esolver2->options[LIS_EOPTIONS_SUBSPACE] = 1;
  esolver2->options[LIS_EOPTIONS_MAXITER] = emaxiter;
  esolver2->options[LIS_EOPTIONS_OUTPUT] = esolver->options[LIS_EOPTIONS_OUTPUT];
  esolver2->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN] = tol; 

  /* compute refined eigenpairs */
  for (i=0;i<ss;i++)
    {
      lis_vector_duplicate(A, &esolver->evector[i]); 
      esolver2->ishift = esolver->evalue[i];
      lis_gesolve(A, B, esolver->evector[i], &evalue, esolver2);
      lis_esolver_work_destroy(esolver2); 
      esolver->evalue[i] = evalue;
      esolver->iter[i] = esolver2->iter[0];      
      esolver->resid[i] = esolver2->resid[0];

      if (i==0) 
	{
	  evalue0 = esolver->evalue[0];
	  iter0 = esolver2->iter[0];
	  resid0 = esolver2->resid[0];
	  if( output & LIS_EPRINT_MEM ) 
	    {
	      for (ic=0;ic<iter0+1;ic++)
		{
		  esolver->rhistory[ic] = esolver2->rhistory[ic]; 
		}
	    }
	  esolver->time = esolver2->time;	  
	  esolver->ptime += esolver2->ptime;
	  esolver->itime += esolver2->itime;
	  esolver->p_c_time += esolver2->p_c_time;
	  esolver->p_i_time += esolver2->p_i_time;
	}

      if (output) 
	{
	  lis_printf(comm,"Generalized Lanczos: mode number          = %D\n", i);
#ifdef _COMPLEX
	  lis_printf(comm,"Generalized Lanczos: eigenvalue           = (%e, %e)\n", (double)creal(esolver->evalue[i]), (double)cimag(esolver->evalue[i]));
#else	  
	  lis_printf(comm,"Generalized Lanczos: eigenvalue           = %e\n", (double)esolver->evalue[i]);
#endif
	  lis_printf(comm,"Generalized Lanczos: elapsed time         = %e sec.\n", esolver2->time);	  	  
	  lis_printf(comm,"Generalized Lanczos: number of iterations = %D\n",esolver2->iter[0]);
	  lis_printf(comm,"Generalized Lanczos: relative residual    = %e\n\n",(double)esolver2->resid[0]);
	}
    }
  esolver->evalue[0] = evalue0; 
  esolver->iter[0] = iter0;
  esolver->resid[0] = resid0;

  lis_vector_copy(esolver->evector[0], esolver->x);

  lis_esolver_destroy(esolver2); 

  lis_free(t); 
  lis_free(tq);
  lis_free(tr);
  
  lis_precon_destroy(precon);
  lis_solver_destroy(solver);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}



