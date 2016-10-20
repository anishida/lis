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
 compute eigenvalues of a real symmetric tridiagonal matrix 
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
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SUBSPACE(=%d) is less than 0\n",ss);
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
  LIS_MATRIX A,B;
  LIS_INT ss,ic;
  LIS_SCALAR gshift;
  LIS_INT emaxiter,iter0,qriter;
  LIS_REAL tol,qrerr;
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

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  gshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN];      
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];

  t = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::t");
  tq = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::tq");
  tr = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::tr");
  
  A = esolver->A;
  B = esolver->B;
  r = esolver->work[0];
  v = &esolver->work[1];
  lis_vector_set_all(0.0,v[0]);
  lis_vector_set_all(1.0,r);
  lis_vector_nrm2(r, &nrm2);

  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p none",solver);
  lis_solver_set_optionC(solver);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);
  lis_esolver_get_esolvername(niesolver, esolvername);
  if( output & A->my_rank==0 )
    {
      printf("inner eigensolver     : %s\n", esolvername);
      printf("linear solver         : %s\n", solvername);
      printf("preconditioner        : %s\n", preconname);
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
	lis_vector_scale(1/nrm2, v[j]);
	lis_matvec(A, v[j], r);
      }
      else {
	lis_vector_scale(1/t[(j-2)*ss+j-1], v[j]);
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

  /* compute eigenvalues of a real symmetric tridiagonal matrix 
     T(j) = ST'(j)S^* */
  lis_array_qr(ss,t,tq,tr,&qriter,&qrerr);

  for (i=0;i<ss;i++)
    {
      esolver->evalue[i] = t[i*ss+i];
    }

  if( output & A->my_rank==0 ) 
    {
#ifdef _LONG__LONG
      printf("size of subspace      : %lld\n\n", ss);
#else
      printf("size of subspace      : %d\n\n", ss);
#endif
      printf("approximate eigenvalues in subspace:\n\n");
      for (i=0;i<ss;i++)
	{
#ifdef _LONG__LONG
	  printf("Lanczos: mode number              = %lld\n", i);
#else
	  printf("Lanczos: mode number              = %d\n", i);
#endif
#ifdef _COMPLEX	  
#ifdef _LONG__DOUBLE
	  printf("Lanczos: eigenvalue               = (%Le, %Le)\n", creall(esolver->evalue[i] - gshift), cimagl(esolver->evalue[i] - gshift));
#else
	  printf("Lanczos: eigenvalue               = (%e, %e)\n", creal(esolver->evalue[i] - gshift), cimag(esolver->evalue[i] - gshift));
#endif
#else
#ifdef _LONG__DOUBLE
	  printf("Lanczos: eigenvalue               = %Le\n", esolver->evalue[i] - gshift);
#else
	  printf("Lanczos: eigenvalue               = %e\n", esolver->evalue[i] - gshift);
#endif
#endif	  
	}
      printf("\n");
      printf("compute refined eigenpairs:\n\n");
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
      esolver2->lshift = -(esolver->evalue[i]);
      lis_gesolve(A, B, esolver->evector[i], &evalue, esolver2);
      lis_esolver_work_destroy(esolver2); 
      esolver->evalue[i] = evalue - esolver2->lshift;
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
	  esolver->ptime = esolver2->ptime;
	  esolver->itime = esolver2->itime;
	  esolver->p_c_time = esolver2->p_c_time;
	  esolver->p_i_time = esolver2->p_i_time;
	}

      if (output & A->my_rank==0) 
	{

#ifdef _LONG__LONG
	  printf("Lanczos: mode number          = %lld\n", i);
#else
	  printf("Lanczos: mode number          = %d\n", i);
#endif
#ifdef _COMPLEX
#ifdef _LONG__DOUBLE
	  printf("Lanczos: eigenvalue           = (%Le, %Le)\n", creall(esolver->evalue[i] - gshift), cimagl(esolver->evalue[i] - gshift));
#else
	  printf("Lanczos: eigenvalue           = (%e, %e)\n", creal(esolver->evalue[i] - gshift), cimag(esolver->evalue[i] - gshift));
#endif
#else	  
#ifdef _LONG__DOUBLE
	  printf("Lanczos: eigenvalue           = %Le\n", esolver->evalue[i] - gshift);
#else
	  printf("Lanczos: eigenvalue           = %e\n", esolver->evalue[i] - gshift);
#endif
#endif	  
#ifdef _LONG__LONG
	  printf("Lanczos: number of iterations = %lld\n",esolver2->iter[0]);
#else
	  printf("Lanczos: number of iterations = %d\n",esolver2->iter[0]);
#endif
#ifdef _LONG__DOUBLE
	  printf("Lanczos: relative residual    = %Le\n\n",esolver2->resid[0]);
#else
	  printf("Lanczos: relative residual    = %e\n\n",esolver2->resid[0]);
#endif	  
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

