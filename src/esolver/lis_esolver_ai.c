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
 * Arnoldi Iteration                   *
 ***************************************
 v(0)    = (1,...,1)^T
 v(0)    = v(0)/||v(0)||_2
 ***************************************
 for j=0,1,2,...,m-1
   w = A * v(j)
   for i=0,1,2,...,j-1
     h(i,j) = <w, v(i)>
     w = w - h(i,j) * v(i)
   end for
   h(j+1,j) = ||w||_2
   if h(j+1,j) = 0, stop
   v(j+1) = w / h(j+1,j)
 end for
 compute eigenvalues of an upper Hessenberg matrix 
 H(m) = SH'(m)S^*, where
       (h(0,0)   h(0,1)                            )
       (h(1,0)   h(1,1)                            )
       (  0      h(2,1)                            )
   H = (           0   ...                         )
       (                      h(m-3,m-2) h(m-3,m-1))
       (                      h(m-2,m-2) h(m-2,m-1))                       
       (                   0  h(m-1,m-2)   h(m-1)  )
 compute refined eigenpairs
 ***************************************/

#define NWORK 2
#undef __FUNC__
#define __FUNC__ "lis_eai_check_params"
LIS_INT lis_eai_check_params(LIS_ESOLVER esolver)
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
#define __FUNC__ "lis_eai_malloc_work"
LIS_INT lis_eai_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err,ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

	worklen = NWORK + ss;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_eai_malloc_work::work" );
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
#define __FUNC__ "lis_eai"
LIS_INT lis_eai(LIS_ESOLVER esolver)
{
  LIS_Comm comm;  
  LIS_MATRIX A;
  LIS_INT err;
  LIS_INT ss,ic;
  LIS_INT emaxiter,iter0,hqriter;
  LIS_REAL tol,hqrerr,D;
  double time,time0;    
  LIS_INT i,j;
  LIS_INT output, niesolver;
  LIS_REAL nrm2,resid0; 
  LIS_VECTOR *v,w;
  LIS_SCALAR *h,*hq,*hr,evalue,evalue0;
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

  h = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eai::h");
  hq = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eai::hq");
  hr = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eai::hr");
  
  A = esolver->A;
  w = esolver->work[0];
  v = &esolver->work[1];
  lis_vector_set_all(1.0,v[0]);
  lis_vector_nrm2(v[0],&nrm2);
  lis_vector_scale(1.0/nrm2,v[0]);
  
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

  for (i=0;i<ss*ss;i++) h[i] = 0.0;

  j=-1;
  while (j<ss-1)
    {
      j = j+1;

      /* w = A * v(j) */
      lis_matvec(A, v[j], w);

      /* reorthogonalization */
      for (i=0;i<=j;i++)
	{
	  /* h(i,j) = <v(i), w> */
	  lis_vector_dot(v[i], w, &h[i+j*ss]);
	  /* w = w - h(i,j) * v(i) */
	  lis_vector_axpy(-h[i+j*ss], v[i], w); 
	}

      /* h(j+1,j) = ||w||_2 */
      lis_vector_nrm2(w, (LIS_REAL *)&h[j+1+j*ss]);

      /* convergence check */
      if (fabs(h[j+1+j*ss])<tol) break;

      /* v(j+1) = w / h(j,j-1) */
      lis_vector_scale(1/h[j+1+j*ss],w);
      lis_vector_copy(w,v[j+1]);
      
    }

  /* compute eigenvalues of an upper
     Hessenberg matrix H(j) = SH'(j)S^* */
  time0 = lis_wtime();
  lis_array_qr(ss,h,hq,hr,&hqriter,&hqrerr);
  time = lis_wtime() - time0;
  
  if( output ) 
    {
      lis_printf(comm,"size of subspace      : %D\n\n", ss);
      lis_printf(comm,"Ritz values:\n\n");
    }
  
  i=0;
  while (i<ss) 
    {
      i = i + 1;

      /* eigenvalues on the diagonal of H */
      if (fabs(h[i+(i-1)*ss])<tol)
	{
	  if( output ) 
	    {
	      lis_printf(comm,"Arnoldi: mode number          = %D\n",i-1);
#ifdef _COMPLEX
	      lis_printf(comm,"Arnoldi: Rits value           = (%e, %e)\n",(double)creal(h[i-1+(i-1)*ss]),(double)cimag(h[i-1+(i-1)*ss]));
#else	      
	      lis_printf(comm,"Arnoldi: Ritz value           = %e\n",(double)(h[i-1+(i-1)*ss]));
#endif
	    }
	  esolver->evalue[i-1] = h[i-1+(i-1)*ss];
	}

      /* eigenvalues of 2x2 matrices on the diagonal of H */
      else
	{
	  D = (h[i-1+(i-1)*ss]+h[i+i*ss])*(h[i-1+(i-1)*ss]+h[i+i*ss])
	    - 4*(h[i-1+(i-1)*ss]*h[i+i*ss]-h[i-1+i*ss]*h[i+(i-1)*ss]);
	  if (D<0)
	    {
	      if( output ) 
		{
		  lis_printf(comm,"Arnoldi: mode number          = %D\n",i-1);
		  lis_printf(comm,"Arnoldi: Ritz value           = (%e, %e)\n", (double)((h[i-1+(i-1)*ss]+h[i+i*ss])/2), (double)sqrt(-D)/2);
		  lis_printf(comm,"Arnoldi: mode number          = %D\n",i);
		  lis_printf(comm,"Arnoldi: Ritz value           = (%e, %e)\n", (double)((h[i-1+(i-1)*ss]+h[i+i*ss])/2), (double)-sqrt(-D)/2);
		}
#ifdef _COMPLEX		  
	      esolver->evalue[i-1] = (h[i-1+(i-1)*ss]+h[i+i*ss])/2 + sqrt(-D)/2 * _Complex_I;
	      esolver->evalue[i]   = (h[i-1+(i-1)*ss]+h[i+i*ss])/2 - sqrt(-D)/2 * _Complex_I;
#else
	      esolver->evalue[i-1] = (h[i-1+(i-1)*ss]+h[i+i*ss])/2;
	      esolver->evalue[i]   = (h[i-1+(i-1)*ss]+h[i+i*ss])/2;     
#endif
		  
	      i=i+1;
	    }
	  else
	    {
	      if( output ) 
		{
		  lis_printf(comm,"Arnoldi: mode number          = %D\n",i-1);
#ifdef _COMPLEX
		  lis_printf(comm,"Arnoldi: Ritz value           = (%e, %e)\n",(double)creal(h[i-1+(i-1)*ss]),(double)cimag(h[i-1+(i-1)*ss]));
#else		  
		  lis_printf(comm,"Arnoldi: Ritz value           = %e\n",(double)(h[i-1+(i-1)*ss]));
#endif
		}
	      esolver->evalue[i-1] = h[i-1+(i-1)*ss];
	    }
	}
    }
  lis_printf(comm,"Arnoldi: elapsed time         = %e sec.\n\n", time);
  
  if( rval )
    {
      lis_free(h); 
      lis_free(hq);
      lis_free(hr);
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

  /* compute refined (real) eigenpairs */

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
	  lis_printf(comm,"Arnoldi: mode number          = %D\n", i);
#ifdef _COMPLEX	  
	  lis_printf(comm,"Arnoldi: eigenvalue           = (%e, %e)\n", (double)creal(esolver->evalue[i]),(double)cimag(esolver->evalue[i]));
#else
	  lis_printf(comm,"Arnoldi: eigenvalue           = %e\n", (double)esolver->evalue[i]);
#endif
	  lis_printf(comm,"Arnoldi: elapsed time         = %e sec.\n", esolver2->time);	  
	  lis_printf(comm,"Arnoldi: number of iterations = %D\n",esolver2->iter[0]);
	  lis_printf(comm,"Arnoldi: relative residual    = %e\n\n",(double)esolver2->resid[0]);
	}
    }

  lis_vector_copy(esolver->evector[0], esolver->x);

  lis_esolver_destroy(esolver2); 

  lis_free(h); 
  lis_free(hq);
  lis_free(hr);

  lis_solver_destroy(solver);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
} 

/***************************************
 * Generalized Arnoldi Iteration       *
 ***************************************
 v(0)    = (1,...,1)^T
 v(0)    = v(0)/||v(0)||_2
 ***************************************
 for j=0,1,2,...,m-1
   w = B^-1 * A * v(j)
   for i=0,1,2,...,j-1
     h(i,j) = <w, v(i)>
     w = w - h(i,j) * v(i)
   end for
   h(j+1,j) = ||w||_2
   if h(j+1,j) = 0, stop
   v(j+1) = w / h(j+1,j)
 end for
 compute eigenvalues of an upper Hessenberg matrix 
 H(m) = SH'(m)S^*, where
       (h(0,0)   h(0,1)                            )
       (h(1,0)   h(1,1)                            )
       (  0      h(2,1)                            )
   H = (           0   ...                         )
       (                      h(m-3,m-2) h(m-3,m-1))
       (                      h(m-2,m-2) h(m-2,m-1))                       
       (                   0  h(m-1,m-2)   h(m-1)  )
 compute refined eigenpairs
 ***************************************/

#undef NWORK
#define NWORK 3
#undef __FUNC__
#define __FUNC__ "lis_egai_check_params"
LIS_INT lis_egai_check_params(LIS_ESOLVER esolver)
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
#define __FUNC__ "lis_egai_malloc_work"
LIS_INT lis_egai_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err,ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

	worklen = NWORK + ss;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_egai_malloc_work::work" );
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
#define __FUNC__ "lis_egai"
LIS_INT lis_egai(LIS_ESOLVER esolver)
{
  LIS_Comm comm;  
  LIS_INT err;
  LIS_MATRIX A,B;
  LIS_INT ss,ic;
  LIS_INT emaxiter,iter0,hqriter,iter2;
  LIS_REAL tol,hqrerr,D;
  double time,time0;  
  LIS_INT i,j;
  LIS_INT output, nigesolver;
  LIS_REAL nrm2,resid0; 
  LIS_VECTOR *v,w,y;
  LIS_SCALAR *h,*hq,*hr,evalue,evalue0;
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
  
  h = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_egai::h");
  hq = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_egai::hq");
  hr = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_egai::hr");
  
  A = esolver->A;
  B = esolver->B;  
  w = esolver->work[0];
  y = esolver->work[1];
  v = &esolver->work[2];
  lis_vector_set_all(1.0,v[0]);
  lis_vector_nrm2(v[0],&nrm2);
  lis_vector_scale(1.0/nrm2,v[0]);
  
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

  for (i=0;i<ss*ss;i++) h[i] = 0.0;

  j=-1;
  while (j<ss-1)
    {
      j = j+1;

      /* w = A * v(j) */
      lis_matvec(A, v[j], w);

      /* w = B^-1 * w */
      err = lis_solve_kernel(B, w, y, solver, precon);
      if( err )
	{
	  lis_solver_work_destroy(solver);	  
	  solver->retcode = err;
	  return err;
	}
      lis_solver_get_iter(solver, &iter2);
      lis_vector_copy(y, w);
      
      /* reorthogonalization */
      for (i=0;i<=j;i++)
	{
	  /* h(i,j) = <v(i), w> */
	  lis_vector_dot(v[i], w, &h[i+j*ss]);
	  /* w = w - h(i,j) * v(i) */
	  lis_vector_axpy(-h[i+j*ss], v[i], w); 
	}

      /* h(j+1,j) = ||w||_2 */
      lis_vector_nrm2(w, (LIS_REAL *)&h[j+1+j*ss]);

      /* convergence check */
      if (fabs(h[j+1+j*ss])<tol) break;

      /* v(j+1) = w / h(j,j-1) */
      lis_vector_scale(1/h[j+1+j*ss],w);
      lis_vector_copy(w,v[j+1]);
      
    }

  /* compute eigenvalues of an upper
     Hessenberg matrix H(j) = SH'(j)S^* */
  time0 = lis_wtime();
  lis_array_qr(ss,h,hq,hr,&hqriter,&hqrerr);
  time = lis_wtime() - time0;
  
  if( output ) 
    {
      lis_printf(comm,"size of subspace      : %D\n\n", ss);
      lis_printf(comm,"Ritz values:\n\n");
    }

  i=0;
  while (i<ss) 
    {
      i = i + 1;

      /* eigenvalues on the diagonal of H */
      if (fabs(h[i+(i-1)*ss])<tol)
	{
	  if( output ) 
	    {
	      lis_printf(comm,"Generalized Arnoldi: mode number              = %D\n",i-1);
#ifdef _COMPLEX
	      lis_printf(comm,"Generalized Arnoldi: Ritz value               = (%e, %e)\n",(double)creal(h[i-1+(i-1)*ss]),(double)cimag(h[i-1+(i-1)*ss]));
#else	      
	      lis_printf(comm,"Generalized Arnoldi: Ritz value               = %e\n",(LIS_REAL)(h[i-1+(i-1)*ss]));
#endif
	    }
	  esolver->evalue[i-1] = h[i-1+(i-1)*ss];
	}

      /* eigenvalues of 2x2 matrices on the diagonal of H */
      else
	{
	  D = (h[i-1+(i-1)*ss]+h[i+i*ss])*(h[i-1+(i-1)*ss]+h[i+i*ss])
	    - 4*(h[i-1+(i-1)*ss]*h[i+i*ss]-h[i-1+i*ss]*h[i+(i-1)*ss]);
	  if (D<0)
	    {
	      if( output ) 
		{
		  lis_printf(comm,"Generalized Arnoldi: mode number              = %D\n",i-1);
		  lis_printf(comm,"Generalized Arnoldi: Ritz value               = (%e, %e)\n", (double)((h[i-1+(i-1)*ss]+h[i+i*ss])/2), (double)sqrt(-D)/2);
		  lis_printf(comm,"Generalized Arnoldi: mode number              = %D\n",i);
		  lis_printf(comm,"Generalized Arnoldi: Ritz value               = (%e, %e)\n", (double)((h[i-1+(i-1)*ss]+h[i+i*ss])/2), (double)-sqrt(-D)/2);
		}
#ifdef _COMPLEX		  
	      esolver->evalue[i-1] = (h[i-1+(i-1)*ss]+h[i+i*ss])/2 + sqrt(-D)/2 * _Complex_I;
	      esolver->evalue[i]   = (h[i-1+(i-1)*ss]+h[i+i*ss])/2 - sqrt(-D)/2 * _Complex_I;
#else
	      esolver->evalue[i-1] = (h[i-1+(i-1)*ss]+h[i+i*ss])/2;
	      esolver->evalue[i]   = (h[i-1+(i-1)*ss]+h[i+i*ss])/2;     
#endif
		  
	      i=i+1;
	    }
	  else
	    {
	      if( output ) 
		{
		  lis_printf(comm,"Generalized Arnoldi: mode number              = %D\n",i-1);
#ifdef _COMPLEX
		  lis_printf(comm,"Generalized Arnoldi: Ritz value               = (%e, %e)\n",(double)creal(h[i-1+(i-1)*ss]),(double)cimag(h[i-1+(i-1)*ss]));
#else		  
		  lis_printf(comm,"Generalized Arnoldi: Ritz value               = %e\n",(double)(h[i-1+(i-1)*ss]));
#endif
		}
	      esolver->evalue[i-1] = h[i-1+(i-1)*ss];
	    }
	}
    }
  lis_printf(comm,"Generalized Arnoldi: elapsed time         = %e sec.\n\n", time);	      

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

  /* compute refined (real) eigenpairs */

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

	  lis_printf(comm,"Generalized Arnoldi: mode number          = %D\n", i);
#ifdef _COMPLEX	  
	  lis_printf(comm,"Generalized Arnoldi: eigenvalue           = (%e, %e)\n", (double)creal(esolver->evalue[i]),(double)cimag(esolver->evalue[i]));
#else
	  lis_printf(comm,"Generalized Arnoldi: eigenvalue           = %e\n", (double)esolver->evalue[i]);
#endif
	  lis_printf(comm,"Generalized Arnoldi: elapsed time         = %e sec.\n", esolver2->time);	  
	  lis_printf(comm,"Generalized Arnoldi: number of iterations = %D\n",esolver2->iter[0]);
	  lis_printf(comm,"Generalized Arnoldi: relative residual    = %e\n\n",(double)esolver2->resid[0]);
	}
    }

  lis_vector_copy(esolver->evector[0], esolver->x);

  lis_esolver_destroy(esolver2); 

  lis_free(h); 
  lis_free(hq);
  lis_free(hr);

  lis_precon_destroy(precon);  
  lis_solver_destroy(solver);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
} 

