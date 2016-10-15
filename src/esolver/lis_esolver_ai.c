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
 v(0)    = (0,...,0)^T
 ***************************************
 for j=1,2,...,m
   w = A * v(j)
   for i=1,2,...,j
     h(i,j) = <w, v(i)>
     w = w - h(i,j) * v(i)
   end for
   h(j+1,j) = ||w||_2
   if h(j+1,j) = 0, stop
   u(j+1) = w / h(j+1,j)
 end for
 compute eigenvalues of a real upper Hessenberg matrix 
 H(m) = SH'(m)S^*, where
       (h(1,1)   h(1,2)                            )
       (h(2,1)   h(2,2)                            )
       (  0      h(3,2)                            )
   H = (           0   ...                         )
       (                      h(m-2,m-1)   h(m-2,m))
       (                      h(m-1,m-1)   h(m-1,m))                       
       (                   0   h(m,m-1)      h(m)  )
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
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SUBSPACE(=%d) is less than 0\n",ss);
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
  LIS_MATRIX A,B;
  LIS_INT ss,ic;
  LIS_SCALAR gshift;
  LIS_INT emaxiter,iter0,hqriter;
  LIS_REAL tol,hqrerr,D;
  LIS_INT i,j;
  LIS_INT output, niesolver;
  LIS_REAL nrm2,resid0; 
  LIS_VECTOR *v,w;
  LIS_SCALAR *h,*hq,*hr,evalue,evalue0;
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

  h = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eai::h");
  hq = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eai::hq");
  hr = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eai::hr");
  
  A = esolver->A;
  B = esolver->B;
  w = esolver->work[0];
  v = &esolver->work[1];
  lis_vector_set_all(0.0,v[0]);
  lis_vector_set_all(1.0,w);
  lis_vector_nrm2(w, &nrm2);

  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p none",solver);  
  lis_solver_set_optionC(solver);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);
  lis_esolver_get_esolvername(niesolver, esolvername);
  if( A->my_rank==0 ) printf("inner eigensolver     : %s\n", esolvername);
  if( A->my_rank==0 ) printf("linear solver         : %s\n", solvername);
  if( A->my_rank==0 ) printf("preconditioner        : %s\n", preconname);

  for (i=0;i<ss*ss;i++) h[i] = 0.0;

  j=-1;
  while (j<ss-1)
    {
      j = j+1;
      lis_vector_copy(w, v[j]);

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

      /* v(j+1) = w / h(i+1,j) */
      lis_vector_scale(1/h[j+1+j*ss],w);
      lis_vector_copy(w,v[j+1]);
      
    }

  /* compute eigenvalues of a real upper
     Hessenberg matrix H(j) = SH'(j)S^* */
  lis_array_qr(ss,h,hq,hr,&hqriter,&hqrerr);


  if( A->my_rank==0 ) 
    {
#ifdef _LONG__LONG
      printf("size of subspace      : %lld\n\n", ss);
#else
      printf("size of subspace      : %d\n\n", ss);
#endif
      if( output ) printf("approximate eigenvalues in subspace:\n\n");


      i=0;
      while (i<ss-1) 
	{
	  i = i + 1;
	  if (fabs(h[i+(i-1)*ss])<tol)
	    {
#ifdef _LONG__LONG
	      printf("Arnoldi: mode number              = %lld\n",i-1);
#else
	      printf("Arnoldi: mode number              = %d\n",i-1);
#endif
#ifdef _LONG__DOUBLE
	      printf("Arnoldi: eigenvalue               = %Le\n",(LIS_REAL)(h[i-1+(i-1)*ss] - gshift));
#else
	      printf("Arnoldi: eigenvalue               = %e\n",(LIS_REAL)(h[i-1+(i-1)*ss] - gshift));
#endif
	      esolver->evalue[i-1] = h[i-1+(i-1)*ss] - gshift;
	    }
	  else
	    {
	      D = (h[i-1+(i-1)*ss]-h[i+i*ss]) * (h[i-1+(i-1)*ss]-h[i+i*ss])
		+ 4 * h[i-1+i*ss] * h[i+(i-1)*ss];
	      if (D<0)
		{
#ifdef _LONG__LONG
		  printf("Arnoldi: mode number              = %lld\n",i-1);
#else
		  printf("Arnoldi: mode number              = %d\n",i-1);
#endif
#ifdef _LONG__DOUBLE	      
		  printf("Arnoldi: eigenvalue               = (%Le, %Le)\n", (LIS_REAL)((h[i-1+(i-1)*ss]+h[i+i*ss])/2 - gshift), (LIS_REAL)sqrt(-D)/2);
#else
		  printf("Arnoldi: eigenvalue               = (%e, %e)\n", (LIS_REAL)((h[i-1+(i-1)*ss]+h[i+i*ss])/2 - gshift), (LIS_REAL)sqrt(-D)/2);
#endif
#ifdef _LONG__LONG	      
		  printf("Arnoldi: mode number              = %lld\n",i);
#else
		  printf("Arnoldi: mode number              = %d\n",i);
#endif
#ifdef _LONG__DOUBLE	      	      
		  printf("Arnoldi: eigenvalue               = (%Le, %Le)\n", (LIS_REAL)((h[i-1+(i-1)*ss]+h[i+i*ss])/2 - gshift), (LIS_REAL)-sqrt(-D)/2);
#else
		  printf("Arnoldi: eigenvalue               = (%e, %e)\n", (LIS_REAL)((h[i-1+(i-1)*ss]+h[i+i*ss])/2 - gshift), (LIS_REAL)-sqrt(-D)/2);
#endif
		  
#ifdef _COMPLEX		  
		  esolver->evalue[i-1] = (h[i-1+(i-1)*ss]+h[i+i*ss])/2 + sqrt(-D)/2 * _Complex_I - gshift;
		  esolver->evalue[i]   = (h[i-1+(i-1)*ss]+h[i+i*ss])/2 - sqrt(-D)/2 * _Complex_I - gshift;
#else
		  esolver->evalue[i-1] = (h[i-1+(i-1)*ss]+h[i+i*ss])/2 - gshift;
		  esolver->evalue[i]   = (h[i-1+(i-1)*ss]+h[i+i*ss])/2 - gshift;     
#endif
		  
		  i=i+1;
		}
	      else
		{
#ifdef _LONG__LONG	      	      
		  printf("Arnoldi: mode number              = %lld\n",i-1);
#else
		  printf("Arnoldi: mode number              = %d\n",i-1);
#endif
#ifdef _LONG__DOUBLE	      	      	      
		  printf("Arnoldi: eigenvalue               = %Le\n",(LIS_REAL)(h[i-1+(i-1)*ss] - gshift));
#else
		  printf("Arnoldi: eigenvalue               = %e\n",(LIS_REAL)(h[i-1+(i-1)*ss] - gshift));
#endif	      
		  esolver->evalue[i-1] = h[i-1+(i-1)*ss] - gshift;
		}
	    }
	}
      if (i<ss)
	{
#ifdef _LONG__LONG	            
	  printf("Arnoldi: mode number              = %lld\n",i);
#else
	  printf("Arnoldi: mode number              = %d\n",i);
#endif
#ifdef _LONG__DOUBLE	      	      	      	      
	  printf("Arnoldi: eigenvalue               = %Le\n",(LIS_REAL)(h[i+i*ss] - gshift));
#else
	  printf("Arnoldi: eigenvalue               = %e\n",(LIS_REAL)(h[i+i*ss] - gshift));
#endif	      
	}

      if( output ) printf("\n");
      if( output ) printf("compute refined eigenpairs:\n\n");
  
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
      esolver2->lshift = -(esolver->evalue[i]);
      lis_gesolve(A, B, esolver->evector[i], &evalue, esolver2);
      lis_esolver_work_destroy(esolver2); 
      esolver->evalue[i] = evalue - esolver2->lshift - gshift;
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

      if (A->my_rank==0) 
	{

#ifdef _LONG__LONG
	  if( output ) printf("Arnoldi: mode number          = %lld\n", i);
#else
	  if( output ) printf("Arnoldi: mode number          = %d\n", i);
#endif
#ifdef _COMPLEX	  
#ifdef _LONG__DOUBLE
	  if( output ) printf("Arnoldi: eigenvalue           = (%Le, %Le)\n", creall(esolver->evalue[i]), cimagl(esolver->evalue[i]));
#else
	  if( output ) printf("Arnoldi: eigenvalue           = (%e, %e)\n", creal(esolver->evalue[i]), cimag(esolver->evalue[i]));
#endif
#else
#ifdef _LONG__DOUBLE
	  if( output ) printf("Arnoldi: eigenvalue           = %Le\n", (LIS_REAL)esolver->evalue[i]);
#else
	  if( output ) printf("Arnoldi: eigenvalue           = %e\n", (LIS_REAL)esolver->evalue[i]);
#endif
#endif	  
#ifdef _LONG__LONG
	  if( output ) printf("Arnoldi: number of iterations = %lld\n",esolver2->iter[0]);
#else
	  if( output ) printf("Arnoldi: number of iterations = %d\n",esolver2->iter[0]);
#endif
#ifdef _LONG__DOUBLE
	  if( output ) printf("Arnoldi: relative residual    = %Le\n\n",esolver2->resid[0]);
#else
	  if( output ) printf("Arnoldi: relative residual    = %e\n\n",esolver2->resid[0]);
#endif
	}
    }

  esolver->evalue[0] = evalue0; 
  esolver->iter[0] = iter0;
  esolver->resid[0] = resid0;

  lis_vector_copy(esolver->evector[0], esolver->x);

  lis_esolver_destroy(esolver2); 

  lis_free(h); 
  lis_free(hq);
  lis_free(hr);

  lis_solver_destroy(solver);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
} 

