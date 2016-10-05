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
 * Approximate Inverse Iteration       *
 ***************************************
 x(0)    = (1,...,1)^T
 ***************************************
 for k=1,2,...
   x(k-1)    = x(k-1)/||x(k-1)||_2
   z         = (M - lshift I)^-1 * x(k-1)
   1/evalue  = <x(k-1),z>
   resid     = ||z - 1/evalue * x||_2 / |1/evalue|
   x(k)      = z         
 ***************************************/

#define NWORK 2
#undef __FUNC__
#define __FUNC__ "lis_eaii_check_params"
LIS_INT lis_eaii_check_params(LIS_ESOLVER esolver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_eaii_malloc_work"
LIS_INT lis_eaii_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_eaii_malloc_work::work" );
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
#define __FUNC__ "lis_eaii"
LIS_INT lis_eaii(LIS_ESOLVER esolver)
{
  LIS_MATRIX A;
  LIS_VECTOR x;
  LIS_SCALAR evalue, ievalue;
  LIS_SCALAR lshift;
  LIS_INT emaxiter;
  LIS_REAL tol;
  LIS_INT iter,output;
  LIS_REAL nrm2,resid;
  LIS_VECTOR z,q;
  LIS_SOLVER solver;
  LIS_PRECON precon;
  double time,itime,ptime,p_c_time,p_i_time;
  LIS_INT nsol, precon_type;
  char solvername[128], preconname[128];

  LIS_DEBUG_FUNC_IN;

  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  lshift = esolver->lshift;
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];

  A = esolver->A;
  x = esolver->x;
  if (esolver->options[LIS_EOPTIONS_INITGUESS_ONES] ) 
    {
      lis_vector_set_all(1.0,x);
    }
  evalue = 1.0;
  z = esolver->work[0];
  q = esolver->work[1];
  ptime = 0;

  iter=0;
  ievalue = 1/(evalue);
#ifdef _COMPLEX  
#ifdef _LONG__DOUBLE
  if( output & (A->my_rank==0) ) printf("local shift           : (%Le, %Le)\n", creall(lshift), cimagl(lshift));
#else
  if( output & (A->my_rank==0) ) printf("local shift           : (%e, %e)\n", creal(lshift), cimag(lshift));
#endif
#else
#ifdef _LONG__DOUBLE
  if( output & (A->my_rank==0) ) printf("local shift           : %Le\n", lshift);
#else
  if( output & (A->my_rank==0) ) printf("local shift           : %e\n", lshift);
#endif
#endif  
  if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p ilu",solver);
  lis_solver_set_optionC(solver);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);
  if( output & (A->my_rank==0) ) printf("linear solver         : %s\n", solvername);
  if( output & (A->my_rank==0) ) printf("preconditioner        : %s\n", preconname);

  lis_vector_set_all(1.0,q);
  lis_solve(A, q, x, solver);
  lis_precon_create(solver, &precon);
  solver->precon = precon;

  while (iter<emaxiter)
    {
      iter = iter+1;

      /* x = x / ||x||_2 */
      lis_vector_nrm2(x, &nrm2);
      lis_vector_scale(1/nrm2, x);

      /* z = (M - lshift I)^-1 * x */
      time = lis_wtime();
      lis_psolve(solver, x, z);
      ptime += lis_wtime() - time;

      /* 1/evalue = <x,z> */
      lis_vector_dot(x, z, &ievalue); 

      /* resid = ||z - 1/evalue * x||_2 / |1/evalue| */
      lis_vector_axpyz(-ievalue,x,z,q); 
      lis_vector_nrm2(q, &resid); 
      resid = fabs(resid/ievalue);

      if( output )
	{
	  if( output & LIS_EPRINT_MEM ) esolver->rhistory[iter] = resid;
	  if( output & LIS_EPRINT_OUT && A->my_rank==0 ) lis_print_rhistory(iter,resid);
	}

      /* x = z */
      lis_vector_copy(z,x);

      /* convergence check */      
      if( tol >= resid ) 
	{
	  esolver->retcode    = LIS_SUCCESS;
	  esolver->iter[0]    = iter;
	  esolver->resid[0]   = resid;
	  esolver->evalue[0]  = 1/ievalue;
	  lis_vector_nrm2(x, &nrm2);
	  lis_vector_scale(1/nrm2, x);
	  if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
	  lis_precon_destroy(precon);
	  lis_solver_destroy(solver);
	  LIS_DEBUG_FUNC_OUT;
	  return LIS_SUCCESS;
	}
    }
  esolver->retcode    = LIS_MAXITER;
  esolver->iter[0]    = iter;
  esolver->resid[0]   = resid;
  esolver->evalue[0]  = 1/ievalue;
  lis_vector_nrm2(x, &nrm2);
  lis_vector_scale(1/nrm2, x);

  lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
  esolver->itime = solver->itime;
  esolver->ptime = ptime;
  esolver->p_i_time = solver->p_i_time;
  esolver->p_c_time = solver->p_c_time;

  if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
  lis_precon_destroy(precon);
  lis_solver_destroy(solver);
  LIS_DEBUG_FUNC_OUT;
  return LIS_MAXITER;
}
