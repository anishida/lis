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
 * Inverse Iteration                   *
 ***************************************
 v    = (1,...,1)^T
 ***************************************
 for k=1,2,...
   v     = v / ||v||_2
   y     = A^-1 * v
   theta = <v,y>
   resid = ||y - <v,y> * v||_2 / |<v,y>|
   v     = y
   if resid < tol then stop
 lambda  = 1 / theta
 x       = v / ||v||_2
 ***************************************/

#define NWORK 2
#undef __FUNC__
#define __FUNC__ "lis_eii_check_params"
LIS_INT lis_eii_check_params(LIS_ESOLVER esolver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_eii_malloc_work"
LIS_INT lis_eii_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_eii_malloc_work::work" );
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
#define __FUNC__ "lis_eii"
LIS_INT lis_eii(LIS_ESOLVER esolver)
{
  LIS_Comm comm;  
  LIS_MATRIX A;
  LIS_VECTOR v;
  LIS_SCALAR theta;  
  LIS_SCALAR oshift,ishift;
  LIS_INT emaxiter;
  LIS_REAL tol;
  LIS_INT iter,iter2,output;
  LIS_REAL nrm2,resid;
  LIS_VECTOR y,q;
  LIS_SOLVER solver;
  double time,itime,ptime,p_c_time,p_i_time;

  LIS_INT err;
  LIS_PRECON precon;
  LIS_INT nsol, precon_type;
  char solvername[128], preconname[128];

  LIS_DEBUG_FUNC_IN;

  comm = LIS_COMM_WORLD;

  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN];
#ifdef _COMPLEX
  oshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN] + esolver->params[LIS_EPARAMS_SHIFT_IM - LIS_EOPTIONS_LEN] * _Complex_I;
#else
  oshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN];
#endif	
  ishift = esolver->ishift;
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];

  A = esolver->A;
  v = esolver->x;
  if (esolver->options[LIS_EOPTIONS_INITGUESS_ONES] ) 
    {
      lis_vector_set_all(1.0,v);
    }
  y = esolver->work[0];
  q = esolver->work[1];

  if ( esolver->ishift != 0.0 ) oshift = ishift;
  if ( oshift != 0.0 ) lis_matrix_shift_diagonal(A, oshift);

  if( output )
    {
#ifdef _COMPLEX
      lis_printf(comm,"shift                 : (%e, %e)\n", (double)creal(oshift), (double)cimag(oshift));      
#else
      lis_printf(comm,"shift                 : %e\n", (double)oshift);
#endif
    }
  
  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p none",solver);
  err = lis_solver_set_optionC(solver);
  CHKERR(err);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);
  if( output )
    {
      lis_printf(comm,"linear solver         : %s\n", solvername);
      lis_printf(comm,"preconditioner        : %s\n", preconname);
    }

  /* create preconditioner */
  solver->A = A;
  err = lis_precon_create(solver, &precon);
  if( err )
    {
      lis_solver_work_destroy(solver);
      solver->retcode = err;
      return err;
    }

  iter=0;
  while (iter<emaxiter)
    {
      iter = iter+1;

      /* v = v / ||v||_2 */
      lis_vector_nrm2(v, &nrm2);
      lis_vector_scale(1.0/nrm2, v);

      /* y = A^-1 * v */
      err = lis_solve_kernel(A, v, y, solver, precon);
      if( err )
	{
	  lis_solver_work_destroy(solver);	  
	  solver->retcode = err;
	  return err;
	}
      lis_solver_get_iter(solver, &iter2);

      /* theta = <v,y> */
      lis_vector_dot(v, y, &theta);

      /* resid = ||y - theta * v||_2 / |theta| */
      lis_vector_axpyz(-theta, v, y, q); 
      lis_vector_nrm2(q, &resid); 
      resid = resid / fabs(theta);

      /* v = y */
      lis_vector_copy(y, v);

      /* convergence check */
      lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
      esolver->ptime += solver->ptime;
      esolver->itime += solver->itime;
      esolver->p_c_time += solver->p_c_time;
      esolver->p_i_time += solver->p_i_time;

      if( output )
	{
	  if( output & LIS_EPRINT_MEM ) esolver->rhistory[iter] = resid;
	  if( output & LIS_EPRINT_OUT ) lis_print_rhistory(comm,iter,resid);
	}

      if( tol >= resid ) 
	{
	  esolver->retcode    = LIS_SUCCESS;
	  esolver->iter[0]    = iter;
	  esolver->resid[0]   = resid;
	  esolver->evalue[0]  = 1.0/theta + oshift;
	  lis_vector_nrm2(v, &nrm2);
	  lis_vector_scale(1.0/nrm2, v);
	  if ( oshift != 0.0 ) lis_matrix_shift_diagonal(A, -oshift);
	  lis_precon_destroy(precon);
	  lis_solver_destroy(solver); 
	  LIS_DEBUG_FUNC_OUT;
	  return LIS_SUCCESS;
	}
    }

  lis_precon_destroy(precon);

  esolver->retcode    = LIS_MAXITER;
  esolver->iter[0]    = iter;
  esolver->resid[0]   = resid;
  esolver->evalue[0]  = 1.0/theta + oshift;
  lis_vector_nrm2(v, &nrm2);
  lis_vector_scale(1.0/nrm2, v);
  if ( oshift != 0.0 ) lis_matrix_shift_diagonal(A, -oshift);
  lis_solver_destroy(solver); 
  LIS_DEBUG_FUNC_OUT;
  return LIS_MAXITER;
}


/***************************************
 * Generalized Inverse Iteration       *
 ***************************************
 v    = (1,...,1)^T
 ***************************************
 for k=1,2,...
   v     = v / ||v||_2
   w     = B * v 
   v     = v / <v,w>^1/2
   w     = w / <v,w>^1/2 
   y     = A^-1 * w
   theta = <w,y>
   resid = ||y - <v,y> * v||_2 / |<v,y>|
   v     = y
   if resid < tol then stop
 lambda  = 1 / theta
 x       = v / ||v||_2
 ***************************************/

#undef NWORK
#define NWORK 3
#undef __FUNC__
#define __FUNC__ "lis_egii_check_params"
LIS_INT lis_egii_check_params(LIS_ESOLVER esolver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_egii_malloc_work"
LIS_INT lis_egii_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_egii_malloc_work::work" );
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
#define __FUNC__ "lis_egii"
LIS_INT lis_egii(LIS_ESOLVER esolver)
{
  LIS_Comm comm;  
  LIS_MATRIX A,B;
  LIS_VECTOR v;
  LIS_SCALAR eta,theta;
  LIS_SCALAR oshift,ishift;
  LIS_INT emaxiter;
  LIS_REAL tol;
  LIS_INT iter,iter2,output;
  LIS_REAL nrm2,resid;
  LIS_VECTOR w,y,q;
  LIS_SOLVER solver;
  double time,itime,ptime,p_c_time,p_i_time;

  LIS_INT err;
  LIS_PRECON precon;
  LIS_INT nsol, precon_type;
  char solvername[128], preconname[128];

  LIS_DEBUG_FUNC_IN;

  comm = LIS_COMM_WORLD;

  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN];
#ifdef _COMPLEX
  oshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN] + esolver->params[LIS_EPARAMS_SHIFT_IM - LIS_EOPTIONS_LEN] * _Complex_I;
#else
  oshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN];
#endif	
  ishift = esolver->ishift;
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];

  A = esolver->A;
  B = esolver->B;  
  v = esolver->x;
  if (esolver->options[LIS_EOPTIONS_INITGUESS_ONES] ) 
    {
      lis_vector_set_all(1.0,v);
    }
  w = esolver->work[0];
  y = esolver->work[1];
  q = esolver->work[2];

  if ( esolver->ishift != 0.0 ) oshift = ishift;
  if ( oshift != 0.0 ) lis_matrix_shift_matrix(A, B, oshift);

  if( output )
    {
#ifdef _COMPLEX
      lis_printf(comm,"shift                 : (%e, %e)\n", (double)creal(oshift), (double)cimag(oshift));
#else
      lis_printf(comm,"shift                 : %e\n", (double)oshift);
#endif
    }
  
  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p none",solver);
  err = lis_solver_set_optionC(solver);
  CHKERR(err);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);
  if( output )
    {
      lis_printf(comm,"linear solver         : %s\n", solvername);
      lis_printf(comm,"preconditioner        : %s\n", preconname);
    }

  /* create preconditioner */
  solver->A = A;
  err = lis_precon_create(solver, &precon);
  if( err )
    {
      lis_solver_work_destroy(solver);
      solver->retcode = err;
      return err;
    }

  iter=0;
  while (iter<emaxiter)
    {
      iter = iter+1;

      /* v = v / ||v||_2 */
      lis_vector_nrm2(v, &nrm2);
      lis_vector_scale(1.0/nrm2, v);

      /* w = B * v */
      lis_matvec(B, v, w);

      /* v = v / <w,v>^1/2, w = w / <w,v>^1/2 */
      lis_vector_dot(w, v, &eta);
      eta = sqrt(eta);
      lis_vector_scale(1.0/eta, v);
      lis_vector_scale(1.0/eta, w);

      /* y = A^-1 * w */
      err = lis_solve_kernel(A, w, y, solver, precon);
      if( err )
	{
	  lis_solver_work_destroy(solver);	  
	  solver->retcode = err;
	  return err;
	}
      lis_solver_get_iter(solver,&iter2);

      /* theta = <w,y> */
      lis_vector_dot(w, y, &theta);

      /* resid = ||y - theta * v||_2 / |theta| */

      lis_vector_axpyz(-theta, v, y, q); 
      lis_vector_nrm2(q, &resid); 
      resid = resid / fabs(theta);

      /* v = y */
      lis_vector_copy(y, v);

      /* convergence check */
      lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
      esolver->ptime += solver->ptime;
      esolver->itime += solver->itime;
      esolver->p_c_time += solver->p_c_time;
      esolver->p_i_time += solver->p_i_time;

      if( output )
	{
	  if( output & LIS_EPRINT_MEM ) esolver->rhistory[iter] = resid;
	  if( output & LIS_EPRINT_OUT ) lis_print_rhistory(comm,iter,resid);
	}

      if( tol >= resid ) 
	{
	  esolver->retcode    = LIS_SUCCESS;
	  esolver->iter[0]    = iter;
	  esolver->resid[0]   = resid;
	  esolver->evalue[0]  = 1.0/theta + oshift;
	  lis_vector_nrm2(v, &nrm2);
	  lis_vector_scale(1.0/nrm2, v);
	  if ( oshift != 0.0 ) lis_matrix_shift_matrix(A, B, -oshift);
	  lis_precon_destroy(precon);
	  lis_solver_destroy(solver); 
	  LIS_DEBUG_FUNC_OUT;
	  return LIS_SUCCESS;
	}
    }

  lis_precon_destroy(precon);

  esolver->retcode    = LIS_MAXITER;
  esolver->iter[0]    = iter;
  esolver->resid[0]   = resid;
  esolver->evalue[0]  = 1.0/theta + oshift;
  lis_vector_nrm2(v, &nrm2);
  lis_vector_scale(1.0/nrm2, v);
  if ( oshift != 0.0 ) lis_matrix_shift_matrix(A, B, -oshift);
  lis_solver_destroy(solver); 
  LIS_DEBUG_FUNC_OUT;
  return LIS_MAXITER;
}

