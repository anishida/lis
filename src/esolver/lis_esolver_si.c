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
 * Subspace Iteration                  *
 ***************************************
 QR factorization V * R = Z for the starting matrix Z
 for k=1,2,...
   if Power Iteration
     R = A * V
   if Inverse Iteration
     R = A^-1 * V
   Z = V*R
   resid     = ||Z - V * R||_2
   QR factorization V * R = Z
 ***************************************/

#define NWORK 4
#undef __FUNC__
#define __FUNC__ "lis_esi_check_params"
LIS_INT lis_esi_check_params(LIS_ESOLVER esolver)
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
#define __FUNC__ "lis_esi_malloc_work"
LIS_INT lis_esi_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err,ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

	worklen = NWORK + ss;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_esi_malloc_work::work" );
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
#define __FUNC__ "lis_esi"
LIS_INT lis_esi(LIS_ESOLVER esolver)
{
  LIS_Comm comm;  
  LIS_MATRIX A;
  LIS_VECTOR x, Ax;
  LIS_SCALAR xAx, xx;
  LIS_SCALAR oshift;
  LIS_INT ss;
  LIS_INT emaxiter;
  LIS_REAL tol;
  LIS_INT j,k;
  LIS_SCALAR theta;
  LIS_INT iter,giter,output,niesolver;
  LIS_REAL nrm2,resid;
  LIS_SCALAR dot;
  LIS_VECTOR *v,r,q;
  LIS_SOLVER solver;
  LIS_PRECON precon;
  double time,itime,ptime,p_c_time,p_i_time;
  double etime0,etime;
  LIS_INT err;
  LIS_INT nsol, precon_type;
  char solvername[128], preconname[128], esolvername[128];

  LIS_DEBUG_FUNC_IN;

  comm = LIS_COMM_WORLD;

  A = esolver->A;
  x = esolver->x;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN];
#ifdef _COMPLEX
  oshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN] + esolver->params[LIS_EPARAMS_SHIFT_IM - LIS_EOPTIONS_LEN] * _Complex_I;
#else
  oshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN];
#endif	
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];

  r = esolver->work[0];
  q = esolver->work[1];
  v = &esolver->work[2];
  Ax = esolver->work[3];
  lis_vector_set_all(1.0,r);
  lis_vector_nrm2(r, &nrm2);
  lis_vector_scale(1.0/nrm2,r);

  if ( oshift != 0.0 ) lis_matrix_shift_diagonal(A, oshift);  

  if( output )
    {
#ifdef _COMPLEX
      lis_printf(comm,"shift                 : (%e, %e)\n", (double)creal(oshift), (double)cimag(oshift));      
#else
      lis_printf(comm,"shift                 : %e\n", (double)oshift);
#endif
    }
  
  lis_esolver_get_esolvername(niesolver, esolvername);
  if( output ) lis_printf(comm,"inner eigensolver     : %s\n", esolvername);

  switch ( niesolver )
    {
    case LIS_ESOLVER_II:

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
      break;

    }

  if( output ) 
    {
      lis_printf(comm,"size of subspace      : %D\n\n", ss);
      lis_printf(comm,"compute eigenpairs in subspace:\n\n");
    }

  giter=0;
  j=0;
  while (j<ss)
    {
      etime0 = lis_wtime();
      lis_vector_duplicate(A,&esolver->evector[j]); 
      j = j+1;
      lis_vector_copy(r, v[j]);

      if (niesolver==LIS_ESOLVER_II )
	{
	  /* create preconditioner */
	  solver->A = A;
	  err = lis_precon_create(solver, &precon);
	  if( err )
	    {
	      lis_solver_work_destroy(solver);
	      solver->retcode = err;
	      return err;
	    }
	}

      ptime = 0;
      iter = 0;
      while (iter<emaxiter)
	{
	  iter = iter+1;
	  giter = giter+1;

	  /* QR factorization V*R = Z for starting vector Z */
	  for (k=1;k<j;k++)
	    {
	      lis_vector_dot(v[j], v[k], &dot); 
	      lis_vector_axpy(-dot, v[k], v[j]);
	    }

	  /* kernel */
	  switch( niesolver )
	    {
	    case LIS_ESOLVER_PI:

	      /* R = A * V */
	      lis_matvec(A,v[j],r); 

	      break;

	    case LIS_ESOLVER_II:

	      /* R = A^-1 * V */
	      err = lis_solve_kernel(A, v[j], r, solver, precon);
	      if( err )
		{
		  lis_solver_work_destroy(solver);	  
		  solver->retcode = err;
		  return err;
		}

	      break;

	    }

	  /* elapsed time of linear solver */
	  if ( j==1 &&  niesolver==LIS_ESOLVER_II )
	    {
	      lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
	      esolver->ptime += solver->ptime;
	      esolver->itime += solver->itime;
	      esolver->p_c_time += solver->p_c_time;
	      esolver->p_i_time += solver->p_i_time;
	    }

	  /* Z = V*R */
	  lis_vector_nrm2(r, &nrm2);
	  lis_vector_dot(v[j],r,&theta);

	  /* resid = ||Z - V*R||_2 */
	  lis_vector_axpyz(-theta,v[j],r,q);
	  lis_vector_nrm2(q, &resid);
	  resid = resid / fabs(theta);

	  lis_vector_scale(1.0/nrm2,r);
	  lis_vector_copy(r, v[j]);

	  /* convergence check */
	  if ( j==1 ) 
	    {
	      if( output & LIS_PRINT_MEM ) esolver->rhistory[iter] = resid; 
	      esolver->iter[j-1] = iter;
	      /* esolver->resid = resid; */
	    }


	  if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,resid);
	  if (tol>resid) break;
	}

      switch ( niesolver )
	{
	case LIS_ESOLVER_PI:
	  esolver->evalue[j-1] = theta + oshift;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	case LIS_ESOLVER_II:
	  esolver->evalue[j-1] = 1/theta + oshift;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	}

      lis_vector_copy(v[j], esolver->evector[j-1]);
      etime = lis_wtime() - etime0;

      if (output & (ss>1))
	{
	  lis_printf(comm,"Subspace: mode number          = %D\n", j-1);
#ifdef _COMPLEX	  
	  lis_printf(comm,"Subspace: eigenvalue           = (%e, %e)\n", (double)creal(esolver->evalue[j-1]), (double)cimag(esolver->evalue[j-1]));
#else
	  lis_printf(comm,"Subspace: eigenvalue           = %e\n", (double)esolver->evalue[j-1]);
#endif
	  lis_printf(comm,"Subspace: elapsed time         = %e sec.\n", etime);	  
	  lis_printf(comm,"Subspace: number of iterations = %D\n",iter);
	  lis_printf(comm,"Subspace: relative residual    = %e\n\n",(double)resid);
	}
    }
  
  if ( oshift != 0.0 ) lis_matrix_shift_diagonal(A, -oshift);  
  lis_vector_copy(esolver->evector[0], esolver->x);

  switch ( niesolver )
    {
    case LIS_ESOLVER_II:
      lis_precon_destroy(precon);
      lis_solver_destroy(solver);
      break;
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

/***************************************
 * Generalized Subspace Iteration      *
 ***************************************
 QR factorization V * R = Z for the starting matrix Z
 for k=1,2,...
   if Power Iteration
     R = B^-1 * A * V
   if Inverse Iteration
     R = A^-1 * B * V
   Z = V * R
   resid = ||Z - V * R||_2
   QR factorization V * R = Z
 ***************************************/

#undef NWORK
#define NWORK 5
#undef __FUNC__
#define __FUNC__ "lis_egsi_check_params"
LIS_INT lis_egsi_check_params(LIS_ESOLVER esolver)
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
#define __FUNC__ "lis_egsi_malloc_work"
LIS_INT lis_egsi_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err,ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

	worklen = NWORK + ss;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_egsi_malloc_work::work" );
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
#define __FUNC__ "lis_egsi"
LIS_INT lis_egsi(LIS_ESOLVER esolver)
{
  LIS_Comm comm;  
  LIS_MATRIX A, B;
  LIS_VECTOR x, Ax;
  LIS_SCALAR xAx, xx;
  LIS_SCALAR oshift;
  LIS_INT ss;
  LIS_INT emaxiter;
  LIS_REAL tol;
  LIS_INT j,k;
  LIS_SCALAR theta,eta;
  LIS_INT iter,giter,output,nigesolver;
  LIS_REAL nrm2,resid;
  LIS_SCALAR dot;
  LIS_VECTOR *v,w,y,q;
  LIS_SOLVER solver;
  LIS_PRECON precon;
  double time,itime,ptime,p_c_time,p_i_time;
  double etime0,etime;  
  LIS_INT err;
  LIS_INT nsol, precon_type;
  char solvername[128], preconname[128], esolvername[128];

  LIS_DEBUG_FUNC_IN;

  comm = LIS_COMM_WORLD;

  A = esolver->A;
  B = esolver->B;
  x = esolver->x;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN];
#ifdef _COMPLEX
  oshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN] + esolver->params[LIS_EPARAMS_SHIFT_IM - LIS_EOPTIONS_LEN] * _Complex_I;
#else
  oshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN];
#endif	
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  nigesolver = esolver->options[LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER];

  y = esolver->work[0];
  w = esolver->work[1];
  q = esolver->work[2];
  v = &esolver->work[3];
  Ax = esolver->work[4];
  lis_vector_set_all(1.0,y);
  lis_vector_nrm2(y,&nrm2);
  lis_vector_scale(1.0/nrm2,y);

  if ( oshift != 0.0 ) lis_matrix_shift_matrix(A, B, oshift);

  if( output )
    {
#ifdef _COMPLEX
      lis_printf(comm,"shift                 : (%e, %e)\n", (double)creal(oshift), (double)cimag(oshift));      
#else
      lis_printf(comm,"shift                 : %e\n", (double)oshift);
#endif
    }
  
  lis_esolver_get_esolvername(nigesolver, esolvername);
  if( output ) lis_printf(comm,"inner eigensolver     : %s\n", esolvername);

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

  if( output ) 
    {
      lis_printf(comm,"size of subspace      : %D\n\n", ss);
      lis_printf(comm,"compute eigenpairs in subspace:\n\n");
    }

  giter=0;
  j=0;
  while (j<ss)
    {
      etime0 = lis_wtime();
      lis_vector_duplicate(A,&esolver->evector[j]); 
      j = j+1;
      lis_vector_copy(y, v[j]);

      switch ( nigesolver )
	{
	case LIS_ESOLVER_GII:
	  {
	    /* create preconditioner */
	    solver->A = A;
	    err = lis_precon_create(solver, &precon);
	    if( err )
	      {
		lis_solver_work_destroy(solver);
		solver->retcode = err;
		return err;
	      }
	  }
	case LIS_ESOLVER_GPI:
	  {
	    /* create preconditioner */
	    solver->A = B;
	    err = lis_precon_create(solver, &precon);
	    if( err )
	      {
		lis_solver_work_destroy(solver);
		solver->retcode = err;
		return err;
	      }
	  }
	}

      ptime = 0;
      iter = 0;
      while (iter<emaxiter)
	{
	  iter = iter+1;
	  giter = giter+1;

	  /* QR factorization V * R = Z for starting vector Z */
	  for (k=1;k<j;k++)
	    {
	      lis_vector_dot(v[j], v[k], &dot);
	      lis_vector_axpy(-dot, v[k], v[j]);
	    }

	  /* kernel */
	  switch( nigesolver )
	    {
	    case LIS_ESOLVER_GPI:

	      /* y = A * v */	      
	      lis_matvec(A, v[j], w); 

	      /* v = v / <v,w>^1/2, w = w / <v,w>^1/2 */	      
	      lis_vector_dot(v[j], w, &eta);
	      eta = sqrt(eta);
	      lis_vector_scale(1.0/eta, v[j]);
	      lis_vector_scale(1.0/eta, w);

	      /* y = B^-1 * w */
	      err = lis_solve_kernel(B, w, y, solver, precon);
	      if( err )
		{
		  lis_solver_work_destroy(solver);	  
		  solver->retcode = err;
		  return err;
		}

	      break;

	    case LIS_ESOLVER_GII:

	      /* w = B * v */
	      lis_matvec(B, v[j], w);
	      
	      /* v = v / <v,w>^1/2, w = w / <v,w>^1/2 */
	      lis_vector_dot(v[j], w, &eta);
	      eta = sqrt(eta);
	      lis_vector_scale(1.0/eta, v[j]);
	      lis_vector_scale(1.0/eta, w);

	      /* y = A^-1 * w */
	      err = lis_solve_kernel(A, w, y, solver, precon);
	      if( err )
		{
		  lis_solver_work_destroy(solver);	  
		  solver->retcode = err;
		  return err;
		}
	      break;

	    }

	  /* elapsed time of linear solver */
	  if ( j==1 )
	    {
	      lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
	      esolver->ptime += solver->ptime;
	      esolver->itime += solver->itime;
	      esolver->p_c_time += solver->p_c_time;
	      esolver->p_i_time += solver->p_i_time;
	    }

	  /* theta = <w,y> */
	  lis_vector_dot(w, y, &theta);

	  /* resid = ||y - theta * v||_2 / |theta| */
	  lis_vector_axpyz(-theta, v[j], y, q);
	  lis_vector_nrm2(q, &resid);
	  resid = resid / fabs(theta);

	  /* y = y / ||y||_2 */
	  lis_vector_nrm2(y, &nrm2);
	  lis_vector_scale(1.0/nrm2, y);
	  
	  /* v = y */
	  lis_vector_copy(y, v[j]);

	  /* convergence check */
	  if ( j==1 ) 
	    {
	      if( output & LIS_PRINT_MEM ) esolver->rhistory[iter] = resid; 
	      esolver->iter[j-1] = iter;
	      /* esolver->resid = resid; */
	    }


	  if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,resid);
	  if (tol>resid) break;
	}

      switch ( nigesolver )
	{
	case LIS_ESOLVER_GPI:
	  esolver->evalue[j-1] = theta + oshift;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	case LIS_ESOLVER_GII:
	  esolver->evalue[j-1] = 1.0/theta + oshift;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	}

      lis_vector_copy(v[j], esolver->evector[j-1]);
      etime = lis_wtime() - etime0;

      if (output & (ss>1))
	{
	  lis_printf(comm,"Generalized Subspace: mode number          = %D\n", j-1);
#ifdef _COMPLEX	  
	  lis_printf(comm,"Generalized Subspace: eigenvalue           = (%e, %e)\n", (double)creal(esolver->evalue[j-1]), (double)cimag(esolver->evalue[j-1]));
#else
	  lis_printf(comm,"Generalized Subspace: eigenvalue           = %e\n", (double)esolver->evalue[j-1]);
#endif
	  lis_printf(comm,"Generalized Subspace: elapsed time         = %e sec.\n", etime);	  	  
	  lis_printf(comm,"Generalized Subspace: number of iterations = %D\n",iter);
	  lis_printf(comm,"Generalized Subspace: relative residual    = %e\n\n",(double)resid);
	}
    }
  
  if ( oshift != 0.0 ) lis_matrix_shift_matrix(A, B, -oshift);  
  lis_vector_copy(esolver->evector[0], esolver->x);

  lis_precon_destroy(precon);
  lis_solver_destroy(solver);
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}
