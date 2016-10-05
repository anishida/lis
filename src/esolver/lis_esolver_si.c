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
 QR factorization VR = Z for the starting matrix Z
 for k=1,2,...
   if Power Iteration
     R = A * V
   if Inverse Iteration
     R = (A - lshift I)^-1 * V
   if Approximate Inverse Iteration
     R = (M - lshift I)^-1 * V
   if Rayleigh Quotient Iteration
     R = (A - mu I)^-1 * V
   R = V*R
   resid     = ||Z - VR||_2
   QR factorization VR = Z
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
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SUBSPACE(=%d) is less than 0\n",ss);
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
  LIS_MATRIX A;
  LIS_VECTOR x, Ax;
  LIS_SCALAR xAx, xx, mu;
  LIS_SCALAR lshift;
  LIS_INT ss;
  LIS_INT emaxiter;
  LIS_REAL tol;
  LIS_INT j,k;
  LIS_SCALAR dotvr;
  LIS_INT iter,giter,output,niesolver;
  LIS_REAL nrm2,resid;
  LIS_SCALAR dot;
  LIS_VECTOR *v,r,q;
  LIS_SOLVER solver;
  LIS_PRECON precon;
  double time,itime,ptime,p_c_time,p_i_time;
  LIS_INT err;
  LIS_INT nsol, precon_type;
  char solvername[128], preconname[128], esolvername[128];

  LIS_DEBUG_FUNC_IN;

  A = esolver->A;
  x = esolver->x;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  lshift = esolver->lshift;
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];

  r = esolver->work[0];
  q = esolver->work[1];
  v = &esolver->work[2];
  Ax = esolver->work[3];
  lis_vector_set_all(1.0,r);
  lis_vector_nrm2(r, &nrm2);
  lis_vector_scale(1/nrm2,r);

  lis_esolver_get_esolvername(niesolver, esolvername);
  if( A->my_rank==0 ) printf("inner eigensolver     : %s\n", esolvername);

  switch ( niesolver )
    {
    case LIS_ESOLVER_II:

      lis_solver_create(&solver);
      lis_solver_set_option("-i bicg -p none",solver);
      lis_solver_set_optionC(solver);
      lis_solver_get_solver(solver, &nsol);
      lis_solver_get_precon(solver, &precon_type);
      lis_solver_get_solvername(nsol, solvername);
      lis_solver_get_preconname(precon_type, preconname);
      if( A->my_rank==0 ) printf("linear solver         : %s\n", solvername);
      if( A->my_rank==0 ) printf("preconditioner        : %s\n", preconname);
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
      break;

    case LIS_ESOLVER_AII:

      lis_solver_create(&solver);
      lis_solver_set_option("-i bicg -p none",solver);
      lis_solver_set_optionC(solver);
      lis_solver_get_solver(solver, &nsol);
      lis_solver_get_precon(solver, &precon_type);
      lis_solver_get_solvername(nsol, solvername);
      lis_solver_get_preconname(precon_type, preconname);
      if( A->my_rank==0 ) printf("linear solver         : %s\n", solvername);
      if( A->my_rank==0 ) printf("preconditioner        : %s\n", preconname);
      if( A->my_rank==0 ) printf("\n");
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
      lis_vector_set_all(1.0,q);
      lis_solve(A, q, x, solver);
      lis_precon_create(solver, &precon);
      solver->precon = precon;
      break;

    case LIS_ESOLVER_RQI:

      lis_solver_create(&solver);
      lis_solver_set_option("-p none -maxiter 10",solver);
      lis_solver_set_optionC(solver);
      lis_solver_get_solver(solver, &nsol);
      lis_solver_get_precon(solver, &precon_type);
      lis_solver_get_solvername(nsol, solvername);
      lis_solver_get_preconname(precon_type, preconname);
      if( A->my_rank==0 ) printf("linear solver         : %s\n", solvername);
      if( A->my_rank==0 ) printf("preconditioner        : %s\n", preconname);
      if( A->my_rank==0 ) printf("\n");
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
      break;

    }

  giter=0;
  j=0;

  if( A->my_rank==0 ) 
    {
#ifdef _LONG__LONG
      printf("size of subspace      : %lld\n\n", ss);
#else
      printf("size of subspace      : %d\n\n", ss);
#endif
      if( output ) printf("compute eigenpairs in subspace:\n\n");
    }

  while (j<ss)
    {
      lis_vector_duplicate(A,&esolver->evector[j]); 
      j = j+1;
      lis_vector_copy(r, v[j]);

      if (niesolver==LIS_ESOLVER_II || niesolver==LIS_ESOLVER_RQI)
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

      if (niesolver==LIS_ESOLVER_RQI)
	{
	  lis_vector_nrm2(x, &nrm2);
	  lis_vector_scale(1/nrm2, x);
	  lis_matvec(A, x, Ax);
	  lis_vector_dot(x, Ax, &xAx);
	  lis_vector_dot(x, x, &xx);
	  mu = xAx / xx;
	}

      iter = 0;
      ptime = 0;

      while (iter<emaxiter)
	{
	  iter = iter+1;
	  giter = giter+1;

	  /* QR factorization VR = Z for starting vector Z */
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

	      /* R = (A - lshift I)^-1 * V */
	      lis_solve_kernel(A, v[j], r, solver, precon);

	      break;

	    case LIS_ESOLVER_AII:

 	      /* R = (M - lshift I)^-1 * V */
	      time = lis_wtime();
	      lis_psolve(solver, v[j], r); 
	      ptime += lis_wtime() - time;

	      break;

	    case LIS_ESOLVER_RQI:

	      /* R = (A - mu I)^-1 * V */
	      lis_vector_nrm2(v[j], &nrm2);
	      lis_vector_scale(1/nrm2, v[j]);
	      lis_matrix_shift_diagonal(A, -mu);
	      lis_solve_kernel(A, v[j], r, solver, precon);
	      lis_matrix_shift_diagonal(A, mu);

	      break;
	    }

	  /* elapsed time of linear solver */
	  if ( j==1 && ( niesolver==LIS_ESOLVER_II || niesolver==LIS_ESOLVER_AII || niesolver==LIS_ESOLVER_RQI ))
	    {
	      lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
	      esolver->ptime += solver->ptime;
	      esolver->itime += solver->itime;
	      esolver->p_c_time += solver->p_c_time;
	      esolver->p_i_time += solver->p_i_time;
	    }

	  /* R = V*R */
	  lis_vector_nrm2(r, &nrm2);
	  lis_vector_dot(v[j],r,&dotvr);
	  mu = mu + 1/dotvr;

	  /* resid = ||Z - VR||_2 */
	  lis_vector_axpyz(-dotvr,v[j],r,q);
	  lis_vector_nrm2(q, &resid);
	  resid = fabs(resid / dotvr);

	  lis_vector_scale(1/nrm2,r);
	  lis_vector_copy(r, v[j]);

	  /* convergence check */
	  if ( j==1 ) 
	    {
	      if( output & LIS_PRINT_MEM ) esolver->rhistory[iter] = resid; 
	      esolver->iter[j-1] = iter;
	      /* esolver->resid = resid; */
	    }


	  if( output & LIS_PRINT_OUT && A->my_rank==0 ) lis_print_rhistory(iter,resid);
	  if (tol>resid) break;
	}

      switch ( niesolver )
	{
	case LIS_ESOLVER_PI:
	  esolver->evalue[j-1] = dotvr;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	case LIS_ESOLVER_II:
	  esolver->evalue[j-1] = 1/dotvr;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	case LIS_ESOLVER_AII:
	  esolver->evalue[j-1] = 1/dotvr;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	case LIS_ESOLVER_RQI:
	  esolver->evalue[j-1] = mu;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	}

      lis_vector_copy(v[j], esolver->evector[j-1]);  

      /* elapsed time of linear solver for Approximate Inverse */
      if ( niesolver==LIS_ESOLVER_AII )
	{
	  lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
	  esolver->ptime = ptime;
	  esolver->itime = solver->itime;
	  esolver->p_c_time = solver->p_c_time;
	  esolver->p_i_time = solver->p_i_time;
	}
      
      if (A->my_rank==0 && ss>1)
	{
#ifdef _LONG__LONG
	  if( output ) printf("Subspace: mode number          = %lld\n", j-1);
#else
	  if( output ) printf("Subspace: mode number          = %d\n", j-1);
#endif
#ifdef _COMPLEX	  
#ifdef _LONG__DOUBLE
	  if( output ) printf("Subspace: eigenvalue           = (Le, %Le)\n", creall(esolver->evalue[j-1]), cimagl(esolver->evalue[j-1]));
#else
	  if( output ) printf("Subspace: eigenvalue           = (%e, %e)\n", creal(esolver->evalue[j-1]), cimag(esolver->evalue[j-1]));
#endif
#else
#ifdef _LONG__DOUBLE
	  if( output ) printf("Subspace: eigenvalue           = %Le\n", esolver->evalue[j-1]);
#else
	  if( output ) printf("Subspace: eigenvalue           = %e\n", esolver->evalue[j-1]);
#endif
#endif	  
#ifdef _LONG__LONG
	  if( output ) printf("Subspace: number of iterations = %lld\n",iter);
#else
	  if( output ) printf("Subspace: number of iterations = %d\n",iter);
#endif
#ifdef _LONG__DOUBLE
	  if( output ) printf("Subspace: relative residual    = %Le\n\n",resid);
#else
	  if( output ) printf("Subspace: relative residual    = %e\n\n",resid);
#endif	  
	}
    }
  
  lis_vector_copy(esolver->evector[0], esolver->x);

  switch ( niesolver )
    {
    case LIS_ESOLVER_II:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_precon_destroy(precon);
      lis_solver_destroy(solver);
      break;
    case LIS_ESOLVER_AII:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_precon_destroy(precon);
      lis_solver_destroy(solver);
      break;
    case LIS_ESOLVER_RQI:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_precon_destroy(precon);
      lis_solver_destroy(solver);
      break;
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_esi_quad"
LIS_INT lis_esi_quad(LIS_ESOLVER esolver)
{
  LIS_MATRIX A;
  LIS_VECTOR x, Ax;
  LIS_SCALAR xAx, xx, mu;
  LIS_REAL lshift;
  LIS_INT ss;
  LIS_INT emaxiter;
  LIS_REAL tol;
  LIS_INT j,k;
  LIS_SCALAR dotvr;
  LIS_INT iter,giter,output,niesolver;
  LIS_REAL nrm2,resid;
  LIS_QUAD_PTR qdot_vv, qdot_vr;
  LIS_VECTOR *v,r,q;
  LIS_SOLVER solver;
  LIS_PRECON precon;
  double time,itime,ptime,p_c_time,p_i_time;
  LIS_INT err;
  LIS_INT nsol, precon_type;
  char solvername[128], preconname[128];

  LIS_DEBUG_FUNC_IN;

  A = esolver->A;
  x = esolver->x;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  lshift = esolver->lshift;
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];

  r = esolver->work[0];
  q = esolver->work[1];
  v = &esolver->work[2];
  Ax = esolver->work[3];

  LIS_QUAD_SCALAR_MALLOC(qdot_vv,0,1);
  LIS_QUAD_SCALAR_MALLOC(qdot_vr,1,1);

  lis_vector_set_all(1.0,r);
  lis_vector_nrm2(r, &nrm2);
  lis_vector_scale(1/nrm2,r);
	  
  switch ( niesolver )
    {
    case LIS_ESOLVER_II:

      lis_solver_create(&solver);
      lis_solver_set_option("-i bicg -p none -precision quad",solver);
      lis_solver_set_optionC(solver);
      lis_solver_get_solver(solver, &nsol);
      lis_solver_get_precon(solver, &precon_type);
      lis_solver_get_solvername(nsol, solvername);
      lis_solver_get_preconname(precon_type, preconname);
      if( A->my_rank==0 ) printf("linear solver         : %s\n", solvername);
      if( A->my_rank==0 ) printf("preconditioner        : %s\n", preconname);
      if( A->my_rank==0 ) printf("\n");
#ifdef _LONG__DOUBLE
      if( A->my_rank==0 ) printf("local shift           : %Le\n", lshift);
#else
      if( A->my_rank==0 ) printf("local shift           : %e\n", lshift);
#endif
      if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
      break;

    case LIS_ESOLVER_AII:

      lis_solver_create(&solver);
      lis_solver_set_option("-i bicg -p none -precision quad",solver);
      lis_solver_set_optionC(solver);
      lis_solver_get_solver(solver, &nsol);
      lis_solver_get_precon(solver, &precon_type);
      lis_solver_get_solvername(nsol, solvername);
      lis_solver_get_preconname(precon_type, preconname);
      if( A->my_rank==0 ) printf("linear solver         : %s\n", solvername);
      if( A->my_rank==0 ) printf("preconditioner        : %s\n", preconname);
      if( A->my_rank==0 ) printf("\n");
#ifdef _LONG__DOUBLE
      if( A->my_rank==0 ) printf("local shift           : %Le\n", lshift);
#else
      if( A->my_rank==0 ) printf("local shift           : %e\n", lshift);
#endif
      if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
      lis_vector_set_all(1.0,q);
      lis_solve(A, q, x, solver);
      lis_precon_create(solver, &precon);
      solver->precon = precon;
      break;

    case LIS_ESOLVER_RQI:

      lis_solver_create(&solver);
      lis_solver_set_option("-p none -precision quad -maxiter 10",solver);
      lis_solver_set_optionC(solver);
      lis_solver_get_solver(solver, &nsol);
      lis_solver_get_precon(solver, &precon_type);
      lis_solver_get_solvername(nsol, solvername);
      lis_solver_get_preconname(precon_type, preconname);
      if( A->my_rank==0 ) printf("linear solver         : %s\n", solvername);
      if( A->my_rank==0 ) printf("preconditioner        : %s\n", preconname);
      if( A->my_rank==0 ) printf("\n");
#ifdef _LONG__DOUBLE
      if( A->my_rank==0 ) printf("local shift           : %Le\n", lshift);
#else
      if( A->my_rank==0 ) printf("local shift           : %e\n", lshift);
#endif
      if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
      break;

    }

  giter=0;
  j=0;
  while (j<ss)
    {
      lis_vector_duplicate(A,&esolver->evector[j]); 
      j = j+1;
      lis_vector_copy(r, v[j]);

      if (niesolver==LIS_ESOLVER_II || niesolver==LIS_ESOLVER_RQI)
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

      if (niesolver==LIS_ESOLVER_RQI)
	{
	  lis_vector_nrm2(x, &nrm2);
	  lis_vector_scale(1/nrm2, x);
	  lis_matvec(A, x, Ax);
	  lis_vector_dot(x, Ax, &xAx);
	  lis_vector_dot(x, x, &xx);
	  mu = xAx / xx;
	}

      iter = 0;
      while (iter<emaxiter)
	{
	  iter = iter+1;
	  giter = giter+1;

	  /* QR factorization VR = Z for starting vector Z */
	  for (k=1;k<j;k++)
	    { 
	      lis_vector_dotex_mmm(v[j], v[k], &qdot_vv);
	      lis_quad_minus((LIS_QUAD *)qdot_vv.hi);
	      lis_vector_axpyex_mmm(qdot_vv,v[k],v[j]);
	    }

	  /* kernel */
	  switch( niesolver )
	    {
	    case LIS_ESOLVER_PI:

	      /* R = A * V */
	      lis_matvec(A,v[j],r); 

	      break;

	    case LIS_ESOLVER_II:

	      /* R = (A - lshift I)^-1 * V */
	      lis_solve_kernel(A, v[j], r, solver, precon);

	      break;

	    case LIS_ESOLVER_AII:

 	      /* R = (M - lshift I)^-1 * V */
	      lis_psolve(solver, v[j], r); 

	      break;

	    case LIS_ESOLVER_RQI:

	      /* R = (A - mu I)^-1 * V */
	      lis_vector_nrm2(v[j], &nrm2);
	      lis_vector_scale(1/nrm2, v[j]);
	      lis_matrix_shift_diagonal(A, -mu);
	      lis_solve_kernel(A, v[j], r, solver, precon);
	      lis_matrix_shift_diagonal(A, mu);

	      break;
	    }

	  /* elapsed time of linear solver */
	  if ( j==1 && ( niesolver==LIS_ESOLVER_II || niesolver==LIS_ESOLVER_AII || niesolver==LIS_ESOLVER_RQI ))
	    {
	      lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
	      esolver->ptime += solver->ptime;
	      esolver->itime += solver->itime;
	      esolver->p_c_time += solver->p_c_time;
	      esolver->p_i_time += solver->p_i_time;
	    }

	  /* R = V*R */
	  lis_vector_nrm2(r, &nrm2);
	  lis_vector_dotex_mmm(v[j], r, &qdot_vr);
	  lis_quad_minus((LIS_QUAD *)qdot_vr.hi);
	  lis_vector_axpyzex_mmmm(qdot_vr,v[j],r,q);
	  lis_quad_minus((LIS_QUAD *)qdot_vr.hi);	  
	  dotvr = qdot_vr.hi[0];
	  mu = mu + 1/dotvr;

	  /* resid = ||Z - VR||_2 */
	  lis_vector_nrm2(q, &resid);
	  resid = fabs(resid / dotvr);
	  lis_vector_scale(1/nrm2,r);
	  lis_vector_copy(r, v[j]);

	  /* convergence check */
	  if ( j==1 ) 
	    {
	      if( output & LIS_PRINT_MEM ) esolver->rhistory[iter] = resid; 
	      esolver->iter[j-1] = iter;
	      /* esolver->resid = resid;*/
	    }
	  if( output & LIS_PRINT_OUT && A->my_rank==0 ) lis_print_rhistory(iter,resid);
	  if (tol>resid) break;
	}

      switch ( niesolver )
	{
	case LIS_ESOLVER_PI:
  	  esolver->evalue[j-1] = dotvr;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	case LIS_ESOLVER_II:
	  esolver->evalue[j-1] = 1/dotvr;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	case LIS_ESOLVER_AII:
	  esolver->evalue[j-1] = 1/dotvr;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	case LIS_ESOLVER_RQI:
	  esolver->evalue[j-1] = mu;
	  esolver->resid[j-1] = resid;
	  esolver->iter[j-1] = iter;
	  break;
	}
      lis_vector_copy(v[j], esolver->evector[j-1]);  

      if (A->my_rank==0 && ss>1)
	{
#ifdef _LONG__LONG
	  if( output ) printf("Subspace: mode number          = %lld\n", j-1);
#else
	  if( output ) printf("Subspace: mode number          = %d\n", j-1);
#endif
#ifdef _COMPLEX	  
#ifdef _LONG__DOUBLE
	  if( output ) printf("Subspace: eigenvalue           = (%Le, %Le)\n", creall(esolver->evalue[j-1]), cimagl(esolver->evalue[j-1]));
#else
	  if( output ) printf("Subspace: eigenvalue           = (%e, %e)\n", creal(esolver->evalue[j-1]), cimag(esolver->evalue[j-1]));
#endif
#else
#ifdef _LONG__DOUBLE
	  if( output ) printf("Subspace: eigenvalue           = %Le\n", esolver->evalue[j-1]);
#else
	  if( output ) printf("Subspace: eigenvalue           = %e\n", esolver->evalue[j-1]);
#endif
#endif	  
#ifdef _LONG__LONG
	  if( output ) printf("Subspace: number of iterations = %lld\n",iter);
#else
	  if( output ) printf("Subspace: number of iterations = %d\n",iter);
#endif
#ifdef _LONG__DOUBLE
	  if( output ) printf("Subspace: relative residual    = %Le\n\n",resid);
#else
	  if( output ) printf("Subspace: relative residual    = %e\n\n",resid);
#endif
	}
    }
  
  lis_vector_copy(esolver->evector[0], esolver->x);

  switch ( niesolver )
    {
    case LIS_ESOLVER_II:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_precon_destroy(precon);
      lis_solver_destroy(solver);
      break;
    case LIS_ESOLVER_AII:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_precon_destroy(precon);
      lis_solver_destroy(solver);
      break;
    case LIS_ESOLVER_RQI:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_precon_destroy(precon);
      lis_solver_destroy(solver);
      break;
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}
#endif

