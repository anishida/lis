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
  LIS_REAL nrm2,resid0,cert_resid;
  LIS_SCALAR dot;
  LIS_VECTOR *v,r;
  LIS_SCALAR *t,*tq,*tr,evalue,evalue0;
  LIS_SOLVER solver;
  LIS_ESOLVER esolver2;
  char esolvername[128],solvername[128],preconname[128];
  LIS_INT nsol,precon_type;
  LIS_INT rval,final_status;

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

      /*
       * Reorthogonalize the residual that generates v(j+1).
       * Include the current vector as well to remove round-off
       * left after the alpha(j) recurrence update.
       */
      for (k=1;k<=j;k++)
	{
	  lis_vector_dot(v[k], r, &dot);
	  lis_vector_axpy(-dot, v[k], r);
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
  final_status = LIS_SUCCESS;
  for (i=0;i<ss;i++)
    {
      lis_vector_duplicate(A, &esolver->evector[i]); 
      esolver2->ishift = esolver->evalue[i];
      lis_esolve(A, esolver->evector[i], &evalue, esolver2);
      lis_esolver_work_destroy(esolver2); 
      esolver->evalue[i] = evalue;
      esolver->iter[i] = esolver2->iter[0];

      err = lis_esolver_standard_residual(
              A, esolver->evector[i], esolver->evalue[i],
              r, v[0], &cert_resid);
      if( err )
        {
          esolver->retcode = err;
          lis_esolver_destroy(esolver2);
          lis_free(t);
          lis_free(tq);
          lis_free(tr);
          lis_solver_destroy(solver);
          LIS_DEBUG_FUNC_OUT;
          return err;
        }

      esolver->resid[i] = cert_resid;
      if( cert_resid>tol && final_status==LIS_SUCCESS )
        final_status = LIS_MAXITER;

      if (i==0) 
	{
	  evalue0 = esolver->evalue[0];
	  iter0 = esolver2->iter[0];
	  resid0 = esolver->resid[0];
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
	  lis_printf(comm,"Lanczos: relative residual    = %e\n\n",(double)esolver->resid[i]);
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
  esolver->retcode = final_status;
  return final_status;
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

	/*
	 * The untransformed generalized Lanczos recurrence keeps both
	 * w_j = B*v_j and v_j for full reorthogonalization.  They must
	 * not alias.  The transformed path uses one Lanczos basis.
	 */
	if( esolver->transform==LIS_EIGEN_TRANSFORM_AUTO )
	{
		worklen = NWORK + 2*ss;
	}
	else
	{
		worklen = NWORK + ss;
	}
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
#define __FUNC__ "lis_egli_projected_jacobi"
static LIS_INT lis_egli_projected_jacobi(
	LIS_INT n,
	LIS_SCALAR *a,
	LIS_SCALAR *eval,
	LIS_SCALAR *evec,
	LIS_INT *niter,
	LIS_REAL *error)
{
#ifdef _COMPLEX
	LIS_SETERR(
		LIS_ERR_NOT_IMPLEMENTED,
		"transformed generalized Lanczos projected "
		"eigenvectors are not yet implemented for complex arithmetic\n");
	return LIS_ERR_NOT_IMPLEMENTED;
#else
	LIS_INT i,j,p,q,sweep;
	LIS_INT max_sweeps;
	LIS_REAL app,aqq,apq;
	LIS_REAL akp,akq;
	LIS_REAL vip,viq;
	LIS_REAL tau,t,c,s;
	LIS_REAL off,maxoff,scale,threshold;

	LIS_DEBUG_FUNC_IN;

	if( n<=0 || a==NULL || eval==NULL || evec==NULL )
	{
		LIS_SETERR(
			LIS_ERR_ILL_ARG,
			"invalid projected eigensystem arguments\n");
		return LIS_ERR_ILL_ARG;
	}

	for(i=0;i<n*n;i++)
	{
		evec[i] = 0.0;
	}

	for(i=0;i<n;i++)
	{
		evec[i*n+i] = 1.0;
	}

	scale = 0.0;
	for(i=0;i<n;i++)
	{
		if( fabs(a[i*n+i]) > scale )
		{
			scale = fabs(a[i*n+i]);
		}
	}

	if( scale < 1.0 )
	{
		scale = 1.0;
	}

	threshold = 1.0e-14*scale;
	max_sweeps = 50;

	maxoff = 0.0;

	for(sweep=0;sweep<max_sweeps;sweep++)
	{
		maxoff = 0.0;

		for(p=0;p<n-1;p++)
		{
			for(q=p+1;q<n;q++)
			{
				apq = a[p*n+q];
				off = fabs(apq);

				if( off > maxoff )
				{
					maxoff = off;
				}

				if( off <= threshold )
				{
					continue;
				}

				app = a[p*n+p];
				aqq = a[q*n+q];

				tau = (aqq-app)/(2.0*apq);

				if( tau>=0.0 )
				{
					t = 1.0/
						(tau+sqrt(1.0+tau*tau));
				}
				else
				{
					t = -1.0/
						(-tau+sqrt(1.0+tau*tau));
				}

				c = 1.0/sqrt(1.0+t*t);
				s = t*c;

				for(j=0;j<n;j++)
				{
					if( j==p || j==q )
					{
						continue;
					}

					akp = a[j*n+p];
					akq = a[j*n+q];

					a[j*n+p] = c*akp-s*akq;
					a[p*n+j] = a[j*n+p];

					a[j*n+q] = s*akp+c*akq;
					a[q*n+j] = a[j*n+q];
				}

				a[p*n+p] = app-t*apq;
				a[q*n+q] = aqq+t*apq;
				a[p*n+q] = 0.0;
				a[q*n+p] = 0.0;

				for(j=0;j<n;j++)
				{
					vip = evec[j*n+p];
					viq = evec[j*n+q];

					evec[j*n+p] = c*vip-s*viq;
					evec[j*n+q] = s*vip+c*viq;
				}
			}
		}

		if( maxoff <= threshold )
		{
			break;
		}
	}

	for(i=0;i<n;i++)
	{
		eval[i] = a[i*n+i];
	}

	if( niter )
	{
		*niter = sweep+1;
	}

	if( error )
	{
		*error = maxoff;
	}

	if( maxoff > threshold )
	{
		LIS_SETERR(
			LIS_MAXITER,
			"projected Jacobi eigensolver did not converge\n");
		return LIS_MAXITER;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}


static LIS_INT
lis_egli_solve_with_retry(
        LIS_MATRIX S,
        LIS_VECTOR b,
        LIS_VECTOR x,
        LIS_SOLVER solver,
        LIS_PRECON *precon,
        LIS_REAL acceptance_tol)
{
        LIS_INT err,status;

        status =
                lis_esolver_inner_solve_attempt(
                        S,b,x,solver,*precon);
        if( status==LIS_SUCCESS )
        {
                return LIS_SUCCESS;
        }

        if( !lis_esolver_inner_solve_retryable(status) )
        {
                return status;
        }

        err = lis_esolver_rebuild_precon(
                S,solver,precon);
        if( err )
        {
                return err;
        }

        lis_vector_set_all(0.0,x);

        status =
                lis_esolver_inner_solve_attempt(
                        S,b,x,solver,*precon);

        return lis_esolver_inner_solve_certify(
                S,b,x,solver,
                status,
                acceptance_tol);
}


static LIS_INT
lis_egli_replace_shifted_operator(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_SCALAR shift,
        LIS_SOLVER solver,
        LIS_PRECON *precon,
        LIS_MATRIX *S,
        LIS_MATRIX *Sop)
{
        LIS_INT err;
        LIS_MATRIX trial_operator;

        trial_operator = NULL;

        err = lis_matrix_create_operator(
                1.0,A,-shift,B,
                &trial_operator);
        if( err )
        {
                return err;
        }

        err = lis_esolver_rebuild_precon(
                trial_operator,
                solver,
                precon);
        if( err )
        {
                /*
                 * The rebuild helper has detached work/preconditioner
                 * state tied to the previous operator.  Keep solver->A
                 * valid for cleanup before discarding the failed trial.
                 */
                solver->A = *S;
                lis_matrix_destroy(trial_operator);
                return err;
        }

        if( *Sop )
        {
                lis_matrix_destroy(*Sop);
        }

        *Sop = trial_operator;
        *S = trial_operator;

        return LIS_SUCCESS;
}


static LIS_INT
lis_egli_first_solve_with_recovery(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_VECTOR b,
        LIS_VECTOR x,
        LIS_SOLVER solver,
        LIS_PRECON *precon,
        LIS_MATRIX *S,
        LIS_MATRIX *Sop,
        LIS_SCALAR requested_shift,
        LIS_REAL acceptance_tol,
        LIS_SCALAR *active_shift)
{
        LIS_INT err,status,attempt;
        LIS_SCALAR perturbation;
        LIS_SCALAR trial_shift;

        *active_shift = requested_shift;

        status =
                lis_esolver_inner_solve_attempt(
                        *S,b,x,solver,*precon);
        if( status==LIS_SUCCESS )
        {
                return LIS_SUCCESS;
        }

        if( !lis_esolver_inner_solve_retryable(status) )
        {
                return status;
        }

        err = lis_esolver_rebuild_precon(
                *S,solver,precon);
        if( err )
        {
                return err;
        }

        lis_vector_set_all(0.0,x);

        status =
                lis_esolver_inner_solve_attempt(
                        *S,b,x,solver,*precon);

                status =
                        lis_esolver_inner_solve_certify(
                                *S,b,x,solver,
                                status,
                                acceptance_tol);
        if( status==LIS_SUCCESS )
        {
                return LIS_SUCCESS;
        }

        if( !lis_esolver_inner_solve_retryable(status) )
        {
                return status;
        }

        perturbation =
                lis_esolver_shift_perturbation(
                        requested_shift);
        if( perturbation==0.0 )
        {
                return status;
        }

        for(attempt=0;attempt<2;attempt++)
        {
                trial_shift =
                        attempt==0 ?
                        requested_shift+perturbation :
                        requested_shift-perturbation;

                err = lis_egli_replace_shifted_operator(
                        A,B,trial_shift,
                        solver,precon,
                        S,Sop);
                if( err )
                {
                        return err;
                }

                lis_vector_set_all(0.0,x);

                status =
                        lis_esolver_inner_solve_attempt(
                                *S,b,x,solver,*precon);

                status =
                        lis_esolver_inner_solve_certify(
                                *S,b,x,solver,
                                status,
                                acceptance_tol);

                if( status==LIS_SUCCESS )
                {
                        *active_shift = trial_shift;
                        return LIS_SUCCESS;
                }

                if( !lis_esolver_inner_solve_retryable(status) )
                {
                        return status;
                }
        }

        return status;
}


#undef __FUNC__
#define __FUNC__ "lis_egli_transformed"
static LIS_INT lis_egli_transformed(LIS_ESOLVER esolver)
{
	LIS_Comm comm;
	LIS_MATRIX A,B,M,S,Sop;
	LIS_VECTOR q,r,z,*v;
	LIS_SOLVER solver;
	LIS_PRECON precon;
	LIS_ESOLVER esolver2;
	LIS_SCALAR *t,*tq,*tr;
	LIS_SCALAR oshift,ishift,active_shift;
	LIS_SCALAR alpha,dot,tmp,nu;
	LIS_SCALAR evalue,evalue0;
	LIS_REAL beta,beta_prev,mag,maxmag;
	LIS_REAL tol,projerr;
	LIS_REAL resid,resid0,nrm2;
	LIS_INT transform,ss,nev,output,rval;
	LIS_INT nigesolver,emaxiter;
	LIS_INT i,j,k,imax,ic;
	LIS_INT projiter,iter0;
	LIS_INT err,inner_err;
	LIS_INT nsol,precon_type;
	LIS_INT final_status;
	char solvername[128],preconname[128],esolvername[128];
	double time,time0;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A = esolver->A;
	B = esolver->B;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
        nev = esolver->nev>0 ? esolver->nev : ss;
	emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
	tol = esolver->params[
		LIS_EPARAMS_RESID-LIS_EOPTIONS_LEN];

	output = esolver->options[LIS_EOPTIONS_OUTPUT];
	rval = esolver->options[LIS_EOPTIONS_RVAL];

	nigesolver =
		esolver->options[
			LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER];

	transform = esolver->transform;

#ifdef _COMPLEX
	oshift =
		esolver->params[
			LIS_EPARAMS_SHIFT-LIS_EOPTIONS_LEN] +
		esolver->params[
			LIS_EPARAMS_SHIFT_IM-LIS_EOPTIONS_LEN]
			*_Complex_I;
#else
	oshift =
		esolver->params[
			LIS_EPARAMS_SHIFT-LIS_EOPTIONS_LEN];
#endif

	ishift = esolver->ishift;

	if( ishift!=0.0 )
	{
		oshift = ishift;
	}

	active_shift = oshift;

	if( ss<1 )
	{
		LIS_SETERR1(
			LIS_ERR_ILL_ARG,
			"Parameter LIS_EOPTIONS_SUBSPACE is %D "
			"(Set greater than 0 for transformed Lanczos)\n",
			ss);
		return LIS_ERR_ILL_ARG;
	}

	if( nev<1 || nev>ss )
	{
	        LIS_SETERR2(
	                LIS_ERR_ILL_ARG,
	                "Parameter -nev is %D "
	                "(Set between 1 and subspace size %D)\n",
	                nev,
	                ss);
	        return LIS_ERR_ILL_ARG;
	}

	if( transform==LIS_EIGEN_TRANSFORM_BUCKLING )
	{
		M = A;
	}
	else
	{
		M = B;
	}

	S = A;
	Sop = NULL;
	solver = NULL;
	precon = NULL;
	esolver2 = NULL;

	t = NULL;
	tq = NULL;
	tr = NULL;

	t = (LIS_SCALAR *)lis_malloc(
		ss*ss*sizeof(LIS_SCALAR),
		"lis_egli_transformed::t");

	tq = (LIS_SCALAR *)lis_malloc(
		ss*ss*sizeof(LIS_SCALAR),
		"lis_egli_transformed::tq");

	tr = (LIS_SCALAR *)lis_malloc(
		ss*ss*sizeof(LIS_SCALAR),
		"lis_egli_transformed::tr");

	if( t==NULL || tq==NULL || tr==NULL )
	{
		LIS_SETERR_MEM(ss*ss*sizeof(LIS_SCALAR));
		err = LIS_ERR_OUT_OF_MEMORY;
		goto cleanup;
	}

	for(i=0;i<ss*ss;i++)
	{
		t[i] = 0.0;
		tq[i] = 0.0;
		tr[i] = 0.0;
	}

	if( oshift!=0.0 )
	{
		err = lis_matrix_create_operator(
			1.0,A,-oshift,B,&Sop);

		if( err ) goto cleanup;

		S = Sop;
	}

	err = lis_solver_create(&solver);
	if( err ) goto cleanup;

	lis_solver_set_option("-i bicg -p none",solver);

	err = lis_solver_set_optionC(solver);
	if( err ) goto cleanup;

	lis_solver_get_solver(solver,&nsol);
	lis_solver_get_precon(solver,&precon_type);

	lis_solver_get_solvername(
		nsol,solvername);

	lis_solver_get_preconname(
		precon_type,preconname);

	lis_esolver_get_esolvername(
		nigesolver,esolvername);

	if( output )
	{
#ifdef _COMPLEX
		lis_printf(
			comm,
			"shift                 : (%e, %e)\n",
			(double)creal(oshift),
			(double)cimag(oshift));
#else
		lis_printf(
			comm,
			"shift                 : %e\n",
			(double)oshift);
#endif

		if( transform==
		    LIS_EIGEN_TRANSFORM_SHIFT_INVERT )
		{
			lis_printf(
				comm,
				"spectral transform    : shift-invert\n");
		}
		else
		{
			lis_printf(
				comm,
				"spectral transform    : buckling\n");
		}

		lis_printf(
			comm,
			"linear solver         : %s\n",
			solvername);

		lis_printf(
			comm,
			"preconditioner        : %s\n",
			preconname);

		lis_printf(
			comm,
			"size of subspace      : %D\n",
			ss);
	}

	solver->A = S;

	err = lis_precon_create(
		solver,&precon);

	if( err ) goto cleanup;

	q = esolver->work[0];
	r = esolver->work[1];
	z = esolver->work[2];
	v = &esolver->work[3];

	lis_vector_set_all(0.0,v[0]);
	lis_vector_set_all(1.0,q);

	/*
	 * Normalize the starting vector in the M inner product.
	 */
	err = lis_matvec(M,q,r);
	if( err ) goto cleanup;

	err = lis_vector_dot(q,r,&dot);
	if( err ) goto cleanup;

#ifdef _COMPLEX
	beta = sqrt(cabs(dot));
#else
	beta = sqrt(fabs(dot));
#endif

	if( !isfinite((double)beta) || beta<=0.0 )
	{
		LIS_SETERR(
			LIS_BREAKDOWN,
			"invalid transformed Lanczos starting norm\n");
		err = LIS_BREAKDOWN;
		goto cleanup;
	}

	lis_vector_copy(q,v[1]);
	lis_vector_scale(1.0/beta,v[1]);

	beta_prev = 0.0;

	time0 = lis_wtime();

	for(j=1;j<=ss;j++)
	{
		/*
		 * r = M v_j
		 */
		err = lis_matvec(M,v[j],r);
		if( err ) goto cleanup;

		/*
		 * q = (A-sigma B)^-1 M v_j
		 *
		 * Only the first transformed solve may adopt a bounded
		 * perturbed shift.  Once the Lanczos basis is started, all
		 * remaining solves use the same transformed operator.
		 */
		if( j==1 )
		{
			err = lis_egli_first_solve_with_recovery(
				A,B,r,q,solver,&precon,
				&S,&Sop,oshift,tol,&active_shift);
		}
		else
		{
			err = lis_egli_solve_with_retry(
				S,r,q,solver,&precon,tol);
		}

		if( err ) goto cleanup;

		esolver->ptime += solver->ptime;
		esolver->itime += solver->itime;
		esolver->p_c_time += solver->p_c_time;
		esolver->p_i_time += solver->p_i_time;

		/*
		 * alpha_j = <v_j,Tv_j>_M
		 */
		err = lis_vector_dot(r,q,&alpha);
		if( err ) goto cleanup;

		t[(j-1)*ss+(j-1)] = alpha;

		if( j>1 )
		{
			lis_vector_axpy(
				-beta_prev,
				v[j-1],
				q);
		}

		lis_vector_axpy(
			-alpha,
			v[j],
			q);

		/*
		 * One full M-orthogonalization pass.
		 */
		err = lis_matvec(M,q,r);
		if( err ) goto cleanup;

		for(k=1;k<=j;k++)
		{
			err = lis_vector_dot(
				v[k],r,&dot);

			if( err ) goto cleanup;

			lis_vector_axpy(
				-dot,
				v[k],
				q);
		}

		if( j==ss )
		{
			break;
		}

		err = lis_matvec(M,q,r);
		if( err ) goto cleanup;

		err = lis_vector_dot(q,r,&dot);
		if( err ) goto cleanup;

#ifdef _COMPLEX
		beta = sqrt(cabs(dot));
#else
		beta = sqrt(fabs(dot));
#endif

		if( !isfinite((double)beta) ||
		    beta<=tol )
		{
			LIS_SETERR(
				LIS_BREAKDOWN,
				"transformed Lanczos invariant subspace "
				"closed before requested subspace size\n");
			err = LIS_BREAKDOWN;
			goto cleanup;
		}

		t[(j-1)*ss+j] = (LIS_SCALAR)beta;
		t[j*ss+j-1] = (LIS_SCALAR)beta;

		lis_vector_copy(
			q,
			v[j+1]);

		lis_vector_scale(
			1.0/beta,
			v[j+1]);

		beta_prev = beta;
	}

	time = lis_wtime()-time0;

	/*
	 * Compute the complete eigensystem of the projected
	 * symmetric tridiagonal matrix.
	 */
	for(i=0;i<ss*ss;i++)
	{
		tq[i] = t[i];
	}

	err = lis_egli_projected_jacobi(
		ss,
		tq,
		esolver->evalue,
		tr,
		&projiter,
		&projerr);

	if( err ) goto cleanup;

	/*
	 * Sort transformed eigenpairs by decreasing |nu|.
	 */
	for(i=0;i<ss-1;i++)
	{
		imax = i;

#ifdef _COMPLEX
		maxmag = cabs(esolver->evalue[i]);
#else
		maxmag = fabs(esolver->evalue[i]);
#endif

		for(j=i+1;j<ss;j++)
		{
#ifdef _COMPLEX
			mag = cabs(esolver->evalue[j]);
#else
			mag = fabs(esolver->evalue[j]);
#endif

			if( mag>maxmag )
			{
				maxmag = mag;
				imax = j;
			}
		}

		if( imax!=i )
		{
			tmp = esolver->evalue[i];
			esolver->evalue[i] =
				esolver->evalue[imax];
			esolver->evalue[imax] = tmp;

			for(k=0;k<ss;k++)
			{
				tmp = tr[k*ss+i];
				tr[k*ss+i] =
					tr[k*ss+imax];
				tr[k*ss+imax] = tmp;
			}
		}
	}

	/*
	 * Map nu to the original generalized eigenvalue lambda.
	 */
	for(i=0;i<ss;i++)
	{
		nu = esolver->evalue[i];

		if( transform==
		    LIS_EIGEN_TRANSFORM_SHIFT_INVERT )
		{
#ifdef _COMPLEX
			mag = cabs(nu);
#else
			mag = fabs(nu);
#endif

			if( mag<=1.0e-30 )
			{
				LIS_SETERR(
					LIS_BREAKDOWN,
					"zero shift-invert Ritz value\n");
				err = LIS_BREAKDOWN;
				goto cleanup;
			}

			esolver->evalue[i] =
				active_shift+1.0/nu;
		}
		else
		{
#ifdef _COMPLEX
			mag = cabs(nu-1.0);
#else
			mag = fabs(nu-1.0);
#endif

			if( mag<=1.0e-30 )
			{
				LIS_SETERR(
					LIS_BREAKDOWN,
					"buckling Ritz value is too close to one\n");
				err = LIS_BREAKDOWN;
				goto cleanup;
			}

			esolver->evalue[i] =
				active_shift*nu/(nu-1.0);
		}

		esolver->iter[i] = ss;
		esolver->resid[i] = 0.0;
	}

	if( output )
	{
		lis_printf(
			comm,
			"projected eigensolver iterations : %D\n",
			projiter);

		lis_printf(
			comm,
			"projected eigensolver error      : %e\n\n",
			(double)projerr);

		lis_printf(
			comm,
			"transformed Lanczos Ritz values:\n\n");

		for(i=0;i<ss;i++)
		{
#ifdef _COMPLEX
			lis_printf(
				comm,
				"Generalized Lanczos: mode number = %D, "
				"eigenvalue = (%e, %e)\n",
				i,
				(double)creal(esolver->evalue[i]),
				(double)cimag(esolver->evalue[i]));
#else
			lis_printf(
				comm,
				"Generalized Lanczos: mode number = %D, "
				"eigenvalue = %e\n",
				i,
				(double)esolver->evalue[i]);
#endif
		}

		lis_printf(
			comm,
			"Generalized Lanczos: elapsed time = %e sec.\n\n",
			time);
	}

	/*
	 * -rval keeps the established GLI behaviour: return the
	 * Ritz values without constructing eigenvectors.
	 */
	if( rval )
	{
		esolver->retcode = LIS_SUCCESS;
		err = LIS_SUCCESS;
		goto cleanup;
	}

	/*
	 * Build physical Ritz vectors:
	 *
	 *     x_i = V y_i
	 *
	 * where columns of tr are projected eigenvectors.
	 */
	final_status = LIS_SUCCESS;
	evalue0 = 0.0;
	iter0 = ss;
	resid0 = 0.0;

	err = lis_esolver_create(&esolver2);
	if( err ) goto cleanup;

	esolver2->options[LIS_EOPTIONS_ESOLVER] =
		nigesolver;

	esolver2->options[LIS_EOPTIONS_SUBSPACE] =
		1;

	esolver2->options[LIS_EOPTIONS_MAXITER] =
		emaxiter;

	esolver2->options[LIS_EOPTIONS_OUTPUT] =
		esolver->options[LIS_EOPTIONS_OUTPUT];

	esolver2->options[LIS_EOPTIONS_INITGUESS_ONES] =
		0;

	esolver2->params[
		LIS_EPARAMS_RESID-LIS_EOPTIONS_LEN] =
		tol;

	esolver2->transform =
		LIS_EIGEN_TRANSFORM_AUTO;

	if( output )
	{
		lis_printf(
			comm,
			"inner eigensolver     : %s\n\n",
			esolvername);
	}

	for(i=0;i<nev;i++)
	{
		err = lis_vector_duplicate(
			A,
			&esolver->evector[i]);

		if( err ) goto cleanup;

		lis_vector_set_all(
			0.0,
			esolver->evector[i]);

		for(k=0;k<ss;k++)
		{
			lis_vector_axpy(
				tr[k*ss+i],
				v[k+1],
				esolver->evector[i]);
		}

		err = lis_vector_nrm2(
			esolver->evector[i],
			&nrm2);

		if( err ) goto cleanup;

		if( !isfinite((double)nrm2) ||
		    nrm2<=0.0 )
		{
			LIS_SETERR(
				LIS_BREAKDOWN,
				"invalid transformed generalized Lanczos Ritz vector norm\n");
			err = LIS_BREAKDOWN;
			goto cleanup;
		}

		lis_vector_scale(
			1.0/nrm2,
			esolver->evector[i]);

		err = lis_esolver_generalized_residual(
			A,
			B,
			esolver->evector[i],
			esolver->evalue[i],
			q,
			r,
			z,
			&resid);

		if( err ) goto cleanup;


		/*
		 * Refine only when the physical generalized residual
		 * is not already converged.
		 */
		if( resid>tol )
		{
			esolver2->ishift =
				esolver->evalue[i];

			/*
			 * Keep the inner inverse-iteration system away from the
			 * nearly singular Ritz pole.  Use one percent of the
			 * nearest distinct Ritz-value spacing and move toward
			 * the original spectral target.
			 */
#ifndef _COMPLEX
			{
			        LIS_REAL gap;
			        LIS_REAL d;
			        LIS_REAL delta;
			        LIS_INT jj;

			        gap = 0.0;

			        for(jj=0;jj<ss;jj++)
			        {
			                if( jj==i )
			                        continue;

			                d = fabs(
			                        (double)(
			                                esolver->evalue[jj] -
			                                esolver->evalue[i]));

			                if( d>0.0 &&
			                    (gap==0.0 || d<gap) )
			                {
			                        gap = d;
			                }
			        }

			        if( gap>0.0 )
			        {
			                delta = 0.01*gap;

			                if( oshift<esolver->evalue[i] )
			                        delta = -delta;

			                esolver2->ishift +=
			                        (LIS_SCALAR)delta;
			        }
			}
#endif


			inner_err = lis_gesolve(
				A,
				B,
				esolver->evector[i],
				&evalue,
				esolver2);

			if( inner_err!=LIS_SUCCESS &&
			    inner_err!=LIS_MAXITER )
			{
				err = inner_err;
				goto cleanup;
			}

			esolver->evalue[i] =
				evalue;

			esolver->iter[i] =
				esolver2->iter[0];

			err = lis_esolver_generalized_residual(
				A,
				B,
				esolver->evector[i],
				esolver->evalue[i],
				q,
				r,
				z,
				&resid);

			if( err ) goto cleanup;

			if( inner_err==LIS_MAXITER &&
			    resid>tol )
			{
				final_status = LIS_MAXITER;
			}

			if( i==0 )
			{
				esolver->ptime +=
					esolver2->ptime;

				esolver->itime +=
					esolver2->itime;

				esolver->p_c_time +=
					esolver2->p_c_time;

				esolver->p_i_time +=
					esolver2->p_i_time;

				if( output & LIS_EPRINT_MEM )
				{
					for(ic=0;
					    ic<esolver2->iter[0]+1;
					    ic++)
					{
						esolver->rhistory[ic] =
							esolver2->rhistory[ic];
					}
				}
			}

			lis_esolver_work_destroy(
				esolver2);
		}

		esolver->resid[i] = resid;

		if( resid>tol )
		{
			final_status = LIS_MAXITER;
		}

		if( output )
		{
			lis_printf(
				comm,
				"Generalized Lanczos: mode number          = %D\n",
				i);

#ifdef _COMPLEX
			lis_printf(
				comm,
				"Generalized Lanczos: eigenvalue           = (%e, %e)\n",
				(double)creal(esolver->evalue[i]),
				(double)cimag(esolver->evalue[i]));
#else
			lis_printf(
				comm,
				"Generalized Lanczos: eigenvalue           = %e\n",
				(double)esolver->evalue[i]);
#endif

			lis_printf(
				comm,
				"Generalized Lanczos: true residual        = %e\n\n",
				(double)esolver->resid[i]);
		}

		if( i==0 )
		{
			evalue0 = esolver->evalue[0];
			iter0 = esolver->iter[0];
			resid0 = esolver->resid[0];
		}
	}

	esolver->evalue[0] = evalue0;
	esolver->iter[0] = iter0;
	esolver->resid[0] = resid0;

	lis_vector_copy(
		esolver->evector[0],
		esolver->x);

	esolver->time = time;
	esolver->retcode = final_status;

	err = final_status;

cleanup:

	if( err!=LIS_SUCCESS )
	{
		esolver->retcode = err;
	}

	if( esolver2 )
	{
		lis_esolver_destroy(esolver2);
	}

	if( precon )
	{
		lis_precon_destroy(precon);
	}

	if( solver )
	{
		lis_solver_destroy(solver);
	}

	if( Sop )
	{
		lis_matrix_destroy(Sop);
	}

	if( t ) lis_free(t);
	if( tq ) lis_free(tq);
	if( tr ) lis_free(tr);

	LIS_DEBUG_FUNC_OUT;
	return err;
}


#undef __FUNC__
#define __FUNC__ "lis_egli"
LIS_INT lis_egli(LIS_ESOLVER esolver)
{
  LIS_Comm comm;

  if( esolver->transform != LIS_EIGEN_TRANSFORM_AUTO )
    {
      return lis_egli_transformed(esolver);
    }
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
  v = &esolver->work[2+ss];
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

      /*
       * Reorthogonalize the dual residual.  With w(k)=B*v(k),
       * subtracting <v(k),r>*w(k) makes q=B^-1*r orthogonal
       * to the existing basis in the B inner product.
       */
      for (k=1;k<=j;k++)
	{
	  lis_vector_dot(v[k], r, &dot);
	  lis_vector_axpy(-dot, w[k], r);
	}

      /* solve B * q = r */
      err = lis_solve_kernel(B, r, q, solver, precon);
      if( err )
	{
	  lis_solver_work_destroy(solver);	  
	  solver->retcode = err;
	  lis_precon_destroy(precon);
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



