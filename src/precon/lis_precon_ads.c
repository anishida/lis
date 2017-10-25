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
#include <string.h>
#include <stdarg.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_precon_create
 * lis_psolve
 * lis_psolveh
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_precon_create_adds"
LIS_INT lis_precon_create_adds(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	i,j;
	LIS_INT	precon_type,worklen;
	LIS_INT	err;
	LIS_VECTOR *work;

	LIS_DEBUG_FUNC_IN;

	precon_type = solver->options[LIS_OPTIONS_PRECON];
	worklen     = 2;
	work        = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_precon_create_adds::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_OUT_OF_MEMORY;
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
	precon->worklen = worklen;
	precon->work    = work;

	err = lis_precon_create_xxx[precon_type](solver,precon);
	if( err )
	{
		lis_precon_destroy(precon);
		return err;
	}

    precon->A       = solver->A;
    precon->is_copy = LIS_FALSE;

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_adds"
LIS_INT lis_psolve_adds(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_INT i,k,n,np,iter,ptype;
	LIS_SCALAR *b,*x,*w,*r,*rl;
	LIS_VECTOR W,R;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	n     = precon->A->n;
	np    = precon->A->np;
	W     = precon->work[0];
	R     = precon->work[1];
	b     = B->value;
	x     = X->value;
	w     = W->value;
	r     = R->value;
	rl    = R->value_lo;
	iter  = solver->options[LIS_OPTIONS_ADDS_ITER];
	ptype = solver->options[LIS_OPTIONS_PRECON];

	#ifdef USE_QUAD_PRECISION
	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
	#endif
		lis_vector_set_all(0.0,X);
		lis_vector_copy(B,R);
		for(k=0;k<iter+1;k++)
		{
			for(i=n;i<np;i++)
			{
				r[i] = 0.0;
			}

			lis_psolve_xxx[ptype](solver,R,W);
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<n;i++)
			{
				x[i] += w[i];
			}
		
			if(k!=iter)
			{
				lis_matvec(precon->A,X,R);
				#ifdef _OPENMP
				#pragma omp parallel for private(i)
				#endif
				for(i=0;i<n;i++)
				{
					r[i] = b[i] - r[i];
				}
			}
		}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_set_allex_nm(0.0,X);
			lis_vector_copyex_mm(B,R);
			for(k=0;k<iter+1;k++)
			{
				for(i=n;i<np;i++)
				{
					r[i] = 0.0;
					rl[i] = 0.0;
				}

				lis_psolve_xxx[ptype](solver,R,W);
				for(i=0;i<n;i++)
				{
					#ifndef USE_SSE2
						LIS_QUAD_ADD(X->value[i],X->value_lo[i],X->value[i],X->value_lo[i],W->value[i],W->value_lo[i]);
					#else
						LIS_QUAD_ADD_SSE2(X->value[i],X->value_lo[i],X->value[i],X->value_lo[i],W->value[i],W->value_lo[i]);
					#endif
	/*				x[i] += w[i];*/
				}
			
				if(k==iter) break;

				lis_matvec(precon->A,X,R);
				for(i=0;i<n;i++)
				{
					#ifndef USE_SSE2
						LIS_QUAD_ADD(R->value[i],R->value_lo[i],B->value[i],B->value_lo[i],-R->value[i],-R->value_lo[i]);
					#else
						LIS_QUAD_ADD_SSE2(R->value[i],R->value_lo[i],B->value[i],B->value_lo[i],-R->value[i],-R->value_lo[i]);
					#endif
	/*				r[i] = b[i] - r[i];*/
				}
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_adds"
LIS_INT lis_psolveh_adds(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_INT i,k,n,np,iter,ptype;
	LIS_SCALAR *b,*x,*w,*r;
	LIS_VECTOR W,R;
	LIS_PRECON precon;


	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	n     = precon->A->n;
	np    = precon->A->np;
	W     = precon->work[0];
	R     = precon->work[1];
	b     = B->value;
	x     = X->value;
	w     = W->value;
	r     = R->value;
	iter  = solver->options[LIS_OPTIONS_ADDS_ITER];
	ptype = solver->options[LIS_OPTIONS_PRECON];

	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
		lis_vector_set_all(0.0,X);
		lis_vector_copy(B,R);
		for(k=0;k<iter+1;k++)
		{
			for(i=n;i<np;i++)
			{
				r[i] = 0.0;
			}

			lis_psolveh_xxx[ptype](solver,R,W);
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<n;i++)
			{
				x[i] += w[i];
			}
		
			if(k!=iter)
			{
				lis_matvech(precon->A,X,R);
				#ifdef _OPENMP
				#pragma omp parallel for private(i)
				#endif
				for(i=0;i<n;i++)
				{
					r[i] = b[i] - r[i];
				}
			}
		}
	}
	else
	{
		lis_vector_set_all(0.0,X);
		lis_vector_copy(B,R);
		for(k=0;k<iter+1;k++)
		{
			for(i=n;i<np;i++)
			{
				r[i] = 0.0;
			}

			lis_psolveh_xxx[ptype](solver,R,W);
			for(i=0;i<n;i++)
			{
				x[i] += w[i];
			}
		
			if(k==iter) break;

			X->precision = LIS_PRECISION_DEFAULT;
			lis_matvech(precon->A,X,R);
			X->precision = LIS_PRECISION_QUAD;
			for(i=0;i<n;i++)
			{
				r[i] = b[i] - r[i];
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

