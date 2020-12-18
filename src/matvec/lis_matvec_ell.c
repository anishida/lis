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
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

void lis_matvec_ell(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,jj,is,ie;
	LIS_INT n,maxnzr,nprocs,my_rank;

	n      = A->n;
	if( A->is_splited )
	{
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0; i<n; i++)
		{
			y[i] = A->D->value[i]*x[i];
		}
		for(j=0;j<A->L->maxnzr;j++)
		{
			jj = j*n;
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=0;i<n;i++)
			{
				y[i] += A->L->value[jj + i] * x[A->L->index[jj + i]];
			}
		}
		for(j=0;j<A->U->maxnzr;j++)
		{
			jj = j*n;
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=0;i<n;i++)
			{
				y[i] += A->U->value[jj + i] * x[A->U->index[jj + i]];
			}
		}
	}
	else
	{
		maxnzr = A->maxnzr;
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif
		#ifdef _OPENMP
		#pragma omp parallel private(i,j,jj,is,ie,my_rank)
		#endif
		{
			#ifdef _OPENMP
				my_rank = omp_get_thread_num();
			#else
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=is;i<ie;i++)
			{
				y[i] = 0.0;
			}
			for(j=0;j<maxnzr;j++)
			{
				jj = j*n;
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=is;i<ie;i++)
				{
					y[i] += A->value[jj + i] * x[A->index[jj + i]];
				}
			}
		}
	}
}

void lis_matvech_ell(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,jj;
	LIS_INT n,np,maxnzr;
	#ifdef _OPENMP
		LIS_INT k,is,ie,nprocs;
		LIS_SCALAR t;
		LIS_SCALAR *w;
	#endif

	n      = A->n;
	np     = A->np;
	if( A->is_splited )
	{
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0; i<n; i++)
		{
			y[i] = conj(A->D->value[i])*x[i];
		}
		for(j=0;j<A->L->maxnzr;j++)
		{
			jj = j*n;
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=0;i<n;i++)
			{
				y[A->L->index[jj + i]] += conj(A->L->value[jj + i]) * x[i];
			}
		}
		for(j=0;j<A->U->maxnzr;j++)
		{
			jj = j*n;
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=0;i<n;i++)
			{
				y[A->U->index[jj + i]] += conj(A->U->value[jj + i]) * x[i];
			}
		}
	}
	else
	{
		#ifdef _OPENMP
			maxnzr = A->maxnzr;
			nprocs = omp_get_max_threads();
			w = (LIS_SCALAR *)lis_malloc( nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_ell::w" );
			#pragma omp parallel private(i,j,t,jj,k,is,ie)
			{
				k = omp_get_thread_num();
				LIS_GET_ISIE(k,nprocs,n,is,ie);
				#pragma omp for
				for(j=0;j<nprocs;j++)
				{
					memset( &w[j*np], 0, np*sizeof(LIS_SCALAR) );
				}
				for(j=0;j<maxnzr;j++)
				{
					jj = j*n;
					#ifdef USE_VEC_COMP
					#pragma cdir nodep
					#pragma _NEC ivdep
					#endif
					for(i=is;i<ie;i++)
					{
						w[k*np + A->index[jj + i]] += conj(A->value[jj + i]) * x[i];
					}
				}
				#pragma omp barrier
				#pragma omp for 
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=0;i<np;i++)
				{
					t = 0.0;
					for(j=0;j<nprocs;j++)
					{
						t += w[j*np+i];
					}
					y[i] = t;
				}
			}
			lis_free(w);
		#else
			maxnzr = A->maxnzr;
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=0; i<n; i++)
			{
				y[i] = 0.0;
			}
			for(j=0;j<maxnzr;j++)
			{
				jj = j*n;
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=0;i<n;i++)
				{
					y[A->index[jj + i]] += conj(A->value[jj + i]) * x[i];
				}
			}
		#endif
	}
}
