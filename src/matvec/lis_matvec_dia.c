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

void lis_matvec_dia(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,is,ie,js,je,jj,ii;
	LIS_INT n,np,nnd,k;
	LIS_INT my_rank,nprocs;

	if( A->is_splited )
	{
		n      = A->n;
		np     = A->np;
		#ifdef _OPENMP
		#pragma omp parallel private(i,j,k,is,ie,jj,js,je,ii,my_rank)
		#endif
		{
			#ifdef _OPENMP
				nprocs  = omp_get_max_threads();
				my_rank = omp_get_thread_num();
			#else
				nprocs  = 1;
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie)

			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=is;i<ie;i++)
			{ 
				y[i] = A->D->value[i] * x[i];
			}
			for(j=0;j<A->L->nnd;j++)
			{
				jj = A->L->index[j];
				js = _max(is,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?ie:_min(ie,np-jj);
				#else
					je = _min(ie,n-jj);
				#endif
				k  = is*A->L->nnd + j*(ie-is);
				ii = js-is;
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=js;i<je;i++)
				{ 
					y[i] += A->L->value[k + ii] * x[jj+i];
					ii++;
				}
			}
			for(j=0;j<A->U->nnd;j++)
			{
				jj = A->U->index[j];
				js = _max(is,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?ie:_min(ie,np-jj);
				#else
					je = _min(ie,n-jj);
				#endif
				k  = is*A->U->nnd + j*(ie-is);
				ii = js-is;
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=js;i<je;i++)
				{ 
					y[i] += A->U->value[k + ii] * x[jj+i];
					ii++;
				}
			}
		}
	}
	else
	{
		n      = A->n;
		np     = A->np;
		nnd    = A->nnd;

		#ifdef _OPENMP
		#pragma omp parallel private(i,j,k,is,ie,jj,js,je,ii,my_rank)
		#endif
		{
			#ifdef _OPENMP
				nprocs  = omp_get_max_threads();
				my_rank = omp_get_thread_num();
			#else
				nprocs  = 1;
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie)

			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=is; i<ie; i++)
			{
				y[i] = 0.0;
			}
			for(j=0;j<nnd;j++)
			{
				jj = A->index[j];
				js = _max(is,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?ie:_min(ie,np-jj);
				#else
					je = _min(ie,n-jj);
				#endif
				k  = is*nnd + j*(ie-is);
				ii = js-is;
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=js;i<je;i++)
				{ 
					y[i] += A->value[k + ii] * x[jj+i];
					ii++;
				}
			}
		}
	}
}

void lis_matvech_dia(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,is,ie,js,je,jj,ii;
	LIS_INT n,nnd,np,k;
	LIS_INT my_rank,nprocs;
	#ifdef _OPENMP
		LIS_SCALAR t,*w;
	#endif

	if( A->is_splited )
	{
		n      = A->n;
		#ifdef _OPENMP
		#pragma omp parallel private(i,j,k,is,ie,jj,js,je,ii,my_rank)
		#endif
		{
			#ifdef _OPENMP
				nprocs  = omp_get_max_threads();
				my_rank = omp_get_thread_num();
			#else
				nprocs  = 1;
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie)

			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=is; i<ie; i++)
			{
				y[i] = 0.0;
			}
			for(j=0;j<A->L->nnd;j++)
			{
				jj = A->L->index[j];
				js = _max(is,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?ie:_min(ie,np-jj);
				#else
					je = _min(ie,n-jj);
				#endif
				k  = is*A->L->nnd + j*(ie-is);
				ii = js-is;
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=js;i<je;i++)
				{ 
					y[jj+i] += conj(A->L->value[k + ii]) * x[i];
					ii++;
				}
			}
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=is;i<ie;i++)
			{ 
				y[i] += conj(A->D->value[i]) * x[i];
			}
			for(j=0;j<A->U->nnd;j++)
			{
				jj = A->U->index[j];
				js = _max(is,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?ie:_min(ie,np-jj);
				#else
					je = _min(ie,n-jj);
				#endif
				k  = is*A->U->nnd + j*(ie-is);
				ii = js-is;
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=js;i<je;i++)
				{ 
					y[jj+i] += conj(A->U->value[k + ii]) * x[i];
					ii++;
				}
			}
		}
	}
	else
	{
		n      = A->n;
		np     = A->np;
		nnd    = A->nnd;

		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			w = (LIS_SCALAR *)lis_malloc( nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_dia::w" );
			#pragma omp parallel private(i,j,k,is,ie,jj,js,je,ii,my_rank)
			{
				nprocs  = omp_get_max_threads();
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie)

				#pragma omp for
				for(j=0;j<nprocs;j++)
				{
					memset( &w[j*np], 0, np*sizeof(LIS_SCALAR) );
				}
				for(j=0;j<nnd;j++)
				{
					jj = A->index[j];
					js = _max(is,-jj);
					#ifdef USE_MPI
						je = jj<=(np-n)?ie:_min(ie,np-jj);
					#else
						je = _min(ie,n-jj);
					#endif
					k  = is*nnd + j*(ie-is);
					ii = js-is;
					#ifdef USE_VEC_COMP
					#pragma cdir nodep
					#pragma _NEC ivdep
					#endif
					for(i=js;i<je;i++)
					{ 
						w[my_rank*np + jj+i] += conj(A->value[k + ii]) * x[i];
						ii++;
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
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=0;i<np;i++)
			{
				y[i] = 0.0;
			}
			for(j=0;j<nnd;j++)
			{
				jj = A->index[j];
				js = _max(0,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?n:_min(n,np-jj);
				#else
					je = _min(n,n-jj);
				#endif
				k  = j*n;
				ii = js;
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=js;i<je;i++)
				{ 
					y[jj+i] += conj(A->value[k + ii]) * x[i];
					ii++;
				}
			}
		#endif
	}
}
