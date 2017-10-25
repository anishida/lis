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
#include <string.h>
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

void lis_matvec_dns(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,is,ie;
	LIS_INT np,n,nprocs,my_rank;

	n = A->n;
	np = A->np;
	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
	#else
		nprocs  = 1;
	#endif
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,is,ie,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		for(i=is;i<ie;i++)
		{
			y[i] = 0;
		}
		for(j=0;j<np;j++)
		{
			for(i=is;i<ie;i++)
			{
				y[i] += A->value[j*n+i] * x[j];
			}
		}
	}
}

void lis_matvech_dns(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j;
	LIS_INT np,n;
	LIS_SCALAR t;
	#ifdef _OPENMP
		LIS_INT is,ie,nprocs,my_rank;
		LIS_SCALAR *w;
	#endif

	n  = A->n;
	np = A->np;
	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		w = (LIS_SCALAR *)lis_malloc( nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_dns::w" );
		#pragma omp parallel private(i,j,t,is,ie,my_rank)
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			memset( &w[my_rank*np], 0, np*sizeof(LIS_SCALAR) );
			for(j=0;j<np;j++)
			{
				t = 0.0;
				for(i=is;i<ie;i++)
				{
					t += conj(A->value[j*n+i]) * x[i];
				}
				w[my_rank*np + j] = t;
			}
			#pragma omp barrier
			#pragma omp for 
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
		for(j=0;j<np;j++)
		{
			t = 0.0;
			for(i=0;i<n;i++)
			{
				t += conj(A->value[j*n+i]) * x[i];
			}
			y[j] = t;
		}
	#endif
}
