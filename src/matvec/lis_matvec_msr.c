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

void lis_matvec_msr(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT n;
	LIS_SCALAR t;

	n      = A->n;
	if( A->is_splited )
	{
		n    = A->n;
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j,js,je,t,jj)
		#endif
		for(i=0; i<n; i++)
		{
			t = A->D->value[i] * x[i];
			js = A->L->index[i];
			je = A->L->index[i+1];
			for(j=js;j<je;j++)
			{
				jj = A->L->index[j];
				t += A->L->value[j] * x[jj];
			}
			js = A->U->index[i];
			je = A->U->index[i+1];
			for(j=js;j<je;j++)
			{
				jj = A->U->index[j];
				t += A->U->value[j] * x[jj];
			}
			y[i] = t;
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j,js,je,t,jj)
		#endif
		for(i=0; i<n; i++)
		{
			t = A->value[i] * x[i];
			js = A->index[i];
			je = A->index[i+1];
			for(j=js;j<je;j++)
			{
				jj = A->index[j];
				t += A->value[j] * x[jj];
			}
			y[i] = t;
		}
	}
}

void lis_matvech_msr(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT n,np;
	LIS_SCALAR t;
	#ifdef _OPENMP
		LIS_INT k,nprocs;
		LIS_SCALAR *w;
	#endif

	n    = A->n;
	np   = A->np;

	if( A->is_splited )
	{
		for(i=0; i<n; i++)
		{
			y[i] = conj(A->D->value[i]) * x[i];
		}
		for(i=0; i<n; i++)
		{
			t = x[i];
			js = A->L->index[i];
			je = A->L->index[i+1];
			for(j=js;j<je;j++)
			{
				jj = A->L->index[j];
				y[jj] += conj(A->L->value[j]) * t;
			}
			js = A->U->index[i];
			je = A->U->index[i+1];
			for(j=js;j<je;j++)
			{
				jj = A->U->index[j];
				y[jj] += conj(A->U->value[j]) * t;
			}
		}
	}
	else
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			w = (LIS_SCALAR *)lis_malloc( nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_msr::w" );
			#pragma omp parallel private(i,j,js,je,t,jj,k)
			{
				k = omp_get_thread_num();
				#pragma omp for
				for(j=0;j<nprocs;j++)
				{
					memset( &w[j*np], 0, np*sizeof(LIS_SCALAR) );
				}
				#pragma omp for 
				for(i=0; i<n; i++)
				{
					js = A->index[i];
					je = A->index[i+1];
					t = x[i];
					for(j=js;j<je;j++)
					{
						jj  = k*np+A->index[j];
						w[jj] += conj(A->value[j]) * t;
					}
					w[k*np+i] += conj(A->value[i]) * x[i];
				}
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
			for(i=0; i<n; i++)
			{
				y[i] = conj(A->value[i]) * x[i];
			}
			for(i=0; i<n; i++)
			{
				t = x[i];
				js = A->index[i];
				je = A->index[i+1];
				for(j=js;j<je;j++)
				{
					jj = A->index[j];
					y[jj] += conj(A->value[j]) * t;
				}
			}
		#endif
	}
}
