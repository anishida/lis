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
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

void lis_matvec_coo(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k;
	LIS_INT n,nnz;

	n    = A->n;
	nnz  = A->nnz;
	if( A->is_splited )
	{
		for(i=0; i<n; i++)
		{
			y[i] = A->D->value[i] * x[i];
		}
		for(k=0; k<A->L->nnz; k++)
		{
			i = A->L->row[k];
			j = A->L->col[k];
			y[i] += A->L->value[k] * x[j];
		}
		for(k=0; k<A->U->nnz; k++)
		{
			i = A->U->row[k];
			j = A->U->col[k];
			y[i] += A->U->value[k] * x[j];
		}
	}
	else
	{
		for(i=0; i<n; i++)
		{
			y[i] = 0.0;
		}
		for(k=0; k<nnz; k++)
		{
			i = A->row[k];
			j = A->col[k];
			y[i] += A->value[k] * x[j];
		}
	}
}

void lis_matvech_coo(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k;
	LIS_INT n,nnz;

	n    = A->n;
	nnz  = A->nnz;
	if( A->is_splited )
	{
		for(i=0;i<n;i++)
		{
			y[i] = conj(A->D->value[i]) * x[i];
		}
		for(k=0;k<A->L->nnz;k++)
		{
			i = A->L->row[k];
			j = A->L->col[k];
			y[j] += conj(A->L->value[k]) * x[i];
		}
		for(k=0;k<A->U->nnz;k++)
		{
			i = A->U->row[k];
			j = A->U->col[k];
			y[j] += conj(A->U->value[k]) * x[i];
		}
	}
	else
	{
		for(i=0; i<n; i++)
		{
			y[i] = 0.0;
		}
		for(k=0; k<nnz; k++)
		{
			i = A->row[k];
			j = A->col[k];
			y[j] += A->value[k] * x[i];
		}
	}
}
