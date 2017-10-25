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

#undef __FUNC__
#define __FUNC__ "lis_matvec_vbr"
void lis_matvec_vbr(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k;
	LIS_INT bi,bj,bc,bn;
	LIS_INT nr;
	LIS_INT n;
	LIS_SCALAR t;

	n   = A->n;
	nr  = A->nr;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,i,j,k,t,bn)
		#endif
		for(bi=0; bi<nr; bi++)
		{
			bn = A->D->bns[bi];
			k  = A->L->row[bi];
			for(i=0;i<bn;i++)
			{
				t = 0.0;
				for(j=0;j<bn;j++)
				{
					t += A->D->v_value[bi][i*bn+j] * x[k+j];
				}
				y[k+i] = t;
			}
/*			lis_array_matvec(bn,A->D->v_value[i],&x[k],&y[k],LIS_INS_VALUE);*/
		}
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j,k,bi,bj,bc)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			for(bc=A->L->bptr[bi];bc<A->L->bptr[bi+1];bc++)
			{
				bj   = A->L->bindex[bc];
				k    = A->L->ptr[bc];
				for(j=A->L->col[bj];j<A->L->col[bj+1];j++)
				{
					for(i=A->L->row[bi];i<A->L->row[bi+1];i++)
					{
						y[i] += A->L->value[k] * x[j];
						k++;
					}
				}
			}
			for(bc=A->U->bptr[bi];bc<A->U->bptr[bi+1];bc++)
			{
				bj   = A->U->bindex[bc];
				k    = A->U->ptr[bc];
				for(j=A->U->col[bj];j<A->U->col[bj+1];j++)
				{
					for(i=A->U->row[bi];i<A->U->row[bi+1];i++)
					{
						y[i] += A->U->value[k] * x[j];
						k++;
					}
				}
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0; i<n; i++)
		{
			y[i] = 0.0;
		}
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j,k,bi,bj,bc)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			for(bc=A->bptr[bi];bc<A->bptr[bi+1];bc++)
			{
				bj   = A->bindex[bc];
				k    = A->ptr[bc];
				for(j=A->col[bj];j<A->col[bj+1];j++)
				{
					for(i=A->row[bi];i<A->row[bi+1];i++)
					{
						y[i] += A->value[k] * x[j];
						k++;
					}
				}
			}
		}
	}
}

#undef __FUNC__
#define __FUNC__ "lis_matvech_vbr"
void lis_matvech_vbr(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k;
	LIS_INT bi,bj,bc,bn;
	LIS_INT nr;
	LIS_INT n,np;
	#ifdef _OPENMP
		LIS_INT nprocs,my_rank;
		LIS_SCALAR t;
		LIS_SCALAR *w;
	#endif

	n   = A->n;
	np  = A->np;
	nr  = A->nr;

  	if( A->is_splited )
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			w = (LIS_SCALAR *)lis_malloc( nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_vbr::w" );
			#pragma omp parallel private(bi,bc,bj,i,j,k,bn,my_rank,t)
			{
				my_rank = omp_get_thread_num();

				#pragma omp for
				for(j=0;j<nprocs;j++)
				{
					memset( &w[j*np], 0, np*sizeof(LIS_SCALAR) );
				}
				#pragma omp for 
				for(bi=0;bi<nr;bi++)
				{
					bn = A->D->bns[bi];
					k  = A->L->row[bi];
					for(i=0;i<bn;i++)
					{
						t = 0.0;
						for(j=0;j<bn;j++)
						{
							t += conj(A->D->v_value[bi][j*bn+i]) * x[k+j];
						}
						w[my_rank*np + k+i] += t;
					}
					for(bc=A->L->bptr[bi];bc<A->L->bptr[bi+1];bc++)
					{
						bj   = A->L->bindex[bc];
						k    = A->L->ptr[bc];
						for(j=A->L->col[bj];j<A->L->col[bj+1];j++)
						{
							for(i=A->L->row[bi];i<A->L->row[bi+1];i++)
							{
								w[my_rank*np + j] += conj(A->L->value[k]) * x[i];
								k++;
							}
						}
					}
					for(bc=A->U->bptr[bi];bc<A->U->bptr[bi+1];bc++)
					{
						bj   = A->U->bindex[bc];
						k    = A->U->ptr[bc];
						for(j=A->U->col[bj];j<A->U->col[bj+1];j++)
						{
							for(i=A->U->row[bi];i<A->U->row[bi+1];i++)
							{
								w[my_rank*np + j] += conj(A->U->value[k]) * x[i];
								k++;
							}
						}
					}
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
			for(i=0; i<nr; i++)
			{
				bn = A->D->bns[i];
				k  = A->L->row[i];
				lis_array_matvec(bn,A->D->v_value[i],&x[k],&y[k],LIS_INS_VALUE);
			}
			for(bi=0;bi<nr;bi++)
			{
				for(bc=A->L->bptr[bi];bc<A->L->bptr[bi+1];bc++)
				{
					bj   = A->L->bindex[bc];
					k    = A->L->ptr[bc];
					for(j=A->L->col[bj];j<A->L->col[bj+1];j++)
					{
						for(i=A->L->row[bi];i<A->L->row[bi+1];i++)
						{
							y[j] += conj(A->L->value[k]) * x[i];
							k++;
						}
					}
				}
				for(bc=A->U->bptr[bi];bc<A->U->bptr[bi+1];bc++)
				{
					bj   = A->U->bindex[bc];
					k    = A->U->ptr[bc];
					for(j=A->U->col[bj];j<A->U->col[bj+1];j++)
					{
						for(i=A->U->row[bi];i<A->U->row[bi+1];i++)
						{
							y[j] += conj(A->U->value[k]) * x[i];
							k++;
						}
					}
				}
			}
		#endif
	}
	else
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			w = (LIS_SCALAR *)lis_malloc( nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_vbr::w" );
			#pragma omp parallel private(bi,bc,bj,i,j,k,my_rank)
			{
				my_rank = omp_get_thread_num();

				#pragma omp for
				for(j=0;j<nprocs;j++)
				{
					memset( &w[j*np], 0, np*sizeof(LIS_SCALAR) );
				}
				#pragma omp for 
				for(bi=0;bi<nr;bi++)
				{
					for(bc=A->bptr[bi];bc<A->bptr[bi+1];bc++)
					{
						bj   = A->bindex[bc];
						k    = A->ptr[bc];
						for(j=A->col[bj];j<A->col[bj+1];j++)
						{
							for(i=A->row[bi];i<A->row[bi+1];i++)
							{
								w[my_rank*np + j] += conj(A->value[k]) * x[i];
								k++;
							}
						}
					}
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
			for(i=0; i<n; i++)
			{
				y[i] = 0.0;
			}
			for(bi=0;bi<nr;bi++)
			{
				for(bc=A->bptr[bi];bc<A->bptr[bi+1];bc++)
				{
					bj   = A->bindex[bc];
					k    = A->ptr[bc];
					for(j=A->col[bj];j<A->col[bj+1];j++)
					{
						for(i=A->row[bi];i<A->row[bi+1];i++)
						{
							y[j] += conj(A->value[k]) * x[i];
							k++;
						}
					}
				}
			}
		#endif
	}
}
