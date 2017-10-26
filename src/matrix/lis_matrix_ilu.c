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
#define __FUNC__ "lis_matrix_ilu_create"
LIS_INT lis_matrix_ilu_create(LIS_INT n, LIS_INT bs, LIS_MATRIX_ILU *A)
{
	LIS_INT i;
	LIS_INT *nnz;
	LIS_INT	**index;

	LIS_DEBUG_FUNC_IN;

	*A      = NULL;
	nnz     = NULL;
	index   = NULL;

	*A = (LIS_MATRIX_ILU)lis_malloc( sizeof(struct LIS_MATRIX_ILU_STRUCT),"lis_matrix_ilu_create::A" );
	if( NULL==*A )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_ILU_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	memset(*A,0,sizeof(struct LIS_MATRIX_ILU_STRUCT));

	nnz = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_ilu_create::nnz" );
	if( nnz==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	index = (LIS_INT **)lis_malloc( n*sizeof(LIS_INT *),"lis_matrix_ilu_create::index" );
	if( index==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT *));
		return LIS_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0;i<n;i++)
	{
		nnz[i] = 0;
		index[i] = NULL;
	}

	(*A)->n      = n;
	(*A)->bs     = bs;
	(*A)->nnz    = nnz;
	(*A)->index  = index;
	(*A)->nnz_ma = NULL;
	(*A)->value  = NULL;
	(*A)->values = NULL;
	(*A)->bsz    = NULL;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_ilu_setCR"
LIS_INT lis_matrix_ilu_setCR(LIS_MATRIX_ILU A)
{
	LIS_INT n;
	LIS_SCALAR **value;

	LIS_DEBUG_FUNC_IN;

	n = A->n;
	value = (LIS_SCALAR **)lis_malloc( n*sizeof(LIS_SCALAR *),"lis_matrix_ilu_setCR::value" );
	if( value==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	A->value = value;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_ilu_setVR"
LIS_INT lis_matrix_ilu_setVR(LIS_MATRIX_ILU A)
{
	LIS_INT n;
	LIS_INT *bsz;
	LIS_SCALAR ***values;

	LIS_DEBUG_FUNC_IN;

	n = A->n;
	bsz = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_ilu_setVR::bsz" );
	if( bsz==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	values = (LIS_SCALAR ***)lis_malloc( n*sizeof(LIS_SCALAR **),"lis_matrix_ilu_setVR::values" );
	if( values==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR *));
		return LIS_OUT_OF_MEMORY;
	}

	A->bsz    = bsz;
	A->values = values;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_ilu_destroy"
LIS_INT lis_matrix_ilu_destroy(LIS_MATRIX_ILU A)
{
	LIS_INT i,j;

	LIS_DEBUG_FUNC_IN;

	if( lis_is_malloc(A) )
	{
		if( A->bsz )
		{
			for(i=0;i<A->n;i++)
			{
				free(A->index[i]);
				for(j=0;j<A->nnz[i];j++)
				{
					free(A->values[i][j]);
				}
				if( A->nnz[i]>0 ) free(A->values[i]);
			}
			lis_free2(5,A->bsz,A->nnz,A->index,A->values,A->nnz_ma);
		}
		else
		{
			for(i=0;i<A->n;i++)
			{
				if( A->nnz[i]>0 )
				{
					free(A->index[i]);
					free(A->value[i]);
				}
			}
			lis_free2(4,A->nnz,A->index,A->value,A->nnz_ma);
		}
		lis_free(A);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_ilu_create"
LIS_INT lis_matrix_ilu_premalloc(LIS_INT nnzrow, LIS_MATRIX_ILU A)
{
	LIS_INT i,n;
	LIS_INT *nnz_ma;

	LIS_DEBUG_FUNC_IN;

	n = A->n;

	nnz_ma = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_ilu_premalloc::nnz_ma" );
	if( nnz_ma==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0;i<n;i++)
	{
		nnz_ma[i] = nnzrow;
		A->index[i] = (LIS_INT *)malloc( nnzrow*sizeof(LIS_INT) );
		A->value[i] = (LIS_SCALAR *)malloc( nnzrow*sizeof(LIS_SCALAR) );
	}
	for(i=0;i<n;i++)
	{
		if( A->index[i]==NULL )
		{
			LIS_SETERR_MEM(nnzrow*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		if( A->value[i]==NULL )
		{
			LIS_SETERR_MEM(nnzrow*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
	}

	A->nnz_ma   = nnz_ma;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_ilu_realloc"
LIS_INT lis_matrix_ilu_realloc(LIS_INT row, LIS_INT nnz, LIS_MATRIX_ILU A)
{

	LIS_DEBUG_FUNC_IN;


	A->index[row] = (LIS_INT *)realloc(A->index[row],nnz*sizeof(LIS_INT));
	if( A->index[row]==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	A->value[row] = (LIS_SCALAR *)realloc(A->value[row],nnz*sizeof(LIS_SCALAR));
	if( A->value[row]==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matvech_ilu"
LIS_INT lis_matvech_ilu(LIS_MATRIX A, LIS_MATRIX_ILU LU, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT i,j,jj,n;
	LIS_SCALAR t,*x;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_INT	j0,j1;
		LIS_QUAD_PD tt;
	#endif

	LIS_DEBUG_FUNC_IN;

	n = LU->n;
	x = X->value;

	#ifdef USE_QUAD_PRECISION
	if( X->precision==LIS_PRECISION_DEFAULT )
	#endif
	{
		#ifdef USE_MPI
			LIS_MATVEC_SENDRECV;
		#endif
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j,jj,t)
		#endif
		for(i=0;i<n;i++)
		{
			t = 0.0;
			for(j=0;j<LU->nnz[i];j++)
			{
				jj = LU->index[i][j];
				t += conj(LU->value[i][j]) * X->value[jj];
			}
			Y->value[i] = t;
		}
	}
	#ifdef USE_QUAD_PRECISION
	else
	{
		#ifdef USE_MPI
			lis_send_recv_mp(A->commtable,X);
		#endif
		#ifndef USE_FMA2_SSE2
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,jj,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,jj,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			for(i=0;i<n;i++)
			{
				Y->value[i] = Y->value_lo[i] = 0.0;
				for(j=0;j<LU->nnz[i];j++)
				{
					jj = LU->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(Y->value[i],Y->value_lo[i],Y->value[i],Y->value_lo[i],X->value[jj],X->value_lo[jj],LU->value[i][j]);
					#else
						LIS_QUAD_FMAD_SSE2(Y->value[i],Y->value_lo[i],Y->value[i],Y->value_lo[i],X->value[jj],X->value_lo[jj],LU->value[i][j]);
					#endif
				}
			}
		#else
			#ifdef _OPENMP
			#ifndef USE_SSE2
				#pragma omp parallel for private(i,j,j0,j1,tt,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel for private(i,j,j0,j1,tt,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			#endif
			for(i=0;i<n;i++)
			{
				tt.hi[0] = tt.hi[1] = tt.lo[0] = tt.lo[1] = 0.0;
				for(j=0;j<LU->nnz[i]-1;j+=2)
				{
					j0 = LU->index[i][j];
					j1 = LU->index[i][j+1];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD2_SSE2_LDSD(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],X->value[j0],X->value_lo[j0],X->value[j1],X->value_lo[j1],LU->value[i][j]);
					#endif
				}
				#ifdef USE_SSE2
					LIS_QUAD_ADD_SSE2(Y->value[i],Y->value_lo[i],tt.hi[0],tt.lo[0],tt.hi[1],tt.lo[1]);
				#endif
				for(;j<LU->nnz[i];j++)
				{
					j0 = LU->index[i][j];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD_SSE2(Y->value[i],Y->value_lo[i],Y->value[i],Y->value_lo[i],X->value[j0],X->value_lo[j0],LU->value[i][j]);
					#endif
				}
			}
		#endif
	}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matvec_ilu"
LIS_INT lis_matvec_ilu(LIS_MATRIX A, LIS_MATRIX_ILU LU, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT i,j,jj,n,np;
	LIS_SCALAR *x;
	#ifdef _OPENMP
		LIS_INT nprocs,k;
		LIS_SCALAR t,*w;
	#endif
	#ifdef USE_QUAD_PRECISION
		LIS_INT j0,j1;
		#ifdef _OPENMP
				LIS_SCALAR *ww,*wwl;
		#endif
	#endif
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	np = A->np;
	n  = LU->n;
	x  = X->value;

	#ifdef USE_QUAD_PRECISION
	if( X->precision==LIS_PRECISION_DEFAULT )
	#endif
	{
		#ifdef USE_MPI
			LIS_MATVEC_SENDRECV;
		#endif
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			w = (LIS_SCALAR *)lis_malloc( nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_csr::w" );
			#pragma omp parallel private(i,j,k,jj,t)
			{
				k = omp_get_thread_num();
				#pragma omp for
				for(j=0;j<nprocs;j++)
				{
					memset( &w[j*np], 0, np*sizeof(LIS_SCALAR) );
				}
				#pragma omp for 
				for(i=0;i<n;i++)
				{
					for(j=0;j<LU->nnz[i];j++)
					{
						jj = k*np + LU->index[i][j];
						w[jj] += LU->value[i][j] * X->value[i];
					}
				}
				#pragma omp for 
				for(i=0;i<np;i++)
				{
					t = 0.0;
					for(j=0;j<nprocs;j++)
					{
						t += w[j*np+i];
					}
					Y->value[i] = t;
				}
			}
			lis_free(w);
		#else
			for(i=0;i<np;i++)
			{
				Y->value[i] = 0.0;
			}
			for(i=0;i<n;i++)
			{
				for(j=0;j<LU->nnz[i];j++)
				{
					jj = LU->index[i][j];
					Y->value[jj] += LU->value[i][j] * X->value[i];
				}
			}
		#endif
	}
	#ifdef USE_QUAD_PRECISION
	else
	{
		#ifdef USE_MPI
			lis_send_recv_mp(A->commtable,X);
		#endif
		#ifdef _OPENMP
			#ifndef USE_FMA2_SSE2
				nprocs = omp_get_max_threads();
				ww  = (LIS_SCALAR *)lis_malloc( 2*nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_csr_mp::ww" );
				wwl = &ww[nprocs*np];
				#ifndef USE_SSE2
					#pragma omp parallel private(i,j,jj,k,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
				#else
					#pragma omp parallel private(i,j,jj,k,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
				#endif
				{
					k = omp_get_thread_num();
					#pragma omp for
					for(j=0;j<nprocs;j++)
					{
						memset( &ww[j*np], 0, np*sizeof(LIS_SCALAR) );
						memset( &wwl[j*np], 0, np*sizeof(LIS_SCALAR) );
					}
					#pragma omp for 
					for(i=0;i<n;i++)
					{
						for(j=0;j<LU->nnz[i];j++)
						{
							jj  = k*np + LU->index[i][j];
							#ifndef USE_SSE2
							LIS_QUAD_FMAD(ww[jj],wwl[jj],ww[jj],wwl[jj],X->value[i],X->value_lo[i],LU->value[i][j]);
							#else
								LIS_QUAD_FMAD_SSE2(ww[jj],wwl[jj],ww[jj],wwl[jj],X->value[i],X->value_lo[i],LU->value[i][j]);
							#endif
						}
					}
					#pragma omp for 
					for(i=0;i<np;i++)
					{
						Y->value[i] = Y->value_lo[i] = 0.0;
						for(j=0;j<nprocs;j++)
						{
							#ifndef USE_SSE2
								LIS_QUAD_ADD(Y->value[i],Y->value_lo[i],Y->value[i],Y->value_lo[i],ww[j*np+i],wwl[j*np+i]);
							#else
								LIS_QUAD_ADD_SSE2(Y->value[i],Y->value_lo[i],Y->value[i],Y->value_lo[i],ww[j*np+i],wwl[j*np+i]);
							#endif
						}
					}
				}
				lis_free(ww);
			#else
				nprocs = omp_get_max_threads();
				ww  = (LIS_SCALAR *)lis_malloc( 2*nprocs*np*sizeof(LIS_SCALAR), "lis_matvech_csr_mp2::ww" );
				wwl = &ww[nprocs*np];
				#pragma omp parallel private(i,j,j0,j1,k,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
				{
					k = omp_get_thread_num();
					#pragma omp for
					for(j=0;j<nprocs;j++)
					{
						memset( &ww[j*np], 0, np*sizeof(LIS_SCALAR) );
						memset( &wwl[j*np], 0, np*sizeof(LIS_SCALAR) );
					}
					#pragma omp for
					for(i=0; i<n; i++)
					{
						for(j=0;j<LU->nnz[i]-1;j+=2)
						{
							j0  = k*np + LU->index[i][j];
							j1  = k*np + LU->index[i][j+1];
							#ifdef USE_SSE2
								LIS_QUAD_FMAD2_SSE2_STSD(ww[j0],wwl[j0],ww[j1],wwl[j1],ww[j0],wwl[j0],ww[j1],wwl[j1],X->value[i],X->value_lo[i],X->value[i],X->value_lo[i],LU->value[i][j]);
							#endif
						}
						for(;j<LU->nnz[i];j++)
						{
							j0  = LU->index[i][j];
							#ifdef USE_SSE2
								LIS_QUAD_FMAD_SSE2(ww[j0],wwl[j0],ww[j0],wwl[j0],X->value[i],X->value_lo[i],LU->value[i][j]);
							#endif
						}
					}
					#pragma omp for 
					for(i=0;i<np;i++)
					{
						Y->value[i] = Y->value_lo[i] = 0.0;
						for(j=0;j<nprocs;j++)
						{
							#ifdef USE_SSE2
								LIS_QUAD_ADD_SSE2(Y->value[i],Y->value_lo[i],Y->value[i],Y->value_lo[i],ww[j*np+i],wwl[j*np+i]);
							#endif
						}
					}
				}
				lis_free(ww);
			#endif
		#else
			#ifndef USE_FMA2_SSE2
				for(i=0;i<np;i++)
				{
					Y->value[i]    = 0.0;
					Y->value_lo[i] = 0.0;
				}
				for(i=0;i<n;i++)
				{
					for(j=0;j<LU->nnz[i];j++)
					{
						jj  = LU->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(Y->value[jj],Y->value_lo[jj],Y->value[jj],Y->value_lo[jj],X->value[i],X->value_lo[i],LU->value[i][j]);
						#else
							LIS_QUAD_FMAD_SSE2(Y->value[jj],Y->value_lo[jj],Y->value[jj],Y->value_lo[jj],X->value[i],X->value_lo[i],LU->value[i][j]);
						#endif
					}
				}
			#else
				for(i=0; i<np; i++)
				{
					Y->value[i]  = 0.0;
					Y->value_lo[i] = 0.0;
				}
				for(i=0; i<n; i++)
				{
					for(j=0;j<LU->nnz[i]-1;j+=2)
					{
						j0  = LU->index[i][j];
						j1  = LU->index[i][j+1];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD2_SSE2_STSD(Y->value[j0],Y->value_lo[j0],Y->value[j1],Y->value_lo[j1],Y->value[j0],Y->value_lo[j0],Y->value[j1],Y->value_lo[j1],X->value[i],X->value_lo[i],X->value[i],X->value_lo[i],LU->value[i][j]);
						#endif
					}
					for(;j<LU->nnz[i];j++)
					{
						j0  = LU->index[i][j];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD_SSE2(Y->value[j0],Y->value_lo[j0],Y->value[j0],Y->value_lo[j0],X->value[i],X->value_lo[i],LU->value[i][j]);
						#endif
					}
				}
			#endif
		#endif
	}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
