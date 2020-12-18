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

#ifndef USE_OVERLAP
void lis_matvec_jad(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k,is,ie,js,je;
	LIS_INT n,maxnzr;
	LIS_INT nprocs,my_rank;
	LIS_SCALAR *w;

	n      = A->n;
	w      = A->work;
	if( A->is_splited )
	{
		n      = A->n;
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif

		#ifdef _OPENMP
		#pragma omp parallel private(i,j,k,is,ie,js,je,my_rank)
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
			for(i=is; i<ie; i++)
			{
				y[i] = A->D->value[i]*x[i];
				w[i] = 0.0;
			}
			for(j=0;j<A->L->maxnzr;j++)
			{
				k  = is;
				js = A->L->ptr[my_rank*(A->L->maxnzr+1) + j];
				je = A->L->ptr[my_rank*(A->L->maxnzr+1) + j+1];
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=js;i<je;i++)
				{
					w[k] += A->L->value[i] * x[A->L->index[i]];
					k++;
				}
			}
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=is; i<ie; i++)
			{
				y[A->L->row[i]]  += w[i];
			}
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=is; i<ie; i++)
			{
				w[i] = 0.0;
			}
			for(j=0;j<A->U->maxnzr;j++)
			{
				k  = is;
				js = A->U->ptr[my_rank*(A->U->maxnzr+1) + j];
				je = A->U->ptr[my_rank*(A->U->maxnzr+1) + j+1];
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=js;i<je;i++)
				{
					w[k] += A->U->value[i] * x[A->U->index[i]];
					k++;
				}
			}
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=is; i<ie; i++)
			{
				y[A->U->row[i]] += w[i];
			}
		}
	}
	else
	{
		n      = A->n;
		maxnzr = A->maxnzr;
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif

		#ifdef _OPENMP
		#pragma omp parallel private(i,j,k,is,ie,js,je,my_rank)
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
			for(i=is; i<ie; i++)
			{
				w[i] = 0.0;
			}
			for(j=0;j<maxnzr;j++)
			{
				k  = is;
				js = A->ptr[my_rank*(maxnzr+1) + j];
				je = A->ptr[my_rank*(maxnzr+1) + j+1];
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=js;i<je;i++)
				{
					w[k] += A->value[i] * x[A->index[i]];
					k++;
				}
			}
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(i=is; i<ie; i++)
			{
				y[A->row[i]] = w[i];
			}
		}
	}
}
#else
void lis_matvec_jad(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k,is,ie,js,je,jj,kk;
	LIS_INT n,maxnzr,maxnzr2;
	LIS_INT nprocs,my_rank;
	LIS_SCALAR t;
	LIS_SCALAR *w;
	LIS_INT neib,inum,neibpetot,pad;
	LIS_SCALAR *ws,*wr;
	LIS_INT	*iw,err;
	LIS_COMMTABLE commtable;

	n         = A->n;
	w         = A->work;
	commtable = A->commtable;
	neibpetot = commtable->neibpetot;
	ws        = commtable->ws;
	wr        = commtable->wr;
	pad       = commtable->pad;

	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->export_ptr[neib];
		inum = commtable->export_ptr[neib+1] - is;
		for(i=is;i<is+inum;i++)
		{
			ws[i] = x[commtable->export_index[i]];
		}
		MPI_Isend(&ws[is],inum,LIS_MPI_SCALAR,commtable->neibpe[neib],0,commtable->comm,&commtable->req1[neib]);
	}

	if( A->is_splited )
	{
		n      = A->n;
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif

		#pragma omp parallel private(i,j,k,is,ie,js,je,my_rank)
		{
			#ifdef _OPENMP
				my_rank = omp_get_thread_num();
			#else
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				y[i] = A->D->value[i]*x[i];
				w[i] = 0.0;
			}
			for(j=0;j<A->L->maxnzr;j++)
			{
				k  = is;
				js = A->L->ptr[my_rank*(A->L->maxnzr+1) + j];
				je = A->L->ptr[my_rank*(A->L->maxnzr+1) + j+1];
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=js;i<je;i++)
				{
					w[k] += A->L->value[i] * x[A->L->index[i]];
					k++;
				}
			}
			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				y[A->L->row[i]]  += w[i];
			}
			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				w[i] = 0.0;
			}
			for(j=0;j<A->U->maxnzr;j++)
			{
				k  = is;
				js = A->U->ptr[my_rank*(A->U->maxnzr+1) + j];
				je = A->U->ptr[my_rank*(A->U->maxnzr+1) + j+1];
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=js;i<je;i++)
				{
					w[k] += A->U->value[i] * x[A->U->index[i]];
					k++;
				}
			}
			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				y[A->U->row[i]] += w[i];
			}
		}
	}
	else
	{
		n       = A->n;
		maxnzr  = A->maxnzr;
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif

		#pragma omp parallel private(i,j,k,is,ie,js,je,my_rank)
		{
			#ifdef _OPENMP
				my_rank = omp_get_thread_num();
			#else
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				w[i] = 0.0;
			}
			for(j=0;j<maxnzr;j++)
			{
				k  = is;
				js = A->ptr[my_rank*(maxnzr+1) + j];
				je = A->ptr[my_rank*(maxnzr+1) + j+1];
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=js;i<je;i++)
				{
					w[k] += A->value[i] * x[A->index[i]];
					k++;
				}
			}
			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				y[A->row[i]] = w[i];
			}
		}
	}
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->import_ptr[neib];
		inum = commtable->import_ptr[neib+1] - is;
		MPI_Irecv(&wr[is],inum,LIS_MPI_SCALAR,commtable->neibpe[neib],0,commtable->comm,&commtable->req2[neib]);
	}
	MPI_Waitall(neibpetot, commtable->req2, commtable->sta2);

	k = commtable->import_index[0] + pad;
	for(i=commtable->import_ptr[0];i<commtable->import_ptr[neibpetot];i++)
	{
		x[k++] = wr[i];
	}

	MPI_Waitall(neibpetot, commtable->req1, commtable->sta1);

	if( A->is_splited )
	{
		n      = A->n;
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif

		#pragma omp parallel private(i,j,k,is,ie,js,je,my_rank)
		{
			#ifdef _OPENMP
				my_rank = omp_get_thread_num();
			#else
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				y[i] = A->D->value[i]*x[i];
				w[i] = 0.0;
			}
			for(j=0;j<A->L->maxnzr;j++)
			{
				k  = is;
				js = A->L->ptr[my_rank*(A->L->maxnzr+1) + j];
				je = A->L->ptr[my_rank*(A->L->maxnzr+1) + j+1];
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=js;i<je;i++)
				{
					w[k] += A->L->value[i] * x[A->L->index[i]];
					k++;
				}
			}
			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				y[A->L->row[i]]  += w[i];
			}
			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				w[i] = 0.0;
			}
			for(j=0;j<A->U->maxnzr;j++)
			{
				k  = is;
				js = A->U->ptr[my_rank*(A->U->maxnzr+1) + j];
				je = A->U->ptr[my_rank*(A->U->maxnzr+1) + j+1];
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=js;i<je;i++)
				{
					w[k] += A->U->value[i] * x[A->U->index[i]];
					k++;
				}
			}
			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				y[A->U->row[i]] += w[i];
			}
		}
	}
	else
	{
		maxnzr2 = A->U->maxnzr;
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif

		#pragma omp parallel private(i,j,k,is,ie,js,je,my_rank)
		{
			#ifdef _OPENMP
				my_rank = omp_get_thread_num();
			#else
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				w[i] = 0.0;
			}
			for(j=0;j<maxnzr2;j++)
			{
				k  = is;
				js = A->U->ptr[my_rank*(maxnzr2+1) + j];
				je = A->U->ptr[my_rank*(maxnzr2+1) + j+1];
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=js;i<je;i++)
				{
					w[k] += A->U->value[i] * x[A->U->index[i]];
					k++;
				}
			}
			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=is; i<ie; i++)
			{
				y[A->U->row[i]] += w[i];
			}
		}
	}
}
#endif

void lis_matvech_jad(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k,js,je,jj;
	LIS_INT n,np,maxnzr;
	#ifdef _OPENMP
		LIS_INT is,ie,my_rank,nprocs;
		LIS_SCALAR t;
	#endif
	LIS_SCALAR *w;

	n      = A->n;
	np     = A->np;
	maxnzr = A->maxnzr;
	w      = A->work;

	if( A->is_splited )
	{
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0; i<n; i++)
		{
			y[i] = conj(A->D->value[i]) * x[i];
		}
		for(i=0;i<A->L->maxnzr;i++)
		{
			k  = 0;
			js = A->L->ptr[i];
			je = A->L->ptr[i+1];
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(j=js;j<je;j++)
			{
				jj = A->L->index[j];
				y[jj] += conj(A->L->value[j]) * x[A->L->row[k]];
				k++;
			}
		}
		for(i=0;i<A->U->maxnzr;i++)
		{
			k  = 0;
			js = A->U->ptr[i];
			je = A->U->ptr[i+1];
			#ifdef USE_VEC_COMP
			#pragma cdir nodep
			#pragma _NEC ivdep
			#endif
			for(j=js;j<je;j++)
			{
				jj = A->U->index[j];
				y[jj] += conj(A->U->value[j]) * x[A->U->row[k]];
				k++;
			}
		}
	}
	else
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			w = (LIS_SCALAR *)lis_malloc( nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_jad::w" );
			#pragma omp parallel private(i,j,t,is,ie,js,je,jj,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				memset( &w[my_rank*np], 0, np*sizeof(LIS_SCALAR) );

				for(j=0;j<maxnzr;j++)
				{
					k  = is;
					js = A->ptr[my_rank*(maxnzr+1) + j];
					je = A->ptr[my_rank*(maxnzr+1) + j+1];
					#ifdef USE_VEC_COMP
					#pragma cdir nodep
					#pragma _NEC ivdep
					#endif
					for(i=js;i<je;i++)
					{
						w[my_rank*np + A->index[i]] += conj(A->value[i]) * x[A->row[k]];
						k++;
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
			for(i=0; i<np; i++)
			{
				y[i] = 0.0;
			}
			for(i=0;i<maxnzr;i++)
			{
				k  = 0;
				js = A->ptr[i];
				je = A->ptr[i+1];
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(j=js;j<je;j++)
				{
					jj = A->index[j];
					y[jj] += conj(A->value[j]) * x[A->row[k]];
					k++;
				}
			}
		#endif
	}
}

void lis_matvec_jad_u4_1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT	i,j;
	LIS_INT	n,np,maxnzr;
	LIS_SCALAR *yy;
	LIS_INT	is0,is1,is2,is3;
	LIS_INT	ie0,ie1,ie2,ie3;

	n      = A->n;
	np     = A->np;
	maxnzr = A->maxnzr;
	yy     = A->work;
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0; i<np; i++)
	{
		yy[i] = 0.0;
	}
	for(j=3;j<maxnzr;j+=4)
	{
		is0 = A->ptr[j-3];
		is1 = A->ptr[j-2];
		is2 = A->ptr[j-1];
		is3 = A->ptr[j-0];
		ie3 = A->ptr[j+1-0] - A->ptr[j-0];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie3-0;i+=1)
		{
			yy[i] += A->value[is0+i]*x[A->index[is0+i]] 
			      +  A->value[is1+i]*x[A->index[is1+i]]
			      +  A->value[is2+i]*x[A->index[is2+i]]
			      +  A->value[is3+i]*x[A->index[is3+i]];
		}
		ie2 = A->ptr[j+1-1] - A->ptr[j-1];
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			yy[i] += A->value[is0+i]*x[A->index[is0+i]] 
			      +  A->value[is1+i]*x[A->index[is1+i]]
			      +  A->value[is2+i]*x[A->index[is2+i]];
		}
		ie1 = A->ptr[j+1-2] - A->ptr[j-2];
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			yy[i] += A->value[is0+i]*x[A->index[is0+i]] 
			      +  A->value[is1+i]*x[A->index[is1+i]];
		}
		ie0 = A->ptr[j+1-3] - A->ptr[j-3];
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			yy[i] += A->value[is0+i]*x[A->index[is0+i]];
		}
	}
	for(j=j-1;j<maxnzr;j+=3)
	{
		is0 = A->ptr[j-2];
		is1 = A->ptr[j-1];
		is2 = A->ptr[j-0];
		ie2 = A->ptr[j+1-0] - A->ptr[j-0];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie2-0;i+=1)
		{
			yy[i] += A->value[is0+i]*x[A->index[is0+i]] 
			      +  A->value[is1+i]*x[A->index[is1+i]]
			      +  A->value[is2+i]*x[A->index[is2+i]];
		}
		ie1 = A->ptr[j+1-1] - A->ptr[j-1];
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			yy[i] += A->value[is0+i]*x[A->index[is0+i]] 
			      +  A->value[is1+i]*x[A->index[is1+i]];
		}
		ie0 = A->ptr[j+1-2] - A->ptr[j-2];
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			yy[i] += A->value[is0+i]*x[A->index[is0+i]];
		}
	}
	for(j=j-1;j<maxnzr;j+=2)
	{
		is0 = A->ptr[j-1];
		is1 = A->ptr[j-0];
		ie1 = A->ptr[j+1-0] - A->ptr[j-0];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie1-0;i+=1)
		{
			yy[i] += A->value[is0+i]*x[A->index[is0+i]] 
			      +  A->value[is1+i]*x[A->index[is1+i]];
		}
		ie0 = A->ptr[j+1-1] - A->ptr[j-1];
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			yy[i] += A->value[is0+i]*x[A->index[is0+i]];
		}
	}
	for(j=j-1;j<maxnzr;j+=1)
	{
		is0 = A->ptr[j-0];
		ie0 = A->ptr[j+1-0] - A->ptr[j-0];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie0-0;i+=1)
		{
			yy[i] += A->value[is0+i]*x[A->index[is0+i]];
		}
	}
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0;i<n;i++)
	{
		y[A->row[i]] = yy[i];
	}
}

void lis_matvec_jad_u5_1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT	i,j;
	LIS_INT	n,np,maxnzr;
	LIS_SCALAR *yy;
	LIS_INT	is0,is1,is2,is3,is4;
	LIS_INT	ie0,ie1,ie2,ie3,ie4;
	LIS_INT	*jj0,*jj1,*jj2,*jj3,*jj4;
	LIS_SCALAR *vv0,*vv1,*vv2,*vv3,*vv4;
	LIS_INT	j00;
	LIS_INT	j10;
	LIS_INT	j20;
	LIS_INT	j30;
	LIS_INT	j40;

	n      = A->n;
	np     = A->np;
	maxnzr = A->maxnzr;
	yy     = A->work;
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0; i<np; i++)
	{
		yy[i] = 0.0;
	}
	for(j=4;j<maxnzr;j+=5)
	{
		yy   = A->work;
		is0 = A->ptr[j-4];
		ie0 = A->ptr[j+1-4] - A->ptr[j-4];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-3];
		ie1 = A->ptr[j+1-3] - A->ptr[j-3];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-2];
		ie2 = A->ptr[j+1-2] - A->ptr[j-2];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-1];
		ie3 = A->ptr[j+1-1] - A->ptr[j-1];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];
		is4 = A->ptr[j-0];
		ie4 = A->ptr[j+1-0] - A->ptr[j-0];
		jj4 = &A->index[is4];
		vv4 = &A->value[is4];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie4-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=4)
	{
		yy   = A->work;
		is0 = A->ptr[j-3];
		ie0 = A->ptr[j+1-3] - A->ptr[j-3];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-2];
		ie1 = A->ptr[j+1-2] - A->ptr[j-2];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-1];
		ie2 = A->ptr[j+1-1] - A->ptr[j-1];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-0];
		ie3 = A->ptr[j+1-0] - A->ptr[j-0];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=3)
	{
		yy   = A->work;
		is0 = A->ptr[j-2];
		ie0 = A->ptr[j+1-2] - A->ptr[j-2];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-1];
		ie1 = A->ptr[j+1-1] - A->ptr[j-1];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-0];
		ie2 = A->ptr[j+1-0] - A->ptr[j-0];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=2)
	{
		yy   = A->work;
		is0 = A->ptr[j-1];
		ie0 = A->ptr[j+1-1] - A->ptr[j-1];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-0];
		ie1 = A->ptr[j+1-0] - A->ptr[j-0];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=1)
	{
		yy   = A->work;
		is0 = A->ptr[j-0];
		ie0 = A->ptr[j+1-0] - A->ptr[j-0];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	yy = A->work;
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0;i<n;i++)
	{
		y[A->row[i]] = yy[i];
	}
}

void lis_matvec_jad_u6_1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT	i,j;
	LIS_INT	n,np,maxnzr;
	LIS_SCALAR *yy;
	LIS_INT	is0,is1,is2,is3,is4,is5;
	LIS_INT	ie0,ie1,ie2,ie3,ie4,ie5;
	LIS_INT	*jj0,*jj1,*jj2,*jj3,*jj4,*jj5;
	LIS_SCALAR *vv0,*vv1,*vv2,*vv3,*vv4,*vv5;
	LIS_INT	j00;
	LIS_INT	j10;
	LIS_INT	j20;
	LIS_INT	j30;
	LIS_INT	j40;
	LIS_INT	j50;

	n      = A->n;
	np     = A->np;
	maxnzr = A->maxnzr;
	yy     = A->work;
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0; i<np; i++)
	{
		yy[i] = 0.0;
	}
	for(j=5;j<maxnzr;j+=6)
	{
		yy   = A->work;
		is0 = A->ptr[j-5];
		ie0 = A->ptr[j+1-5] - A->ptr[j-5];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-4];
		ie1 = A->ptr[j+1-4] - A->ptr[j-4];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-3];
		ie2 = A->ptr[j+1-3] - A->ptr[j-3];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-2];
		ie3 = A->ptr[j+1-2] - A->ptr[j-2];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];
		is4 = A->ptr[j-1];
		ie4 = A->ptr[j+1-1] - A->ptr[j-1];
		jj4 = &A->index[is4];
		vv4 = &A->value[is4];
		is5 = A->ptr[j-0];
		ie5 = A->ptr[j+1-0] - A->ptr[j-0];
		jj5 = &A->index[is5];
		vv5 = &A->value[is5];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie5-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			j50 = jj5[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40] + vv5[i+0]*x[j50];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie4-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=5)
	{
		yy   = A->work;
		is0 = A->ptr[j-4];
		ie0 = A->ptr[j+1-4] - A->ptr[j-4];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-3];
		ie1 = A->ptr[j+1-3] - A->ptr[j-3];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-2];
		ie2 = A->ptr[j+1-2] - A->ptr[j-2];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-1];
		ie3 = A->ptr[j+1-1] - A->ptr[j-1];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];
		is4 = A->ptr[j-0];
		ie4 = A->ptr[j+1-0] - A->ptr[j-0];
		jj4 = &A->index[is4];
		vv4 = &A->value[is4];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie4-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=4)
	{
		yy   = A->work;
		is0 = A->ptr[j-3];
		ie0 = A->ptr[j+1-3] - A->ptr[j-3];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-2];
		ie1 = A->ptr[j+1-2] - A->ptr[j-2];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-1];
		ie2 = A->ptr[j+1-1] - A->ptr[j-1];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-0];
		ie3 = A->ptr[j+1-0] - A->ptr[j-0];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=3)
	{
		yy   = A->work;
		is0 = A->ptr[j-2];
		ie0 = A->ptr[j+1-2] - A->ptr[j-2];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-1];
		ie1 = A->ptr[j+1-1] - A->ptr[j-1];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-0];
		ie2 = A->ptr[j+1-0] - A->ptr[j-0];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=2)
	{
		yy   = A->work;
		is0 = A->ptr[j-1];
		ie0 = A->ptr[j+1-1] - A->ptr[j-1];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-0];
		ie1 = A->ptr[j+1-0] - A->ptr[j-0];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=1)
	{
		yy   = A->work;
		is0 = A->ptr[j-0];
		ie0 = A->ptr[j+1-0] - A->ptr[j-0];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	yy = A->work;
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0;i<n;i++)
	{
		y[A->row[i]] = yy[i];
	}
}

void lis_matvec_jad_u7_1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT	i,j;
	LIS_INT	n,np,maxnzr;
	LIS_SCALAR *yy;
	LIS_INT	is0,is1,is2,is3,is4,is5,is6;
	LIS_INT	ie0,ie1,ie2,ie3,ie4,ie5,ie6;
	LIS_INT	*jj0,*jj1,*jj2,*jj3,*jj4,*jj5,*jj6;
	LIS_SCALAR *vv0,*vv1,*vv2,*vv3,*vv4,*vv5,*vv6;
	LIS_INT	j00;
	LIS_INT	j10;
	LIS_INT	j20;
	LIS_INT	j30;
	LIS_INT	j40;
	LIS_INT	j50;
	LIS_INT	j60;

	n      = A->n;
	np     = A->np;
	maxnzr = A->maxnzr;
	yy     = A->work;
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0; i<np; i++)
	{
		yy[i] = 0.0;
	}
	for(j=6;j<maxnzr;j+=7)
	{
		yy   = A->work;
		is0 = A->ptr[j-6];
		ie0 = A->ptr[j+1-6] - A->ptr[j-6];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-5];
		ie1 = A->ptr[j+1-5] - A->ptr[j-5];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-4];
		ie2 = A->ptr[j+1-4] - A->ptr[j-4];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-3];
		ie3 = A->ptr[j+1-3] - A->ptr[j-3];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];
		is4 = A->ptr[j-2];
		ie4 = A->ptr[j+1-2] - A->ptr[j-2];
		jj4 = &A->index[is4];
		vv4 = &A->value[is4];
		is5 = A->ptr[j-1];
		ie5 = A->ptr[j+1-1] - A->ptr[j-1];
		jj5 = &A->index[is5];
		vv5 = &A->value[is5];
		is6 = A->ptr[j-0];
		ie6 = A->ptr[j+1-0] - A->ptr[j-0];
		jj6 = &A->index[is6];
		vv6 = &A->value[is6];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie6-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			j50 = jj5[i+0];
			j60 = jj6[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40] + vv5[i+0]*x[j50] + vv6[i+0]*x[j60];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie5-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			j50 = jj5[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40] + vv5[i+0]*x[j50];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie4-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=6)
	{
		yy   = A->work;
		is0 = A->ptr[j-5];
		ie0 = A->ptr[j+1-5] - A->ptr[j-5];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-4];
		ie1 = A->ptr[j+1-4] - A->ptr[j-4];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-3];
		ie2 = A->ptr[j+1-3] - A->ptr[j-3];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-2];
		ie3 = A->ptr[j+1-2] - A->ptr[j-2];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];
		is4 = A->ptr[j-1];
		ie4 = A->ptr[j+1-1] - A->ptr[j-1];
		jj4 = &A->index[is4];
		vv4 = &A->value[is4];
		is5 = A->ptr[j-0];
		ie5 = A->ptr[j+1-0] - A->ptr[j-0];
		jj5 = &A->index[is5];
		vv5 = &A->value[is5];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie5-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			j50 = jj5[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40] + vv5[i+0]*x[j50];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie4-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=5)
	{
		yy   = A->work;
		is0 = A->ptr[j-4];
		ie0 = A->ptr[j+1-4] - A->ptr[j-4];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-3];
		ie1 = A->ptr[j+1-3] - A->ptr[j-3];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-2];
		ie2 = A->ptr[j+1-2] - A->ptr[j-2];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-1];
		ie3 = A->ptr[j+1-1] - A->ptr[j-1];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];
		is4 = A->ptr[j-0];
		ie4 = A->ptr[j+1-0] - A->ptr[j-0];
		jj4 = &A->index[is4];
		vv4 = &A->value[is4];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie4-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=4)
	{
		yy   = A->work;
		is0 = A->ptr[j-3];
		ie0 = A->ptr[j+1-3] - A->ptr[j-3];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-2];
		ie1 = A->ptr[j+1-2] - A->ptr[j-2];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-1];
		ie2 = A->ptr[j+1-1] - A->ptr[j-1];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-0];
		ie3 = A->ptr[j+1-0] - A->ptr[j-0];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=3)
	{
		yy   = A->work;
		is0 = A->ptr[j-2];
		ie0 = A->ptr[j+1-2] - A->ptr[j-2];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-1];
		ie1 = A->ptr[j+1-1] - A->ptr[j-1];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-0];
		ie2 = A->ptr[j+1-0] - A->ptr[j-0];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=2)
	{
		yy   = A->work;
		is0 = A->ptr[j-1];
		ie0 = A->ptr[j+1-1] - A->ptr[j-1];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-0];
		ie1 = A->ptr[j+1-0] - A->ptr[j-0];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=1)
	{
		yy   = A->work;
		is0 = A->ptr[j-0];
		ie0 = A->ptr[j+1-0] - A->ptr[j-0];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	yy = A->work;
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0;i<n;i++)
	{
		y[A->row[i]] = yy[i];
	}
}

void lis_matvec_jad_u8_1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT	i,j;
	LIS_INT	n,np,maxnzr;
	LIS_SCALAR *yy;
	LIS_INT	is0,is1,is2,is3,is4,is5,is6,is7;
	LIS_INT	ie0,ie1,ie2,ie3,ie4,ie5,ie6,ie7;
	LIS_INT	*jj0,*jj1,*jj2,*jj3,*jj4,*jj5,*jj6,*jj7;
	LIS_SCALAR *vv0,*vv1,*vv2,*vv3,*vv4,*vv5,*vv6,*vv7;
	LIS_INT	j00;
	LIS_INT	j10;
	LIS_INT	j20;
	LIS_INT	j30;
	LIS_INT	j40;
	LIS_INT	j50;
	LIS_INT	j60;
	LIS_INT	j70;

	n      = A->n;
	np     = A->np;
	maxnzr = A->maxnzr;
	yy     = A->work;
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0; i<np; i++)
	{
		yy[i] = 0.0;
	}
	for(j=7;j<maxnzr;j+=8)
	{
		yy   = A->work;
		is0 = A->ptr[j-7];
		ie0 = A->ptr[j+1-7] - A->ptr[j-7];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-6];
		ie1 = A->ptr[j+1-6] - A->ptr[j-6];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-5];
		ie2 = A->ptr[j+1-5] - A->ptr[j-5];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-4];
		ie3 = A->ptr[j+1-4] - A->ptr[j-4];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];
		is4 = A->ptr[j-3];
		ie4 = A->ptr[j+1-3] - A->ptr[j-3];
		jj4 = &A->index[is4];
		vv4 = &A->value[is4];
		is5 = A->ptr[j-2];
		ie5 = A->ptr[j+1-2] - A->ptr[j-2];
		jj5 = &A->index[is5];
		vv5 = &A->value[is5];
		is6 = A->ptr[j-1];
		ie6 = A->ptr[j+1-1] - A->ptr[j-1];
		jj6 = &A->index[is6];
		vv6 = &A->value[is6];
		is7 = A->ptr[j-0];
		ie7 = A->ptr[j+1-0] - A->ptr[j-0];
		jj7 = &A->index[is7];
		vv7 = &A->value[is7];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie7-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			j50 = jj5[i+0];
			j60 = jj6[i+0];
			j70 = jj7[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40] + vv5[i+0]*x[j50] + vv6[i+0]*x[j60] + vv7[i+0]*x[j70];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie6-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			j50 = jj5[i+0];
			j60 = jj6[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40] + vv5[i+0]*x[j50] + vv6[i+0]*x[j60];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie5-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			j50 = jj5[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40] + vv5[i+0]*x[j50];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie4-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=7)
	{
		yy   = A->work;
		is0 = A->ptr[j-6];
		ie0 = A->ptr[j+1-6] - A->ptr[j-6];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-5];
		ie1 = A->ptr[j+1-5] - A->ptr[j-5];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-4];
		ie2 = A->ptr[j+1-4] - A->ptr[j-4];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-3];
		ie3 = A->ptr[j+1-3] - A->ptr[j-3];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];
		is4 = A->ptr[j-2];
		ie4 = A->ptr[j+1-2] - A->ptr[j-2];
		jj4 = &A->index[is4];
		vv4 = &A->value[is4];
		is5 = A->ptr[j-1];
		ie5 = A->ptr[j+1-1] - A->ptr[j-1];
		jj5 = &A->index[is5];
		vv5 = &A->value[is5];
		is6 = A->ptr[j-0];
		ie6 = A->ptr[j+1-0] - A->ptr[j-0];
		jj6 = &A->index[is6];
		vv6 = &A->value[is6];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie6-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			j50 = jj5[i+0];
			j60 = jj6[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40] + vv5[i+0]*x[j50] + vv6[i+0]*x[j60];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie5-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			j50 = jj5[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40] + vv5[i+0]*x[j50];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie4-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=6)
	{
		yy   = A->work;
		is0 = A->ptr[j-5];
		ie0 = A->ptr[j+1-5] - A->ptr[j-5];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-4];
		ie1 = A->ptr[j+1-4] - A->ptr[j-4];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-3];
		ie2 = A->ptr[j+1-3] - A->ptr[j-3];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-2];
		ie3 = A->ptr[j+1-2] - A->ptr[j-2];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];
		is4 = A->ptr[j-1];
		ie4 = A->ptr[j+1-1] - A->ptr[j-1];
		jj4 = &A->index[is4];
		vv4 = &A->value[is4];
		is5 = A->ptr[j-0];
		ie5 = A->ptr[j+1-0] - A->ptr[j-0];
		jj5 = &A->index[is5];
		vv5 = &A->value[is5];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie5-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			j50 = jj5[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40] + vv5[i+0]*x[j50];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie4-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=5)
	{
		yy   = A->work;
		is0 = A->ptr[j-4];
		ie0 = A->ptr[j+1-4] - A->ptr[j-4];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-3];
		ie1 = A->ptr[j+1-3] - A->ptr[j-3];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-2];
		ie2 = A->ptr[j+1-2] - A->ptr[j-2];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-1];
		ie3 = A->ptr[j+1-1] - A->ptr[j-1];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];
		is4 = A->ptr[j-0];
		ie4 = A->ptr[j+1-0] - A->ptr[j-0];
		jj4 = &A->index[is4];
		vv4 = &A->value[is4];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie4-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			j40 = jj4[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30] + vv4[i+0]*x[j40];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=4)
	{
		yy   = A->work;
		is0 = A->ptr[j-3];
		ie0 = A->ptr[j+1-3] - A->ptr[j-3];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-2];
		ie1 = A->ptr[j+1-2] - A->ptr[j-2];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-1];
		ie2 = A->ptr[j+1-1] - A->ptr[j-1];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];
		is3 = A->ptr[j-0];
		ie3 = A->ptr[j+1-0] - A->ptr[j-0];
		jj3 = &A->index[is3];
		vv3 = &A->value[is3];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie3-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			j30 = jj3[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20] + vv3[i+0]*x[j30];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=3)
	{
		yy   = A->work;
		is0 = A->ptr[j-2];
		ie0 = A->ptr[j+1-2] - A->ptr[j-2];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-1];
		ie1 = A->ptr[j+1-1] - A->ptr[j-1];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];
		is2 = A->ptr[j-0];
		ie2 = A->ptr[j+1-0] - A->ptr[j-0];
		jj2 = &A->index[is2];
		vv2 = &A->value[is2];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie2-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			j20 = jj2[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10] + vv2[i+0]*x[j20];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=2)
	{
		yy   = A->work;
		is0 = A->ptr[j-1];
		ie0 = A->ptr[j+1-1] - A->ptr[j-1];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];
		is1 = A->ptr[j-0];
		ie1 = A->ptr[j+1-0] - A->ptr[j-0];
		jj1 = &A->index[is1];
		vv1 = &A->value[is1];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie1-0;i+=1)
		{
			j00 = jj0[i+0];
			j10 = jj1[i+0];
			yy[i+0] += vv0[i+0]*x[j00] + vv1[i+0]*x[j10];
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	for(j=j-1;j<maxnzr;j+=1)
	{
		yy   = A->work;
		is0 = A->ptr[j-0];
		ie0 = A->ptr[j+1-0] - A->ptr[j-0];
		jj0 = &A->index[is0];
		vv0 = &A->value[is0];

		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=0;i<ie0-0;i+=1)
		{
			j00 = jj0[i+0];
			yy[i+0] += vv0[i+0]*x[j00];
		}
	}
	yy = A->work;
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0;i<n;i++)
	{
		y[A->row[i]] = yy[i];
	}
}
