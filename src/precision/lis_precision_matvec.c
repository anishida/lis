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


#ifdef USE_QUAD_PRECISION
void lis_matvec_csr_mp(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT	i,j,n;
	LIS_INT	is,ie,j0;
	LIS_INT	*jj0;
	LIS_SCALAR *vv0;
	LIS_SCALAR *x,*y,*xl,*yl;
	LIS_QUAD_DECLAR;


	n     = A->n;
	x     = X->value;
	y     = Y->value;
	xl    = X->value_lo;
	yl    = Y->value_lo;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#ifndef USE_SSE2
			#pragma omp parallel for private(i,j,is,ie,j0,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		#else
			#pragma omp parallel for private(i,j,is,ie,j0,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
		#endif
		#endif
		for(i=0;i<n;i++)
		{
			#ifndef USE_SSE2
				LIS_QUAD_MULD(y[i],yl[i],x[i],xl[i],A->D->value[i]);
			#else
				LIS_QUAD_MULD_SSE2(y[i],yl[i],x[i],xl[i],A->D->value[i]);
			#endif
			is = A->L->ptr[i];
			ie = A->L->ptr[i+1];
			for(j=is;j<ie-0;j+=1)
			{
				j0 = A->L->index[j+0];
				#ifndef USE_SSE2
					LIS_QUAD_FMAD(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->L->value[j]);
				#else
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->L->value[j]);
				#endif
			}
			is = A->U->ptr[i];
			ie = A->U->ptr[i+1];
			for(j=is;j<ie-0;j+=1)
			{
				j0 = A->U->index[j+0];
				#ifndef USE_SSE2
					LIS_QUAD_FMAD(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->U->value[j]);
				#else
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->U->value[j]);
				#endif
			}
		}
	}
	else
	{
		jj0 = A->index;
		vv0 = A->value;
		#ifdef _OPENMP
		#ifndef USE_SSE2
			#pragma omp parallel for private(i,j,is,ie,j0,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		#else
			#pragma omp parallel for private(i,j,is,ie,j0,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
		#endif
		#endif
		for(i=0;i<n;i++)
		{
			y[i] = yl[i] = 0.0;

			is = A->ptr[i];
			ie = A->ptr[i+1];
			for(j=is;j<ie-0;j+=1)
			{
				j0 = jj0[j+0];
				#ifndef USE_SSE2
					LIS_QUAD_FMAD(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],vv0[j]);
				#else
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],vv0[j]);
				#endif
			}
		}
	}
}

void lis_matvec_csr_mp2(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT i,j,n;
	LIS_INT	is,ie;
	LIS_INT	j0,j1;
	LIS_INT	*jj0;
	LIS_SCALAR *vv0;
	LIS_SCALAR *x,*y,*xl,*yl;
	LIS_QUAD_PD tt;
	LIS_QUAD_DECLAR;

	n     = A->n;
	x     = X->value;
	y     = Y->value;
	xl    = X->value_lo;
	yl    = Y->value_lo;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#ifndef USE_SSE2
			#pragma omp parallel for private(i,j,is,ie,j0,j1,tt,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		#else
			#pragma omp parallel for private(i,j,is,ie,j0,j1,tt,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
		#endif
		#endif
		for(i=0;i<n;i++)
		{
			#ifndef USE_SSE2
				LIS_QUAD_MULD(y[i],yl[i],x[i],xl[i],A->D->value[i]);
			#else
				LIS_QUAD_MULD_SSE2(y[i],yl[i],x[i],xl[i],A->D->value[i]);
			#endif

			tt.hi[0] = tt.hi[1] = tt.lo[0] = tt.lo[1] = 0.0;
			is = A->L->ptr[i];
			ie = A->L->ptr[i+1];
			for(j=is;j<ie-1;j+=2)
			{
				j0 = A->L->index[j+0];
				j1 = A->L->index[j+1];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD2_SSE2_LDSD(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],x[j0],xl[j0],x[j1],xl[j1],A->L->value[j]);
				#endif
			}
			for(;j<ie;j++)
			{
				j0 = A->L->index[j+0];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->L->value[j]);
				#endif
			}
			is = A->U->ptr[i];
			ie = A->U->ptr[i+1];
			for(j=is;j<ie-1;j+=2)
			{
				j0 = A->U->index[j+0];
				j1 = A->U->index[j+1];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD2_SSE2_LDSD(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],x[j0],xl[j0],x[j1],xl[j1],A->U->value[j]);
				#endif
			}
			for(;j<ie;j++)
			{
				j0 = A->U->index[j+0];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->U->value[j]);
				#endif
			}
			#ifdef USE_SSE2
				LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],tt.hi[0],tt.lo[0]);
				LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],tt.hi[1],tt.lo[1]);
			#endif
		}
	}
	else
	{
		jj0 = A->index;
		vv0 = A->value;
		#ifdef _OPENMP
		#ifndef USE_SSE2
			#pragma omp parallel for private(i,j,is,ie,j0,j1,tt,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		#else
			#pragma omp parallel for private(i,j,is,ie,j0,j1,tt,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
		#endif
		#endif
		for(i=0;i<n;i++)
		{
			tt.hi[0] = tt.hi[1] = tt.lo[0] = tt.lo[1] = 0.0;

			is = A->ptr[i];
			ie = A->ptr[i+1];
			for(j=is;j<ie-1;j+=2)
			{
				j0 = jj0[j+0];
				j1 = jj0[j+1];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD2_SSE2_LDSD(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],x[j0],xl[j0],x[j1],xl[j1],vv0[j]);
				#endif
			}
			#ifdef USE_SSE2
				LIS_QUAD_ADD_SSE2(y[i],yl[i],tt.hi[0],tt.lo[0],tt.hi[1],tt.lo[1]);
			#endif
			for(;j<ie;j++)
			{
				j0 = jj0[j+0];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],vv0[j]);
				#endif
			}
		}
	}
}

void lis_matvech_csr_mp(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT	i,j,js,je,jj;
	LIS_INT	n,np;
	LIS_QUAD_PTR tt0;
	LIS_SCALAR *x,*y,*xl,*yl;
	#ifdef _OPENMP
		LIS_INT k,nprocs;
		LIS_SCALAR *ww,*wwl;
	#endif
	LIS_QUAD_DECLAR;

	n    = A->n;
	np   = A->np;
	x     = X->value;
	y     = Y->value;
	xl    = X->value_lo;
	yl    = Y->value_lo;
	tt0.hi = &X->work[0];
	tt0.lo = &X->work[1];
	if( A->is_splited )
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			ww  = (LIS_SCALAR *)lis_malloc( 2*nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_csr_mp::ww" );
			wwl = &ww[nprocs*np];
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,js,je,jj,k,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,js,je,jj,k,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
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
				for(i=0; i<n; i++)
				{
					js = A->L->ptr[i];
					je = A->L->ptr[i+1];
					for(j=js;j<je;j++)
					{
						jj  = k*np+A->L->index[j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->L->value[j]);
						#else
							LIS_QUAD_FMAD_SSE2(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->L->value[j]);
						#endif
					}
					js = A->U->ptr[i];
					je = A->U->ptr[i+1];
					for(j=js;j<je;j++)
					{
						jj  = k*np+A->U->index[j];
						#ifndef USE_SSE2
						LIS_QUAD_FMAD(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->U->value[j]);
						#else
							LIS_QUAD_FMAD_SSE2(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->U->value[j]);
						#endif
					}
				}
				#pragma omp for 
				for(i=0;i<np;i++)
				{
					#ifndef USE_SSE2
					LIS_QUAD_MULD(y[i],yl[i],x[i],xl[i],A->D->value[i]);
					#else
						LIS_QUAD_MULD_SSE2(y[i],yl[i],x[i],xl[i],A->D->value[i]);
					#endif
					for(j=0;j<nprocs;j++)
					{
						#ifndef USE_SSE2
							LIS_QUAD_ADD(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#else
							LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#endif
					}
				}
			}
			lis_free(ww);
		#else
			for(i=0; i<np; i++)
			{
				#ifndef USE_SSE2
				LIS_QUAD_MULD(y[i],yl[i],x[i],xl[i],A->D->value[i]);
				#else
				LIS_QUAD_MULD_SSE2(y[i],yl[i],x[i],xl[i],A->D->value[i]);
				#endif
			}
			for(i=0; i<n; i++)
			{
				js = A->L->ptr[i];
				je = A->L->ptr[i+1];
				for(j=js;j<je;j++)
				{
					jj  = A->L->index[j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(y[jj],yl[jj],y[jj],yl[jj],x[i],xl[i],A->L->value[j]);
					#else
						LIS_QUAD_FMAD_SSE2(y[jj],yl[jj],y[jj],yl[jj],x[i],xl[i],A->L->value[j]);
					#endif
				}
				js = A->U->ptr[i];
				je = A->U->ptr[i+1];
				for(j=js;j<je;j++)
				{
					jj  = A->U->index[j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(y[jj],yl[jj],y[jj],yl[jj],x[i],xl[i],A->U->value[j]);
					#else
						LIS_QUAD_FMAD_SSE2(y[jj],yl[jj],y[jj],yl[jj],x[i],xl[i],A->U->value[j]);
					#endif
				}
			}
		#endif
	}
	else
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			ww  = (LIS_SCALAR *)lis_malloc( 2*nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_csr_mp::ww" );
			wwl = &ww[nprocs*np];
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,js,je,jj,k,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,js,je,jj,k,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
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
				for(i=0; i<n; i++)
				{
					js = A->ptr[i];
					je = A->ptr[i+1];
					for(j=js;j<je;j++)
					{
						jj  = k*np+A->index[j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->value[j]);
						#else
							LIS_QUAD_FMAD_SSE2(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->value[j]);
						#endif
					}
				}
				#pragma omp for 
				for(i=0;i<np;i++)
				{
					y[i] = yl[i] = 0.0;
					for(j=0;j<nprocs;j++)
					{
						#ifndef USE_SSE2
							LIS_QUAD_ADD(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#else
							LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#endif
					}
				}
			}
			lis_free(ww);
		#else
			for(i=0; i<np; i++)
			{
				y[i]  = 0.0;
				yl[i] = 0.0;
			}
			for(i=0; i<n; i++)
			{
				js = A->ptr[i];
				je = A->ptr[i+1];
				tt0.hi[0] = x[i];
				tt0.lo[0] = xl[i];
				for(j=js;j<je;j++)
				{
					jj  = A->index[j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(y[jj],yl[jj],y[jj],yl[jj],tt0.hi[0],tt0.lo[0],A->value[j]);
					#else
						LIS_QUAD_FMAD_SSE2(y[jj],yl[jj],y[jj],yl[jj],tt0.hi[0],tt0.lo[0],A->value[j]);
					#endif
				}
			}
		#endif
	}
}

void lis_matvech_csr_mp2(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT i,j,js,je,j0,j1;
	LIS_INT	n,np;
	LIS_QUAD_PTR tt0;
	LIS_SCALAR *x,*y,*xl,*yl;
	#ifdef _OPENMP
		LIS_INT k,nprocs;
		LIS_SCALAR *ww,*wwl;
	#endif
	LIS_QUAD_DECLAR;

	n    = A->n;
	np   = A->np;
	x     = X->value;
	y     = Y->value;
	xl    = X->value_lo;
	yl    = Y->value_lo;
	tt0.hi = &X->work[0];
	tt0.lo = &X->work[2];
	if( A->is_splited )
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			ww  = (LIS_SCALAR *)lis_malloc( 2*nprocs*np*sizeof(LIS_SCALAR), "lis_matvech_csr_mp2::ww" );
			wwl = &ww[nprocs*np];
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,js,je,j0,j1,k,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,js,je,j0,j1,k,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			{
				k = omp_get_thread_num();
				#pragma omp for
				for(i=0; i<np; i++)
				{
					#ifndef USE_SSE2
					LIS_QUAD_MULD(ww[i],wwl[i],x[i],xl[i],A->D->value[i]);
					#else
					LIS_QUAD_MULD_SSE2(ww[i],wwl[i],x[i],xl[i],A->D->value[i]);
					#endif
				}
				#pragma omp for 
				for(i=0; i<n; i++)
				{
					js = A->L->ptr[i];
					je = A->L->ptr[i+1];
					for(j=js;j<je-1;j+=2)
					{
						j0  = k*np + A->L->index[j];
						j1  = k*np + A->L->index[j+1];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD2_SSE2_STSD(ww[j0],wwl[j0],ww[j1],wwl[j1],ww[j0],wwl[j0],ww[j1],wwl[j1],x[i],xl[i],x[i],xl[i],A->L->value[j]);
						#endif
					}
					for(;j<je;j++)
					{
						j0  = A->L->index[j];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD_SSE2(ww[j0],wwl[j0],ww[j0],wwl[j0],x[i],xl[i],A->L->value[j]);
						#endif
					}
					js = A->U->ptr[i];
					je = A->U->ptr[i+1];
					for(j=js;j<je-1;j+=2)
					{
						j0  = k*np + A->U->index[j];
						j1  = k*np + A->U->index[j+1];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD2_SSE2_STSD(ww[j0],wwl[j0],ww[j1],wwl[j1],ww[j0],wwl[j0],ww[j1],wwl[j1],x[i],xl[i],x[i],xl[i],A->U->value[j]);
						#endif
					}
					for(;j<je;j++)
					{
						j0  = A->U->index[j];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD_SSE2(ww[j0],wwl[j0],ww[j0],wwl[j0],x[i],xl[i],A->U->value[j]);
						#endif
					}
				}
				#pragma omp for 
				for(i=0;i<np;i++)
				{
					y[i] = yl[i] = 0.0;
					for(j=0;j<nprocs;j++)
					{
						#ifdef USE_SSE2
							LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#endif
					}
				}
			}
			lis_free(ww);
		#else
			for(i=0; i<np; i++)
			{
				#ifndef USE_SSE2
			  	LIS_QUAD_MULD(y[i],yl[i],x[i],xl[i],A->D->value[i]);
				#else
				LIS_QUAD_MULD_SSE2(y[i],yl[i],x[i],xl[i],A->D->value[i]);
				#endif
			}
			for(i=0; i<n; i++)
			{
				js = A->L->ptr[i];
				je = A->L->ptr[i+1];
				for(j=js;j<je-1;j+=2)
				{
					j0  = A->L->index[j];
					j1  = A->L->index[j+1];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD2_SSE2_STSD(y[j0],yl[j0],y[j1],yl[j1],y[j0],yl[j0],y[j1],yl[j1],x[i],xl[i],x[i],xl[i],A->L->value[j]);
					#endif
				}
				for(;j<je;j++)
				{
					j0  = A->L->index[j];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD_SSE2(y[j0],yl[j0],y[j0],yl[j0],tt0.hi[0],tt0.lo[0],A->L->value[j]);
					#endif
				}
				js = A->U->ptr[i];
				je = A->U->ptr[i+1];
				for(j=js;j<je-1;j+=2)
				{
					j0  = A->U->index[j];
					j1  = A->U->index[j+1];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD2_SSE2_STSD(y[j0],yl[j0],y[j1],yl[j1],y[j0],yl[j0],y[j1],yl[j1],x[i],xl[i],x[i],xl[i],A->U->value[j]);
					#endif
				}
				for(;j<je;j++)
				{
					j0  = A->U->index[j];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD_SSE2(y[j0],yl[j0],y[j0],yl[j0],tt0.hi[0],tt0.lo[0],A->U->value[j]);
					#endif
				}
			}
		#endif
	}
	else
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			ww  = (LIS_SCALAR *)lis_malloc( 2*nprocs*np*sizeof(LIS_SCALAR), "lis_matvech_csr_mp2::ww" );
			wwl = &ww[nprocs*np];
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,js,je,j0,j1,k,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,js,je,j0,j1,k,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
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
				for(i=0; i<n; i++)
				{
					js = A->ptr[i];
					je = A->ptr[i+1];
					for(j=js;j<je-1;j+=2)
					{
						j0  = k*np + A->index[j];
						j1  = k*np + A->index[j+1];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD2_SSE2_STSD(ww[j0],wwl[j0],ww[j1],wwl[j1],ww[j0],wwl[j0],ww[j1],wwl[j1],x[i],xl[i],x[i],xl[i],A->value[j]);
						#endif
					}
					for(;j<je;j++)
					{
						j0  = A->index[j];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD_SSE2(ww[j0],wwl[j0],ww[j0],wwl[j0],x[i],xl[i],A->value[j]);
						#endif
					}
				}
				#pragma omp for 
				for(i=0;i<np;i++)
				{
					y[i] = yl[i] = 0.0;
					for(j=0;j<nprocs;j++)
					{
						#ifdef USE_SSE2
							LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#endif
					}
				}
			}
			lis_free(ww);
		#else
			for(i=0; i<np; i++)
			{
				y[i]  = 0.0;
				yl[i] = 0.0;
			}
			for(i=0; i<n; i++)
			{
				js = A->ptr[i];
				je = A->ptr[i+1];
				for(j=js;j<je-1;j+=2)
				{
					j0  = A->index[j];
					j1  = A->index[j+1];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD2_SSE2_STSD(y[j0],yl[j0],y[j1],yl[j1],y[j0],yl[j0],y[j1],yl[j1],x[i],xl[i],x[i],xl[i],A->value[j]);
					#endif
				}
				for(;j<je;j++)
				{
					j0  = A->index[j];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD_SSE2(y[j0],yl[j0],y[j0],yl[j0],x[i],xl[i],A->value[j]);
					#endif
				}
			}
		#endif
	}
}
#endif

#ifdef USE_QUAD_PRECISION
void lis_matvec_csc_mp(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT	i,j,js,je,jj;
	LIS_INT	n,np;
	LIS_QUAD_PTR tt0;
	LIS_SCALAR *x,*y,*xl,*yl;
	#ifdef _OPENMP
		LIS_INT k,nprocs;
		LIS_SCALAR *ww,*wwl;
	#endif
	LIS_QUAD_DECLAR;

	n    = A->n;
	np   = A->np;
	x     = X->value;
	y     = Y->value;
	xl    = X->value_lo;
	yl    = Y->value_lo;
	tt0.hi = &X->work[0];
	tt0.lo = &X->work[1];
	if( A->is_splited )
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			ww  = (LIS_SCALAR *)lis_malloc( 2*nprocs*np*sizeof(LIS_SCALAR),"lis_matvec_csr_mp::ww" );
			wwl = &ww[nprocs*np];
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,js,je,jj,k,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,js,je,jj,k,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
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
				for(i=0; i<np; i++)
				{
					js = A->L->ptr[i];
					je = A->L->ptr[i+1];
					for(j=js;j<je;j++)
					{
						jj  = k*np+A->L->index[j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->L->value[j]);
						#else
							LIS_QUAD_FMAD_SSE2(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->L->value[j]);
						#endif
					}
					js = A->U->ptr[i];
					je = A->U->ptr[i+1];
					for(j=js;j<je;j++)
					{
						jj  = k*np+A->U->index[j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->U->value[j]);
						#else
							LIS_QUAD_FMAD_SSE2(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->U->value[j]);
						#endif
					}
				}
				#pragma omp for 
				for(i=0;i<n;i++)
				{
					#ifndef USE_SSE2
						LIS_QUAD_MULD(y[i],yl[i],x[i],xl[i],A->D->value[i]);
					#else
						LIS_QUAD_MULD_SSE2(y[i],yl[i],x[i],xl[i],A->D->value[i]);
					#endif
					for(j=0;j<nprocs;j++)
					{
						#ifndef USE_SSE2
							LIS_QUAD_ADD(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#else
							LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#endif
					}
				}
			}
			lis_free(ww);
		#else
			for(i=0; i<n; i++)
			{
				#ifndef USE_SSE2
					LIS_QUAD_MULD(y[i],yl[i],x[i],xl[i],A->D->value[i]);
				#else
					LIS_QUAD_MULD_SSE2(y[i],yl[i],x[i],xl[i],A->D->value[i]);
				#endif
			}
			for(i=0; i<np; i++)
			{
				js = A->L->ptr[i];
				je = A->L->ptr[i+1];
				for(j=js;j<je;j++)
				{
					jj  = A->L->index[j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(y[jj],yl[jj],y[jj],yl[jj],x[i],xl[i],A->L->value[j]);
					#else
						LIS_QUAD_FMAD_SSE2(y[jj],yl[jj],y[jj],yl[jj],x[i],xl[i],A->L->value[j]);
					#endif
				}
				js = A->U->ptr[i];
				je = A->U->ptr[i+1];
				for(j=js;j<je;j++)
				{
					jj  = A->U->index[j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(y[jj],yl[jj],y[jj],yl[jj],x[i],xl[i],A->U->value[j]);
					#else
						LIS_QUAD_FMAD_SSE2(y[jj],yl[jj],y[jj],yl[jj],x[i],xl[i],A->U->value[j]);
					#endif
				}
			}
		#endif
	}
	else
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			ww  = (LIS_SCALAR *)lis_malloc( 2*nprocs*np*sizeof(LIS_SCALAR),"lis_matvec_csr_mp::ww" );
			wwl = &ww[nprocs*np];
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,js,je,jj,k,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,js,je,jj,k,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
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
				for(i=0; i<np; i++)
				{
					js = A->ptr[i];
					je = A->ptr[i+1];
					for(j=js;j<je;j++)
					{
						jj  = k*np+A->L->index[j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->L->value[j]);
						#else
							LIS_QUAD_FMAD_SSE2(ww[jj],wwl[jj],ww[jj],wwl[jj],x[i],xl[i],A->L->value[j]);
						#endif
					}
				}
				#pragma omp for 
				for(i=0;i<np;i++)
				{
					y[i] = yl[i] = 0.0;
					for(j=0;j<nprocs;j++)
					{
						#ifndef USE_SSE2
							LIS_QUAD_ADD(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#else
							LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#endif
					}
				}
			}
			lis_free(ww);
		#else
			for(i=0; i<n; i++)
			{
				y[i]  = 0.0;
				yl[i] = 0.0;
			}
			for(i=0; i<np; i++)
			{
				js = A->ptr[i];
				je = A->ptr[i+1];
				tt0.hi[0] = x[i];
				tt0.lo[0] = xl[i];
				for(j=js;j<je;j++)
				{
					jj  = A->index[j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(y[jj],yl[jj],y[jj],yl[jj],tt0.hi[0],tt0.lo[0],A->value[j]);
					#else
						LIS_QUAD_FMAD_SSE2(y[jj],yl[jj],y[jj],yl[jj],tt0.hi[0],tt0.lo[0],A->value[j]);
					#endif
				}
			}
		#endif
	}
}

void lis_matvec_csc_mp2(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT	i,j,js,je,j0,j1;
	LIS_INT	n,np;
	LIS_QUAD_PTR tt0;
	LIS_SCALAR *x,*y,*xl,*yl;
	#ifdef _OPENMP
		LIS_INT k,nprocs;
		LIS_SCALAR *ww,*wwl;
	#endif
	LIS_QUAD_DECLAR;

	n    = A->n;
	np   = A->np;
	x     = X->value;
	y     = Y->value;
	xl    = X->value_lo;
	yl    = Y->value_lo;
	tt0.hi = &X->work[0];
	tt0.lo = &X->work[2];
	if( A->is_splited )
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			ww  = (LIS_SCALAR *)lis_malloc( 2*nprocs*np*sizeof(LIS_SCALAR), "lis_matvec_csr_mp2::ww" );
			wwl = &ww[nprocs*np];
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,js,je,j0,j1,k,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,js,je,j0,j1,k,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
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
				for(i=0; i<np; i++)
				{
					js = A->ptr[i];
					je = A->ptr[i+1];
					for(j=js;j<je-1;j+=2)
					{
						j0  = k*np + A->index[j];
						j1  = k*np + A->index[j+1];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD2_SSE2_STSD(ww[j0],wwl[j0],ww[j1],wwl[j1],ww[j0],wwl[j0],ww[j1],wwl[j1],x[i],xl[i],x[i],xl[i],A->value[j]);
						#endif
					}
					for(;j<je;j++)
					{
						j0  = A->index[j];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD_SSE2(ww[j0],wwl[j0],ww[j0],wwl[j0],x[i],xl[i],A->value[j]);
						#endif
					}
				}
				#pragma omp for 
				for(i=0;i<n;i++)
				{
					y[i] = yl[i] = 0.0;
					for(j=0;j<nprocs;j++)
					{
						#ifdef USE_SSE2
							LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#endif
					}
				}
			}
			lis_free(ww);
		#else
			for(i=0; i<n; i++)
			{
				#ifndef USE_SSE2
					LIS_QUAD_MULD(y[i],yl[i],x[i],xl[i],A->D->value[i]);
				#else
					LIS_QUAD_MULD_SSE2(y[i],yl[i],x[i],xl[i],A->D->value[i]);
				#endif
			}
			for(i=0; i<np; i++)
			{
				js = A->L->ptr[i];
				je = A->L->ptr[i+1];
				for(j=js;j<je-1;j+=2)
				{
					j0  = A->L->index[j];
					j1  = A->L->index[j+1];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD2_SSE2_STSD(y[j0],yl[j0],y[j1],yl[j1],y[j0],yl[j0],y[j1],yl[j1],x[i],xl[i],x[i],xl[i],A->L->value[j]);
					#endif
				}
				for(;j<je;j++)
				{
					j0  = A->L->index[j];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD_SSE2(y[j0],yl[j0],y[j0],yl[j0],tt0.hi[0],tt0.lo[0],A->L->value[j]);
					#endif
				}
				js = A->U->ptr[i];
				je = A->U->ptr[i+1];
				for(j=js;j<je-1;j+=2)
				{
					j0  = A->U->index[j];
					j1  = A->U->index[j+1];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD2_SSE2_STSD(y[j0],yl[j0],y[j1],yl[j1],y[j0],yl[j0],y[j1],yl[j1],x[i],xl[i],x[i],xl[i],A->U->value[j]);
					#endif
				}
				for(;j<je;j++)
				{
					j0  = A->U->index[j];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD_SSE2(y[j0],yl[j0],y[j0],yl[j0],tt0.hi[0],tt0.lo[0],A->U->value[j]);
					#endif
				}
			}
		#endif
	}
	else
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			ww  = (LIS_SCALAR *)lis_malloc( 2*nprocs*np*sizeof(LIS_SCALAR), "lis_matvec_csr_mp2::ww" );
			wwl = &ww[nprocs*np];
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,js,je,j0,j1,k,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,js,je,j0,j1,k,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
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
				for(i=0; i<np; i++)
				{
					js = A->ptr[i];
					je = A->ptr[i+1];
					for(j=js;j<je-1;j+=2)
					{
						j0  = k*np + A->index[j];
						j1  = k*np + A->index[j+1];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD2_SSE2_STSD(ww[j0],wwl[j0],ww[j1],wwl[j1],ww[j0],wwl[j0],ww[j1],wwl[j1],x[i],xl[i],x[i],xl[i],A->value[j]);
						#endif
					}
					for(;j<je;j++)
					{
						j0  = A->index[j];
						#ifdef USE_SSE2
							LIS_QUAD_FMAD_SSE2(ww[j0],wwl[j0],ww[j0],wwl[j0],x[i],xl[i],A->value[j]);
						#endif
					}
				}
				#pragma omp for 
				for(i=0;i<n;i++)
				{
					y[i] = yl[i] = 0.0;
					for(j=0;j<nprocs;j++)
					{
						#ifdef USE_SSE2
							LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],ww[j*np+i],wwl[j*np+i]);
						#endif
					}
				}
			}
			lis_free(ww);
		#else
			for(i=0; i<n; i++)
			{
				y[i]  = 0.0;
				yl[i] = 0.0;
			}
			for(i=0; i<np; i++)
			{
				js = A->ptr[i];
				je = A->ptr[i+1];
				for(j=js;j<je-1;j+=2)
				{
					j0  = A->index[j];
					j1  = A->index[j+1];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD2_SSE2_STSD(y[j0],yl[j0],y[j1],yl[j1],y[j0],yl[j0],y[j1],yl[j1],x[i],xl[i],x[i],xl[i],A->value[j]);
					#endif
				}
				for(;j<je;j++)
				{
					j0  = A->index[j];
					#ifdef USE_SSE2
						LIS_QUAD_FMAD_SSE2(y[j0],yl[j0],y[j0],yl[j0],x[i],xl[i],A->value[j]);
					#endif
				}
			}
		#endif
	}
}

void lis_matvech_csc_mp(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT	i,j,np;
	LIS_INT	is,ie,j0;
	LIS_INT	*jj0;
	LIS_SCALAR *vv0;
	LIS_SCALAR *x,*y,*xl,*yl;
	LIS_QUAD_DECLAR;


	np    = A->np;
	x     = X->value;
	y     = Y->value;
	xl    = X->value_lo;
	yl    = Y->value_lo;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#ifndef USE_SSE2
			#pragma omp parallel for private(i,j,is,ie,j0,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		#else
			#pragma omp parallel for private(i,j,is,ie,j0,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
		#endif
		#endif
		for(i=0;i<np;i++)
		{
			#ifndef USE_SSE2
		  		LIS_QUAD_MULD(y[i],yl[i],x[i],xl[i],A->D->value[i]);
			#else
				LIS_QUAD_MULD_SSE2(y[i],yl[i],x[i],xl[i],A->D->value[i]);
			#endif
			is = A->L->ptr[i];
			ie = A->L->ptr[i+1];
			for(j=is;j<ie-0;j+=1)
			{
				j0 = A->L->index[j+0];
				#ifndef USE_SSE2
					LIS_QUAD_FMAD(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->L->value[j]);
				#else
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->L->value[j]);
				#endif
			}
			is = A->U->ptr[i];
			ie = A->U->ptr[i+1];
			for(j=is;j<ie-0;j+=1)
			{
				j0 = A->U->index[j+0];
				#ifndef USE_SSE2
					LIS_QUAD_FMAD(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->U->value[j]);
				#else
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->U->value[j]);
				#endif
			}
		}
	}
	else
	{
		jj0 = A->index;
		vv0 = A->value;
		#ifdef _OPENMP
		#ifndef USE_SSE2
			#pragma omp parallel for private(i,j,is,ie,j0,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		#else
			#pragma omp parallel for private(i,j,is,ie,j0,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
		#endif
		#endif
		for(i=0;i<np;i++)
		{
			y[i] = yl[i] = 0.0;

			is = A->ptr[i];
			ie = A->ptr[i+1];
			for(j=is;j<ie-0;j+=1)
			{
				j0 = jj0[j+0];
				#ifndef USE_SSE2
					LIS_QUAD_FMAD(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],vv0[j]);
				#else
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],vv0[j]);
				#endif
			}
		}
	}
}

void lis_matvech_csc_mp2(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT	i,j,np;
	LIS_INT	is,ie;
	LIS_INT	j0,j1;
	LIS_INT	*jj0;
	LIS_SCALAR *vv0;
	LIS_SCALAR *x,*y,*xl,*yl;
	LIS_QUAD_PD tt;
	LIS_QUAD_DECLAR;

	np    = A->np;
	x     = X->value;
	y     = Y->value;
	xl    = X->value_lo;
	yl    = Y->value_lo;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#ifndef USE_SSE2
			#pragma omp parallel for private(i,j,is,ie,j0,j1,tt,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		#else
			#pragma omp parallel for private(i,j,is,ie,j0,j1,tt,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
		#endif
		#endif
		for(i=0;i<np;i++)
		{
			#ifndef USE_SSE2
				LIS_QUAD_MULD(y[i],yl[i],x[i],xl[i],A->D->value[i]);
			#else
				LIS_QUAD_MULD_SSE2(y[i],yl[i],x[i],xl[i],A->D->value[i]);
			#endif

			tt.hi[0] = tt.hi[1] = tt.lo[0] = tt.lo[1] = 0.0;
			is = A->L->ptr[i];
			ie = A->L->ptr[i+1];
			for(j=is;j<ie-1;j+=2)
			{
				j0 = A->L->index[j+0];
				j1 = A->L->index[j+1];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD2_SSE2_LDSD(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],x[j0],xl[j0],x[j1],xl[j1],A->L->value[j]);
				#endif
			}
			for(;j<ie;j++)
			{
				j0 = A->L->index[j+0];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->L->value[j]);
				#endif
			}
			is = A->U->ptr[i];
			ie = A->U->ptr[i+1];
			for(j=is;j<ie-1;j+=2)
			{
				j0 = A->U->index[j+0];
				j1 = A->U->index[j+1];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD2_SSE2_LDSD(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],x[j0],xl[j0],x[j1],xl[j1],A->U->value[j]);
				#endif
			}
			for(;j<ie;j++)
			{
				j0 = A->U->index[j+0];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],A->U->value[j]);
				#endif
			}
			#ifdef USE_SSE2
				LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],tt.hi[0],tt.lo[0]);
				LIS_QUAD_ADD_SSE2(y[i],yl[i],y[i],yl[i],tt.hi[1],tt.lo[1]);
			#endif
		}
	}
	else
	{
		jj0 = A->index;
		vv0 = A->value;
		#ifdef _OPENMP
		#ifndef USE_SSE2
			#pragma omp parallel for private(i,j,is,ie,j0,j1,tt,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		#else
			#pragma omp parallel for private(i,j,is,ie,j0,j1,tt,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
		#endif
		#endif
		for(i=0;i<np;i++)
		{
			tt.hi[0] = tt.hi[1] = tt.lo[0] = tt.lo[1] = 0.0;

			is = A->ptr[i];
			ie = A->ptr[i+1];
			for(j=is;j<ie-1;j+=2)
			{
				j0 = jj0[j+0];
				j1 = jj0[j+1];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD2_SSE2_LDSD(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],x[j0],xl[j0],x[j1],xl[j1],vv0[j]);
				#endif
			}
			#ifdef USE_SSE2
				LIS_QUAD_ADD_SSE2(y[i],yl[i],tt.hi[0],tt.lo[0],tt.hi[1],tt.lo[1]);
			#endif
			for(;j<ie;j++)
			{
				j0 = jj0[j+0];
				#ifdef USE_SSE2
					LIS_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],vv0[j]);
				#endif
			}
		}
	}
}
#endif
