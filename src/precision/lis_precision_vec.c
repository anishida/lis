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
#undef __FUNC__
#define __FUNC__ "lis_quad_malloc"
LIS_INT lis_quad_malloc(LIS_QUAD_PTR *a, LIS_INT n)
{
	double *ah, *al;
	
	LIS_DEBUG_FUNC_IN;

	ah = (double *)lis_malloc(2*n*sizeof(double),"lis_quad_malloc::ah");
	al = &ah[n];
	a->hi = ah;
	a->lo = al;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_quad_free"
LIS_INT lis_quad_free(LIS_QUAD_PTR *a)
{
	LIS_DEBUG_FUNC_IN;

	lis_free(a->hi);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_vector_axpyex_mmm"
LIS_INT lis_vector_axpyex_mmm(LIS_QUAD_PTR alpha, LIS_VECTOR vx, LIS_VECTOR vy)
{
	LIS_INT i,n,is,ie,nprocs,my_rank;
	LIS_QUAD_PTR aa;
	LIS_SCALAR *x,*xl,*y,*yl;
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	n   = vx->n;
	x   = vx->value;
	y   = vy->value;
	xl  = vx->value_lo;
	yl  = vy->value_lo;
	aa.hi = &vx->work[4];
	aa.lo = &vx->work[6];
	#ifndef USE_FMA2_SSE2
	    #pragma cdir nodep
	    #pragma _NEC ivdep
		#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		for(i=0; i<n; i++)
		{
			LIS_QUAD_FMA(y[i],yl[i],y[i],yl[i],alpha.hi[0],alpha.lo[0],x[i],xl[i]);
		}
	#else
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif
		aa.hi[0] = aa.hi[1] = alpha.hi[0];
		aa.lo[0] = aa.lo[1] = alpha.lo[0];

		#ifdef _OPENMP
		#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,is,ie,my_rank)
		#endif
		{
			#ifdef _OPENMP
				my_rank = omp_get_thread_num();
			#else
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			for(i=is;i<ie-1;i+=2)
			{
				LIS_QUAD_FMA2_SSE2(y[i],yl[i],y[i],yl[i],aa.hi[0],aa.lo[0],x[i],xl[i]);
			}
			for(;i<ie;i++)
			{
				LIS_QUAD_FMA_SSE2(y[i],yl[i],y[i],yl[i],alpha.hi[0],alpha.lo[0],x[i],xl[i]);
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_axpyzex_mmmm"
LIS_INT lis_vector_axpyzex_mmmm(LIS_QUAD_PTR alpha, LIS_VECTOR vx, LIS_VECTOR vy, LIS_VECTOR vz)
{
	LIS_INT i,n,is,ie,nprocs,my_rank;
	LIS_QUAD_PTR aa;
	LIS_SCALAR *x,*y,*z;
	LIS_SCALAR *xl,*yl,*zl;
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	n    = vx->n;
	x    = vx->value;
	y    = vy->value;
	z    = vz->value;
	xl   = vx->value_lo;
	yl   = vy->value_lo;
	zl   = vz->value_lo;
	aa.hi = &vx->work[4];
	aa.lo = &vx->work[6];
	#ifndef USE_FMA2_SSE2
	    #pragma cdir nodep
	    #pragma _NEC ivdep	
		#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		for(i=0; i<n; i++)
		{
			LIS_QUAD_FMA(z[i],zl[i],y[i],yl[i],alpha.hi[0],alpha.lo[0],x[i],xl[i]);
		}
	#else
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif
		aa.hi[0] = aa.hi[1] = alpha.hi[0];
		aa.lo[0] = aa.lo[1] = alpha.lo[0];

		#ifdef _OPENMP
		#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,is,ie,my_rank)
		#endif
		{
			#ifdef _OPENMP
				my_rank = omp_get_thread_num();
			#else
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			for(i=is;i<ie-1;i+=2)
			{
				LIS_QUAD_FMA2_SSE2(z[i],zl[i],y[i],yl[i],aa.hi[0],aa.lo[0],x[i],xl[i]);
			}
			for(;i<ie;i++)
			{
				LIS_QUAD_FMA_SSE2(z[i],zl[i],y[i],yl[i],alpha.hi[0],alpha.lo[0],x[i],xl[i]);
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_xpayex_mmm"
LIS_INT lis_vector_xpayex_mmm(LIS_VECTOR vx, LIS_QUAD_PTR alpha, LIS_VECTOR vy)
{
	LIS_INT i,n,is,ie,nprocs,my_rank;
	LIS_QUAD_PTR aa;
	LIS_SCALAR *x,*y,*xl,*yl;
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	n   = vx->n;
	x   = vx->value;
	y   = vy->value;
	xl  = vx->value_lo;
	yl  = vy->value_lo;
	aa.hi = &vx->work[4];
	aa.lo = &vx->work[6];
	#ifndef USE_FMA2_SSE2
	    #pragma cdir nodep
	    #pragma _NEC ivdep	
		#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
		for(i=0; i<n; i++)
		{
			LIS_QUAD_FMA(y[i],yl[i],x[i],xl[i],alpha.hi[0],alpha.lo[0],y[i],yl[i]);
		}
	#else
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif
		aa.hi[0] = aa.hi[1] = alpha.hi[0];
		aa.lo[0] = aa.lo[1] = alpha.lo[0];

		#ifdef _OPENMP
		#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,is,ie,my_rank)
		#endif
		{
			#ifdef _OPENMP
				my_rank = omp_get_thread_num();
			#else
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			for(i=is;i<ie-1;i+=2)
			{
				LIS_QUAD_FMA2_SSE2(y[i],yl[i],x[i],xl[i],aa.hi[0],aa.lo[0],y[i],yl[i]);
			}
			for(;i<ie;i++)
			{
				LIS_QUAD_FMA_SSE2(y[i],yl[i],x[i],xl[i],alpha.hi[0],alpha.lo[0],y[i],yl[i]);
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_dotex_mmm"
LIS_INT lis_vector_dotex_mmm(LIS_VECTOR vx, LIS_VECTOR vy, LIS_QUAD_PTR *val)
{
	LIS_INT i,n;
	LIS_SCALAR *x,*y,*xl,*yl;
	LIS_QUAD_PTR dotm2,dotm,tmpm;
	#ifdef _OPENMP
		LIS_INT is,ie,nprocs,my_rank;
		LIS_SCALAR *gt;
	#endif
	#ifdef USE_MPI
		MPI_Comm comm;
	#endif
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	n  = vx->n;
	x  = vx->value;
	y  = vy->value;
	xl = vx->value_lo;
	yl = vy->value_lo;
	dotm2.hi = &vx->work[0];
	dotm2.lo = &vx->work[2];
	dotm.hi = &vx->work[8];
	dotm.lo = &vx->work[9];
	tmpm.hi = &vx->work[10];
	tmpm.lo = &vx->work[11];
	#ifndef NO_ERROR_CHECK
		if( n!=vy->n )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"length of vector x and y is not equal\n");
			return LIS_ERR_ILL_ARG;
		}
	#endif

	#ifdef USE_MPI
		comm   = vx->comm;
	#endif
	#ifdef _OPENMP
		gt     = lis_vec_tmp;
		nprocs = omp_get_max_threads();
		#ifndef USE_SSE2
			#pragma omp parallel private(i,p1,p2,tq,bhi,blo,chi,clo,sh,th,sl,tl,eh,el,is,ie,my_rank)
		#else
			#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,is,ie,my_rank)
		#endif
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			#ifndef USE_FMA2_SSE2
				gt[my_rank*LIS_VEC_TMP_PADD] = gt[my_rank*LIS_VEC_TMP_PADD+1] = 0.0;
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=is;i<ie;i++)
				{
					LIS_QUAD_FMA(gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+1],gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+1],y[i],yl[i],x[i],xl[i]);
				}
			#else
				gt[my_rank*LIS_VEC_TMP_PADD  ] = gt[my_rank*LIS_VEC_TMP_PADD+1] = 0.0;
				gt[my_rank*LIS_VEC_TMP_PADD+2] = gt[my_rank*LIS_VEC_TMP_PADD+3] = 0.0;
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=is;i<ie-1;i+=2)
				{
					LIS_QUAD_FMA2_SSE2(gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+2],gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+2],y[i],yl[i],x[i],xl[i]);
				}
				LIS_QUAD_ADD_SSE2(gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+1],gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+2],gt[my_rank*LIS_VEC_TMP_PADD+1],gt[my_rank*LIS_VEC_TMP_PADD+3]);
				for(;i<ie;i++)
				{
					LIS_QUAD_FMA_SSE2(gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+1],gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+1],y[i],yl[i],x[i],xl[i]);
				}
			#endif
		}
		dotm.hi[0] = dotm.lo[0] = 0.0;
		for(i=0;i<nprocs;i++)
		{
			#ifndef USE_SSE2
				LIS_QUAD_ADD(dotm.hi[0],dotm.lo[0],dotm.hi[0],dotm.lo[0],gt[i*LIS_VEC_TMP_PADD],gt[i*LIS_VEC_TMP_PADD+1]);
			#else
				LIS_QUAD_ADD_SSE2(dotm.hi[0],dotm.lo[0],dotm.hi[0],dotm.lo[0],gt[i*LIS_VEC_TMP_PADD],gt[i*LIS_VEC_TMP_PADD+1]);
			#endif
		}
	#else
		#ifndef USE_FMA2_SSE2
			dotm.hi[0] = dotm.lo[0] = 0.0;
			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=0;i<n;i++)
			{
				LIS_QUAD_FMA(dotm.hi[0],dotm.lo[0],dotm.hi[0],dotm.lo[0],y[i],yl[i],x[i],xl[i]);
			}
		#else
			dotm2.hi[0] = dotm2.hi[1] = 0.0;
			dotm2.lo[0] = dotm2.lo[1] = 0.0;
			for(i=0;i<n-1;i+=2)
			{
				LIS_QUAD_FMA2_SSE2(dotm2.hi[0],dotm2.lo[0],dotm2.hi[0],dotm2.lo[0],y[i],yl[i],x[i],xl[i]);
			}
			LIS_QUAD_ADD_SSE2(dotm.hi[0],dotm.lo[0],dotm2.hi[0],dotm2.lo[0],dotm2.hi[1],dotm2.lo[1]);
			for(;i<n;i++)
			{
				LIS_QUAD_FMA_SSE2(dotm.hi[0],dotm.lo[0],dotm.hi[0],dotm.lo[0],y[i],yl[i],x[i],xl[i]);
			}
		#endif
	#endif
	#ifdef USE_MPI
		MPI_Allreduce(dotm.hi,tmpm.hi,1,LIS_MPI_MSCALAR,LIS_MPI_MSUM,comm);
		val->hi[0] = tmpm.hi[0];
		val->lo[0] = tmpm.lo[0];
	#else
		val->hi[0] = dotm.hi[0];
		val->lo[0] = dotm.lo[0];
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_nrm2ex_mm"
LIS_INT lis_vector_nrm2ex_mm(LIS_VECTOR vx, LIS_QUAD_PTR *val)
{
	LIS_INT i,n;
	LIS_SCALAR *x,*xl;
	LIS_QUAD_PTR dotm2,dotm,tmpm;
	#ifdef _OPENMP
		LIS_INT is,ie,nprocs,my_rank;
		LIS_SCALAR *gt;
	#endif
	#ifdef USE_MPI
		MPI_Comm comm;
	#endif
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	n  = vx->n;
	x  = vx->value;
	xl = vx->value_lo;
	dotm2.hi = &vx->work[0];
	dotm2.lo = &vx->work[2];
	dotm.hi = &vx->work[8];
	dotm.lo = &vx->work[9];
	tmpm.hi = &vx->work[10];
	tmpm.lo = &vx->work[11];
	#ifdef USE_MPI
		comm   = vx->comm;
	#endif
	#ifdef _OPENMP
		gt     = lis_vec_tmp;
		nprocs = omp_get_max_threads();
		#ifndef USE_SSE2
			#pragma omp parallel private(i,is,ie,my_rank,p1,p2,tq,bhi,blo,chi,clo,sh,eh,sl,el,th,tl)
		#else
			#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,is,ie,my_rank)
		#endif
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			#ifndef USE_FMA2_SSE2
				gt[my_rank*LIS_VEC_TMP_PADD] = gt[my_rank*LIS_VEC_TMP_PADD+1] = 0.0;
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=is;i<ie;i++)
				{
					LIS_QUAD_FSA(gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+1],gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+1],x[i],xl[i]);
				}
			#else
				gt[my_rank*LIS_VEC_TMP_PADD  ] = gt[my_rank*LIS_VEC_TMP_PADD+1] = 0.0;
				gt[my_rank*LIS_VEC_TMP_PADD+2] = gt[my_rank*LIS_VEC_TMP_PADD+3] = 0.0;
				#ifdef USE_VEC_COMP
				#pragma cdir nodep
				#pragma _NEC ivdep
				#endif
				for(i=is;i<ie-1;i+=2)
				{
					LIS_QUAD_FSA2_SSE2(gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+2],gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+2],x[i],xl[i]);
				}
				LIS_QUAD_ADD_SSE2(gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+1],gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+2],gt[my_rank*LIS_VEC_TMP_PADD+1],gt[my_rank*LIS_VEC_TMP_PADD+3]);
				for(;i<ie;i++)
				{
					LIS_QUAD_FSA_SSE2(gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+1],gt[my_rank*LIS_VEC_TMP_PADD],gt[my_rank*LIS_VEC_TMP_PADD+1],x[i],xl[i]);
				}
			#endif
		}
		dotm.hi[0] = dotm.lo[0] = 0.0;
		for(i=0;i<nprocs;i++)
		{
			#ifndef USE_SSE2
				LIS_QUAD_ADD(dotm.hi[0],dotm.lo[0],dotm.hi[0],dotm.lo[0],gt[i*LIS_VEC_TMP_PADD],gt[i*LIS_VEC_TMP_PADD+1]);
			#else
				LIS_QUAD_ADD_SSE2(dotm.hi[0],dotm.lo[0],dotm.hi[0],dotm.lo[0],gt[i*LIS_VEC_TMP_PADD],gt[i*LIS_VEC_TMP_PADD+1]);
			#endif
		}
	#else
		#ifndef USE_FMA2_SSE2
			dotm.hi[0] = dotm.lo[0] = 0.0;
			#pragma cdir nodep
			#pragma _NEC ivdep
			for(i=0;i<n;i++)
			{
				LIS_QUAD_FSA(dotm.hi[0],dotm.lo[0],dotm.hi[0],dotm.lo[0],x[i],xl[i]);
			}
		#else
			dotm2.hi[0] = dotm2.hi[1] = 0.0;
			dotm2.lo[0] = dotm2.lo[1] = 0.0;
			for(i=0;i<n-1;i+=2)
			{
				LIS_QUAD_FSA2_SSE2(dotm2.hi[0],dotm2.lo[0],dotm2.hi[0],dotm2.lo[0],x[i],xl[i]);
			}
			LIS_QUAD_ADD_SSE2(dotm.hi[0],dotm.lo[0],dotm2.hi[0],dotm2.lo[0],dotm2.hi[1],dotm2.lo[1]);
			for(;i<n;i++)
			{
				LIS_QUAD_FSA_SSE2(dotm.hi[0],dotm.lo[0],dotm.hi[0],dotm.lo[0],x[i],xl[i]);
			}
		#endif
	#endif
	#ifdef USE_MPI
		MPI_Allreduce(dotm.hi,tmpm.hi,1,LIS_MPI_MSCALAR,LIS_MPI_MSUM,comm);
		#ifndef USE_SSE2
			LIS_QUAD_SQRT(val->hi[0],val->lo[0],tmpm.hi[0],tmpm.lo[0]);
		#else
			LIS_QUAD_SQRT_SSE2(val->hi[0],val->lo[0],tmpm.hi[0],tmpm.lo[0]);
		#endif
	#else
		#ifndef USE_SSE2
			LIS_QUAD_SQRT(val->hi[0],val->lo[0],dotm.hi[0],dotm.lo[0]);
		#else
			LIS_QUAD_SQRT_SSE2(val->hi[0],val->lo[0],dotm.hi[0],dotm.lo[0]);
		#endif
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_copyex_mm"
LIS_INT lis_vector_copyex_mm(LIS_VECTOR vx, LIS_VECTOR vy)
{
	LIS_INT i,n;
	LIS_SCALAR *x,*xl,*y,*yl;

	LIS_DEBUG_FUNC_IN;

	n = vx->n;
	x   = vx->value;
	y   = vy->value;
	xl  = vx->value_lo;
	yl  = vy->value_lo;
	#ifdef USE_VEC_COMP
    #pragma cdir nodep
    #pragma _NEC ivdep	
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
		y[i]  = x[i];
		yl[i] = xl[i];
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_copyex_nm"
LIS_INT lis_vector_copyex_nm(LIS_VECTOR vx, LIS_VECTOR vy)
{
	LIS_INT i,n;
	LIS_SCALAR *x,*y,*yl;

	LIS_DEBUG_FUNC_IN;

	n = vx->n;
	x   = vx->value;
	y   = vy->value;
	yl  = vy->value_lo;
	#ifdef USE_VEC_COMP
    #pragma cdir nodep
    #pragma _NEC ivdep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
		y[i]  = x[i];
		yl[i] = 0.0;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_copyex_mn"
LIS_INT lis_vector_copyex_mn(LIS_VECTOR vx, LIS_VECTOR vy)
{
	LIS_INT i,n;
	LIS_SCALAR *x,*y;

	LIS_DEBUG_FUNC_IN;

	n = vx->n;
	x   = vx->value;
	y   = vy->value;
	#ifdef USE_VEC_COMP
    #pragma cdir nodep
    #pragma _NEC ivdep	
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
		y[i]  = x[i];
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_scaleex_nm"
LIS_INT lis_vector_scaleex_nm(LIS_SCALAR alpha, LIS_VECTOR vx)
{
	LIS_INT i,n,is,ie,nprocs,my_rank;
	LIS_SCALAR *aa;
	LIS_SCALAR *x,*xl;
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	n  = vx->n;
	x  = vx->value;
	xl = vx->value_lo;
	aa = vx->work;
	#ifndef USE_FMA2_SSE2
	    #pragma cdir nodep
	    #pragma _NEC ivdep	
		#ifndef USE_SSE2
			#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,eh,sl,el,th,tl)
		#else
			#pragma omp parallel for private(i,bh,ch,sh,th,bl,sl,tl,p1,p2,t0,t1,t2,is,ie,my_rank)
		#endif
		for(i=0; i<n; i++)
		{
			LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],alpha);
		}
	#else
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif
		aa[0] = aa[1] = alpha;
		#ifdef _OPENMP
		#ifndef USE_SSE2
			#pragma omp parallel private(i,is,ie,my_rank,p1,p2,tq,bhi,blo,chi,clo,sh,eh,sl,el,th,tl)
		#else
			#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,is,ie,my_rank)
		#endif
		#endif
		{
			#ifdef _OPENMP
				my_rank = omp_get_thread_num();
			#else
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			for(i=is;i<ie-1;i+=2)
			{
				LIS_QUAD_MULD2_SSE2(x[i],xl[i],x[i],xl[i],aa[0]);
			}
			for(;i<ie;i++)
			{
				LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],alpha);
			}
		}
	#endif
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_scaleex_mm"
LIS_INT lis_vector_scaleex_mm(LIS_QUAD_PTR alpha, LIS_VECTOR vx)
{
	LIS_INT i,n,is,ie,nprocs,my_rank;
	LIS_QUAD_PTR aa;
	LIS_SCALAR *x,*xl;
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	n  = vx->n;
	x  = vx->value;
	xl = vx->value_lo;
	aa.hi = &vx->work[0];
	aa.lo = &vx->work[2];
	#ifndef USE_FMA2_SSE2
	    #pragma cdir nodep
	    #pragma _NEC ivdep		    	
		#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,eh,sl,el,th,tl)
		for(i=0; i<n; i++)
		{
			LIS_QUAD_MUL(x[i],xl[i],x[i],xl[i],alpha.hi[0],alpha.lo[0]);
		}
	#else
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
		#else
			nprocs = 1;
		#endif
		aa.hi[0] = aa.hi[1] = alpha.hi[0];
		aa.lo[0] = aa.lo[1] = alpha.lo[0];
		#ifdef _OPENMP
		#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,is,ie,my_rank)
		#endif
		{
			#ifdef _OPENMP
				my_rank = omp_get_thread_num();
			#else
				my_rank = 0;
			#endif
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			for(i=is;i<ie-1;i+=2)
			{
				LIS_QUAD_MUL2_SSE2(x[i],xl[i],x[i],xl[i],aa.hi[0],aa.lo[0]);
			}
			for(;i<ie;i++)
			{
				LIS_QUAD_MUL_SSE2(x[i],xl[i],x[i],xl[i],aa.hi[0],aa.lo[0]);
			}
		}
	#endif
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_allex_nm"
LIS_INT lis_vector_set_allex_nm(LIS_SCALAR alpha, LIS_VECTOR vx)
{
	LIS_INT i,n;
	LIS_SCALAR *x,*xl;

	LIS_DEBUG_FUNC_IN;

	n   = vx->n;
	x   = vx->value;
	xl  = vx->value_lo;
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
		x[i]  = alpha;
		xl[i] = 0.0;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_reciprocalex_m"
LIS_INT lis_vector_reciprocalex_m(LIS_VECTOR vx)
{
	LIS_INT i,n;
	LIS_SCALAR *x,*xl;
	LIS_SCALAR one_hi,one_lo;
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	n   = vx->n;
	x   = vx->value;
	xl  = vx->value_lo;
	one_hi = 1.0;
	one_lo = 0.0;
	#ifdef USE_VEC_COMP
    #pragma cdir nodep
    #pragma _NEC ivdep		
	#endif
	#ifdef _OPENMP
	#ifndef USE_SSE2
		#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,eh,sl,el,th,tl)
	#else
		#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
	#endif
	#endif
	for(i=0; i<n; i++)
	{
		#ifndef USE_SSE2
			LIS_QUAD_DIV(x[i],xl[i],one_hi,one_lo,x[i],xl[i]);
		#else
			LIS_QUAD_DIV_SSE2(x[i],xl[i],one_hi,one_lo,x[i],xl[i]);
		#endif
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#ifdef USE_MPI
#undef __FUNC__
#define __FUNC__ "lis_mpi_msum"
void lis_mpi_msum(LIS_QUAD *invec, LIS_QUAD *inoutvec, LIS_INT *len, MPI_Datatype *datatype)
{
	LIS_INT	i;
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	for(i=0;i<*len;i++)
	{
		#ifndef USE_SSE2
			LIS_QUAD_ADD(inoutvec[i].hi,inoutvec[i].lo,inoutvec[i].hi,inoutvec[i].lo,invec[i].hi,invec[i].lo);
		#else
			LIS_QUAD_ADD_SSE2(inoutvec[i].hi,inoutvec[i].lo,inoutvec[i].hi,inoutvec[i].lo,invec[i].hi,invec[i].lo);
		#endif
	}

	LIS_DEBUG_FUNC_OUT;
}

#undef __FUNC__
#define __FUNC__ "lis_send_recv_mp"
LIS_INT lis_send_recv_mp(LIS_COMMTABLE commtable, LIS_VECTOR X)
{
	LIS_INT neib,i,is,inum,neibpetot,k,pad;
	LIS_SCALAR *ws,*wr;
	LIS_SCALAR *x,*xl;

	LIS_DEBUG_FUNC_IN;

	neibpetot = commtable->neibpetot;
	ws        = commtable->ws;
	wr        = commtable->wr;
	pad       = commtable->pad;
	x         = X->value;
	xl        = X->value_lo;

	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->export_ptr[neib];
		inum = commtable->export_ptr[neib+1] - is;
		for(i=0;i<inum;i++)
		{
			ws[is*2+i]      = x[commtable->export_index[is+i]];
			ws[is*2+inum+i] = xl[commtable->export_index[is+i]];
		}
		MPI_Isend(&ws[is*2],inum*2,LIS_MPI_SCALAR,commtable->neibpe[neib],0,commtable->comm,&commtable->req1[neib]);
	}
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->import_ptr[neib];
		inum = commtable->import_ptr[neib+1] - is;
		MPI_Irecv(&wr[is*2],inum*2,LIS_MPI_SCALAR,commtable->neibpe[neib],0,commtable->comm,&commtable->req2[neib]);
	}
	MPI_Waitall(neibpetot, commtable->req2, commtable->sta2);


	k  = commtable->import_index[0]+pad;
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->import_ptr[neib];
		inum = commtable->import_ptr[neib+1] - is;
		for(i=0;i<inum;i++)
		{
			x[k]    = wr[is*2+i];
			xl[k++] = wr[is*2+inum+i];
		}
	}

	MPI_Waitall(neibpetot, commtable->req1, commtable->sta1);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_reduce_mp"
LIS_INT lis_reduce_mp(LIS_COMMTABLE commtable, LIS_VECTOR X)
{
	LIS_INT neib,i,is,inum,neibpetot,pad;
	LIS_SCALAR *x,*xl;
	LIS_SCALAR *ws,*wr;
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	neibpetot = commtable->neibpetot;
	ws        = commtable->ws;
	wr        = commtable->wr;
	pad       = commtable->pad;
	x         = X->value;
	xl        = X->value_lo;

	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->import_ptr[neib];
		inum = commtable->import_ptr[neib+1] - is;
		for(i=0;i<inum;i++)
		{
			wr[is*2+i]      = x[commtable->import_index[is+i]+pad];
			wr[is*2+inum+i] = xl[commtable->import_index[is+i]+pad];
		}
		MPI_Isend(&wr[is*2],inum*2,LIS_MPI_SCALAR,commtable->neibpe[neib],0,commtable->comm,&commtable->req1[neib]);
	}
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->export_ptr[neib];
		inum = commtable->export_ptr[neib+1] - is;
		MPI_Irecv(&ws[is*2],inum*2,LIS_MPI_SCALAR,commtable->neibpe[neib],0,commtable->comm,&commtable->req2[neib]);
	}
	MPI_Waitall(neibpetot, commtable->req2, commtable->sta2);
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->export_ptr[neib];
		inum = commtable->export_ptr[neib+1] - is;
		for(i=0;i<inum;i++)
		{
			/*x[commtable->export_index[i]] += ws[i];*/
			#ifndef USE_SSE2
				LIS_QUAD_ADD(x[commtable->export_index[is+i]],xl[commtable->export_index[is+i]],x[commtable->export_index[is+i]],xl[commtable->export_index[is+i]],ws[is*2+i],ws[is*2+inum+i]);
			#else
				LIS_QUAD_ADD_SSE2(x[commtable->export_index[is+i]],xl[commtable->export_index[is+i]],x[commtable->export_index[is+i]],xl[commtable->export_index[is+i]],ws[is*2+i],ws[is*2+inum+i]);
			#endif
		}
	}
	MPI_Waitall(neibpetot, commtable->req1, commtable->sta1);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#endif

#endif
