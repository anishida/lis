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
#include <string.h>
#include <math.h>
#ifdef USE_QUAD_PRECISION
#ifdef _WIN32
	#include <float.h>
#endif
#endif
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
#define __FUNC__ "lis_quad_x87_fpu_init"
void lis_quad_x87_fpu_init(LIS_UNSIGNED_INT *cw_old)
{
#ifdef HAS_X87_FPU
#ifdef _WIN32
#ifndef _WIN64
	LIS_UNSIGNED_INT cw = _control87(0, 0);
	_control87(0x00010000, 0x00030000);
	*cw_old = cw;
	cw = _control87(0, 0);
#endif
#else
	LIS_INT cw,cw_new;
	asm volatile ("fnstcw %0":"=m" (cw));
	cw_new = (cw & ~0x0300) | 0x0200;
	asm volatile ("fldcw %0": :"m" (cw_new));
	*cw_old = cw;
#endif	
#else
	*cw_old = 0;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_quad_x87_fpu_finalize"
void lis_quad_x87_fpu_finalize(LIS_UNSIGNED_INT cw)
{
#ifdef HAS_X87_FPU
#ifdef _WIN32
#ifndef _WIN64  
    _control87(cw, 0xFFFFFFFF);
#endif
#else
	asm volatile ("fldcw %0": :"m" (cw));
#endif	
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_quad_minus"
void  lis_quad_minus(LIS_QUAD *a)
{
	a->hi = -a->hi;
	a->lo = -a->lo;
}

#undef __FUNC__
#define __FUNC__ "lis_quad_zero"
void  lis_quad_zero(LIS_QUAD *a)
{
	a->hi = 0.0;
	a->lo = 0.0;
}

#undef __FUNC__
#define __FUNC__ "lis_quad_one"
void  lis_quad_one(LIS_QUAD *a)
{
	a->hi = 1.0;
	a->lo = 0.0;
}

#undef __FUNC__
#define __FUNC__ "lis_quad_max"
void  lis_quad_max(LIS_QUAD *a, LIS_QUAD *b, LIS_QUAD *c)
{
	if( (a->hi > b->hi) || ((a->hi==b->hi) && (a->lo > b->lo)) )
	{
		c->hi = a->hi;
		c->lo = a->lo;
	}
	else
	{
		c->hi = b->hi;
		c->lo = b->lo;
	}
}

#undef __FUNC__
#define __FUNC__ "lis_quad_min"
void  lis_quad_min(LIS_QUAD *a, LIS_QUAD *b, LIS_QUAD *c)
{
	if( (a->hi < b->hi) || ((a->hi==b->hi) && (a->lo < b->lo)) )
	{
		c->hi = a->hi;
		c->lo = a->lo;
	}
	else
	{
		c->hi = b->hi;
		c->lo = b->lo;
	}
}



#undef __FUNC__
#define __FUNC__ "lis_quad_add"
void lis_quad_add(LIS_QUAD *a, const LIS_QUAD *b, const LIS_QUAD *c)
{
	LIS_QUAD_DECLAR;

	#ifndef USE_SSE2
		LIS_QUAD_ADD(a->hi,a->lo,b->hi,b->lo,c->hi,c->lo);
	#else
		LIS_QUAD_ADD_SSE2(a->hi,a->lo,b->hi,b->lo,c->hi,c->lo);
	#endif
}

#undef __FUNC__
#define __FUNC__ "lis_quad_sub"
void lis_quad_sub(LIS_QUAD *a, const LIS_QUAD *b, const LIS_QUAD *c)
{
	LIS_QUAD_DECLAR;

	#ifndef USE_SSE2
		LIS_QUAD_ADD(a->hi,a->lo,b->hi,b->lo,-c->hi,-c->lo);
	#else
		LIS_QUAD_ADD_SSE2(a->hi,a->lo,b->hi,b->lo,-c->hi,-c->lo);
	#endif
}


#undef __FUNC__
#define __FUNC__ "lis_quad_mul"
void lis_quad_mul(LIS_QUAD *a, const LIS_QUAD *b, const LIS_QUAD *c)
{
	LIS_QUAD_DECLAR;

	#ifndef USE_SSE2
		LIS_QUAD_MUL(a->hi,a->lo,b->hi,b->lo,c->hi,c->lo);
	#else
		LIS_QUAD_MUL_SSE2(a->hi,a->lo,b->hi,b->lo,c->hi,c->lo);
	#endif
}

#undef __FUNC__
#define __FUNC__ "lis_quad_mul_dd_d"
 void lis_quad_mul_dd_d(LIS_QUAD *a, const LIS_QUAD *b, const double c)
{
	LIS_QUAD_DECLAR;

	#ifndef USE_SSE2
		LIS_QUAD_MULD(a->hi,a->lo,b->hi,b->lo,c);
	#else
		LIS_QUAD_MULD_SSE2(a->hi,a->lo,b->hi,b->lo,c);
	#endif
}

#undef __FUNC__
#define __FUNC__ "lis_quad_sqr"
void lis_quad_sqr(LIS_QUAD *a, const LIS_QUAD *b)
{
	LIS_QUAD_DECLAR;

	#ifndef USE_SSE2
		LIS_QUAD_SQR(a->hi,a->lo,b->hi,b->lo);
	#else
		LIS_QUAD_SQR_SSE2(a->hi,a->lo,b->hi,b->lo);
	#endif
}

#undef __FUNC__
#define __FUNC__ "lis_quad_div"
void lis_quad_div(LIS_QUAD *a, const LIS_QUAD *b, const LIS_QUAD *c)
{
	LIS_QUAD_DECLAR;

	#ifndef USE_SSE2
		LIS_QUAD_DIV(a->hi,a->lo,b->hi,b->lo,c->hi,c->lo);
	#else
		LIS_QUAD_DIV_SSE2(a->hi,a->lo,b->hi,b->lo,c->hi,c->lo);
	#endif
}

#undef __FUNC__
#define __FUNC__ "lis_quad_sqrt"
LIS_INT lis_quad_sqrt(LIS_QUAD *a, const LIS_QUAD *b)
{
	LIS_QUAD_DECLAR;

	#ifndef USE_SSE2
		LIS_QUAD_SQRT(a->hi,a->lo,b->hi,b->lo);
	#else
		LIS_QUAD_SQRT_SSE2(a->hi,a->lo,b->hi,b->lo);
	#endif
	return LIS_SUCCESS;
}
#endif

