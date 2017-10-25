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
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_array_swap_f
 * lis_array_copy_f
 * lis_array_axpy_f
 * lis_array_xpay_f
 * lis_array_axpyz_f
 * lis_array_scale_f
 * lis_array_pmul_f
 * lis_array_pdiv_f
 * lis_array_set_all_f
 * lis_array_abs_f
 * lis_array_reciprocal_f
 * lis_array_conjugate_f
 * lis_array_shift_f
 * lis_array_dot_f
 * lis_array_nhdot_f
 * lis_array_nrm1_f
 * lis_array_nrm2_f
 * lis_array_nrmi_f
 * lis_array_sum_f
 * lis_array_matvec_f
 * lis_array_matvech_f
 * lis_array_matvec_ns_f
 * lis_array_matmat_f
 * lis_array_matmat_ns_f
 * lis_array_ge_f
 * lis_array_solve_f
 * lis_array_cgs_f
 * lis_array_mgs_f
 * lis_array_qr_f
 ************************************************/

#ifdef USE_FORTRAN

#undef __FUNC__
#define __FUNC__ "lis_array_swap_f"
void lis_array_swap_f(LIS_INT *n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_swap(*n,x,y);
	if( *ierr )	return;
	
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_copy_f"
void lis_array_copy_f(LIS_INT *n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_copy(*n,x,y);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_axpy_f"
void lis_array_axpy_f(LIS_INT *n, LIS_SCALAR *alpha, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_axpy(*n,*alpha,x,y);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_xpay_f"
void lis_array_xpay_f(LIS_INT *n, LIS_SCALAR *x, LIS_SCALAR *alpha, LIS_SCALAR *y, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_xpay(*n,x,*alpha,y);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_axpyz_f"
void lis_array_axpyz_f(LIS_INT *n, LIS_SCALAR *alpha, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_axpyz(*n,*alpha,x,y,z);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_scale_f"
void lis_array_scale_f(LIS_INT *n, LIS_SCALAR *alpha, LIS_SCALAR *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_scale(*n,*alpha,x);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_pmul_f"
void lis_array_pmul_f(LIS_INT *n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_pmul(*n,x,y,z);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_pdiv_f"
void lis_array_pdiv_f(LIS_INT *n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_pdiv(*n,x,y,z);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_set_all_f"
void lis_array_set_all_f(LIS_INT *n, LIS_SCALAR *alpha, LIS_SCALAR *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_set_all(*n,*alpha,x);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_abs_f"
void lis_array_abs_f(LIS_INT *n, LIS_SCALAR *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_abs(*n,x);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_reciprocal_f"
void lis_array_reciprocal_f(LIS_INT *n, LIS_SCALAR *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_reciprocal(*n,x);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_conjugate_f"
void lis_array_conjugate_f(LIS_INT *n, LIS_SCALAR *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_conjugate(*n,x);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_shift_f"
void lis_array_shift_f(LIS_INT *n, LIS_SCALAR *t, LIS_SCALAR *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_shift(*n,*t,x);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_dot_f"
void lis_array_dot_f(LIS_INT *n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_dot(*n,x,y,value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_nhdot_f"
void lis_array_nhdot_f(LIS_INT *n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_nhdot(*n,x,y,value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_nrm1_f"
void lis_array_nrm1_f(LIS_INT *n, LIS_SCALAR *x, LIS_REAL *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_nrm1(*n,x,value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_nrm2_f"
void lis_array_nrm2_f(LIS_INT *n, LIS_SCALAR *x, LIS_REAL *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_nrm2(*n,x,value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_nrmi_f"
void lis_array_nrmi_f(LIS_INT *n, LIS_SCALAR *x, LIS_REAL *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_nrmi(*n,x,value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_sum_f"
void lis_array_sum_f(LIS_INT *n, LIS_SCALAR *x, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_sum(*n,x,value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_matvec_f"
void lis_array_matvec_f(LIS_INT *n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT *op, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_matvec(*n,a,x,y,*op);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_matvech_f"
void lis_array_matvech_f(LIS_INT *n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT *op, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_matvech(*n,a,x,y,*op);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_matvec_ns_f"
 void lis_array_matvec_ns_f(LIS_INT *m, LIS_INT *n, LIS_SCALAR *a, LIS_INT *lda, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT *op, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_matvec_ns(*m,*n,a,*lda,x,y,*op);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_matmat_f"
void lis_array_matmat_f(LIS_INT *n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c, LIS_INT *op, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_matmat(*n,a,b,c,*op);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_matmat_ns_f"
 void lis_array_matmat_ns_f(LIS_INT *m, LIS_INT *n, LIS_INT *k, LIS_SCALAR *a, LIS_INT *lda, LIS_SCALAR *b, LIS_INT *ldb, LIS_SCALAR *c, LIS_INT *ldc, LIS_INT *op, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_matmat_ns(*m,*n,*k,a,*lda,b,*ldb,c,*ldc,*op);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_ge_f"
void lis_array_ge_f(LIS_INT *n, LIS_SCALAR *a, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_ge(*n,a);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_solve_f"
 void lis_array_solve_f(LIS_INT *n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *x, LIS_SCALAR *w, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_solve(*n,a,b,x,w);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_cgs_f"
 void lis_array_cgs_f(LIS_INT *n, LIS_SCALAR *a, LIS_SCALAR *q, LIS_SCALAR *r, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_cgs(*n,a,q,r);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_mgs_f"
 void lis_array_mgs_f(LIS_INT *n, LIS_SCALAR *a, LIS_SCALAR *q, LIS_SCALAR *r, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_mgs(*n,a,q,r);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_array_qr_f"
 void lis_array_qr_f(LIS_INT *n, LIS_SCALAR *a, LIS_SCALAR *q, LIS_SCALAR *r, LIS_INT *qriter, LIS_REAL *qrerr, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_array_qr(*n,a,q,r,qriter,qrerr);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#endif


