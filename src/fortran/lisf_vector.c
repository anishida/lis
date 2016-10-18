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
 * lis_vector_create_f
 * lis_vector_duplicate_f
 ************************************************/

#ifdef USE_FORTRAN

#undef __FUNC__
#define __FUNC__ "lis_vector_create_f"
void lis_vector_create_f(LIS_Comm_f *comm, LIS_VECTOR_F *vec, LIS_INT *ierr)
{
	LIS_VECTOR v;
	LIS_Comm c_comm;

	LIS_DEBUG_FUNC_IN;

	#ifdef USE_MPI
		if( *comm==lis_comm_world_f )
		{
			c_comm = MPI_COMM_WORLD;
		}
		else
		{
			c_comm = MPI_Comm_f2c(*comm);
		}
	#else
		c_comm = *comm;
	#endif
	*ierr = lis_vector_create(c_comm,&v);
	if( *ierr )	return;

	v->origin = LIS_ORIGIN_1;

	*vec = LIS_P2V(v);
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_size_f"
void lis_vector_set_size_f(LIS_VECTOR_F *vec, LIS_INT *local_n, LIS_INT *global_n, LIS_INT *ierr)
{
	LIS_VECTOR v;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_set_size((LIS_VECTOR)LIS_V2P(vec),*local_n,*global_n);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_vector_psd_reset_scale_f"
void lis_vector_psd_reset_scale_f(LIS_VECTOR_F *vec, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_psd_reset_scale((LIS_VECTOR)LIS_V2P(vec));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}


#undef __FUNC__
#define __FUNC__ "lis_vector_duplicate_f"
void lis_vector_duplicate_f(LIS_VECTOR_F *vin, LIS_VECTOR_F *vout, LIS_INT *ierr)
{
	LIS_VECTOR v;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_duplicate((LIS_VECTOR)LIS_V2P(vin),&v);
	if( *ierr )	return;
	*vout = LIS_P2V(v);

	LIS_DEBUG_FUNC_OUT;
	return;
}

/*
#undef __FUNC__
#define __FUNC__ "lis_vector_duplicateM_f"
void lis_vector_duplicateM_f(LIS_MATRIX_F *Amat, LIS_VECTOR_F *vout, LIS_INT *ierr)
{
	LIS_VECTOR v;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_duplicateM((LIS_MATRIX)LIS_V2P(Amat),&v);
	if( *ierr )	return;
	*vout = LIS_P2V(v);

	LIS_DEBUG_FUNC_OUT;
	return;
}
*/

#undef __FUNC__
#define __FUNC__ "lis_vector_destroy_f"
void lis_vector_destroy_f(LIS_VECTOR_F *vec, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_destroy((LIS_VECTOR)LIS_V2P(vec));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_range_f"
void lis_vector_get_range_f(LIS_VECTOR_F *vec, LIS_INT *is, LIS_INT *ie, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_get_range((LIS_VECTOR)LIS_V2P(vec),is,ie);
	if( *ierr )	return;
	(*is)++;
	(*ie)++;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_size_f"
void lis_vector_get_size_f(LIS_VECTOR_F *vec, LIS_INT *local_n, LIS_INT *global_n, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_get_size((LIS_VECTOR)LIS_V2P(vec),local_n,global_n);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_value_f"
void lis_vector_set_value_f(LIS_INT *flag, LIS_INT *i, LIS_SCALAR *value, LIS_VECTOR_F *v, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_set_value(*flag,*i,*value,((LIS_VECTOR)LIS_V2P(v)));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_values_f"
void lis_vector_set_values_f(LIS_INT *flag, LIS_INT *count, LIS_INT *index, LIS_SCALAR *values, LIS_VECTOR_F *v, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_set_values(*flag,*count,index,values,((LIS_VECTOR)LIS_V2P(v)));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_values2_f"
void lis_vector_set_values2_f(LIS_INT *flag, LIS_INT *start, LIS_INT *count, LIS_SCALAR *values, LIS_VECTOR_F *v, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_set_values2(*flag,*start,*count,values,((LIS_VECTOR)LIS_V2P(v)));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_value_f"
void lis_vector_get_value_f(LIS_VECTOR_F *v, LIS_INT *i, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_get_value((LIS_VECTOR)LIS_V2P(v),*i,value);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_values_f"
void lis_vector_get_values_f(LIS_VECTOR_F *v, LIS_INT *start, LIS_INT *count, LIS_SCALAR *values, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_get_values((LIS_VECTOR)LIS_V2P(v),*start,*count,values);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_scatter_f"
void lis_vector_scatter_f(LIS_SCALAR *values, LIS_VECTOR_F *v, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_scatter(values,(LIS_VECTOR)LIS_V2P(v));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_gather_f"
void lis_vector_gather_f(LIS_VECTOR_F *v, LIS_SCALAR *values, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_gather((LIS_VECTOR)LIS_V2P(v),values);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_print_f"
void lis_vector_print_f(LIS_VECTOR_F *v)
{
	LIS_DEBUG_FUNC_IN;

	lis_vector_print((LIS_VECTOR)LIS_V2P(v));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_swap_f"
void lis_vector_swap_f(LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_swap((LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_copy_f"
void lis_vector_copy_f(LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_copy((LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_axpy_f"
void lis_vector_axpy_f(LIS_SCALAR *alpha, LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_axpy(*alpha,(LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_xpay_f"
void lis_vector_xpay_f(LIS_VECTOR_F *x, LIS_SCALAR *alpha, LIS_VECTOR_F *y, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_xpay((LIS_VECTOR)LIS_V2P(x),*alpha,(LIS_VECTOR)LIS_V2P(y));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_axpyz_f"
void lis_vector_axpyz_f(LIS_SCALAR *alpha, LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_VECTOR_F *z, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_axpyz(*alpha,(LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y),(LIS_VECTOR)LIS_V2P(z));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_scale_f"
void lis_vector_scale_f(LIS_SCALAR *alpha, LIS_VECTOR_F *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_scale(*alpha,(LIS_VECTOR)LIS_V2P(x));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_pmul_f"
void lis_vector_pmul_f(LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_VECTOR_F *z, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_pmul((LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y),(LIS_VECTOR)LIS_V2P(z));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_pdiv_f"
void lis_vector_pdiv_f(LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_VECTOR_F *z, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_pdiv((LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y),(LIS_VECTOR)LIS_V2P(z));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_all_f"
void lis_vector_set_all_f(LIS_SCALAR *alpha, LIS_VECTOR_F *v, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_set_all(*alpha,((LIS_VECTOR)LIS_V2P(v)));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_abs_f"
void lis_vector_abs_f(LIS_VECTOR_F *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_abs((LIS_VECTOR)LIS_V2P(x));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_reciprocal_f"
void lis_vector_reciprocal_f(LIS_VECTOR_F *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_reciprocal((LIS_VECTOR)LIS_V2P(x));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_conjugate_f"
void lis_vector_conjugate_f(LIS_VECTOR_F *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_conjugate((LIS_VECTOR)LIS_V2P(x));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_shift_f"
void lis_vector_shift_f(LIS_SCALAR *alpha, LIS_VECTOR_F *x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_shift(*alpha, (LIS_VECTOR)LIS_V2P(x));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_dot_f"
void lis_vector_dot_f(LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_dot((LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y),value);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_nhdot_f"
void lis_vector_nhdot_f(LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_nhdot((LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y),value);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_nrm2_f"
void lis_vector_nrm2_f(LIS_VECTOR_F *x, LIS_REAL *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_nrm2((LIS_VECTOR)LIS_V2P(x),value);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_nrm1_f"
void lis_vector_nrm1_f(LIS_VECTOR_F *x, LIS_REAL *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_nrm1((LIS_VECTOR)LIS_V2P(x),value);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_nrmi_f"
void lis_vector_nrmi_f(LIS_VECTOR_F *x, LIS_REAL *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_nrmi((LIS_VECTOR)LIS_V2P(x),value);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_sum_f"
void lis_vector_sum_f(LIS_VECTOR_F *x, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_sum((LIS_VECTOR)LIS_V2P(x),value);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_is_null_f"
void lis_vector_is_null_f(LIS_VECTOR_F *vec, LIS_INT *ierr)
{
	LIS_VECTOR v;

	LIS_DEBUG_FUNC_IN;

	v = (LIS_VECTOR)LIS_V2P(vec);
	if( !lis_is_malloc(v) )
	{
		*ierr = LIS_TRUE;
	}
	else
	{
		if( v->status==LIS_VECTOR_NULL )
		{
			*ierr = LIS_TRUE;
		}
		else
		{
			*ierr = LIS_FALSE;
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return;
}

#endif
