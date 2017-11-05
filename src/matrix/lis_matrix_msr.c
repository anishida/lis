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
 * function                    | SOM |
 *-----------------------------+-----+
 * lis_matrix_set              | o   |
 * lis_matrix_setDLU           | o   |
 * lis_matrix_malloc           | o   |
 * lis_matrix_elements_copy    | o   |
 * lis_matrix_transpose        | o   |
 * lis_matrix_split            | o   |
 * lis_matrix_merge            | o   |
 *-----------------------------+-----+-----+
 * function                    |merge|split|
 *-----------------------------+-----+-----|
 * lis_matrix_convert          | o   |     |
 * lis_matrix_copy             | o   | o   |
 * lis_matrix_get_diagonal     | o   | o   |
 * lis_matrix_shift_diagonal   | o   | o   |
 * lis_matrix_scale            | o   | o   |
 * lis_matrix_scale_symm       | o   | o   |
 * lis_matrix_normf            | o   | o   |
 * lis_matrix_sort             | o   | o   |
 * lis_matrix_solve            | xxx | o   |
 * lis_matrix_solveh           | xxx | o   |
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_msr"
LIS_INT lis_matrix_set_msr(LIS_INT nnz, LIS_INT ndz, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

#if 0
	err = lis_matrix_check(A,LIS_MATRIX_CHECK_SET);
	if( err ) return err;
#else
	if(lis_matrix_is_assembled(A)) {
	  LIS_DEBUG_FUNC_OUT;
	  return LIS_SUCCESS;
	}
	else {
	  err = lis_matrix_check(A,LIS_MATRIX_CHECK_SET);
	  if( err ) return err;
	}
#endif

	A->index       = index;
	A->value       = value;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_MSR;
	A->nnz         = nnz;
	A->ndz         = ndz;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_setDLU_msr"
LIS_INT lis_matrix_setDLU_msr(LIS_INT lnnz, LIS_INT unnz, LIS_INT lndz, LIS_INT undz, LIS_SCALAR *diag, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A)
{
	LIS_INT err;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;

#if 0
	err = lis_matrix_check(A,LIS_MATRIX_CHECK_SET);
	if( err ) return err;
#else
	if(lis_matrix_is_assembled(A))  return LIS_SUCCESS;
	else {
	  err = lis_matrix_check(A,LIS_MATRIX_CHECK_SET);
	  if( err ) return err;
	}
#endif

	A->L = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_msr::A->L");
	if( A->L==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	A->U = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_msr::A->U");
	if( A->U==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		lis_matrix_DLU_destroy(A);
		return LIS_OUT_OF_MEMORY;
	}
	err = lis_matrix_diag_create(A->n,0,A->comm,&D);
	if( err )
	{
		lis_matrix_DLU_destroy(A);
		return err;
	}

	lis_free(D->value);
	D->value       = diag;
	A->D           = D;
	A->L->nnz      = lnnz;
	A->L->ndz      = lndz;
	A->L->index    = lindex;
	A->L->value    = lvalue;
	A->U->nnz      = unnz;
	A->U->ndz      = undz;
	A->U->index    = uindex;
	A->U->value    = uvalue;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_MSR;
	A->is_splited  = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_msr"
LIS_INT lis_matrix_malloc_msr(LIS_INT n, LIS_INT nnz, LIS_INT ndz, LIS_INT **index, LIS_SCALAR **value)
{
	LIS_DEBUG_FUNC_IN;

	*index   = NULL;
	*value   = NULL;

	*index = (LIS_INT *)lis_malloc( (nnz+ndz+1)*sizeof(LIS_INT),"lis_matrix_malloc_msr::index" );
	if( *index==NULL )
	{
		LIS_SETERR_MEM((nnz+ndz+1)*sizeof(LIS_INT));
		lis_free2(2,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*value = (LIS_SCALAR *)lis_malloc( (nnz+ndz+1)*sizeof(LIS_SCALAR),"lis_matrix_malloc_msr::value" );
	if( *value==NULL )
	{
		LIS_SETERR_MEM((nnz+ndz+1)*sizeof(LIS_SCALAR));
		lis_free2(2,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_elements_copy_msr"
LIS_INT lis_matrix_elements_copy_msr(LIS_INT n, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_index, LIS_SCALAR *o_value)
{
	LIS_INT i,j;

	LIS_DEBUG_FUNC_IN;

	#ifdef _OPENMP
	#pragma omp parallel private(i,j)
	#endif
	{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n+1;i++)
		{
			o_index[i] = index[i];
			o_value[i] = value[i];
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n;i++)
		{
			for(j=index[i];j<index[i+1];j++)
			{
				o_value[j]   = value[j];
				o_index[j]   = index[j];
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_msr"
LIS_INT lis_matrix_copy_msr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT i,n,nnz,ndz,lnnz,unnz,lndz,undz;
	LIS_INT *index;
	LIS_INT *lindex;
	LIS_INT *uindex;
	LIS_SCALAR *value,*lvalue,*uvalue,*diag;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;

	if( Ain->is_splited )
	{
		lnnz     = Ain->L->nnz;
		unnz     = Ain->U->nnz;
		lndz     = Ain->L->ndz;
		undz     = Ain->U->ndz;
		lindex   = NULL;
		uindex   = NULL;
		diag     = NULL;

		err = lis_matrix_malloc_msr(n,lnnz,lndz,&lindex,&lvalue);
		if( err )
		{
			return err;
		}
		err = lis_matrix_malloc_msr(n,unnz,undz,&uindex,&uvalue);
		if( err )
		{
			lis_free2(5,diag,uindex,lindex,uvalue,lvalue);
			return err;
		}
		diag = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_matrix_copy_msr::diag");
		if( diag==NULL )
		{
			lis_free2(5,diag,uindex,lindex,uvalue,lvalue);
			return err;
		}

		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<n;i++)
		{
			diag[i] = Ain->D->value[i];
		}
		lis_matrix_elements_copy_msr(n,Ain->L->index,Ain->L->value,lindex,lvalue);
		lis_matrix_elements_copy_msr(n,Ain->U->index,Ain->U->value,uindex,uvalue);

		err = lis_matrix_setDLU_msr(lnnz,unnz,lndz,undz,diag,lindex,lvalue,uindex,uvalue,Aout);
		if( err )
		{
			lis_free2(5,diag,uindex,lindex,uvalue,lvalue);
			return err;
		}
	}
	if( !Ain->is_splited || (Ain->is_splited && Ain->is_save) )
	{
		index   = NULL;
		value   = NULL;
		nnz     = Ain->nnz;
		ndz     = Ain->ndz;
		err = lis_matrix_malloc_msr(n,nnz,ndz,&index,&value);
		if( err )
		{
			return err;
		}

		lis_matrix_elements_copy_msr(n,Ain->index,Ain->value,index,value);

		err = lis_matrix_set_msr(nnz,ndz,index,value,Aout);
		if( err )
		{
			lis_free2(2,index,value);
			return err;
		}
	}

	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal_msr"
LIS_INT lis_matrix_get_diagonal_msr(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i;
	LIS_INT n;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0; i<n; i++)
		{
			d[i] = A->D->value[i];
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0; i<n; i++)
		{
			d[i] = A->value[i];
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_msr"
LIS_INT lis_matrix_shift_diagonal_msr(LIS_MATRIX A, LIS_SCALAR sigma)
{
	LIS_INT i;
	LIS_INT n;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0; i<n; i++)
		{
			A->D->value[i] -= sigma;
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0; i<n; i++)
		{
			A->value[i] -= sigma;
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_msr"
LIS_INT lis_matrix_scale_msr(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j;
	LIS_INT n;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			A->D->value[i] = 1.0;
			for(j=A->L->index[i];j<A->L->index[i+1];j++)
			{
				A->L->value[j] *= d[i];
			}
			for(j=A->U->index[i];j<A->U->index[i+1];j++)
			{
				A->U->value[j] *= d[i];
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			A->value[i] = 1.0;
			for(j=A->index[i];j<A->index[i+1];j++)
			{
				A->value[j] *= d[i];
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_symm_msr"
LIS_INT lis_matrix_scale_symm_msr(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j;
	LIS_INT n;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			A->D->value[i] = 1.0;
			for(j=A->L->index[i];j<A->L->index[i+1];j++)
			{
				A->L->value[j] = A->L->value[j]*d[i]*d[A->L->index[j]];
			}
			for(j=A->U->index[i];j<A->U->index[i+1];j++)
			{
				A->U->value[j] = A->U->value[j]*d[i]*d[A->U->index[j]];
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			A->value[i] = 1.0;
			for(j=A->index[i];j<A->index[i+1];j++)
			{
				A->value[j] = A->value[j]*d[i]*d[A->index[j]];
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_normf_msr"
LIS_INT lis_matrix_normf_msr(LIS_MATRIX A, LIS_SCALAR *nrm)
{
	LIS_INT i,j;
	LIS_INT n;
	LIS_SCALAR sum;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	sum  = (LIS_SCALAR)0;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for reduction(+:sum) private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			sum += A->D->value[i]*A->D->value[i];
			for(j=A->L->index[i];j<A->L->index[i+1];j++)
			{
				sum += A->L->value[j]*A->L->value[j];
			}
			for(j=A->U->index[i];j<A->U->index[i+1];j++)
			{
				sum += A->U->value[j]*A->U->value[j];
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for reduction(+:sum) private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			sum += A->value[i]*A->value[i];
			for(j=A->index[i];j<A->index[i+1];j++)
			{
				sum += A->value[j]*A->value[j];
			}
		}
	}
	*nrm = sqrt(sum);
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_transpose_msr"
LIS_INT lis_matrix_transpose_msr(LIS_MATRIX Ain, LIS_MATRIX *Aout)
{

	LIS_DEBUG_FUNC_IN;

/*	err = lis_matrix_convert_msr2csc(Ain,Aout);*/
	(*Aout)->matrix_type = LIS_MATRIX_MSR;
	(*Aout)->status      = LIS_MATRIX_MSR;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_msr"
LIS_INT lis_matrix_split_msr(LIS_MATRIX A)
{
	LIS_INT i,j,n;
	LIS_INT lnnz,unnz;
	LIS_INT lndz,undz;
	LIS_INT err;
	LIS_INT *lindex,*uindex;
	LIS_SCALAR *lvalue,*uvalue;
	#ifdef _OPENMP
		LIS_INT kl,ku;
		LIS_INT *liw,*uiw;
	#endif
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;

	n        = A->n;
	lnnz     = 0;
	unnz     = 0;
	lndz     = n;
	undz     = n;
	D        = NULL;
	lindex   = NULL;
	lvalue   = NULL;
	uindex   = NULL;
	uvalue   = NULL;

	#ifdef _OPENMP
		liw = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_split_msr::liw");
		if( liw==NULL )
		{
			LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		uiw = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_split_msr::uiw");
		if( uiw==NULL )
		{
			LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
			lis_free(liw);
			return LIS_OUT_OF_MEMORY;
		}
		#pragma omp parallel for private(i)
		for(i=0;i<n+1;i++)
		{
			liw[i] = 0;
			uiw[i] = 0;
		}
		#pragma omp parallel for private(i,j)
		for(i=0;i<n;i++)
		{
			for(j=A->index[i];j<A->index[i+1];j++)
			{
				if( A->index[j]<i )
				{
					liw[i+1]++;
				}
				else if( A->index[j]>i )
				{
					uiw[i+1]++;
				}
			}
		}
		liw[0] = n+1;
		uiw[0] = n+1;
		for(i=0;i<n;i++)
		{
			liw[i+1] += liw[i];
			uiw[i+1] += uiw[i];
		}
		lnnz = liw[n];
		unnz = uiw[n];
	#else
		for(i=0;i<n;i++)
		{
			for(j=A->index[i];j<A->index[i+1];j++)
			{
				if( A->index[j]<i )
				{
					lnnz++;
				}
				else if( A->index[j]>i )
				{
					unnz++;
				}
			}
		}
	#endif

	err = lis_matrix_LU_create(A);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_msr(n,lnnz,lndz,&lindex,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_msr(n,unnz,undz,&uindex,&uvalue);
	if( err )
	{
		lis_free2(4,lindex,lvalue,uindex,uvalue);
		return err;
	}
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		lis_free2(4,lindex,lvalue,uindex,uvalue);
		return err;
	}

	#ifdef _OPENMP
		#pragma omp parallel for private(i)
		for(i=0;i<n+1;i++)
		{
			lindex[i] = liw[i];
			uindex[i] = uiw[i];
		}
		#pragma omp parallel for private(i,j,kl,ku)
		for(i=0;i<n;i++)
		{
			kl = lindex[i];
			ku = uindex[i];
			D->value[i] = A->value[i];
			for(j=A->index[i];j<A->index[i+1];j++)
			{
				if( A->index[j]<i )
				{
					lindex[kl]   = A->index[j];
					lvalue[kl]   = A->value[j];
					kl++;
				}
				else if( A->index[j]>i )
				{
					uindex[ku]   = A->index[j];
					uvalue[ku]   = A->value[j];
					ku++;
				}
			}
		}
		lis_free2(2,liw,uiw);
	#else
		lnnz = n+1;
		unnz = n+1;
		lindex[0] = n+1;
		uindex[0] = n+1;
		for(i=0;i<n;i++)
		{
			D->value[i] = A->value[i];
			for(j=A->index[i];j<A->index[i+1];j++)
			{
				if( A->index[j]<i )
				{
					lindex[lnnz]   = A->index[j];
					lvalue[lnnz]   = A->value[j];
					lnnz++;
				}
				else if( A->index[j]>i )
				{
					uindex[unnz]   = A->index[j];
					uvalue[unnz]   = A->value[j];
					unnz++;
				}
			}
			lindex[i+1] = lnnz;
			uindex[i+1] = unnz;
		}
	#endif
	A->L->nnz     = lnnz - (n+1);
	A->L->ndz     = lndz;
	A->L->index   = lindex;
	A->L->value   = lvalue;
	A->U->nnz     = unnz - (n+1);
	A->U->ndz     = undz;
	A->U->index   = uindex;
	A->U->value   = uvalue;
	A->D          = D;
	A->is_splited = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_merge_msr"
LIS_INT lis_matrix_merge_msr(LIS_MATRIX A)
{
	LIS_INT i,j,n,is;
	LIS_INT nnz,ndz;
	LIS_INT err;
	LIS_INT *index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;


	n       = A->n;
	nnz     = 0;
	ndz     = 0;
	is      = A->is;
	index   = NULL;
	value   = NULL;
	nnz     = A->L->nnz + A->U->nnz + n;

	err = lis_matrix_malloc_msr(n,nnz,ndz,&index,&value);
	if( err )
	{
		return err;
	}

	nnz      = n+1;
	index[0] = n+1;
	if( A->matrix_type==LIS_MATRIX_MSR )
	{
		for(i=0;i<n;i++)
		{
			value[i]   = A->D->value[i];
			for(j=A->L->index[i];j<A->L->index[i+1];j++)
			{
				index[nnz]   = A->L->index[j];
				value[nnz]   = A->L->value[j];
				nnz++;
			}
			for(j=A->U->index[i];j<A->U->index[i+1];j++)
			{
				index[nnz]   = A->U->index[j];
				value[nnz]   = A->U->value[j];
				nnz++;
			}
			index[i+1] = nnz;
		}
	}
	else
	{
		for(i=0;i<n;i++)
		{
			value[i]   = A->D->value[i];
			for(j=A->U->index[i];j<A->U->index[i+1];j++)
			{
				index[nnz]   = A->U->index[j];
				value[nnz]   = A->U->value[j];
				nnz++;
			}
			for(j=A->L->index[i];j<A->L->index[i+1];j++)
			{
				index[nnz]   = A->L->index[j];
				value[nnz]   = A->L->value[j];
				nnz++;
			}
			index[i+1] = nnz;
		}
	}

	A->nnz        = nnz;
	A->ndz        = ndz;
	A->value      = value;
	A->index      = index;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_sort_msr"
LIS_INT lis_matrix_sort_msr(LIS_MATRIX A)
{
	LIS_INT i,n;

	LIS_DEBUG_FUNC_IN;

	if( !A->is_sorted )
	{
		n = A->n;
		if( A->is_splited )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<n;i++)
			{
				lis_sort_id(A->L->ptr[i],A->L->ptr[i+1]-1,A->L->index,A->L->value);
				lis_sort_id(A->U->ptr[i],A->U->ptr[i+1]-1,A->U->index,A->U->value);
			}
		}
		else
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<n;i++)
			{
				lis_sort_id(A->ptr[i],A->ptr[i+1]-1,A->index,A->value);
			}
		}
		A->is_sorted = LIS_TRUE;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solve_msr"
LIS_INT lis_matrix_solve_msr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,n;
	LIS_SCALAR t;
	LIS_SCALAR *b,*x;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	b  = B->value;
	x  = X->value;

	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=0;i<n;i++)
		{
			t = b[i];
			for(j=A->L->index[i];j<A->L->index[i+1];j++)
			{
				t -= A->L->value[j] * x[A->L->index[j]];
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			t = b[i];
			for(j=A->U->index[i];j<A->U->index[i+1];j++)
			{
				t -= A->U->value[j] * x[A->U->index[j]];
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<n;i++)
		{
			t = b[i];
			for(j=A->L->index[i];j<A->L->index[i+1];j++)
			{

				t -= A->L->value[j] * x[A->L->index[j]];
			}
			x[i]   = t * A->WD->value[i];
		}
		for(i=n-1;i>=0;i--)
		{
			t = 0.0;
			for(j=A->U->index[i];j<A->U->index[i+1];j++)
			{
				if( A->U->index[j]>=n ) continue; 
				t += A->U->value[j] * x[A->U->index[j]];
			}
			x[i]  -= t * A->WD->value[i];
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solveh_msr"
LIS_INT lis_matrix_solveh_msr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,n;
	LIS_SCALAR t;
	LIS_SCALAR *x;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	x  = X->value;

	lis_vector_copy(B,X);
	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=n-1;i>=0;i--)
		{
			x[i]   = x[i] * conj(A->WD->value[i]);
			for(j=A->U->index[i];j<A->U->index[i+1];j++)
			{
				x[A->U->index[j]] -= conj(A->U->value[j]) * x[i];
			}
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=0;i<n;i++)
		{
			x[i]   = x[i] * conj(A->WD->value[i]);
			for(j=A->L->index[i];j<A->L->index[i+1];j++)
			{
				x[A->L->index[j]] -= conj(A->L->value[j]) * x[i];
			}
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<n;i++)
		{
			t   = x[i] * conj(A->WD->value[i]);
			for(j=A->U->index[i];j<A->U->index[i+1];j++)
			{
				x[A->U->index[j]] -= conj(A->U->value[j]) * t;
			}
		}
		for(i=n-1;i>=0;i--)
		{
			t    = x[i] * conj(A->WD->value[i]);
			x[i] = t;
			for(j=A->L->index[i];j<A->L->index[i+1];j++)
			{
				x[A->L->index[j]] -= conj(A->L->value[j]) * t;
			}
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_csr2msr"
LIS_INT lis_matrix_convert_csr2msr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,k,jj;
	LIS_INT err;
	LIS_INT n,nnz,ndz;
	LIS_INT count;
	LIS_INT *iw;
	LIS_INT *index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	nnz		= Ain->nnz;

	iw      = NULL;
	index   = NULL;
	value   = NULL;

	iw = (LIS_INT *)lis_malloc( (n+1)*sizeof(LIS_INT),"lis_matrix_convert_csr2msr::iw" );
	if( iw==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
		return LIS_ERR_OUT_OF_MEMORY;
	}

	/* check ndz */
	for(i=0;i<n+1;i++) iw[i] = 0;
	count = 0;
	#ifdef _OPENMP
	#pragma omp parallel private(i,j)
	#endif
	{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n;i++)
		{
			iw[i+1] = 0;
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				if( i==Ain->index[j] )
				{
					iw[i+1] = 1;
				}
			}
		}
		#ifdef _OPENMP
		#pragma omp for reduction(+:count)
		#endif
		for(i=0;i<n;i++)
		{
			count += iw[i+1];
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n;i++)
		{
			iw[i+1] = Ain->ptr[i+1]-Ain->ptr[i]-iw[i+1];
		}
	}
	ndz = n - count;

	err = lis_matrix_malloc_msr(n,nnz,ndz,&index,&value);
	if( err )
	{
		lis_free2(3,index,value,iw);
		return err;
	}

	/* convert msr */
	iw[0] = n+1;
	for(i=0;i<n;i++)
	{
		iw[i+1] = iw[i+1] + iw[i];
	}
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k)
	#endif
	{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n+1;i++)
		{
			index[i] = iw[i];
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n;i++)
		{
			k = index[i];
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				jj = Ain->index[j];
				if( jj==i )
				{
					value[i]   = Ain->value[j];
				}
				else
				{
					value[k]   = Ain->value[j];
					index[k]   = Ain->index[j];
					k++;
				}
			}
		}
	}

	err = lis_matrix_set_msr(nnz,ndz,index,value,Aout);
	if( err )
	{
		lis_free2(3,index,value,iw);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_free(iw);
		lis_matrix_storage_destroy(Aout);
		return err;
	}

	lis_free(iw);
	LIS_DEBUG_FUNC_OUT;

	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_msr2csr"
LIS_INT lis_matrix_convert_msr2csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,k;
	LIS_INT err;
	LIS_INT n,nnz;
	LIS_INT *ptr,*index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	nnz     = Ain->nnz;

	ptr     = NULL;
	index   = NULL;
	value   = NULL;

	err = lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);
	if( err )
	{
		return err;
	}

	/* convert csr */
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0;i<n;i++)
	{
		ptr[i+1] = Ain->index[i+1] - Ain->index[i];
		if( Ain->value[i]!=0.0 )
		{
			ptr[i+1]++;
		}
	}
	ptr[0] = 0;
	for(i=0;i<n;i++)
	{
		ptr[i+1] += ptr[i];
	}

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k)
	#endif
	for(i=0;i<n;i++)
	{
		k = ptr[i];
		if( Ain->value[i]!=(LIS_SCALAR)0.0 )
		{
			value[k]   = Ain->value[i];
			index[k]   = i;
			k++;
		}
		for(j=Ain->index[i];j<Ain->index[i+1];j++)
		{
			value[k]   = Ain->value[j];
			index[k]   = Ain->index[j];
			k++;
		}
	}

	err = lis_matrix_set_csr(nnz,ptr,index,value,Aout);
	if( err )
	{
		lis_free2(3,ptr,index,value);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

