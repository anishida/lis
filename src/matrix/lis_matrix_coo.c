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
#define __FUNC__ "lis_matrix_set_coo"
LIS_INT lis_matrix_set_coo(LIS_INT nnz, LIS_INT *row, LIS_INT *col, LIS_SCALAR *value, LIS_MATRIX A)
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

	A->row         = row;
	A->col         = col;
	A->value       = value;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_COO;
	A->nnz         = nnz;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_setDLU_coo"
LIS_INT lis_matrix_setDLU_coo(LIS_INT lnnz, LIS_INT unnz, LIS_SCALAR *diag, LIS_INT *lrow, LIS_INT *lcol, LIS_SCALAR *lvalue, LIS_INT *urow, LIS_INT *ucol, LIS_SCALAR *uvalue, LIS_MATRIX A)
{
	LIS_INT	err;
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

	A->L = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_coo::A->L");
	if( A->L==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	A->U = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_coo::A->U");
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
	A->L->row      = lrow;
	A->L->col      = lcol;
	A->L->value    = lvalue;
	A->U->nnz      = unnz;
	A->U->row      = urow;
	A->U->col      = ucol;
	A->U->value    = uvalue;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_COO;
	A->is_splited  = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_coo"
LIS_INT lis_matrix_malloc_coo(LIS_INT nnz, LIS_INT **row, LIS_INT **col, LIS_SCALAR **value)
{
	LIS_DEBUG_FUNC_IN;

	*row     = NULL;
	*col     = NULL;
	*value   = NULL;

	*row = (LIS_INT *)lis_malloc( nnz*sizeof(LIS_INT),"lis_matrix_malloc_coo::row" );
	if( *row==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		lis_free2(3,*row,*col,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*col = (LIS_INT *)lis_malloc( nnz*sizeof(LIS_INT),"lis_matrix_malloc_coo::col" );
	if( *col==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		lis_free2(3,*row,*col,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_malloc_coo::value" );
	if( *value==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		lis_free2(3,*row,*col,*value);
		return LIS_OUT_OF_MEMORY;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_elements_copy_coo"
LIS_INT lis_matrix_elements_copy_coo(LIS_INT nnz, LIS_INT *row, LIS_INT *col, LIS_SCALAR *value, LIS_INT *o_row, LIS_INT *o_col, LIS_SCALAR *o_value)
{
	LIS_INT	i;

	LIS_DEBUG_FUNC_IN;


	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0;i<nnz;i++)
	{
		o_row[i]     = row[i];
		o_col[i]     = col[i];
		o_value[i]   = value[i];
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_coo"
LIS_INT lis_matrix_copy_coo(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT i,n,nnz,lnnz,unnz;
	LIS_INT *row,*col;
	LIS_INT *lrow,*lcol;
	LIS_INT *urow,*ucol;
	LIS_SCALAR *value,*lvalue,*uvalue,*diag;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;

	if( Ain->is_splited )
	{
		lnnz     = Ain->L->nnz;
		unnz     = Ain->U->nnz;
		lrow     = NULL;
		lcol     = NULL;
		lvalue   = NULL;
		urow     = NULL;
		ucol     = NULL;
		uvalue   = NULL;
		diag     = NULL;

		err = lis_matrix_malloc_coo(lnnz,&lrow,&lcol,&lvalue);
		if( err )
		{
			return err;
		}
		err = lis_matrix_malloc_coo(unnz,&urow,&ucol,&uvalue);
		if( err )
		{
			lis_free2(7,diag,urow,lcol,urow,lcol,uvalue,lvalue);
			return err;
		}
		diag = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_matrix_copy_coo::diag");
		if( diag==NULL )
		{
			lis_free2(7,diag,urow,lcol,urow,lcol,uvalue,lvalue);
			return err;
		}

		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<n;i++)
		{
			diag[i] = Ain->D->value[i];
		}
		lis_matrix_elements_copy_coo(lnnz,Ain->L->row,Ain->L->col,Ain->L->value,lrow,lcol,lvalue);
		lis_matrix_elements_copy_coo(unnz,Ain->U->row,Ain->U->col,Ain->U->value,urow,ucol,uvalue);

		err = lis_matrix_setDLU_coo(lnnz,unnz,diag,lrow,lcol,lvalue,urow,ucol,uvalue,Aout);
		if( err )
		{
			lis_free2(7,diag,urow,lcol,urow,lcol,uvalue,lvalue);
			return err;
		}
	}
	if( !Ain->is_splited || (Ain->is_splited && Ain->is_save) )
	{
		row     = NULL;
		col     = NULL;
		value   = NULL;
		nnz     = Ain->nnz;
		err = lis_matrix_malloc_coo(nnz,&row,&col,&value);
		if( err )
		{
			return err;
		}

		lis_matrix_elements_copy_coo(nnz,Ain->row,Ain->col,Ain->value,row,col,value);

		err = lis_matrix_set_coo(nnz,row,col,value,Aout);
		if( err )
		{
			lis_free2(3,row,col,value);
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
#define __FUNC__ "lis_matrix_get_diagonal_coo"
LIS_INT lis_matrix_get_diagonal_coo(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i;
	LIS_INT n,nnz;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	nnz  = A->nnz;
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
		for(i=0; i<nnz;i++)
		{
			if( A->row[i]==A->col[i] )
			{
				d[A->row[i]] = A->value[i];
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_coo"
LIS_INT lis_matrix_shift_diagonal_coo(LIS_MATRIX A, LIS_SCALAR sigma)
{
	LIS_INT i;
	LIS_INT n,nnz;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	nnz  = A->nnz;
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
		for(i=0; i<nnz;i++)
		{
			if( A->row[i]==A->col[i] )
			{
				A->value[i] -= sigma;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_coo"
LIS_INT lis_matrix_scale_coo(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j;
	LIS_INT n,nnz;

	LIS_DEBUG_FUNC_IN;

	n = A->n;
	if( A->is_splited )
	{
		for(j=0; j<A->L->nnz; j++)
		{
			i = A->L->row[j];
			A->L->value[j] *= d[i];
		}
		for(i=0;i<n;i++) A->D->value[i] = 1.0;
		for(j=0; j<A->U->nnz; j++)
		{
			i = A->U->row[j];
			A->U->value[j] *= d[i];
		}
	}
	else
	{
		nnz = A->nnz;
		for(j=0; j<nnz; j++)
		{
			i = A->row[j];
			A->value[j] *= d[i];
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_symm_coo"
LIS_INT lis_matrix_scale_symm_coo(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j,k;
	LIS_INT n,nnz;

	LIS_DEBUG_FUNC_IN;

	n = A->n;
	if( A->is_splited )
	{
		for(k=0; k<A->L->nnz; k++)
		{
			i = A->L->row[k];
			j = A->L->row[k];
			A->L->value[k] *= d[i]*d[j];
		}
		for(i=0;i<n;i++) A->D->value[i] = 1.0;
		for(k=0; k<A->U->nnz; k++)
		{
			i = A->U->row[k];
			j = A->U->row[k];
			A->U->value[k] *= d[i]*d[j];
		}
	}
	else
	{
		nnz = A->nnz;
		for(k=0; k<nnz; k++)
		{
			i = A->row[k];
			j = A->row[k];
			A->value[k] *= d[i]*d[j];
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_normf_coo"
LIS_INT lis_matrix_normf_coo(LIS_MATRIX A, LIS_SCALAR *nrm)
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
#define __FUNC__ "lis_matrix_transpose_coo"
LIS_INT lis_matrix_transpose_coo(LIS_MATRIX Ain, LIS_MATRIX *Aout)
{

	LIS_DEBUG_FUNC_IN;

/*	err = lis_matrix_convert_coo2csc(Ain,Aout);*/
	(*Aout)->matrix_type = LIS_MATRIX_COO;
	(*Aout)->status      = LIS_MATRIX_COO;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_coo"
LIS_INT lis_matrix_split_coo(LIS_MATRIX A)
{
	LIS_INT i,nnz;
	LIS_INT nnzl,nnzu;
	LIS_INT err;
	LIS_INT *lrow,*lcol,*urow,*ucol;
	LIS_SCALAR *lvalue,*uvalue;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;

	nnz      = A->nnz;
	nnzl     = 0;
	nnzu     = 0;
	D        = NULL;
	lrow     = NULL;
	lcol     = NULL;
	lvalue   = NULL;
	urow     = NULL;
	ucol     = NULL;
	uvalue   = NULL;

	for(i=0;i<nnz;i++)
	{
		if( A->col[i]<A->row[i] )
		{
			nnzl++;
		}
		else if( A->col[i]>A->row[i] )
		{
			nnzu++;
		}
	}

	err = lis_matrix_LU_create(A);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_coo(nnzl,&lrow,&lcol,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_coo(nnzu,&urow,&ucol,&uvalue);
	if( err )
	{
		lis_free2(6,lrow,lcol,lvalue,urow,ucol,uvalue);
		return err;
	}
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		lis_free2(6,lrow,lcol,lvalue,urow,ucol,uvalue);
		return err;
	}

	nnzl = 0;
	nnzu = 0;
	for(i=0;i<nnz;i++)
	{
		if( A->col[i]<A->row[i] )
		{
			lrow[nnzl]   = A->row[i];
			lcol[nnzl]   = A->col[i];
			lvalue[nnzl] = A->value[i];
			nnzl++;
		}
		else if( A->col[i]>A->row[i] )
		{
			urow[nnzu]   = A->row[i];
			ucol[nnzu]   = A->col[i];
			uvalue[nnzu] = A->value[i];
			nnzu++;
		}
		else
		{
			D->value[A->row[i]] = A->value[i];
		}
	}
	A->L->nnz     = nnzl;
	A->L->row     = lrow;
	A->L->col     = lcol;
	A->L->value   = lvalue;
	A->U->nnz     = nnzu;
	A->U->row     = urow;
	A->U->col     = ucol;
	A->U->value   = uvalue;
	A->D          = D;
	A->is_splited = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_merge_coo"
LIS_INT lis_matrix_merge_coo(LIS_MATRIX A)
{
	LIS_INT i;
	LIS_INT nnz,nnzl,nnzu;
	LIS_INT err;
	LIS_INT *row,*col;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;


	nnzl    = A->L->nnz;
	nnzu    = A->U->nnz;
	row     = NULL;
	col     = NULL;
	value   = NULL;
	nnz     = A->L->nnz + A->U->nnz + A->D->n;

	err = lis_matrix_malloc_coo(nnz,&row,&col,&value);
	if( err )
	{
		return err;
	}

	nnz    = 0;
	for(i=0;i<nnzl;i++)
	{
		row[nnz]   = A->L->row[i];
		col[nnz]   = A->L->col[i];
		value[nnz] = A->L->value[i];
		nnz++;
	}
	for(i=0;i<A->D->n;i++)
	{
		row[nnz]   = i;
		col[nnz]   = i;
		value[nnz] = A->D->value[i];
		nnz++;
	}
	for(i=0;i<nnzu;i++)
	{
		row[nnz]   = A->U->row[i];
		col[nnz]   = A->U->col[i];
		value[nnz] = A->U->value[i];
		nnz++;
	}

	A->nnz        = nnz;
	A->row        = row;
	A->col        = col;
	A->value      = value;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_sort_coo"
LIS_INT lis_matrix_sort_coo(LIS_MATRIX A)
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
#define __FUNC__ "lis_matrix_convert_csr2coo"
LIS_INT lis_matrix_convert_csr2coo(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,k;
	LIS_INT err;
	LIS_INT n,nnz;
	LIS_INT *row,*col;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	nnz		= Ain->nnz;

	row     = NULL;
	col     = NULL;
	value   = NULL;

	err = lis_matrix_malloc_coo(nnz,&row,&col,&value);
	if( err )
	{
		return err;
	}

	/* convert coo */
	k    = 0;
	for(i=0;i<n;i++)
	{
		for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
		{
			row[k]     = i;
			col[k]     = Ain->index[j];
			value[k]   = Ain->value[j];
			k++;
		}
	}

	err = lis_matrix_set_coo(nnz,row,col,value,Aout);
	if( err )
	{
		lis_free2(3,row,col,value);
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

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_coo2csr"
LIS_INT lis_matrix_convert_coo2csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j;
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
	lis_sort_iid(0,nnz-1,Ain->row,Ain->col,Ain->value);
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0;i<n+1;i++)
	{
		ptr[i] = 0;
	}
	for(i=0;i<nnz;i++)
	{
		ptr[Ain->row[i]+1]++;
	}
	for(i=0;i<n;i++)
	{
		ptr[i+1] += ptr[i];
	}
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j)
	#endif
	for(i=0;i<n;i++)
	{
		for(j=ptr[i];j<ptr[i+1];j++)
		{
			value[j]   = Ain->value[j];
			index[j]   = Ain->col[j];
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
