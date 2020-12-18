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
#define __FUNC__ "lis_matrix_set_jad"
LIS_INT lis_matrix_set_jad(LIS_INT nnz, LIS_INT maxnzr, LIS_INT *perm, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A)
{
	LIS_INT i,n;
	LIS_INT *col;
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

	n = A->n;
	col = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_set_jad::col");
	if( col==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<n;i++)
	{
		col[perm[i]] = i;
	}

	A->col         = col;
	A->row         = perm;
	A->ptr         = ptr;
	A->index       = index;
	A->value       = value;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_JAD;
	A->nnz         = nnz;
	A->maxnzr      = maxnzr;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_setDLU_jad"
LIS_INT lis_matrix_setDLU_jad(LIS_INT lnnz, LIS_INT unnz, LIS_INT lmaxnzr, LIS_INT umaxnzr, LIS_SCALAR *diag, LIS_INT *lperm, LIS_INT *lptr, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uperm, LIS_INT *uptr, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A)
{
	LIS_INT n,i,err;
	LIS_INT *lcol,*ucol;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;

	n = A->n;

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

	A->L = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_jad::A->L");
	if( A->L==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	A->U = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_jad::A->U");
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
	lcol = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_setDLU_jad::lcol");
	if( lcol==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		lis_matrix_DLU_destroy(A);
		return LIS_OUT_OF_MEMORY;
	}
	ucol = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_setDLU_jad::ucol");
	if( ucol==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		lis_matrix_DLU_destroy(A);
		lis_free(lcol);
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<n;i++)
	{
		lcol[lperm[i]] = i;
		ucol[uperm[i]] = i;
	}

	lis_free(D->value);
	D->value       = diag;
	A->D           = D;
	A->L->nnz      = lnnz;
	A->L->maxnzr   = lmaxnzr;
	A->L->col      = lcol;
	A->L->row      = lperm;
	A->L->ptr      = lptr;
	A->L->index    = lindex;
	A->L->value    = lvalue;
	A->U->nnz      = unnz;
	A->U->maxnzr   = umaxnzr;
	A->U->col      = ucol;
	A->U->row      = uperm;
	A->U->ptr      = uptr;
	A->U->index    = uindex;
	A->U->value    = uvalue;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_JAD;
	A->is_splited  = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_jad"
LIS_INT lis_matrix_malloc_jad(LIS_INT n, LIS_INT nnz, LIS_INT maxnzr, LIS_INT **perm, LIS_INT **ptr, LIS_INT **index, LIS_SCALAR **value)
{
	LIS_INT	nprocs;

	LIS_DEBUG_FUNC_IN;

	*perm    = NULL;
	*ptr     = NULL;
	*index   = NULL;
	*value   = NULL;
	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
	#else
		nprocs  = 1;
	#endif

	*perm = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_malloc_jad::perm" );
	if( *perm==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		lis_free2(4,*perm,*ptr,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*ptr = (LIS_INT *)lis_malloc( nprocs*(maxnzr+1)*sizeof(LIS_INT),"lis_matrix_malloc_jad::ptr" );
	if( *ptr==NULL )
	{
		LIS_SETERR_MEM(nprocs*(maxnzr+1)*sizeof(LIS_INT));
		lis_free2(4,*perm,*ptr,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*index = (LIS_INT *)lis_malloc( nnz*sizeof(LIS_INT),"lis_matrix_malloc_jad::index" );
	if( *index==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		lis_free2(4,*perm,*ptr,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_malloc_jad::value" );
	if( *value==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		lis_free2(4,*perm,*ptr,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_elements_copy_jad"
LIS_INT lis_matrix_elements_copy_jad(LIS_INT n, LIS_INT maxnzr, LIS_INT *perm, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_perm, LIS_INT *o_ptr, LIS_INT *o_index, LIS_SCALAR *o_value)
{
	LIS_INT i,j,is,ie;
	LIS_INT nprocs,my_rank;

	LIS_DEBUG_FUNC_IN;

	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
	#else
		nprocs = 1;
	#endif

	#ifdef _OPENMP
	#pragma omp parallel private(i,j,is,ie,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		for(j=0;j<maxnzr+1;j++)
		{
			o_ptr[my_rank*(maxnzr+1) + j] = ptr[my_rank*(maxnzr+1) + j];
		}
		for(i=is;i<ie;i++)
		{
			o_perm[i] = perm[i];
		}
		for(j=0;j<maxnzr;j++)
		{
			for(i=ptr[my_rank*(maxnzr+1) + j];i<ptr[my_rank*(maxnzr+1) + j+1];i++)
			{
				o_value[i]   = value[i];
				o_index[i]   = index[i];
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_jad"
LIS_INT lis_matrix_copy_jad(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT i,n,nnz,lnnz,unnz,maxnzr,lmaxnzr,umaxnzr;
	LIS_INT *perm,*ptr,*index;
	LIS_INT *lperm,*lptr,*lindex;
	LIS_INT *uperm,*uptr,*uindex;
	LIS_SCALAR *value,*lvalue,*uvalue,*diag;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;

	if( Ain->is_splited )
	{
		lmaxnzr  = Ain->L->maxnzr;
		umaxnzr  = Ain->U->maxnzr;
		lnnz     = Ain->L->nnz;
		unnz     = Ain->U->nnz;
		lperm    = NULL;
		lptr     = NULL;
		lindex   = NULL;
		uperm    = NULL;
		uptr     = NULL;
		uindex   = NULL;
		diag     = NULL;

		err = lis_matrix_malloc_jad(n,lnnz,lmaxnzr,&lperm,&lptr,&lindex,&lvalue);
		if( err )
		{
			return err;
		}
		err = lis_matrix_malloc_jad(n,unnz,umaxnzr,&uperm,&uptr,&uindex,&uvalue);
		if( err )
		{
			lis_free2(9,diag,uperm,lperm,uptr,lptr,uindex,lindex,uvalue,lvalue);
			return err;
		}
		diag = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_matrix_copy_jad::diag");
		if( diag==NULL )
		{
			lis_free2(9,diag,uperm,lperm,uptr,lptr,uindex,lindex,uvalue,lvalue);
			return err;
		}

		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<n;i++)
		{
			diag[i] = Ain->D->value[i];
		}
		lis_matrix_elements_copy_jad(n,lmaxnzr,Ain->L->row,Ain->L->ptr,Ain->L->index,Ain->L->value,lperm,lptr,lindex,lvalue);
		lis_matrix_elements_copy_jad(n,umaxnzr,Ain->U->row,Ain->U->ptr,Ain->U->index,Ain->U->value,uperm,uptr,uindex,uvalue);

		err = lis_matrix_setDLU_jad(lnnz,unnz,lmaxnzr,umaxnzr,diag,lperm,lptr,lindex,lvalue,uperm,uptr,uindex,uvalue,Aout);
		if( err )
		{
			lis_free2(9,diag,uperm,lperm,uptr,lptr,uindex,lindex,uvalue,lvalue);
			return err;
		}
	}
	if( !Ain->is_splited || (Ain->is_splited && Ain->is_save) )
	{
		perm    = NULL;
		ptr     = NULL;
		index   = NULL;
		value   = NULL;
		maxnzr  = Ain->maxnzr;
		nnz     = Ain->nnz;
		err = lis_matrix_malloc_jad(n,nnz,maxnzr,&perm,&ptr,&index,&value);
		if( err )
		{
			return err;
		}

		lis_matrix_elements_copy_jad(n,maxnzr,Ain->row,Ain->ptr,Ain->index,Ain->value,perm,ptr,index,value);

		err = lis_matrix_set_jad(nnz,maxnzr,perm,ptr,index,value,Aout);
		if( err )
		{
			lis_free2(4,perm,ptr,index,value);
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
#define __FUNC__ "lis_matrix_get_diagonal_jad"
LIS_INT lis_matrix_get_diagonal_jad(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j,k,l;
	LIS_INT n,maxnzr;
	#ifdef _OPENMP
		LIS_INT is,ie,js,je;
		LIS_INT nprocs,my_rank;
	#endif

	LIS_DEBUG_FUNC_IN;

	n      = A->n;
	maxnzr = A->maxnzr;
	k      = n;
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
			nprocs = omp_get_max_threads();
			#pragma omp parallel private(i,j,k,l,is,ie,js,je,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				k = ie-is;
				for(j=0;j<maxnzr;j++)
				{
					l  = is;
					js = A->ptr[my_rank*(maxnzr+1) + j];
					je = A->ptr[my_rank*(maxnzr+1) + j+1];
					for(i=js;i<je;i++)
					{
						if( A->row[l]==A->index[i] )
						{
							d[A->row[l]] = A->value[i];
							k--;
							if( k==0 ) goto get_diag_end;
						}
						l++;
					}
				}
				get_diag_end:
				;
			}
		#else
			for(j=0;j<maxnzr;j++)
			{
				l = 0;
				for(i=A->ptr[j];i<A->ptr[j+1];i++)
				{
					if( A->row[l]==A->index[i] )
					{
						d[A->row[l]] = A->value[i];
						k--;
						if( k==0 ) return LIS_SUCCESS;
					}
					l++;
				}
			}
		#endif
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_jad"
LIS_INT lis_matrix_shift_diagonal_jad(LIS_MATRIX A, LIS_SCALAR sigma)
{
	LIS_INT i,j,k,l;
	LIS_INT n,maxnzr;
	#ifdef _OPENMP
		LIS_INT is,ie,js,je;
		LIS_INT nprocs,my_rank;
	#endif

	LIS_DEBUG_FUNC_IN;

	n      = A->n;
	maxnzr = A->maxnzr;
	k      = n;
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
			nprocs = omp_get_max_threads();
			#pragma omp parallel private(i,j,k,l,is,ie,js,je,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				k = ie-is;
				for(j=0;j<maxnzr;j++)
				{
					l  = is;
					js = A->ptr[my_rank*(maxnzr+1) + j];
					je = A->ptr[my_rank*(maxnzr+1) + j+1];
					for(i=js;i<je;i++)
					{
						if( A->row[l]==A->index[i] )
						{
							A->value[i] -= sigma;
							k--;
							if( k==0 ) goto get_diag_end;
						}
						l++;
					}
				}
				get_diag_end:
				;
			}
		#else
			for(j=0;j<maxnzr;j++)
			{
				l = 0;
				for(i=A->ptr[j];i<A->ptr[j+1];i++)
				{
					if( A->row[l]==A->index[i] )
					{
						A->value[i] -= sigma;
						k--;
						if( k==0 ) return LIS_SUCCESS;
					}
					l++;
				}
			}
		#endif
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_jad"
LIS_INT lis_matrix_scale_jad(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j,k,is,ie,js,je;
	LIS_INT n,maxnzr;
	LIS_INT nprocs,my_rank;

	LIS_DEBUG_FUNC_IN;

	n = A->n;
	if( A->is_splited )
	{
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
			for(i=is;i<ie;i++)
			{
				A->D->value[i] = 1.0;
			}
			for(j=0;j<A->L->maxnzr;j++)
			{
				k  = is;
				js = A->L->ptr[my_rank*(A->L->maxnzr+1) + j];
				je = A->L->ptr[my_rank*(A->L->maxnzr+1) + j+1];
				for(i=js;i<je;i++)
				{
					A->L->value[i] *= d[A->L->row[k]];
					k++;
				}
			}
			for(j=0;j<A->U->maxnzr;j++)
			{
				k  = is;
				js = A->U->ptr[my_rank*(A->U->maxnzr+1) + j];
				je = A->U->ptr[my_rank*(A->U->maxnzr+1) + j+1];
				for(i=js;i<je;i++)
				{
					A->U->value[i] *= d[A->U->row[k]];
					k++;
				}
			}
		}
	}
	else
	{
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
			for(j=0;j<maxnzr;j++)
			{
				k  = is;
				js = A->ptr[my_rank*(maxnzr+1) + j];
				je = A->ptr[my_rank*(maxnzr+1) + j+1];
				for(i=js;i<je;i++)
				{
					A->value[i] *= d[A->row[k]];
					k++;
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_symm_jad"
LIS_INT lis_matrix_scale_symm_jad(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j,k,is,ie,js,je;
	LIS_INT n,maxnzr;
	LIS_INT nprocs,my_rank;

	LIS_DEBUG_FUNC_IN;

	n = A->n;
	if( A->is_splited )
	{
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
			for(i=is;i<ie;i++)
			{
				A->D->value[i] = 1.0;
			}
			for(j=0;j<A->L->maxnzr;j++)
			{
				k  = is;
				js = A->L->ptr[my_rank*(A->L->maxnzr+1) + j];
				je = A->L->ptr[my_rank*(A->L->maxnzr+1) + j+1];
				for(i=js;i<je;i++)
				{
					A->L->value[i] *= d[A->L->row[k]]*d[A->L->index[i]];
					k++;
				}
			}
			for(j=0;j<A->U->maxnzr;j++)
			{
				k  = is;
				js = A->U->ptr[my_rank*(A->U->maxnzr+1) + j];
				je = A->U->ptr[my_rank*(A->U->maxnzr+1) + j+1];
				for(i=js;i<je;i++)
				{
					A->U->value[i] *= d[A->U->row[k]]*d[A->U->index[i]];
					k++;
				}
			}
		}
	}
	else
	{
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
			for(j=0;j<maxnzr;j++)
			{
				k  = is;
				js = A->ptr[my_rank*(maxnzr+1) + j];
				je = A->ptr[my_rank*(maxnzr+1) + j+1];
				for(i=js;i<je;i++)
				{
					A->value[i] *= d[A->row[k]]*d[A->index[i]];
					k++;
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_normf_jad"
LIS_INT lis_matrix_normf_jad(LIS_MATRIX A, LIS_SCALAR *nrm)
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
#define __FUNC__ "lis_matrix_transpose_jad"
LIS_INT lis_matrix_transpose_jad(LIS_MATRIX Ain, LIS_MATRIX *Aout)
{

	LIS_DEBUG_FUNC_IN;

/*	err = lis_matrix_convert_jad2csc(Ain,Aout);*/
	(*Aout)->matrix_type = LIS_MATRIX_JAD;
	(*Aout)->status      = LIS_MATRIX_JAD;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_jad"
LIS_INT lis_matrix_split_jad(LIS_MATRIX A)
{
	LIS_INT i,j,k,kk,n,maxnzr;
	LIS_INT lnnz,unnz,lmaxnzr,umaxnzr;
	LIS_INT err;
	LIS_INT *liw,*uiw,*liw2,*uiw2;
	LIS_INT *lperm,*lptr,*lindex;
	LIS_INT *uperm,*uptr,*uindex;
	#ifdef _OPENMP
		LIS_INT *iw;
		LIS_INT my_rank,nprocs,is,ie,js,je;
	#endif
	LIS_SCALAR *lvalue,*uvalue;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;

	n        = A->n;
	maxnzr   = A->maxnzr;
	lmaxnzr  = 0;
	umaxnzr  = 0;
	lnnz     = 0;
	unnz     = 0;
	liw      = NULL;
	uiw      = NULL;
	liw2     = NULL;
	uiw2     = NULL;
	D        = NULL;
	lperm    = NULL;
	lptr     = NULL;
	lindex   = NULL;
	lvalue   = NULL;
	uperm    = NULL;
	uptr     = NULL;
	uindex   = NULL;
	uvalue   = NULL;

	liw = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_split_jad::liw");
	if( liw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	uiw = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_split_jad::uiw");
	if( uiw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		lis_free2(12,liw,uiw,liw2,uiw2,lperm,lptr,lindex,lvalue,uperm,uptr,uindex,uvalue);
		return LIS_ERR_OUT_OF_MEMORY;
	}
	liw2 = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_split_jad::liw2");
	if( liw2==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		lis_free2(12,liw,uiw,liw2,uiw2,lperm,lptr,lindex,lvalue,uperm,uptr,uindex,uvalue);
		return LIS_ERR_OUT_OF_MEMORY;
	}
	uiw2 = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_split_jad::uiw2");
	if( uiw2==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		lis_free2(12,liw,uiw,liw2,uiw2,lperm,lptr,lindex,lvalue,uperm,uptr,uindex,uvalue);
		return LIS_ERR_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		iw = (LIS_INT *)lis_malloc((nprocs+1)*LIS_VEC_TMP_PADD*sizeof(LIS_INT),"lis_matrix_split_jad::iw");
		if( iw==NULL )
		{
			LIS_SETERR_MEM((nprocs+1)*LIS_VEC_TMP_PADD*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		#pragma omp parallel private(i,j,k,is,ie,js,je,my_rank)
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			memset(&liw[is],0,(ie-is)*sizeof(LIS_INT));
			memset(&uiw[is],0,(ie-is)*sizeof(LIS_INT));
			iw[my_rank*LIS_VEC_TMP_PADD]   = 0;
			iw[my_rank*LIS_VEC_TMP_PADD+1] = 0;
			iw[my_rank*LIS_VEC_TMP_PADD+2] = 0;
			iw[my_rank*LIS_VEC_TMP_PADD+3] = 0;
			for(j=0;j<maxnzr;j++)
			{
				k = is;
				js = A->ptr[my_rank*(maxnzr+1) + j];
				je = A->ptr[my_rank*(maxnzr+1) + j+1];
				for(i=js;i<je;i++)
				{
					if( A->index[i]<A->row[k] )
					{
						iw[my_rank*LIS_VEC_TMP_PADD+2]++;
						liw[k]++;
					}
					else if( A->index[i]>A->row[k] )
					{
						iw[my_rank*LIS_VEC_TMP_PADD+3]++;
						uiw[k]++;
					}
					k++;
				}
			}
			for(i=is;i<ie;i++)
			{
				if( iw[my_rank*LIS_VEC_TMP_PADD]<liw[i]   ) iw[my_rank*LIS_VEC_TMP_PADD]   = liw[i];
				if( iw[my_rank*LIS_VEC_TMP_PADD+1]<uiw[i] ) iw[my_rank*LIS_VEC_TMP_PADD+1] = uiw[i];
			}
		}
		iw[4] = 0;
		iw[5] = 0;
		for(i=0;i<nprocs;i++)
		{
			if( iw[i*LIS_VEC_TMP_PADD]>lmaxnzr   ) lmaxnzr = iw[i*LIS_VEC_TMP_PADD];
			if( iw[i*LIS_VEC_TMP_PADD+1]>umaxnzr ) umaxnzr = iw[i*LIS_VEC_TMP_PADD+1];
			iw[(i+1)*LIS_VEC_TMP_PADD+4] = iw[i*LIS_VEC_TMP_PADD+4] + iw[i*LIS_VEC_TMP_PADD+2];
			iw[(i+1)*LIS_VEC_TMP_PADD+5] = iw[i*LIS_VEC_TMP_PADD+5] + iw[i*LIS_VEC_TMP_PADD+3];
		}
		lnnz = iw[nprocs*LIS_VEC_TMP_PADD+4];
		unnz = iw[nprocs*LIS_VEC_TMP_PADD+5];
	#else
		memset(liw,0,n*sizeof(LIS_INT));
		memset(uiw,0,n*sizeof(LIS_INT));
		for(j=0;j<maxnzr;j++)
		{
			k = 0;
			for(i=A->ptr[j];i<A->ptr[j+1];i++)
			{
				if( A->index[i]<A->row[k] )
				{
					lnnz++;
					liw[k]++;
				}
				else if( A->index[i]>A->row[k] )
				{
					unnz++;
					uiw[k]++;
				}
				k++;
			}
		}
		for(i=0;i<n;i++)
		{
			if( lmaxnzr<liw[i] ) lmaxnzr = liw[i];
			if( umaxnzr<uiw[i] ) umaxnzr = uiw[i];
		}
	#endif

	err = lis_matrix_LU_create(A);
	if( err )
	{
		lis_free2(12,liw,uiw,liw2,uiw2,lperm,lptr,lindex,lvalue,uperm,uptr,uindex,uvalue);
		return err;
	}
	err = lis_matrix_malloc_jad(n,lnnz,lmaxnzr,&lperm,&lptr,&lindex,&lvalue);
	if( err )
	{
		lis_free2(12,liw,uiw,liw2,uiw2,lperm,lptr,lindex,lvalue,uperm,uptr,uindex,uvalue);
		return err;
	}
	err = lis_matrix_malloc_jad(n,unnz,umaxnzr,&uperm,&uptr,&uindex,&uvalue);
	if( err )
	{
		lis_free2(12,liw,uiw,liw2,uiw2,lperm,lptr,lindex,lvalue,uperm,uptr,uindex,uvalue);
		return err;
	}
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		lis_free2(12,liw,uiw,liw2,uiw2,lperm,lptr,lindex,lvalue,uperm,uptr,uindex,uvalue);
		return err;
	}

	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		#pragma omp parallel private(i,j,k,is,ie,js,je,my_rank)
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			memset(&lptr[my_rank*(lmaxnzr+1)],0,(lmaxnzr+1)*sizeof(LIS_INT));
			memset(&uptr[my_rank*(umaxnzr+1)],0,(umaxnzr+1)*sizeof(LIS_INT));
			for(i=is;i<ie;i++)
			{
				lperm[i] = A->row[i];
				uperm[i] = A->row[i];
				for(j=0;j<liw[i];j++)
				{
					lptr[my_rank*(lmaxnzr+1) + j+1]++;
				}
				for(j=0;j<uiw[i];j++)
				{
					uptr[my_rank*(umaxnzr+1) + j+1]++;
				}
			}
			lis_sortr_ii(is,ie-1,liw,lperm);
			lis_sortr_ii(is,ie-1,uiw,uperm);
			lptr[my_rank*(lmaxnzr+1)] = iw[my_rank*LIS_VEC_TMP_PADD+4];
			uptr[my_rank*(umaxnzr+1)] = iw[my_rank*LIS_VEC_TMP_PADD+5];
			for(j=0;j<lmaxnzr;j++)
			{
				lptr[my_rank*(lmaxnzr+1) + j+1] += lptr[my_rank*(lmaxnzr+1) + j];
			}
			for(j=0;j<umaxnzr;j++)
			{
				uptr[my_rank*(umaxnzr+1) + j+1] += uptr[my_rank*(umaxnzr+1) + j];
			}

			for(i=is;i<ie;i++)
			{
				liw[i]         = 0;
				uiw[i]         = 0;
				liw2[lperm[i]] = i;
				uiw2[uperm[i]] = i;
			}
			for(j=0;j<maxnzr;j++)
			{
				k = is;
				for(i=A->ptr[my_rank*(maxnzr+1) + j];i<A->ptr[my_rank*(maxnzr+1) + j+1];i++)
				{
					if( A->index[i]<A->row[k] )
					{
						kk              = lptr[my_rank*(lmaxnzr+1) + liw[A->row[k]]] + liw2[A->row[k]] - is;
						liw[A->row[k]]++;
						lindex[kk]      = A->index[i];
						lvalue[kk]      = A->value[i];
					}
					else if( A->index[i]>A->row[k] )
					{
						kk              = uptr[my_rank*(umaxnzr+1) + uiw[A->row[k]]] + uiw2[A->row[k]] - is;
						uiw[A->row[k]]++;
						uindex[kk]      = A->index[i];
						uvalue[kk]      = A->value[i];
					}
					else
					{
						D->value[A->row[k]] = A->value[i];
					}
					k++;
				}
			}
		}
		lis_free(iw);
	#else
		memset(lptr,0,(lmaxnzr+1)*sizeof(LIS_INT));
		memset(uptr,0,(umaxnzr+1)*sizeof(LIS_INT));
		for(i=0;i<n;i++)
		{
			lperm[i] = A->row[i];
			uperm[i] = A->row[i];
			for(j=0;j<liw[i];j++)
			{
				lptr[j+1]++;
			}
			for(j=0;j<uiw[i];j++)
			{
				uptr[j+1]++;
			}
		}
		lis_sortr_ii(0,n-1,liw,lperm);
		lis_sortr_ii(0,n-1,uiw,uperm);
		for(j=0;j<lmaxnzr;j++)
		{
			lptr[j+1] += lptr[j];
		}
		for(j=0;j<umaxnzr;j++)
		{
			uptr[j+1] += uptr[j];
		}

		for(i=0;i<n;i++)
		{
			liw[i]         = 0;
			uiw[i]         = 0;
			liw2[lperm[i]] = i;
			uiw2[uperm[i]] = i;
		}
		for(j=0;j<maxnzr;j++)
		{
			k = 0;
			for(i=A->ptr[j];i<A->ptr[j+1];i++)
			{
				if( A->index[i]<A->row[k] )
				{
					kk              = lptr[liw[A->row[k]]] + liw2[A->row[k]];
					liw[A->row[k]]++;
					lindex[kk]      = A->index[i];
					lvalue[kk]      = A->value[i];
				}
				else if( A->index[i]>A->row[k] )
				{
					kk              = uptr[uiw[A->row[k]]] + uiw2[A->row[k]];
					uiw[A->row[k]]++;
					uindex[kk]      = A->index[i];
					uvalue[kk]      = A->value[i];
				}
				else
				{
					D->value[A->row[k]] = A->value[i];
				}
				k++;
			}
		}
	#endif
	A->L->nnz     = lnnz;
	A->L->maxnzr  = lmaxnzr;
	A->L->col     = liw2;
	A->L->row     = lperm;
	A->L->ptr     = lptr;
	A->L->index   = lindex;
	A->L->value   = lvalue;
	A->U->nnz     = unnz;
	A->U->maxnzr  = umaxnzr;
	A->U->col     = uiw2;
	A->U->row     = uperm;
	A->U->ptr     = uptr;
	A->U->index   = uindex;
	A->U->value   = uvalue;
	A->D          = D;
	A->is_splited = LIS_TRUE;

	lis_free2(2,liw,uiw);
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_merge_jad"
LIS_INT lis_matrix_merge_jad(LIS_MATRIX A)
{
	LIS_INT i,j,k,kk,ie;
	LIS_INT n,nnz,maxnzr;
	LIS_INT err;
	LIS_INT *perm,*ptr,*index,*iw,*iw2;
	LIS_SCALAR *value;
	#ifdef _OPENMP
		LIS_INT is,js,je,nprocs,my_rank;
		LIS_INT *iw3;
	#endif

	LIS_DEBUG_FUNC_IN;


	n       = A->n;
	perm    = NULL;
	ptr     = NULL;
	index   = NULL;
	value   = NULL;
	iw      = NULL;
	iw2     = NULL;
	nnz     = A->L->nnz + A->U->nnz + n;

	iw = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_merge_jad::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw2 = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_merge_jad::iw2");
	if( iw2==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		lis_free2(2,iw,iw2);
		return LIS_OUT_OF_MEMORY;
	}
	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		iw3 = (LIS_INT *)lis_malloc((nprocs+1)*LIS_VEC_TMP_PADD*sizeof(LIS_INT),"lis_matrix_merge_jad::iw3");
		if( iw3==NULL )
		{
			LIS_SETERR_MEM((nprocs+1)*LIS_VEC_TMP_PADD*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		#pragma omp parallel private(i,j,k,is,ie,js,je,my_rank)
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			iw3[my_rank*LIS_VEC_TMP_PADD]   = 0;
			iw3[my_rank*LIS_VEC_TMP_PADD+1] = (ie-is);
			for(i=is;i<ie;i++)
			{
				iw[i] = 1;
			}
			for(j=0;j<A->L->maxnzr;j++)
			{
				k  = is;
				js = A->L->ptr[my_rank*(A->L->maxnzr+1) + j];
				je = A->L->ptr[my_rank*(A->L->maxnzr+1) + j+1];
				for(i=js;i<je;i++)
				{
					iw3[my_rank*LIS_VEC_TMP_PADD+1]++;
					iw[A->L->row[k++]]++;
				}
			}
			for(j=0;j<A->U->maxnzr;j++)
			{
				k  = is;
				js = A->U->ptr[my_rank*(A->U->maxnzr+1) + j];
				je = A->U->ptr[my_rank*(A->U->maxnzr+1) + j+1];
				for(i=js;i<je;i++)
				{
					iw3[my_rank*LIS_VEC_TMP_PADD+1]++;
					iw[A->U->row[k++]]++;
				}
			}
			for(i=is;i<ie;i++)
			{
				if( iw3[my_rank*LIS_VEC_TMP_PADD]<iw[i] ) iw3[my_rank*LIS_VEC_TMP_PADD] = iw[i];
			}
		}
		maxnzr = 0;
		iw3[2]  = 0;
		for(i=0;i<nprocs;i++)
		{
			if( iw3[i*LIS_VEC_TMP_PADD]>maxnzr ) maxnzr = iw3[i*LIS_VEC_TMP_PADD];
			iw3[(i+1)*LIS_VEC_TMP_PADD+2] = iw3[i*LIS_VEC_TMP_PADD+2] + iw3[i*LIS_VEC_TMP_PADD+1];
		}
	#else
		for(i=0;i<n;i++)
		{
			iw[i] = 1;
		}
		for(j=0;j<A->L->maxnzr;j++)
		{
			ie = A->L->ptr[j+1] - A->L->ptr[j];
			for(i=0;i<ie;i++)
			{
				iw[A->L->row[i]]++;
			}
		}
		for(j=0;j<A->U->maxnzr;j++)
		{
			ie = A->U->ptr[j+1] - A->U->ptr[j];
			for(i=0;i<ie;i++)
			{
				iw[A->U->row[i]]++;
			}
		}
		maxnzr = 0;
		for(i=0;i<n;i++)
		{
			if( maxnzr<iw[i] ) maxnzr = iw[i];
		}
	#endif

	err = lis_matrix_malloc_jad(n,nnz,maxnzr,&perm,&ptr,&index,&value);
	if( err )
	{
		lis_free2(2,iw,iw2);
		return err;
	}

	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		#pragma omp parallel private(i,j,k,is,ie,js,je,my_rank)
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			memset(&ptr[my_rank*(maxnzr+1)],0,(maxnzr+1)*sizeof(LIS_INT));
			for(i=is;i<ie;i++)
			{
				perm[i] = i;
				for(j=0;j<iw[i];j++)
				{
					ptr[my_rank*(maxnzr+1) + j+1]++;
				}
			}
			lis_sortr_ii(is,ie-1,iw,perm);
			ptr[my_rank*(maxnzr+1)] = iw3[my_rank*LIS_VEC_TMP_PADD+2];
			for(j=0;j<maxnzr;j++)
			{
				ptr[my_rank*(maxnzr+1) + j+1] += ptr[my_rank*(maxnzr+1) + j];
			}

			for(i=is;i<ie;i++)
			{
				iw[i]        = 0;
				iw2[perm[i]] = i;
			}
			for(j=0;j<A->L->maxnzr;j++)
			{
				k = is;
				for(i=A->L->ptr[my_rank*(A->L->maxnzr+1) + j];i<A->L->ptr[my_rank*(A->L->maxnzr+1) + j+1];i++)
				{
					kk                 = ptr[my_rank*(maxnzr+1) + iw[A->L->row[k]]] + iw2[A->L->row[k]] - is;
					iw[A->L->row[k]]++;
					index[kk]          = A->L->index[i];
					value[kk]          = A->L->value[i];
					k++;
				}
			}
			for(i=is;i<ie;i++)
			{
				kk         = ptr[my_rank*(maxnzr+1) + iw[i]] + iw2[i] - is;
				iw[i]++;
				index[kk]  = i;
				value[kk]  = A->D->value[i];
			}
			for(j=0;j<A->U->maxnzr;j++)
			{
				k = is;
				for(i=A->U->ptr[my_rank*(A->U->maxnzr+1) + j];i<A->U->ptr[my_rank*(A->U->maxnzr+1) + j+1];i++)
				{
					kk                 = ptr[my_rank*(maxnzr+1) + iw[A->U->row[k]]] + iw2[A->U->row[k]] - is;
					iw[A->U->row[k]]++;
					index[kk]          = A->U->index[i];
					value[kk]          = A->U->value[i];
					k++;
				}
			}
		}
	#else
		memset(ptr,0,(maxnzr+1)*sizeof(LIS_INT));
		for(i=0;i<n;i++)
		{
			perm[i] = i;
			for(j=0;j<iw[i];j++)
			{
				ptr[j+1]++;
			}
		}
		lis_sortr_ii(0,n-1,iw,perm);
		for(j=0;j<maxnzr;j++)
		{
			ptr[j+1] += ptr[j];
		}

		for(i=0;i<n;i++)
		{
			iw[i]        = 0;
			iw2[perm[i]] = i;
		}
		for(j=0;j<A->L->maxnzr;j++)
		{
			k = 0;
			for(i=A->L->ptr[j];i<A->L->ptr[j+1];i++)
			{
				kk                 = ptr[iw[A->L->row[k]]] + iw2[A->L->row[k]];
				iw[A->L->row[k]]++;
				index[kk]          = A->L->index[i];
				value[kk]          = A->L->value[i];
				k++;
			}
		}
		for(i=0;i<n;i++)
		{
			kk         = ptr[iw[i]] + iw2[i];
			iw[i]++;
			index[kk]  = i;
			value[kk]  = A->D->value[i];
		}
		for(j=0;j<A->U->maxnzr;j++)
		{
			k = 0;
			for(i=A->U->ptr[j];i<A->U->ptr[j+1];i++)
			{
				kk                 = ptr[iw[A->U->row[k]]] + iw2[A->U->row[k]];
				iw[A->U->row[k]]++;
				index[kk]          = A->U->index[i];
				value[kk]          = A->U->value[i];
				k++;
			}
		}
	#endif

	A->nnz        = nnz;
	A->row        = perm;
	A->ptr        = ptr;
	A->value      = value;
	A->index      = index;

	lis_free2(2,iw,iw2);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_sort_jad"
LIS_INT lis_matrix_sort_jad(LIS_MATRIX A)
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
#define __FUNC__ "lis_matrix_solve_jad"
LIS_INT lis_matrix_solve_jad(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,k,l,n;
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
			k = A->L->col[i];
			l = 0;
			j = A->L->ptr[l++] + k;
			t = b[i];
			while( j<A->L->ptr[l] && l<=A->L->maxnzr )
			{
				t -= A->L->value[j] * x[A->L->index[j]];
				j = A->L->ptr[l++] + k;
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			k = A->U->col[i];
			l = 0;
			j = A->U->ptr[l++] + k;
			t = b[i];
			while( j<A->U->ptr[l] && l<=A->U->maxnzr )
			{
				t -= A->U->value[j] * x[A->U->index[j]];
				j = A->U->ptr[l++] + k;
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<n;i++)
		{
			k = A->L->col[i];
			l = 0;
			j = A->L->ptr[l++] + k;
			t = b[i];
			while( j<A->L->ptr[l] && l<=A->L->maxnzr )
			{
				t -= A->L->value[j] * x[A->L->index[j]];
				j = A->L->ptr[l++] + k;
			}
			x[i]   = t * A->WD->value[i];
		}
		for(i=n-1;i>=0;i--)
		{
			k = A->U->col[i];
			l = 0;
			j = A->U->ptr[l++] + k;
			t = 0.0;
			while( j<A->U->ptr[l] && l<=A->U->maxnzr )
			{
				t += A->U->value[j] * x[A->U->index[j]];
				j = A->U->ptr[l++] + k;
			}
			x[i]  -= t * A->WD->value[i];
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solveh_jad"
LIS_INT lis_matrix_solveh_jad(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,k,l,n;
	LIS_SCALAR t;
	LIS_SCALAR *x;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	x  = X->value;

	lis_vector_copy(B,X);
	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=0;i<n;i++)
		{
			k = A->U->col[i];
			l = 0;
			j = A->U->ptr[l++] + k;
			x[i]   = x[i] * conj(A->WD->value[i]);
			while( j<A->U->ptr[l] && l<=A->U->maxnzr )
			{
				x[A->U->index[j]] -= conj(A->U->value[j]) * x[i];
				j = A->U->ptr[l++] + k;
			}
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			k = A->L->col[i];
			l = 0;
			j = A->L->ptr[l++] + k;
			x[i]   = x[i] * conj(A->WD->value[i]);
			while( j<A->L->ptr[l] && l<=A->L->maxnzr )
			{
				x[A->L->index[j]] -= conj(A->L->value[j]) * x[i];
				j = A->L->ptr[l++] + k;
			}
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<n;i++)
		{
			k = A->U->col[i];
			l = 0;
			j = A->U->ptr[l++] + k;
			t = x[i] * conj(A->WD->value[i]);
			while( j<A->U->ptr[l] && l<=A->U->maxnzr )
			{
				x[A->U->index[j]] -= conj(A->U->value[j]) * t;
				j = A->U->ptr[l++] + k;
			}
		}
		for(i=n-1;i>=0;i--)
		{
			k = A->L->col[i];
			l = 0;
			j = A->L->ptr[l++] + k;
			x[i] = x[i] * conj(A->WD->value[i]);
			t    = x[i];
			while( j<A->L->ptr[l] && l<=A->L->maxnzr )
			{
				x[A->L->index[j]] -= conj(A->L->value[j]) * t;
				j = A->L->ptr[l++] + k;
			}
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#ifndef USE_OVERLAP
#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_csr2jad"
LIS_INT lis_matrix_convert_csr2jad(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,l,js,je;
	LIS_INT err;
	LIS_INT n,nnz,maxnzr,nprocs,my_rank;
	LIS_INT is,ie;
	LIS_INT *iw,*maxnzrpe,*nnzpe;
	LIS_INT *perm,*ptr,*index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	nnz	= Ain->nnz;
	my_rank = Ain->my_rank;
	is      = Ain->is;
	ie      = Ain->ie;
	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
	#else
		nprocs  = 1;
	#endif

	perm     = NULL;
	ptr      = NULL;
	index    = NULL;
	value    = NULL;
	iw       = NULL;
	maxnzrpe = NULL;
	nnzpe    = NULL;

	iw = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_convert_csr2jad::iw" );
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	maxnzrpe = (LIS_INT *)lis_malloc( nprocs*sizeof(LIS_INT),"lis_matrix_convert_csr2jad::maxnzrpe" );
	if( maxnzrpe==NULL )
	{
		LIS_SETERR_MEM(nprocs*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	nnzpe = (LIS_INT *)lis_malloc( (nprocs+1)*sizeof(LIS_INT),"lis_matrix_convert_csr2jad::nnzpe" );
	if( nnzpe==NULL )
	{
		LIS_SETERR_MEM((nprocs+1)*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
	#pragma omp parallel private(i,is,ie,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		maxnzrpe[my_rank] = 0;
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(i=is;i<ie;i++)
		{
			iw[i] = Ain->ptr[i+1] - Ain->ptr[i];
			if( iw[i] > maxnzrpe[my_rank] ) maxnzrpe[my_rank] = iw[i];
		}
		nnzpe[my_rank+1] = Ain->ptr[ie] - Ain->ptr[is];
	}
	maxnzr    = 0;
	nnzpe[0]  = 0;
	for(i=0;i<nprocs;i++)
	{
		if( maxnzrpe[i] > maxnzr ) maxnzr = maxnzrpe[i];
		nnzpe[i+1] += nnzpe[i];
	}

	err = lis_matrix_malloc_jad(n,nnz,maxnzr,&perm,&ptr,&index,&value);
	if( err )
	{
		return err;
	}

	/* convert jad */
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,is,ie,js,je,l,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		memset(&ptr[my_rank*(maxnzr+1)],0,(maxnzr+1)*sizeof(LIS_INT));
		#if 1
			for(i=is;i<ie;i++)
			{
				perm[i] = i;
				for(j=0;j<iw[i];j++)
				{
					ptr[my_rank*(maxnzr+1) + j+1]++;
				}
			}
			lis_sortr_ii(is,ie-1,iw,perm);
		#else
			lis_sort_jad(is,ie,maxnzr,iw,perm);
		
			j = 0;
			for(i=ie-1;i>=is;i--)
			{
				for(;j<iw[i];j++)
				{
					ptr[my_rank*(maxnzr+1) + j+1] = i-is+1;
				}
				if( iw[i]==maxnzr ) break;
			}
			
		#endif

		ptr[my_rank*(maxnzr+1)] = nnzpe[my_rank];
		for(j=0;j<maxnzr;j++)
		{
			ptr[my_rank*(maxnzr+1) + j+1] += ptr[my_rank*(maxnzr+1) + j];
		}



		#ifndef USE_VEC_COMP
			for(i=is;i<ie;i++)
			{
				js = Ain->ptr[perm[i]];
				je = Ain->ptr[perm[i]+1];
				for(j=js;j<je;j++)
				{
					l          = ptr[my_rank*(maxnzr+1) + j-js]+i-is;
					value[l]   = Ain->value[j];
					index[l]   = Ain->index[j];
				}
			}
		#else
			for(j=0;j<maxnzr;j++)
			{
				js = ptr[my_rank*(maxnzr+1) + j];
				je = ptr[my_rank*(maxnzr+1) + j+1];
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=js;i<je;i++)
				{
					l          = Ain->ptr[perm[is+(i-js)]] + j;
					value[i]   = Ain->value[l];
					index[i]   = Ain->index[l];
				}
			}
		#endif
	}

	err = lis_matrix_set_jad(nnz,maxnzr,perm,ptr,index,value,Aout);
	if( err )
	{
		lis_free2(7,perm,ptr,index,value,iw,maxnzrpe,nnzpe);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_free2(2,iw,nnzpe);
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	lis_free2(3,iw,nnzpe,maxnzrpe);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#else
#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_csr2jad"
LIS_INT lis_matrix_convert_csr2jad(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,jj,k,kk,l,js,je;
	LIS_INT err;
	LIS_INT np,n,nnz,nnz2,maxnzr,maxnzr2,nprocs,my_rank;
	LIS_INT is,ie,pe;
	LIS_INT *iw,*maxnzrpe,*nnzpe;
	LIS_INT *iw2,*maxnzrpe2,*nnzpe2;
	LIS_INT *perm,*ptr,*index;
	LIS_INT *perm2,*ptr2,*index2;
	LIS_SCALAR *value;
	LIS_SCALAR *value2;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	np      = Ain->np;
	nnz		= Ain->nnz;
	my_rank = Ain->my_rank;
	is      = Ain->is;
	ie      = Ain->ie;
	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
	#else
		nprocs  = 1;
	#endif

	perm     = NULL;
	ptr      = NULL;
	index    = NULL;
	value    = NULL;
	iw       = NULL;
	iw2      = NULL;
	maxnzrpe = NULL;
	nnzpe    = NULL;

	lis_matrix_split2_csr(Ain);

	iw = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_convert_csr2jad::iw" );
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw2 = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_convert_csr2jad::iw2" );
	if( iw2==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	maxnzrpe = (LIS_INT *)lis_malloc( nprocs*sizeof(LIS_INT),"lis_matrix_convert_csr2jad::maxnzrpe" );
	if( maxnzrpe==NULL )
	{
		LIS_SETERR_MEM(nprocs*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	maxnzrpe2 = (LIS_INT *)lis_malloc( nprocs*sizeof(LIS_INT),"lis_matrix_convert_csr2jad::maxnzrpe2" );
	if( maxnzrpe2==NULL )
	{
		LIS_SETERR_MEM(nprocs*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	nnzpe = (LIS_INT *)lis_malloc( (nprocs+1)*sizeof(LIS_INT),"lis_matrix_convert_csr2jad::nnzpe" );
	if( nnzpe==NULL )
	{
		LIS_SETERR_MEM((nprocs+1)*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	nnzpe2 = (LIS_INT *)lis_malloc( (nprocs+1)*sizeof(LIS_INT),"lis_matrix_convert_csr2jad::nnzpe2" );
	if( nnzpe2==NULL )
	{
		LIS_SETERR_MEM((nprocs+1)*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k,is,ie,jj,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		maxnzrpe[my_rank]  = 0;
		maxnzrpe2[my_rank] = 0;
		for(i=is;i<ie;i++)
		{
			iw[i]  = Ain->L->ptr[i+1] - Ain->L->ptr[i];
			iw2[i] = Ain->U->ptr[i+1] - Ain->U->ptr[i];
			if( iw[i]  > maxnzrpe[my_rank]  ) maxnzrpe[my_rank]  = iw[i];
			if( iw2[i] > maxnzrpe2[my_rank] ) maxnzrpe2[my_rank] = iw2[i];
		}
		nnzpe[my_rank+1]  = Ain->L->ptr[ie] - Ain->L->ptr[is];
		nnzpe2[my_rank+1] = Ain->U->ptr[ie] - Ain->U->ptr[is];
	}
	maxnzr     = 0;
	maxnzr2    = 0;
	nnzpe[0]   = 0;
	nnzpe2[0]  = 0;
	for(i=0;i<nprocs;i++)
	{
		if( maxnzrpe[i] >  maxnzr  ) maxnzr  = maxnzrpe[i];
		if( maxnzrpe2[i] > maxnzr2 ) maxnzr2 = maxnzrpe2[i];
		nnzpe[i+1]  += nnzpe[i];
		nnzpe2[i+1] += nnzpe2[i];
	}
	nnz  = nnzpe[nprocs];
	nnz2 = nnzpe2[nprocs];

	err = lis_matrix_malloc_jad(n,nnz,maxnzr,&perm,&ptr,&index,&value);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_jad(n,nnz2,maxnzr2,&perm2,&ptr2,&index2,&value2);
	if( err )
	{
		return err;
	}
	err = lis_matrix_LU_create(Aout);
	if( err )
	{
		lis_free2(7,perm,ptr,index,value,iw,maxnzrpe,nnzpe);
		return err;
	}

	/* convert jad */
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k,is,ie,jj,kk,js,je,l,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		memset(&ptr[my_rank*(maxnzr+1)],0,(maxnzr+1)*sizeof(LIS_INT));
		memset(&ptr2[my_rank*(maxnzr2+1)],0,(maxnzr2+1)*sizeof(LIS_INT));
		#if 0
			for(i=is;i<ie;i++)
			{
				perm[i] = i;
				perm2[i] = i;
				for(j=0;j<iw[i];j++)
				{
					ptr[my_rank*(maxnzr+1) + j+1]++;
				}
				for(j=0;j<iw2[i];j++)
				{
					ptr2[my_rank*(maxnzr2+1) + j+1]++;
				}
			}
			lis_sortr_ii(is,ie-1,iw,perm);
			lis_sortr_ii(is,ie-1,iw2,perm2);
		#else
			lis_sort_jad(is,ie,maxnzr,iw,perm);
			lis_sort_jad(is,ie,maxnzr2,iw2,perm2);
			
			j = 0;
			for(i=ie-1;i>=is;i--)
			{
				for(;j<iw[i];j++)
				{
					ptr[my_rank*(maxnzr+1) + j+1] = i-is+1;
				}
				if( iw[i]==maxnzr ) break;
			}
			j = 0;
			for(i=ie-1;i>=is;i--)
			{
				for(;j<iw2[i];j++)
				{
					ptr2[my_rank*(maxnzr2+1) + j+1] = i-is+1;
				}
				if( iw2[i]==maxnzr2 ) break;
			}
		#endif

		ptr[my_rank*(maxnzr+1)] = nnzpe[my_rank];
		ptr2[my_rank*(maxnzr2+1)] = nnzpe2[my_rank];
		for(j=0;j<maxnzr;j++)
		{
			ptr[my_rank*(maxnzr+1) + j+1] += ptr[my_rank*(maxnzr+1) + j];
		}
		for(j=0;j<maxnzr2;j++)
		{
			ptr2[my_rank*(maxnzr2+1) + j+1] += ptr2[my_rank*(maxnzr2+1) + j];
		}
		#if 0
			for(i=is;i<ie;i++)
			{
				kk  = 0;
				js = Ain->L->ptr[perm[i]];
				je = Ain->L->ptr[perm[i]+1];
				for(j=js;j<je;j++)
				{
					l          = ptr[my_rank*(maxnzr+1) + kk]+i-is;
					value[l]   = Ain->L->value[j];
					index[l]   = Ain->L->index[j];
					kk++;
				}
				kk  = 0;
				js = Ain->U->ptr[perm2[i]];
				je = Ain->U->ptr[perm2[i]+1];
				for(j=js;j<je;j++)
				{
					l          = ptr2[my_rank*(maxnzr2+1) + kk]+i-is;
					value2[l]   = Ain->U->value[j];
					index2[l]   = Ain->U->index[j];
					kk++;
				}
			}
		#else
			for(j=0;j<maxnzr;j++)
			{
				js = ptr[my_rank*(maxnzr+1) + j];
				je = ptr[my_rank*(maxnzr+1) + j+1];
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=js;i<je;i++)
				{
					l  = Ain->L->ptr[perm[i-js]] + j;
					value[i]   = Ain->L->value[l];
					index[i]   = Ain->L->index[l];
				}
			}
			for(j=0;j<maxnzr2;j++)
			{
				js = ptr2[my_rank*(maxnzr2+1) + j];
				je = ptr2[my_rank*(maxnzr2+1) + j+1];
				#pragma cdir nodep
				#pragma _NEC ivdep
				for(i=js;i<je;i++)
				{
					l  = Ain->U->ptr[perm2[i-js]] + j;
					value2[i]   = Ain->U->value[l];
					index2[i]   = Ain->U->index[l];
				}
			}
		#endif
	}

	err = lis_matrix_set_jad(nnz,maxnzr,perm,ptr,index,value,Aout);
	if( err )
	{
		lis_free2(7,perm,ptr,index,value,iw,maxnzrpe,nnzpe);
		return err;
	}
	Aout->U->maxnzr = maxnzr2;
	Aout->U->row    = perm2;
	Aout->U->ptr    = ptr2;
	Aout->U->index  = index2;
	Aout->U->value  = value2;
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_free2(2,iw,nnzpe);
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	Aout->work = (LIS_SCALAR *)lis_malloc(np*sizeof(LIS_SCALAR));
	lis_free2(6,iw,nnzpe,maxnzrpe,iw2,nnzpe2,maxnzrpe2);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#endif

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_jad2csr"
LIS_INT lis_matrix_convert_jad2csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,jj,k,is,ie;
	LIS_INT err;
	LIS_INT n,nnz,maxnzr;
	LIS_INT *iw;
	LIS_INT *ptr,*index;
	LIS_SCALAR *value;
	#ifdef _OPENMP
		LIS_INT	nprocs,my_rank;
	#endif

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	nnz     = Ain->nnz;
	maxnzr  = Ain->maxnzr;
	is      = Ain->is;
	ie      = Ain->ie;

	ptr     = NULL;
	index   = NULL;
	value   = NULL;
	iw      = NULL;

	iw = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_convert_jad2csr::iw" );
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	err = lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);
	if( err )
	{
		lis_free2(4,ptr,index,value,iw);
		return err;
	}

	/* convert csr */
	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
		ptr[0] = 0;
		#pragma omp parallel private(i,j,k,is,ie,jj,my_rank)
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

			for(i=is;i<ie;i++)
			{
				ptr[i+1] = 0;
			}
			for(j=0;j<maxnzr;j++)
			{
				k = is;
				for(i=Ain->ptr[my_rank*(maxnzr+1) + j];i<Ain->ptr[my_rank*(maxnzr+1) + j+1];i++)
				{ 
					ptr[Ain->row[k]+1]++;
					k++;
				}
			}
			#pragma omp barrier
			#pragma omp single
			for(i=0;i<n;i++)
			{
				ptr[i+1] += ptr[i];
			}
			for(i=is;i<ie;i++)
			{
				iw[i]    = ptr[i];
			}
			for(j=0;j<maxnzr;j++)
			{
				jj = is;
				for(i=Ain->ptr[my_rank*(maxnzr+1) + j];i<Ain->ptr[my_rank*(maxnzr+1) + j+1];i++)
				{
					k          = iw[Ain->row[jj]]++;
					value[k]   = Ain->value[i];
					index[k]   = Ain->index[i];
					jj++;
				}
			}
		}
	#else
		for(i=0;i<n+1;i++)
		{
			ptr[i] = 0;
		}
		for(j=0;j<maxnzr;j++)
		{
			k = 0;
			for(i=Ain->ptr[j];i<Ain->ptr[j+1];i++)
			{ 
				ptr[Ain->row[k]+1]++;
				k++;
			}
		}
		for(i=0;i<n;i++)
		{
			ptr[i+1] += ptr[i];
		}
		for(i=0;i<n;i++)
		{
			iw[i]    = ptr[i];
		}
		for(j=0;j<maxnzr;j++)
		{
			jj = 0;
			for(i=Ain->ptr[j];i<Ain->ptr[j+1];i++)
			{
				k          = iw[Ain->row[jj]]++;
				value[k]   = Ain->value[i];
				index[k]   = Ain->index[i];
				jj++;
			}
		}
	#endif

	err = lis_matrix_set_csr(nnz,ptr,index,value,Aout);
	if( err )
	{
		lis_free2(4,ptr,index,value,iw);
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
#define __FUNC__ "lis_vector_sort_jad_order"
LIS_INT lis_vector_sort_jad_order(LIS_MATRIX A, LIS_VECTOR v)
{
	LIS_INT i,n,np;
	LIS_SCALAR *t;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	np = A->np;

	t  = (LIS_SCALAR *)lis_malloc(np*sizeof(LIS_SCALAR),"lis_vector_sort_jad_order::t");
	if( t==NULL )
	{
		LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0;i<n;i++)
	{
		t[i] = v->value[A->row[i]];
	}
	lis_free(v->value);
	v->value = t;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_sort_ascending_order"
LIS_INT lis_vector_sort_ascending_order(LIS_MATRIX A, LIS_VECTOR v)
{
	LIS_INT i,n,np;
	LIS_SCALAR *t;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	np = A->np;

	t  = (LIS_SCALAR *)lis_malloc(np*sizeof(LIS_SCALAR),"lis_vector_sort_ascending_order::t");
	if( t==NULL )
	{
		LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0;i<n;i++)
	{
		t[A->row[i]] = v->value[i];
	}
	lis_free(v->value);
	v->value = t;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

