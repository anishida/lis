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
#define __FUNC__ "lis_matrix_set_ell"
LIS_INT lis_matrix_set_ell(LIS_INT maxnzr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A)
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
	A->status      = -LIS_MATRIX_ELL;
	A->maxnzr      = maxnzr;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_setDLU_ell"
LIS_INT lis_matrix_setDLU_ell(LIS_INT lmaxnzr, LIS_INT umaxnzr, LIS_SCALAR *diag, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A)
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

	A->L = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_ell::A->L");
	if( A->L==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	A->U = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_ell::A->U");
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
	A->L->maxnzr   = lmaxnzr;
	A->L->index    = lindex;
	A->L->value    = lvalue;
	A->U->maxnzr   = umaxnzr;
	A->U->index    = uindex;
	A->U->value    = uvalue;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_ELL;
	A->is_splited  = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_ell"
LIS_INT lis_matrix_malloc_ell(LIS_INT n, LIS_INT maxnzr, LIS_INT **index, LIS_SCALAR **value)
{
	LIS_DEBUG_FUNC_IN;

	*index   = NULL;
	*value   = NULL;

	*index = (LIS_INT *)lis_malloc( n*maxnzr*sizeof(LIS_INT),"lis_matrix_malloc_ell::index" );
	if( *index==NULL )
	{
		LIS_SETERR_MEM(n*maxnzr*sizeof(LIS_INT));
		lis_free2(2,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*value = (LIS_SCALAR *)lis_malloc( n*maxnzr*sizeof(LIS_SCALAR),"lis_matrix_malloc_ell::value" );
	if( *value==NULL )
	{
		LIS_SETERR_MEM(n*maxnzr*sizeof(LIS_SCALAR));
		lis_free2(2,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_elements_copy_ell"
LIS_INT lis_matrix_elements_copy_ell(LIS_INT n, LIS_INT maxnzr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_index, LIS_SCALAR *o_value)
{
	LIS_INT	i,j;
	#ifdef _OPENMP
		LIS_INT	is,ie,my_rank,nprocs;
	#endif

	LIS_DEBUG_FUNC_IN;

	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
		#pragma omp parallel private(i,j,is,ie,my_rank)
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			for(j=0;j<maxnzr;j++)
			{
				for(i=is;i<ie;i++)
				{
					o_value[j*n+i] = value[j*n+i];
					o_index[j*n+i] = index[j*n+i];
				}
			}
		}
	#else
		for(j=0;j<maxnzr;j++)
		{
			for(i=0;i<n;i++)
			{
				o_value[j*n+i] = value[j*n+i];
				o_index[j*n+i] = index[j*n+i];
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_ell"
LIS_INT lis_matrix_copy_ell(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT i,n,maxnzr,lmaxnzr,umaxnzr;
	LIS_INT *index;
	LIS_INT *lindex;
	LIS_INT *uindex;
	LIS_SCALAR *value,*lvalue,*uvalue,*diag;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;

	if( Ain->is_splited )
	{
		lmaxnzr  = Ain->L->maxnzr;
		umaxnzr  = Ain->U->maxnzr;
		lindex   = NULL;
		uindex   = NULL;
		diag     = NULL;

		err = lis_matrix_malloc_ell(n,lmaxnzr,&lindex,&lvalue);
		if( err )
		{
			return err;
		}
		err = lis_matrix_malloc_ell(n,umaxnzr,&uindex,&uvalue);
		if( err )
		{
			lis_free2(5,diag,uindex,lindex,uvalue,lvalue);
			return err;
		}
		diag = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_matrix_copy_ell::diag");
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
		lis_matrix_elements_copy_ell(n,lmaxnzr,Ain->L->index,Ain->L->value,lindex,lvalue);
		lis_matrix_elements_copy_ell(n,umaxnzr,Ain->U->index,Ain->U->value,uindex,uvalue);

		err = lis_matrix_setDLU_ell(lmaxnzr,umaxnzr,diag,lindex,lvalue,uindex,uvalue,Aout);
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
		maxnzr  = Ain->maxnzr;

		err = lis_matrix_malloc_ell(n,maxnzr,&index,&value);
		if( err )
		{
			return err;
		}

		lis_matrix_elements_copy_ell(n,maxnzr,Ain->index,Ain->value,index,value);

		err = lis_matrix_set_ell(maxnzr,index,value,Aout);
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
#define __FUNC__ "lis_matrix_split_ell"
LIS_INT lis_matrix_split_ell(LIS_MATRIX A)
{
	LIS_INT i,j,n,maxnzr,my_rank,nprocs,is,ie;
	LIS_INT lmaxnzr,umaxnzr,lcount,ucount;
	LIS_INT err;
	#ifdef _OPENMP
		LIS_INT *iw;
	#endif
	LIS_INT *lindex,*uindex;
	LIS_SCALAR *lvalue,*uvalue;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;

	n        = A->n;
	maxnzr   = A->maxnzr;
	lmaxnzr  = 0;
	umaxnzr  = 0;
	D        = NULL;
	lindex   = NULL;
	lvalue   = NULL;
	uindex   = NULL;
	uvalue   = NULL;

	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
	#else
		nprocs  = 1;
	#endif
	#ifdef _OPENMP
		iw = (LIS_INT *)lis_malloc(nprocs*LIS_VEC_TMP_PADD*sizeof(LIS_INT),"lis_matrix_split_ell::iw");
		if( iw==NULL )
		{
			LIS_SETERR_MEM(nprocs*LIS_VEC_TMP_PADD*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		#pragma omp parallel private(i,j,is,ie,lcount,ucount,my_rank)
		{
			my_rank = omp_get_thread_num();
			iw[my_rank*LIS_VEC_TMP_PADD]   = 0;
			iw[my_rank*LIS_VEC_TMP_PADD+1] = 0;
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			for(i=is;i<ie;i++)
			{
				lcount = 0;
				ucount = 0;
				for(j=0;j<maxnzr;j++)
				{
					if( A->index[j*n+i]<i )
					{
						lcount++;
					}
					else if( A->index[j*n+i]>i )
					{
						ucount++;
					}
				}
				if( lcount>iw[my_rank*LIS_VEC_TMP_PADD]   ) iw[my_rank*LIS_VEC_TMP_PADD]   = lcount;
				if( ucount>iw[my_rank*LIS_VEC_TMP_PADD+1] ) iw[my_rank*LIS_VEC_TMP_PADD+1] = ucount;
			}
		}
		for(i=0;i<nprocs;i++)
		{
			if( iw[i*LIS_VEC_TMP_PADD]>lmaxnzr   ) lmaxnzr = iw[i*LIS_VEC_TMP_PADD];
			if( iw[i*LIS_VEC_TMP_PADD+1]>umaxnzr ) umaxnzr = iw[i*LIS_VEC_TMP_PADD+1];
		}
		lis_free(iw);
	#else
		for(i=0;i<n;i++)
		{
			lcount = 0;
			ucount = 0;
			for(j=0;j<maxnzr;j++)
			{
				if( A->index[j*n+i]<i )
				{
					lcount++;
				}
				else if( A->index[j*n+i]>i )
				{
					ucount++;
				}
			}
			if( lcount>lmaxnzr ) lmaxnzr = lcount;
			if( ucount>umaxnzr ) umaxnzr = ucount;
		}
	#endif

	err = lis_matrix_LU_create(A);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_ell(n,lmaxnzr,&lindex,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_ell(n,umaxnzr,&uindex,&uvalue);
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
	#pragma omp parallel private(i,j,is,ie,lcount,ucount,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
		for(j=0;j<lmaxnzr;j++)
		{
			for(i=is;i<ie;i++)
			{
				lvalue[j*n + i]   = 0.0;
				lindex[j*n + i]   = i;
				D->value[i]       = 0.0;
			}
		}
		for(j=0;j<umaxnzr;j++)
		{
			for(i=is;i<ie;i++)
			{
				uvalue[j*n + i]   = 0.0;
				uindex[j*n + i]   = i;
			}
		}
		for(i=is;i<ie;i++)
		{
			lcount = 0;
			ucount = 0;
			for(j=0;j<maxnzr;j++)
			{
				if( A->index[j*n+i]<i )
				{
					lindex[lcount*n+i] = A->index[j*n+i];
					lvalue[lcount*n+i] = A->value[j*n+i];
					lcount++;
				}
				else if( A->index[j*n+i]>i )
				{
					uindex[ucount*n+i] = A->index[j*n+i];
					uvalue[ucount*n+i] = A->value[j*n+i];
					ucount++;
				}
				else
				{
					if( A->value[j*n+i]!=0.0 ) D->value[i] = A->value[j*n+i];
				}
			}
		}
	}
	A->L->maxnzr  = lmaxnzr;
	A->L->index   = lindex;
	A->L->value   = lvalue;
	A->U->maxnzr  = umaxnzr;
	A->U->index   = uindex;
	A->U->value   = uvalue;
	A->D          = D;
	A->is_splited = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_merge_ell"
LIS_INT lis_matrix_merge_ell(LIS_MATRIX A)
{
	LIS_INT i,j,n;
	LIS_INT maxnzr,lmaxnzr,umaxnzr,count;
	LIS_INT err;
	LIS_INT *index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;


	n       = A->n;
	maxnzr  = 0;
	lmaxnzr = A->L->maxnzr;
	umaxnzr = A->U->maxnzr;
	index   = NULL;
	value   = NULL;

	for(i=0;i<n;i++)
	{
		count = 0;
		for(j=0;j<lmaxnzr;j++)
		{
			if( A->L->index[j*n+i]<i )
			{
				count++;
			}
		}
		for(j=0;j<umaxnzr;j++)
		{
			if( A->U->index[j*n+i]>i )
			{
				count++;
			}
		}
		count++;
		if( count>maxnzr ) maxnzr = count;
	}

	err = lis_matrix_malloc_ell(n,maxnzr,&index,&value);
	if( err )
	{
		return err;
	}

	for(j=0;j<maxnzr;j++)
	{
		for(i=0;i<n;i++)
		{
			value[j*n + i]   = 0.0;
			index[j*n + i]   = i;
		}
	}
	for(i=0;i<n;i++)
	{
		count = 0;
		for(j=0;j<lmaxnzr;j++)
		{
			if( A->L->index[j*n+i]<i )
			{
				index[count*n+i] = A->L->index[j*n+i];
				value[count*n+i] = A->L->value[j*n+i];
				count++;
			}
		}
		index[count*n+i] = i;
		value[count*n+i] = A->D->value[i];
		count++;
		for(j=0;j<umaxnzr;j++)
		{
			if( A->U->index[j*n+i]>i )
			{
				index[count*n+i] = A->U->index[j*n+i];
				value[count*n+i] = A->U->value[j*n+i];
				count++;
			}
		}
	}

	A->maxnzr     = maxnzr;
	A->value      = value;
	A->index      = index;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal_ell"
LIS_INT lis_matrix_get_diagonal_ell(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j;
	LIS_INT n,maxnzr;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			d[i] = A->D->value[i];
		}
	}
	else
	{
		maxnzr = A->maxnzr;
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			d[i] = (LIS_SCALAR)0.0;
			for(j=0;j<maxnzr;j++)
			{
				if( i==A->index[j*n+i] )
				{
					d[i] = A->value[j*n+i];
					break;
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_ell"
LIS_INT lis_matrix_shift_diagonal_ell(LIS_MATRIX A, LIS_SCALAR sigma)
{
	LIS_INT i,j;
	LIS_INT n,maxnzr;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			A->D->value[i] -= sigma;
		}
	}
	else
	{
		maxnzr = A->maxnzr;
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			for(j=0;j<maxnzr;j++)
			{
				if( i==A->index[j*n+i] )
				{
					A->value[j*n+i] -= sigma;
					break;
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_ell"
LIS_INT lis_matrix_scale_ell(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j,is,ie;
	LIS_INT n,maxnzr,nprocs,my_rank;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	if( A->is_splited )
	{
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
			for(i=is;i<ie;i++)
			{
				A->D->value[i] = 1.0;
			}
			for(j=0;j<A->L->maxnzr;j++)
			{
				for(i=is;i<ie;i++)
				{
					A->L->value[j*n + i] *= d[i];
				}
			}
			for(j=0;j<A->U->maxnzr;j++)
			{
				for(i=is;i<ie;i++)
				{
					A->U->value[j*n + i] *= d[i];
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
		#pragma omp parallel private(i,j,is,ie,my_rank)
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
				for(i=is;i<ie;i++)
				{
					A->value[j*n + i] *= d[i];
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_symm_ell"
LIS_INT lis_matrix_scale_symm_ell(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j,is,ie;
	LIS_INT n,maxnzr,nprocs,my_rank;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	if( A->is_splited )
	{
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
			for(i=is;i<ie;i++)
			{
				A->D->value[i] = 1.0;
			}
			for(j=0;j<A->L->maxnzr;j++)
			{
				for(i=is;i<ie;i++)
				{
					A->L->value[j*n + i] *= d[i]*d[A->L->index[j*n + i]];
				}
			}
			for(j=0;j<A->U->maxnzr;j++)
			{
				for(i=is;i<ie;i++)
				{
					A->U->value[j*n + i] *= d[i]*d[A->U->index[j*n + i]];
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
		#pragma omp parallel private(i,j,is,ie,my_rank)
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
				for(i=is;i<ie;i++)
				{
					A->value[j*n + i] *= d[i]*d[A->index[j*n + i]];;
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solve_ell"
LIS_INT lis_matrix_solve_ell(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,n;
	LIS_SCALAR t;
	LIS_SCALAR *b,*x;

	LIS_DEBUG_FUNC_IN;

	n       = A->n;
	b       = B->value;
	x       = X->value;

	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=0;i<n;i++)
		{
			t = b[i];
			for(j=0;j<A->L->maxnzr;j++)
			{
				t -= A->L->value[j*n + i] * x[A->L->index[j*n + i]];
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			t = b[i];
			for(j=0;j<A->U->maxnzr;j++)
			{
				t -= A->U->value[j*n + i] * x[A->U->index[j*n + i]];
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<n;i++)
		{
			t = b[i];
			for(j=0;j<A->L->maxnzr;j++)
			{
				t -= A->L->value[j*n + i] * x[A->L->index[j*n + i]];
			}
			x[i]   = t * A->WD->value[i];
		}
		for(i=n-1;i>=0;i--)
		{
			t = 0.0;
			for(j=0;j<A->U->maxnzr;j++)
			{
				if( A->U->index[j*n + i]>=n ) continue; 
				t += A->U->value[j*n + i] * x[A->U->index[j*n + i]];
			}
			x[i]  -= t * A->WD->value[i];
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solveh_ell"
LIS_INT lis_matrix_solveh_ell(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,n;
	LIS_SCALAR t;
	LIS_SCALAR *x;

	LIS_DEBUG_FUNC_IN;

	n       = A->n;
	x       = X->value;

	lis_vector_copy(B,X);
	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=0;i<n;i++)
		{
			x[i]   = x[i] * conj(A->WD->value[i]);
			for(j=0;j<A->U->maxnzr;j++)
			{
				x[A->U->index[j*n + i]] -= conj(A->U->value[j*n + i]) * x[i];
			}
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			x[i]   = x[i] * conj(A->WD->value[i]);
			for(j=0;j<A->L->maxnzr;j++)
			{
				x[A->L->index[j*n +i]] -= conj(A->L->value[j*n + i]) * x[i];
			}
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<n;i++)
		{
			t   = x[i] * conj(A->WD->value[i]);
			for(j=0;j<A->U->maxnzr;j++)
			{
				x[A->U->index[j*n + i]] -= conj(A->U->value[j*n + i]) * t;
			}
		}
		for(i=n-1;i>=0;i--)
		{
			t    = x[i] * conj(A->WD->value[i]);
			x[i] = t;
			for(j=0;j<A->L->maxnzr;j++)
			{
				x[A->L->index[j*n + i]] -= conj(A->L->value[j*n + i]) * t;
			}
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_csr2ell"
LIS_INT lis_matrix_convert_csr2ell(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,k;
	LIS_INT err;
	LIS_INT n,maxnzr,nprocs,my_rank;
	LIS_INT is,ie,count;
	LIS_INT *iw;
	LIS_INT *index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;

	index   = NULL;
	value   = NULL;
	iw      = NULL;


	/* check maxnzr */
	#ifdef _OPENMP
		#define PADD 32
		nprocs  = omp_get_max_threads();
		iw = (LIS_INT *)lis_malloc( nprocs*PADD*sizeof(LIS_INT),"lis_matrix_convert_csr2ell::iw" );
		if( iw==NULL )
		{
			LIS_SETERR_MEM(nprocs*PADD*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		#pragma omp parallel private(i,j,k,is,ie,my_rank)
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

			iw[my_rank*PADD] = 0;
			for(i=is;i<ie;i++)
			{
				k = Ain->ptr[i+1] - Ain->ptr[i];
				if( k > iw[my_rank*PADD] ) iw[my_rank*PADD] = k;
			}
		}
		maxnzr = 0;
		for(i=0;i<nprocs;i++)
		{
			if( iw[i*PADD] > maxnzr ) maxnzr = iw[i*PADD]; 
		}
		lis_free(iw);
	#else
		maxnzr  = 0;
		for(i=0;i<n;i++)
		{
			count = Ain->ptr[i+1] - Ain->ptr[i];
			if( count > maxnzr ) maxnzr = count;
		}
	#endif


	err = lis_matrix_malloc_ell(n,maxnzr,&index,&value);
	if( err )
	{
		return err;
	}

	/* convert ell */
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k,is,ie,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
			nprocs  = omp_get_max_threads();
		#else
			my_rank = 0;
			nprocs  = 1;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		for(j=0;j<maxnzr;j++)
		{
			for(i=is;i<ie;i++)
			{
				value[j*n + i]   = 0.0;
				index[j*n + i]   = i;
			}
		}
		for(i=is;i<ie;i++)
		{
			k=0;
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				value[k*n + i]   = Ain->value[j];
				index[k*n + i]   = Ain->index[j];
				k++;
			}
		}
	}

	err = lis_matrix_set_ell(maxnzr,index,value,Aout);
	if( err )
	{
		lis_free2(2,index,value);
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
#define __FUNC__ "lis_matrix_convert_ell2csr"
LIS_INT lis_matrix_convert_ell2csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,k;
	LIS_INT err;
	LIS_INT n,nnz,maxnzr,is,ie,nprocs,my_rank;
	LIS_INT *iw;
	LIS_INT *ptr,*index;
	LIS_SCALAR *value;

	n       = Ain->n;
	maxnzr  = Ain->maxnzr;
	is      = Ain->is;
	ie      = Ain->ie;

	ptr     = NULL;
	index   = NULL;
	value   = NULL;
	iw      = NULL;

	iw = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_convert_ell2csr::iw" );
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	ptr = (LIS_INT *)lis_malloc( (n+1)*sizeof(LIS_INT),"lis_matrix_convert_ell2csr::ptr" );
	if( ptr==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
		lis_free2(4,ptr,index,value,iw);
		return LIS_OUT_OF_MEMORY;
	}

	/* check nnz */
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k,is,ie,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
			nprocs  = omp_get_max_threads();
		#else
			my_rank = 0;
			nprocs  = 1;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
		memset(&iw[is],0,(ie-is)*sizeof(LIS_INT));

		for(j=0;j<maxnzr;j++)
		{
			for(i=is;i<ie;i++)
			{ 
				if( Ain->value[j*n + i]!=(LIS_SCALAR)0.0 )
				{
					iw[i]++;
				}
			}
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n+1;i++)
		{
			ptr[i] = 0;
		}
		#ifdef _OPENMP
		#pragma omp single
		#endif
		for(i=0;i<n;i++)
		{
			ptr[i+1] = ptr[i] + iw[i];
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n;i++)
		{
			iw[i]    = ptr[i];
		}
	}
	nnz = ptr[n];

	index = (LIS_INT *)lis_malloc( nnz*sizeof(LIS_INT),"lis_matrix_convert_ell2csr::index" );
	if( index==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		lis_free2(4,ptr,index,value,iw);
		return LIS_OUT_OF_MEMORY;
	}
	value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_convert_ell2csr::value" );
	if( value==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		lis_free2(4,ptr,index,value,iw);
		return LIS_OUT_OF_MEMORY;
	}

	/* convert csr */
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k,is,ie,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
			nprocs  = omp_get_max_threads();
		#else
			my_rank = 0;
			nprocs  = 1;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		for(j=0;j<maxnzr;j++)
		{
			for(i=is;i<ie;i++)
			{ 
				if( Ain->value[j*n + i]!=(LIS_SCALAR)0.0 )
				{
					k        = iw[i]++;
					value[k] = Ain->value[j*n + i];
					index[k] = Ain->index[j*n + i];
				}
			}
		}
	}

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

