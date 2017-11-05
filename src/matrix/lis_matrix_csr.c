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
#define __FUNC__ "lis_matrix_set_csr"
LIS_INT lis_matrix_set_csr(LIS_INT nnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A)
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

	A->ptr         = ptr;
	A->index       = index;
	A->value       = value;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_CSR;
	A->nnz         = nnz;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_setDLU_csr"
LIS_INT lis_matrix_setDLU_csr(LIS_INT nnzl, LIS_INT nnzu, LIS_SCALAR *diag, LIS_INT *lptr, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uptr, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A)
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

	A->L = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT), "lis_matrix_setDLU_csr::A->L");
	if( A->L==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	A->U = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT), "lis_matrix_setDLU_csr::A->U");
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
	A->L->nnz      = nnzl;
	A->L->ptr      = lptr;
	A->L->index    = lindex;
	A->L->value    = lvalue;
	A->U->nnz      = nnzu;
	A->U->ptr      = uptr;
	A->U->index    = uindex;
	A->U->value    = uvalue;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_CSR;
	A->is_splited  = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_csr"
LIS_INT lis_matrix_malloc_csr(LIS_INT n, LIS_INT nnz, LIS_INT **ptr, LIS_INT **index, LIS_SCALAR **value)
{
	LIS_DEBUG_FUNC_IN;

	*ptr     = NULL;
	*index   = NULL;
	*value   = NULL;

	*ptr = (LIS_INT *)lis_malloc( (n+1)*sizeof(LIS_INT),"lis_matrix_malloc_csr::ptr" );
	if( *ptr==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
		lis_free2(3,*ptr,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*index = (LIS_INT *)lis_malloc( nnz*sizeof(LIS_INT),"lis_matrix_malloc_csr::index" );
	if( *index==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		lis_free2(3,*ptr,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_malloc_csr::value" );
	if( *value==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		lis_free2(3,*ptr,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_psd_set_value_csr"
LIS_INT lis_matrix_psd_set_value_csr(LIS_INT flag, LIS_INT i, LIS_INT j, LIS_SCALAR value, LIS_MATRIX A)
{
    LIS_INT n,gn,is,ie,k,k0,k1,j_global_test,j_global;
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

    n   = A->n;
	gn  = A->gn;
    is  = A->is;
    ie  = A->ie;

	if( A->origin )
	{
		i--;
		j--;
	}

	if( i<0  || j<0 )
	{
		k = 0;
		if( A->origin )
		{
			i++;
			j++;
			k++;
		}
		LIS_SETERR3(LIS_ERR_ILL_ARG,"i(=%D) or j(=%D) are less than %D\n",i,j,k);
		return LIS_ERR_ILL_ARG;
	}
	if( i>=gn || j>=gn )
	{
		if( A->origin )
		{
			i++;
			j++;
		}
		LIS_SETERR3(LIS_ERR_ILL_ARG,"i(=%D) or j(=%D) are larger than global n=(%D)\n",i,j,gn);
		return LIS_ERR_ILL_ARG;
	}


    k0 = A->ptr[i-is];
    k1 = A->ptr[i-is+1];

    for(k=k0;k<k1;k++)
    {
        j_global_test=A->index[k]+is;
        if ( j_global_test>=is && j_global_test<ie ) {
            j_global=j_global_test;
        } else {
            j_global=A->l2g_map[A->index[k]-n];
        }

        if (j_global==j) {
            if(flag==LIS_INS_VALUE)
            {
                A->value[k]=value;
            }
            else 
            {
                A->value[k]+=value;
            }
            break;
        }
    }

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_elements_copy_csr"
LIS_INT lis_matrix_elements_copy_csr(LIS_INT n, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_ptr, LIS_INT *o_index, LIS_SCALAR *o_value)
{
	LIS_INT	i,j;

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
			o_ptr[i] = ptr[i];
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n;i++)
		{
			for(j=ptr[i];j<ptr[i+1];j++)
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
#define __FUNC__ "lis_matrix_copy_csr"
LIS_INT lis_matrix_copy_csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT i,n,nnz,lnnz,unnz;
	LIS_INT *ptr,*index;
	LIS_INT *lptr,*lindex;
	LIS_INT *uptr,*uindex;
	LIS_SCALAR *value,*lvalue,*uvalue,*diag;

	LIS_DEBUG_FUNC_IN;


	n       = Ain->n;

	if( Ain->is_splited )
	{
		lnnz     = Ain->L->nnz;
		unnz     = Ain->U->nnz;
		lptr     = NULL;
		lindex   = NULL;
		uptr     = NULL;
		uindex   = NULL;
		diag     = NULL;

		err = lis_matrix_malloc_csr(n,lnnz,&lptr,&lindex,&lvalue);
		if( err )
		{
			return err;
		}
		err = lis_matrix_malloc_csr(n,unnz,&uptr,&uindex,&uvalue);
		if( err )
		{
			lis_free2(7,diag,uptr,lptr,uindex,lindex,uvalue,lvalue);
			return err;
		}
		diag = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_matrix_copy_csr::diag");
		if( diag==NULL )
		{
			lis_free2(7,diag,uptr,lptr,uindex,lindex,uvalue,lvalue);
			return err;
		}

		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<n;i++)
		{
			diag[i] = Ain->D->value[i];
		}
		lis_matrix_elements_copy_csr(n,Ain->L->ptr,Ain->L->index,Ain->L->value,lptr,lindex,lvalue);
		lis_matrix_elements_copy_csr(n,Ain->U->ptr,Ain->U->index,Ain->U->value,uptr,uindex,uvalue);

		err = lis_matrix_setDLU_csr(lnnz,unnz,diag,lptr,lindex,lvalue,uptr,uindex,uvalue,Aout);
		if( err )
		{
			lis_free2(7,diag,uptr,lptr,uindex,lindex,uvalue,lvalue);
			return err;
		}
	}
	if( !Ain->is_splited || (Ain->is_splited && Ain->is_save) )
	{
		ptr     = NULL;
		index   = NULL;
		value   = NULL;
		nnz     = Ain->nnz;
		err = lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);
		if( err )
		{
			return err;
		}

		lis_matrix_elements_copy_csr(n,Ain->ptr,Ain->index,Ain->value,ptr,index,value);

		err = lis_matrix_set_csr(nnz,ptr,index,value,Aout);
		if( err )
		{
			lis_free2(3,ptr,index,value);
			return err;
		}
	}
	if( Ain->matrix_type==LIS_MATRIX_CSC )
	{
		Aout->matrix_type = LIS_MATRIX_CSC;
		Aout->status = -LIS_MATRIX_CSC;
		err = lis_matrix_assemble(Aout);
	}
	else
	{
		err = lis_matrix_assemble(Aout);
	}
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copyDLU_csr"
LIS_INT lis_matrix_copyDLU_csr(LIS_MATRIX Ain, LIS_MATRIX_DIAG *D, LIS_MATRIX *L, LIS_MATRIX *U)
{
	LIS_INT err;
	LIS_INT i,n,np,lnnz,unnz;
	LIS_INT *lptr,*lindex;
	LIS_INT *uptr,*uindex;
	LIS_SCALAR *lvalue,*uvalue,*diag;

	LIS_DEBUG_FUNC_IN;
	
	*D = NULL;
	*L = NULL;
	*U = NULL;

	err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	n       = Ain->n;
	np      = Ain->np;

	err = lis_matrix_duplicate(Ain,L);
	if( err )
	{
		return err;
	}
	err = lis_matrix_duplicate(Ain,U);
	if( err )
	{
		lis_matrix_destroy(*L);
		return err;
	}
	err = lis_matrix_diag_duplicateM(Ain,D);
	if( err )
	{
		lis_matrix_destroy(*L);
		lis_matrix_destroy(*U);
		return err;
	}
	lis_free((*D)->value);

	if( Ain->is_splited )
	{
	}
	lnnz     = Ain->L->nnz;
	unnz     = Ain->U->nnz;
	lptr     = NULL;
	lindex   = NULL;
	uptr     = NULL;
	uindex   = NULL;
	diag     = NULL;

	err = lis_matrix_malloc_csr(n,lnnz,&lptr,&lindex,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_csr(n,unnz,&uptr,&uindex,&uvalue);
	if( err )
	{
		lis_free2(7,diag,uptr,lptr,uindex,lindex,uvalue,lvalue);
		return err;
	}
	diag = (LIS_SCALAR *)lis_malloc(np*sizeof(LIS_SCALAR),"lis_matrix_copyDLU_csr::diag");
	if( diag==NULL )
	{
		lis_free2(7,diag,uptr,lptr,uindex,lindex,uvalue,lvalue);
		return err;
	}

	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0;i<n;i++)
	{
		diag[i] = Ain->D->value[i];
	}
	lis_matrix_elements_copy_csr(n,Ain->L->ptr,Ain->L->index,Ain->L->value,lptr,lindex,lvalue);
	lis_matrix_elements_copy_csr(n,Ain->U->ptr,Ain->U->index,Ain->U->value,uptr,uindex,uvalue);

	(*D)->value = diag;
	err = lis_matrix_set_csr(lnnz,lptr,lindex,lvalue,*L);
	if( err )
	{
		lis_free2(7,diag,uptr,lptr,uindex,lindex,uvalue,lvalue);
		return err;
	}
	err = lis_matrix_set_csr(unnz,uptr,uindex,uvalue,*U);
	if( err )
	{
		lis_free2(7,diag,uptr,lptr,uindex,lindex,uvalue,lvalue);
		return err;
	}

	err = lis_matrix_assemble(*L);
	if( err )
	{
		return err;
	}
	err = lis_matrix_assemble(*U);
	if( err )
	{
		return err;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal_csr"
LIS_INT lis_matrix_get_diagonal_csr(LIS_MATRIX A, LIS_SCALAR d[])
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
			d[i] = A->D->value[i];
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			d[i] = (LIS_SCALAR)0.0;
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( i==A->index[j] )
				{
					d[i] = A->value[j];
					break;
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_csr"
LIS_INT lis_matrix_shift_diagonal_csr(LIS_MATRIX A, LIS_SCALAR sigma)
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
			A->D->value[i] -= sigma;
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( i==A->index[j] )
				{
					A->value[j] -= sigma;
					break;
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_csr"
LIS_INT lis_matrix_scale_csr(LIS_MATRIX A, LIS_SCALAR d[])
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
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				A->L->value[j] *= d[i];
			}
			for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
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
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				A->value[j] *= d[i];
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_symm_csr"
LIS_INT lis_matrix_scale_symm_csr(LIS_MATRIX A, LIS_SCALAR d[])
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
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				A->L->value[j] = A->L->value[j]*d[i]*d[A->L->index[j]];
			}
			for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
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
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				A->value[j] = A->value[j]*d[i]*d[A->index[j]];
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_normf_csr"
LIS_INT lis_matrix_normf_csr(LIS_MATRIX A, LIS_SCALAR *nrm)
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
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				sum += A->L->value[j]*A->L->value[j];
			}
			for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
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
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
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
#define __FUNC__ "lis_matrix_transpose_csr"
LIS_INT lis_matrix_transpose_csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_convert_csr2csc(Ain,Aout);
	if( err )
	{
		return err;
	}
	Aout->matrix_type = LIS_MATRIX_CSR;
	Aout->status      = LIS_MATRIX_CSR;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_csr"
LIS_INT lis_matrix_split_csr(LIS_MATRIX A)
{
	LIS_INT i,j,n;
	LIS_INT nnzl,nnzu;
	LIS_INT err;
	LIS_INT *lptr,*lindex,*uptr,*uindex;
	LIS_SCALAR *lvalue,*uvalue;
	LIS_MATRIX_DIAG	D;
	#ifdef _OPENMP
		LIS_INT kl,ku;
		LIS_INT *liw,*uiw;
	#endif

	LIS_DEBUG_FUNC_IN;

	n        = A->n;
	nnzl     = 0;
	nnzu     = 0;
	D        = NULL;
	lptr     = NULL;
	lindex   = NULL;
	lvalue   = NULL;
	uptr     = NULL;
	uindex   = NULL;
	uvalue   = NULL;

	#ifdef _OPENMP
		liw = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_split_csr::liw");
		if( liw==NULL )
		{
			LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		uiw = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_split_csr::uiw");
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
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
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
		for(i=0;i<n;i++)
		{
			liw[i+1] += liw[i];
			uiw[i+1] += uiw[i];
		}
		nnzl = liw[n];
		nnzu = uiw[n];
	#else
		for(i=0;i<n;i++)
		{
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<i )
				{
					nnzl++;
				}
				else if( A->index[j]>i )
				{
					nnzu++;
				}
			}
		}
	#endif

	err = lis_matrix_LU_create(A);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_csr(n,nnzl,&lptr,&lindex,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_csr(n,nnzu,&uptr,&uindex,&uvalue);
	if( err )
	{
		lis_free2(6,lptr,lindex,lvalue,uptr,uindex,uvalue);
		return err;
	}
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		lis_free2(6,lptr,lindex,lvalue,uptr,uindex,uvalue);
		return err;
	}

	#ifdef _OPENMP
		#pragma omp parallel for private(i)
		for(i=0;i<n+1;i++)
		{
			lptr[i] = liw[i];
			uptr[i] = uiw[i];
		}
		#pragma omp parallel for private(i,j,kl,ku)
		for(i=0;i<n;i++)
		{
			kl = lptr[i];
			ku = uptr[i];
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
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
				else
				{
					D->value[i] = A->value[j];
				}
			}
		}
		lis_free2(2,liw,uiw);
	#else
		nnzl = 0;
		nnzu = 0;
		lptr[0] = 0;
		uptr[0] = 0;
		for(i=0;i<n;i++)
		{
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<i )
				{
					lindex[nnzl]   = A->index[j];
					lvalue[nnzl]   = A->value[j];
					nnzl++;
				}
				else if( A->index[j]>i )
				{
					uindex[nnzu]   = A->index[j];
					uvalue[nnzu]   = A->value[j];
					nnzu++;
				}
				else
				{
					D->value[i] = A->value[j];
				}
			}
			lptr[i+1] = nnzl;
			uptr[i+1] = nnzu;
		}
	#endif
	A->L->nnz     = nnzl;
	A->L->ptr     = lptr;
	A->L->index   = lindex;
	A->L->value   = lvalue;
	A->U->nnz     = nnzu;
	A->U->ptr     = uptr;
	A->U->index   = uindex;
	A->U->value   = uvalue;
	A->D          = D;
	A->is_splited = LIS_TRUE;
	
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_matrix_split_create_csr"
LIS_INT lis_matrix_split_create_csr(LIS_MATRIX A)
{
	LIS_INT i,j,n;
	LIS_INT nnzl,nnzu;
	LIS_INT err;
	LIS_INT *lptr,*lindex,*uptr,*uindex;
	LIS_SCALAR *lvalue,*uvalue;
	LIS_MATRIX_DIAG	D;
	#ifdef _OPENMP
		LIS_INT kl,ku;
		LIS_INT *liw,*uiw;
	#endif

	LIS_DEBUG_FUNC_IN;

	n        = A->n;
	nnzl     = 0;
	nnzu     = 0;
	D        = NULL;
	lptr     = NULL;
	lindex   = NULL;
	lvalue   = NULL;
	uptr     = NULL;
	uindex   = NULL;
	uvalue   = NULL;

	#ifdef _OPENMP
		liw = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_split_create_csr::liw");
		if( liw==NULL )
		{
			LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		uiw = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_split_create_csr::uiw");
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
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
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
		for(i=0;i<n;i++)
		{
			liw[i+1] += liw[i];
			uiw[i+1] += uiw[i];
		}
		nnzl = liw[n];
		nnzu = uiw[n];
	#else
		for(i=0;i<n;i++)
		{
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<i )
				{
					nnzl++;
				}
				else if( A->index[j]>i )
				{
					nnzu++;
				}
			}
		}
	#endif

	err = lis_matrix_LU_create(A);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_csr(n,nnzl,&lptr,&lindex,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_csr(n,nnzu,&uptr,&uindex,&uvalue);
	if( err )
	{
		lis_free2(6,lptr,lindex,lvalue,uptr,uindex,uvalue);
		return err;
	}
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		lis_free2(6,lptr,lindex,lvalue,uptr,uindex,uvalue);
		return err;
	}

	#ifdef _OPENMP
		#pragma omp parallel for private(i)
		for(i=0;i<n+1;i++)
		{
			lptr[i] = liw[i];
			uptr[i] = uiw[i];
		}
		#pragma omp parallel for private(i,j,kl,ku)
		for(i=0;i<n;i++)
		{
			kl = lptr[i];
			ku = uptr[i];
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<i )
				{
					lindex[kl]   = A->index[j];
					kl++;
				}
				else if( A->index[j]>i )
				{
					uindex[ku]   = A->index[j];
					ku++;
				}
			}
		}
		lis_free2(2,liw,uiw);
	#else
		nnzl = 0;
		nnzu = 0;
		lptr[0] = 0;
		uptr[0] = 0;
		for(i=0;i<n;i++)
		{
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<i )
				{
					lindex[nnzl]   = A->index[j];
					nnzl++;
				}
				else if( A->index[j]>i )
				{
					uindex[nnzu]   = A->index[j];
					nnzu++;
				}
			}
			lptr[i+1] = nnzl;
			uptr[i+1] = nnzu;
		}
	#endif
	A->L->nnz     = nnzl;
	A->L->ptr     = lptr;
	A->L->index   = lindex;
	A->L->value   = lvalue;
	A->U->nnz     = nnzu;
	A->U->ptr     = uptr;
	A->U->index   = uindex;
	A->U->value   = uvalue;
	A->D          = D;
	A->is_splited = LIS_TRUE;
	
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_matrix_split_update_csr"
LIS_INT lis_matrix_split_update_csr(LIS_MATRIX A)
{
	LIS_INT i,j,n;
	LIS_INT nnzl,nnzu;
	LIS_INT err;
/*    LIS_INT *lptr,*uptr;*/
	#ifdef _OPENMP
		LIS_INT kl,ku;
		LIS_INT *liw,*uiw;
	#endif

	LIS_DEBUG_FUNC_IN;

	n        = A->n;
	nnzl     = 0;
	nnzu     = 0;
/*    lptr     = NULL;*/
/*    uptr     = NULL;*/

    #ifdef _OPENMP
        liw = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_split_update_csr::liw");
        if( liw==NULL )
        {
            LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
            return LIS_OUT_OF_MEMORY;
        }
        uiw = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_split_update_csr::uiw");
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
            for(j=A->ptr[i];j<A->ptr[i+1];j++)
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
        for(i=0;i<n;i++)
        {
            liw[i+1] += liw[i];
            uiw[i+1] += uiw[i];
        }
    #endif

	#ifdef _OPENMP
        /*NEH I have made the changes here that I think will be necessary, */
        /*NEH but as we don't use threads, this hasn't actually been checked at all . . . */
/*	  #pragma omp parallel for private(i)*/
/*        for(i=0;i<n+1;i++)*/
/*        {*/
/*            lptr[i] = liw[i];*/
/*            uptr[i] = uiw[i];*/
/*        }*/
		#pragma omp parallel for private(i,j,kl,ku)
		for(i=0;i<n;i++)
		{
/*            kl = lptr[i];*/
/*            ku = uptr[i];*/
			kl = A->L->ptr[i];
			ku = A->U->ptr[i];
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<i )
				{
					A->L->value[kl]   = A->value[j];
					kl++;
				}
				else if( A->index[j]>i )
				{
					A->U->value[ku]   = A->value[j];
					ku++;
				}
				else
				{
					A->D->value[i] = A->value[j];
				}
			}
		}
		lis_free2(2,liw,uiw);
	#else
		nnzl = 0;
		nnzu = 0;
		for(i=0;i<n;i++)
		{
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<i )
				{
					A->L->value[nnzl]   = A->value[j];
					nnzl++;
				}
				else if( A->index[j]>i )
				{
					A->U->value[nnzu]   = A->value[j];
					nnzu++;
				}
				else
				{
					A->D->value[i] = A->value[j];
				}
			}
		}
	#endif
	
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_merge_csr"
LIS_INT lis_matrix_merge_csr(LIS_MATRIX A)
{
	LIS_INT i,j,n;
	LIS_INT nnz;
	LIS_INT err;
	LIS_INT *ptr,*index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;


	n       = A->n;
	nnz     = 0;
	ptr     = NULL;
	index   = NULL;
	value   = NULL;
	nnz     = A->L->nnz + A->U->nnz + n;

	err = lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);
	if( err )
	{
		return err;
	}

	nnz    = 0;
	ptr[0] = 0;
	for(i=0;i<n;i++)
	{
		for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
		{
			index[nnz]   = A->L->index[j];
			value[nnz]   = A->L->value[j];
			nnz++;
		}
		index[nnz]   = i;
		value[nnz]   = A->D->value[i];
		nnz++;
		for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
		{
			index[nnz]   = A->U->index[j];
			value[nnz]   = A->U->value[j];
			nnz++;
		}
		ptr[i+1] = nnz;
	}

	A->nnz        = nnz;
	A->ptr        = ptr;
	A->value      = value;
	A->index      = index;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split2_csr"
LIS_INT lis_matrix_split2_csr(LIS_MATRIX A)
{
	LIS_INT i,j,n;
	LIS_INT nnzl,nnzu;
	LIS_INT err;
	LIS_INT *lptr,*lindex,*uptr,*uindex;
	LIS_SCALAR *lvalue,*uvalue;
	#ifdef _OPENMP
		LIS_INT kl,ku;
		LIS_INT *liw,*uiw;
	#endif

	LIS_DEBUG_FUNC_IN;

	n        = A->n;
	nnzl     = 0;
	nnzu     = 0;
	lptr     = NULL;
	lindex   = NULL;
	lvalue   = NULL;
	uptr     = NULL;
	uindex   = NULL;
	uvalue   = NULL;

	#ifdef _OPENMP
		liw = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_split2_csr::liw");
		if( liw==NULL )
		{
			LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		uiw = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_split2_csr::uiw");
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
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<n )
				{
					liw[i+1]++;
				}
				else
				{
					uiw[i+1]++;
				}
			}
		}
		for(i=0;i<n;i++)
		{
			liw[i+1] += liw[i];
			uiw[i+1] += uiw[i];
		}
		nnzl = liw[n];
		nnzu = uiw[n];
	#else
		for(i=0;i<n;i++)
		{
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<n )
				{
					nnzl++;
				}
				else
				{
					nnzu++;
				}
			}
		}
	#endif

	err = lis_matrix_LU_create(A);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_csr(n,nnzl,&lptr,&lindex,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_csr(n,nnzu,&uptr,&uindex,&uvalue);
	if( err )
	{
		lis_free2(6,lptr,lindex,lvalue,uptr,uindex,uvalue);
		return err;
	}

	#ifdef _OPENMP
		#pragma omp parallel for private(i)
		for(i=0;i<n+1;i++)
		{
			lptr[i] = liw[i];
			uptr[i] = uiw[i];
		}
		#pragma omp parallel for private(i,j,kl,ku)
		for(i=0;i<n;i++)
		{
			kl = lptr[i];
			ku = uptr[i];
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<n )
				{
					lindex[kl]   = A->index[j];
					lvalue[kl]   = A->value[j];
					kl++;
				}
				else
				{
					uindex[ku]   = A->index[j];
					uvalue[ku]   = A->value[j];
					ku++;
				}
			}
		}
		lis_free2(2,liw,uiw);
	#else
		nnzl = 0;
		nnzu = 0;
		lptr[0] = 0;
		uptr[0] = 0;
		for(i=0;i<n;i++)
		{
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]<n )
				{
					lindex[nnzl]   = A->index[j];
					lvalue[nnzl]   = A->value[j];
					nnzl++;
				}
				else
				{
					uindex[nnzu]   = A->index[j];
					uvalue[nnzu]   = A->value[j];
					nnzu++;
				}
			}
			lptr[i+1] = nnzl;
			uptr[i+1] = nnzu;
		}
	#endif
	A->L->nnz     = nnzl;
	A->L->ptr     = lptr;
	A->L->index   = lindex;
	A->L->value   = lvalue;
	A->U->nnz     = nnzu;
	A->U->ptr     = uptr;
	A->U->index   = uindex;
	A->U->value   = uvalue;
	A->is_splited = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_sort_csr"
LIS_INT lis_matrix_sort_csr(LIS_MATRIX A)
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
#define __FUNC__ "lis_matrix_solve_csr"
LIS_INT lis_matrix_solve_csr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,n,jj;
	LIS_SCALAR t;
	LIS_SCALAR *b,*x;
	#ifdef _OPENMP
		LIS_INT is,ie,my_rank,nprocs;
	#endif
	#ifdef USE_QUAD_PRECISION
		LIS_QUAD w1,w2;
		LIS_SCALAR *xl;
	#endif
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	b  = B->value;
	x  = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif

	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=0;i<n;i++)
		{
			t = b[i];
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
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
			for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
			{
				t -= A->U->value[j] * x[A->U->index[j]];
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_SSOR:
	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			#ifdef _OPENMP
				nprocs = omp_get_max_threads();
				#pragma omp parallel private(i,j,jj,t,is,ie,my_rank)
				{
					my_rank = omp_get_thread_num();
					LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

					for(i=is;i<ie;i++)
					{
						t = b[i];
						for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
						{
							jj = A->L->index[j];
							if( jj<is ) continue;
							t -= A->L->value[j] * x[jj];
						}
						x[i]   = t * A->WD->value[i];
					}
					for(i=ie-1;i>=is;i--)
					{
						t = 0.0;
						for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
						{
							jj = A->U->index[j];
							if( jj<is || jj>=ie ) continue; 
							t += A->U->value[j] * x[jj];
						}
						x[i]  -= t * A->WD->value[i];
					}
				}
			#else
				for(i=0;i<n;i++)
				{
					t = b[i];
					for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
					{

						t -= A->L->value[j] * x[A->L->index[j]];
					}
					x[i]   = t * A->WD->value[i];
				}
				for(i=n-1;i>=0;i--)
				{
					t = 0.0;
					for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
					{
						if( A->U->index[j]>=n ) continue; 
						t += A->U->value[j] * x[A->U->index[j]];
					}
					x[i]  -= t * A->WD->value[i];
				}
			#endif
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			#ifdef _OPENMP
				nprocs = omp_get_max_threads();
				#ifndef USE_SSE2
					#pragma omp parallel private(i,j,jj,is,ie,w1,my_rank,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
				#else
					#pragma omp parallel private(i,j,jj,is,ie,w1,my_rank,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
				#endif
				{
					my_rank = omp_get_thread_num();
					LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
					for(i=is;i<ie;i++)
					{
						w1.hi = B->value[i];
						w1.lo = B->value_lo[i];
						/* t = b[i]; */
						for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
						{
							jj = A->L->index[j];
							if( jj<is || jj>=ie ) continue;
							#ifndef USE_SSE2
								LIS_QUAD_FMAD(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-A->L->value[j]);
							#else
								LIS_QUAD_FMAD_SSE2(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-A->L->value[j]);
							#endif
							/* t -= A->L->value[j] * x[A->L->index[j]]; */
						}
						#ifndef USE_SSE2
							LIS_QUAD_MULD(x[i],xl[i],w1.hi,w1.lo,A->WD->value[i]);
						#else
							LIS_QUAD_MULD_SSE2(x[i],xl[i],w1.hi,w1.lo,A->WD->value[i]);
						#endif
						/* x[i]   = t * A->WD->value[i]; */
					}
					for(i=ie-1;i>=is;i--)
					{
						w1.hi = 0.0;
						w1.lo = 0.0;
						/* t = 0.0; */
						for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
						{
							jj = A->U->index[j];
							if( jj<is || jj>=ie ) continue;
							#ifndef USE_SSE2
								LIS_QUAD_FMAD(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],A->U->value[j]);
							#else
								LIS_QUAD_FMAD_SSE2(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],A->U->value[j]);
							#endif
							/* t += A->U->value[j] * x[A->U->index[j]]; */
						}
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],w1.hi,w1.lo,-A->WD->value[i]);
						#else
							LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],w1.hi,w1.lo,-A->WD->value[i]);
						#endif
						/* x[i]  -= t * A->WD->value[i]; */
					}
				}
			#else
				for(i=0;i<n;i++)
				{
					w1.hi = B->value[i];
					w1.lo = B->value_lo[i];
					/* t = b[i]; */
					for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
					{
						jj = A->L->index[j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-A->L->value[j]);
						#else
							LIS_QUAD_FMAD_SSE2(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-A->L->value[j]);
						#endif
						/* t -= A->L->value[j] * x[A->L->index[j]]; */
					}
					#ifndef USE_SSE2
						LIS_QUAD_MULD(x[i],xl[i],w1.hi,w1.lo,A->WD->value[i]);
					#else
						LIS_QUAD_MULD_SSE2(x[i],xl[i],w1.hi,w1.lo,A->WD->value[i]);
					#endif
					/* x[i]   = t * A->WD->value[i]; */
				}
				for(i=n-1;i>=0;i--)
				{
					w1.hi = 0.0;
					w1.lo = 0.0;
					/* t = 0.0; */
					for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
					{
						jj = A->U->index[j];
						if( jj>=n ) continue; 
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],A->U->value[j]);
						#else
							LIS_QUAD_FMAD_SSE2(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],A->U->value[j]);
						#endif
						/* t += A->U->value[j] * x[A->U->index[j]]; */
					}
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],w1.hi,w1.lo,-A->WD->value[i]);
					#else
						LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],w1.hi,w1.lo,-A->WD->value[i]);
					#endif
					/* x[i]  -= t * A->WD->value[i]; */
				}
			#endif
		}
	#endif
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solveh_csr"
LIS_INT lis_matrix_solveh_csr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,jj,n;
	LIS_SCALAR t;
	LIS_SCALAR *x;
	#ifdef _OPENMP
		LIS_INT is,ie,my_rank,nprocs;
	#endif
	#ifdef USE_QUAD_PRECISION
		LIS_QUAD w1,w2;
		LIS_SCALAR *xl;
	#endif
	LIS_QUAD_DECLAR;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	x  = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,X);
		}
	#endif
	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=0;i<n;i++)
		{
			x[i]   = x[i] * conj(A->WD->value[i]);
			for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
			{
				x[A->U->index[j]] -= conj(A->U->value[j]) * x[i];
			}
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			x[i]   = x[i] * conj(A->WD->value[i]);
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				x[A->L->index[j]] -= conj(A->L->value[j]) * x[i];
			}
		}
		break;
	case LIS_MATRIX_SSOR:
	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			#ifdef _OPENMP
				nprocs = omp_get_max_threads();
				#pragma omp parallel private(i,j,jj,t,is,ie,my_rank)
				{
					my_rank = omp_get_thread_num();
					LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
					for(i=is;i<ie;i++)
					{
						t   = x[i] * conj(A->WD->value[i]);
						for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
						{
							jj = A->U->index[j];
							if( jj<is || jj>=ie ) continue;
							x[jj] -= conj(A->U->value[j]) * t;
						}
					}
					for(i=ie-1;i>=is;i--)
					{
						t    = x[i] * A->WD->value[i];
						x[i] = t;
						for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
						{
							jj = A->L->index[j];
							if( jj<is ) continue;
							x[jj] -= conj(A->L->value[j]) * t;
						}
					}
				}
			#else
				for(i=0;i<n;i++)
				{
					t   = x[i] * conj(A->WD->value[i]);
					for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
					{
						x[A->U->index[j]] -= conj(A->U->value[j]) * t;
					}
				}
				for(i=n-1;i>=0;i--)
				{
					t    = x[i] * conj(A->WD->value[i]);
					x[i] = t;
					for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
					{
						x[A->L->index[j]] -= conj(A->L->value[j]) * t;
					}
				}
			#endif
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			#ifdef _OPENMP
				nprocs = omp_get_max_threads();
				#ifndef USE_SSE2
					#pragma omp parallel private(i,j,jj,is,ie,w1,my_rank,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
				#else
					#pragma omp parallel private(i,j,jj,is,ie,w1,my_rank,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
				#endif
				{
					my_rank = omp_get_thread_num();
					LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
					for(i=is;i<ie;i++)
					{
						#ifndef USE_SSE2
							LIS_QUAD_MULD(w1.hi,w1.lo,x[i],xl[i],conj(A->WD->value[i]));
						#else
							LIS_QUAD_MULD_SSE2(w1.hi,w1.lo,x[i],xl[i],conj(A->WD->value[i]));
						#endif
						/* t   = x[i] * A->WD->value[i]; */
						for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
						{
							jj = A->U->index[j];
							if( jj<is || jj>=ie ) continue;
							#ifndef USE_SSE2
								LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],w1.hi,w1.lo,-conj(A->U->value[j]));
							#else
								LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],w1.hi,w1.lo,-conj(A->U->value[j]));
							#endif
							/* x[A->U->index[j]] -= A->U->value[j] * t; */
						}
					}
					for(i=ie-1;i>=is;i--)
					{
						#ifndef USE_SSE2
							LIS_QUAD_MULD(w1.hi,w1.lo,x[i],xl[i],conj(A->WD->value[i]));
						#else
							LIS_QUAD_MULD_SSE2(w1.hi,w1.lo,x[i],xl[i],conj(A->WD->value[i]));
						#endif
						x[i]  = w1.hi;
						xl[i] = w1.lo;
						/* t    = x[i] * conj(A->WD->value[i]); */
						/* x[i] = t; */
						for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
						{
							jj = A->L->index[j];
							if( jj<is || jj>=ie ) continue;
							#ifndef USE_SSE2
								LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],w1.hi,w1.lo,-conj(A->L->value[j]));
							#else
								LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],w1.hi,w1.lo,-conj(A->L->value[j]));
							#endif
							/* x[A->L->index[j]] -= conj(A->L->value[j]) * t; */
						}
					}
				}
			#else
				for(i=0;i<n;i++)
				{
					#ifndef USE_SSE2
						LIS_QUAD_MULD(w1.hi,w1.lo,x[i],xl[i],conj(A->WD->value[i]));
					#else
						LIS_QUAD_MULD_SSE2(w1.hi,w1.lo,x[i],xl[i],conj(A->WD->value[i]));
					#endif
					/* t   = x[i] * A->WD->value[i]; */
					for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
					{
						jj = A->U->index[j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],w1.hi,w1.lo,-conj(A->U->value[j]));
						#else
							LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],w1.hi,w1.lo,-conj(A->U->value[j]));
						#endif
						/* x[A->U->index[j]] -= conj(A->U->value[j]) * t; */
					}
				}
				for(i=n-1;i>=0;i--)
				{
					#ifndef USE_SSE2
						LIS_QUAD_MULD(w1.hi,w1.lo,x[i],xl[i],conj(A->WD->value[i]));
					#else
						LIS_QUAD_MULD_SSE2(w1.hi,w1.lo,x[i],xl[i],conj(A->WD->value[i]));
					#endif
					x[i]  = w1.hi;
					xl[i] = w1.lo;
					/* t    = x[i] * conj(A->WD->value[i]); */
					/* x[i] = t; */
					for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
					{
						jj = A->L->index[j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],w1.hi,w1.lo,-conj(A->L->value[j]));
						#else
							LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],w1.hi,w1.lo,-conj(A->L->value[j]));
						#endif
						/* x[A->L->index[j]] -= conj(A->L->value[j]) * t; */
					}
				}
			#endif
		}
	#endif
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

LIS_INT lis_matrix_ordering_mc21(LIS_MATRIX A, LIS_INT *iperm)
{
	LIS_INT n,numnz,jord,i,j,k,ii,kk,in1,in2,j1;
	LIS_INT *pr,*cv,*arp,*out;

	n = A->n;
	/*
	iperm = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_ordering_mc21:iperm");
	if( iperm==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	*/
	pr = (LIS_INT *)lis_malloc(4*n*sizeof(LIS_INT),"lis_matrix_ordering_mc21:pr");
	if( pr==NULL )
	{
		LIS_SETERR_MEM(4*n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	cv  = pr + n;
	arp = pr + 2*n;
	out = pr + 3*n;

	for(i=0;i<n;i++)
	{
		arp[i]   = A->ptr[i+1] - A->ptr[i] - 1;
		cv[i]    = -1;
		iperm[i] = -1;
	}
	numnz = 0;

	for(jord=0;jord<n;)
	{
		j = jord;
		pr[j] = -1;
		for(k=0;k<=jord;k++)
		{
			in1 = arp[j];
			if( in1>=0 )
			{
				in2 = A->ptr[j+1] - 1;
				in1 = in2 - in1;
				for(ii=in1;ii<=in2;ii++)
				{
					i = A->index[ii];
					if( iperm[i]==-1 ) goto mc21_80;
				}
				arp[j] = -1;
			}
			out[j] = A->ptr[j+1] - A->ptr[j] - 1;
			for(kk=0;kk<jord;kk++)
			{
				in1 = out[j];
				if( in1>=0 )
				{
					in2 = A->ptr[j+1] - 1;
					in1 = in2 - in1;
					for(ii=in1;ii<=in2;ii++)
					{
						i = A->index[ii];
						if( cv[i]!=jord ) break;
					}
					if( cv[i]!=jord )
					{
						j1 = j;
						j = iperm[i];
						cv[i] = jord;
						pr[j] = j1;
						out[j1] = in2 - ii - 1;
						break;
					}
				}
				j = pr[j];
				if( j==-1 ) goto mc21_100;
			}
		}
mc21_80:
		iperm[i] = j;
		arp[j] = in2 - ii - 1;
		numnz = numnz + 1;
		for(k=0;k<jord;k++)
		{
			j = pr[j];
			if( j==-1 ) break;
			ii = A->ptr[j+1] - out[j] - 2;
			i = A->index[ii];
			iperm[i] = j;
		}
mc21_100:
	jord++;
	}

	if( numnz!=n )
	{
		for(i=0;i<n;i++)
		{
			arp[i] = 0;
		}
		k = 0;
		for(i=0;i<n;i++)
		{
			if( iperm[i]!=0 )
			{
				arp[j] = 1;
			}
			else
			{
				k = k + 1;
				out[k] = 1;
			}
		}
		k = 0;
		for(i=0;i<n;i++)
		{
			if( arp[i]!=0 ) continue;
			k = k + 1;
			iperm[out[k]] = i;
		}
	}
	lis_free(pr);
	return LIS_SUCCESS;
}
