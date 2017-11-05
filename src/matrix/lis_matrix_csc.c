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
#define __FUNC__ "lis_matrix_set_csc"
LIS_INT lis_matrix_set_csc(LIS_INT nnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A)
{
	LIS_INT	err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_set_csr(nnz,ptr,index,value,A);
	if( err ) return err;

	A->status      = -LIS_MATRIX_CSC;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_setDLU_csc"
LIS_INT lis_matrix_setDLU_csc(LIS_INT nnzl, LIS_INT nnzu, LIS_SCALAR *diag, LIS_INT *lptr, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uptr, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_setDLU_csr(nnzl,nnzu,diag,lptr,lindex,lvalue,uptr,uindex,uvalue,A);
	if( err ) return err;

	A->status      = -LIS_MATRIX_CSC;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_csc"
LIS_INT lis_matrix_malloc_csc(LIS_INT n, LIS_INT nnz, LIS_INT **ptr, LIS_INT **index, LIS_SCALAR **value)
{
	LIS_INT err;
	
	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_malloc_csr(n,nnz,ptr,index,value);

	LIS_DEBUG_FUNC_OUT;
	return err;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_elements_copy_csc"
LIS_INT lis_matrix_elements_copy_csc(LIS_INT np, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_ptr, LIS_INT *o_index, LIS_SCALAR *o_value)
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
		for(i=0;i<np+1;i++)
		{
			o_ptr[i] = ptr[i];
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<np;i++)
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
#define __FUNC__ "lis_matrix_copy_csc"
LIS_INT lis_matrix_copy_csc(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT i,n,np,nnz,lnnz,unnz;
	LIS_INT *ptr,*index;
	LIS_INT *lptr,*lindex;
	LIS_INT *uptr,*uindex;
	LIS_SCALAR *value,*lvalue,*uvalue,*diag;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	np      = Ain->np;

	if( Ain->is_splited )
	{
		lnnz     = Ain->L->nnz;
		unnz     = Ain->U->nnz;
		lptr     = NULL;
		lindex   = NULL;
		uptr     = NULL;
		uindex   = NULL;
		diag     = NULL;

		err = lis_matrix_malloc_csr(np,lnnz,&lptr,&lindex,&lvalue);
		if( err )
		{
			return err;
		}
		err = lis_matrix_malloc_csr(np,unnz,&uptr,&uindex,&uvalue);
		if( err )
		{
			lis_free2(7,diag,uptr,lptr,uindex,lindex,uvalue,lvalue);
			return err;
		}
		diag = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_matrix_copy_csc::diag");
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
		lis_matrix_elements_copy_csr(np,Ain->L->ptr,Ain->L->index,Ain->L->value,lptr,lindex,lvalue);
		lis_matrix_elements_copy_csr(np,Ain->U->ptr,Ain->U->index,Ain->U->value,uptr,uindex,uvalue);

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
		err = lis_matrix_malloc_csr(np,nnz,&ptr,&index,&value);
		if( err )
		{
			return err;
		}

		lis_matrix_elements_copy_csr(np,Ain->ptr,Ain->index,Ain->value,ptr,index,value);

		err = lis_matrix_set_csr(nnz,ptr,index,value,Aout);
		if( err )
		{
			lis_free2(3,ptr,index,value);
			return err;
		}
	}
	Aout->status = -LIS_MATRIX_CSC;
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
#define __FUNC__ "lis_matrix_get_diagonal_csc"
LIS_INT lis_matrix_get_diagonal_csc(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j;
	LIS_INT n,np;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	np   = A->np;
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
		for(i=0; i<np; i++)
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
#define __FUNC__ "lis_matrix_shift_diagonal_csc"
LIS_INT lis_matrix_shift_diagonal_csc(LIS_MATRIX A, LIS_SCALAR sigma)
{
	LIS_INT i,j;
	LIS_INT n,np;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	np   = A->np;
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
		for(i=0; i<np; i++)
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
#define __FUNC__ "lis_matrix_scale_csc"
LIS_INT lis_matrix_scale_csc(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j;
	LIS_INT np;

	LIS_DEBUG_FUNC_IN;

	np   = A->np;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<np; i++)
		{
			A->D->value[i] = 1.0;
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				A->L->value[j] *= d[A->L->index[j]];
			}
			for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
			{
				A->U->value[j] *= d[A->U->index[j]];
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<np; i++)
		{
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				A->value[j] *= d[A->index[j]];
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_symm_csc"
LIS_INT lis_matrix_scale_symm_csc(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j;
	LIS_INT np;

	LIS_DEBUG_FUNC_IN;

	np   = A->np;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0; i<np; i++)
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
		for(i=0; i<np; i++)
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
#define __FUNC__ "lis_matrix_normf_csc"
LIS_INT lis_matrix_normf_csc(LIS_MATRIX A, LIS_SCALAR *nrm)
{
	LIS_INT i,j;
	LIS_INT np;
	LIS_SCALAR sum;

	LIS_DEBUG_FUNC_IN;

	np   = A->np;
	sum  = (LIS_SCALAR)0;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for reduction(+:sum) private(i,j)
		#endif
		for(i=0; i<np; i++)
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
		for(i=0; i<np; i++)
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
#define __FUNC__ "lis_matrix_transpose_csc"
LIS_INT lis_matrix_transpose_csc(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_convert_csr2csc(Ain,Aout);

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_csc"
LIS_INT lis_matrix_split_csc(LIS_MATRIX A)
{
	LIS_INT i,j,np;
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

	np       = A->np;
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
		liw = (LIS_INT *)lis_malloc((np+1)*sizeof(LIS_INT),"lis_matrix_split_csc::liw");
		if( liw==NULL )
		{
			LIS_SETERR_MEM((np+1)*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		uiw = (LIS_INT *)lis_malloc((np+1)*sizeof(LIS_INT),"lis_matrix_split_csc::uiw");
		if( uiw==NULL )
		{
			LIS_SETERR_MEM((np+1)*sizeof(LIS_INT));
			lis_free(liw);
			return LIS_OUT_OF_MEMORY;
		}
		#pragma omp parallel for private(i)
		for(i=0;i<np+1;i++)
		{
			liw[i] = 0;
			uiw[i] = 0;
		}
		#pragma omp parallel for private(i,j)
		for(i=0;i<np;i++)
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
		for(i=0;i<np;i++)
		{
			liw[i+1] += liw[i];
			uiw[i+1] += uiw[i];
		}
		nnzl = liw[np];
		nnzu = uiw[np];
	#else
		for(i=0;i<np;i++)
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
	err = lis_matrix_malloc_csr(np,nnzl,&lptr,&lindex,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_csr(np,nnzu,&uptr,&uindex,&uvalue);
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
		for(i=0;i<np+1;i++)
		{
			lptr[i] = liw[i];
			uptr[i] = uiw[i];
		}
		#pragma omp parallel for private(i,j,kl,ku)
		for(i=0;i<np;i++)
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
		for(i=0;i<np;i++)
		{
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				if( A->index[j]>i )
				{
					lindex[nnzl]   = A->index[j];
					lvalue[nnzl]   = A->value[j];
					nnzl++;
				}
				else if( A->index[j]<i )
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


#undef __FUNC__
#define __FUNC__ "lis_matrix_merge_csc"
LIS_INT lis_matrix_merge_csc(LIS_MATRIX A)
{
	LIS_INT i,j,n,np;
	LIS_INT nnz;
	LIS_INT err;
	LIS_INT *ptr,*index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;


	n       = A->n;
	np      = A->np;
	ptr     = NULL;
	index   = NULL;
	value   = NULL;
	nnz     = A->L->nnz + A->U->nnz + n;

	err = lis_matrix_malloc_csr(np,nnz,&ptr,&index,&value);
	if( err )
	{
		return err;
	}

	nnz    = 0;
	ptr[0] = 0;
	for(i=0;i<np;i++)
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
#define __FUNC__ "lis_matrix_sort_csc"
LIS_INT lis_matrix_sort_csc(LIS_MATRIX A)
{
	LIS_INT i,np;

	LIS_DEBUG_FUNC_IN;

	if( !A->is_sorted )
	{
		np = A->np;
		if( A->is_splited )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<np;i++)
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
			for(i=0;i<np;i++)
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
#define __FUNC__ "lis_matrix_solve_csc"
LIS_INT lis_matrix_solve_csc(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,n,np;
	LIS_SCALAR t;
	LIS_SCALAR *x;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	np = A->np;
	x  = X->value;

	lis_vector_copy(B,X);
	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=0;i<np;i++)
		{
			x[i]   = x[i] * A->WD->value[i];
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				x[A->L->index[j]] -= A->L->value[j] * x[i];
			}
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			x[i]   = x[i] * A->WD->value[i];
			for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
			{
				x[A->U->index[j]] -= A->U->value[j] * x[i];
			}
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<np;i++)
		{
			t   = x[i] * A->WD->value[i];
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				x[A->L->index[j]] -= A->L->value[j] * t;
			}
		}
		for(i=np-1;i>=0;i--)
		{
			t    = x[i] * A->WD->value[i];
			x[i] = t;
			for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
			{
				x[A->U->index[j]] -= A->U->value[j] * t;
			}
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solveh_csc"
LIS_INT lis_matrix_solveh_csc(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,np;
	LIS_SCALAR t;
	LIS_SCALAR *b,*x;

	LIS_DEBUG_FUNC_IN;

	np = A->np;
	b  = B->value;
	x  = X->value;

	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=0;i<np;i++)
		{
			t = b[i];
			for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
			{
				t -= conj(A->U->value[j]) * x[A->U->index[j]];
			}
			x[i]   = t * conj(A->WD->value[i]);
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=np-1;i>=0;i--)
		{
			t = b[i];
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				t -= conj(A->L->value[j]) * x[A->L->index[j]];
			}
			x[i]   = t * conj(A->WD->value[i]);
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<np;i++)
		{
			t = b[i];
			for(j=A->U->ptr[i];j<A->U->ptr[i+1];j++)
			{
				t -= conj(A->U->value[j]) * x[A->U->index[j]];
			}
			x[i]   = t * conj(A->WD->value[i]);
		}
		for(i=np-1;i>=0;i--)
		{
			t = 0.0;
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				t += conj(A->L->value[j]) * x[A->L->index[j]];
			}
			x[i]  -= t * conj(A->WD->value[i]);
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_csr2csc"
LIS_INT lis_matrix_convert_csr2csc(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,jj,k,l;
	LIS_INT err;
	LIS_INT n,nnz,np;
	LIS_INT *iw,*iw2;
	LIS_INT *ptr,*index;
	LIS_SCALAR *value;
	#ifdef _OPENMP
		LIS_INT my_rank,nprocs;
	#endif

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	np      = Ain->np;
	nnz     = Ain->nnz;
	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
	#endif

	ptr     = NULL;
	index   = NULL;
	value   = NULL;
	iw      = NULL;
	iw2     = NULL;

	ptr = (LIS_INT *)lis_malloc( (np+1)*sizeof(LIS_INT),"lis_matrix_convert_csr2csc::ptr" );
	if( ptr==NULL )
	{
		LIS_SETERR_MEM((np+1)*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	index = (LIS_INT *)lis_malloc( nnz*sizeof(LIS_INT),"lis_matrix_convert_csr2csc::index" );
	if( index==NULL )
	{
		lis_free2(5,ptr,index,value,iw,iw2);
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_convert_csr2csc::value" );
	if( value==NULL )
	{
		lis_free2(5,ptr,index,value,iw,iw2);
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	iw = (LIS_INT *)lis_malloc( (np+1)*sizeof(LIS_INT),"lis_matrix_convert_csr2csc::iw" );
	if( iw==NULL )
	{
		lis_free2(5,ptr,index,value,iw,iw2);
		LIS_SETERR_MEM((np+1)*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	#ifdef _OPENMP
		iw2 = (LIS_INT *)lis_malloc( nprocs*np*sizeof(LIS_INT),"lis_matrix_convert_csr2csc::iw2" );
		if( iw2==NULL )
		{
			lis_free2(5,ptr,index,value,iw,iw2);
			LIS_SETERR_MEM(nprocs*np*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
	#endif

	/* convert csc */
	#ifdef _OPENMP
		#pragma omp parallel private(i,j,jj,k,l,my_rank)
		{
			my_rank  = omp_get_thread_num();
			#pragma omp for
			for(i=0;i<nprocs;i++)
			{
				memset( &iw2[i*np], 0, np*sizeof(LIS_INT) );
			}
			#pragma omp for
			for(i=0;i<=np;i++)
			{
				ptr[i] = 0;
			}
			#pragma omp for
			for(i=0;i<n;i++)
			{
				for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
				{
					jj = Ain->index[j];
					iw2[my_rank*np + jj]++;
				}
			}
			#pragma omp for
			for(j=0;j<np;j++)
			{
				k = 0;
				for(i=0;i<nprocs;i++)
				{
					k += iw2[i*np+j];
				}
				iw[j] = k;
			}
			#pragma omp single
			for(j=0;j<np;j++)
			{
				ptr[j+1] = ptr[j] + iw[j];
			}
			#pragma omp for
			for(j=0;j<np;j++)
			{
				k = ptr[j];
				for(i=0;i<nprocs;i++)
				{
					l = iw2[i*np+j];
					iw2[i*np+j] = k;
					k = l + k;
				}
			}
			#pragma omp for
			for(i=0;i<np;i++)
			{
				for(j=ptr[i];j<ptr[i+1];j++)
				{
					value[j] = 0.0;
					index[j] = 0;
				}
			}
			#pragma omp for
			for(i=0;i<n;i++)
			{
				for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
				{
					k        = Ain->index[j];
					l        = iw2[my_rank*np+k];
					value[l] = Ain->value[j];
					index[l] = i;
					iw2[my_rank*np+k]++;
				}
			}
		}
	#else
		for(i=0;i<np+1;i++) iw[i] = 0;
		for(i=0;i<n;i++)
		{
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				jj = Ain->index[j];
				iw[jj]++;
			}
		}
		ptr[0] = 0;
		for(i=0;i<np;i++)
		{
			ptr[i+1] = ptr[i] + iw[i];
			iw[i]    = ptr[i];
		}

		for(i=0;i<n;i++)
		{
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				k        = Ain->index[j];
				l        = iw[k];
				value[l] = Ain->value[j];
				index[l] = i;
				iw[k]++;
			}
		}
	#endif

	err = lis_matrix_set_csc(nnz,ptr,index,value,Aout);
	if( err )
	{
		lis_free2(5,ptr,index,value,iw,iw2);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_free2(2,iw,iw2);
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	lis_free2(2,iw,iw2);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_csc2csr"
LIS_INT lis_matrix_convert_csc2csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,jj,k,l;
	LIS_INT err;
	LIS_INT n,np,nnz;
	LIS_INT *iw,*iw2;
	LIS_INT *ptr,*index;
	LIS_SCALAR *value;
	#ifdef _OPENMP
		LIS_INT nprocs,my_rank;
	#endif

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	np      = Ain->np;
	nnz     = Ain->nnz;
	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
	#endif

	ptr     = NULL;
	index   = NULL;
	value   = NULL;
	iw      = NULL;
	iw2     = NULL;

	err = lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);
	if( err )
	{
		return err;
	}
	iw = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_convert_csc2csr::iw" );
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		lis_free2(5,ptr,index,value,iw,iw2);
		return LIS_OUT_OF_MEMORY;
	}
	#ifdef _OPENMP
		iw2 = (LIS_INT *)lis_malloc( nprocs*n*sizeof(LIS_INT),"lis_matrix_convert_csc2csr::iw2" );
		if( iw2==NULL )
		{
			lis_free2(5,ptr,index,value,iw,iw2);
			LIS_SETERR_MEM(nprocs*n*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
	#endif

	/* convert csr */
	#ifdef _OPENMP
		#pragma omp parallel private(i,j,k,jj,l,my_rank)
		{
			my_rank  = omp_get_thread_num();
			#pragma omp for
			for(i=0;i<nprocs;i++)
			{
				memset( &iw2[i*n], 0, n*sizeof(LIS_INT) );
			}
			#pragma omp for
			for(i=0;i<=n;i++)
			{
				ptr[i] = 0;
			}
			#pragma omp for
			for(i=0;i<np;i++)
			{
				for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
				{
					jj = Ain->index[j];
					iw2[my_rank*n + jj]++;
				}
			}
			#pragma omp for
			for(j=0;j<n;j++)
			{
				k = 0;
				for(i=0;i<nprocs;i++)
				{
					k += iw2[i*n+j];
				}
				iw[j] = k;
			}
			#pragma omp single
			for(j=0;j<n;j++)
			{
				ptr[j+1] = ptr[j] + iw[j];
			}
			#pragma omp for
			for(j=0;j<n;j++)
			{
				k = ptr[j];
				for(i=0;i<nprocs;i++)
				{
					l = iw2[i*n+j];
					iw2[i*n+j] = k;
					k = l + k;
				}
			}
			#pragma omp for
			for(i=0;i<n;i++)
			{
				for(j=ptr[i];j<ptr[i+1];j++)
				{
					value[j] = 0.0;
					index[j] = 0;
				}
			}
			#pragma omp for
			for(i=0;i<np;i++)
			{
				for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
				{
					k        = Ain->index[j];
					l        = iw2[my_rank*n+k];
					value[l] = Ain->value[j];
					index[l] = i;
					iw2[my_rank*n+k]++;
				}
			}
		}
	#else
		for(i=0;i<n;i++) iw[i] = 0;
		for(i=0;i<np;i++)
		{
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				jj = Ain->index[j];
				iw[jj]++;
			}
		}
		ptr[0] = 0;
		for(i=0;i<n;i++)
		{
			ptr[i+1] = ptr[i] + iw[i];
			iw[i]    = ptr[i];
		}

		for(i=0;i<np;i++)
		{
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				k        = Ain->index[j];
				l        = iw[k];
				value[l] = Ain->value[j];
				index[l] = i;
				iw[k]++;
			}
		}
	#endif

	err = lis_matrix_set_csr(nnz,ptr,index,value,Aout);
	if( err )
	{
		lis_free2(5,ptr,index,value,iw,iw2);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_free2(2,iw,iw2);
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	lis_free2(2,iw,iw2);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
