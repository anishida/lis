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
#define __FUNC__ "lis_matrix_set_dia"
LIS_INT lis_matrix_set_dia(LIS_INT nnd, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A)
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
	A->status      = -LIS_MATRIX_DIA;
	A->nnd         = nnd;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_setDLU_dia"
LIS_INT lis_matrix_setDLU_dia(LIS_INT lnnd, LIS_INT unnd, LIS_SCALAR *diag, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A)
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

	A->L = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_dia::A->L");
	if( A->L==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	A->U = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_dia::A->U");
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
	A->L->nnd      = lnnd;
	A->L->index    = lindex;
	A->L->value    = lvalue;
	A->U->nnd      = unnd;
	A->U->index    = uindex;
	A->U->value    = uvalue;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_DIA;
	A->is_splited  = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_dia"
LIS_INT lis_matrix_malloc_dia(LIS_INT n, LIS_INT nnd, LIS_INT **index, LIS_SCALAR **value)
{
	LIS_DEBUG_FUNC_IN;

	*index   = NULL;
	*value   = NULL;

	*index = (LIS_INT *)lis_malloc( n*nnd*sizeof(LIS_INT),"lis_matrix_malloc_dia::index" );
	if( *index==NULL )
	{
		LIS_SETERR_MEM(n*nnd*sizeof(LIS_INT));
		lis_free2(2,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*value = (LIS_SCALAR *)lis_malloc( n*nnd*sizeof(LIS_SCALAR),"lis_matrix_malloc_dia::value" );
	if( *value==NULL )
	{
		LIS_SETERR_MEM(n*nnd*sizeof(LIS_SCALAR));
		lis_free2(2,*index,*value);
		return LIS_OUT_OF_MEMORY;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_elements_copy_dia"
LIS_INT lis_matrix_elements_copy_dia(LIS_INT n, LIS_INT nnd, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_index, LIS_SCALAR *o_value)
{
	LIS_INT is,ie;
	LIS_INT nprocs,my_rank;

	LIS_DEBUG_FUNC_IN;

	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
	#else
		nprocs = 1;
	#endif

	memcpy(o_index,index,nnd*sizeof(LIS_INT));
	#ifdef _OPENMP
	#pragma omp parallel private(is,ie,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie)

		memcpy(&o_value[is*nnd],&value[is*nnd],(ie-is)*nnd*sizeof(LIS_SCALAR));
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_dia"
LIS_INT lis_matrix_copy_dia(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT i,n,nnd,lnnd,unnd;
	LIS_INT *index;
	LIS_INT *lindex;
	LIS_INT *uindex;
	LIS_SCALAR *value,*lvalue,*uvalue,*diag;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;

	if( Ain->is_splited )
	{
		lnnd     = Ain->L->nnd;
		unnd     = Ain->U->nnd;
		lindex   = NULL;
		uindex   = NULL;
		diag     = NULL;

		err = lis_matrix_malloc_dia(n,lnnd,&lindex,&lvalue);
		if( err )
		{
			return err;
		}
		err = lis_matrix_malloc_dia(n,unnd,&uindex,&uvalue);
		if( err )
		{
			lis_free2(5,diag,uindex,lindex,uvalue,lvalue);
			return err;
		}
		diag = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_matrix_copy_dia::diag");
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
		lis_matrix_elements_copy_dia(n,lnnd,Ain->L->index,Ain->L->value,lindex,lvalue);
		lis_matrix_elements_copy_dia(n,unnd,Ain->U->index,Ain->U->value,uindex,uvalue);

		err = lis_matrix_setDLU_dia(lnnd,unnd,diag,lindex,lvalue,uindex,uvalue,Aout);
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
		nnd     = Ain->nnd;
		err = lis_matrix_malloc_dia(n,nnd,&index,&value);
		if( err )
		{
			return err;
		}

		lis_matrix_elements_copy_dia(n,nnd,Ain->index,Ain->value,index,value);

		err = lis_matrix_set_dia(nnd,index,value,Aout);
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
#define __FUNC__ "lis_matrix_get_diagonal_dia"
LIS_INT lis_matrix_get_diagonal_dia(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j;
	LIS_INT n,nnd;
	#ifdef _OPENMP
		LIS_INT is,ie,my_rank,nprocs;
	#endif

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	nnd  = A->nnd;
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
			nprocs  = omp_get_max_threads();
			for(j=0;j<nnd;j++)
			{
				if( A->index[j]==0 ) break;
			}
			#pragma omp parallel private(is,ie,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				memcpy(&d[is],&A->value[is*nnd+j*(ie-is)],(ie-is)*sizeof(LIS_SCALAR));
			}
		#else
			for(j=0;j<nnd;j++)
			{
				if( A->index[j]==0 ) break;
			}
			for(i=0;i<n;i++)
			{
				d[i] = A->value[j*n+i];
			}
		#endif
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_dia"
LIS_INT lis_matrix_shift_diagonal_dia(LIS_MATRIX A, LIS_SCALAR sigma)
{
	LIS_INT i,j,k;
	LIS_INT n,nnd;
	#ifdef _OPENMP
		LIS_INT is,ie,my_rank,nprocs;
	#endif

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	nnd  = A->nnd;
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
			nprocs  = omp_get_max_threads();
			for(j=0;j<nnd;j++)
			{
				if( A->index[j]==0 ) break;
			}
			#pragma omp parallel private(is,ie,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				/*
				memcpy(&d[is],&A->value[is*nnd+j*(ie-is)],(ie-is)*sizeof(LIS_SCALAR));
				*/
				for (k=is;k<ie;k++)
				  {
				    A->value[is*nnd+j*k] -= sigma;
				  }
			}
		#else
			for(j=0;j<nnd;j++)
			{
				if( A->index[j]==0 ) break;
			}
			for(i=0;i<n;i++)
			{
				A->value[j*n+i] -= sigma;
			}
		#endif
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_dia"
LIS_INT lis_matrix_scale_dia(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT n,np,nnd;
	#ifdef _OPENMP
		LIS_INT k,is,ie,ii;
		LIS_INT my_rank,nprocs;
	#endif


	LIS_DEBUG_FUNC_IN;

	n      = A->n;
	np     = A->np;
	if( A->is_splited )
	{
		#ifdef _OPENMP
			#pragma omp parallel private(i,j,k,is,ie,jj,js,je,ii,my_rank)
			{
				nprocs  = omp_get_max_threads();
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie)
				for(i=is;i<ie;i++)
				{
					A->D->value[i] = 1.0;
				}
				for(j=0;j<A->L->nnd;j++)
				{
					jj = A->L->index[j];
					js = _max(is,-jj);
					#ifdef USE_MPI
						je = jj<=(np-n)?ie:_min(ie,np-jj);
					#else
						je = _min(ie,n-jj);
					#endif
					k  = is*A->L->nnd + j*(ie-is);
					ii = js-is;
					for(i=js;i<je;i++)
					{
						A->L->value[k + ii++] *= d[i];
					}
				}
				for(j=0;j<A->U->nnd;j++)
				{
					jj = A->U->index[j];
					js = _max(is,-jj);
					#ifdef USE_MPI
						je = jj<=(np-n)?ie:_min(ie,np-jj);
					#else
						je = _min(ie,n-jj);
					#endif
					k  = is*A->U->nnd + j*(ie-is);
					ii = js-is;
					for(i=js;i<je;i++)
					{
						A->U->value[k + ii++] *= d[i];
					}
				}
			}
		#else
			for(i=0;i<n;i++)
			{
				A->D->value[i] = 1.0;
			}
			for(j=0;j<A->L->nnd;j++)
			{
				jj = A->L->index[j];
				js = _max(0,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?n:_min(n,np-jj);
				#else
					je = _min(n,n-jj);
				#endif
				for(i=js;i<je;i++)
				{
					A->L->value[j*n + i] *= d[i];
				}
			}
			for(j=0;j<A->U->nnd;j++)
			{
				jj = A->U->index[j];
				js = _max(0,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?n:_min(n,np-jj);
				#else
					je = _min(n,n-jj);
				#endif
				for(i=js;i<je;i++)
				{
					A->U->value[j*n + i] *= d[i];
				}
			}
		#endif
	}
	else
	{
		nnd    = A->nnd;
		#ifdef _OPENMP
			#pragma omp parallel private(i,j,k,is,ie,jj,js,je,ii,my_rank)
			{
				nprocs  = omp_get_max_threads();
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie)
				for(j=0;j<nnd;j++)
				{
					jj = A->index[j];
					js = _max(is,-jj);
					#ifdef USE_MPI
						je = jj<=(np-n)?ie:_min(ie,np-jj);
					#else
						je = _min(ie,n-jj);
					#endif
					k  = is*nnd + j*(ie-is);
					ii = js-is;
					for(i=js;i<je;i++)
					{
						A->value[k + ii] *= d[i];
						ii++;
					}
				}
			}
		#else
			for(j=0;j<nnd;j++)
			{
				jj = A->index[j];
				js = _max(0,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?n:_min(n,np-jj);
				#else
					je = _min(n,n-jj);
				#endif
				for(i=js;i<je;i++)
				{
					A->value[j*n + i] *= d[i];
				}
			}
		#endif
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_symm_dia"
LIS_INT lis_matrix_scale_symm_dia(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT n,np,nnd;
	#ifdef _OPENMP
		LIS_INT k,is,ie,ii;
		LIS_INT my_rank,nprocs;
	#endif

	LIS_DEBUG_FUNC_IN;

	n      = A->n;
	np     = A->np;
	nnd  = A->nnd;
	if( A->is_splited )
	{
		#ifdef _OPENMP
			#pragma omp parallel private(i,j,k,is,ie,jj,js,je,ii,my_rank)
			{
				nprocs  = omp_get_max_threads();
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie)
				for(i=is;i<ie;i++)
				{
					A->D->value[i] = 1.0;
				}
				for(j=0;j<A->L->nnd;j++)
				{
					jj = A->L->index[j];
					js = _max(is,-jj);
					#ifdef USE_MPI
						je = jj<=(np-n)?ie:_min(ie,np-jj);
					#else
						je = _min(ie,n-jj);
					#endif
					k  = is*A->L->nnd + j*(ie-is);
					ii = js-is;
					for(i=js;i<je;i++)
					{
						A->L->value[k + ii++] *= d[i]*d[i+A->L->index[j]];;
					}
				}
				for(j=0;j<A->U->nnd;j++)
				{
					jj = A->U->index[j];
					js = _max(is,-jj);
					#ifdef USE_MPI
						je = jj<=(np-n)?ie:_min(ie,np-jj);
					#else
						je = _min(ie,n-jj);
					#endif
					k  = is*A->U->nnd + j*(ie-is);
					ii = js-is;
					for(i=js;i<je;i++)
					{
						A->U->value[k + ii++] *= d[i]*d[i+A->U->index[j]];;
					}
				}
			}
		#else
			for(i=0;i<n;i++)
			{
				A->D->value[i] = 1.0;
			}
			for(j=0;j<A->L->nnd;j++)
			{
				jj = A->L->index[j];
				js = _max(0,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?n:_min(n,np-jj);
				#else
					je = _min(n,n-jj);
				#endif
				for(i=js;i<je;i++)
				{
					A->L->value[j*n + i] *= d[i]*d[i+A->L->index[j]];;
				}
			}
			for(j=0;j<A->U->nnd;j++)
			{
				jj = A->U->index[j];
				js = _max(0,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?n:_min(n,np-jj);
				#else
					je = _min(n,n-jj);
				#endif
				for(i=js;i<je;i++)
				{
					A->U->value[j*n + i] *= d[i]*d[i+A->U->index[j]];;
				}
			}
		#endif
	}
	else
	{
		nnd    = A->nnd;
		#ifdef _OPENMP
			#pragma omp parallel private(i,j,k,is,ie,jj,js,je,ii,my_rank)
			{
				nprocs  = omp_get_max_threads();
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie)
				for(j=0;j<nnd;j++)
				{
					jj = A->index[j];
					js = _max(is,-jj);
					#ifdef USE_MPI
						je = jj<=(np-n)?ie:_min(ie,np-jj);
					#else
						je = _min(ie,n-jj);
					#endif
					k  = is*nnd + j*(ie-is);
					ii = js-is;
					for(i=js;i<je;i++)
					{
						A->value[k + ii] *= d[i]*d[i+A->index[j]];
						ii++;
					}
				}
			}
		#else
			for(j=0;j<nnd;j++)
			{
				jj = A->index[j];
				js = _max(0,-jj);
				#ifdef USE_MPI
					je = jj<=(np-n)?n:_min(n,np-jj);
				#else
					je = _min(n,n-jj);
				#endif
				for(i=js;i<je;i++)
				{
					A->value[j*n + i] *= d[i]*d[i+A->index[j]];
				}
			}
		#endif
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_normf_dia"
LIS_INT lis_matrix_normf_dia(LIS_MATRIX A, LIS_SCALAR *nrm)
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
#define __FUNC__ "lis_matrix_transpose_dia"
LIS_INT lis_matrix_transpose_dia(LIS_MATRIX Ain, LIS_MATRIX *Aout)
{

	LIS_DEBUG_FUNC_IN;

/*	err = lis_matrix_convert_dia2csc(Ain,Aout);*/
	(*Aout)->matrix_type = LIS_MATRIX_DIA;
	(*Aout)->status      = LIS_MATRIX_DIA;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_dia"
LIS_INT lis_matrix_split_dia(LIS_MATRIX A)
{
	LIS_INT i,j,n,nnd;
	LIS_INT lnnd,unnd;
	#ifdef _OPENMP
		LIS_INT kl,ku,nprocs,my_rank,is,ie;
	#endif
	LIS_INT err;
	LIS_INT *lindex,*uindex;
	LIS_SCALAR *lvalue,*uvalue;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;

	n        = A->n;
	nnd      = A->nnd;
	lnnd     = 0;
	unnd     = 0;
	D        = NULL;
	lindex   = NULL;
	lvalue   = NULL;
	uindex   = NULL;
	uvalue   = NULL;

	for(j=0;j<nnd;j++)
	{
		if( A->index[j]<0 )
		{
			lnnd++;
		}
		else if( A->index[j]>0 )
		{
			unnd++;
		}
	}

	err = lis_matrix_LU_create(A);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_dia(n,lnnd,&lindex,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_dia(n,unnd,&uindex,&uvalue);
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
		kl = 0;
		ku = 0;
		nprocs  = omp_get_max_threads();
		for(j=0;j<nnd;j++)
		{
			if( A->index[j]<0 )
			{
				lindex[kl]   = A->index[j];
				#pragma omp parallel private(i,is,ie,my_rank)
				{
					my_rank = omp_get_thread_num();
					LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
					memcpy(&lvalue[is*lnnd+kl*(ie-is)],&A->value[is*nnd+j*(ie-is)],(ie-is)*sizeof(LIS_SCALAR));
				}
				kl++;
			}
			else if( A->index[j]>0 )
			{
				uindex[ku]   = A->index[j];
				#pragma omp parallel private(i,is,ie,my_rank)
				{
					my_rank = omp_get_thread_num();
					LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
					memcpy(&uvalue[is*unnd+ku*(ie-is)],&A->value[is*nnd+j*(ie-is)],(ie-is)*sizeof(LIS_SCALAR));
				}
				ku++;
			}
			else
			{
				#pragma omp parallel private(i,is,ie,my_rank)
				{
					my_rank = omp_get_thread_num();
					LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
					for(i=is;i<ie;i++)
					{
						D->value[i] = A->value[is*nnd+j*(ie-is)+i-is];
					}
				}
			}
		}
	#else
		lnnd = 0;
		unnd = 0;
		for(j=0;j<nnd;j++)
		{
			if( A->index[j]<0 )
			{
				lindex[lnnd]   = A->index[j];
				for(i=0;i<n;i++)
				{
					lvalue[lnnd*n+i] = A->value[j*n+i];
				}
				lnnd++;
			}
			else if( A->index[j]>0 )
			{
				uindex[unnd]   = A->index[j];
				for(i=0;i<n;i++)
				{
					uvalue[unnd*n+i] = A->value[j*n+i];
				}
				unnd++;
			}
			else
			{
				for(i=0;i<n;i++)
				{
					D->value[i] = A->value[j*n+i];
				}
			}
		}
	#endif
	A->L->nnd     = lnnd;
	A->L->index   = lindex;
	A->L->value   = lvalue;
	A->U->nnd     = unnd;
	A->U->index   = uindex;
	A->U->value   = uvalue;
	A->D          = D;
	A->is_splited = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_merge_dia"
LIS_INT lis_matrix_merge_dia(LIS_MATRIX A)
{
        LIS_INT i,j,k,n,is;
	#ifdef _OPENMP
		LIS_INT nprocs,my_rank,ie;
	#endif
	LIS_INT nnd;
	LIS_INT err;
	LIS_INT *index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;


	n       = A->n;
	is      = A->is;
	index   = NULL;
	value   = NULL;
	nnd     = A->L->nnd + A->U->nnd + 1;

	err = lis_matrix_malloc_dia(n,nnd,&index,&value);
	if( err )
	{
		return err;
	}

	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
		k = 0;
		for(j=0;j<A->L->nnd;j++)
		{
			index[k]   = A->L->index[j];
			#pragma omp parallel private(is,ie,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				memcpy(&value[is*nnd+k*(ie-is)],&A->L->value[is*A->L->nnd+j*(ie-is)],(ie-is)*sizeof(LIS_SCALAR));
			}
			k++;
		}
		index[k]   = 0;
		#pragma omp parallel private(is,ie,my_rank)
		{
			my_rank = omp_get_thread_num();
			LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
			memcpy(&value[is*nnd+k*(ie-is)],&A->D->value[is],(ie-is)*sizeof(LIS_SCALAR));
		}
		k++;
		for(j=0;j<A->U->nnd;j++)
		{
			index[k]   = A->U->index[j];
			#pragma omp parallel private(is,ie,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				memcpy(&value[is*nnd+k*(ie-is)],&A->U->value[is*A->U->nnd+j*(ie-is)],(ie-is)*sizeof(LIS_SCALAR));
			}
			k++;
		}
	#else
		k = 0;
		for(j=0;j<A->L->nnd;j++)
		{
			index[k]   = A->L->index[j];
			for(i=0;i<n;i++)
			{
				value[k*n+i] = A->L->value[j*n+i];
			}
			k++;
		}
		index[k]   = 0;
		for(i=0;i<n;i++)
		{
			value[k*n+i] = A->D->value[i];
		}
		k++;
		for(j=0;j<A->U->nnd;j++)
		{
			index[k]   = A->U->index[j];
			for(i=0;i<n;i++)
			{
				value[k*n+i] = A->U->value[j*n+i];
			}
			k++;
		}
	#endif

	A->nnd        = nnd;
	A->value      = value;
	A->index      = index;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_sort_dia"
LIS_INT lis_matrix_sort_dia(LIS_MATRIX A)
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
#define __FUNC__ "lis_matrix_solve_dia"
LIS_INT lis_matrix_solve_dia(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
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
			for(j=0;j<A->L->nnd;j++)
			{
				if( i+A->L->index[j] >= 0 ) t -= A->L->value[j*n + i] * x[i + A->L->index[j]];
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			t = b[i];
			for(j=0;j<A->U->nnd;j++)
			{
				if( i+A->U->index[j] < n ) t -= A->U->value[j*n + i] * x[i + A->U->index[j]];
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<n;i++)
		{
			t = b[i];
			for(j=0;j<A->L->nnd;j++)
			{
				if( i+A->L->index[j] >= 0 ) t -= A->L->value[j*n + i] * x[i + A->L->index[j]];
			}
			x[i]   = t * A->WD->value[i];
		}
		for(i=n-1;i>=0;i--)
		{
			t = 0.0;
			for(j=0;j<A->U->nnd;j++)
			{
				if( i+A->U->index[j] < n ) t += A->U->value[j*n + i] * x[i + A->U->index[j]];
			}
			x[i]  -= t * A->WD->value[i];
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solveh_dia"
LIS_INT lis_matrix_solveh_dia(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
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
		for(i=0;i<n;i++)
		{
			x[i]  = x[i] * conj(A->WD->value[i]);
			for(j=0;j<A->U->nnd;j++)
			{
				if( i+A->U->index[j] < n ) x[i + A->U->index[j]] -= conj(A->U->value[j*n + i]) * x[i];
			}
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			x[i]  = x[i] * conj(A->WD->value[i]);
			for(j=0;j<A->L->nnd;j++)
			{
				if( i+A->L->index[j] >= 0 ) x[i + A->L->index[j]] -= conj(A->L->value[j*n + i]) * x[i];
			}
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<n;i++)
		{
			t  = x[i] * conj(A->WD->value[i]);
			for(j=0;j<A->U->nnd;j++)
			{
				if( i+A->U->index[j] < n ) x[i + A->U->index[j]] -= conj(A->U->value[j*n + i]) * t;
			}
		}
		for(i=n-1;i>=0;i--)
		{
			x[i]  = x[i] * conj(A->WD->value[i]);
			t     = x[i];
			for(j=0;j<A->L->nnd;j++)
			{
				if( i+A->L->index[j] >= 0 ) x[i + A->L->index[j]] -= conj(A->L->value[j*n + i]) * t;
			}
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_csr2dia"
LIS_INT lis_matrix_convert_csr2dia(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,jj,k;
	LIS_INT err;
	LIS_INT n,nnz,nnd,nprocs,my_rank;
	LIS_INT is,ie;
	LIS_INT *iw;
	LIS_INT *index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	nnz		= Ain->nnz;

	index   = NULL;
	value   = NULL;
	iw      = NULL;

	iw = (LIS_INT *)lis_malloc( nnz*sizeof(LIS_INT),"lis_matrix_convert_csr2dia::iw" );
	if( iw==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		return LIS_ERR_OUT_OF_MEMORY;
	}

	lis_matrix_sort_csr(Ain);

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j)
	#endif
	for(i=0;i<n;i++)
	{
		for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
		{
			iw[j] = Ain->index[j] - i;
		}
	}
	lis_sort_i(0,nnz-1,iw);
	nnd = 1;
	k   = iw[0];
	for(i=1;i<nnz;i++)
	{
		if( k!=iw[i] )
		{
			k = iw[i];
			nnd++;
		}
	}

	err = lis_matrix_malloc_dia(n,nnd,&index,&value);
	if( err )
	{
		lis_free(iw);
		return err;
	}

	/* convert dia */
	k        = iw[0];
	index[0] = k;
	j        = 1;
	for(i=1;i<nnz;i++)
	{
		if( k!=iw[i] )
		{
			k = iw[i];
			index[j] = k;
			j++;
		}
	}
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k,is,ie,jj,my_rank)
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
		memset(&value[is*nnd],0,(ie-is)*nnd*sizeof(LIS_SCALAR));

		for(i=is;i<ie;i++)
		{
			k = 0;
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				jj=Ain->index[j] - i;
				while( jj!=index[k] ) k++;
				value[is*nnd + k*(ie-is) + i-is] = Ain->value[j];
			}
		}
	}

	err = lis_matrix_set_dia(nnd,index,value,Aout);
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
#define __FUNC__ "lis_matrix_convert_dia2csr"
LIS_INT lis_matrix_convert_dia2csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,jj,k,l;
	LIS_INT err,js,je;
	LIS_INT n,np,nnz,nnd,is,ie,nprocs,my_rank;
	LIS_INT *iw;
	LIS_INT *ptr,*index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	np      = Ain->np;
	nnd     = Ain->nnd;
	is      = Ain->is;
	ie      = Ain->ie;

	ptr     = NULL;
	index   = NULL;
	value   = NULL;
	iw      = NULL;

	iw = (LIS_INT *)lis_malloc( (n+1)*sizeof(LIS_INT),"lis_matrix_convert_dia2csr::iw" );
	if( iw==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
		return LIS_ERR_OUT_OF_MEMORY;
	}

	iw[0]  = 0;
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k,is,ie,js,je,jj,my_rank)
	#endif
	{
		#ifdef _OPENMP
			nprocs  = omp_get_max_threads();
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
			nprocs  = 1;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
		memset(&iw[is+1],0,(ie-is)*sizeof(LIS_INT));
		k = ie-is;

		for(j=0;j<nnd;j++)
		{
			jj = Ain->index[j];
			js = _max(is,-jj) -is;
			#ifdef USE_MPI
				je = jj<=(np-n)?ie:_min(ie,np-jj);
			#else
				je = _min(ie,n-jj);
			#endif
			je -= is;
			for(i=js;i<je;i++)
			{ 
				if( Ain->value[is*nnd + j*k + i]!=(LIS_SCALAR)0.0 )
				{
					iw[i+is+1]++;
				}
			}
		}
	}
	for(i=0;i<n;i++)
	{
		iw[i+1] += iw[i];
	}
	nnz = iw[n];

	err = lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);
	if( err )
	{
		lis_free2(4,ptr,index,value,iw);
		return err;
	}

	/* convert csr */
	ptr[0] = 0;
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k,l,is,ie,js,je,jj,my_rank)
	#endif
	{
		#ifdef _OPENMP
			nprocs  = omp_get_max_threads();
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
			nprocs  = 1;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
		l = ie-is;

		for(i=is;i<ie;i++)
		{
			ptr[i+1] = iw[i+1];
		}
		for(j=0;j<nnd;j++)
		{
			jj = Ain->index[j];
			js = _max(is,-jj) - is;
			#ifdef USE_MPI
				je = jj<=(np-n)?ie:_min(ie,np-jj);
			#else
				je = _min(ie,n-jj);
			#endif
			je -= is;
			for(i=js;i<je;i++)
			{ 
				if( Ain->value[is*nnd + j*l + i]!=(LIS_SCALAR)0.0 )
				{
					k        = iw[i+is]++;
					value[k] = Ain->value[is*nnd + j*l + i];
					index[k] = i+is + jj;
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

