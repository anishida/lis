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
 * lis_matrix_set
 * lis_matrix_malloc
 * lis_matrix_copy
 * lis_matrix_convert
 * lis_matrix_get_diagonal
 * lis_matrix_shift_diagonal
 * lis_matrix_scale
 * lis_matrix_scale_symm
 * lis_matrix_normf
 * lis_matrix_transpose
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_bsr"
LIS_INT lis_matrix_set_bsr(LIS_INT bnr, LIS_INT bnc, LIS_INT bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_MATRIX A)
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
	A->bptr        = bptr;
	A->bindex      = bindex;
	A->value       = value;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_BSR;
	A->is_block    = LIS_TRUE;
	A->bnnz        = bnnz;
	A->nr          = (A->n-1)/bnr+1;
	A->nc          = (A->gn-1)/bnc+1;
	if( A->n==A->np )
	{
		A->nc      = 1 + (A->n - 1)/bnc;
		A->pad     = (bnc - A->n%bnc)%bnc;
	}
	else
	{
		A->nc      = 2 + (A->n - 1)/bnc + (A->np - A->n - 1)/bnc;
		A->pad     = (bnc - A->n%bnc)%bnc + (bnc - (A->np-A->n)%bnc)%bnc;
	}
	A->bnr         = bnr;
	A->bnc         = bnc;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_setDLU_bsr"
LIS_INT lis_matrix_setDLU_bsr(LIS_INT bnr, LIS_INT bnc, LIS_INT lbnnz, LIS_INT ubnnz, LIS_MATRIX_DIAG D, LIS_INT *lbptr, LIS_INT *lbindex, LIS_SCALAR *lvalue, LIS_INT *ubptr, LIS_INT *ubindex, LIS_SCALAR *uvalue, LIS_MATRIX A)
{
	LIS_INT	err;

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

	A->L = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_bsr::A->L");
	if( A->L==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	A->U = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_setDLU_bsr::A->U");
	if( A->U==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		lis_matrix_DLU_destroy(A);
		return LIS_OUT_OF_MEMORY;
	}

	A->D           = D;
	A->L->bnnz     = lbnnz;
	A->L->bptr     = lbptr;
	A->L->bindex   = lbindex;
	A->L->value    = lvalue;
	A->U->bnnz     = ubnnz;
	A->U->bptr     = ubptr;
	A->U->bindex   = ubindex;
	A->U->value    = uvalue;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_BSR;
	A->is_splited  = LIS_TRUE;
	A->is_block    = LIS_TRUE;
	A->nr          = (A->n-1)/bnr+1;
	A->nc          = (A->gn-1)/bnc+1;
	if( A->n==A->np )
	{
		A->nc      = 1 + (A->n - 1)/bnc;
		A->pad     = (bnc - A->n%bnc)%bnc;
	}
	else
	{
		A->nc      = 2 + (A->n - 1)/bnc + (A->np - A->n - 1)/bnc;
		A->pad     = (bnc - A->n%bnc)%bnc + (bnc - (A->np-A->n)%bnc)%bnc;
	}
	A->bnr         = bnr;
	A->bnc         = bnc;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_bsr"
LIS_INT lis_matrix_malloc_bsr(LIS_INT n, LIS_INT bnr, LIS_INT bnc, LIS_INT bnnz, LIS_INT **bptr, LIS_INT **bindex, LIS_SCALAR **value)
{
	LIS_INT	nr;

	LIS_DEBUG_FUNC_IN;

	nr        = 1 + (n -1)/bnr;
	*bptr     = NULL;
	*bindex   = NULL;
	*value    = NULL;

	*bptr = (LIS_INT *)lis_malloc( (nr+1)*sizeof(LIS_INT),"lis_matrix_malloc_bsr::bptr" );
	if( *bptr==NULL )
	{
		LIS_SETERR_MEM((nr+1)*sizeof(LIS_INT));
		lis_free2(3,*bptr,*bindex,*value);
		return LIS_FAILS;
	}
	*bindex = (LIS_INT *)lis_malloc( bnnz*sizeof(LIS_INT),"lis_matrix_malloc_bsr::bindex" );
	if( *bindex==NULL )
	{
		LIS_SETERR_MEM(bnnz*sizeof(LIS_INT));
		lis_free2(3,*bptr,*bindex,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*value = (LIS_SCALAR *)lis_malloc( bnnz*bnr*bnc*sizeof(LIS_SCALAR),"lis_matrix_malloc_bsr::value" );
	if( *value==NULL )
	{
		LIS_SETERR_MEM(bnnz*bnr*bnc*sizeof(LIS_SCALAR));
		lis_free2(3,*bptr,*bindex,*value);
		return LIS_OUT_OF_MEMORY;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_elements_copy_bsr"
LIS_INT lis_matrix_elements_copy_bsr(LIS_INT n, LIS_INT bnr, LIS_INT bnc, LIS_INT bnnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_ptr, LIS_INT *o_index, LIS_SCALAR *o_value)
{
	LIS_INT	i,j,k;
	LIS_INT	nr,bs;

	LIS_DEBUG_FUNC_IN;

	nr  = 1 + (n - 1)/bnr;
	bs  = bnr*bnc;
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k)
	#endif
	{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<nr+1;i++)
		{
			o_ptr[i] = ptr[i];
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<nr;i++)
		{
			for(j=ptr[i];j<ptr[i+1];j++)
			{
				for(k=0;k<bs;k++)
				{
					o_value[j*bs+k]   = value[j*bs+k];
				}
				o_index[j]   = index[j];
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_bsr"
LIS_INT lis_matrix_copy_bsr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT n,bnnz,bnr,bnc;
	LIS_INT lbnnz,ubnnz;
	LIS_INT *bptr,*bindex;
	LIS_INT *lbptr,*lbindex;
	LIS_INT *ubptr,*ubindex;
	LIS_SCALAR	*value,*lvalue,*uvalue;
	LIS_MATRIX_DIAG D;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	bnnz    = Ain->bnnz;
	bnr     = Ain->bnr;
	bnc     = Ain->bnc;
    
	if( Ain->is_splited )
	{
		lbnnz    = Ain->L->bnnz;
		ubnnz    = Ain->U->bnnz;
		lbptr    = NULL;
		lbindex  = NULL;
		lvalue   = NULL;
		ubptr    = NULL;
		ubindex  = NULL;
		uvalue   = NULL;
		D        = NULL;

		err = lis_matrix_malloc_bsr(n,bnr,bnc,lbnnz,&lbptr,&lbindex,&lvalue);
		if( err )
		{
			return err;
		}
		err = lis_matrix_malloc_bsr(n,bnr,bnc,ubnnz,&ubptr,&ubindex,&uvalue);
		if( err )
		{
			lis_free2(6,ubptr,lbptr,ubindex,lbindex,uvalue,lvalue);
			return err;
		}
		err = lis_matrix_diag_duplicateM(Ain,&D);
		if( err )
		{
			lis_free2(6,ubptr,lbptr,ubindex,lbindex,uvalue,lvalue);
			return err;
		}

		lis_matrix_diag_copy(Ain->D,D);
		lis_matrix_elements_copy_bsr(n,bnr,bnc,lbnnz,Ain->L->bptr,Ain->L->bindex,Ain->L->value,lbptr,lbindex,lvalue);
		lis_matrix_elements_copy_bsr(n,bnr,bnc,ubnnz,Ain->U->bptr,Ain->U->bindex,Ain->U->value,ubptr,ubindex,uvalue);

		err = lis_matrix_setDLU_bsr(bnr,bnc,lbnnz,ubnnz,D,lbptr,lbindex,lvalue,ubptr,ubindex,uvalue,Aout);
		if( err )
		{
			lis_free2(6,ubptr,lbptr,ubindex,lbindex,uvalue,lvalue);
			return err;
		}
	}
	if( !Ain->is_splited || (Ain->is_splited && Ain->is_save) )
	{
		bptr    = NULL;
		bindex  = NULL;
		value   = NULL;

		err = lis_matrix_malloc_bsr(n,bnr,bnc,bnnz,&bptr,&bindex,&value);
		if( err )
		{
			return err;
		}

		lis_matrix_elements_copy_bsr(n,bnr,bnc,bnnz,Ain->bptr,Ain->bindex,Ain->value,bptr,bindex,value);

		err = lis_matrix_set_bsr(bnr,bnc,bnnz,bptr,bindex,value,Aout);
		if( err )
		{
			lis_free2(3,bptr,bindex,value);
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
#define __FUNC__ "lis_matrix_convert_csr2bsr"
LIS_INT lis_matrix_convert_csr2bsr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,k,n,bnr,bnc;
	LIS_INT ii,jj,kk,pad;
	LIS_INT bnnz,bj,nr,nc,jpos,nnz,ij,kv,bi;
	LIS_INT err;
	LIS_INT np,nprocs,my_rank;
	LIS_INT *iw,*iw2;
	LIS_INT *bptr,*bindex;
	LIS_SCALAR	*value;

	LIS_DEBUG_FUNC_IN;

	bnr   = Aout->conv_bnr;
	bnc   = Aout->conv_bnc;

	n       = Ain->n;
	np      = Ain->np;
	nr      = 1 + (n - 1)/bnr;
	pad     = (bnc - n%bnc)%bnc;
	if( n==np )
	{
		nc      = 1 + (n - 1)/bnc;
	}
	else
	{
		nc      = 2 + (n - 1)/bnc + (pad + np - n - 1)/bnc;
	}

	bptr    = NULL;
	bindex  = NULL;
	value   = NULL;
	iw      = NULL;
	iw2     = NULL;

	bptr = (LIS_INT *)lis_malloc( (nr+1)*sizeof(LIS_INT),"lis_matrix_convert_csr2bsr::bptr" );
	if( bptr==NULL )
	{
		LIS_SETERR_MEM((nr+1)*sizeof(LIS_INT));
		lis_free2(5,bptr,bindex,value,iw,iw2);
		return LIS_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
	#else
		nprocs = 1;
	#endif
	iw    = (LIS_INT *)lis_malloc( nprocs*nc*sizeof(LIS_INT),"lis_matrix_convert_csr2bsr::iw" );
	iw2   = (LIS_INT *)lis_malloc( nprocs*nc*sizeof(LIS_INT),"lis_matrix_convert_csr2bsr::iw2" );

	#ifdef _OPENMP
	#pragma omp parallel private(i,k,ii,j,bj,kk,ij,jj,my_rank,kv,jpos)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		memset(&iw[my_rank*nc],0,nc*sizeof(LIS_INT));

		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<nr;i++)
		{
			k = 0;
			kk   = bnr*i;
			jj   = 0;
			for(ii=0;ii+kk<n&&ii<bnr;ii++)
			{
				for(j=Ain->ptr[kk+ii];j<Ain->ptr[kk+ii+1];j++)
				{
					#ifdef USE_MPI
						bj = Ain->index[j];
						bj = (bj<n)?bj/bnc:(bj+pad)/bnc;
					#else
						bj   = Ain->index[j]/bnc;
					#endif
					jpos = iw[my_rank*nc + bj];
					if( jpos==0 )
					{
						iw[my_rank*nc + bj] = 1;
						iw2[my_rank*nc + jj] = bj;
						jj++;
					}
				}
			}
			for(bj=0;bj<jj;bj++)
			{
				k++;
				ii = iw2[my_rank*nc + bj];
				iw[my_rank*nc + ii]=0;
			}
			bptr[i+1] = k;
		}
	}

	bptr[0] = 0;
	for(i=0;i<nr;i++)
	{
		bptr[i+1] += bptr[i];
	}
	bnnz = bptr[nr];
	nnz  = bnnz*bnr*bnc;


	bindex = (LIS_INT *)lis_malloc( bnnz*sizeof(LIS_INT),"lis_matrix_convert_csr2bsr::bindex" );
	if( bindex==NULL )
	{
		LIS_SETERR_MEM((nr+1)*sizeof(LIS_INT));
		lis_free2(5,bptr,bindex,value,iw,iw2);
		return LIS_OUT_OF_MEMORY;
	}
	value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_convert_csr2bsr::value" );
	if( value==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		lis_free2(5,bptr,bindex,value,iw,iw2);
		return LIS_OUT_OF_MEMORY;
	}


	/* convert bsr */
	#ifdef _OPENMP
	#pragma omp parallel private(bi,i,ii,k,j,bj,jpos,kv,kk,ij,jj,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		memset(&iw[my_rank*nc],0,nc*sizeof(LIS_INT));

		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(bi=0;bi<nr;bi++)
		{
			i  = bi*bnr;
			ii = 0;
			kk = bptr[bi];
			while( i+ii<n && ii<=bnr-1 )
			{
				for( k=Ain->ptr[i+ii];k<Ain->ptr[i+ii+1];k++)
				{
					#ifdef USE_MPI
						j    = (Ain->index[k]<n)?Ain->index[k]:Ain->index[k]+pad;
					#else
						j    = Ain->index[k];
					#endif
					bj   = j/bnc;
					j    = j%bnc;
					jpos = iw[my_rank*nc + bj];
					if( jpos==0 )
					{
						kv                  = kk * bnr * bnc;
						iw[my_rank*nc + bj] = kv+1;
						bindex[kk]          = bj;
						for(jj=0;jj<bnr*bnc;jj++) value[kv+jj] = 0.0;
						ij = j*bnr + ii;
						value[kv+ij]   = Ain->value[k];
						kk = kk+1;
					}
					else
					{
						ij = j*bnr + ii;
						value[jpos+ij-1]   = Ain->value[k];
					}
				}
				ii = ii+1;
			}
			for(j=bptr[bi];j<bptr[bi+1];j++)
			{
				iw[my_rank*nc + bindex[j]] = 0;
			}
		}
	}
	lis_free2(2,iw,iw2);

	err = lis_matrix_set_bsr(bnr,bnc,bnnz,bptr,bindex,value,Aout);
	if( err )
	{
		lis_free2(3,bptr,bindex,value);
		return err;
	}
	Aout->pad_comm = pad;
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	#ifdef USE_MPI
		Aout->commtable->pad = pad;
		MPI_Barrier(Ain->comm);
	#endif
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_bsr2csr"
LIS_INT lis_matrix_convert_bsr2csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,k;
	LIS_INT nr,bnr,bnc,bs,bi,bj;
	LIS_INT err;
	LIS_INT n,nnz;
	LIS_INT *ptr,*index;
	LIS_SCALAR	*value;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	nr      = Ain->nr;
	bnr     = Ain->bnr;
	bnc     = Ain->bnc;
	bs      = bnr*bnc;

	ptr     = NULL;
	index   = NULL;
	value   = NULL;

	ptr = (LIS_INT *)lis_malloc( (n+1)*sizeof(LIS_INT),"lis_matrix_convert_bsr2csr::ptr" );
	if( ptr==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	/* check nnz */
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,bi,bj)
	#endif
	{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n+1;i++)
		{
			ptr[i] = 0;
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(bi=0;bi<nr;bi++)
		{
			for(bj=Ain->bptr[bi];bj<Ain->bptr[bi+1];bj++)
			{
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						if( Ain->value[bj*bs + j*bnr + i] != (LIS_SCALAR)0.0 )
						{
							ptr[bi*bnr+i+1]++;
						}
					}
				}
			}
		}
	}
	for(i=0;i<n;i++)
	{
		ptr[i+1] += ptr[i];
	}
	nnz = ptr[n];

	index = (LIS_INT *)lis_malloc( nnz*sizeof(LIS_INT),"lis_matrix_convert_bsr2csr::index" );
	if( index==NULL )
	{
		lis_free2(3,ptr,index,value);
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_convert_bsr2csr::value" );
	if( value==NULL )
	{
		lis_free2(3,ptr,index,value);
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	/* convert csr */
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,bi,bj,k)
	#endif
	for(bi=0;bi<nr;bi++)
	{
		for(i=0;i<bnr;i++)
		{
			if( bi*bnr+i==n ) break;
			k = ptr[bi*bnr+i];
			for(bj=Ain->bptr[bi];bj<Ain->bptr[bi+1];bj++)
			{
				for(j=0;j<bnc;j++)
				{
					if( Ain->value[bj*bs + j*bnr + i] != (LIS_SCALAR)0.0 )
					{
						value[k]   = Ain->value[bj*bs + j*bnr + i];
						#ifdef USE_MPI
							index[k]   = Ain->bindex[bj]*bnc + j;
							if( index[k]>=n ) index[k] -= Ain->pad_comm;
						#else
							index[k]   = Ain->bindex[bj]*bnc + j;
						#endif
						k++;
					}
				}
			}
		}
	}

	err = lis_matrix_set_csr(nnz,ptr,index,value,Aout);
	if( err )
	{
		lis_free2(3,ptr,index,value);
		return err;
	}
	Aout->pad      = 0;
	Aout->pad_comm = 0;
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	#ifdef USE_MPI
		Aout->commtable->pad = 0;
	#endif
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal_bsr"
LIS_INT lis_matrix_get_diagonal_bsr(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j,k,bi,bj,bjj,nr;
	LIS_INT bnr,bnc,bs;
	LIS_INT n;

	LIS_DEBUG_FUNC_IN;

	n   = A->n;
	nr  = A->nr;
	bnr = A->bnr;
	bnc = A->bnc;
	bs  = bnr*bnc;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0;i<nr;i++)
		{
			for(j=0;j<bnr;j++)
			{
				d[i*bnr+j] = A->D->value[i*bs+j*bnr+j];
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj,bjj,i,j,k)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			k = 0;
			i = bi*bnr;
			for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
			{
				bjj = A->bindex[bj];
				if( i>=bjj*bnc && i<(bjj+1)*bnc )
				{
					for(j=i%bnc;j<bnc&&k<bnr&&i<n;j++)
					{
						d[i] = A->value[bj*bs + j*bnr + k];
						i++;
						k++;
					}
				}
				if( k==bnr ) break;
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_bsr"
LIS_INT lis_matrix_shift_diagonal_bsr(LIS_MATRIX A, LIS_SCALAR sigma)
{
	LIS_INT i,j,k,bi,bj,bjj,nr;
	LIS_INT bnr,bnc,bs;
	LIS_INT n;

	LIS_DEBUG_FUNC_IN;

	n   = A->n;
	nr  = A->nr;
	bnr = A->bnr;
	bnc = A->bnc;
	bs  = bnr*bnc;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		for(i=0;i<nr;i++)
		{
			for(j=0;j<bnr;j++)
			{
				A->D->value[i*bs+j*bnr+j] -= sigma;
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj,bjj,i,j,k)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			k = 0;
			i = bi*bnr;
			for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
			{
				bjj = A->bindex[bj];
				if( i>=bjj*bnc && i<(bjj+1)*bnc )
				{
					for(j=i%bnc;j<bnc&&k<bnr&&i<n;j++)
					{
						A->value[bj*bs + j*bnr + k] -= sigma;
						i++;
						k++;
					}
				}
				if( k==bnr ) break;
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_bsr"
LIS_INT lis_matrix_scale_bsr(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j;
	LIS_INT bi,bj,bs;
	LIS_INT nr;
	LIS_INT bnr,bnc;

	LIS_DEBUG_FUNC_IN;

	bnr  = A->bnr;
	bnc  = A->bnc;
	nr   = A->nr;
	bs   = A->bnr*A->bnc;

	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj,i,j)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			for(bj=A->L->bptr[bi];bj<A->L->bptr[bi+1];bj++)
			{
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						A->L->value[bj*bs+j*bnr+i] *= d[bi*bnr+i];
					}
				}
			}
			for(bj=A->U->bptr[bi];bj<A->U->bptr[bi+1];bj++)
			{
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						A->U->value[bj*bs+j*bnr+i] *= d[bi*bnr+i];
					}
				}
			}
			for(j=0;j<bnc;j++)
			{
				for(i=0;i<bnr;i++)
				{
					A->D->value[bi*bs+j*bnr+i] *= d[bi*bnr+i];
				}
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj,i,j)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
			{
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						A->value[bj*bs+j*bnr+i] *= d[bi*bnr+i];
					}
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_symm_bsr"
LIS_INT lis_matrix_scale_symm_bsr(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i,j;
	LIS_INT bi,bj,bjj,bs;
	LIS_INT nr;
	LIS_INT bnr,bnc;

	LIS_DEBUG_FUNC_IN;

	bnr  = A->bnr;
	bnc  = A->bnc;
	nr   = A->nr;
	bs   = A->bnr*A->bnc;

	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj,i,j)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			for(bj=A->L->bptr[bi];bj<A->L->bptr[bi+1];bj++)
			{
				bjj = A->L->bindex[bj];
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						A->L->value[bj*bs+j*bnr+i] *= d[bi*bnr+i]*d[bjj*bnc+j];
					}
				}
			}
			for(bj=A->U->bptr[bi];bj<A->U->bptr[bi+1];bj++)
			{
				bjj = A->U->bindex[bj];
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						A->U->value[bj*bs+j*bnr+i] *= d[bi*bnr+i]*d[bjj*bnc+j];
					}
				}
			}
			for(j=0;j<bnc;j++)
			{
				for(i=0;i<bnr;i++)
				{
					A->D->value[bi*bs+j*bnr+i] *= d[bi*bnr+i]*d[bi*bnr+i];
				}
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj,bjj,i,j)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
			{
				bjj = A->bindex[bj];
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						A->value[bj*bs+j*bnr+i] *= d[bi*bnr+i]*d[bjj*bnc+j];
					}
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_bscale_bsr"
LIS_INT lis_matrix_bscale_bsr(LIS_MATRIX A, LIS_MATRIX_DIAG D)
{
	LIS_INT bi,bj;
	LIS_INT nr;
	LIS_INT bnr;
	LIS_SCALAR *d,b[16];

	LIS_DEBUG_FUNC_IN;

	bnr  = A->bnr;
	nr   = A->nr;
	d    = D->value;

	if( bnr==1 )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			A->D->value[bi] = 1.0;
			for(bj=A->L->bptr[bi];bj<A->L->bptr[bi+1];bj++)
			{
				A->L->value[bj] *= d[bi];
			}
			for(bj=A->U->bptr[bi];bj<A->U->bptr[bi+1];bj++)
			{
				A->U->value[bj] *= d[bi];
			}
		}
	}
	else if( bnr==2 )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			A->D->value[4*bi]   = 1.0;
			A->D->value[4*bi+1] = 0.0;
			A->D->value[4*bi+2] = 0.0;
			A->D->value[4*bi+3] = 1.0;
			for(bj=A->L->bptr[bi];bj<A->L->bptr[bi+1];bj++)
			{
				b[0] = d[4*bi+0] * A->L->value[4*bj+0] + d[4*bi+2] * A->L->value[4*bj+1];
				b[1] = d[4*bi+1] * A->L->value[4*bj+0] + d[4*bi+3] * A->L->value[4*bj+1];
				b[2] = d[4*bi+0] * A->L->value[4*bj+2] + d[4*bi+2] * A->L->value[4*bj+3];
				b[3] = d[4*bi+1] * A->L->value[4*bj+2] + d[4*bi+3] * A->L->value[4*bj+3];
				memcpy(&A->L->value[4*bj],b,4*sizeof(LIS_SCALAR));
			}
			for(bj=A->U->bptr[bi];bj<A->U->bptr[bi+1];bj++)
			{
				b[0] = d[4*bi+0] * A->U->value[4*bj+0] + d[4*bi+2] * A->U->value[4*bj+1];
				b[1] = d[4*bi+1] * A->U->value[4*bj+0] + d[4*bi+3] * A->U->value[4*bj+1];
				b[2] = d[4*bi+0] * A->U->value[4*bj+2] + d[4*bi+2] * A->U->value[4*bj+3];
				b[3] = d[4*bi+1] * A->U->value[4*bj+2] + d[4*bi+3] * A->U->value[4*bj+3];
				memcpy(&A->U->value[4*bj],b,4*sizeof(LIS_SCALAR));
			}
		}
	}
	else if( bnr==3 )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			A->D->value[9*bi]   = 1.0;
			A->D->value[9*bi+1] = 0.0;
			A->D->value[9*bi+2] = 0.0;
			A->D->value[9*bi+3] = 0.0;
			A->D->value[9*bi+4] = 1.0;
			A->D->value[9*bi+5] = 0.0;
			A->D->value[9*bi+6] = 0.0;
			A->D->value[9*bi+7] = 0.0;
			A->D->value[9*bi+8] = 1.0;
			for(bj=A->L->bptr[bi];bj<A->L->bptr[bi+1];bj++)
			{
				b[0] = d[9*bi+0] * A->L->value[9*bj+0] + d[9*bi+3] * A->L->value[9*bj+1] + d[9*bi+6] * A->L->value[9*bj+2];
				b[1] = d[9*bi+1] * A->L->value[9*bj+0] + d[9*bi+4] * A->L->value[9*bj+1] + d[9*bi+7] * A->L->value[9*bj+2];
				b[2] = d[9*bi+2] * A->L->value[9*bj+0] + d[9*bi+5] * A->L->value[9*bj+1] + d[9*bi+8] * A->L->value[9*bj+2];
				b[3] = d[9*bi+0] * A->L->value[9*bj+3] + d[9*bi+3] * A->L->value[9*bj+4] + d[9*bi+6] * A->L->value[9*bj+5];
				b[4] = d[9*bi+1] * A->L->value[9*bj+3] + d[9*bi+4] * A->L->value[9*bj+4] + d[9*bi+7] * A->L->value[9*bj+5];
				b[5] = d[9*bi+2] * A->L->value[9*bj+3] + d[9*bi+5] * A->L->value[9*bj+4] + d[9*bi+8] * A->L->value[9*bj+5];
				b[6] = d[9*bi+0] * A->L->value[9*bj+6] + d[9*bi+3] * A->L->value[9*bj+7] + d[9*bi+6] * A->L->value[9*bj+8];
				b[7] = d[9*bi+1] * A->L->value[9*bj+6] + d[9*bi+4] * A->L->value[9*bj+7] + d[9*bi+7] * A->L->value[9*bj+8];
				b[8] = d[9*bi+2] * A->L->value[9*bj+6] + d[9*bi+5] * A->L->value[9*bj+7] + d[9*bi+8] * A->L->value[9*bj+8];
				memcpy(&A->L->value[9*bj],b,8*sizeof(LIS_SCALAR));
			}
			for(bj=A->U->bptr[bi];bj<A->U->bptr[bi+1];bj++)
			{
				b[0] = d[9*bi+0] * A->U->value[9*bj+0] + d[9*bi+3] * A->U->value[9*bj+1] + d[9*bi+6] * A->U->value[9*bj+2];
				b[1] = d[9*bi+1] * A->U->value[9*bj+0] + d[9*bi+4] * A->U->value[9*bj+1] + d[9*bi+7] * A->U->value[9*bj+2];
				b[2] = d[9*bi+2] * A->U->value[9*bj+0] + d[9*bi+5] * A->U->value[9*bj+1] + d[9*bi+8] * A->U->value[9*bj+2];
				b[3] = d[9*bi+0] * A->U->value[9*bj+3] + d[9*bi+3] * A->U->value[9*bj+4] + d[9*bi+6] * A->U->value[9*bj+5];
				b[4] = d[9*bi+1] * A->U->value[9*bj+3] + d[9*bi+4] * A->U->value[9*bj+4] + d[9*bi+7] * A->U->value[9*bj+5];
				b[5] = d[9*bi+2] * A->U->value[9*bj+3] + d[9*bi+5] * A->U->value[9*bj+4] + d[9*bi+8] * A->U->value[9*bj+5];
				b[6] = d[9*bi+0] * A->U->value[9*bj+6] + d[9*bi+3] * A->U->value[9*bj+7] + d[9*bi+6] * A->U->value[9*bj+8];
				b[7] = d[9*bi+1] * A->U->value[9*bj+6] + d[9*bi+4] * A->U->value[9*bj+7] + d[9*bi+7] * A->U->value[9*bj+8];
				b[8] = d[9*bi+2] * A->U->value[9*bj+6] + d[9*bi+5] * A->U->value[9*bj+7] + d[9*bi+8] * A->U->value[9*bj+8];
				memcpy(&A->U->value[9*bj],b,8*sizeof(LIS_SCALAR));
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_normf_bsr"
LIS_INT lis_matrix_normf_bsr(LIS_MATRIX A, LIS_SCALAR *nrm)
{
	LIS_INT j;
	LIS_INT bi,bj,bs;
	LIS_INT nr;
	LIS_INT bnr,bnc;
	LIS_SCALAR sum;

	LIS_DEBUG_FUNC_IN;

	bnr  = A->bnr;
	bnc  = A->bnc;
	bs   = bnr*bnc;
	sum  = (LIS_SCALAR)0;

	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for reduction(+:sum) private(bi,bj,j)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			for(bj=A->L->bptr[bi];bj<A->L->bptr[bi+1];bj++)
			{
				for(j=0;j<bs;j++)
				{
					sum += A->L->value[bj+j]*A->L->value[bj+j];
				}
			}
			for(bj=A->U->bptr[bi];bj<A->U->bptr[bi+1];bj++)
			{
				for(j=0;j<bs;j++)
				{
					sum += A->U->value[bj+j]*A->U->value[bj+j];
				}
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for reduction(+:sum) private(bi,bj,j)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
			{
				for(j=0;j<bs;j++)
				{
					sum += A->value[bj+j]*A->value[bj+j];
				}
			}
		}
	}
	*nrm = sqrt(sum);
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_bsr"
LIS_INT lis_matrix_split_bsr(LIS_MATRIX A)
{
	LIS_INT i,j,n;
	LIS_INT bnr,bnc,nr,nc,bs;
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
	bnr      = A->bnr;
	bnc      = A->bnc;
	nr       = A->nr;
	nc       = A->nc;
	bs       = A->bnr*A->bnc;
	nnzl     = 0;
	nnzu     = 0;
	D        = NULL;
	lptr     = NULL;
	lindex   = NULL;
	lvalue   = NULL;
	uptr     = NULL;
	uindex   = NULL;
	uvalue   = NULL;

	if( bnr!=bnc )
	{
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
	}
	#ifdef _OPENMP
		liw = (LIS_INT *)lis_malloc((nr+1)*sizeof(LIS_INT),"lis_matrix_split_bsr::liw");
		if( liw==NULL )
		{
			LIS_SETERR_MEM((nr+1)*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		uiw = (LIS_INT *)lis_malloc((nr+1)*sizeof(LIS_INT),"lis_matrix_split_bsr::uiw");
		if( uiw==NULL )
		{
			LIS_SETERR_MEM((nr+1)*sizeof(LIS_INT));
			lis_free(liw);
			return LIS_OUT_OF_MEMORY;
		}
		#pragma omp parallel for private(i)
		for(i=0;i<nr+1;i++)
		{
			liw[i] = 0;
			uiw[i] = 0;
		}
		#pragma omp parallel for private(i,j)
		for(i=0;i<nr;i++)
		{
			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				if( A->bindex[j]<i )
				{
					liw[i+1]++;
				}
				else if( A->bindex[j]>i )
				{
					uiw[i+1]++;
				}
			}
		}
		for(i=0;i<nr;i++)
		{
			liw[i+1] += liw[i];
			uiw[i+1] += uiw[i];
		}
		nnzl = liw[nr];
		nnzu = uiw[nr];
	#else
		for(i=0;i<nr;i++)
		{
			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				if( A->bindex[j]<i )
				{
					nnzl++;
				}
				else if( A->bindex[j]>i )
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
	err = lis_matrix_malloc_bsr(n,bnr,bnc,nnzl,&lptr,&lindex,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_bsr(n,bnr,bnc,nnzu,&uptr,&uindex,&uvalue);
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
		for(i=0;i<nr+1;i++)
		{
			lptr[i] = liw[i];
			uptr[i] = uiw[i];
		}
		#pragma omp parallel for private(i,j,kl,ku)
		for(i=0;i<nr;i++)
		{
			kl = lptr[i];
			ku = uptr[i];
			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				if( A->bindex[j]<i )
				{
					lindex[kl]   = A->bindex[j];
					memcpy(&lvalue[bs*kl],&A->value[bs*j],bs*sizeof(LIS_SCALAR));;
					kl++;
				}
				else if( A->bindex[j]>i )
				{
					uindex[ku]   = A->bindex[j];
					memcpy(&uvalue[bs*ku],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
					ku++;
				}
				else
				{
					memcpy(&D->value[bs*i],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
				}
			}
		}
		lis_free2(2,liw,uiw);
	#else
		nnzl = 0;
		nnzu = 0;
		lptr[0] = 0;
		uptr[0] = 0;
		for(i=0;i<nr;i++)
		{
			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				if( A->bindex[j]<i )
				{
					lindex[nnzl]   = A->bindex[j];
					memcpy(&lvalue[bs*nnzl],&A->value[bs*j],bs*sizeof(LIS_SCALAR));;
					nnzl++;
				}
				else if( A->bindex[j]>i )
				{
					uindex[nnzu]   = A->bindex[j];
					memcpy(&uvalue[bs*nnzu],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
					nnzu++;
				}
				else
				{
					memcpy(&D->value[bs*i],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
				}
			}
			lptr[i+1] = nnzl;
			uptr[i+1] = nnzu;
		}
	#endif
	A->L->bnr     = bnr;
	A->L->bnc     = bnc;
	A->L->nr      = nr;
	A->L->nc      = nc;
	A->L->bnnz    = nnzl;
	A->L->bptr    = lptr;
	A->L->bindex  = lindex;
	A->L->value   = lvalue;
	A->U->bnr     = bnr;
	A->U->bnc     = bnc;
	A->U->nr      = nr;
	A->U->nc      = nc;
	A->U->bnnz    = nnzu;
	A->U->bptr    = uptr;
	A->U->bindex  = uindex;
	A->U->value   = uvalue;
	A->D          = D;
	A->is_splited = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_merge_bsr"
LIS_INT lis_matrix_merge_bsr(LIS_MATRIX A)
{
	LIS_INT i,j,n,nr;
	LIS_INT bnnz,bnr,bnc,bs;
	LIS_INT err;
	LIS_INT *bptr,*bindex;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;


	n       = A->n;
	nr      = A->nr;
	bnr     = A->bnr;
	bnc     = A->bnc;
	bs      = bnr*bnc;
	bptr    = NULL;
	bindex  = NULL;
	value   = NULL;
	bnnz    = A->L->bnnz + A->U->bnnz + nr;

	err = lis_matrix_malloc_bsr(n,bnr,bnc,bnnz,&bptr,&bindex,&value);
	if( err )
	{
		return err;
	}

	bnnz    = 0;
	bptr[0] = 0;
	for(i=0;i<nr;i++)
	{
		for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
		{
			bindex[bnnz]   = A->L->bindex[j];
			memcpy(&value[bs*bnnz],&A->L->value[bs*j],bs*sizeof(LIS_SCALAR));;
			bnnz++;
		}
		bindex[bnnz] = i;
		memcpy(&value[bs*bnnz],&A->D->value[bs*i],bs*sizeof(LIS_SCALAR));;
		bnnz++;
		for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
		{
			bindex[bnnz]   = A->U->bindex[j];
			memcpy(&value[bs*bnnz],&A->U->value[bs*j],bs*sizeof(LIS_SCALAR));;
			bnnz++;
		}
		bptr[i+1] = bnnz;
	}

	A->bnnz       = bnnz;
	A->bptr       = bptr;
	A->value      = value;
	A->bindex      = bindex;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solve_bsr"
LIS_INT lis_matrix_solve_bsr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,k,ii,jj,nr,bnr,bnc,bs;
	LIS_SCALAR t0,t1,t2;
	LIS_SCALAR *b,*x,*w;

	LIS_DEBUG_FUNC_IN;

	nr  = A->nr;
	bnr = A->bnr;
	bnc = A->bnc;
	bs  = A->bnr*A->bnc;
	b   = B->value;
	x   = X->value;

	
	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		switch(bnr)
		{
		case 1:
			for(i=0;i<nr;i++)
			{
				t0 = b[i];
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj  = A->L->bindex[j];
					t0 -= A->L->value[j] * x[jj];
				}
				x[i] = A->WD->value[i] * t0;
			}
			break;
		case 2:
			for(i=0;i<nr;i++)
			{
				t0 = b[i*2];
				t1 = b[i*2+1];
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj  = A->L->bindex[j];
					t0 -= A->L->value[j*4+0] * x[jj*2+0];
					t1 -= A->L->value[j*4+1] * x[jj*2+0];
					t0 -= A->L->value[j*4+2] * x[jj*2+1];
					t1 -= A->L->value[j*4+3] * x[jj*2+1];
				}
				x[i*2+0] = A->WD->value[4*i+0] * t0 + A->WD->value[4*i+2] * t1;
				x[i*2+1] = A->WD->value[4*i+1] * t0 + A->WD->value[4*i+3] * t1;
			}
			break;
		case 3:
			for(i=0;i<nr;i++)
			{
				t0 = b[i*3];
				t1 = b[i*3+1];
				t2 = b[i*3+2];
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj  = A->L->bindex[j];
					t0 -= A->L->value[j*9+0] * x[jj*3+0];
					t1 -= A->L->value[j*9+1] * x[jj*3+0];
					t2 -= A->L->value[j*9+2] * x[jj*3+0];
					t0 -= A->L->value[j*9+3] * x[jj*3+1];
					t1 -= A->L->value[j*9+4] * x[jj*3+1];
					t2 -= A->L->value[j*9+5] * x[jj*3+1];
					t0 -= A->L->value[j*9+6] * x[jj*3+2];
					t1 -= A->L->value[j*9+7] * x[jj*3+2];
					t2 -= A->L->value[j*9+8] * x[jj*3+2];
				}
				x[i*3+0] = A->WD->value[9*i+0] * t0 + A->WD->value[9*i+3] * t1 + A->WD->value[9*i+6] * t2;
				x[i*3+1] = A->WD->value[9*i+1] * t0 + A->WD->value[9*i+4] * t1 + A->WD->value[9*i+7] * t2;
				x[i*3+2] = A->WD->value[9*i+2] * t0 + A->WD->value[9*i+5] * t1 + A->WD->value[9*i+8] * t2;
			}
			break;
		default:
			w = (LIS_SCALAR *)lis_malloc(bnr*sizeof(LIS_SCALAR),"lis_matrix_solve_bsr::w");
			for(i=0;i<nr;i++)
			{
				for(j=0;j<bnr;j++)
				{
					w[j] = b[i*bnr+j];
				}
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					k   = A->L->bindex[j] * bnc;
					for(ii=0;ii<bnr;ii++)
					{
						t0   = w[ii];
						for(jj=0;jj<bnc;jj++)
						{
							t0 -= A->L->value[j*bs + jj*bnr+ii] * x[k + jj];
						}
						w[ii] = t0;
					}
				}
				for(ii=0;ii<bnr;ii++)
				{
					t0 = 0.0;
					for(jj=0;jj<bnc;jj++)
					{
						t0 += A->WD->value[i*bs + jj*bnr+ii] * w[jj];
					}
					x[i*bnr+ii] = t0;
				}
			}
			lis_free(w);
			break;
		}
		break;
	case LIS_MATRIX_UPPER:
		switch(bnr)
		{
		case 1:
			for(i=nr-1;i>=0;i--)
			{
				t0 = b[i];
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj  = A->U->bindex[j];
					t0 -= A->U->value[j] * x[jj];
				}
				x[i] = A->WD->value[i] * t0;
			}
			break;
		case 2:
			for(i=nr-1;i>=0;i--)
			{
				t0 = b[i*2];
				t1 = b[i*2+1];
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj  = A->U->bindex[j];
					t0 -= A->U->value[j*4+0] * x[jj*2+0];
					t1 -= A->U->value[j*4+1] * x[jj*2+0];
					t0 -= A->U->value[j*4+2] * x[jj*2+1];
					t1 -= A->U->value[j*4+3] * x[jj*2+1];
				}
				x[i*2+0] = A->WD->value[4*i+0] * t0 + A->WD->value[4*i+2] * t1;
				x[i*2+1] = A->WD->value[4*i+1] * t0 + A->WD->value[4*i+3] * t1;
			}
			break;
		case 3:
			for(i=nr-1;i>=0;i--)
			{
				t0 = b[i*3];
				t1 = b[i*3+1];
				t2 = b[i*3+2];
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj  = A->U->bindex[j];
					t0 -= A->U->value[j*9+0] * x[jj*3+0];
					t1 -= A->U->value[j*9+1] * x[jj*3+0];
					t2 -= A->U->value[j*9+2] * x[jj*3+0];
					t0 -= A->U->value[j*9+3] * x[jj*3+1];
					t1 -= A->U->value[j*9+4] * x[jj*3+1];
					t2 -= A->U->value[j*9+5] * x[jj*3+1];
					t0 -= A->U->value[j*9+6] * x[jj*3+2];
					t1 -= A->U->value[j*9+7] * x[jj*3+2];
					t2 -= A->U->value[j*9+8] * x[jj*3+2];
				}
				x[i*3+0] = A->WD->value[9*i+0] * t0 + A->WD->value[9*i+3] * t1 + A->WD->value[9*i+6] * t2;
				x[i*3+1] = A->WD->value[9*i+1] * t0 + A->WD->value[9*i+4] * t1 + A->WD->value[9*i+7] * t2;
				x[i*3+2] = A->WD->value[9*i+2] * t0 + A->WD->value[9*i+5] * t1 + A->WD->value[9*i+8] * t2;
			}
			break;
		default:
			w = (LIS_SCALAR *)lis_malloc(bnr*sizeof(LIS_SCALAR),"lis_matrix_solve_bsr::w");
			for(i=nr-1;i>=0;i--)
			{
				for(j=0;j<bnr;j++)
				{
					w[j] = b[i*bnr+j];
				}
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					k   = A->U->bindex[j] * bnc;
					for(ii=0;ii<bnr;ii++)
					{
						t0   = w[ii];
						for(jj=0;jj<bnc;jj++)
						{
							t0 -= A->U->value[j*bs + jj*bnr+ii] * x[k + jj];
						}
						w[ii] = t0;
					}
				}
				for(ii=0;ii<bnr;ii++)
				{
					t0 = 0.0;
					for(jj=0;jj<bnc;jj++)
					{
						t0 += A->WD->value[i*bs + jj*bnr+ii] * w[jj];
					}
					x[i*bnr+ii] = t0;
				}
			}
			lis_free(w);
			break;
		}
		break;
	case LIS_MATRIX_SSOR:
		switch(bnr)
		{
		case 1:
			for(i=0;i<nr;i++)
			{
				t0 = b[i];
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj  = A->L->bindex[j];
					t0 -= A->L->value[j] * x[jj];
				}
				x[i] = A->WD->value[i] * t0;
			}
			for(i=nr-1;i>=0;i--)
			{
				t0 = 0.0;
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj  = A->U->bindex[j];
					t0 += A->U->value[j] * x[jj];
				}
				x[i] -= A->WD->value[i] * t0;
			}
			break;
		case 2:
			for(i=0;i<nr;i++)
			{
				t0 = b[i*2];
				t1 = b[i*2+1];
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj  = A->L->bindex[j];
					t0 -= A->L->value[j*4+0] * x[jj*2+0];
					t1 -= A->L->value[j*4+1] * x[jj*2+0];
					t0 -= A->L->value[j*4+2] * x[jj*2+1];
					t1 -= A->L->value[j*4+3] * x[jj*2+1];
				}
				x[i*2+0] = A->WD->value[4*i+0] * t0 + A->WD->value[4*i+2] * t1;
				x[i*2+1] = A->WD->value[4*i+1] * t0 + A->WD->value[4*i+3] * t1;
			}
			for(i=nr-1;i>=0;i--)
			{
				t0 = 0.0;
				t1 = 0.0;
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj  = A->U->bindex[j];
					t0 += A->U->value[j*4+0] * x[jj*2+0];
					t1 += A->U->value[j*4+1] * x[jj*2+0];
					t0 += A->U->value[j*4+2] * x[jj*2+1];
					t1 += A->U->value[j*4+3] * x[jj*2+1];
				}
				x[i*2+0] -= A->WD->value[4*i+0] * t0 + A->WD->value[4*i+2] * t1;
				x[i*2+1] -= A->WD->value[4*i+1] * t0 + A->WD->value[4*i+3] * t1;
			}
			break;
		case 3:
			for(i=0;i<nr;i++)
			{
				t0 = b[i*bnr];
				t1 = b[i*bnr+1];
				t2 = b[i*bnr+2];
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj  = A->L->bindex[j];
					t0 -= A->L->value[j*9+0] * x[jj*3+0];
					t1 -= A->L->value[j*9+1] * x[jj*3+0];
					t2 -= A->L->value[j*9+2] * x[jj*3+0];
					t0 -= A->L->value[j*9+3] * x[jj*3+1];
					t1 -= A->L->value[j*9+4] * x[jj*3+1];
					t2 -= A->L->value[j*9+5] * x[jj*3+1];
					t0 -= A->L->value[j*9+6] * x[jj*3+2];
					t1 -= A->L->value[j*9+7] * x[jj*3+2];
					t2 -= A->L->value[j*9+8] * x[jj*3+2];
				}
				x[i*bnr+0] = A->WD->value[9*i+0] * t0 + A->WD->value[9*i+3] * t1 + A->WD->value[9*i+6] * t2;
				x[i*bnr+1] = A->WD->value[9*i+1] * t0 + A->WD->value[9*i+4] * t1 + A->WD->value[9*i+7] * t2;
				x[i*bnr+2] = A->WD->value[9*i+2] * t0 + A->WD->value[9*i+5] * t1 + A->WD->value[9*i+8] * t2;
			}
			for(i=nr-1;i>=0;i--)
			{
				t0 = 0.0;
				t1 = 0.0;
				t2 = 0.0;
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj  = A->U->bindex[j];
					t0 += A->U->value[j*9+0] * x[jj*3+0];
					t1 += A->U->value[j*9+1] * x[jj*3+0];
					t2 += A->U->value[j*9+2] * x[jj*3+0];
					t0 += A->U->value[j*9+3] * x[jj*3+1];
					t1 += A->U->value[j*9+4] * x[jj*3+1];
					t2 += A->U->value[j*9+5] * x[jj*3+1];
					t0 += A->U->value[j*9+6] * x[jj*3+2];
					t1 += A->U->value[j*9+7] * x[jj*3+2];
					t2 += A->U->value[j*9+8] * x[jj*3+2];
				}
				x[i*3+0] -= A->WD->value[9*i+0] * t0 + A->WD->value[9*i+3] * t1 + A->WD->value[9*i+6] * t2;
				x[i*3+1] -= A->WD->value[9*i+1] * t0 + A->WD->value[9*i+4] * t1 + A->WD->value[9*i+7] * t2;
				x[i*3+2] -= A->WD->value[9*i+2] * t0 + A->WD->value[9*i+5] * t1 + A->WD->value[9*i+8] * t2;
			}
			break;
		default:
			w = (LIS_SCALAR *)lis_malloc(bnr*sizeof(LIS_SCALAR),"lis_matrix_solve_bsr::w");
			for(i=0;i<nr;i++)
			{
				for(j=0;j<bnr;j++)
				{
					w[j] = b[i*bnr+j];
				}
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					k   = A->L->bindex[j] * bnc;
					for(ii=0;ii<bnr;ii++)
					{
						t0   = w[ii];
						for(jj=0;jj<bnc;jj++)
						{
							t0 -= A->L->value[j*bs + jj*bnr+ii] * x[k + jj];
						}
						w[ii] = t0;
					}
				}
				for(ii=0;ii<bnr;ii++)
				{
					t0 = 0.0;
					for(jj=0;jj<bnc;jj++)
					{
						t0 += A->WD->value[i*bs + jj*bnr+ii] * w[jj];
					}
					x[i*bnr+ii] = t0;
				}
			}
			for(i=nr-1;i>=0;i--)
			{
				for(j=0;j<bnr;j++)
				{
					w[j] = 0.0;
				}
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					k   = A->U->bindex[j] * bnc;
					for(ii=0;ii<bnr;ii++)
					{
						t0   = w[ii];
						for(jj=0;jj<bnc;jj++)
						{
							t0 += A->U->value[j*bs + jj*bnr+ii] * x[k + jj];
						}
						w[ii] = t0;
					}
				}
				for(ii=0;ii<bnr;ii++)
				{
					t0 = 0.0;
					for(jj=0;jj<bnc;jj++)
					{
						t0 += A->WD->value[i*bs + jj*bnr+ii] * w[jj];
					}
					x[i*bnr+ii] -= t0;
				}
			}
			lis_free(w);
			break;
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solveh_bsr"
LIS_INT lis_matrix_solveh_bsr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag)
{
	LIS_INT i,j,k,ii,jj,nr,bnr,bnc,bs;
	LIS_SCALAR t0,t1,t2;
	LIS_SCALAR *x,*w;

	LIS_DEBUG_FUNC_IN;

	nr  = A->nr;
	bnr = A->bnr;
	bnc = A->bnc;
	bs  = A->bnr*A->bnc;
	x   = X->value;

	lis_vector_copy(B,X);
	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		switch(bnr)
		{
		case 1:
			for(i=0;i<nr;i++)
			{
				x[i] = x[i] * conj(A->WD->value[i]);
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj     = A->U->bindex[j];
					x[jj] -= conj(A->U->value[j]) * x[i];
				}
			}
			break;
		case 2:
			for(i=0;i<nr;i++)
			{
				t0 = conj(A->WD->value[4*i+0]) * x[i*2] + conj(A->WD->value[4*i+1]) * x[i*2+1];
				t1 = conj(A->WD->value[4*i+2]) * x[i*2] + conj(A->WD->value[4*i+3]) * x[i*2+1];
				x[i*2+0] = t0;
				x[i*2+1] = t1;
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj  = A->U->bindex[j];
					x[jj*2+0] -= conj(A->U->value[j*4+0]) * t0 + A->U->value[j*4+1] * t1;
					x[jj*2+1] -= conj(A->U->value[j*4+2]) * t0 + A->U->value[j*4+3] * t1;
				}
			}
			break;
		case 3:
			for(i=0;i<nr;i++)
			{
				t0 = conj(A->WD->value[9*i+0]) * x[i*3] + conj(A->WD->value[9*i+1]) * x[i*3+1] + conj(A->WD->value[9*i+2]) * x[i*3+2];
				t1 = conj(A->WD->value[9*i+3]) * x[i*3] + conj(A->WD->value[9*i+4]) * x[i*3+1] + conj(A->WD->value[9*i+5]) * x[i*3+2];
				t2 = conj(A->WD->value[9*i+6]) * x[i*3] + conj(A->WD->value[9*i+7]) * x[i*3+1] + conj(A->WD->value[9*i+8]) * x[i*3+2];
				x[i*3]   = t0;
				x[i*3+1] = t1;
				x[i*3+2] = t2;
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj  = A->U->bindex[j];
					x[jj*3+0] -= conj(A->U->value[j*9+0]) * t0 + A->U->value[j*9+1] * t1 + A->U->value[j*9+2] * t2;
					x[jj*3+1] -= conj(A->U->value[j*9+3]) * t0 + A->U->value[j*9+4] * t1 + A->U->value[j*9+5] * t2;
					x[jj*3+2] -= conj(A->U->value[j*9+6]) * t0 + A->U->value[j*9+7] * t1 + A->U->value[j*9+8] * t2;
				}
			}
			break;
		default:
			w = (LIS_SCALAR *)lis_malloc(bnc*sizeof(LIS_SCALAR),"lis_matrix_solve_bsr::w");
			for(i=0;i<nr;i++)
			{
				for(jj=0;jj<bnc;jj++)
				{
					t0 = 0.0;
					for(ii=0;ii<bnr;ii++)
					{
						t0 += conj(A->WD->value[i*bs + jj*bnr+ii]) * x[i*bnr + ii];
					}
					w[jj] = t0;
				}
				memcpy(&x[i*bnr],w,bnr*sizeof(LIS_SCALAR));
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					k   = A->U->bindex[j] * bnc;
					for(jj=0;jj<bnc;jj++)
					{
						t0 = 0.0;
						for(ii=0;ii<bnr;ii++)
						{
							t0 += conj(A->U->value[j*bs + jj*bnr+ii]) * w[ii];
						}
						x[k + jj] -= t0;
					}
				}
			}
			lis_free(w);
			break;
		}
		break;
	case LIS_MATRIX_UPPER:
		switch(bnr)
		{
		case 1:
			for(i=nr-1;i>=0;i--)
			{
				x[i] = x[i] * A->WD->value[i];
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj     = A->L->bindex[j];
					x[jj] -= conj(A->L->value[j]) * x[i];
				}
			}
			break;
		case 2:
			for(i=nr-1;i>=0;i--)
			{
				t0 = conj(A->WD->value[4*i+0]) * x[i*2] + conj(A->WD->value[4*i+1]) * x[i*2+1];
				t1 = conj(A->WD->value[4*i+2]) * x[i*2] + conj(A->WD->value[4*i+3]) * x[i*2+1];
				x[i*2+0] = t0;
				x[i*2+1] = t1;
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj  = A->L->bindex[j];
					x[jj*2+0] -= conj(A->L->value[j*4+0]) * t0 + A->L->value[j*4+1] * t1;
					x[jj*2+1] -= conj(A->L->value[j*4+2]) * t0 + A->L->value[j*4+3] * t1;
				}
			}
			break;
		case 3:
			for(i=nr-1;i>=0;i--)
			{
				t0 = conj(A->WD->value[9*i+0]) * x[i*3] + conj(A->WD->value[9*i+1]) * x[i*3+1] + conj(A->WD->value[9*i+2]) * x[i*3+2];
				t1 = conj(A->WD->value[9*i+3]) * x[i*3] + conj(A->WD->value[9*i+4]) * x[i*3+1] + conj(A->WD->value[9*i+5]) * x[i*3+2];
				t2 = conj(A->WD->value[9*i+6]) * x[i*3] + conj(A->WD->value[9*i+7]) * x[i*3+1] + conj(A->WD->value[9*i+8]) * x[i*3+2];
				x[i*3]   = t0;
				x[i*3+1] = t1;
				x[i*3+2] = t2;
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj  = A->L->bindex[j];
					x[jj*3+0] -= conj(A->L->value[j*9+0]) * t0 + conj(A->L->value[j*9+1]) * t1 + conj(A->L->value[j*9+2]) * t2;
					x[jj*3+1] -= conj(A->L->value[j*9+3]) * t0 + conj(A->L->value[j*9+4]) * t1 + conj(A->L->value[j*9+5]) * t2;
					x[jj*3+2] -= conj(A->L->value[j*9+6]) * t0 + conj(A->L->value[j*9+7]) * t1 + conj(A->L->value[j*9+8]) * t2;
				}
			}
			break;
		default:
			w = (LIS_SCALAR *)lis_malloc(bnr*sizeof(LIS_SCALAR),"lis_matrix_solve_bsr::w");
			for(i=nr-1;i>=0;i--)
			{
				for(jj=0;jj<bnc;jj++)
				{
					t0 = 0.0;
					for(ii=0;ii<bnr;ii++)
					{
						t0 += conj(A->WD->value[i*bs + jj*bnr+ii]) * x[i*bnr + ii];
					}
					w[jj] = t0;
				}
				memcpy(&x[i*bnr],w,bnr*sizeof(LIS_SCALAR));
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					k   = A->L->bindex[j] * bnc;
					for(jj=0;jj<bnc;jj++)
					{
						t0 = 0.0;
						for(ii=0;ii<bnr;ii++)
						{
							t0 += conj(A->L->value[j*bs + jj*bnr+ii]) * w[ii];
						}
						x[k + jj] -= t0;
					}
				}
			}
			lis_free(w);
			break;
		}
		break;
	case LIS_MATRIX_SSOR:
		switch(bnr)
		{
		case 1:
			for(i=0;i<nr;i++)
			{
				t0 = x[i] * conj(A->WD->value[i]);
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj     = A->U->bindex[j];
					x[jj] -= conj(A->U->value[j]) * t0;
				}
			}
			for(i=nr-1;i>=0;i--)
			{
				t0   = x[i] * conj(A->WD->value[i]);
				x[i] = t0;
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj     = A->L->bindex[j];
					x[jj] -= conj(A->L->value[j]) * t0;
				}
			}
			break;
		case 2:
			for(i=0;i<nr;i++)
			{
				t0 = conj(A->WD->value[4*i+0]) * x[i*2] + conj(A->WD->value[4*i+1]) * x[i*2+1];
				t1 = conj(A->WD->value[4*i+2]) * x[i*2] + conj(A->WD->value[4*i+3]) * x[i*2+1];
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj  = A->U->bindex[j];
					x[jj*2+0] -= conj(A->U->value[j*4+0]) * t0 + conj(A->U->value[j*4+1]) * t1;
					x[jj*2+1] -= conj(A->U->value[j*4+2]) * t0 + conj(A->U->value[j*4+3]) * t1;
				}
			}
			for(i=nr-1;i>=0;i--)
			{
				t0 = conj(A->WD->value[4*i+0]) * x[i*2] + conj(A->WD->value[4*i+1]) * x[i*2+1];
				t1 = conj(A->WD->value[4*i+2]) * x[i*2] + conj(A->WD->value[4*i+3]) * x[i*2+1];
				x[i*2+0] = t0;
				x[i*2+1] = t1;
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj  = A->L->bindex[j];
					x[jj*2+0] -= conj(A->L->value[j*4+0]) * t0 + conj(A->L->value[j*4+1]) * t1;
					x[jj*2+1] -= conj(A->L->value[j*4+2]) * t0 + conj(A->L->value[j*4+3]) * t1;
				}
			}
			break;
		case 3:
			for(i=0;i<nr;i++)
			{
				t0 = conj(A->WD->value[9*i+0]) * x[i*3] + conj(A->WD->value[9*i+1]) * x[i*3+1] + conj(A->WD->value[9*i+2]) * x[i*3+2];
				t1 = conj(A->WD->value[9*i+3]) * x[i*3] + conj(A->WD->value[9*i+4]) * x[i*3+1] + conj(A->WD->value[9*i+5]) * x[i*3+2];
				t2 = conj(A->WD->value[9*i+6]) * x[i*3] + conj(A->WD->value[9*i+7]) * x[i*3+1] + conj(A->WD->value[9*i+8]) * x[i*3+2];
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					jj  = A->U->bindex[j];
					x[jj*3+0] -= conj(A->U->value[j*9+0]) * t0 + conj(A->U->value[j*9+1]) * t1 + conj(A->U->value[j*9+2]) * t2;
					x[jj*3+1] -= conj(A->U->value[j*9+3]) * t0 + conj(A->U->value[j*9+4]) * t1 + conj(A->U->value[j*9+5]) * t2;
					x[jj*3+2] -= conj(A->U->value[j*9+6]) * t0 + conj(A->U->value[j*9+7]) * t1 + conj(A->U->value[j*9+8]) * t2;
				}
			}
			for(i=nr-1;i>=0;i--)
			{
				t0 = conj(A->WD->value[9*i+0]) * x[i*3] + conj(A->WD->value[9*i+1]) * x[i*3+1] + conj(A->WD->value[9*i+2]) * x[i*3+2];
				t1 = conj(A->WD->value[9*i+3]) * x[i*3] + conj(A->WD->value[9*i+4]) * x[i*3+1] + conj(A->WD->value[9*i+5]) * x[i*3+2];
				t2 = conj(A->WD->value[9*i+6]) * x[i*3] + conj(A->WD->value[9*i+7]) * x[i*3+1] + conj(A->WD->value[9*i+8]) * x[i*3+2];
				x[i*3]   = t0;
				x[i*3+1] = t1;
				x[i*3+2] = t2;
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					jj  = A->L->bindex[j];
					x[jj*3+0] -= conj(A->L->value[j*9+0]) * t0 + conj(A->L->value[j*9+1]) * t1 + conj(A->L->value[j*9+2]) * t2;
					x[jj*3+1] -= conj(A->L->value[j*9+3]) * t0 + conj(A->L->value[j*9+4]) * t1 + conj(A->L->value[j*9+5]) * t2;
					x[jj*3+2] -= conj(A->L->value[j*9+6]) * t0 + conj(A->L->value[j*9+7]) * t1 + conj(A->L->value[j*9+8]) * t2;
				}
			}
			break;
		default:
			w = (LIS_SCALAR *)lis_malloc(bnc*sizeof(LIS_SCALAR),"lis_matrix_solve_bsr::w");
			for(i=0;i<nr;i++)
			{
				for(jj=0;jj<bnc;jj++)
				{
					t0 = 0.0;
					for(ii=0;ii<bnr;ii++)
					{
						t0 += conj(A->WD->value[i*bs + jj*bnr+ii]) * x[i*bnr + ii];
					}
					w[jj] = t0;
				}
				for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
				{
					k   = A->U->bindex[j] * bnc;
					for(jj=0;jj<bnc;jj++)
					{
						t0 = 0.0;
						for(ii=0;ii<bnr;ii++)
						{
							t0 += conj(A->U->value[j*bs + jj*bnr+ii]) * w[ii];
						}
						x[k + jj] -= t0;
					}
				}
			}
			for(i=nr-1;i>=0;i--)
			{
				for(jj=0;jj<bnc;jj++)
				{
					t0 = 0.0;
					for(ii=0;ii<bnr;ii++)
					{
						t0 += conj(A->WD->value[i*bs + jj*bnr+ii]) * x[i*bnr + ii];
					}
					w[jj] = t0;
				}
				memcpy(&x[i*bnr],w,bnr*sizeof(LIS_SCALAR));
				for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
				{
					k   = A->L->bindex[j] * bnc;
					for(jj=0;jj<bnc;jj++)
					{
						t0 = 0.0;
						for(ii=0;ii<bnr;ii++)
						{
							t0 += conj(A->L->value[j*bs + jj*bnr+ii]) * w[ii];
						}
						x[k + jj] -= t0;
					}
				}
			}
			lis_free(w);
			break;
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_sort_bsr"
LIS_INT lis_matrix_sort_bsr(LIS_MATRIX A)
{
	LIS_INT i,nr,bnr,bs;

	LIS_DEBUG_FUNC_IN;

	if( !A->is_sorted )
	{
		nr  = A->nr;
		bnr = A->bnr;
		bs  = bnr*bnr;
		if( A->is_splited )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<nr;i++)
			{
				lis_sort_id_block(A->L->bptr[i],A->L->bptr[i+1]-1,A->L->bindex,A->L->value,bs);
				lis_sort_id_block(A->U->bptr[i],A->U->bptr[i+1]-1,A->U->bindex,A->U->value,bs);
			}
		}
		else
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<nr;i++)
			{
				lis_sort_id_block(A->bptr[i],A->bptr[i+1]-1,A->bindex,A->value,bs);
			}
		}
		A->is_sorted = LIS_TRUE;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
