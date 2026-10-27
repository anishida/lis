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
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_matrix_set_blocksize
 * lis_matrix_convert
 * lis_matrix_convert_self
 * lis_matrix_copy
 * lis_matrix_copyDLU
 * lis_matrix_axpy
 * lis_matrix_xpay
 * lis_matrix_axpyz
 * lis_matrix_scale
 * lis_matrix_get_diagonal
 * lis_matrix_shift_diagonal
 * lis_matrix_shift_matrix
 * lis_matrix_split
 * lis_matrix_merge
 * lis_matrix_solve
 * lis_matrix_solveh
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_blocksize"
LIS_INT lis_matrix_set_blocksize(LIS_MATRIX A, LIS_INT bnr, LIS_INT bnc, LIS_INT row[], LIS_INT col[])
{
	LIS_INT i,n;
	LIS_INT err;
	LIS_INT *conv_row,*conv_col;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;

	if( bnr<=0 || bnc<=0 )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"bnr=%D <= 0 or bnc=%D <= 0\n",bnr,bnc);
		return LIS_ERR_ILL_ARG;
	}
	if( (row==NULL && col!=NULL) || (row!=NULL && col==NULL) )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"either row[]=%x or col[]=%x is NULL\n",row,col);
		return LIS_ERR_ILL_ARG;
	}
	if( row==NULL )
	{
		A->conv_bnr = bnr;
		A->conv_bnc = bnc;
	}
	else
	{
		n = A->n;
		conv_row = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_set_blocksize::conv_row");
		if( conv_row==NULL )
		{
			LIS_SETERR_MEM(n*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
		conv_col = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_set_blocksize::conv_col");
		if( conv_col==NULL )
		{
			LIS_SETERR_MEM(n*sizeof(LIS_INT));
			lis_free(conv_row);
			return LIS_OUT_OF_MEMORY;
		}
		for(i=0;i<n;i++)
		{
			conv_row[i] = row[i];
			conv_col[i] = col[i];
		}
		A->conv_row = conv_row;
		A->conv_col = conv_col;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_user2csr"
static LIS_INT lis_matrix_convert_user2csr(
	LIS_MATRIX Ain,
	LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT i,n;
	LIS_INT row_nnz;
	LIS_INT total_nnz;
	LIS_INT capacity;
	LIS_INT *ptr;
	LIS_INT *index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	n = Ain->n;

	ptr   = NULL;
	index = NULL;
	value = NULL;

	ptr = (LIS_INT *)lis_malloc(
		(n+1)*sizeof(LIS_INT),
		"lis_matrix_convert_user2csr::ptr");

	if( ptr==NULL )
	{
		LIS_SETERR_MEM(
			(n+1)*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	ptr[0] = 0;

	for(i=0;i<n;i++)
	{
		row_nnz = 0;

		err = lis_matrix_get_user_row(
			Ain,
			i,
			&row_nnz,
			NULL,
			NULL);

		if( err )
		{
			lis_free(ptr);
			return err;
		}

		ptr[i+1] = ptr[i] + row_nnz;
	}

	total_nnz = ptr[n];

	if( total_nnz>0 )
	{
		index = (LIS_INT *)lis_malloc(
			total_nnz*sizeof(LIS_INT),
			"lis_matrix_convert_user2csr::index");

		if( index==NULL )
		{
			LIS_SETERR_MEM(
				total_nnz*sizeof(LIS_INT));
			lis_free(ptr);
			return LIS_OUT_OF_MEMORY;
		}

		value = (LIS_SCALAR *)lis_malloc(
			total_nnz*sizeof(LIS_SCALAR),
			"lis_matrix_convert_user2csr::value");

		if( value==NULL )
		{
			LIS_SETERR_MEM(
				total_nnz*sizeof(LIS_SCALAR));
			lis_free2(2,ptr,index);
			return LIS_OUT_OF_MEMORY;
		}
	}

	for(i=0;i<n;i++)
	{
		row_nnz = ptr[i+1] - ptr[i];

		if( row_nnz==0 )
		{
			continue;
		}

		capacity = row_nnz;

		err = lis_matrix_get_user_row(
			Ain,
			i,
			&capacity,
			&index[ptr[i]],
			&value[ptr[i]]);

		if( err )
		{
			lis_free2(
				3,
				ptr,
				index,
				value);

			return err;
		}

		if( capacity!=row_nnz )
		{
			lis_free2(
				3,
				ptr,
				index,
				value);

			LIS_SETERR(
				LIS_ERR_ILL_ARG,
				"LIS_MATRIX_USER row size changed during conversion\n");

			return LIS_ERR_ILL_ARG;
		}
	}

	err = lis_matrix_set_csr(
		total_nnz,
		ptr,
		index,
		value,
		Aout);

	if( err )
	{
		lis_free2(
			3,
			ptr,
			index,
			value);

		return err;
	}

	/*
	 * USER row callbacks expose global column indices.
	 * Normal matrix assembly performs the standard LIS
	 * global-to-local mapping for distributed matrices.
	 */
	err = lis_matrix_assemble(Aout);
	if( err ) return err;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#undef __FUNC__
#define __FUNC__ "lis_matrix_convert"
LIS_INT lis_matrix_convert(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT istmp;
	LIS_INT convert_matrix_type;
	LIS_MATRIX Atmp,Atmp2;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;
	err = lis_matrix_check(Aout,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;

	err = lis_matrix_merge(Ain);
	if( err ) return err;

	convert_matrix_type = Aout->matrix_type;

	if( Ain->matrix_type==convert_matrix_type && !Ain->is_block )
	{
		err = lis_matrix_copy(Ain,Aout);
		return err;
	}
	if( Ain->matrix_type!=LIS_MATRIX_CSR )
	{
		istmp = LIS_TRUE;
		switch( Ain->matrix_type )
		{
		case LIS_MATRIX_RCO:
			switch( convert_matrix_type )
			{
			case LIS_MATRIX_CSR:
				err = lis_matrix_convert_rco2csr(Ain,Aout);
				LIS_DEBUG_FUNC_OUT;
				return err;
				break;
			case LIS_MATRIX_BSR:
				err = lis_matrix_convert_rco2bsr(Ain,Aout);
				LIS_DEBUG_FUNC_OUT;
				return err;
			case LIS_MATRIX_CSC:
				err = lis_matrix_convert_rco2csc(Ain,Aout);
				LIS_DEBUG_FUNC_OUT;
				return err;
				break;
			default:
				err = lis_matrix_duplicate(Ain,&Atmp);
				if( err ) return err;
				err = lis_matrix_convert_rco2csr(Ain,Atmp);
				break;
			}
			break;
		case LIS_MATRIX_CSC:
			switch( convert_matrix_type )
			{
			case LIS_MATRIX_BSC:
				err = lis_matrix_convert_csc2bsc(Ain,Aout);
				LIS_DEBUG_FUNC_OUT;
				return err;
			default:
				err = lis_matrix_duplicate(Ain,&Atmp);
				if( err ) return err;
				err = lis_matrix_convert_csc2csr(Ain,Atmp);
				break;
			}
			break;
		case LIS_MATRIX_MSR:
			err = lis_matrix_duplicate(Ain,&Atmp);
			if( err ) return err;
			err = lis_matrix_convert_msr2csr(Ain,Atmp);
			break;
		case LIS_MATRIX_DIA:
			err = lis_matrix_duplicate(Ain,&Atmp);
			if( err ) return err;
			err = lis_matrix_convert_dia2csr(Ain,Atmp);
			break;
		case LIS_MATRIX_ELL:
			err = lis_matrix_duplicate(Ain,&Atmp);
			if( err ) return err;
			err = lis_matrix_convert_ell2csr(Ain,Atmp);
			break;
		case LIS_MATRIX_JAD:
			err = lis_matrix_duplicate(Ain,&Atmp);
			if( err ) return err;
			err = lis_matrix_convert_jad2csr(Ain,Atmp);
			break;
		case LIS_MATRIX_BSR:
			err = lis_matrix_duplicate(Ain,&Atmp);
			if( err ) return err;
			err = lis_matrix_convert_bsr2csr(Ain,Atmp);
			break;
		case LIS_MATRIX_BSC:
			err = lis_matrix_duplicate(Ain,&Atmp);
			if( err ) return err;
			err = lis_matrix_convert_bsc2csr(Ain,Atmp);
			break;
		case LIS_MATRIX_VBR:
			err = lis_matrix_duplicate(Ain,&Atmp);
			if( err ) return err;
			err = lis_matrix_convert_vbr2csr(Ain,Atmp);
			break;
		case LIS_MATRIX_DNS:
			err = lis_matrix_duplicate(Ain,&Atmp);
			if( err ) return err;
			err = lis_matrix_convert_dns2csr(Ain,Atmp);
			break;
		case LIS_MATRIX_COO:
			err = lis_matrix_duplicate(Ain,&Atmp);
			if( err ) return err;
			err = lis_matrix_convert_coo2csr(Ain,Atmp);
			break;
		case LIS_MATRIX_USER:
			err = lis_matrix_duplicate(Ain,&Atmp);
			if( err ) return err;
			err = lis_matrix_convert_user2csr(Ain,Atmp);
			break;
		default:
			LIS_SETERR_IMP;
			err = LIS_ERR_NOT_IMPLEMENTED;
		}
		if( err )
		{
			return err;
		}
		if( convert_matrix_type== LIS_MATRIX_CSR )
		{
			lis_matrix_storage_destroy(Aout);
			lis_matrix_DLU_destroy(Aout);
			lis_matrix_diag_destroy(Aout->WD);
			if( Aout->l2g_map ) lis_free( Aout->l2g_map );
			if( Aout->commtable ) lis_commtable_destroy( Aout->commtable );
			if( Aout->ranges ) lis_free( Aout->ranges );
			lis_matrix_copy_struct(Atmp,Aout);
			lis_free(Atmp);
			return LIS_SUCCESS;
		}
	}
	else
	{
		istmp = LIS_FALSE;
		Atmp = Ain;
	}

	switch( convert_matrix_type )
	{
	case LIS_MATRIX_BSR:
		err = lis_matrix_convert_csr2bsr(Atmp,Aout);
		break;
	case LIS_MATRIX_CSC:
		err = lis_matrix_convert_csr2csc(Atmp,Aout); 
		break;
	case LIS_MATRIX_MSR:
		err = lis_matrix_convert_csr2msr(Atmp,Aout);
		break;
	case LIS_MATRIX_ELL:
		err = lis_matrix_convert_csr2ell(Atmp,Aout);
		break;
	case LIS_MATRIX_DIA:
		err = lis_matrix_convert_csr2dia(Atmp,Aout);
		break;
	case LIS_MATRIX_JAD:
		err = lis_matrix_convert_csr2jad(Atmp,Aout);
		break;
	case LIS_MATRIX_BSC:
		err = lis_matrix_duplicate(Atmp,&Atmp2);
		if( err ) return err;
		err = lis_matrix_convert_csr2csc(Atmp,Atmp2);
		if( err ) return err;
		if( Atmp!=Ain )
		{
			lis_matrix_destroy(Atmp);
		}
		Atmp = Atmp2;
		istmp = LIS_TRUE;
		err = lis_matrix_convert_csc2bsc(Atmp,Aout);
		break;
	case LIS_MATRIX_VBR:
		err = lis_matrix_convert_csr2vbr(Atmp,Aout);
		break;
	case LIS_MATRIX_DNS:
		err = lis_matrix_convert_csr2dns(Atmp,Aout);
		break;
	case LIS_MATRIX_COO:
		err = lis_matrix_convert_csr2coo(Atmp,Aout);
		break;
	default:
		LIS_SETERR_IMP;
		err = LIS_ERR_NOT_IMPLEMENTED;
	}

	if( istmp ) lis_matrix_destroy(Atmp);
	if( err )
	{
		return err;
	}

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_self"
LIS_INT lis_matrix_convert_self(LIS_SOLVER solver)
{
	LIS_INT err;
	LIS_INT storage,block;
	LIS_MATRIX A,B;

	LIS_DEBUG_FUNC_IN;

	A = solver->A;
	storage     = solver->options[LIS_OPTIONS_STORAGE];
	block       = solver->options[LIS_OPTIONS_STORAGE_BLOCK];

	if( storage>0 && A->matrix_type!=storage )
	{
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		err = lis_matrix_set_blocksize(B,block,block,NULL,NULL);
		if( err )
		{
			lis_matrix_destroy(B);
			return err;
		}
		err = lis_matrix_set_type(B,storage);
		if( err )
		{
			lis_matrix_destroy(B);
			return err;
		}
		err = lis_matrix_convert(A,B);
		if( err )
		{
			lis_matrix_destroy(B);
			return err;
		}
		lis_matrix_storage_destroy(A);
		lis_matrix_DLU_destroy(A);
		lis_matrix_diag_destroy(A->WD);
		if( A->l2g_map ) lis_free( A->l2g_map );
		if( A->commtable ) lis_commtable_destroy( A->commtable );
		if( A->ranges ) lis_free( A->ranges );
		err = lis_matrix_copy_struct(B,A);
		if( err ) return err;
		lis_free(B);
		if( A->matrix_type==LIS_MATRIX_JAD )
		{
			A->work = (LIS_SCALAR *)lis_malloc(A->n*sizeof(LIS_SCALAR),"lis_precon_create_bjacobi::A->work");
			if( A->work==NULL )
			{
				LIS_SETERR_MEM(A->n*sizeof(LIS_SCALAR));
				return LIS_OUT_OF_MEMORY;
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy"
LIS_INT lis_matrix_copy(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;
	err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;

	switch( Ain->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_matrix_copy_csr(Ain,Aout);
		break;
	case LIS_MATRIX_CSC:
		err = lis_matrix_copy_csc(Ain,Aout);
		break;
	case LIS_MATRIX_MSR:
		err = lis_matrix_copy_msr(Ain,Aout);
		break;
	case LIS_MATRIX_DIA:
		err = lis_matrix_copy_dia(Ain,Aout);
		break;
	case LIS_MATRIX_ELL:
		err = lis_matrix_copy_ell(Ain,Aout);
		break;
	case LIS_MATRIX_JAD:
		err = lis_matrix_copy_jad(Ain,Aout);
		break;
	case LIS_MATRIX_BSR:
		err = lis_matrix_copy_bsr(Ain,Aout);
		break;
	case LIS_MATRIX_VBR:
		err = lis_matrix_copy_vbr(Ain,Aout);
		break;
	case LIS_MATRIX_DNS:
		err = lis_matrix_copy_dns(Ain,Aout);
		break;
	case LIS_MATRIX_COO:
		err = lis_matrix_copy_coo(Ain,Aout);
		break;
	case LIS_MATRIX_BSC:
		err = lis_matrix_copy_bsc(Ain,Aout);
		break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copyDLU"
LIS_INT lis_matrix_copyDLU(LIS_MATRIX Ain, LIS_MATRIX_DIAG *D, LIS_MATRIX *L, LIS_MATRIX *U)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	switch( Ain->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_matrix_copyDLU_csr(Ain,D,L,U);
		break;
/*
	case LIS_MATRIX_CSC:
		err = lis_matrix_copy_csc(Ain,Aout);
		break;
	case LIS_MATRIX_MSR:
		err = lis_matrix_copy_msr(Ain,Aout);
		break;
	case LIS_MATRIX_DIA:
		err = lis_matrix_copy_dia(Ain,Aout);
		break;
	case LIS_MATRIX_ELL:
		err = lis_matrix_copy_ell(Ain,Aout);
		break;
	case LIS_MATRIX_JAD:
		err = lis_matrix_copy_jad(Ain,Aout);
		break;
	case LIS_MATRIX_BSR:
		err = lis_matrix_copy_bsr(Ain,Aout);
		break;
	case LIS_MATRIX_BSC:
		err = lis_matrix_copy_bsc(Ain,Aout);
		break;
	case LIS_MATRIX_VBR:
		err = lis_matrix_copy_vbr(Ain,Aout);
		break;
	case LIS_MATRIX_DNS:
		err = lis_matrix_copy_dns(Ain,Aout);
		break;
	case LIS_MATRIX_COO:
		err = lis_matrix_copy_coo(Ain,Aout);
		break;
*/
		default:
			LIS_SETERR_IMP;
			*D = NULL;
			*L = NULL;
			*U = NULL;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_axpy"
LIS_INT lis_matrix_axpy(LIS_SCALAR alpha, LIS_MATRIX A, LIS_MATRIX B)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	err = lis_matrix_check(B,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;
	
	switch( A->matrix_type )
	{
	case LIS_MATRIX_DNS:
		err = lis_matrix_axpy_dns(alpha,A,B);
		break;

		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_xpay"
LIS_INT lis_matrix_xpay(LIS_SCALAR alpha, LIS_MATRIX A, LIS_MATRIX B)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	err = lis_matrix_check(B,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;
	
	switch( A->matrix_type )
	{
	case LIS_MATRIX_DNS:
		err = lis_matrix_xpay_dns(alpha,A,B);
		break;

		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_axpyz"
LIS_INT lis_matrix_axpyz(LIS_SCALAR alpha, LIS_MATRIX A, LIS_MATRIX B, LIS_MATRIX C)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	err = lis_matrix_check(B,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;
	
	err = lis_matrix_check(C,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;
	
	switch( A->matrix_type )
	{
	case LIS_MATRIX_DNS:
		err = lis_matrix_axpyz_dns(alpha,A,B,C);
		break;

		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_prepare_symm_diag_scaling"
LIS_INT lis_matrix_prepare_symm_diag_scaling(LIS_VECTOR D)
{
	LIS_INT i,n,invalid;
	LIS_REAL mag;
	LIS_SCALAR *d;
#ifdef USE_MPI
	LIS_INT invalid_global;
#endif

	n       = D->n;
	d       = D->value;
	invalid = 0;

	for(i=0;i<n;i++)
	{
		mag = fabs(d[i]);

		if( mag==0.0 || mag!=mag || mag>LIS_SCALAR_MAX )
		{
			invalid = 1;
			break;
		}
	}

#ifdef USE_MPI
	invalid_global = 0;
	MPI_Allreduce(
		&invalid,
		&invalid_global,
		1,
		LIS_MPI_INT,
		MPI_SUM,
		D->comm);

	invalid = invalid_global;
#endif

	if( invalid )
	{
		LIS_SETERR(
			LIS_BREAKDOWN,
			"symmetric diagonal scaling requires finite nonzero diagonal entries\n");

		return LIS_BREAKDOWN;
	}

#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	for(i=0;i<n;i++)
	{
		d[i] = 1.0 / sqrt(fabs(d[i]));
	}

	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_scale"
LIS_INT lis_matrix_scale(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR D, LIS_INT action)
{
	LIS_INT i,n,np,err;
	LIS_SCALAR *b,*d;

	n  = A->n;
	np = A->np;
	b  = B->value;
	d  = D->value;

	err = lis_matrix_get_diagonal(A,D);
	if( err ) return err;
	if( action==LIS_SCALE_SYMM_DIAG )
	{
		err = lis_matrix_prepare_symm_diag_scaling(D);
		if( err ) return err;
#ifdef USE_MPI
		if( A->np>D->np )
		{
			D->value = (LIS_SCALAR *)lis_realloc(D->value,A->np*sizeof(LIS_SCALAR));
			if( D->value==NULL )
			{
				LIS_SETERR_MEM(A->np*sizeof(LIS_SCALAR));
				return LIS_OUT_OF_MEMORY;
			}
			d = D->value;
		}
		lis_send_recv(A->commtable,d);
#endif

		switch( A->matrix_type )
		{
		case LIS_MATRIX_CSR:
			lis_matrix_scale_symm_csr(A, d);
			break;
		case LIS_MATRIX_CSC:
			lis_matrix_scale_symm_csc(A, d);
			break;
		case LIS_MATRIX_MSR:
			lis_matrix_scale_symm_msr(A, d);
			break;
		case LIS_MATRIX_DIA:
			lis_matrix_scale_symm_dia(A, d);
			break;
		case LIS_MATRIX_ELL:
			lis_matrix_scale_symm_ell(A, d);
			break;
		case LIS_MATRIX_JAD:
			lis_matrix_scale_symm_jad(A, d);
			break;
		case LIS_MATRIX_BSR:
			lis_matrix_scale_symm_bsr(A, d);
			break;
		case LIS_MATRIX_BSC:
			lis_matrix_scale_symm_bsc(A, d);
			break;
		case LIS_MATRIX_DNS:
			lis_matrix_scale_symm_dns(A, d);
			break;
		case LIS_MATRIX_COO:
			lis_matrix_scale_symm_coo(A, d);
			break;
		case LIS_MATRIX_VBR:
			lis_matrix_scale_symm_vbr(A, d);
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
	else if (action==LIS_SCALE_JACOBI)
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0; i<n; i++)
		{
			d[i] = 1.0 / d[i];
		}
		switch( A->matrix_type )
		{
		case LIS_MATRIX_CSR:
			lis_matrix_scale_csr(A, d);
			break;
		case LIS_MATRIX_CSC:
			lis_matrix_scale_csc(A, d);
			break;
		case LIS_MATRIX_MSR:
			lis_matrix_scale_msr(A, d);
			break;
		case LIS_MATRIX_DIA:
			lis_matrix_scale_dia(A, d);
			break;
		case LIS_MATRIX_ELL:
			lis_matrix_scale_ell(A, d);
			break;
		case LIS_MATRIX_JAD:
			lis_matrix_scale_jad(A, d);
			break;
		case LIS_MATRIX_BSR:
			lis_matrix_scale_bsr(A, d);
			break;
		case LIS_MATRIX_BSC:
			lis_matrix_scale_bsc(A, d);
			break;
		case LIS_MATRIX_DNS:
			lis_matrix_scale_dns(A, d);
			break;
		case LIS_MATRIX_COO:
			lis_matrix_scale_coo(A, d);
			break;
		case LIS_MATRIX_VBR:
			lis_matrix_scale_vbr(A, d);
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}

	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
		b[i] = b[i]*d[i];
	}
	A->is_scaled = LIS_TRUE;
	B->is_scaled = LIS_TRUE;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_inf_norm"
LIS_INT lis_matrix_get_inf_norm(LIS_MATRIX A, LIS_REAL *value)
{
	LIS_MATRIX C;
	LIS_PRECON_ROW_SOURCE row;
	LIS_INT err;
	LIS_INT i,k;
	LIS_REAL local_max,global_max,row_sum;

	LIS_DEBUG_FUNC_IN;

	if( A==NULL || value==NULL )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"invalid matrix infinity-norm arguments\n");
		return LIS_ERR_ILL_ARG;
	}

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	/*
	 * CSR and row-backed USER matrices can be traversed directly.
	 * Other explicit storage formats use LIS' existing conversion path.
	 * Matrix operators do not expose entries, so their exact row-sum norm
	 * cannot in general be reconstructed from matvec alone.
	 */
	if( A->matrix_type!=LIS_MATRIX_CSR && A->matrix_type!=LIS_MATRIX_USER )
	{
		if( A->matrix_type==LIS_MATRIX_OPERATOR )
		{
			LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
				"exact infinity norm requires explicit or row-backed matrix entries\n");
			return LIS_ERR_NOT_IMPLEMENTED;
		}

		C = NULL;
		err = lis_matrix_duplicate(A,&C);
		if( err ) return err;

		err = lis_matrix_set_type(C,LIS_MATRIX_CSR);
		if( !err ) err = lis_matrix_convert(A,C);
		if( !err ) err = lis_matrix_get_inf_norm(C,value);

		lis_matrix_destroy(C);
		return err;
	}

	err = lis_precon_row_source_init(A,&row);
	if( err ) return err;

	local_max = (LIS_REAL)0.0;
	for(i=0;i<A->n;i++)
	{
		err = lis_precon_row_source_get(&row,i);
		if( err )
		{
			lis_precon_row_source_destroy(&row);
			return err;
		}

		row_sum = (LIS_REAL)0.0;
		for(k=0;k<row.nnz;k++)
		{
			row_sum += (LIS_REAL)fabs(row.value[k]);
		}

		if( row_sum>local_max ) local_max = row_sum;
	}

	lis_precon_row_source_destroy(&row);

#ifdef USE_MPI
	err = MPI_Allreduce(
		&local_max,&global_max,1,LIS_MPI_REAL,MPI_MAX,A->comm);
	if( err!=MPI_SUCCESS )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"MPI_Allreduce failed for matrix infinity norm\n");
		return LIS_ERR_ILL_ARG;
	}
#else
	global_max = local_max;
#endif

	*value = global_max;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_matrix_psd_reset_scale"
LIS_INT lis_matrix_psd_reset_scale(LIS_MATRIX A)
{
	A->is_scaled=LIS_FALSE;
	return LIS_SUCCESS;
}


static LIS_INT
lis_matrix_get_diagonal_user_rows(
	LIS_MATRIX A,
	LIS_SCALAR *d)
{
	LIS_PRECON_ROW_SOURCE row;
	LIS_INT err;
	LIS_INT i,k;
	LIS_INT global_row;
	LIS_SCALAR diagonal;

	err = lis_precon_row_source_init(A,&row);
	if( err ) return err;

	for(i=0;i<A->n;i++)
	{
		err = lis_precon_row_source_get(&row,i);
		if( err )
		{
			lis_precon_row_source_destroy(&row);
			return err;
		}

		global_row = A->is+i;
		diagonal = (LIS_SCALAR)0.0;

		for(k=0;k<row.nnz;k++)
		{
			if( row.index[k]==global_row )
			{
				diagonal += row.value[k];
			}
		}

		d[i] = diagonal;
	}

	lis_precon_row_source_destroy(&row);

	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal"
LIS_INT lis_matrix_get_diagonal(LIS_MATRIX A, LIS_VECTOR D)
{
	LIS_INT err;
	LIS_SCALAR *d;

	LIS_DEBUG_FUNC_IN;

	d = D->value;
	switch( A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		lis_matrix_get_diagonal_csr(A, d);
		break;
	case LIS_MATRIX_CSC:
		lis_matrix_get_diagonal_csc(A, d);
		break;
	case LIS_MATRIX_MSR:
		lis_matrix_get_diagonal_msr(A, d);
		break;
	case LIS_MATRIX_DIA:
		lis_matrix_get_diagonal_dia(A, d);
		break;
	case LIS_MATRIX_ELL:
		lis_matrix_get_diagonal_ell(A, d);
		break;
	case LIS_MATRIX_JAD:
		lis_matrix_get_diagonal_jad(A, d);
		break;
	case LIS_MATRIX_BSR:
		lis_matrix_get_diagonal_bsr(A, d);
		break;
	case LIS_MATRIX_BSC:
		lis_matrix_get_diagonal_bsc(A, d);
		break;
	case LIS_MATRIX_DNS:
		lis_matrix_get_diagonal_dns(A, d);
		break;
	case LIS_MATRIX_COO:
		lis_matrix_get_diagonal_coo(A, d);
		break;
	case LIS_MATRIX_VBR:
		lis_matrix_get_diagonal_vbr(A, d);
		break;
	case LIS_MATRIX_USER:
	        if( A->user_get_diagonal!=NULL )
	        {
	                err = A->user_get_diagonal(A->user_data,d);
	        }
	        else if( A->user_get_row!=NULL )
	        {
	                err = lis_matrix_get_diagonal_user_rows(A,d);
	        }
	        else
	        {
	                LIS_SETERR(
	                        LIS_ERR_NOT_IMPLEMENTED,
	                        "user diagonal callback and row callback are not set\n");
	                return LIS_ERR_NOT_IMPLEMENTED;
	        }
	        if( err ) return err;
	        break;

	case LIS_MATRIX_OPERATOR:
	        if( A->operator_scale!=NULL )
	        {
	                if( A->operator_A==NULL )
	                {
	                        LIS_SETERR(
	                                LIS_ERR_ILL_ARG,
	                                "scaled matrix operator is not initialized\n");
	                        return LIS_ERR_ILL_ARG;
	                }

	                err = lis_matrix_get_diagonal(
	                        A->operator_A,D);
	                if( err ) return err;

	                err = lis_vector_pmul(
	                        D,A->operator_scale,D);
	                if( err ) return err;

	                err = lis_vector_pmul(
	                        D,A->operator_scale,D);
	                if( err ) return err;

	                break;
	        }

	        if( A->operator_A==NULL || A->operator_B==NULL ||
	            A->operator_work==NULL )
	        {
	                LIS_SETERR(LIS_ERR_ILL_ARG,
	                           "matrix operator is not initialized\n");
	                return LIS_ERR_ILL_ARG;
	        }

	        err = lis_matrix_get_diagonal(
	                A->operator_A,A->operator_work);
	        if( err ) return err;

	        err = lis_matrix_get_diagonal(A->operator_B,D);
	        if( err ) return err;

	        err = lis_vector_scale(
	                A->operator_alpha,A->operator_work);
	        if( err ) return err;

	        err = lis_vector_scale(A->operator_beta,D);
	        if( err ) return err;

	        err = lis_vector_axpy(
	                (LIS_SCALAR)1.0,A->operator_work,D);
	        if( err ) return err;
	        break;

	default:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
		break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal"
LIS_INT lis_matrix_shift_diagonal(LIS_MATRIX A, LIS_SCALAR sigma)
{

	LIS_DEBUG_FUNC_IN;

	switch( A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		lis_matrix_shift_diagonal_csr(A, sigma);
		break;
	case LIS_MATRIX_CSC:
		lis_matrix_shift_diagonal_csc(A, sigma);
		break;
	case LIS_MATRIX_MSR:
		lis_matrix_shift_diagonal_msr(A, sigma);
		break;
	case LIS_MATRIX_DIA:
		lis_matrix_shift_diagonal_dia(A, sigma);
		break;
	case LIS_MATRIX_ELL:
		lis_matrix_shift_diagonal_ell(A, sigma);
		break;
	case LIS_MATRIX_JAD:
		lis_matrix_shift_diagonal_jad(A, sigma);
		break;
	case LIS_MATRIX_BSR:
		lis_matrix_shift_diagonal_bsr(A, sigma);
		break;
	case LIS_MATRIX_BSC:
		lis_matrix_shift_diagonal_bsc(A, sigma);
		break;
	case LIS_MATRIX_DNS:
		lis_matrix_shift_diagonal_dns(A, sigma);
		break;
	case LIS_MATRIX_COO:
		lis_matrix_shift_diagonal_coo(A, sigma);
		break;
	case LIS_MATRIX_VBR:
		lis_matrix_shift_diagonal_vbr(A, sigma);
		break;

	default:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
		break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_matrix"
LIS_INT lis_matrix_shift_matrix(LIS_MATRIX A, LIS_MATRIX B, LIS_SCALAR sigma)
{
	LIS_INT err;
	LIS_MATRIX Atmp,Btmp;
  
	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_duplicate(A,&Atmp);
	err = lis_matrix_duplicate(B,&Btmp);	
	if( err ) return err;
	lis_matrix_set_type(Atmp,LIS_MATRIX_DNS);
	err = lis_matrix_convert(A,Atmp);
	if( err ) return err;
	lis_matrix_set_type(Btmp,LIS_MATRIX_DNS);
	err = lis_matrix_convert(B,Btmp);
	if( err ) return err;
	lis_matrix_axpy_dns(-sigma,Btmp,Atmp);
	lis_matrix_convert(Atmp,A);
	lis_matrix_destroy(Atmp);	
	lis_matrix_destroy(Btmp);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split"
LIS_INT lis_matrix_split(LIS_MATRIX A)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	if( A->is_splited )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	switch( A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_matrix_split_csr(A);
		break;
	case LIS_MATRIX_CSC:
		err = lis_matrix_split_csc(A);
		break;
	case LIS_MATRIX_BSR:
		err = lis_matrix_split_bsr(A);
		break;
	case LIS_MATRIX_MSR:
		err = lis_matrix_split_msr(A);
		break;
	case LIS_MATRIX_ELL:
		err = lis_matrix_split_ell(A);
		break;
	case LIS_MATRIX_DIA:
		err = lis_matrix_split_dia(A);
		break;
	case LIS_MATRIX_JAD:
		err = lis_matrix_split_jad(A);
		break;
	case LIS_MATRIX_BSC:
		err = lis_matrix_split_bsc(A);
		break;
	case LIS_MATRIX_DNS:
		err = lis_matrix_split_dns(A);
		break;
	case LIS_MATRIX_COO:
		err = lis_matrix_split_coo(A);
		break;
	case LIS_MATRIX_VBR:
		err = lis_matrix_split_vbr(A);
		break;
	default:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
		break;
	}

	if( err ) return err;
	/*
       	if( !A->is_save )
       	{
       		lis_matrix_storage_destroy(A);
       	}
	*/
	A->is_splited = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_matrix_split_create"
LIS_INT lis_matrix_split_create(LIS_MATRIX A)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

    /* don't try to split if it is already split */
	if( A->is_splited )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}

	switch( A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_matrix_split_create_csr(A);
		break;
    case LIS_MATRIX_CSC:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_BSR:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_MSR:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_ELL:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_DIA:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_JAD:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_BSC:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_DNS:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_COO:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_VBR:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
	default:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
	}

	if( err ) return err;

	A->is_splited = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_matrix_split_update"
LIS_INT lis_matrix_split_update(LIS_MATRIX A)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

    /* don't try to perform "split_update" if the matrix is not already split */
    if( !A->is_splited )
    {
        return LIS_ERR_ILL_ARG;
    }

	switch( A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_matrix_split_update_csr(A);
		break;
    case LIS_MATRIX_CSC:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_BSR:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_MSR:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_ELL:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_DIA:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_JAD:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_BSC:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_DNS:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_COO:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
    case LIS_MATRIX_VBR:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
	default:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
	}

	if( err ) return err;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_merge"
LIS_INT lis_matrix_merge(LIS_MATRIX A)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	if( !A->is_splited || (A->is_save && A->is_splited) )
	{
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	switch( A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_matrix_merge_csr(A);
		break;
	case LIS_MATRIX_CSC:
		err = lis_matrix_merge_csc(A);
		break;
	case LIS_MATRIX_MSR:
		err = lis_matrix_merge_msr(A);
		break;
	case LIS_MATRIX_BSR:
		err = lis_matrix_merge_bsr(A);
		break;
	case LIS_MATRIX_ELL:
		err = lis_matrix_merge_ell(A);
		break;
	case LIS_MATRIX_JAD:
		err = lis_matrix_merge_jad(A);
		break;
	case LIS_MATRIX_DIA:
		err = lis_matrix_merge_dia(A);
		break;
	case LIS_MATRIX_BSC:
		err = lis_matrix_merge_bsc(A);
		break;
	case LIS_MATRIX_DNS:
		err = lis_matrix_merge_dns(A);
		break;
	case LIS_MATRIX_COO:
		err = lis_matrix_merge_coo(A);
		break;
	case LIS_MATRIX_VBR:
		err = lis_matrix_merge_vbr(A);
		break;
	default:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
		break;
	}

	if( err ) return err;
	if( !A->is_save )
	{
		lis_matrix_DLU_destroy(A);
		A->is_splited = LIS_FALSE;
	}


	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solve"
LIS_INT lis_matrix_solve(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT flag)
{
	LIS_DEBUG_FUNC_IN;

	if( !A->is_splited ) lis_matrix_split(A);

	switch( A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		lis_matrix_solve_csr(A,b,x,flag);
		break;
	case LIS_MATRIX_BSR:
		lis_matrix_solve_bsr(A,b,x,flag);
		break;
	case LIS_MATRIX_CSC:
		lis_matrix_solve_csc(A,b,x,flag);
		break;
	case LIS_MATRIX_MSR:
		lis_matrix_solve_msr(A,b,x,flag);
		break;
	case LIS_MATRIX_ELL:
		lis_matrix_solve_ell(A,b,x,flag);
		break;
	case LIS_MATRIX_JAD:
		lis_matrix_solve_jad(A,b,x,flag);
		break;
	case LIS_MATRIX_DIA:
		lis_matrix_solve_dia(A,b,x,flag);
		break;
	case LIS_MATRIX_DNS:
		lis_matrix_solve_dns(A,b,x,flag);
		break;
	case LIS_MATRIX_BSC:
		lis_matrix_solve_bsc(A,b,x,flag);
		break;
	case LIS_MATRIX_VBR:
		lis_matrix_solve_vbr(A,b,x,flag);
		break;
	default:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solveh"
LIS_INT lis_matrix_solveh(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT flag)
{
	LIS_DEBUG_FUNC_IN;

	if( !A->is_splited ) lis_matrix_split(A);

	switch( A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		lis_matrix_solveh_csr(A,b,x,flag);
		break;
	case LIS_MATRIX_BSR:
		lis_matrix_solveh_bsr(A,b,x,flag);
		break;
	case LIS_MATRIX_CSC:
		lis_matrix_solveh_csc(A,b,x,flag);
		break;
	case LIS_MATRIX_MSR:
		lis_matrix_solveh_msr(A,b,x,flag);
		break;
	case LIS_MATRIX_ELL:
		lis_matrix_solveh_ell(A,b,x,flag);
		break;
	case LIS_MATRIX_JAD:
		lis_matrix_solveh_jad(A,b,x,flag);
		break;
	case LIS_MATRIX_DIA:
		lis_matrix_solveh_dia(A,b,x,flag);
		break;
	case LIS_MATRIX_DNS:
		lis_matrix_solveh_dns(A,b,x,flag);
		break;
	case LIS_MATRIX_BSC:
		lis_matrix_solveh_bsc(A,b,x,flag);
		break;
	case LIS_MATRIX_VBR:
		lis_matrix_solveh_vbr(A,b,x,flag);
		break;
	default:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}





