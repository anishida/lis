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
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_matrix_init
 * lis_matrix_check
 * lis_matrix_create
 * lis_matrix_destroy
 * lis_matrix_duplicate
 * lis_matrix_malloc
 * lis_matrix_assemble
 * lis_matrix_is_assembled
 * lis_matrix_get_range
 * lis_matrix_get_size
 * lis_matrix_get_nnz
 * lis_matrix_get_type
 * lis_matrix_get_option
 * lis_matrix_set_type
 * lis_matrix_set_value
 * lis_matrix_set_option
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_matrix_init"
LIS_INT lis_matrix_init(LIS_MATRIX *Amat)
{

	LIS_DEBUG_FUNC_IN;

	memset(*Amat,0,sizeof(struct LIS_MATRIX_STRUCT));

	(*Amat)->label       = LIS_LABEL_MATRIX;
	(*Amat)->matrix_type = LIS_MATRIX_CSR;
	(*Amat)->status      = LIS_MATRIX_DECIDING_SIZE;
/*	(*Amat)->status      = LIS_MATRIX_NULL;*/
	(*Amat)->w_annz      = LIS_MATRIX_W_ANNZ;
	(*Amat)->conv_bnr    = 2;
	(*Amat)->conv_bnc    = 2;
	(*Amat)->is_destroy  = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_check"
LIS_INT lis_matrix_check(LIS_MATRIX A, LIS_INT level)
{
	LIS_DEBUG_FUNC_IN;

	switch( level )
	{
	case LIS_MATRIX_CHECK_NULL:
		if( !lis_is_malloc(A) )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix A is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	case LIS_MATRIX_CHECK_SIZE:
		if( !lis_is_malloc(A) )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix A is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		if( A->status==LIS_MATRIX_DECIDING_SIZE )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix size is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	case LIS_MATRIX_CHECK_TYPE:
		if( !lis_is_malloc(A) )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix A is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		if( A->status==LIS_MATRIX_DECIDING_SIZE )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix size is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
        /* only want this condition satisfied for the case of non-zero number of DOFs */
        /* on the current process */
		if( A->status==LIS_MATRIX_NULL && A->n>0 )
		{
            /* as the status is actually changed from ...NULL to ...ASSEMBLING in lis_matrix_set_value, */
            /* ...NULL here really means that the components are undefined, not the type */
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix components undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	case LIS_MATRIX_CHECK_NOT_ASSEMBLED:
		if( !lis_is_malloc(A) )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix A is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		if( A->status!=LIS_MATRIX_DECIDING_SIZE && A->status!=LIS_MATRIX_NULL && A->status!=LIS_MATRIX_ASSEMBLING )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix A has already been assembled\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	case LIS_MATRIX_CHECK_SET:
		if( !lis_is_malloc(A) )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix A is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		if( A->status==LIS_MATRIX_DECIDING_SIZE )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix size is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		if( A->status!=LIS_MATRIX_NULL )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix A has already been assembled\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	default:
		if( !lis_is_malloc(A) )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix A is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		if( A->status==LIS_MATRIX_DECIDING_SIZE )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix size is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
/*        if( A->status==LIS_MATRIX_NULL )*/
		if( A->status==LIS_MATRIX_NULL && A->n>0 )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix type is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
/*        if( A->status<=LIS_MATRIX_ASSEMBLING )*/
		if( A->status<=LIS_MATRIX_ASSEMBLING && A->n>0 )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"matrix A is assembling\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_create"
LIS_INT lis_matrix_create(LIS_Comm comm, LIS_MATRIX *Amat)
{
	LIS_INT nprocs,my_rank;  
	int int_nprocs,int_my_rank;

	LIS_DEBUG_FUNC_IN;

	*Amat = NULL;


	*Amat = (LIS_MATRIX)lis_malloc( sizeof(struct LIS_MATRIX_STRUCT),"lis_matrix_create::Amat" );
	if( NULL==*Amat )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_matrix_init(Amat);

	#ifdef USE_MPI
#if _DEBUG
	printf("c_comm = %d f_comm = %d\n",MPI_COMM_WORLD,comm);
#endif
		MPI_Comm_size(comm,&int_nprocs);
		MPI_Comm_rank(comm,&int_my_rank);
		nprocs = int_nprocs;
		my_rank = int_my_rank;
	#else
		nprocs  = 1;
		my_rank = 0;
	#endif
	(*Amat)->comm        = comm;
	(*Amat)->nprocs      = nprocs;
	(*Amat)->my_rank     = my_rank;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_size"
LIS_INT lis_matrix_set_size(LIS_MATRIX Amat, LIS_INT local_n, LIS_INT global_n)
{
	LIS_INT nprocs,my_rank;
	int int_nprocs,int_my_rank;
	LIS_INT is,ie;
	LIS_INT err;
	LIS_INT *ranges;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(Amat,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;

	if( global_n>0 && local_n>global_n )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%D) is larger than global n(=%D)\n",local_n,global_n);
		return LIS_ERR_ILL_ARG;
	}
	if( local_n<0 || global_n<0 )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%D) or global n(=%D) are less than 0\n",local_n,global_n);
		return LIS_ERR_ILL_ARG;
	}
    /* actually want to allow this case for MPI, but not otherwise */
    #ifndef USE_MPI
    if( local_n==0 && global_n==0 )
    {
        LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%D) and global n(=%D) are 0\n",local_n,global_n);
        return LIS_ERR_ILL_ARG;
    }
    #endif
	#ifdef USE_MPI
	MPI_Comm_size(Amat->comm,&int_nprocs);
	nprocs = int_nprocs;
    /* change the logic to something a little more direct, i.e., that will not */
    /* throw an error if local_n=0 . . . */
	if( global_n>0 && global_n<nprocs )
	  {
	    LIS_SETERR2(LIS_ERR_ILL_ARG,"global n(=%D) is smaller than nprocs(=%D)\n",global_n,nprocs);
	    return LIS_ERR_ILL_ARG;
	  }
	#endif

	err = lis_ranges_create(Amat->comm,&local_n,&global_n,&ranges,&is,&ie,&nprocs,&my_rank);
	if( err )
	{
		return err;
	}

	Amat->status      = LIS_MATRIX_NULL;
	Amat->ranges      = ranges;
	Amat->n           = local_n;
	Amat->gn          = global_n;
	Amat->np          = local_n;
	Amat->my_rank     = my_rank;
	Amat->nprocs      = nprocs;
	Amat->is          = is;
	Amat->ie          = ie;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_LU_create"
LIS_INT lis_matrix_LU_create(LIS_MATRIX A)
{
	LIS_MATRIX_CORE		L,U;

	LIS_DEBUG_FUNC_IN;

	L = NULL;
	U = NULL;

	L = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_LU_create::L");
	if( L==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	U = (LIS_MATRIX_CORE)lis_calloc(sizeof(struct LIS_MATRIX_CORE_STRUCT),"lis_matrix_LU_create::U");
	if( U==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_CORE_STRUCT));
		lis_free(L);
		return LIS_OUT_OF_MEMORY;
	}

	A->L = L;
	A->U = U;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_LU_destroy"
LIS_INT lis_matrix_LU_destroy(LIS_MATRIX_CORE Amat)
{
	LIS_DEBUG_FUNC_IN;

	if( Amat )
	{
		if( Amat->ptr ) lis_free( Amat->ptr );
		if( Amat->row ) lis_free( Amat->row );
		if( Amat->col ) lis_free( Amat->col );
		if( Amat->index) lis_free( Amat->index );
		if( Amat->bptr ) lis_free( Amat->bptr );
		if( Amat->bindex ) lis_free( Amat->bindex );
		if( Amat->value ) lis_free( Amat->value );
		if( Amat->work ) lis_free( Amat->work );
		lis_free(Amat);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_DLU_destroy"
LIS_INT lis_matrix_DLU_destroy(LIS_MATRIX Amat)
{
	LIS_DEBUG_FUNC_IN;

	if( Amat )
	{
		if( Amat->D ) lis_matrix_diag_destroy(Amat->D);
		if( Amat->L ) lis_matrix_LU_destroy(Amat->L);
		if( Amat->U ) lis_matrix_LU_destroy(Amat->U);
		Amat->D = NULL;
		Amat->L = NULL;
		Amat->U = NULL;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_storage_destroy"
LIS_INT lis_matrix_storage_destroy(LIS_MATRIX Amat)
{

	LIS_DEBUG_FUNC_IN;

	if( Amat )
	{
		if( Amat->is_destroy )
		{
			if( Amat->ptr ) lis_free( Amat->ptr );
			if( Amat->row ) lis_free( Amat->row );
			if( Amat->col ) lis_free( Amat->col );
			if( Amat->index) lis_free( Amat->index );
			if( Amat->bptr ) lis_free( Amat->bptr );
			if( Amat->bindex ) lis_free( Amat->bindex );
			if( Amat->value ) lis_free( Amat->value );
			if( Amat->work ) lis_free( Amat->work );
			if( Amat->conv_row ) lis_free( Amat->conv_row );
			if( Amat->conv_col ) lis_free( Amat->conv_col );
			if( Amat->w_nnz ) lis_free( Amat->w_nnz );
			if( Amat->w_row )
			{
#if 0
				for(i=0;i<Amat->n;i++)
				{
					lis_free(Amat->w_index[i]);
					lis_free(Amat->w_value[i]);
				}
#else
				lis_free_mat(Amat);
#endif
				lis_free2(3,Amat->w_row,Amat->w_index,Amat->w_value);
			}
		}
		Amat->row       = NULL;
		Amat->col       = NULL;
		Amat->ptr       = NULL;
		Amat->index     = NULL;
		Amat->bptr      = NULL;
		Amat->bindex    = NULL;
		Amat->value     = NULL;
		Amat->work      = NULL;
		Amat->conv_row  = NULL;
		Amat->conv_col  = NULL;
		Amat->w_nnz     = NULL;
		Amat->w_row     = NULL;
		Amat->w_index   = NULL;
		Amat->w_value   = NULL;
		Amat->v_value   = NULL;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_destroy"
LIS_INT lis_matrix_destroy(LIS_MATRIX Amat)
{
	LIS_DEBUG_FUNC_IN;

	if( lis_is_malloc(Amat) )
	{
		if( !Amat->is_fallocated ) lis_matrix_storage_destroy(Amat);
		lis_matrix_DLU_destroy(Amat);
		lis_matrix_diag_destroy(Amat->WD);
		if( Amat->l2g_map ) lis_free( Amat->l2g_map );
		if( Amat->commtable ) lis_commtable_destroy( Amat->commtable );
		if( Amat->ranges ) lis_free( Amat->ranges );
		lis_free(Amat);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_duplicate"
LIS_INT lis_matrix_duplicate(LIS_MATRIX Ain, LIS_MATRIX *Aout)
{
	LIS_INT err;
	LIS_INT nprocs;
	#ifdef USE_MPI
		LIS_INT i;
		LIS_INT *ranges;
		LIS_INT *l2g_map;
	#endif

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	nprocs = Ain->nprocs;
	*Aout = NULL;

	*Aout = (LIS_MATRIX)lis_malloc( sizeof(struct LIS_MATRIX_STRUCT),"lis_matrix_duplicate::Aout" );
	if( NULL==*Aout )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_matrix_init(Aout);

	#ifdef USE_MPI
		l2g_map = (LIS_INT *)lis_malloc( (Ain->np-Ain->n)*sizeof(LIS_INT),"lis_matrix_duplicate::l2g_map" );
		if( l2g_map==NULL )
		{
			LIS_SETERR_MEM((Ain->np-Ain->n)*sizeof(LIS_INT));
			lis_matrix_destroy(*Aout);
			*Aout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		memcpy(l2g_map,Ain->l2g_map,(Ain->np-Ain->n)*sizeof(LIS_INT));
		ranges = (LIS_INT *)lis_malloc( (nprocs+1)*sizeof(LIS_INT),"lis_matrix_duplicate::range" );
		if( ranges==NULL )
		{
			LIS_SETERR_MEM((nprocs+1)*sizeof(LIS_INT));
			lis_free(l2g_map);
			lis_matrix_destroy(*Aout);
			*Aout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		for(i=0;i<nprocs+1;i++) ranges[i] = Ain->ranges[i];
		(*Aout)->ranges      = ranges;
		(*Aout)->l2g_map     = l2g_map;
	#else
		(*Aout)->ranges      = NULL;
	#endif

	(*Aout)->status      = LIS_MATRIX_NULL;
	(*Aout)->is_block    = Ain->is_block;
	(*Aout)->n           = Ain->n;
	(*Aout)->gn          = Ain->gn;
	(*Aout)->np          = Ain->np;
	(*Aout)->comm        = Ain->comm;
	(*Aout)->my_rank     = Ain->my_rank;
	(*Aout)->nprocs      = Ain->nprocs;
	(*Aout)->is          = Ain->is;
	(*Aout)->ie          = Ain->ie;
	(*Aout)->origin      = Ain->origin;
	(*Aout)->is_destroy  = Ain->is_destroy;

#ifdef USE_MPI
	err = lis_commtable_duplicateM(Ain,Aout);
	if( err )
	{
		lis_matrix_destroy(*Aout);
		*Aout = NULL;
		return err;
	}
	(*Aout)->is_comm   = LIS_TRUE;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_struct"
LIS_INT lis_matrix_copy_struct(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_DEBUG_FUNC_IN;

	memcpy(Aout,Ain,sizeof(struct LIS_MATRIX_STRUCT));

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_range"
LIS_INT lis_matrix_get_range(LIS_MATRIX A, LIS_INT *is, LIS_INT *ie)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_SIZE);
	if( err ) return err;

	*is = A->is;
	*ie = A->ie;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_size"
LIS_INT lis_matrix_get_size(LIS_MATRIX A, LIS_INT *local_n, LIS_INT *global_n)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_SIZE);
	if( err ) return err;

	*local_n  = A->n;
	*global_n = A->gn;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_nnz"
LIS_INT lis_matrix_get_nnz(LIS_MATRIX A, LIS_INT *nnz)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_SIZE);
	if( err ) return err;

	*nnz = A->nnz;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_assemble"
LIS_INT lis_matrix_assemble(LIS_MATRIX A)
{
	LIS_INT err;
	LIS_INT matrix_type;
	LIS_MATRIX B;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_TYPE);
	if( err ) return err;

	matrix_type = A->matrix_type;
	if( A->status==LIS_MATRIX_ASSEMBLING )
	{
		A->matrix_type = LIS_MATRIX_RCO;
		A->status      = LIS_MATRIX_RCO;
		#ifdef USE_MPI
			err = lis_matrix_g2l(A);
			if( err ) return err;
			err = lis_commtable_create(A);
			if( err ) return err;
			A->is_comm = LIS_TRUE;
		#endif
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,matrix_type);
		err = lis_matrix_convert(A,B);
		if( err ) return err;
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
			A->work = (LIS_SCALAR *)lis_malloc(A->n*sizeof(LIS_SCALAR),"lis_matrix_assemble::A->work");
			if( A->work==NULL )
			{
				LIS_SETERR_MEM(A->n*sizeof(LIS_SCALAR));
				return LIS_OUT_OF_MEMORY;
			}
		}

		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	else if( A->status<0 && A->n>0 )
	{
		A->status      = -A->status;
		A->matrix_type = A->status;
		if( A->matrix_type==LIS_MATRIX_JAD )
		{
			A->work = (LIS_SCALAR *)lis_malloc(A->n*sizeof(LIS_SCALAR),"lis_matrix_assemble::A->work");
			if( A->work==NULL )
			{
				LIS_SETERR_MEM(A->n*sizeof(LIS_SCALAR));
				return LIS_OUT_OF_MEMORY;
			}
		}
	}
	#ifdef USE_MPI
		if( !A->is_pmat )
		{
			if( A->l2g_map==NULL )
			{
				err = lis_matrix_g2l(A);
				if( err ) return err;
			}
			if( A->commtable==NULL )			
			{
				err = lis_commtable_create(A);
				if( err ) return err;
			}
		}
	#endif
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_is_assembled"
LIS_INT lis_matrix_is_assembled(LIS_MATRIX A)
{
  if( A->status!=LIS_MATRIX_NULL ) return !LIS_SUCCESS;
  else                             return  LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_value"
LIS_INT lis_matrix_set_value(LIS_INT flag, LIS_INT i, LIS_INT j, LIS_SCALAR value, LIS_MATRIX A)
{
	LIS_INT n,gn,is,k;
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_NOT_ASSEMBLED);
	if( err ) return err;

	n   = A->n;
	gn  = A->gn;
	is  = A->is;
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

	if( A->status==LIS_MATRIX_NULL )
	{
		if( A->w_nnz==NULL )
		{
			A->w_nnz = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_set_value::A->w_nnz");
			if( A->w_nnz==NULL )
			{
				LIS_SETERR_MEM(n*sizeof(LIS_INT));
				return LIS_OUT_OF_MEMORY;
			}
			for(k=0;k<n;k++) A->w_nnz[k] = A->w_annz;
		}
		err = lis_matrix_malloc_rco(n,A->w_nnz,&A->w_row,&A->w_index,&A->w_value);
		if( err )
		{
			lis_free(A->w_nnz);
			return err;
		}
		A->status  = LIS_MATRIX_ASSEMBLING;
		A->is_copy = LIS_TRUE;
	}
	if( A->w_nnz[i-is]==A->w_row[i-is] )
	{
		A->w_nnz[i-is] += A->w_annz;
		err = lis_matrix_realloc_rco(i-is,A->w_nnz[i-is],&A->w_index,&A->w_value);
		if( err )
		{
			for(k=0;k<n;k++)
			{
				lis_free(A->w_index[k]);
				lis_free(A->w_value[k]);
			}
			lis_free2(4,A->w_nnz,A->w_row,A->w_index,A->w_value);
			return err;
		}
	}
	for(k=0;k<A->w_row[i-is];k++)
	{
		if( A->w_index[i-is][k]==j ) break;
	}
	if( k<A->w_row[i-is] )
	{
		if( flag==LIS_INS_VALUE )
		{
			A->w_value[i-is][k] = value;
		}
		else
		{
			A->w_value[i-is][k] += value;
		}
	}
	else
	{
		k                   = A->w_row[i-is]++;
		A->w_index[i-is][k] = j;
		A->w_value[i-is][k] = value;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_matrix_psd_set_value"
LIS_INT lis_matrix_psd_set_value(LIS_INT flag, LIS_INT i, LIS_INT j, LIS_SCALAR value, LIS_MATRIX A)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

    /* presuming that the matrix is assembled */
    /* still, if not, default defined below */
    switch(A->status) {
        case LIS_MATRIX_CSR:
            err=lis_matrix_psd_set_value_csr(flag,i,j,value,A);
            if (err) return err;
            break;
        case LIS_MATRIX_CSC:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
        case LIS_MATRIX_MSR:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
        case LIS_MATRIX_DIA:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
        case LIS_MATRIX_ELL:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
        case LIS_MATRIX_JAD:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
        case LIS_MATRIX_BSR:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
        case LIS_MATRIX_BSC:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
        case LIS_MATRIX_VBR:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
        case LIS_MATRIX_COO:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
        case LIS_MATRIX_DNS:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
        default:
            err = LIS_ERR_NOT_IMPLEMENTED;
            return err;
    }

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_values"
LIS_INT lis_matrix_set_values(LIS_INT flag, LIS_INT n, LIS_SCALAR value[], LIS_MATRIX A)
{
  LIS_INT i,j;
  LIS_DEBUG_FUNC_IN;

  for (i=0;i<n;i++)
    {
      for (j=0;j<n;j++)
	{
	  lis_matrix_set_value(flag, i, j, value[i*n+j], A);
	}
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_type"
LIS_INT lis_matrix_set_type(LIS_MATRIX A, LIS_INT matrix_type)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_NOT_ASSEMBLED);
	if( err ) return err;
	if( matrix_type < LIS_MATRIX_CSR || matrix_type > LIS_MATRIX_DNS )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"matrix_type is %D (Set between 1 to %D)\n", matrix_type, LIS_MATRIX_DNS);
		return LIS_ERR_ILL_ARG;
	}

	A->matrix_type = matrix_type;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_type"
LIS_INT lis_matrix_get_type(LIS_MATRIX A, LIS_INT *matrix_type)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;

	*matrix_type = A->matrix_type;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_destroyflag"
LIS_INT lis_matrix_set_destroyflag(LIS_MATRIX A, LIS_INT flag)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;

	if( flag )
	{
		A->is_destroy = LIS_TRUE;
	}
	else
	{
		A->is_destroy = LIS_FALSE;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_destroyflag"
LIS_INT lis_matrix_get_destroyflag(LIS_MATRIX A, LIS_INT *flag)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;

	*flag = A->is_destroy;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc"
LIS_INT lis_matrix_malloc(LIS_MATRIX A, LIS_INT nnz_row, LIS_INT nnz[])
{
	LIS_INT n,k;
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_NOT_ASSEMBLED);
	if( err ) return err;

	n   = A->n;
	if( A->w_nnz==NULL )
	{
		A->w_nnz = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_malloc::A->w_nnz");
		if( A->w_nnz==NULL )
		{
			LIS_SETERR_MEM(n*sizeof(LIS_INT));
			return LIS_OUT_OF_MEMORY;
		}
	}
	if( nnz==NULL )
	{
		A->w_annz = nnz_row;
		for(k=0;k<n;k++) A->w_nnz[k] = nnz_row;
	}
	else
	{
		for(k=0;k<n;k++) A->w_nnz[k] = nnz[k];
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_unset"
LIS_INT lis_matrix_unset(LIS_MATRIX A)
{
	LIS_INT err;

 	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_SIZE);
	if( err ) return err;
	
	if( A->is_copy )
	{
		lis_matrix_storage_destroy(A);
	}
	A->row     = NULL;
	A->col     = NULL;
	A->ptr     = NULL;
	A->index   = NULL;
	A->bptr    = NULL;
	A->bindex  = NULL;
	A->value   = NULL;
	A->is_copy = LIS_FALSE;
	A->status  = LIS_MATRIX_NULL;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

