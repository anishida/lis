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


#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_init"
LIS_INT lis_matrix_diag_init(LIS_MATRIX_DIAG *D)
{
	LIS_DEBUG_FUNC_IN;

	memset(*D,0,sizeof(struct LIS_MATRIX_DIAG_STRUCT));
	(*D)->bn = 1;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_check"
LIS_INT lis_matrix_diag_check(LIS_MATRIX_DIAG D, LIS_INT level)
{
	LIS_DEBUG_FUNC_IN;

	switch( level )
	{
	case LIS_MATRIX_CHECK_NULL:
		if( D==NULL )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"diagonal matrix D is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	default:
		if( D==NULL )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"diagonal matrix D is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_create"
LIS_INT lis_matrix_diag_create(LIS_INT local_n, LIS_INT global_n, LIS_Comm comm, LIS_MATRIX_DIAG *D)
{
	int nprocs,my_rank;
	LIS_INT is,ie;
	LIS_INT *ranges;
	#ifdef USE_MPI
		LIS_INT i;
	#endif

	LIS_DEBUG_FUNC_IN;

	*D = NULL;

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
	if( local_n==0 && global_n==0 )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%D) and global n(=%D) are 0\n",local_n,global_n);
		return LIS_ERR_ILL_ARG;
	}

	*D = (LIS_MATRIX_DIAG)lis_malloc( sizeof(struct LIS_MATRIX_DIAG_STRUCT),"lis_matrix_diag_create::D" );
	if( NULL==*D )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_DIAG_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_matrix_diag_init(D);

	#ifdef USE_MPI
		MPI_Comm_size(comm,&nprocs);
		MPI_Comm_rank(comm,&my_rank);
		ranges = (LIS_INT *)lis_malloc( (nprocs+1)*sizeof(LIS_INT),"lis_matrix_diag_create::ranges" );
		if( ranges==NULL )
		{
			LIS_SETERR_MEM((nprocs+1)*sizeof(LIS_INT));
			lis_matrix_diag_destroy(*D);
			*D = NULL;
			return LIS_OUT_OF_MEMORY;
		}
	#else
		nprocs  = 1;
		my_rank = 0;
		ranges  = NULL;
	#endif

	#ifdef USE_MPI
		MPI_Allreduce(&local_n,&i,1,LIS_MPI_INT,MPI_SUM,comm);
		if( i==0 )
	#else
		if( local_n==0 )
	#endif
	{
		#ifdef USE_MPI
			LIS_GET_ISIE(my_rank,nprocs,global_n,is,ie);
			local_n = ie-is;
			MPI_Allgather(&ie,1,LIS_MPI_INT,&ranges[1],1,LIS_MPI_INT,comm);
			ranges[0] = 0;
		#else
			local_n = global_n;
			is      = 0;
			ie      = global_n;
		#endif
	}
	else
	{
		#ifdef USE_MPI
			MPI_Allgather(&local_n,1,LIS_MPI_INT,&ranges[1],1,LIS_MPI_INT,comm);
			ranges[0] = 0;
			for(i=0;i<nprocs;i++)
			{
				ranges[i+1] += ranges[i];
			}
			global_n = ranges[nprocs];
			is       = ranges[my_rank];
			ie       = ranges[my_rank+1];
		#else
			global_n = local_n;
			is       = 0;
			ie       = local_n;
		#endif
	}
	(*D)->ranges      = ranges;

	(*D)->value = (LIS_SCALAR *)lis_malloc( local_n*sizeof(LIS_SCALAR),"lis_matrix_diag_create::D->value" );
	if( NULL==(*D)->value )
	{
		LIS_SETERR_MEM(local_n*sizeof(LIS_SCALAR));
		lis_matrix_diag_destroy(*D);
		*D = NULL;
		return LIS_OUT_OF_MEMORY;
	}

	(*D)->n           = local_n;
	(*D)->nr          = local_n;
	(*D)->gn          = global_n;
	(*D)->np          = local_n;
	(*D)->comm        = comm;
	(*D)->my_rank     = my_rank;
	(*D)->nprocs      = nprocs;
	(*D)->is          = is;
	(*D)->ie          = ie;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_destroy"
LIS_INT lis_matrix_diag_destroy(LIS_MATRIX_DIAG D)
{
	LIS_INT i;
	LIS_DEBUG_FUNC_IN;

	if( D )
	{
		if( D->value ) lis_free( D->value );
		if( D->work ) lis_free( D->work );
		if( D->bns )
		{
			for(i=0;i<D->nr;i++)
			{
				if( D->v_value[i] ) free(D->v_value[i]);
			}
			lis_free2(2,D->bns,D->v_value);
		}
		if( D->ptr ) lis_free( D->ptr );
		if( D->ranges ) lis_free( D->ranges );
		lis_free(D);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_duplicate"
LIS_INT lis_matrix_diag_duplicate(LIS_MATRIX_DIAG Din, LIS_MATRIX_DIAG *Dout)
{
	LIS_INT err,nr,bnmax,t;
	LIS_INT nprocs;
	LIS_INT i;
	#ifdef USE_MPI
		LIS_INT *ranges;
	#endif

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_diag_check(Din,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	nprocs = Din->nprocs;
	*Dout = NULL;

	*Dout = (LIS_MATRIX_DIAG)lis_malloc( sizeof(struct LIS_MATRIX_DIAG_STRUCT),"lis_matrix_diag_duplicate::Dout" );
	if( NULL==*Dout )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_DIAG_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_matrix_diag_init(Dout);

	#ifdef USE_MPI
		ranges = (LIS_INT *)lis_malloc( (nprocs+1)*sizeof(LIS_INT),"lis_matrix_diag_duplicate::ranges" );
		if( ranges==NULL )
		{
			LIS_SETERR_MEM((nprocs+1)*sizeof(LIS_INT));
			lis_matrix_diag_destroy(*Dout);
			*Dout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		for(i=0;i<nprocs+1;i++) ranges[i] = Din->ranges[i];
		(*Dout)->ranges      = ranges;
	#else
		(*Dout)->ranges      = NULL;
	#endif

	if( Din->bns==NULL )
	{
		(*Dout)->value = (LIS_SCALAR *)lis_malloc( Din->bn*Din->bn*Din->nr*sizeof(LIS_SCALAR),"lis_matrix_diag_duplicate::Dout->value" );
		if( NULL==(*Dout)->value )
		{
			LIS_SETERR_MEM(Din->bn*Din->bn*Din->nr*sizeof(LIS_SCALAR));
			lis_matrix_diag_destroy(*Dout);
			*Dout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		(*Dout)->bn = Din->bn;
	}
	else
	{
		nr = Din->nr;
		(*Dout)->bns = (LIS_INT *)lis_malloc( nr*sizeof(LIS_INT),"lis_matrix_diag_duplicate::Dout->bns" );
		if( NULL==(*Dout)->bns )
		{
			LIS_SETERR_MEM(nr*sizeof(LIS_INT));
			lis_matrix_diag_destroy(*Dout);
			*Dout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		(*Dout)->v_value = (LIS_SCALAR **)lis_malloc( nr*sizeof(LIS_SCALAR *),"lis_matrix_diag_duplicate::Dout->value" );
		if( NULL==(*Dout)->v_value )
		{
			LIS_SETERR_MEM(nr*sizeof(LIS_SCALAR *));
			lis_matrix_diag_destroy(*Dout);
			*Dout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		bnmax = 0;
		for(i=0;i<nr;i++)
		{
			t                 = Din->bns[i]; 
			bnmax             = _max(bnmax,t);
			(*Dout)->bns[i]   = t;
			(*Dout)->v_value[i] = (LIS_SCALAR *)malloc( t*t*sizeof(LIS_SCALAR));
		}
		(*Dout)->bn = bnmax;
		(*Dout)->nr = nr;
	}

	(*Dout)->n           = Din->n;
	(*Dout)->nr          = Din->nr;
	(*Dout)->gn          = Din->gn;
	(*Dout)->np          = Din->np;
	(*Dout)->comm        = Din->comm;
	(*Dout)->my_rank     = Din->my_rank;
	(*Dout)->nprocs      = Din->nprocs;
	(*Dout)->is          = Din->is;
	(*Dout)->ie          = Din->ie;
	(*Dout)->origin      = Din->origin;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_duplicate"
LIS_INT lis_matrix_diag_duplicateM(LIS_MATRIX Ain, LIS_MATRIX_DIAG *Dout)
{
	LIS_INT nr,err;
	LIS_INT nprocs;
	LIS_INT i,k,t,bnmax;
	#ifdef USE_MPI
		LIS_INT *ranges;
	#endif

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	nprocs = Ain->nprocs;
	*Dout = NULL;

	*Dout = (LIS_MATRIX_DIAG)lis_malloc( sizeof(struct LIS_MATRIX_DIAG_STRUCT),"lis_matrix_diag_duplicateM::Dout" );
	if( NULL==*Dout )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_DIAG_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_matrix_diag_init(Dout);

	#ifdef USE_MPI
		ranges = (LIS_INT *)lis_malloc( (nprocs+1)*sizeof(LIS_INT),"lis_matrix_diag_duplicateM::ranges" );
		if( ranges==NULL )
		{
			LIS_SETERR_MEM((nprocs+1)*sizeof(LIS_INT));
			lis_matrix_diag_destroy(*Dout);
			*Dout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		for(i=0;i<nprocs+1;i++) ranges[i] = Ain->ranges[i];
		(*Dout)->ranges      = ranges;
	#else
		(*Dout)->ranges      = NULL;
	#endif

	switch(Ain->matrix_type)
	{
	case LIS_MATRIX_BSR:
	case LIS_MATRIX_BSC:
		k = Ain->nr*Ain->bnr*Ain->bnc;
		(*Dout)->value = (LIS_SCALAR *)lis_malloc( k*sizeof(LIS_SCALAR),"lis_matrix_diag_duplicateM::Dout->value" );
		if( NULL==(*Dout)->value )
		{
			LIS_SETERR_MEM(k*sizeof(LIS_SCALAR));
			lis_matrix_diag_destroy(*Dout);
			*Dout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		(*Dout)->bn = Ain->bnr;
		(*Dout)->nr = Ain->nr;
		break;
	case LIS_MATRIX_VBR:
		nr = Ain->nr;
		(*Dout)->bns = (LIS_INT *)lis_malloc( nr*sizeof(LIS_INT),"lis_matrix_diag_duplicateM::Dout->bns" );
		if( NULL==(*Dout)->bns )
		{
			LIS_SETERR_MEM(nr*sizeof(LIS_INT));
			lis_matrix_diag_destroy(*Dout);
			*Dout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		(*Dout)->v_value = (LIS_SCALAR **)lis_malloc( nr*sizeof(LIS_SCALAR *),"lis_matrix_diag_duplicateM::Dout->value" );
		if( NULL==(*Dout)->v_value )
		{
			LIS_SETERR_MEM(nr*sizeof(LIS_SCALAR *));
			lis_matrix_diag_destroy(*Dout);
			*Dout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		bnmax = 0;
		for(i=0;i<nr;i++)
		{
			t                 = Ain->row[i+1] - Ain->row[i]; 
			bnmax             = _max(bnmax,t);
			(*Dout)->bns[i]   = t;
			(*Dout)->v_value[i] = (LIS_SCALAR *)malloc( t*t*sizeof(LIS_SCALAR));
		}
		(*Dout)->bn = bnmax;
		(*Dout)->nr = nr;
		break;
	default:
		(*Dout)->value = (LIS_SCALAR *)lis_malloc( Ain->np*sizeof(LIS_SCALAR),"lis_matrix_diag_duplicateM::Dout->value" );
		if( NULL==(*Dout)->value )
		{
			LIS_SETERR_MEM(Ain->np*sizeof(LIS_SCALAR));
			lis_matrix_diag_destroy(*Dout);
			*Dout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		(*Dout)->nr          = Ain->n;
		break;
	}

	(*Dout)->n           = Ain->n;
	(*Dout)->gn          = Ain->gn;
	(*Dout)->np          = Ain->np;
	(*Dout)->comm        = Ain->comm;
	(*Dout)->my_rank     = Ain->my_rank;
	(*Dout)->nprocs      = Ain->nprocs;
	(*Dout)->is          = Ain->is;
	(*Dout)->ie          = Ain->ie;
	(*Dout)->origin      = Ain->origin;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_malloc"
LIS_INT lis_matrix_diag_mallocM(LIS_MATRIX A, LIS_SCALAR **diag)
{
	LIS_INT err;
	LIS_INT k;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	switch(A->matrix_type)
	{
	case LIS_MATRIX_BSR:
		k = A->nr*A->bnr*A->bnc;
		*diag = (LIS_SCALAR *)lis_malloc( k*sizeof(LIS_SCALAR),"lis_matrix_diag_mallocM::diag" );
		if( NULL==*diag )
		{
			LIS_SETERR_MEM(k*sizeof(LIS_SCALAR));
			*diag = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		break;
	default:
		*diag = (LIS_SCALAR *)lis_malloc( A->n*sizeof(LIS_SCALAR),"lis_matrix_diag_mallocM::diag" );
		if( NULL==*diag )
		{
			LIS_SETERR_MEM(A->n*sizeof(LIS_SCALAR));
			*diag = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_get_range"
LIS_INT lis_matrix_diag_get_range(LIS_MATRIX_DIAG D, LIS_INT *is, LIS_INT *ie)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_diag_check(D,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;

	*is = D->is;
	*ie = D->ie;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_get_size"
LIS_INT lis_matrix_diag_get_size(LIS_MATRIX_DIAG D, LIS_INT *local_n, LIS_INT *global_n)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_diag_check(D,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;

	*local_n  = D->n;
	*global_n = D->gn;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_set_blocksize"
LIS_INT lis_matrix_diag_set_blocksize(LIS_MATRIX_DIAG D, LIS_INT bn, LIS_INT *bns)
{
	LIS_INT i,n,nr,bnmax,t;
	LIS_INT err;
	LIS_SCALAR *diag;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_diag_check(D,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;

	n = D->n;
	if( bns==NULL )
	{
		nr   = 1 + (n-1)/bn;
		diag = (LIS_SCALAR *)lis_malloc( bn*bn*nr*sizeof(LIS_SCALAR),"lis_matrix_diag_set_blocksize::diag" );
		if( NULL==diag )
		{
			LIS_SETERR_MEM(bn*bn*nr*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		lis_free(D->value);
		D->value = diag;
		D->nr    = nr;
		D->bn    = bn;
	}
	else
	{
		if( D->bns==NULL )
		{
			lis_free(D->value);
			D->bns = (LIS_INT *)lis_malloc( bn*sizeof(LIS_INT),"lis_matrix_diag_duplicateM::Dout->bns" );
			if( NULL==D->bns )
			{
				LIS_SETERR_MEM(bn*sizeof(LIS_INT));
				lis_matrix_diag_destroy(D);
				D = NULL;
				return LIS_OUT_OF_MEMORY;
			}
			D->v_value = (LIS_SCALAR **)lis_malloc( bn*sizeof(LIS_SCALAR *),"lis_matrix_diag_duplicateM::Dout->value" );
			if( NULL==D->v_value )
			{
				LIS_SETERR_MEM(bn*sizeof(LIS_SCALAR *));
				lis_matrix_diag_destroy(D);
				D = NULL;
				return LIS_OUT_OF_MEMORY;
			}
			bnmax = 0;
			for(i=0;i<bn;i++)
			{
				t                 = bns[i]; 
				bnmax             = _max(bnmax,t);
				D->bns[i]   = t;
				D->v_value[i] = (LIS_SCALAR *)malloc( t*t*sizeof(LIS_SCALAR));
			}
			D->bn = bnmax;
			D->nr = bn;
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_copy"
LIS_INT lis_matrix_diag_copy(LIS_MATRIX_DIAG X, LIS_MATRIX_DIAG Y)
{
	LIS_INT i,n,nr,bn;
	LIS_SCALAR *x,*y;

	LIS_DEBUG_FUNC_IN;

	x  = X->value;
	y  = Y->value;
	n  = X->n;
	nr = X->nr;
	if( n!=Y->n )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"length of diagonal matrix X and Y is not equal\n");
		return LIS_ERR_ILL_ARG;
	}
	if( X->bns )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<nr;i++)
		{
			bn = X->bns[i]*X->bns[i];
			memcpy(Y->v_value[i],X->v_value[i],bn*sizeof(LIS_SCALAR));
		}
	}
	else
	{
		bn = X->bn*X->bn;
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0; i<nr; i++)
		{
			memcpy(&y[i*bn],&x[i*bn],bn*sizeof(LIS_SCALAR));
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_scale"
LIS_INT lis_matrix_diag_scale(LIS_SCALAR alpha, LIS_MATRIX_DIAG D)
{
	LIS_INT i,j,nr,bn;
	LIS_SCALAR *d;

	LIS_DEBUG_FUNC_IN;

	d  = D->value;
	nr = D->nr;
	bn = D->bn;
	if( D->bns )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0; i<nr; i++)
		{
			bn = D->bns[i]*D->bns[i];
			for(j=0;j<bn;j++)
			{
				D->v_value[i][j] *= alpha;
			}
		}
	}
	else
	{
		if( bn==1 )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				d[i] = alpha * d[i];
			}
		}
		else if( bn==2 )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				d[4*i]   = alpha * d[4*i];
				d[4*i+1] = alpha * d[4*i+1];
				d[4*i+2] = alpha * d[4*i+2];
				d[4*i+3] = alpha * d[4*i+3];
			}
		}
		else if( bn==3 )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				d[9*i]   = alpha * d[9*i];
				d[9*i+1] = alpha * d[9*i+1];
				d[9*i+2] = alpha * d[9*i+2];
				d[9*i+3] = alpha * d[9*i+3];
				d[9*i+4] = alpha * d[9*i+4];
				d[9*i+5] = alpha * d[9*i+5];
				d[9*i+6] = alpha * d[9*i+6];
				d[9*i+7] = alpha * d[9*i+7];
				d[9*i+8] = alpha * d[9*i+8];
			}
		}
		else if( bn==4 )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				d[16*i]    = alpha * d[16*i];
				d[16*i+1]  = alpha * d[16*i+1];
				d[16*i+2]  = alpha * d[16*i+2];
				d[16*i+3]  = alpha * d[16*i+3];
				d[16*i+4]  = alpha * d[16*i+4];
				d[16*i+5]  = alpha * d[16*i+5];
				d[16*i+6]  = alpha * d[16*i+6];
				d[16*i+7]  = alpha * d[16*i+7];
				d[16*i+8]  = alpha * d[16*i+8];
				d[16*i+9]  = alpha * d[16*i+9];
				d[16*i+10] = alpha * d[16*i+10];
				d[16*i+11] = alpha * d[16*i+11];
				d[16*i+12] = alpha * d[16*i+12];
				d[16*i+13] = alpha * d[16*i+13];
				d[16*i+14] = alpha * d[16*i+14];
				d[16*i+15] = alpha * d[16*i+15];
			}
		}
		else
		{
			bn = bn*bn;
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				for(j=0;j<bn;j++)
				{
					d[i*bn+j] = alpha * d[i*bn+j];
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_inverse"
LIS_INT lis_matrix_diag_inverse(LIS_MATRIX_DIAG D)
{
	LIS_INT i,k,nr,bn,bs,n;
	LIS_SCALAR *d;

	LIS_DEBUG_FUNC_IN;

	n  = D->n;
	d  = D->value;
	nr = D->nr;
	bn = D->bn;
	bs = D->bn*D->bn;
	if( D->bns )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,bn)
		#endif
		for(i=0; i<nr; i++)
		{
			bn = D->bns[i];
			lis_array_ge(bn,D->v_value[i]);
		}
	}
	else
	{
		if( bn==1 )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				d[i] = 1.0 / d[i];
			}
		}
		else
		{
			k = n%bn;
			if( k!=0 )
			{
				for(i=bn-1;i>=k;i--)
				{
					d[bs*(nr-1)+i*(bn+1)] = 1.0;
				}
			}
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<nr;i++)
			{
				lis_array_ge(bn,&d[i*bs]);
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_matvec"
LIS_INT lis_matrix_diag_matvec(LIS_MATRIX_DIAG D, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT i,nr,bn,bs;
	LIS_SCALAR *d,*x,*y;

	LIS_DEBUG_FUNC_IN;

	d  = D->value;
	x  = X->value;
	y  = Y->value;
	nr = D->nr;
	bn = D->bn;
	bs = D->bn*D->bn;
	if( D->bns )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,bn)
		#endif
		for(i=0; i<nr; i++)
		{
			bn = D->bns[i];
			lis_array_matvec(bn,D->v_value[i],&x[i*bn],&y[i*bn],LIS_INS_VALUE);
		}
	}
	else
	{
		if( bn==1 )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				y[i] = x[i] * d[i];
			}
		}
		else if( bn==2 )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				y[2*i]   = d[4*i]   * x[2*i] + d[4*i+2] * x[2*i+1];
				y[2*i+1] = d[4*i+1] * x[2*i] + d[4*i+3] * x[2*i+1];
			}
		}
		else if( bn==3 )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				y[3*i]   = d[9*i]   * x[3*i] + d[9*i+3] * x[3*i+1] + d[9*i+6] * x[3*i+2];
				y[3*i+1] = d[9*i+1] * x[3*i] + d[9*i+4] * x[3*i+1] + d[9*i+7] * x[3*i+2];
				y[3*i+2] = d[9*i+2] * x[3*i] + d[9*i+5] * x[3*i+1] + d[9*i+8] * x[3*i+2];
			}
		}
		else if( bn==4 )
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				y[4*i]   = d[16*i]   * x[4*i] + d[16*i+4] * x[4*i+1] + d[16*i+8]  * x[4*i+2] + d[16*i+12] * x[4*i+3];
				y[4*i+1] = d[16*i+1] * x[4*i] + d[16*i+5] * x[4*i+1] + d[16*i+9]  * x[4*i+2] + d[16*i+13] * x[4*i+3];
				y[4*i+2] = d[16*i+2] * x[4*i] + d[16*i+6] * x[4*i+1] + d[16*i+10] * x[4*i+2] + d[16*i+14] * x[4*i+3];
				y[4*i+3] = d[16*i+3] * x[4*i] + d[16*i+7] * x[4*i+1] + d[16*i+11] * x[4*i+2] + d[16*i+15] * x[4*i+3];
			}
		}
		else
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				lis_array_matvec(bn,&d[i*bs],&x[i*bn],&y[i*bn],LIS_INS_VALUE);
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_matvech"
LIS_INT lis_matrix_diag_matvech(LIS_MATRIX_DIAG D, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT i,nr,bn,bs;
	LIS_SCALAR *d,*x,*y;

	LIS_DEBUG_FUNC_IN;

	d  = D->value;
	x  = X->value;
	y  = Y->value;
	nr = D->nr;
	bn = D->bn;
	bs = D->bn*D->bn;
	if( D->bns )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,bn)
		#endif
		for(i=0; i<nr; i++)
		{
			bn = D->bns[i];
			lis_array_matvech(bn,D->v_value[i],&x[i*bn],&y[i*bn],LIS_INS_VALUE);
		}
	}
	else
	{
#if 0
		if( bn==1 )
		{
			#pragma omp parallel for private(i)
			for(i=0; i<nr; i++)
			{
				y[i] = x[i] * conj(d[i]);
			}
		}
		else if( bn==2 )
		{
			#pragma omp parallel for private(i)
			for(i=0; i<nr; i++)
			{
				y[2*i]   = conj(d[4*i])   * x[2*i] + conj(d[4*i+1]) * x[2*i+1];
				y[2*i+1] = conj(d[4*i+2]) * x[2*i] + conj(d[4*i+3]) * x[2*i+1];
			}
		}
		else if( bn==3 )
		{
			#pragma omp parallel for private(i)
			for(i=0; i<nr; i++)
			{
				y[3*i]   = conj(d[9*i])   * x[3*i] + conj(d[9*i+1]) * x[3*i+1] + conj(d[9*i+2]) * x[3*i+2];
				y[3*i+1] = conj(d[9*i+3]) * x[3*i] + conj(d[9*i+4]) * x[3*i+1] + conj(d[9*i+5]) * x[3*i+2];
				y[3*i+2] = conj(d[9*i+6]) * x[3*i] + conj(d[9*i+7]) * x[3*i+1] + conj(d[9*i+8]) * x[3*i+2];
			}
		}
		else if( bn==4 )
		{
			#pragma omp parallel for private(i)
			for(i=0; i<nr; i++)
			{
				y[4*i]   = conj(d[16*i]     * x[4*i] + conj(d[16*i+ 1]) * x[4*i+1] + conj(d[16*i+ 2]) * x[4*i+2] + conj(d[16*i+ 3]) * x[4*i+3];
				y[4*i+1] = conj(d[16*i+ 4]) * x[4*i] + conj(d[16*i+ 5]) * x[4*i+1] + conj(d[16*i+ 6]) * x[4*i+2] + conj(d[16*i+ 7]) * x[4*i+3];
				y[4*i+2] = conj(d[16*i+ 8]) * x[4*i] + conj(d[16*i+ 9]) * x[4*i+1] + conj(d[16*i+10]) * x[4*i+2] + conj(d[16*i+11]) * x[4*i+3];
				y[4*i+3] = conj(d[16*i+12]) * x[4*i] + conj(d[16*i+13]) * x[4*i+1] + conj(d[16*i+14]) * x[4*i+2] + conj(d[16*i+15]) * x[4*i+3];
			}
		}
		else
#endif
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<nr; i++)
			{
				lis_array_matvech(bn,&d[i*bs],&x[i*bn],&y[i*bn],LIS_INS_VALUE);
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_diag_print"
LIS_INT lis_matrix_diag_print(LIS_MATRIX_DIAG D)
{
	LIS_INT	k,i,j,nr,bn,nn;

	LIS_DEBUG_FUNC_IN;

	nr = D->nr;
	if( D->bns )
	{
		nn = 0;
		for(k=0; k<nr; k++)
		{
			bn = D->bns[k];
			for(i=0;i<bn;i++)
			{
#ifdef _LONG__LONG
				printf("%4lld (", nn+i);
#else
				printf("%4d (", nn+i);
#endif
				for(j=0; j<bn; j++)
				{
#ifdef _COMPLEX				  
					printf("%6.2f %6.2f ", (double)creal(D->v_value[k][j*bn+i]), (double)cimag(D->v_value[k][j*bn+i]));
#else
					printf("%6.2f ", (double)D->v_value[k][j*bn+i]);
#endif					
				}
				printf(")\n");
			}
			nn += bn;
		}
	}
	else
	{
		bn = D->bn;
		nn = D->bn*D->bn;
		for(k=0; k<nr; k++)
		{
			for(i=0;i<bn;i++)
			{
#ifdef _LONG__LONG
				printf("%4lld (", k*nn+i);
#else
				printf("%4d (", k*nn+i);
#endif
				for(j=0; j<bn; j++)
				{
#ifdef _COMPLEX
					printf("%6.2f %6.2f ", (double)creal(D->value[k*nn + j*bn+i]), (double)cimag(D->value[k*nn + j*bn+i]));
#else				  
					printf("%6.2f ", (double)D->value[k*nn + j*bn+i]);
#endif					
				}
				printf(")\n");
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
