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


/*
 * This subroutine is made based on ITSOL.
 *
 * http://www-users.cs.umn.edu/~saad/software/ITSOL/
 *
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
#include <math.h>
#include <memory.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_precon_create
 * lis_precon_psd_create
 * lis_precon_psd_update
 * lis_psolve
 * lis_psolveh
 ************************************************/

LIS_INT lis_symbolic_fact_csr(LIS_SOLVER solver, LIS_PRECON precon);
LIS_INT lis_numerical_fact_csr(LIS_SOLVER solver, LIS_PRECON precon);
LIS_INT lis_symbolic_fact_bsr(LIS_SOLVER solver, LIS_PRECON precon);
LIS_INT lis_numerical_fact_bsr(LIS_SOLVER solver, LIS_PRECON precon);
LIS_INT lis_symbolic_fact_vbr(LIS_SOLVER solver, LIS_PRECON precon);
LIS_INT lis_numerical_fact_vbr(LIS_SOLVER solver, LIS_PRECON precon);

#undef __FUNC__
#define __FUNC__ "lis_precon_create_iluk"
LIS_INT lis_precon_create_iluk(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	storage,err;
	LIS_MATRIX A,B;

	LIS_DEBUG_FUNC_IN;

	storage     = solver->options[LIS_OPTIONS_STORAGE];

	if( storage==LIS_MATRIX_BSR || storage==LIS_MATRIX_VBR )
	{
		if( solver->A->matrix_type!=storage )
		{
			err = lis_matrix_convert_self(solver);
			if( err ) return err;
		}
	}


	switch( solver->A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_symbolic_fact_csr(solver,precon);
		if( err ) return err;
		err = lis_numerical_fact_csr(solver,precon);
		if( err ) return err;
		lis_psolve_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolve_iluk_csr;
		lis_psolveh_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolveh_iluk_csr;
		precon->is_copy = LIS_TRUE;
		break;
	case LIS_MATRIX_BSR:
		err = lis_symbolic_fact_bsr(solver,precon);
		if( err ) return err;
		err = lis_numerical_fact_bsr(solver,precon);
		if( err ) return err;
		lis_psolve_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolve_iluk_bsr;
		lis_psolveh_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolveh_iluk_bsr;
		break;
	case LIS_MATRIX_VBR:
		err = lis_symbolic_fact_vbr(solver,precon);
		if( err ) return err;
		err = lis_numerical_fact_vbr(solver,precon);
		if( err ) return err;
		lis_psolve_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolve_iluk_vbr;
/*		lis_psolveh_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolveh_iluk_vbr;*/
		break;
	default:
		A = solver->A;
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,LIS_MATRIX_CSR);
		err = lis_matrix_convert(A,B);
		if( err ) return err;
		solver->A = B;
		err = lis_symbolic_fact_csr(solver,precon);
		if( err ) return err;
		err = lis_numerical_fact_csr(solver,precon);
		if( err ) return err;
		lis_psolve_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolve_iluk_csr;
		lis_psolveh_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolveh_iluk_csr;
		lis_matrix_destroy(B);
		solver->A = A;
		precon->is_copy = LIS_TRUE;
		break;
	}

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_precon_psd_create_iluk"
LIS_INT lis_precon_psd_create_iluk(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	storage,err;
	LIS_MATRIX A,B;

	LIS_DEBUG_FUNC_IN;

	storage     = solver->options[LIS_OPTIONS_STORAGE];

	if( storage==LIS_MATRIX_BSR || storage==LIS_MATRIX_VBR )
	{
		if( solver->A->matrix_type!=storage )
		{
			err = lis_matrix_convert_self(solver);
			if( err ) return err;
		}
	}

	switch( solver->A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_symbolic_fact_csr(solver,precon);
		if( err ) return err;
        lis_psolve_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolve_iluk_csr;
        lis_psolveh_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolveh_iluk_csr;
        precon->is_copy = LIS_TRUE;
		break;
	case LIS_MATRIX_BSR:
/*        err = lis_symbolic_fact_bsr(solver,precon);*/
/*        if( err ) return err;*/
/*        err = lis_numerical_fact_bsr(solver,precon);*/
/*        if( err ) return err;*/
/*        lis_psolve_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolve_iluk_bsr;*/
/*        lis_psolveh_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolveh_iluk_bsr;*/
/*        break;*/
        err = LIS_ERR_NOT_IMPLEMENTED;
		return err;
	case LIS_MATRIX_VBR:
/*        err = lis_symbolic_fact_vbr(solver,precon);*/
/*        if( err ) return err;*/
/*        err = lis_numerical_fact_vbr(solver,precon);*/
/*        if( err ) return err;*/
/*        lis_psolve_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolve_iluk_vbr;*/
/*		lis_psolveh_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolveh_iluk_vbr;*/
/*        break;*/
        err = LIS_ERR_NOT_IMPLEMENTED;
		return err;
	default:
		A = solver->A;
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,LIS_MATRIX_CSR);
		err = lis_matrix_convert(A,B);
		if( err ) return err;
		solver->A = B;
		err = lis_symbolic_fact_csr(solver,precon);
		if( err ) return err;
		lis_psolve_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolve_iluk_csr;
		lis_psolveh_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolveh_iluk_csr;
		lis_matrix_destroy(B);
		solver->A = A;
		precon->is_copy = LIS_TRUE;
		break;
	}

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_precon_psd_update_iluk"
LIS_INT lis_precon_psd_update_iluk(LIS_SOLVER solver, LIS_PRECON precon)
{
/*    LIS_INT	storage,err;*/
	LIS_INT	err;
	LIS_MATRIX A;

	LIS_DEBUG_FUNC_IN;

	switch( solver->A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_numerical_fact_csr(solver,precon);
		if( err ) return err;
		break;
	case LIS_MATRIX_BSR:
/*        err = lis_symbolic_fact_bsr(solver,precon);*/
/*        if( err ) return err;*/
/*        err = lis_numerical_fact_bsr(solver,precon);*/
/*        if( err ) return err;*/
/*        lis_psolve_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolve_iluk_bsr;*/
/*        lis_psolveh_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolveh_iluk_bsr;*/
/*        break;*/
        err = LIS_ERR_NOT_IMPLEMENTED;
		return err;
	case LIS_MATRIX_VBR:
/*        err = lis_symbolic_fact_vbr(solver,precon);*/
/*        if( err ) return err;*/
/*        err = lis_numerical_fact_vbr(solver,precon);*/
/*        if( err ) return err;*/
/*        lis_psolve_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolve_iluk_vbr;*/
/*		lis_psolveh_xxx[LIS_PRECON_TYPE_ILU]  = lis_psolveh_iluk_vbr;*/
/*        break;*/
        err = LIS_ERR_NOT_IMPLEMENTED;
		return err;
	default:
/*        nothing should ever really happen here, as any undetermined matrix_type*/
/*        should have been handled in "psd_create"*/
		break;
	}

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_symbolic_fact_csr"
LIS_INT lis_symbolic_fact_csr(LIS_SOLVER solver, LIS_PRECON precon)
{
#ifdef _OPENMP
	LIS_INT	err;
	LIS_INT	i,j,k;
	LIS_INT	n,levfill;
	LIS_INT	col,ip,it,jpiv,incl,incu,jmin,kmin;
	LIS_INT	*levls,*jbuf,*iw,**ulvl;
	LIS_INT	is,ie,my_rank,nprocs;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	levfill = solver->options[LIS_OPTIONS_FILL];
	nprocs = omp_get_max_threads();

	L      = NULL;
	U      = NULL;


	err = lis_matrix_ilu_create(n,1,&L);
	if( err ) return err;
	err = lis_matrix_ilu_create(n,1,&U);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(L);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(U);
	if( err ) return err;
	err = lis_vector_duplicate(A,&D);
	if( err )
	{
		return err;
	}


	ulvl   = (LIS_INT **)lis_malloc(n*sizeof(LIS_INT *),"lis_symbolic_fact_csr::ulvl");
	if( ulvl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	levls   = (LIS_INT *)lis_malloc(nprocs*n*sizeof(LIS_INT),"lis_symbolic_fact_csr::levls");
	if( levls==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	jbuf   = (LIS_INT *)lis_malloc(nprocs*n*sizeof(LIS_INT),"lis_symbolic_fact_csr::jbuf");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_symbolic_fact_csr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	#pragma omp parallel private(i,j,k,is,ie,my_rank,incl,incu,col,jpiv,it,ip,kmin,jmin)
	{
		my_rank  = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		for(i=is;i<ie;i++) iw[i]=-1;

		for(i=is;i<ie;i++)
		{
			incl = 0;
			incu = i;

			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				col = A->index[j];
				#ifdef USE_MPI
					if( col>=n ) continue;
				#endif
				if( col>=is && col<ie )
				{
					if( col < i )
					{
						jbuf[my_rank*n+incl] = col;
						levls[my_rank*n+incl] = 0;
						iw[col] = incl++;
					}
					else if( col > i )
					{
						jbuf[my_rank*n+incu] = col;
						levls[my_rank*n+incu] = 0;
						iw[col] = incu++;
					}
				}
			}

			jpiv = -1;
			while( ++jpiv < incl )
			{
				k = jbuf[my_rank*n+jpiv];
				kmin = k;
				jmin = jpiv;
				for(j=jpiv+1;j<incl;j++)
				{
					if( jbuf[my_rank*n+j]<kmin )
					{
						kmin = jbuf[my_rank*n+j];
						jmin = j;
					}
				}

				if( jmin!=jpiv )
				{
					jbuf[my_rank*n+jpiv] = kmin;
					jbuf[my_rank*n+jmin] = k;
					iw[kmin] = jpiv;
					iw[k] = jmin;
					j = levls[my_rank*n+jpiv];
					levls[my_rank*n+jpiv] = levls[my_rank*n+jmin];
					levls[my_rank*n+jmin] = j;
					k = kmin;
				}

				for(j=0;j<U->nnz[k];j++)
				{
					col = U->index[k][j];
					it = ulvl[k][j] + levls[my_rank*n+jpiv]+1;
					if( it > levfill ) continue;
					ip = iw[col];
					if( ip==-1 )
					{
						if( col < i )
						{
							jbuf[my_rank*n+incl] = col;
							levls[my_rank*n+incl] = it;
							iw[col] = incl++;
						}
						else if( col > i )
						{
							jbuf[my_rank*n+incu] = col;
							levls[my_rank*n+incu] = it;
							iw[col] = incu++;
						}
					}
					else
					{
						levls[my_rank*n+ip] = _min(levls[my_rank*n+ip],it);
					}
				}
			}
			for(j=0;j<incl;j++) iw[jbuf[my_rank*n+j]] = -1;
			for(j=i;j<incu;j++) iw[jbuf[my_rank*n+j]] = -1;

			L->nnz[i] = incl;
			if( incl > 0 )
			{
				L->index[i] = (LIS_INT *)malloc(incl*sizeof(LIS_INT));
				L->value[i] = (LIS_SCALAR *)malloc(incl*sizeof(LIS_SCALAR));
				memcpy(L->index[i],&jbuf[my_rank*n],incl*sizeof(LIS_INT));
			}

			k = incu-i;
			U->nnz[i] = k;
			if( k > 0 )
			{
				U->index[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
				U->value[i] = (LIS_SCALAR *)malloc(k*sizeof(LIS_SCALAR));
				ulvl[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
				memcpy(U->index[i],&jbuf[my_rank*n+i],k*sizeof(LIS_INT));
				memcpy(ulvl[i],&levls[my_rank*n+i],k*sizeof(LIS_INT));
			}
		}
	}

	precon->L  = L;
	precon->U  = U;
	precon->D  = D;

	lis_free2(3,levls,jbuf,iw);
	for(i=0;i<n-1;i++)
	{
		if(U->nnz[i]) free(ulvl[i]);
	}
	lis_free(ulvl);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT	err;
	LIS_INT	i,j,k;
	LIS_INT	n,levfill;
	LIS_INT	col,ip,it,jpiv,incl,incu,jmin,kmin;
	LIS_INT	*levls,*jbuf,*iw,**ulvl;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	levfill = solver->options[LIS_OPTIONS_FILL];

	L      = NULL;
	U      = NULL;
	D      = NULL;


	err = lis_matrix_ilu_create(n,1,&L);
	if( err ) return err;
	err = lis_matrix_ilu_create(n,1,&U);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(L);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(U);
	if( err ) return err;
	err = lis_vector_duplicate(A,&D);
	if( err )
	{
		return err;
	}


	ulvl   = (LIS_INT **)lis_malloc(n*sizeof(LIS_INT *),"lis_symbolic_fact_csr::ulvl");
	if( ulvl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	levls   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_symbolic_fact_csr::levls");
	if( levls==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	jbuf   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_symbolic_fact_csr::jbuf");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_symbolic_fact_csr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<n;i++) iw[i]=-1;

	for(i=0;i<n;i++)
	{
		incl = 0;
		incu = i;

		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			col = A->index[j];
			#ifdef USE_MPI
				if( col>=n ) continue;
			#endif
			if( col < i )
			{
				jbuf[incl] = col;
				levls[incl] = 0;
				iw[col] = incl++;
			}
			else if( col > i )
			{
				jbuf[incu] = col;
				levls[incu] = 0;
				iw[col] = incu++;
			}
		}

		jpiv = -1;
		while( ++jpiv < incl )
		{
			k = jbuf[jpiv];
			kmin = k;
			jmin = jpiv;
			for(j=jpiv+1;j<incl;j++)
			{
				if( jbuf[j]<kmin )
				{
					kmin = jbuf[j];
					jmin = j;
				}
			}

			if( jmin!=jpiv )
			{
				jbuf[jpiv] = kmin;
				jbuf[jmin] = k;
				iw[kmin] = jpiv;
				iw[k] = jmin;
				j = levls[jpiv];
				levls[jpiv] = levls[jmin];
				levls[jmin] = j;
				k = kmin;
			}

			for(j=0;j<U->nnz[k];j++)
			{
				col = U->index[k][j];
				it = ulvl[k][j] + levls[jpiv]+1;
				if( it > levfill ) continue;
				ip = iw[col];
				if( ip==-1 )
				{
					if( col < i )
					{
						jbuf[incl] = col;
						levls[incl] = it;
						iw[col] = incl++;
					}
					else if( col > i )
					{
						jbuf[incu] = col;
						levls[incu] = it;
						iw[col] = incu++;
					}
				}
				else
				{
					levls[ip] = _min(levls[ip],it);
				}
			}
		}
		for(j=0;j<incl;j++) iw[jbuf[j]] = -1;
		for(j=i;j<incu;j++) iw[jbuf[j]] = -1;

		L->nnz[i] = incl;
		if( incl > 0 )
		{
			L->index[i] = (LIS_INT *)malloc(incl*sizeof(LIS_INT));
			L->value[i] = (LIS_SCALAR *)malloc(incl*sizeof(LIS_SCALAR));
			memcpy(L->index[i],jbuf,incl*sizeof(LIS_INT));
		}

		k = incu-i;
		U->nnz[i] = k;
		if( k > 0 )
		{
			U->index[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
			U->value[i] = (LIS_SCALAR *)malloc(k*sizeof(LIS_SCALAR));
			ulvl[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
			memcpy(U->index[i],jbuf+i,k*sizeof(LIS_INT));
			memcpy(ulvl[i],levls+i,k*sizeof(LIS_INT));
		}
	}

	precon->L  = L;
	precon->U  = U;
	precon->D  = D;

	lis_free2(3,levls,jbuf,iw);
	for(i=0;i<n-1;i++)
	{
		if(U->nnz[i]) free(ulvl[i]);
	}
	lis_free(ulvl);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_numerical_fact_csr"
LIS_INT lis_numerical_fact_csr(LIS_SOLVER solver, LIS_PRECON precon)
{
#ifdef _OPENMP
	LIS_INT	i,j,k;
	LIS_INT	n;
	LIS_INT	col,jpos,jrow;
	LIS_INT	*jw;
	LIS_INT	is,ie,my_rank,nprocs;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nprocs = omp_get_max_threads();

	L = precon->L;
	U = precon->U;
	D = precon->D;


	jw   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_numerical_fact_csr::jw");
	if( jw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	#pragma omp parallel private(i,j,k,is,ie,my_rank,col,jpos,jrow)
	{
		my_rank  = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		for(i=is;i<ie;i++) jw[i] = -1;


		for(i=is;i<ie;i++)
		{
			for(j=0;j<L->nnz[i];j++)
			{
				col = L->index[i][j];
				jw[col] = j;
				L->value[i][j] = 0;
			}

			jw[i] = i;
			D->value[i] = 0;

			for(j=0;j<U->nnz[i];j++)
			{
				col = U->index[i][j];
				jw[col] = j;
				U->value[i][j] = 0;
			}

			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				col = A->index[j];
				#ifdef USE_MPI
					if( col>=n ) continue;
				#endif
				if( col>=is && col<ie )
				{
					jpos = jw[col];
					if( col<i )
					{
						L->value[i][jpos] = A->value[j];
					}
					else if( col==i )
					{
						D->value[i] = A->value[j];
					}
					else
					{
						U->value[i][jpos] = A->value[j];
					}
				}
			}

			for(j=0;j<L->nnz[i];j++)
			{
				jrow = L->index[i][j];
				L->value[i][j] *= D->value[jrow];

				for(k=0;k<U->nnz[jrow];k++)
				{
					col = U->index[jrow][k];
					jpos = jw[col];
					if( jpos==-1 ) continue;
					if( col<i )
					{
						L->value[i][jpos] -= L->value[i][j] * U->value[jrow][k];
					}
					else if( col==i )
					{
						D->value[i] -= L->value[i][j] * U->value[jrow][k];
					}
					else
					{
						U->value[i][jpos] -= L->value[i][j] * U->value[jrow][k];
					}
				}
			}

			for(j=0;j<L->nnz[i];j++)
			{
				col = L->index[i][j];
				jw[col] = -1;
			}
			jw[i] = -1;
			for(j=0;j<U->nnz[i];j++)
			{
				col = U->index[i][j];
				jw[col] = -1;
			}
			D->value[i] = 1.0 / D->value[i];
		}
	}
	lis_free(jw);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT	i,j,k;
	LIS_INT	n;
	LIS_INT	col,jpos,jrow;
	LIS_INT	*jw;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;

	L = precon->L;
	U = precon->U;
	D = precon->D;


	jw   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_numerical_fact_csr::jw");
	if( jw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}


	for(i=0;i<n;i++) jw[i] = -1;

	for(i=0;i<n;i++)
	{
		for(j=0;j<L->nnz[i];j++)
		{
			col = L->index[i][j];
			jw[col] = j;
			L->value[i][j] = 0;
		}

		jw[i] = i;
		D->value[i] = 0;

		for(j=0;j<U->nnz[i];j++)
		{
			col = U->index[i][j];
			jw[col] = j;
			U->value[i][j] = 0;
		}

		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			col = A->index[j];
			#ifdef USE_MPI
				if( col>=n ) continue;
			#endif
			jpos = jw[col];
			if( col<i )
			{
				L->value[i][jpos] = A->value[j];
			}
			else if( col==i )
			{
				D->value[i] = A->value[j];
			}
			else
			{
				U->value[i][jpos] = A->value[j];
			}
		}

		for(j=0;j<L->nnz[i];j++)
		{
			jrow = L->index[i][j];
			L->value[i][j] *= D->value[jrow];

			for(k=0;k<U->nnz[jrow];k++)
			{
				col = U->index[jrow][k];
				jpos = jw[col];
				if( jpos==-1 ) continue;
				if( col<i )
				{
					L->value[i][jpos] -= L->value[i][j] * U->value[jrow][k];
				}
				else if( col==i )
				{
					D->value[i] -= L->value[i][j] * U->value[jrow][k];
				}
				else
				{
					U->value[i][jpos] -= L->value[i][j] * U->value[jrow][k];
				}
			}
		}

		for(j=0;j<L->nnz[i];j++)
		{
			col = L->index[i][j];
			jw[col] = -1;
		}
		jw[i] = -1;
		for(j=0;j<U->nnz[i];j++)
		{
			col = U->index[i][j];
			jw[col] = -1;
		}
		D->value[i] = 1.0 / D->value[i];
	}
	lis_free(jw);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_iluk_csr"
LIS_INT lis_psolve_iluk_csr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
#ifdef _OPENMP
	LIS_INT i,j,jj,n;
	LIS_INT is,ie,my_rank,nprocs;
	LIS_SCALAR *x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_SCALAR *xl;
	#endif

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif
	n = solver->A->n;
	nprocs = omp_get_max_threads();

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			#pragma omp parallel private(i,j,jj,is,ie,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

				for(i=is;i<ie;i++)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj     = L->index[i][j];
						x[i] -= L->value[i][j] * x[jj];
					}
				}
				for(i=ie-1;i>=is;i--)
				{
					for(j=0;j<U->nnz[i];j++)
					{
						jj = U->index[i][j];
						x[i] -= U->value[i][j] * x[jj];
					}
					x[i] = D->value[i]*x[i];
				}
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,X);
			nprocs = omp_get_max_threads();
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

				for(i=is;i<ie;i++)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj     = L->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-L->value[i][j]);
						#else
							LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-L->value[i][j]);
						#endif
/*						x[i] -= L->value[i][j] * x[jj];*/
					}
				}
				for(i=ie-1;i>=is;i--)
				{
					for(j=0;j<U->nnz[i];j++)
					{
						jj = U->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-U->value[i][j]);
						#else
							LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-U->value[i][j]);
						#endif
/*						x[i] -= U->value[i][j] * x[jj];*/
					}
					#ifndef USE_SSE2
						LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],D->value[i]);
					#else
						LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],D->value[i]);
					#endif
/*					x[i] = D->value[i]*x[i];*/
				}
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT i,j,jj,n;
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_SCALAR *xl;
	#endif

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	b = B->value;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif
	n = solver->A->n;

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			for(i=0; i<n; i++)
			{
				for(j=0;j<L->nnz[i];j++)
				{
					jj     = L->index[i][j];
					x[i] -= L->value[i][j] * x[jj];
				}
			}
			for(i=n-1; i>=0; i--)
			{
				for(j=0;j<U->nnz[i];j++)
				{
					jj = U->index[i][j];
					x[i] -= U->value[i][j] * x[jj];
				}
				x[i] = D->value[i]*x[i];
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copy(B,X);
			for(i=0; i<n; i++)
			{
				for(j=0;j<L->nnz[i];j++)
				{
					jj     = L->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-L->value[i][j]);
					#else
						LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-L->value[i][j]);
					#endif
/*					x[i] -= L->value[i][j] * x[jj];*/
				}
			}
			for(i=n-1; i>=0; i--)
			{
				for(j=0;j<U->nnz[i];j++)
				{
					jj = U->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-U->value[i][j]);
					#else
						LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-U->value[i][j]);
					#endif
/*					x[i] -= U->value[i][j] * x[jj];*/
				}
				#ifndef USE_SSE2
					LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],D->value[i]);
				#else
					LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],D->value[i]);
				#endif
/*				x[i] = D->value[i]*x[i];*/
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_iluk_csr"
LIS_INT lis_psolveh_iluk_csr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
#ifdef _OPENMP
	LIS_INT i,j,jj,n;
	LIS_INT is,ie,my_rank,nprocs;
	LIS_SCALAR *x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_SCALAR *xl;
	#endif

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif
	n = solver->A->n;
	nprocs = omp_get_max_threads();

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			#pragma omp parallel private(i,j,jj,is,ie,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

				for(i=is;i<ie;i++)
				{
					x[i] = D->value[i]*x[i];
					for(j=0;j<U->nnz[i];j++)
					{
						jj     = U->index[i][j];
						x[jj] -= conj(U->value[i][j]) * x[i];
					}
				}
				for(i=ie-1;i>=is;i--)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj     = L->index[i][j];
						x[jj] -= conj(L->value[i][j]) * x[i];
					}
				}
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,X);
			nprocs = omp_get_max_threads();
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

				for(i=is;i<ie;i++)
				{
					#ifndef USE_SSE2
			 			LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],conj(D->value[i]));
					#else
						LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],conj(D->value[i]));
					#endif
/*					x[i] = conj(D->value[i])*x[i];*/
					for(j=0;j<U->nnz[i];j++)
					{
						jj     = U->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(U->value[i][j]));
						#else
							LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(U->value[i][j]));
						#endif
/*						x[jj] -= U->value[i][j] * x[i];*/
					}
				}
				for(i=ie-1;i>=is;i--)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj     = L->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
						#else
							LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
						#endif
/*						x[jj] -= L->value[i][j] * x[i];*/
					}
				}
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT i,j,jj,n;
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_SCALAR *xl;
	#endif


	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	b = B->value;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif
	n = solver->A->n;

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			for(i=0; i<n; i++)
			{
				x[i] = conj(D->value[i])*x[i];
				for(j=0;j<U->nnz[i];j++)
				{
					jj     = U->index[i][j];
					x[jj] -= conj(U->value[i][j]) * x[i];
				}
			}
			for(i=n-1; i>=0; i--)
			{
				for(j=0;j<L->nnz[i];j++)
				{
					jj     = L->index[i][j];
					x[jj] -= conj(L->value[i][j]) * x[i];
				}
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copy(B,X);
			for(i=0; i<n; i++)
			{
				#ifndef USE_SSE2
					LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],conj(D->value[i]));
				#else
					LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],conj(D->value[i]));
				#endif
/*				x[i] = conj(D->value[i]*x[i]);*/
				for(j=0;j<U->nnz[i];j++)
				{
					jj     = U->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(U->value[i][j]));
					#else
						LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(U->value[i][j]));
					#endif
/*					x[jj] -= conj(U->value[i][j]) * x[i];*/
				}
			}
			for(i=n-1; i>=0; i--)
			{
				for(j=0;j<L->nnz[i];j++)
				{
					jj     = L->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(L->value[i][j]));
					#else
						LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(L->value[i][j]));
					#endif
/*					x[jj] -= conj(L->value[i][j]) * x[i];*/
				}
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}


#undef __FUNC__
#define __FUNC__ "lis_symbolic_fact_bsr"
LIS_INT lis_symbolic_fact_bsr(LIS_SOLVER solver, LIS_PRECON precon)
{
#ifdef _OPENMP
	LIS_INT	err;
	LIS_INT	i,j,k,bnr,bs;
	LIS_INT	n,nr,levfill;
	LIS_INT	col,ip,it,jpiv,incl,incu,jmin,kmin;
	LIS_INT	*levls,*jbuf,*iw,**ulvl;
	LIS_INT	is,ie,my_rank,nprocs;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nr     = A->nr;
	bnr    = A->bnr;
	bs     = bnr*bnr;
	levfill = solver->options[LIS_OPTIONS_FILL];
	nprocs = omp_get_max_threads();

	L      = NULL;
	U      = NULL;


	err = lis_matrix_ilu_create(nr,bnr,&L);
	if( err ) return err;
	err = lis_matrix_ilu_create(nr,bnr,&U);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(L);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(U);
	if( err ) return err;
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		return err;
	}


	ulvl   = (LIS_INT **)lis_malloc(nr*sizeof(LIS_INT *),"lis_symbolic_fact_bsr::ulvl");
	if( ulvl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	levls   = (LIS_INT *)lis_malloc(nprocs*nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::levls");
	if( levls==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	jbuf   = (LIS_INT *)lis_malloc(nprocs*nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::jbuf");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	#pragma omp parallel private(i,j,k,is,ie,my_rank,incl,incu,col,jpiv,it,ip,kmin,jmin)
	{
		my_rank  = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,nr,is,ie);

		for(i=is;i<ie;i++) iw[i]=-1;

		for(i=is;i<ie;i++)
		{
			incl = 0;
			incu = i;

			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				col = A->bindex[j];
				#ifdef USE_MPI
					if( col>=nr ) continue;
				#endif
				if( col>=is && col<ie )
				{
					if( col < i )
					{
						jbuf[my_rank*nr+incl] = col;
						levls[my_rank*nr+incl] = 0;
						iw[col] = incl++;
					}
					else if( col > i )
					{
						jbuf[my_rank*nr+incu] = col;
						levls[my_rank*nr+incu] = 0;
						iw[col] = incu++;
					}
				}
			}

			jpiv = -1;
			while( ++jpiv < incl )
			{
				k = jbuf[my_rank*nr+jpiv];
				kmin = k;
				jmin = jpiv;
				for(j=jpiv+1;j<incl;j++)
				{
					if( jbuf[my_rank*nr+j]<kmin )
					{
						kmin = jbuf[my_rank*nr+j];
						jmin = j;
					}
				}

				if( jmin!=jpiv )
				{
					jbuf[my_rank*nr+jpiv] = kmin;
					jbuf[my_rank*nr+jmin] = k;
					iw[kmin] = jpiv;
					iw[k] = jmin;
					j = levls[my_rank*nr+jpiv];
					levls[my_rank*nr+jpiv] = levls[my_rank*nr+jmin];
					levls[my_rank*nr+jmin] = j;
					k = kmin;
				}

				for(j=0;j<U->nnz[k];j++)
				{
					col = U->index[k][j];
					it = ulvl[k][j] + levls[my_rank*nr+jpiv]+1;
					if( it > levfill ) continue;
					ip = iw[col];
					if( ip==-1 )
					{
						if( col < i )
						{
							jbuf[my_rank*nr+incl] = col;
							levls[my_rank*nr+incl] = it;
							iw[col] = incl++;
						}
						else if( col > i )
						{
							jbuf[my_rank*nr+incu] = col;
							levls[my_rank*nr+incu] = it;
							iw[col] = incu++;
						}
					}
					else
					{
						levls[my_rank*nr+ip] = _min(levls[my_rank*nr+ip],it);
					}
				}
			}
			for(j=0;j<incl;j++) iw[jbuf[my_rank*nr+j]] = -1;
			for(j=i;j<incu;j++) iw[jbuf[my_rank*nr+j]] = -1;

			L->nnz[i] = incl;
			if( incl > 0 )
			{
				L->index[i] = (LIS_INT *)malloc(incl*sizeof(LIS_INT));
				L->value[i] = (LIS_SCALAR *)malloc(bs*incl*sizeof(LIS_SCALAR));
				memcpy(L->index[i],&jbuf[my_rank*nr],incl*sizeof(LIS_INT));
			}

			k = incu-i;
			U->nnz[i] = k;
			if( k > 0 )
			{
				U->index[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
				U->value[i] = (LIS_SCALAR *)malloc(bs*k*sizeof(LIS_SCALAR));
				ulvl[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
				memcpy(U->index[i],&jbuf[my_rank*nr+i],k*sizeof(LIS_INT));
				memcpy(ulvl[i],&levls[my_rank*nr+i],k*sizeof(LIS_INT));
			}
		}
	}

	precon->L  = L;
	precon->U  = U;
	precon->WD  = D;

	lis_free2(3,levls,jbuf,iw);
	for(i=0;i<nr-1;i++)
	{
		if(U->nnz[i]) free(ulvl[i]);
	}
	lis_free(ulvl);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT	err;
	LIS_INT	i,j,k,bnr,bs;
	LIS_INT	n,nr,levfill;
	LIS_INT	col,ip,it,jpiv,incl,incu,jmin,kmin;
	LIS_INT	*levls,*jbuf,*iw,**ulvl;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nr     = A->nr;
	bnr    = A->bnr;
	bs     = bnr*bnr;
	levfill = solver->options[LIS_OPTIONS_FILL];

	L      = NULL;
	U      = NULL;


	err = lis_matrix_ilu_create(nr,bnr,&L);
	if( err ) return err;
	err = lis_matrix_ilu_create(nr,bnr,&U);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(L);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(U);
	if( err ) return err;
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		return err;
	}


	ulvl   = (LIS_INT **)lis_malloc(nr*sizeof(LIS_INT *),"lis_symbolic_fact_bsr::ulvl");
	if( ulvl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	levls   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::levls");
	if( levls==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	jbuf   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::jbuf");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<nr;i++) iw[i]=-1;

	for(i=0;i<nr;i++)
	{
		incl = 0;
		incu = i;

		for(j=A->bptr[i];j<A->bptr[i+1];j++)
		{
			col = A->bindex[j];
			#ifdef USE_MPI
				if( col>=nr ) continue;
			#endif
			if( col < i )
			{
				jbuf[incl] = col;
				levls[incl] = 0;
				iw[col] = incl++;
			}
			else if( col > i )
			{
				jbuf[incu] = col;
				levls[incu] = 0;
				iw[col] = incu++;
			}
		}

		jpiv = -1;
		while( ++jpiv < incl )
		{
			k = jbuf[jpiv];
			kmin = k;
			jmin = jpiv;
			for(j=jpiv+1;j<incl;j++)
			{
				if( jbuf[j]<kmin )
				{
					kmin = jbuf[j];
					jmin = j;
				}
			}

			if( jmin!=jpiv )
			{
				jbuf[jpiv] = kmin;
				jbuf[jmin] = k;
				iw[kmin] = jpiv;
				iw[k] = jmin;
				j = levls[jpiv];
				levls[jpiv] = levls[jmin];
				levls[jmin] = j;
				k = kmin;
			}

			for(j=0;j<U->nnz[k];j++)
			{
				col = U->index[k][j];
				it = ulvl[k][j] + levls[jpiv]+1;
				if( it > levfill ) continue;
				ip = iw[col];
				if( ip==-1 )
				{
					if( col < i )
					{
						jbuf[incl] = col;
						levls[incl] = it;
						iw[col] = incl++;
					}
					else if( col > i )
					{
						jbuf[incu] = col;
						levls[incu] = it;
						iw[col] = incu++;
					}
				}
				else
				{
					levls[ip] = _min(levls[ip],it);
				}
			}
		}
		for(j=0;j<incl;j++) iw[jbuf[j]] = -1;
		for(j=i;j<incu;j++) iw[jbuf[j]] = -1;

		L->nnz[i] = incl;
		if( incl > 0 )
		{
			L->index[i] = (LIS_INT *)malloc(incl*sizeof(LIS_INT));
			L->value[i] = (LIS_SCALAR *)malloc(bs*incl*sizeof(LIS_SCALAR));
			memcpy(L->index[i],jbuf,incl*sizeof(LIS_INT));
		}

		k = incu-i;
		U->nnz[i] = k;
		if( k > 0 )
		{
			U->index[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
			U->value[i] = (LIS_SCALAR *)malloc(bs*k*sizeof(LIS_SCALAR));
			ulvl[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
			memcpy(U->index[i],jbuf+i,k*sizeof(LIS_INT));
			memcpy(ulvl[i],levls+i,k*sizeof(LIS_INT));
		}
	}

	precon->L  = L;
	precon->U  = U;
	precon->WD  = D;

	lis_free2(3,levls,jbuf,iw);
	for(i=0;i<nr-1;i++)
	{
		if(U->nnz[i]) free(ulvl[i]);
	}
	lis_free(ulvl);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}



#undef __FUNC__
#define __FUNC__ "lis_numerical_fact_bsr"
LIS_INT lis_numerical_fact_bsr(LIS_SOLVER solver, LIS_PRECON precon)
{
#ifdef _OPENMP
	LIS_INT	i,j,k,bnr,bs;
	LIS_INT	n,nr;
	LIS_INT	col,jpos,jrow;
	LIS_INT	*jw;
	LIS_INT	is,ie,my_rank,nprocs;
	LIS_SCALAR buf[16];
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nr     = A->nr;
	bnr    = A->bnr;
	bs     = bnr*bnr;
	nprocs = omp_get_max_threads();

	L = precon->L;
	U = precon->U;
	D = precon->WD;


	jw   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_numerical_fact_bsr::jw");
	if( jw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	#pragma omp parallel private(i,j,k,is,ie,my_rank,col,jpos,jrow,buf)
	{
		my_rank  = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,nr,is,ie);

		for(i=is;i<ie;i++) jw[i] = -1;


		for(i=is;i<ie;i++)
		{
			for(j=0;j<L->nnz[i];j++)
			{
				col = L->index[i][j];
				jw[col] = j;
				memset(&L->value[i][bs*j],0,bs*sizeof(LIS_SCALAR));
			}

			jw[i] = i;
			memset(&D->value[bs*i],0,bs*sizeof(LIS_SCALAR));

			for(j=0;j<U->nnz[i];j++)
			{
				col = U->index[i][j];
				jw[col] = j;
				memset(&U->value[i][bs*j],0,bs*sizeof(LIS_SCALAR));
			}

			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				col = A->bindex[j];
				#ifdef USE_MPI
					if( col>=nr ) continue;
				#endif
				if( col>=is && col<ie )
				{
					jpos = jw[col];
					if( col<i )
					{
						memcpy(&L->value[i][bs*jpos],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
					}
					else if( col==i )
					{
						memcpy(&D->value[bs*i],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
					}
					else
					{
						memcpy(&U->value[i][bs*jpos],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
					}
				}
			}

			for(j=0;j<L->nnz[i];j++)
			{
				jrow = L->index[i][j];
				lis_array_matmat(bnr,&L->value[i][bs*j],&D->value[bs*jrow],buf,LIS_INS_VALUE);
				memcpy(&L->value[i][bs*j],buf,bs*sizeof(LIS_SCALAR));

				for(k=0;k<U->nnz[jrow];k++)
				{
					col = U->index[jrow][k];
					jpos = jw[col];
					if( jpos==-1 ) continue;
					if( col<i )
					{
						lis_array_matmat(bnr,&L->value[i][bs*j],&U->value[jrow][bs*k],&L->value[i][bs*jpos],LIS_SUB_VALUE);
					}
					else if( col==i )
					{
						lis_array_matmat(bnr,&L->value[i][bs*j],&U->value[jrow][bs*k],&D->value[bs*i],LIS_SUB_VALUE);
					}
					else
					{
						lis_array_matmat(bnr,&L->value[i][bs*j],&U->value[jrow][bs*k],&U->value[i][bs*jpos],LIS_SUB_VALUE);
					}
				}
			}

			for(j=0;j<L->nnz[i];j++)
			{
				col = L->index[i][j];
				jw[col] = -1;
			}
			jw[i] = -1;
			for(j=0;j<U->nnz[i];j++)
			{
				col = U->index[i][j];
				jw[col] = -1;
			}

			if( i==nr-1 )
			{
				switch(bnr)
				{
				case 2:
					if( n%2!=0 )
					{
						D->value[4*(nr-1)+3] = 1.0;
					}
					break;
				case 3:
					if( n%3==1 )
					{
						D->value[9*(nr-1)+4] = 1.0;
						D->value[9*(nr-1)+8] = 1.0;
					}
					else if( n%3==2 )
					{
						D->value[9*(nr-1)+8] = 1.0;
					}
					break;
				}
			}
			lis_array_ge(bnr,&D->value[bs*i]);
		}
	}
	lis_free(jw);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT	i,j,k,bnr,bs;
	LIS_INT	n,nr;
	LIS_INT	col,jpos,jrow;
	LIS_INT	*jw;
	LIS_SCALAR buf[16];
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nr     = A->nr;
	bnr    = A->bnr;
	bs     = bnr*bnr;

	L = precon->L;
	U = precon->U;
	D = precon->WD;


	jw   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_numerical_fact_bsr::jw");
	if( jw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}


	for(i=0;i<nr;i++) jw[i] = -1;

	for(i=0;i<nr;i++)
	{
		for(j=0;j<L->nnz[i];j++)
		{
			col = L->index[i][j];
			jw[col] = j;
			memset(&L->value[i][bs*j],0,bs*sizeof(LIS_SCALAR));
		}

		jw[i] = i;
		memset(&D->value[bs*i],0,bs*sizeof(LIS_SCALAR));

		for(j=0;j<U->nnz[i];j++)
		{
			col = U->index[i][j];
			jw[col] = j;
			memset(&U->value[i][bs*j],0,bs*sizeof(LIS_SCALAR));
		}

		for(j=A->bptr[i];j<A->bptr[i+1];j++)
		{
			col = A->bindex[j];
			#ifdef USE_MPI
				if( col>=nr ) continue;
			#endif
			jpos = jw[col];
			if( col<i )
			{
				memcpy(&L->value[i][bs*jpos],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
			}
			else if( col==i )
			{
				memcpy(&D->value[bs*i],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
			}
			else
			{
				memcpy(&U->value[i][bs*jpos],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
			}
		}

		for(j=0;j<L->nnz[i];j++)
		{
			jrow = L->index[i][j];
			lis_array_matmat(bnr,&L->value[i][bs*j],&D->value[bs*jrow],buf,LIS_INS_VALUE);
			memcpy(&L->value[i][bs*j],buf,bs*sizeof(LIS_SCALAR));

			for(k=0;k<U->nnz[jrow];k++)
			{
				col = U->index[jrow][k];
				jpos = jw[col];
				if( jpos==-1 ) continue;
				if( col<i )
				{
					lis_array_matmat(bnr,&L->value[i][bs*j],&U->value[jrow][bs*k],&L->value[i][bs*jpos],LIS_SUB_VALUE);
				}
				else if( col==i )
				{
					lis_array_matmat(bnr,&L->value[i][bs*j],&U->value[jrow][bs*k],&D->value[bs*i],LIS_SUB_VALUE);
				}
				else
				{
					lis_array_matmat(bnr,&L->value[i][bs*j],&U->value[jrow][bs*k],&U->value[i][bs*jpos],LIS_SUB_VALUE);
				}
			}
		}

		for(j=0;j<L->nnz[i];j++)
		{
			col = L->index[i][j];
			jw[col] = -1;
		}
		jw[i] = -1;
		for(j=0;j<U->nnz[i];j++)
		{
			col = U->index[i][j];
			jw[col] = -1;
		}

		if( i==nr-1 )
		{
			switch(bnr)
			{
			case 2:
				if( n%2!=0 )
				{
					D->value[4*(nr-1)+3] = 1.0;
				}
				break;
			case 3:
				if( n%3==1 )
				{
					D->value[9*(nr-1)+4] = 1.0;
					D->value[9*(nr-1)+8] = 1.0;
				}
				else if( n%3==2 )
				{
					D->value[9*(nr-1)+8] = 1.0;
				}
				break;
			}
		}
		lis_array_ge(bnr,&D->value[bs*i]);
	}
	lis_free(jw);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}


#undef __FUNC__
#define __FUNC__ "lis_psolve_iluk_bsr"
LIS_INT lis_psolve_iluk_bsr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
#ifdef _OPENMP
	LIS_INT i,j,jj,nr,bnr,bs;
	LIS_INT is,ie,my_rank,nprocs;
	LIS_SCALAR w[3];
	LIS_SCALAR *x;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG D;
	LIS_PRECON precon;

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->WD;
	x = X->value;
	nr = solver->A->nr;
	bnr = solver->A->bnr;
	bs  = bnr*bnr;
	nprocs = omp_get_max_threads();

	lis_vector_copy(B,X);
	#pragma omp parallel private(i,j,jj,is,ie,my_rank,w)
	{
		my_rank = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,nr,is,ie);

		for(i=is;i<ie;i++)
		{
			for(j=0;j<L->nnz[i];j++)
			{
				jj     = L->index[i][j];
				lis_array_matvec(bnr,&L->value[i][bs*j],&x[bnr*jj],&x[bnr*i],LIS_SUB_VALUE);
/*				x[bnr*i+0] -= L->value[i][bs*j+0]*x[bnr*jj+0] + L->value[i][bs*j+3]*x[bnr*jj+1] + L->value[i][bs*j+6]*x[bnr*jj+2];
				x[bnr*i+1] -= L->value[i][bs*j+1]*x[bnr*jj+0] + L->value[i][bs*j+4]*x[bnr*jj+1] + L->value[i][bs*j+7]*x[bnr*jj+2];
				x[bnr*i+2] -= L->value[i][bs*j+2]*x[bnr*jj+0] + L->value[i][bs*j+5]*x[bnr*jj+1] + L->value[i][bs*j+8]*x[bnr*jj+2];*/
			}
		}
		for(i=ie-1;i>=is;i--)
		{
			for(j=0;j<U->nnz[i];j++)
			{
				jj = U->index[i][j];
				lis_array_matvec(bnr,&U->value[i][bs*j],&x[bnr*jj],&x[bnr*i],LIS_SUB_VALUE);
/*				x[bnr*i+0] -= U->value[i][bs*j+0]*x[bnr*jj+0] + U->value[i][bs*j+3]*x[bnr*jj+1] + U->value[i][bs*j+6]*x[bnr*jj+2];
				x[bnr*i+1] -= U->value[i][bs*j+1]*x[bnr*jj+0] + U->value[i][bs*j+4]*x[bnr*jj+1] + U->value[i][bs*j+7]*x[bnr*jj+2];
				x[bnr*i+2] -= U->value[i][bs*j+2]*x[bnr*jj+0] + U->value[i][bs*j+5]*x[bnr*jj+1] + U->value[i][bs*j+8]*x[bnr*jj+2];*/
			}
/*			luinv(bnr,&D->value[bs*i],&x[bnr*i],w);
			memcpy(&x[bnr*i],w,bnr*sizeof(LIS_SCALAR));
			w[0] = D->value[bs*i+0]*x[bnr*i+0] + D->value[bs*i+3]*x[bnr*i+1] + D->value[bs*i+6]*x[bnr*i+2];
			w[1] = D->value[bs*i+1]*x[bnr*i+0] + D->value[bs*i+4]*x[bnr*i+1] + D->value[bs*i+7]*x[bnr*i+2];
			w[2] = D->value[bs*i+2]*x[bnr*i+0] + D->value[bs*i+5]*x[bnr*i+1] + D->value[bs*i+8]*x[bnr*i+2];
			x[bnr*i+0] = w[0];
			x[bnr*i+1] = w[1];
			x[bnr*i+2] = w[2];*/
			lis_array_matvec(bnr,&D->value[bs*i],&x[bnr*i],w,LIS_INS_VALUE);
			memcpy(&x[bnr*i],w,bnr*sizeof(LIS_SCALAR));
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT i,j,jj,nr,bnr,bs;
	LIS_SCALAR w[9];
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG D;
	LIS_PRECON precon;

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->WD;
	b = B->value;
	x = X->value;
	nr = solver->A->nr;
	bnr = solver->A->bnr;
	bs  = bnr*bnr;

	lis_vector_copy(B,X);
	for(i=0; i<nr; i++)
	{
		for(j=0;j<L->nnz[i];j++)
		{
			jj     = L->index[i][j];
			lis_array_matvec(bnr,&L->value[i][bs*j],&x[bnr*jj],&x[bnr*i],LIS_SUB_VALUE);
		}
	}
	for(i=nr-1; i>=0; i--)
	{
		for(j=0;j<U->nnz[i];j++)
		{
			jj = U->index[i][j];
			lis_array_matvec(bnr,&U->value[i][bs*j],&x[bnr*jj],&x[bnr*i],LIS_SUB_VALUE);
		}
/*		luinv(bnr,&D->value[bs*i],&x[bnr*i],w);
		memcpy(&x[bnr*i],w,bnr*sizeof(LIS_SCALAR));*/
		lis_array_matvec(bnr,&D->value[bs*i],&x[bnr*i],w,LIS_INS_VALUE);
		memcpy(&x[bnr*i],w,bnr*sizeof(LIS_SCALAR));
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_iluk_bsr"
LIS_INT lis_psolveh_iluk_bsr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
#ifdef _OPENMP
	LIS_INT i,j,jj,nr,bnr,bs;
	LIS_INT is,ie,my_rank,nprocs;
	LIS_SCALAR w[3];
	LIS_SCALAR *x;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG D;
	LIS_PRECON precon;

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->WD;
	x = X->value;
	nr = solver->A->nr;
	bnr = solver->A->bnr;
	bs  = bnr*bnr;
	nprocs = omp_get_max_threads();

	lis_vector_copy(B,X);
	#pragma omp parallel private(i,j,jj,is,ie,my_rank,w)
	{
		my_rank = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,nr,is,ie);

		for(i=is;i<ie;i++)
		{
			for(j=0;j<L->nnz[i];j++)
			{
				jj     = L->index[i][j];
				lis_array_matvec(bnr,&L->value[i][bs*j],&x[bnr*jj],&x[bnr*i],LIS_SUB_VALUE);
/*				x[bnr*i+0] -= L->value[i][bs*j+0]*x[bnr*jj+0] + L->value[i][bs*j+3]*x[bnr*jj+1] + L->value[i][bs*j+6]*x[bnr*jj+2];
				x[bnr*i+1] -= L->value[i][bs*j+1]*x[bnr*jj+0] + L->value[i][bs*j+4]*x[bnr*jj+1] + L->value[i][bs*j+7]*x[bnr*jj+2];
				x[bnr*i+2] -= L->value[i][bs*j+2]*x[bnr*jj+0] + L->value[i][bs*j+5]*x[bnr*jj+1] + L->value[i][bs*j+8]*x[bnr*jj+2];*/
			}
		}
		for(i=ie-1;i>=is;i--)
		{
			for(j=0;j<U->nnz[i];j++)
			{
				jj = U->index[i][j];
				lis_array_matvec(bnr,&U->value[i][bs*j],&x[bnr*jj],&x[bnr*i],LIS_SUB_VALUE);
/*				x[bnr*i+0] -= U->value[i][bs*j+0]*x[bnr*jj+0] + U->value[i][bs*j+3]*x[bnr*jj+1] + U->value[i][bs*j+6]*x[bnr*jj+2];
				x[bnr*i+1] -= U->value[i][bs*j+1]*x[bnr*jj+0] + U->value[i][bs*j+4]*x[bnr*jj+1] + U->value[i][bs*j+7]*x[bnr*jj+2];
				x[bnr*i+2] -= U->value[i][bs*j+2]*x[bnr*jj+0] + U->value[i][bs*j+5]*x[bnr*jj+1] + U->value[i][bs*j+8]*x[bnr*jj+2];*/
			}
/*			luinv(bnr,&D->value[bs*i],&x[bnr*i],w);
			memcpy(&x[bnr*i],w,bnr*sizeof(LIS_SCALAR));
			w[0] = D->value[bs*i+0]*x[bnr*i+0] + D->value[bs*i+3]*x[bnr*i+1] + D->value[bs*i+6]*x[bnr*i+2];
			w[1] = D->value[bs*i+1]*x[bnr*i+0] + D->value[bs*i+4]*x[bnr*i+1] + D->value[bs*i+7]*x[bnr*i+2];
			w[2] = D->value[bs*i+2]*x[bnr*i+0] + D->value[bs*i+5]*x[bnr*i+1] + D->value[bs*i+8]*x[bnr*i+2];
			x[bnr*i+0] = w[0];
			x[bnr*i+1] = w[1];
			x[bnr*i+2] = w[2];*/
			lis_array_matvec(bnr,&D->value[bs*i],&x[bnr*i],w,LIS_INS_VALUE);
			memcpy(&x[bnr*i],w,bnr*sizeof(LIS_SCALAR));
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT i,j,jj,nr,bnr,bs;
	LIS_SCALAR w[9];
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG D;
	LIS_PRECON precon;

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->WD;
	b = B->value;
	x = X->value;
	nr = solver->A->nr;
	bnr = solver->A->bnr;
	bs  = bnr*bnr;

	lis_vector_copy(B,X);
	for(i=0; i<nr; i++)
	{
		lis_array_matvech(bnr,&D->value[bs*i],&x[bnr*i],w,LIS_INS_VALUE);
		memcpy(&x[bnr*i],w,bnr*sizeof(LIS_SCALAR));
		for(j=0;j<U->nnz[i];j++)
		{
			jj = U->index[i][j];
			lis_array_matvech(bnr,&U->value[i][bs*j],&x[bnr*i],&x[bnr*jj],LIS_SUB_VALUE);
		}
	}
	for(i=nr-1; i>=0; i--)
	{
		for(j=0;j<L->nnz[i];j++)
		{
			jj     = L->index[i][j];
			lis_array_matvech(bnr,&L->value[i][bs*j],&x[bnr*i],&x[bnr*jj],LIS_SUB_VALUE);
		}
/*		luinv(bnr,&D->value[bs*i],&x[bnr*i],w);
		memcpy(&x[bnr*i],w,bnr*sizeof(LIS_SCALAR));*/
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_symbolic_fact_vbr"
LIS_INT lis_symbolic_fact_vbr(LIS_SOLVER solver, LIS_PRECON precon)
{
#ifdef _OPENMP
	LIS_INT	err;
	LIS_INT	i,j,k,bnr;
	LIS_INT	n,nr,levfill;
	LIS_INT	col,ip,it,jpiv,incl,incu,jmin,kmin;
	LIS_INT	*levls,*jbuf,*iw,**ulvl;
	LIS_INT	is,ie,my_rank,nprocs;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nr     = A->nr;
	bnr    = A->bnr;
	levfill = solver->options[LIS_OPTIONS_FILL];
	nprocs = omp_get_max_threads();

	L      = NULL;
	U      = NULL;


	err = lis_matrix_ilu_create(nr,bnr,&L);
	if( err ) return err;
	err = lis_matrix_ilu_create(nr,bnr,&U);
	if( err ) return err;
	err = lis_matrix_ilu_setVR(L);
	if( err ) return err;
	err = lis_matrix_ilu_setVR(U);
	if( err ) return err;
	memcpy(L->bsz,A->row,(nr+1)*sizeof(LIS_INT));
	memcpy(U->bsz,A->row,(nr+1)*sizeof(LIS_INT));
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		return err;
	}


	ulvl   = (LIS_INT **)lis_malloc(nr*sizeof(LIS_INT *),"lis_symbolic_fact_bsr::ulvl");
	if( ulvl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	levls   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::levls");
	if( levls==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	jbuf   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::jbuf");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<nr;i++) iw[i]=-1;

	#pragma omp parallel private(i,j,k,is,ie,my_rank,incl,incu,col,jpiv,it,ip,kmin,jmin)
	{
		my_rank  = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,nr,is,ie);
		for(i=0;i<nr;i++)
		{
			incl = 0;
			incu = i;

			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				col = A->bindex[j];
				#ifdef USE_MPI
					if( col>=nr ) continue;
				#endif
				if( col < i )
				{
					jbuf[incl] = col;
					levls[incl] = 0;
					iw[col] = incl++;
				}
				else if( col > i )
				{
					jbuf[incu] = col;
					levls[incu] = 0;
					iw[col] = incu++;
				}
			}

			jpiv = -1;
			while( ++jpiv < incl )
			{
				k = jbuf[jpiv];
				kmin = k;
				jmin = jpiv;
				for(j=jpiv+1;j<incl;j++)
				{
					if( jbuf[j]<kmin )
					{
						kmin = jbuf[j];
						jmin = j;
					}
				}

				if( jmin!=jpiv )
				{
					jbuf[jpiv] = kmin;
					jbuf[jmin] = k;
					iw[kmin] = jpiv;
					iw[k] = jmin;
					j = levls[jpiv];
					levls[jpiv] = levls[jmin];
					levls[jmin] = j;
					k = kmin;
				}

				for(j=0;j<U->nnz[k];j++)
				{
					col = U->index[k][j];
					it = ulvl[k][j] + levls[jpiv]+1;
					if( it > levfill ) continue;
					ip = iw[col];
					if( ip==-1 )
					{
						if( col < i )
						{
							jbuf[incl] = col;
							levls[incl] = it;
							iw[col] = incl++;
						}
						else if( col > i )
						{
							jbuf[incu] = col;
							levls[incu] = it;
							iw[col] = incu++;
						}
					}
					else
					{
						levls[ip] = _min(levls[ip],it);
					}
				}
			}
			for(j=0;j<incl;j++) iw[jbuf[j]] = -1;
			for(j=i;j<incu;j++) iw[jbuf[j]] = -1;

			L->nnz[i] = incl;
			if( incl > 0 )
			{
				L->index[i] = (LIS_INT *)malloc(incl*sizeof(LIS_INT));
				L->values[i] = (LIS_SCALAR **)malloc(incl*sizeof(LIS_SCALAR *));
				memcpy(L->index[i],jbuf,incl*sizeof(LIS_INT));
				/*
				printf("i=%d L_nnz=%d\n",i,incl);
				for(j=0;j<incl;j++)
				{
					printf("(%d,%d) ",j,jbuf[j]);
				}
				printf("\n");
				*/
			}

			k = incu-i;
			U->nnz[i] = k;
			if( k > 0 )
			{
				U->index[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
				U->values[i] = (LIS_SCALAR **)malloc(k*sizeof(LIS_SCALAR *));
				ulvl[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
				memcpy(U->index[i],jbuf+i,k*sizeof(LIS_INT));
				memcpy(ulvl[i],levls+i,k*sizeof(LIS_INT));
				/*
				printf("i=%d U_nnz=%d\n",i,k);
				for(j=i;j<incu;j++)
				{
					printf("(%d,%d) ",j,jbuf[j]);
				}
				printf("\n");
				*/
			}
		}
	}

	precon->L  = L;
	precon->U  = U;
	precon->WD  = D;

	lis_free2(3,levls,jbuf,iw);
	for(i=0;i<nr-1;i++)
	{
		if(U->nnz[i]) free(ulvl[i]);
	}
	lis_free(ulvl);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT	err;
	LIS_INT	i,j,k,bnr,bs;
	LIS_INT	n,nr,levfill;
	LIS_INT	col,ip,it,jpiv,incl,incu,jmin,kmin;
	LIS_INT	*levls,*jbuf,*iw,**ulvl;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nr     = A->nr;
	bnr    = A->bnr;
	bs     = bnr*bnr;
	levfill = solver->options[LIS_OPTIONS_FILL];

	L      = NULL;
	U      = NULL;


	err = lis_matrix_ilu_create(nr,bnr,&L);
	if( err ) return err;
	err = lis_matrix_ilu_create(nr,bnr,&U);
	if( err ) return err;
	err = lis_matrix_ilu_setVR(L);
	if( err ) return err;
	err = lis_matrix_ilu_setVR(U);
	if( err ) return err;
	memcpy(L->bsz,A->row,(nr+1)*sizeof(LIS_INT));
	memcpy(U->bsz,A->row,(nr+1)*sizeof(LIS_INT));
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		return err;
	}


	ulvl   = (LIS_INT **)lis_malloc(nr*sizeof(LIS_INT *),"lis_symbolic_fact_bsr::ulvl");
	if( ulvl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	levls   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::levls");
	if( levls==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	jbuf   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::jbuf");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_symbolic_fact_bsr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<nr;i++) iw[i]=-1;

	for(i=0;i<nr;i++)
	{
		incl = 0;
		incu = i;

		for(j=A->bptr[i];j<A->bptr[i+1];j++)
		{
			col = A->bindex[j];
			#ifdef USE_MPI
				if( col>=nr ) continue;
			#endif
			if( col < i )
			{
				jbuf[incl] = col;
				levls[incl] = 0;
				iw[col] = incl++;
			}
			else if( col > i )
			{
				jbuf[incu] = col;
				levls[incu] = 0;
				iw[col] = incu++;
			}
		}

		jpiv = -1;
		while( ++jpiv < incl )
		{
			k = jbuf[jpiv];
			kmin = k;
			jmin = jpiv;
			for(j=jpiv+1;j<incl;j++)
			{
				if( jbuf[j]<kmin )
				{
					kmin = jbuf[j];
					jmin = j;
				}
			}

			if( jmin!=jpiv )
			{
				jbuf[jpiv] = kmin;
				jbuf[jmin] = k;
				iw[kmin] = jpiv;
				iw[k] = jmin;
				j = levls[jpiv];
				levls[jpiv] = levls[jmin];
				levls[jmin] = j;
				k = kmin;
			}

			for(j=0;j<U->nnz[k];j++)
			{
				col = U->index[k][j];
				it = ulvl[k][j] + levls[jpiv]+1;
				if( it > levfill ) continue;
				ip = iw[col];
				if( ip==-1 )
				{
					if( col < i )
					{
						jbuf[incl] = col;
						levls[incl] = it;
						iw[col] = incl++;
					}
					else if( col > i )
					{
						jbuf[incu] = col;
						levls[incu] = it;
						iw[col] = incu++;
					}
				}
				else
				{
					levls[ip] = _min(levls[ip],it);
				}
			}
		}
		for(j=0;j<incl;j++) iw[jbuf[j]] = -1;
		for(j=i;j<incu;j++) iw[jbuf[j]] = -1;

		L->nnz[i] = incl;
		if( incl > 0 )
		{
			L->index[i] = (LIS_INT *)malloc(incl*sizeof(LIS_INT));
			L->values[i] = (LIS_SCALAR **)malloc(incl*sizeof(LIS_SCALAR *));
			memcpy(L->index[i],jbuf,incl*sizeof(LIS_INT));
			/*
			printf("i=%d L_nnz=%d\n",i,incl);
			for(j=0;j<incl;j++)
			{
				printf("(%d,%d) ",j,jbuf[j]);
			}
			printf("\n");
			*/
		}

		k = incu-i;
		U->nnz[i] = k;
		if( k > 0 )
		{
			U->index[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
			U->values[i] = (LIS_SCALAR **)malloc(k*sizeof(LIS_SCALAR *));
			ulvl[i] = (LIS_INT *)malloc(k*sizeof(LIS_INT));
			memcpy(U->index[i],jbuf+i,k*sizeof(LIS_INT));
			memcpy(ulvl[i],levls+i,k*sizeof(LIS_INT));
			/*
			printf("i=%d U_nnz=%d\n",i,k);
			for(j=i;j<incu;j++)
			{
				printf("(%d,%d) ",j,jbuf[j]);
			}
			printf("\n");
			*/
		}
	}

	precon->L  = L;
	precon->U  = U;
	precon->WD  = D;

	lis_free2(3,levls,jbuf,iw);
	for(i=0;i<nr-1;i++)
	{
		if(U->nnz[i]) free(ulvl[i]);
	}
	lis_free(ulvl);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#define _DIM(bs,i)	(bs[i+1]-bs[i])

#undef __FUNC__
#define __FUNC__ "lis_numerical_fact_vbr"
LIS_INT lis_numerical_fact_vbr(LIS_SOLVER solver, LIS_PRECON precon)
{
#ifdef _OPENMP
	LIS_INT	i,j,k,dim,sz,mm,nn,kk;
	LIS_INT	n,nr;
	LIS_INT	col,jpos,jrow;
	LIS_INT	*jw,*bsz;
	LIS_INT	is,ie,my_rank,nprocs;
	LIS_SCALAR buf[1024];
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nr     = A->nr;
	nprocs = omp_get_max_threads();

	L = precon->L;
	U = precon->U;
	D = precon->WD;

	bsz = A->row;


	jw   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_numerical_fact_bsr::jw");
	if( jw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}


	#pragma omp parallel private(i,j,k,is,ie,my_rank,col,jpos,jrow,buf)
	{
		my_rank  = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,nr,is,ie);

		for(i=0;i<nr;i++) jw[i] = -1;

		for(i=0;i<nr;i++)
		{
			dim = _DIM(bsz,i);
			for(j=0;j<L->nnz[i];j++)
			{
				col = L->index[i][j];
				sz  = dim*_DIM(bsz,col);
				jw[col] = j;
				L->values[i][j] = (LIS_SCALAR *)malloc(sz*sizeof(LIS_SCALAR));
				memset(L->values[i][j],0,sz*sizeof(LIS_SCALAR));
			}

			jw[i] = i;
			memset(D->v_value[i],0,dim*dim*sizeof(LIS_SCALAR));

			for(j=0;j<U->nnz[i];j++)
			{
				col = U->index[i][j];
				sz  = dim*_DIM(bsz,col);
				jw[col] = j;
				U->values[i][j] = (LIS_SCALAR *)malloc(sz*sizeof(LIS_SCALAR));
				memset(U->values[i][j],0,sz*sizeof(LIS_SCALAR));
			}

			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				col = A->bindex[j];
				#ifdef USE_MPI
					if( col>=nr ) continue;
				#endif
				sz   = _DIM(bsz,col);
				jpos = jw[col];
				if( col<i )
				{
					memcpy(L->values[i][jpos],&A->value[A->ptr[j]],dim*sz*sizeof(LIS_SCALAR));
				}
				else if( col==i )
				{
					memcpy(D->v_value[i],&A->value[A->ptr[j]],dim*sz*sizeof(LIS_SCALAR));
				}
				else
				{
					memcpy(U->values[i][jpos],&A->value[A->ptr[j]],dim*sz*sizeof(LIS_SCALAR));
				}
			}

			for(j=0;j<L->nnz[i];j++)
			{
				jrow = L->index[i][j];
				mm   = dim;
				nn   = _DIM(bsz,jrow);
				lis_array_matmat_ns(mm,nn,nn,L->values[i][j],mm,D->v_value[jrow],nn,buf,mm,LIS_INS_VALUE);
				memcpy(L->values[i][j],buf,mm*nn*sizeof(LIS_SCALAR));

				for(k=0;k<U->nnz[jrow];k++)
				{
					col = U->index[jrow][k];
					jpos = jw[col];
					if( jpos==-1 ) continue;
					if( col<i )
					{
						kk = _DIM(bsz,col);
						lis_array_matmat_ns(mm,kk,nn,L->values[i][j],mm,U->values[jrow][k],nn,L->values[i][jpos],mm,LIS_SUB_VALUE);
					}
					else if( col==i )
					{
						lis_array_matmat_ns(mm,mm,nn,L->values[i][j],mm,U->values[jrow][k],nn,D->v_value[i],mm,LIS_SUB_VALUE);
					}
					else
					{
						kk = _DIM(bsz,col);
						lis_array_matmat_ns(mm,kk,nn,L->values[i][j],mm,U->values[jrow][k],nn,U->values[i][jpos],mm,LIS_SUB_VALUE);
					}
				}
			}

			for(j=0;j<L->nnz[i];j++)
			{
				col = L->index[i][j];
				jw[col] = -1;
			}
			jw[i] = -1;
			for(j=0;j<U->nnz[i];j++)
			{
				col = U->index[i][j];
				jw[col] = -1;
			}

			lis_array_ge(dim,D->v_value[i]);
		}
	}
	lis_free(jw);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT	i,j,k,dim,sz,mm,nn,kk;
	LIS_INT	n,nr;
	LIS_INT	col,jpos,jrow;
	LIS_INT	*jw,*bsz;
	LIS_SCALAR buf[1024];
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nr     = A->nr;

	L = precon->L;
	U = precon->U;
	D = precon->WD;

	bsz = A->row;


	jw   = (LIS_INT *)lis_malloc(nr*sizeof(LIS_INT),"lis_numerical_fact_bsr::jw");
	if( jw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}


	for(i=0;i<nr;i++) jw[i] = -1;

	for(i=0;i<nr;i++)
	{
		dim = _DIM(bsz,i);
		for(j=0;j<L->nnz[i];j++)
		{
			col = L->index[i][j];
			sz  = dim*_DIM(bsz,col);
			jw[col] = j;
			L->values[i][j] = (LIS_SCALAR *)malloc(sz*sizeof(LIS_SCALAR));
			memset(L->values[i][j],0,sz*sizeof(LIS_SCALAR));
		}

		jw[i] = i;
		memset(D->v_value[i],0,dim*dim*sizeof(LIS_SCALAR));

		for(j=0;j<U->nnz[i];j++)
		{
			col = U->index[i][j];
			sz  = dim*_DIM(bsz,col);
			jw[col] = j;
			U->values[i][j] = (LIS_SCALAR *)malloc(sz*sizeof(LIS_SCALAR));
			memset(U->values[i][j],0,sz*sizeof(LIS_SCALAR));
		}

		for(j=A->bptr[i];j<A->bptr[i+1];j++)
		{
			col = A->bindex[j];
			#ifdef USE_MPI
				if( col>=nr ) continue;
			#endif
			sz   = _DIM(bsz,col);
			jpos = jw[col];
			if( col<i )
			{
				memcpy(L->values[i][jpos],&A->value[A->ptr[j]],dim*sz*sizeof(LIS_SCALAR));
			}
			else if( col==i )
			{
				memcpy(D->v_value[i],&A->value[A->ptr[j]],dim*sz*sizeof(LIS_SCALAR));
			}
			else
			{
				memcpy(U->values[i][jpos],&A->value[A->ptr[j]],dim*sz*sizeof(LIS_SCALAR));
			}
		}

		for(j=0;j<L->nnz[i];j++)
		{
			jrow = L->index[i][j];
			mm   = dim;
			nn   = _DIM(bsz,jrow);
/*			printf("i=%d j=%d jrow=%d mm=%d nn=%d matmat\n",i,j,jrow,mm,nn);*/
			lis_array_matmat_ns(mm,nn,nn,L->values[i][j],mm,D->v_value[jrow],nn,buf,mm,LIS_INS_VALUE);
			memcpy(L->values[i][j],buf,mm*nn*sizeof(LIS_SCALAR));

			for(k=0;k<U->nnz[jrow];k++)
			{
				col = U->index[jrow][k];
				jpos = jw[col];
				if( jpos==-1 ) continue;
				if( col<i )
				{
					kk = _DIM(bsz,col);
					lis_array_matmat_ns(mm,kk,nn,L->values[i][j],mm,U->values[jrow][k],nn,L->values[i][jpos],mm,LIS_SUB_VALUE);
				}
				else if( col==i )
				{
					lis_array_matmat_ns(mm,mm,nn,L->values[i][j],mm,U->values[jrow][k],nn,D->v_value[i],mm,LIS_SUB_VALUE);
				}
				else
				{
					kk = _DIM(bsz,col);
					lis_array_matmat_ns(mm,kk,nn,L->values[i][j],mm,U->values[jrow][k],nn,U->values[i][jpos],mm,LIS_SUB_VALUE);
				}
			}
		}

		for(j=0;j<L->nnz[i];j++)
		{
			col = L->index[i][j];
			jw[col] = -1;
		}
		jw[i] = -1;
		for(j=0;j<U->nnz[i];j++)
		{
			col = U->index[i][j];
			jw[col] = -1;
		}

		lis_array_ge(dim,D->v_value[i]);
	}
	lis_free(jw);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_iluk_vbr"
LIS_INT lis_psolve_iluk_vbr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
#ifdef _OPENMP
	LIS_INT i,j,jj,nr,bnr,dim,sz,*bsz;
	LIS_INT is,ie,my_rank,nprocs;
	LIS_SCALAR w[1024];
	LIS_SCALAR *x;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG D;
	LIS_PRECON precon;

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->WD;
	x = X->value;
	nr = solver->A->nr;
	bsz = L->bsz;
	nprocs = omp_get_max_threads();

	lis_vector_copy(B,X);
	#pragma omp parallel private(i,j,jj,is,ie,my_rank,w)
	{
		my_rank = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,nr,is,ie);

		for(i=0; i<nr; i++)
		{
			dim = _DIM(bsz,i);
			bnr = bsz[i];
			for(j=0;j<L->nnz[i];j++)
			{
				jj     = L->index[i][j];
				sz     = _DIM(bsz,jj);
				lis_array_matvec_ns(dim,sz,L->values[i][j],dim,&x[bsz[jj]],&x[bnr],LIS_SUB_VALUE);
			}
		}
		for(i=nr-1; i>=0; i--)
		{
			dim = _DIM(bsz,i);
			bnr = bsz[i];
			for(j=0;j<U->nnz[i];j++)
			{
				jj = U->index[i][j];
				sz = _DIM(bsz,jj);
				lis_array_matvec_ns(dim,sz,U->values[i][j],dim,&x[bsz[jj]],&x[bnr],LIS_SUB_VALUE);
			}
			lis_array_matvec_ns(dim,dim,D->v_value[i],dim,&x[bnr],w,LIS_INS_VALUE);
			memcpy(&x[bnr],w,dim*sizeof(LIS_SCALAR));
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT i,j,jj,nr,bnr,dim,sz,*bsz;
	LIS_SCALAR w[1024];
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG D;
	LIS_PRECON precon;

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->WD;
	b = B->value;
	x = X->value;
	nr = solver->A->nr;
	bsz = L->bsz;

	lis_vector_copy(B,X);
	for(i=0; i<nr; i++)
	{
		dim = _DIM(bsz,i);
		bnr = bsz[i];
		for(j=0;j<L->nnz[i];j++)
		{
			jj     = L->index[i][j];
			sz     = _DIM(bsz,jj);
			lis_array_matvec_ns(dim,sz,L->values[i][j],dim,&x[bsz[jj]],&x[bnr],LIS_SUB_VALUE);
		}
	}
	for(i=nr-1; i>=0; i--)
	{
		dim = _DIM(bsz,i);
		bnr = bsz[i];
		for(j=0;j<U->nnz[i];j++)
		{
			jj = U->index[i][j];
			sz = _DIM(bsz,jj);
			lis_array_matvec_ns(dim,sz,U->values[i][j],dim,&x[bsz[jj]],&x[bnr],LIS_SUB_VALUE);
		}
		lis_array_matvec_ns(dim,dim,D->v_value[i],dim,&x[bnr],w,LIS_INS_VALUE);
		memcpy(&x[bnr],w,dim*sizeof(LIS_SCALAR));
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}
