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
 * lis_precon_create
 * lis_psolve
 * lis_psolveh
 ************************************************/

#define BSR_ENABLE 0 
#undef __FUNC__
#define __FUNC__ "lis_precon_create_iluc"
LIS_INT lis_precon_create_iluc(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	storage,block,err;
	LIS_MATRIX A,B;

	LIS_DEBUG_FUNC_IN;

#ifdef BSR_ENABLE
	storage     = solver->options[LIS_OPTIONS_STORAGE];
	block       = solver->options[LIS_OPTIONS_STORAGE_BLOCK];

	if( storage==LIS_MATRIX_BSR )
	{
		if( solver->A->matrix_type!=storage )
		{
			err = lis_matrix_convert_self(solver);
			if( err ) return err;
		}
	}
#endif

	switch( solver->A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_matrix_split(solver->A);
		if( err )
		{
			return err;
		}
		err = lis_precon_create_iluc_csr(solver,precon);
		break;
#ifdef BSR_ENABLE
	case LIS_MATRIX_BSR:
		err = lis_matrix_split(solver->A);
		if( err )
		{
			return err;
		}
		err = lis_precon_create_iluc_bsr(solver,precon);
		lis_psolve_xxx[LIS_PRECON_TYPE_ILUC]  = lis_psolve_iluc_bsr;
		break;
#endif
	default:
		A = solver->A;
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,LIS_MATRIX_CSR);
		err = lis_matrix_convert(A,B);
		if( err ) return err;
		solver->A = B;
		err = lis_matrix_split(solver->A);
		if( err )
		{
			return err;
		}
		err = lis_precon_create_iluc_csr(solver,precon);
		lis_matrix_destroy(B);
		solver->A = A;
		break;
	}

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_precon_create_iluc_csr"
LIS_INT lis_precon_create_iluc_csr(LIS_SOLVER solver, LIS_PRECON precon)
{
#ifdef _OPENMP
	LIS_INT	err;
	LIS_INT	i,j,k,l,ii,jj,kk,ll;
	LIS_INT	is,ie,my_rank,nprocs;
	LIS_INT	n,nnz,lfil,len;
	LIS_INT	cz,cw;
	LIS_INT	*iw,*iw2,*wc,*wl,*iz,*zc,*zl;
	LIS_INT	*index,*ptr;
	LIS_SCALAR val;
	LIS_REAL t,tol,toldd;
	LIS_SCALAR *z,*w,*tmp;
	LIS_SCALAR *value;
	LIS_SCALAR m;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	tol    = solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
	m      = solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
	lfil   = (LIS_INT)((double)A->nnz/(2.0*n))*m;
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

	tmp   = (LIS_SCALAR *)lis_malloc(nprocs*n*sizeof(LIS_SCALAR),"lis_precon_create_iluc_csr::tmp");
	if( tmp==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	w   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_iluc_csr::w");
	if( w==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	z   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_iluc_csr::z");
	if( z==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	iw   = (LIS_INT *)lis_malloc(nprocs*n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(nprocs*n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw2  = (LIS_INT *)lis_malloc( nprocs*(n+1)*sizeof(LIS_INT), "lis_precon_create_iluc_csr::iw2" );
	if( iw2==NULL )
	{
		LIS_SETERR_MEM(nprocs*(n+1)*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	wc   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::wc");
	if( wc==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	wl   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::wl");
	if( wl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iz   = (LIS_INT *)lis_malloc(nprocs*n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::iz");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	zc   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::zc");
	if( wc==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	zl   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::zl");
	if( wl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	nnz = A->L->nnz;
	err = lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);


	#pragma omp parallel private(i,ii,j,jj,k,kk,l,ll,is,ie,my_rank,cz,cw,toldd,t,nnz,len)
	{
		my_rank  = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		
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
		for(i=0;i<n;i++)
		{
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				jj           = A->L->index[j];
/*				if( jj>=is && jj<ie )*/
				{
					iw2[my_rank*n + jj]++;
				}
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
		for(i=0;i<n;i++)
		{
			for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
			{
				k        = A->L->index[j];
/*				if( k>=is && k<ie )*/
				{
					l        = iw2[my_rank*n+k];
					value[l] = A->L->value[j];
					index[l] = i;
					iw2[my_rank*n+k]++;
				}
			}
		}
		

		for(i=is;i<ie;i++)
		{
			D->value[i] = A->D->value[i];
			zl[i]       = -1;
			wl[i]       = -1;
			zc[i]       = 0;
			wc[i]       = 0;
		}


	for(k=is;k<ie;k++)
	{
		cz = 0;
		cw = 0;
		/* z_j = a_kj (j>=k) */
		for(j=A->U->ptr[k];j<A->U->ptr[k+1];j++)
		{
			jj       = A->U->index[j];
			if( jj<is || jj>=ie ) continue;
			#ifdef USE_MPI
				if( jj>n-1 ) continue;
			#endif
			z[jj]    = A->U->value[j];
			iz[my_rank*n+cz++] = jj;
			zc[jj]   = 1;
		}

		/* w_j = a_jk (j>k) */
		for(j=ptr[k];j<ptr[k+1];j++)
		{
			jj       = index[j];
			if( jj<is || jj>=ie ) continue;
			w[jj]    = value[j];
			iw[my_rank*n+cw++] = jj;
			wc[jj]   = 1;
		}

		/* z_j = z_j - l_ki*u_ij */
		i = wl[k];
		while( i>=is )
		{
			ii  = wc[i];
			kk  = zc[i];
			val = L->value[i][ii];
			for(j=kk;j<U->nnz[i];j++)
			{
				jj     = U->index[i][j];
				if( jj==k ) continue;
				if( zc[jj]==1 )
				{
					z[jj]   -= val*U->value[i][j];
				}
				else
				{
					z[jj]    = -val*U->value[i][j];
					iz[my_rank*n+cz++] = jj;
					zc[jj]   = 1;
				}
			}
			ii++;
			wc[i] = ii;
			if( ii<L->nnz[i] )
			{
				kk     = L->index[i][ii];
				ll     = i;
				i      = wl[ll];
				wl[ll] = wl[kk];
				wl[kk] = ll;
			}
			else
			{
				i      = wl[i];
			}
		}

		/* w_j = w_j - l_ji*u_ik */
		i = zl[k];
		while( i>=is )
		{
			ii  = zc[i];
			kk  = wc[i];
			val = U->value[i][ii];
			for(j=kk;j<L->nnz[i];j++)
			{
				jj     = L->index[i][j];
				if( wc[jj]==1 )
				{
					w[jj]   -= val*L->value[i][j];
				}
				else
				{
					w[jj]    = -val*L->value[i][j];
					iw[my_rank*n+cw++] = jj;
					wc[jj]   = 1;
				}
			}
			ii++;
			zc[i] = ii;
			if( ii<U->nnz[i] )
			{
				kk     = U->index[i][ii];
				ll     = i;
				i      = zl[ll];
				zl[ll] = zl[kk];
				zl[kk] = ll;
			}
			else
			{
				i      = zl[i];
			}
		}

		toldd       = fabs(D->value[k])*tol;
		D->value[k] = 1.0/D->value[k];
		t           = D->value[k];
		if( cz<cw )
		{
			for(j=0;j<cz;j++)
			{
				jj = iz[my_rank*n+j];
				if( wc[jj] )
				{
					D->value[jj] -= z[jj]*w[jj]*t;
				}
			}
		}
		else
		{
			for(j=0;j<cw;j++)
			{
				jj = iw[my_rank*n+j];
				if( zc[jj] )
				{
					D->value[jj] -= z[jj]*w[jj]*t;
				}
			}
		}
		/* drop U */
		nnz = 0;
		for(j=0;j<cz;j++)
		{
			jj = iz[my_rank*n+j];
			if( fabs(z[jj])>toldd )
			{
				iz[my_rank*n+nnz++] = jj;
			}
			else
			{
				z[k] = z[k] + fabs(z[jj]);
				zc[jj] = 0;
			}
		}
		len = _min(lfil,nnz);
		for(j=0;j<nnz;j++) tmp[my_rank*n+j] = fabs(z[is+j]);
		lis_sort_di(0,nnz-1,&tmp[my_rank*n],&iz[my_rank*n]);
		lis_sort_i(0,len-1,&iz[my_rank*n]);
		U->nnz[k] = len;
		if( len>0 )
		{
			U->index[k] = (LIS_INT *)malloc(len*sizeof(LIS_INT));
			U->value[k] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jj = iz[my_rank*n+j];
			U->index[k][j] = jj;
			U->value[k][j] = z[jj];
		}
		for(j=len;j<cz;j++) zc[iz[my_rank*n+j]] = 0;
		cz = nnz;

		/* drop L */
		nnz = 0;
		for(j=0;j<cw;j++)
		{
			jj = iw[my_rank*n+j];
			if( fabs(w[jj])>toldd )
			{
				iw[my_rank*n+nnz++] = jj;
			}
			else
			{
				wc[jj] = 0;
			}
		}
		len = _min(lfil,nnz);
		for(j=0;j<nnz;j++) tmp[my_rank*n+j] = fabs(w[is+j]);
		lis_sort_di(0,nnz-1,&tmp[my_rank*n],&iw[my_rank*n]);
		lis_sort_i(0,len-1,&iw[my_rank*n]);
		L->nnz[k] = len;
		if( len>0 )
		{
			L->index[k] = (LIS_INT *)malloc(len*sizeof(LIS_INT));
			L->value[k] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jj = iw[my_rank*n+j];
			L->index[k][j] = jj;
			L->value[k][j] = w[jj]*t;
		}
		for(j=len;j<cw;j++) wc[iw[my_rank*n+j]] = 0;
		cw = nnz;

		/**/
		for(j=0;j<cz;j++) zc[iz[my_rank*n+j]] = 0;
		for(j=0;j<cw;j++) wc[iw[my_rank*n+j]] = 0;
		if(U->nnz[k]>0)
		{
			jj     = iz[my_rank*n+0];
			zc[k]  = 0;
			zl[k]  = zl[jj];
			zl[jj] = k;
		}
		if(L->nnz[k]>0)
		{
			jj     = iw[my_rank*n+0];
			wc[k]  = 0;
			wl[k]  = wl[jj];
			wl[jj] = k;
		}
	}
	}

	precon->L  = L;
	precon->U  = U;
	precon->D  = D;

	lis_free2(13,tmp,w,z,iw,iw2,wc,wl,iz,zc,zl,ptr,index,value);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT	err;
	LIS_INT	i,j,k,l,ii,jj,kk,ll;
	LIS_INT	n,nnz,annz,lfil,len;
	LIS_INT	cz,cw;
	LIS_INT	*iw,*wc,*wl,*iz,*zc,*zl;
	LIS_INT	*index,*ptr;
	LIS_SCALAR gamma,val;
	LIS_REAL t,tol,toldd;
	LIS_SCALAR *z,*w,*tmp;
	LIS_SCALAR *value;
	LIS_SCALAR m;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	tol    = solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
	m      = solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
	gamma  = solver->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN];
	annz   = 10+A->nnz / A->n;
	lfil   = (LIS_INT)(((double)A->nnz/(2.0*n))*m);

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

	tmp   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_iluc_csr::tmp");
	if( tmp==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	w   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_iluc_csr::w");
	if( w==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	z   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_iluc_csr::z");
	if( z==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	iw   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	wc   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::wc");
	if( wc==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	wl   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::wl");
	if( wl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iz   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::iz");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	zc   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::zc");
	if( wc==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	zl   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::zl");
	if( wl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	lis_matrix_sort_csr(A);
	nnz = A->L->nnz;
	err = lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);

	for(i=0;i<n;i++) iw[i] = 0;
	for(i=0;i<n;i++)
	{
		for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
		{
			jj = A->L->index[j];
			iw[jj]++;
		}
	}
	ptr[0] = 0;
	for(i=0;i<n;i++)
	{
		ptr[i+1] = ptr[i] + iw[i];
		iw[i]    = ptr[i];
	}
	for(i=0;i<n;i++)
	{
		for(j=A->L->ptr[i];j<A->L->ptr[i+1];j++)
		{
			k        = A->L->index[j];
			l        = iw[k];
			value[l] = A->L->value[j];
			index[l] = i;
			iw[k]++;
		}
	}


	for(i=0;i<n;i++)
	{
		D->value[i] = gamma * A->D->value[i];
		zl[i]       = -1;
		wl[i]       = -1;
		zc[i]       = 0;
		wc[i]       = 0;
	}


	for(k=0;k<n;k++)
	{
		cz = 0;
		cw = 0;
		/* z_j = a_kj (j>=k) */
		for(j=A->U->ptr[k];j<A->U->ptr[k+1];j++)
		{
			jj       = A->U->index[j];
			#ifdef USE_MPI
				if( jj>n-1 ) continue;
			#endif
			z[jj]    = A->U->value[j];
			iz[cz++] = jj;
			zc[jj]   = 1;
		}

		/* w_j = a_jk (j>k) */
		for(j=ptr[k];j<ptr[k+1];j++)
		{
			jj       = index[j];
			w[jj]    = value[j];
			iw[cw++] = jj;
			wc[jj]   = 1;
		}

		/* z_j = z_j - l_ki*u_ij */
		i = wl[k];
		while( i>=0 )
		{
			ii  = wc[i];
			kk  = zc[i];
			val = L->value[i][ii];
			for(j=kk;j<U->nnz[i];j++)
			{
				jj     = U->index[i][j];
				if( jj==k ) continue;
				if( zc[jj]==1 )
				{
					z[jj]   -= val*U->value[i][j];
				}
				else
				{
					z[jj]    = -val*U->value[i][j];
					iz[cz++] = jj;
					zc[jj]   = 1;
				}
			}
			ii++;
			wc[i] = ii;
			if( ii<L->nnz[i] )
			{
				kk     = L->index[i][ii];
				ll     = i;
				i      = wl[ll];
				wl[ll] = wl[kk];
				wl[kk] = ll;
			}
			else
			{
				i      = wl[i];
			}
		}

		/* w_j = w_j - l_ji*u_ik */
		i = zl[k];
		while( i>=0 )
		{
			ii  = zc[i];
			kk  = wc[i];
			val = U->value[i][ii];
			for(j=kk;j<L->nnz[i];j++)
			{
				jj     = L->index[i][j];
				if( wc[jj]==1 )
				{
					w[jj]   -= val*L->value[i][j];
				}
				else
				{
					w[jj]    = -val*L->value[i][j];
					iw[cw++] = jj;
					wc[jj]   = 1;
				}
			}
			ii++;
			zc[i] = ii;
			if( ii<U->nnz[i] )
			{
				kk     = U->index[i][ii];
				ll     = i;
				i      = zl[ll];
				zl[ll] = zl[kk];
				zl[kk] = ll;
			}
			else
			{
				i      = zl[i];
			}
		}

		toldd       = fabs(D->value[k])*tol;
		D->value[k] = 1.0/D->value[k];
		t           = D->value[k];
		if( cz<cw )
		{
			for(j=0;j<cz;j++)
			{
				jj = iz[j];
				if( wc[jj] )
				{
					D->value[jj] -= z[jj]*w[jj]*t;
				}
			}
		}
		else
		{
			for(j=0;j<cw;j++)
			{
				jj = iw[j];
				if( zc[jj] )
				{
					D->value[jj] -= z[jj]*w[jj]*t;
				}
			}
		}
		/* drop U */
		nnz = 0;
		for(j=0;j<cz;j++)
		{
			jj = iz[j];
			if( fabs(z[jj])>toldd )
			{
				iz[nnz++] = jj;
			}
			else
			{
				z[k] = z[k] + fabs(z[jj]);
				zc[jj] = 0;
			}
		}
		len = _min(lfil,nnz);
		for(j=0;j<nnz;j++) tmp[j] = fabs(z[j]);
		lis_sort_di(0,nnz-1,tmp,iz);
		lis_sort_i(0,len-1,iz);
		U->nnz[k] = len;
		if( len>0 )
		{
			U->index[k] = (LIS_INT *)malloc(len*sizeof(LIS_INT));
			U->value[k] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jj = iz[j];
			U->index[k][j] = jj;
			U->value[k][j] = z[jj];
		}
		for(j=len;j<cz;j++) zc[iz[j]] = 0;
		cz = nnz;

		/* drop L */
		nnz = 0;
		for(j=0;j<cw;j++)
		{
			jj = iw[j];
			if( fabs(w[jj])>toldd )
			{
				iw[nnz++] = jj;
			}
			else
			{
				wc[jj] = 0;
			}
		}
		len = _min(lfil,nnz);
		for(j=0;j<nnz;j++) tmp[j] = fabs(w[j]);
		lis_sort_di(0,nnz-1,tmp,iw);
		lis_sort_i(0,len-1,iw);
		L->nnz[k] = len;
		if( len>0 )
		{
			L->index[k] = (LIS_INT *)malloc(len*sizeof(LIS_INT));
			L->value[k] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jj = iw[j];
			L->index[k][j] = jj;
			L->value[k][j] = w[jj]*t;
		}
		for(j=len;j<cw;j++) wc[iw[j]] = 0;
		cw = nnz;

		/**/
		for(j=0;j<cz;j++) zc[iz[j]] = 0;
		for(j=0;j<cw;j++) wc[iw[j]] = 0;
		if(U->nnz[k]>0)
		{
			jj     = iz[0];
			zc[k]  = 0;
			zl[k]  = zl[jj];
			zl[jj] = k;
		}
		if(L->nnz[k]>0)
		{
			jj     = iw[0];
			wc[k]  = 0;
			wl[k]  = wl[jj];
			wl[jj] = k;
		}
	}

	precon->L  = L;
	precon->U  = U;
	precon->D  = D;

	lis_free2(12,tmp,w,z,iw,wc,wl,iz,zc,zl,ptr,index,value);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_precon_create_iluc_bsr"
LIS_INT lis_precon_create_iluc_bsr(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_Comm comm;  
	LIS_INT err;
	LIS_INT	i,j,k,l,ii,jj,kk,ll,bnr,bs;
	LIS_INT	n,nr,nnz,bnnz,lfil,len;
	LIS_INT	cz,cw;
	LIS_INT	*iw,*wc,*wl,*iz,*zc,*zl;
	LIS_INT	*index,*ptr;
	LIS_SCALAR gamma,val;
	LIS_REAL t,tol,toldd;
	LIS_SCALAR *z,*w,*tmp,wz[9],lu[9];
	LIS_SCALAR *value;
	LIS_SCALAR m,a;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A      = solver->A;
	n      = A->n;
	nr     = A->nr;
	bnr    = A->bnr;
	bs     = bnr*bnr;
	tol    = solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
	m      = solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
	gamma  = solver->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN];
	lfil   = (LIS_INT)(((double)A->bnnz/(2.0*nr))*m);

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

	tmp   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_iluc_csr::tmp");
	if( tmp==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	w   = (LIS_SCALAR *)lis_malloc(bs*nr*sizeof(LIS_SCALAR),"lis_precon_create_iluc_csr::w");
	if( w==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	z   = (LIS_SCALAR *)lis_malloc(bs*nr*sizeof(LIS_SCALAR),"lis_precon_create_iluc_csr::z");
	if( z==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	iw   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	wc   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::wc");
	if( wc==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	wl   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::wl");
	if( wl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iz   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::iz");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	zc   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::zc");
	if( wc==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	zl   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_iluc_csr::zl");
	if( wl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	bnnz = A->L->bnnz;
	err = lis_matrix_malloc_bsr(n,bnr,bnr,bnnz,&ptr,&index,&value);

	for(i=0;i<nr;i++) iw[i] = 0;
	for(i=0;i<nr;i++)
	{
		for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
		{
			jj = A->L->bindex[j];
			iw[jj]++;
		}
	}
	ptr[0] = 0;
	for(i=0;i<nr;i++)
	{
		ptr[i+1] = ptr[i] + iw[i];
		iw[i]    = ptr[i];
	}
	for(i=0;i<nr;i++)
	{
		for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
		{
			k        = A->L->bindex[j];
			l        = iw[k];
			memcpy(&value[bs*l],&A->L->value[bs*j],bs*sizeof(LIS_SCALAR));
			index[l] = i;
			iw[k]++;
		}
	}


	switch(bnr)
	{
	case 1:
		for(i=0;i<nr;i++)
		{
			D->value[i] = gamma*A->D->value[i];
			zl[i]       = -1;
			wl[i]       = -1;
			zc[i]       = 0;
			wc[i]       = 0;
		}
		break;
	case 2:
		for(i=0;i<nr;i++)
		{
			D->value[4*i+0] = gamma*A->D->value[4*i+0];
			D->value[4*i+1] = gamma*A->D->value[4*i+1];
			D->value[4*i+2] = gamma*A->D->value[4*i+2];
			D->value[4*i+3] = gamma*A->D->value[4*i+3];
			zl[i]       = -1;
			wl[i]       = -1;
			zc[i]       = 0;
			wc[i]       = 0;
		}
		if( n%2!=0 )
		{
			D->value[4*(nr-1)+3] = 1.0;
		}
		break;
	case 3:
		for(i=0;i<nr;i++)
		{
			D->value[9*i+0] = gamma*A->D->value[9*i+0];
			D->value[9*i+1] = gamma*A->D->value[9*i+1];
			D->value[9*i+2] = gamma*A->D->value[9*i+2];
			D->value[9*i+3] = gamma*A->D->value[9*i+3];
			D->value[9*i+4] = gamma*A->D->value[9*i+4];
			D->value[9*i+5] = gamma*A->D->value[9*i+5];
			D->value[9*i+6] = gamma*A->D->value[9*i+6];
			D->value[9*i+7] = gamma*A->D->value[9*i+7];
			D->value[9*i+8] = gamma*A->D->value[9*i+8];
			zl[i]       = -1;
			wl[i]       = -1;
			zc[i]       = 0;
			wc[i]       = 0;
		}
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


	for(k=0;k<nr;k++)
	{
		cz = 0;
		cw = 0;
		/* z_j = a_kj (j>=k) */
		for(j=A->U->bptr[k];j<A->U->bptr[k+1];j++)
		{
			jj       = A->U->bindex[j];
			#ifdef USE_MPI
				if( jj>nr-1 ) continue;
			#endif
#if 0
			switch(bnr)
			{
			case 1:
				toldd       = fabs(D->value[k]);
				t           = fabs(A->U->value[j]);
				break;
			case 2:
				toldd       = fabs(D->value[4*k+0]) + fabs(D->value[4*k+1]);
				toldd      += fabs(D->value[4*k+2]) + fabs(D->value[4*k+3]);
/*				toldd       = sqrt(toldd);*/
				t       = fabs(A->U->value[4*j+0]) + fabs(A->U->value[4*j+1]);
				t      += fabs(A->U->value[4*j+2]) + fabs(A->U->value[4*j+3]);
/*				t       = sqrt(t);*/
				break;
			case 3:
				toldd       = fabs(D->value[9*k+0]) + fabs(D->value[9*k+1]) + fabs(D->value[9*k+2]);
				toldd      += fabs(D->value[9*k+3]) + fabs(D->value[9*k+4]) + fabs(D->value[9*k+5]);
				toldd      += fabs(D->value[9*k+6]) + fabs(D->value[9*k+7]) + fabs(D->value[9*k+8]);
/*				toldd       = sqrt(toldd);*/
				t       = fabs(A->U->value[9*j+0]) + fabs(A->U->value[9*j+1]) + fabs(A->U->value[9*j+2]);
				t      += fabs(A->U->value[9*j+3]) + fabs(A->U->value[9*j+4]) + fabs(A->U->value[9*j+5]);
				t      += fabs(A->U->value[9*j+6]) + fabs(A->U->value[9*j+7]) + fabs(A->U->value[9*j+8]);
/*				t       = sqrt(t);*/
				break;
			}
			if( toldd>t*1.0e-10 ) 
			{
				memcpy(&z[bs*jj],&A->U->value[bs*j],bs*sizeof(LIS_SCALAR));
				iz[cz++] = jj;
				zc[jj]   = 1;
			}
#else
			memcpy(&z[bs*jj],&A->U->value[bs*j],bs*sizeof(LIS_SCALAR));
			iz[cz++] = jj;
			zc[jj]   = 1;
#endif
		}

		/* w_j = a_jk (j>k) */
		for(j=ptr[k];j<ptr[k+1];j++)
		{
			jj       = index[j];
#if 0
			switch(bnr)
			{
			case 1:
				toldd       = fabs(D->value[k]);
				t           = fabs(value[j]);
				break;
			case 2:
				toldd       = fabs(D->value[4*k+0]) + fabs(D->value[4*k+1]);
				toldd      += fabs(D->value[4*k+2]) + fabs(D->value[4*k+3]);
/*				toldd       = sqrt(toldd);*/
				t       = fabs(value[4*j+0]) + fabs(value[4*j+1]);
				t      += fabs(value[4*j+2]) + fabs(value[4*j+3]);
/*				t       = sqrt(t);*/
				break;
			case 3:
				toldd       = fabs(D->value[9*k+0]) + fabs(D->value[9*k+1]) + fabs(D->value[9*k+2]);
				toldd      += fabs(D->value[9*k+3]) + fabs(D->value[9*k+4]) + fabs(D->value[9*k+5]);
				toldd      += fabs(D->value[9*k+6]) + fabs(D->value[9*k+7]) + fabs(D->value[9*k+8]);
/*				toldd       = sqrt(toldd);*/
				t       = fabs(value[9*j+0]) + fabs(value[9*j+1]) + fabs(value[9*j+2]);
				t      += fabs(value[9*j+3]) + fabs(value[9*j+4]) + fabs(value[9*j+5]);
				t      += fabs(value[9*j+6]) + fabs(value[9*j+7]) + fabs(value[9*j+8]);
/*				t       = sqrt(t);*/
				break;
			}
			if( toldd>t*1.0e-10 ) 
			{
				memcpy(&w[bs*jj],&value[bs*j],bs*sizeof(LIS_SCALAR));
				iw[cw++] = jj;
				wc[jj]   = 1;
			}
#else
			memcpy(&w[bs*jj],&value[bs*j],bs*sizeof(LIS_SCALAR));
			iw[cw++] = jj;
			wc[jj]   = 1;
#endif
		}

		/* z_j = z_j - l_ki*u_ij */
		i = wl[k];
		while( i>=0 )
		{
			ii  = wc[i];
			kk  = zc[i];
			val = L->value[i][ii];
			for(j=kk;j<U->nnz[i];j++)
			{
				jj     = U->index[i][j];
#if 1
				if( jj==k ) continue;
				if( zc[jj]==1 )
#else
				if( jj==k )
				{
					switch(bnr)
					{
					case 1:
						D->value[jj]   -= val*U->value[i][j];
						break;
					case 2:
						D->value[4*jj+0] -= L->value[i][4*ii+0]*U->value[i][4*j+0];
						D->value[4*jj+0] -= L->value[i][4*ii+2]*U->value[i][4*j+1];
						D->value[4*jj+1] -= L->value[i][4*ii+1]*U->value[i][4*j+0];
						D->value[4*jj+1] -= L->value[i][4*ii+3]*U->value[i][4*j+1];
						D->value[4*jj+2] -= L->value[i][4*ii+0]*U->value[i][4*j+2];
						D->value[4*jj+2] -= L->value[i][4*ii+2]*U->value[i][4*j+3];
						D->value[4*jj+3] -= L->value[i][4*ii+1]*U->value[i][4*j+2];
						D->value[4*jj+3] -= L->value[i][4*ii+3]*U->value[i][4*j+3];
						break;
					case 3:
						D->value[9*jj+0] -= L->value[i][9*ii+0]*U->value[i][9*j+0] + L->value[i][9*ii+3]*U->value[i][9*j+1] + L->value[i][9*ii+6]*U->value[i][9*j+2];
						D->value[9*jj+1] -= L->value[i][9*ii+1]*U->value[i][9*j+0] + L->value[i][9*ii+4]*U->value[i][9*j+1] + L->value[i][9*ii+7]*U->value[i][9*j+2];
						D->value[9*jj+2] -= L->value[i][9*ii+2]*U->value[i][9*j+0] + L->value[i][9*ii+5]*U->value[i][9*j+1] + L->value[i][9*ii+8]*U->value[i][9*j+2];
						D->value[9*jj+3] -= L->value[i][9*ii+0]*U->value[i][9*j+3] + L->value[i][9*ii+3]*U->value[i][9*j+4] + L->value[i][9*ii+6]*U->value[i][9*j+5];
						D->value[9*jj+4] -= L->value[i][9*ii+1]*U->value[i][9*j+3] + L->value[i][9*ii+4]*U->value[i][9*j+4] + L->value[i][9*ii+7]*U->value[i][9*j+5];
						D->value[9*jj+5] -= L->value[i][9*ii+2]*U->value[i][9*j+3] + L->value[i][9*ii+5]*U->value[i][9*j+4] + L->value[i][9*ii+8]*U->value[i][9*j+5];
						D->value[9*jj+6] -= L->value[i][9*ii+0]*U->value[i][9*j+6] + L->value[i][9*ii+3]*U->value[i][9*j+7] + L->value[i][9*ii+6]*U->value[i][9*j+8];
						D->value[9*jj+7] -= L->value[i][9*ii+1]*U->value[i][9*j+6] + L->value[i][9*ii+4]*U->value[i][9*j+7] + L->value[i][9*ii+7]*U->value[i][9*j+8];
						D->value[9*jj+8] -= L->value[i][9*ii+2]*U->value[i][9*j+6] + L->value[i][9*ii+5]*U->value[i][9*j+7] + L->value[i][9*ii+8]*U->value[i][9*j+8];
						break;
					}
				}
				else if( zc[jj]==1 )
#endif
				{
					switch(bnr)
					{
					case 1:
						z[jj]   -= val*U->value[i][j];
						break;
					case 2:
						z[4*jj+0] -= L->value[i][4*ii+0]*U->value[i][4*j+0];
						z[4*jj+0] -= L->value[i][4*ii+2]*U->value[i][4*j+1];
						z[4*jj+1] -= L->value[i][4*ii+1]*U->value[i][4*j+0];
						z[4*jj+1] -= L->value[i][4*ii+3]*U->value[i][4*j+1];
						z[4*jj+2] -= L->value[i][4*ii+0]*U->value[i][4*j+2];
						z[4*jj+2] -= L->value[i][4*ii+2]*U->value[i][4*j+3];
						z[4*jj+3] -= L->value[i][4*ii+1]*U->value[i][4*j+2];
						z[4*jj+3] -= L->value[i][4*ii+3]*U->value[i][4*j+3];
						break;
					case 3:
						z[9*jj+0] -= L->value[i][9*ii+0]*U->value[i][9*j+0] + L->value[i][9*ii+3]*U->value[i][9*j+1] + L->value[i][9*ii+6]*U->value[i][9*j+2];
						z[9*jj+1] -= L->value[i][9*ii+1]*U->value[i][9*j+0] + L->value[i][9*ii+4]*U->value[i][9*j+1] + L->value[i][9*ii+7]*U->value[i][9*j+2];
						z[9*jj+2] -= L->value[i][9*ii+2]*U->value[i][9*j+0] + L->value[i][9*ii+5]*U->value[i][9*j+1] + L->value[i][9*ii+8]*U->value[i][9*j+2];
						z[9*jj+3] -= L->value[i][9*ii+0]*U->value[i][9*j+3] + L->value[i][9*ii+3]*U->value[i][9*j+4] + L->value[i][9*ii+6]*U->value[i][9*j+5];
						z[9*jj+4] -= L->value[i][9*ii+1]*U->value[i][9*j+3] + L->value[i][9*ii+4]*U->value[i][9*j+4] + L->value[i][9*ii+7]*U->value[i][9*j+5];
						z[9*jj+5] -= L->value[i][9*ii+2]*U->value[i][9*j+3] + L->value[i][9*ii+5]*U->value[i][9*j+4] + L->value[i][9*ii+8]*U->value[i][9*j+5];
						z[9*jj+6] -= L->value[i][9*ii+0]*U->value[i][9*j+6] + L->value[i][9*ii+3]*U->value[i][9*j+7] + L->value[i][9*ii+6]*U->value[i][9*j+8];
						z[9*jj+7] -= L->value[i][9*ii+1]*U->value[i][9*j+6] + L->value[i][9*ii+4]*U->value[i][9*j+7] + L->value[i][9*ii+7]*U->value[i][9*j+8];
						z[9*jj+8] -= L->value[i][9*ii+2]*U->value[i][9*j+6] + L->value[i][9*ii+5]*U->value[i][9*j+7] + L->value[i][9*ii+8]*U->value[i][9*j+8];
						break;
					}
				}
				else
				{
					switch(bnr)
					{
					case 1:
						z[jj]   = -val*U->value[i][j];
						break;
					case 2:
						z[4*jj+0] = -L->value[i][4*ii+0]*U->value[i][4*j+0] - L->value[i][4*ii+2]*U->value[i][4*j+1];
						z[4*jj+1] = -L->value[i][4*ii+1]*U->value[i][4*j+0] - L->value[i][4*ii+3]*U->value[i][4*j+1];
						z[4*jj+2] = -L->value[i][4*ii+0]*U->value[i][4*j+2] - L->value[i][4*ii+2]*U->value[i][4*j+3];
						z[4*jj+3] = -L->value[i][4*ii+1]*U->value[i][4*j+2] - L->value[i][4*ii+3]*U->value[i][4*j+3];
						break;
					case 3:
						z[9*jj+0] = -L->value[i][9*ii+0]*U->value[i][9*j+0] - L->value[i][9*ii+3]*U->value[i][9*j+1] - L->value[i][9*ii+6]*U->value[i][9*j+2];
						z[9*jj+1] = -L->value[i][9*ii+1]*U->value[i][9*j+0] - L->value[i][9*ii+4]*U->value[i][9*j+1] - L->value[i][9*ii+7]*U->value[i][9*j+2];
						z[9*jj+2] = -L->value[i][9*ii+2]*U->value[i][9*j+0] - L->value[i][9*ii+5]*U->value[i][9*j+1] - L->value[i][9*ii+8]*U->value[i][9*j+2];
						z[9*jj+3] = -L->value[i][9*ii+0]*U->value[i][9*j+3] - L->value[i][9*ii+3]*U->value[i][9*j+4] - L->value[i][9*ii+6]*U->value[i][9*j+5];
						z[9*jj+4] = -L->value[i][9*ii+1]*U->value[i][9*j+3] - L->value[i][9*ii+4]*U->value[i][9*j+4] - L->value[i][9*ii+7]*U->value[i][9*j+5];
						z[9*jj+5] = -L->value[i][9*ii+2]*U->value[i][9*j+3] - L->value[i][9*ii+5]*U->value[i][9*j+4] - L->value[i][9*ii+8]*U->value[i][9*j+5];
						z[9*jj+6] = -L->value[i][9*ii+0]*U->value[i][9*j+6] - L->value[i][9*ii+3]*U->value[i][9*j+7] - L->value[i][9*ii+6]*U->value[i][9*j+8];
						z[9*jj+7] = -L->value[i][9*ii+1]*U->value[i][9*j+6] - L->value[i][9*ii+4]*U->value[i][9*j+7] - L->value[i][9*ii+7]*U->value[i][9*j+8];
						z[9*jj+8] = -L->value[i][9*ii+2]*U->value[i][9*j+6] - L->value[i][9*ii+5]*U->value[i][9*j+7] - L->value[i][9*ii+8]*U->value[i][9*j+8];
						break;
					}
					iz[cz++] = jj;
					zc[jj]   = 1;
				}
			}
			ii++;
			wc[i] = ii;
			if( ii<L->nnz[i] )
			{
				kk     = L->index[i][ii];
				ll     = i;
				i      = wl[ll];
				wl[ll] = wl[kk];
				wl[kk] = ll;
			}
			else
			{
				i      = wl[i];
			}
		}

		/* w_j = w_j - l_ji*u_ik */
		i = zl[k];
		while( i>=0 )
		{
			ii  = zc[i];
			kk  = wc[i];
			val = U->value[i][ii];
			for(j=kk;j<L->nnz[i];j++)
			{
				jj     = L->index[i][j];
				if( wc[jj]==1 )
				{
					switch(bnr)
					{
					case 1:
						w[jj]   -= val*L->value[i][j];
						break;
					case 2:
						w[4*jj+0] -= U->value[i][4*ii+0]*L->value[i][4*j+0];
						w[4*jj+0] -= U->value[i][4*ii+2]*L->value[i][4*j+1];
						w[4*jj+1] -= U->value[i][4*ii+1]*L->value[i][4*j+0];
						w[4*jj+1] -= U->value[i][4*ii+3]*L->value[i][4*j+1];
						w[4*jj+2] -= U->value[i][4*ii+0]*L->value[i][4*j+2];
						w[4*jj+2] -= U->value[i][4*ii+2]*L->value[i][4*j+3];
						w[4*jj+3] -= U->value[i][4*ii+1]*L->value[i][4*j+2];
						w[4*jj+3] -= U->value[i][4*ii+3]*L->value[i][4*j+3];

/*						w[4*jj+0] -= L->value[i][4*j+0]*U->value[i][4*ii+0];
						w[4*jj+0] -= L->value[i][4*j+2]*U->value[i][4*ii+1];
						w[4*jj+1] -= L->value[i][4*j+1]*U->value[i][4*ii+0];
						w[4*jj+1] -= L->value[i][4*j+3]*U->value[i][4*ii+1];
						w[4*jj+2] -= L->value[i][4*j+0]*U->value[i][4*ii+2];
						w[4*jj+2] -= L->value[i][4*j+2]*U->value[i][4*ii+3];
						w[4*jj+3] -= L->value[i][4*j+1]*U->value[i][4*ii+2];
						w[4*jj+3] -= L->value[i][4*j+3]*U->value[i][4*ii+3];*/
						break;
					case 3:
						w[9*jj+0] -= L->value[i][9*j+0]*U->value[i][9*ii+0] + L->value[i][9*j+3]*U->value[i][9*ii+1] + L->value[i][9*j+6]*U->value[i][9*ii+2];
						w[9*jj+1] -= L->value[i][9*j+1]*U->value[i][9*ii+0] + L->value[i][9*j+4]*U->value[i][9*ii+1] + L->value[i][9*j+7]*U->value[i][9*ii+2];
						w[9*jj+2] -= L->value[i][9*j+2]*U->value[i][9*ii+0] + L->value[i][9*j+5]*U->value[i][9*ii+1] + L->value[i][9*j+8]*U->value[i][9*ii+2];
						w[9*jj+3] -= L->value[i][9*j+0]*U->value[i][9*ii+3] + L->value[i][9*j+3]*U->value[i][9*ii+4] + L->value[i][9*j+6]*U->value[i][9*ii+5];
						w[9*jj+4] -= L->value[i][9*j+1]*U->value[i][9*ii+3] + L->value[i][9*j+4]*U->value[i][9*ii+4] + L->value[i][9*j+7]*U->value[i][9*ii+5];
						w[9*jj+5] -= L->value[i][9*j+2]*U->value[i][9*ii+3] + L->value[i][9*j+5]*U->value[i][9*ii+4] + L->value[i][9*j+8]*U->value[i][9*ii+5];
						w[9*jj+6] -= L->value[i][9*j+0]*U->value[i][9*ii+6] + L->value[i][9*j+3]*U->value[i][9*ii+7] + L->value[i][9*j+6]*U->value[i][9*ii+8];
						w[9*jj+7] -= L->value[i][9*j+1]*U->value[i][9*ii+6] + L->value[i][9*j+4]*U->value[i][9*ii+7] + L->value[i][9*j+7]*U->value[i][9*ii+8];
						w[9*jj+8] -= L->value[i][9*j+2]*U->value[i][9*ii+6] + L->value[i][9*j+5]*U->value[i][9*ii+7] + L->value[i][9*j+8]*U->value[i][9*ii+8];

/*						w[9*jj+0] -= U->value[i][9*ii+0]*L->value[i][9*j+0] + U->value[i][9*ii+3]*L->value[i][9*j+1] + U->value[i][9*ii+6]*L->value[i][9*j+2];
						w[9*jj+1] -= U->value[i][9*ii+1]*L->value[i][9*j+0] + U->value[i][9*ii+4]*L->value[i][9*j+1] + U->value[i][9*ii+7]*L->value[i][9*j+2];
						w[9*jj+2] -= U->value[i][9*ii+2]*L->value[i][9*j+0] + U->value[i][9*ii+5]*L->value[i][9*j+1] + U->value[i][9*ii+8]*L->value[i][9*j+2];
						w[9*jj+3] -= U->value[i][9*ii+0]*L->value[i][9*j+3] + U->value[i][9*ii+3]*L->value[i][9*j+4] + U->value[i][9*ii+6]*L->value[i][9*j+5];
						w[9*jj+4] -= U->value[i][9*ii+1]*L->value[i][9*j+3] + U->value[i][9*ii+4]*L->value[i][9*j+4] + U->value[i][9*ii+7]*L->value[i][9*j+5];
						w[9*jj+5] -= U->value[i][9*ii+2]*L->value[i][9*j+3] + U->value[i][9*ii+5]*L->value[i][9*j+4] + U->value[i][9*ii+8]*L->value[i][9*j+5];
						w[9*jj+6] -= U->value[i][9*ii+0]*L->value[i][9*j+6] + U->value[i][9*ii+3]*L->value[i][9*j+7] + U->value[i][9*ii+6]*L->value[i][9*j+8];
						w[9*jj+7] -= U->value[i][9*ii+1]*L->value[i][9*j+6] + U->value[i][9*ii+4]*L->value[i][9*j+7] + U->value[i][9*ii+7]*L->value[i][9*j+8];
						w[9*jj+8] -= U->value[i][9*ii+2]*L->value[i][9*j+6] + U->value[i][9*ii+5]*L->value[i][9*j+7] + U->value[i][9*ii+8]*L->value[i][9*j+8];*/
						break;
					}
				}
				else
				{
					switch(bnr)
					{
					case 1:
						w[jj]   = -val*L->value[i][j];
						break;
					case 2:
						w[4*jj+0] = -U->value[i][4*ii+0]*L->value[i][4*j+0] - U->value[i][4*ii+2]*L->value[i][4*j+1];
						w[4*jj+1] = -U->value[i][4*ii+1]*L->value[i][4*j+0] - U->value[i][4*ii+3]*L->value[i][4*j+1];
						w[4*jj+2] = -U->value[i][4*ii+0]*L->value[i][4*j+2] - U->value[i][4*ii+2]*L->value[i][4*j+3];
						w[4*jj+3] = -U->value[i][4*ii+1]*L->value[i][4*j+2] - U->value[i][4*ii+3]*L->value[i][4*j+3];

/*						w[4*jj+0] = -L->value[i][4*j+0]*U->value[i][4*ii+0] - L->value[i][4*j+2]*U->value[i][4*ii+1];
						w[4*jj+1] = -L->value[i][4*j+1]*U->value[i][4*ii+0] - L->value[i][4*j+3]*U->value[i][4*ii+1];
						w[4*jj+2] = -L->value[i][4*j+0]*U->value[i][4*ii+2] - L->value[i][4*j+2]*U->value[i][4*ii+3];
						w[4*jj+3] = -L->value[i][4*j+1]*U->value[i][4*ii+2] - L->value[i][4*j+3]*U->value[i][4*ii+3];*/
						break;
					case 3:
						w[9*jj+0] = -L->value[i][9*j+0]*U->value[i][9*ii+0] - L->value[i][9*j+3]*U->value[i][9*ii+1] - L->value[i][9*j+6]*U->value[i][9*ii+2];
						w[9*jj+1] = -L->value[i][9*j+1]*U->value[i][9*ii+0] - L->value[i][9*j+4]*U->value[i][9*ii+1] - L->value[i][9*j+7]*U->value[i][9*ii+2];
						w[9*jj+2] = -L->value[i][9*j+2]*U->value[i][9*ii+0] - L->value[i][9*j+5]*U->value[i][9*ii+1] - L->value[i][9*j+8]*U->value[i][9*ii+2];
						w[9*jj+3] = -L->value[i][9*j+0]*U->value[i][9*ii+3] - L->value[i][9*j+3]*U->value[i][9*ii+4] - L->value[i][9*j+6]*U->value[i][9*ii+5];
						w[9*jj+4] = -L->value[i][9*j+1]*U->value[i][9*ii+3] - L->value[i][9*j+4]*U->value[i][9*ii+4] - L->value[i][9*j+7]*U->value[i][9*ii+5];
						w[9*jj+5] = -L->value[i][9*j+2]*U->value[i][9*ii+3] - L->value[i][9*j+5]*U->value[i][9*ii+4] - L->value[i][9*j+8]*U->value[i][9*ii+5];
						w[9*jj+6] = -L->value[i][9*j+0]*U->value[i][9*ii+6] - L->value[i][9*j+3]*U->value[i][9*ii+7] - L->value[i][9*j+6]*U->value[i][9*ii+8];
						w[9*jj+7] = -L->value[i][9*j+1]*U->value[i][9*ii+6] - L->value[i][9*j+4]*U->value[i][9*ii+7] - L->value[i][9*j+7]*U->value[i][9*ii+8];
						w[9*jj+8] = -L->value[i][9*j+2]*U->value[i][9*ii+6] - L->value[i][9*j+5]*U->value[i][9*ii+7] - L->value[i][9*j+8]*U->value[i][9*ii+8];

/*						w[9*jj+0] = -U->value[i][9*ii+0]*L->value[i][9*j+0] - U->value[i][9*ii+3]*L->value[i][9*j+1] - U->value[i][9*ii+6]*L->value[i][9*j+2];
						w[9*jj+1] = -U->value[i][9*ii+1]*L->value[i][9*j+0] - U->value[i][9*ii+4]*L->value[i][9*j+1] - U->value[i][9*ii+7]*L->value[i][9*j+2];
						w[9*jj+2] = -U->value[i][9*ii+2]*L->value[i][9*j+0] - U->value[i][9*ii+5]*L->value[i][9*j+1] - U->value[i][9*ii+8]*L->value[i][9*j+2];
						w[9*jj+3] = -U->value[i][9*ii+0]*L->value[i][9*j+3] - U->value[i][9*ii+3]*L->value[i][9*j+4] - U->value[i][9*ii+6]*L->value[i][9*j+5];
						w[9*jj+4] = -U->value[i][9*ii+1]*L->value[i][9*j+3] - U->value[i][9*ii+4]*L->value[i][9*j+4] - U->value[i][9*ii+7]*L->value[i][9*j+5];
						w[9*jj+5] = -U->value[i][9*ii+2]*L->value[i][9*j+3] - U->value[i][9*ii+5]*L->value[i][9*j+4] - U->value[i][9*ii+8]*L->value[i][9*j+5];
						w[9*jj+6] = -U->value[i][9*ii+0]*L->value[i][9*j+6] - U->value[i][9*ii+3]*L->value[i][9*j+7] - U->value[i][9*ii+6]*L->value[i][9*j+8];
						w[9*jj+7] = -U->value[i][9*ii+1]*L->value[i][9*j+6] - U->value[i][9*ii+4]*L->value[i][9*j+7] - U->value[i][9*ii+7]*L->value[i][9*j+8];
						w[9*jj+8] = -U->value[i][9*ii+2]*L->value[i][9*j+6] - U->value[i][9*ii+5]*L->value[i][9*j+7] - U->value[i][9*ii+8]*L->value[i][9*j+8];*/
						break;
					}
					iw[cw++] = jj;
					wc[jj]   = 1;
				}
			}
			ii++;
			zc[i] = ii;
			if( ii<U->nnz[i] )
			{
				kk     = U->index[i][ii];
				ll     = i;
				i      = zl[ll];
				zl[ll] = zl[kk];
				zl[kk] = ll;
			}
			else
			{
				i      = zl[i];
			}
		}

		switch(bnr)
		{
		case 1:
			toldd       = fabs(D->value[k])*tol;
			break;
		case 2:
			toldd       = fabs(D->value[4*k+0]) + fabs(D->value[4*k+1]);
			toldd      += fabs(D->value[4*k+2]) + fabs(D->value[4*k+3]);
/*			toldd       = sqrt(toldd)*tol;*/
			toldd       = toldd*tol;
			break;
		case 3:
			toldd       = fabs(D->value[9*k+0]) + fabs(D->value[9*k+1]) + fabs(D->value[9*k+2]);
			toldd      += fabs(D->value[9*k+3]) + fabs(D->value[9*k+4]) + fabs(D->value[9*k+5]);
			toldd      += fabs(D->value[9*k+6]) + fabs(D->value[9*k+7]) + fabs(D->value[9*k+8]);
/*			toldd       = sqrt(toldd)*tol;*/
			toldd       = toldd*tol;
			break;
		}
/*		D->value[k] = 1.0/D->value[k];*/
/*		t           = D->value[k];*/

#if 1
		memcpy(lu,&D->value[bs*k],bs*sizeof(LIS_SCALAR));
		for(kk=0;kk<bnr;kk++)
		{
			lu[kk*bnr+kk] = 1.0 / lu[kk*bnr+kk];
			for(ii=kk+1;ii<bnr;ii++)
			{
				a = lu[kk*bnr+ii] * lu[kk*bnr+kk];
				for(jj=kk+1;jj<bnr;jj++)
				{
					lu[jj*bnr+ii] -= a * lu[jj*bnr+kk];
				}
				lu[kk*bnr+ii] = a;
			}
		}
		memcpy(&D->value[bs*k],lu,bs*sizeof(LIS_SCALAR));
		switch(bnr)
		{
		case 1:
		  lis_printf(comm,"k=%d toldd=%e\n",k,toldd);
		  lis_printf(comm,"k=%d %e\n",k,D->value[k]);
		  break;
		case 2:
		  lis_printf(comm,"k=%d toldd=%e\n",k,toldd);
		  lis_printf(comm,"k=%d %e %e\n",k,D->value[4*k+0],D->value[4*k+2]);
		  lis_printf(comm,"k=%d %e %e\n",k,D->value[4*k+1],D->value[4*k+3]);
		  break;
		case 3:
		  lis_printf(comm,"k=%d toldd=%e\n",k,toldd);
		  lis_printf(comm,"k=%d %e %e %e\n",k,D->value[9*k+0],D->value[9*k+3],D->value[9*k+6]);
		  lis_printf(comm,"k=%d %e %e %e\n",k,D->value[9*k+1],D->value[9*k+4],D->value[9*k+7]);
		  lis_printf(comm,"k=%d %e %e %e\n",k,D->value[9*k+2],D->value[9*k+5],D->value[9*k+8]);
		  break;
		}
#endif
#if 0
		switch(bnr)
		{
		case 1:
			D->value[k] = 1.0/D->value[k];
			break;
		case 2:
			/*
				D[i] = |d0 d2|
				       |d1 d3|
				*/
			t        = 1.0/(D->value[4*k]*D->value[4*k+3] - D->value[4*k+1]*D->value[4*k+2]);

			a        = D->value[4*k];
			lis_printf(comm,"k=%d t=%e\n",k,t);
			lis_printf(comm,"k=%d %e %e\n",k,D->value[4*k+0],D->value[4*k+2]);
			lis_printf(comm,"k=%d %e %e\n",k,D->value[4*k+1],D->value[4*k+3]);
			D->value[4*k]   = t*D->value[4*k+3];
			D->value[4*k+1] = -t*D->value[4*k+1];
			D->value[4*k+2] = -t*D->value[4*k+2];
			D->value[4*k+3] = t*a;
			break;
		case 3:
			/*
				D[i] = |d0 d3 d6|
				    |d1 d4 d7|
					|d2 d5 d8|
				*/
			t       = D->value[9*k]*D->value[9*k+4]*D->value[9*k+8] + D->value[9*k+1]*D->value[9*k+5]*D->value[9*k+6] + D->value[9*k+2]*D->value[9*k+3]*D->value[9*k+7];
			t      -= D->value[9*k]*D->value[9*k+5]*D->value[9*k+7] + D->value[9*k+1]*D->value[9*k+3]*D->value[9*k+8] + D->value[9*k+2]*D->value[9*k+4]*D->value[9*k+6];
			lis_printf(comm,"k=%d t=%e\n",k,t);
			t       = 1.0 / t;
			b[0]    = t*(D->value[9*k+4]*D->value[9*k+8] - D->value[9*k+5]*D->value[9*k+7]);
			b[1]    = t*(D->value[9*k+5]*D->value[9*k+6] - D->value[9*k+3]*D->value[9*k+8]);
			b[2]    = t*(D->value[9*k+3]*D->value[9*k+7] - D->value[9*k+4]*D->value[9*k+6]);
			b[3]    = t*(D->value[9*k+2]*D->value[9*k+7] - D->value[9*k+1]*D->value[9*k+8]);
			b[4]    = t*(D->value[9*k+0]*D->value[9*k+8] - D->value[9*k+2]*D->value[9*k+6]);
			b[5]    = t*(D->value[9*k+1]*D->value[9*k+6] - D->value[9*k+0]*D->value[9*k+7]);
			b[6]    = t*(D->value[9*k+1]*D->value[9*k+5] - D->value[9*k+2]*D->value[9*k+4]);
			b[7]    = t*(D->value[9*k+2]*D->value[9*k+3] - D->value[9*k+0]*D->value[9*k+5]);
			b[8]    = t*(D->value[9*k+0]*D->value[9*k+4] - D->value[9*k+1]*D->value[9*k+3]);
			memcpy(&D->value[9*k],b,9*sizeof(LIS_SCALAR));
/*			lis_printf(comm,"k=%d t=%e\n",k,t);
			lis_printf(comm,"k=%d %e %e %e\n",k,b[0],b[3],b[6]);
			lis_printf(comm,"k=%d %e %e %e\n",k,b[1],b[4],b[7]);
			lis_printf(comm,"k=%d %e %e %e\n",k,b[2],b[5],b[8]);*/
			break;
		}
		memcpy(tt,&D->value[bs*k],bs*sizeof(LIS_SCALAR));
#endif
#if 1
		if( cz<cw )
		{
			for(j=0;j<cz;j++)
			{
				jj = iz[j];
				if( wc[jj] )
				{
#if 1
					for(i=0;i<bnr;i++)
					{
						wz[i] = -w[bs*jj+i]*lu[0];
						for(ii=1;ii<bnr;ii++)
						{
							a = -w[bs*jj+ii*bnr+i];
							for(kk=0;kk<ii-1;kk++)
							{
								a -= wz[kk*bnr+i] * lu[ii*bnr+kk];
							}
							wz[ii*bnr+i] = a * lu[ii*bnr+ii];
						}
					}
					for(i=0;i<bnr;i++)
					{
						for(ii=bnr-1;ii>=0;ii--)
						{
							a = wz[ii*bnr+i];
							for(kk=ii+1;kk<bnr;kk++)
							{
								a -= wz[kk*bnr+i] * lu[ii*bnr+kk];
							}
							wz[ii*bnr+i] = a;
						}
					}
#endif
					switch(bnr)
					{
					case 1:
/*						D->value[jj] -= z[jj]*w[jj];*/
/*						D->value[jj] -= z[jj]*w[jj]*tt[0];*/
						D->value[jj] -= z[jj]*wz[0];
						break;
					case 2:
/*						wz[0] = w[4*jj+0]*tt[0] + w[4*jj+2]*tt[1];
						wz[1] = w[4*jj+1]*tt[0] + w[4*jj+3]*tt[1];
						wz[2] = w[4*jj+0]*tt[2] + w[4*jj+2]*tt[3];
						wz[3] = w[4*jj+1]*tt[2] + w[4*jj+3]*tt[3];*/
						D->value[4*jj+0] -= z[4*jj+0]*wz[0] + z[4*jj+2]*wz[1];
						D->value[4*jj+1] -= z[4*jj+1]*wz[0] + z[4*jj+3]*wz[1];
						D->value[4*jj+2] -= z[4*jj+0]*wz[2] + z[4*jj+2]*wz[3];
						D->value[4*jj+3] -= z[4*jj+1]*wz[2] + z[4*jj+3]*wz[3];

/*						D->value[4*jj+0] -= w[4*jj+0]*z[4*jj+0] + w[4*jj+2]*z[4*jj+1];
						D->value[4*jj+1] -= w[4*jj+1]*z[4*jj+0] + w[4*jj+3]*z[4*jj+1];
						D->value[4*jj+2] -= w[4*jj+0]*z[4*jj+2] + w[4*jj+2]*z[4*jj+3];
						D->value[4*jj+3] -= w[4*jj+1]*z[4*jj+2] + w[4*jj+3]*z[4*jj+3];*/
						break;
					case 3:
/*						wz[0] = tt[0]*z[9*jj+0] + tt[3]*z[9*jj+1] + tt[6]*z[9*jj+2];
						wz[1] = tt[1]*z[9*jj+0] + tt[4]*z[9*jj+1] + tt[7]*z[9*jj+2];
						wz[2] = tt[2]*z[9*jj+0] + tt[5]*z[9*jj+1] + tt[8]*z[9*jj+2];
						wz[3] = tt[0]*z[9*jj+3] + tt[3]*z[9*jj+4] + tt[6]*z[9*jj+5];
						wz[4] = tt[1]*z[9*jj+3] + tt[4]*z[9*jj+4] + tt[7]*z[9*jj+5];
						wz[5] = tt[2]*z[9*jj+3] + tt[5]*z[9*jj+4] + tt[8]*z[9*jj+5];
						wz[6] = tt[0]*z[9*jj+6] + tt[3]*z[9*jj+7] + tt[6]*z[9*jj+8];
						wz[7] = tt[1]*z[9*jj+6] + tt[4]*z[9*jj+7] + tt[7]*z[9*jj+8];
						wz[8] = tt[2]*z[9*jj+6] + tt[5]*z[9*jj+7] + tt[8]*z[9*jj+8];
*/
/*						wz[0] = z[9*jj+0];
						wz[1] = z[9*jj+1];
						wz[2] = z[9*jj+2];
						wz[3] = z[9*jj+3];
						wz[4] = z[9*jj+4];
						wz[5] = z[9*jj+5];
						wz[6] = z[9*jj+6];
						wz[7] = z[9*jj+7];
						wz[8] = z[9*jj+8];*/

						D->value[9*jj+0] -= w[9*jj+0]*wz[0] + w[9*jj+3]*wz[1] + w[9*jj+6]*wz[2];
						D->value[9*jj+1] -= w[9*jj+1]*wz[0] + w[9*jj+4]*wz[1] + w[9*jj+7]*wz[2];
						D->value[9*jj+2] -= w[9*jj+2]*wz[0] + w[9*jj+5]*wz[1] + w[9*jj+8]*wz[2];
						D->value[9*jj+3] -= w[9*jj+0]*wz[3] + w[9*jj+3]*wz[4] + w[9*jj+6]*wz[5];
						D->value[9*jj+4] -= w[9*jj+1]*wz[3] + w[9*jj+4]*wz[4] + w[9*jj+7]*wz[5];
						D->value[9*jj+5] -= w[9*jj+2]*wz[3] + w[9*jj+5]*wz[4] + w[9*jj+8]*wz[5];
						D->value[9*jj+6] -= w[9*jj+0]*wz[6] + w[9*jj+3]*wz[7] + w[9*jj+6]*wz[8];
						D->value[9*jj+7] -= w[9*jj+1]*wz[6] + w[9*jj+4]*wz[7] + w[9*jj+7]*wz[8];
						D->value[9*jj+8] -= w[9*jj+2]*wz[6] + w[9*jj+5]*wz[7] + w[9*jj+8]*wz[8];
						break;
					}
				}
			}
		}
		else
		{
			for(j=0;j<cw;j++)
			{
				jj = iw[j];
				if( zc[jj] )
				{
#if 1
					for(i=0;i<bnr;i++)
					{
						wz[i] = -w[bs*jj+i]*lu[0];
						for(ii=1;ii<bnr;ii++)
						{
							a = -w[bs*jj+ii*bnr+i];
							for(kk=0;kk<ii-1;kk++)
							{
								a -= wz[kk*bnr+i] * lu[ii*bnr+kk];
							}
							wz[ii*bnr+i] = a * lu[ii*bnr+ii];
						}
					}
					for(i=0;i<bnr;i++)
					{
						for(ii=bnr-1;ii>=0;ii--)
						{
							a = wz[ii*bnr+i];
							for(kk=ii+1;kk<bnr;kk++)
							{
								a -= wz[kk*bnr+i] * lu[ii*bnr+kk];
							}
							wz[ii*bnr+i] = a;
						}
					}
#endif
					switch(bnr)
					{
					case 1:
/*						D->value[jj] -= z[jj]*w[jj];
						D->value[jj] -= z[jj]*w[jj]*tt[0];*/
						D->value[jj] -= z[jj]*wz[0];
						break;
					case 2:
/*						wz[0] = w[4*jj+0]*tt[0] + w[4*jj+2]*tt[1];
						wz[1] = w[4*jj+1]*tt[0] + w[4*jj+3]*tt[1];
						wz[2] = w[4*jj+0]*tt[2] + w[4*jj+2]*tt[3];
						wz[3] = w[4*jj+1]*tt[2] + w[4*jj+3]*tt[3];
						*/
						D->value[4*jj+0] -= z[4*jj+0]*wz[0] + z[4*jj+2]*wz[1];
						D->value[4*jj+1] -= z[4*jj+1]*wz[0] + z[4*jj+3]*wz[1];
						D->value[4*jj+2] -= z[4*jj+0]*wz[2] + z[4*jj+2]*wz[3];
						D->value[4*jj+3] -= z[4*jj+1]*wz[2] + z[4*jj+3]*wz[3];

/*						D->value[4*jj+0] -= w[4*jj+0]*z[4*jj+0] + w[4*jj+2]*z[4*jj+1];
						D->value[4*jj+1] -= w[4*jj+1]*z[4*jj+0] + w[4*jj+3]*z[4*jj+1];
						D->value[4*jj+2] -= w[4*jj+0]*z[4*jj+2] + w[4*jj+2]*z[4*jj+3];
						D->value[4*jj+3] -= w[4*jj+1]*z[4*jj+2] + w[4*jj+3]*z[4*jj+3];*/
						break;
					case 3:
/*						wz[0] = tt[0]*z[9*jj+0] + tt[3]*z[9*jj+1] + tt[6]*z[9*jj+2];
						wz[1] = tt[1]*z[9*jj+0] + tt[4]*z[9*jj+1] + tt[7]*z[9*jj+2];
						wz[2] = tt[2]*z[9*jj+0] + tt[5]*z[9*jj+1] + tt[8]*z[9*jj+2];
						wz[3] = tt[0]*z[9*jj+3] + tt[3]*z[9*jj+4] + tt[6]*z[9*jj+5];
						wz[4] = tt[1]*z[9*jj+3] + tt[4]*z[9*jj+4] + tt[7]*z[9*jj+5];
						wz[5] = tt[2]*z[9*jj+3] + tt[5]*z[9*jj+4] + tt[8]*z[9*jj+5];
						wz[6] = tt[0]*z[9*jj+6] + tt[3]*z[9*jj+7] + tt[6]*z[9*jj+8];
						wz[7] = tt[1]*z[9*jj+6] + tt[4]*z[9*jj+7] + tt[7]*z[9*jj+8];
						wz[8] = tt[2]*z[9*jj+6] + tt[5]*z[9*jj+7] + tt[8]*z[9*jj+8];
*/
/*						wz[0] = z[9*jj+0];
						wz[1] = z[9*jj+1];
						wz[2] = z[9*jj+2];
						wz[3] = z[9*jj+3];
						wz[4] = z[9*jj+4];
						wz[5] = z[9*jj+5];
						wz[6] = z[9*jj+6];
						wz[7] = z[9*jj+7];
						wz[8] = z[9*jj+8];*/

						D->value[9*jj+0] -= w[9*jj+0]*wz[0] + w[9*jj+3]*wz[1] + w[9*jj+6]*wz[2];
						D->value[9*jj+1] -= w[9*jj+1]*wz[0] + w[9*jj+4]*wz[1] + w[9*jj+7]*wz[2];
						D->value[9*jj+2] -= w[9*jj+2]*wz[0] + w[9*jj+5]*wz[1] + w[9*jj+8]*wz[2];
						D->value[9*jj+3] -= w[9*jj+0]*wz[3] + w[9*jj+3]*wz[4] + w[9*jj+6]*wz[5];
						D->value[9*jj+4] -= w[9*jj+1]*wz[3] + w[9*jj+4]*wz[4] + w[9*jj+7]*wz[5];
						D->value[9*jj+5] -= w[9*jj+2]*wz[3] + w[9*jj+5]*wz[4] + w[9*jj+8]*wz[5];
						D->value[9*jj+6] -= w[9*jj+0]*wz[6] + w[9*jj+3]*wz[7] + w[9*jj+6]*wz[8];
						D->value[9*jj+7] -= w[9*jj+1]*wz[6] + w[9*jj+4]*wz[7] + w[9*jj+7]*wz[8];
						D->value[9*jj+8] -= w[9*jj+2]*wz[6] + w[9*jj+5]*wz[7] + w[9*jj+8]*wz[8];
						break;
					}
				}
			}
		}
#endif
		/* drop U */
		nnz = 0;
		for(j=0;j<cz;j++)
		{
			jj = iz[j];
			switch(bnr)
			{
			case 1:
				t       = fabs(z[jj]);
				break;
			case 2:
				t       = fabs(z[4*jj+0]) + fabs(z[4*jj+1]);
				t      += fabs(z[4*jj+2]) + fabs(z[4*jj+3]);
/*				t       = sqrt(t);*/
				break;
			case 3:
				t       = fabs(z[9*jj+0]) + fabs(z[9*jj+1]) + fabs(z[9*jj+2]);
				t      += fabs(z[9*jj+3]) + fabs(z[9*jj+4]) + fabs(z[9*jj+5]);
				t      += fabs(z[9*jj+6]) + fabs(z[9*jj+7]) + fabs(z[9*jj+8]);
/*				t       = sqrt(t);*/
				break;
			}
			if( t>toldd )
			{
				iz[nnz++]    = jj;
				tmp[jj] = t;
			}
			else
			{
				switch(bnr)
				{
				case 1:
					z[k] = z[k] + fabs(z[jj]);
					break;
				case 2:
					z[4*k+0] += fabs(z[4*jj+0]);
					z[4*k+1] += fabs(z[4*jj+1]);
					z[4*k+1] += fabs(z[4*jj+2]);
					z[4*k+1] += fabs(z[4*jj+3]);
					break;
				case 3:
					z[9*k+0] += fabs(z[9*jj+0]);
					z[9*k+1] += fabs(z[9*jj+1]);
					z[9*k+2] += fabs(z[9*jj+2]);
					z[9*k+3] += fabs(z[9*jj+3]);
					z[9*k+4] += fabs(z[9*jj+4]);
					z[9*k+5] += fabs(z[9*jj+5]);
					z[9*k+6] += fabs(z[9*jj+6]);
					z[9*k+7] += fabs(z[9*jj+7]);
					z[9*k+8] += fabs(z[9*jj+8]);
					break;
				}
				zc[jj] = 0;
				tmp[jj] = 0.0;
			}
		}
		len = _min(lfil,nnz);
		lis_sort_di(0,nnz-1,tmp,iz);
		lis_sort_i(0,len-1,iz);
		U->nnz[k] = len;
		if( len>0 )
		{
			U->index[k] = (LIS_INT *)malloc(len*sizeof(LIS_INT));
			U->value[k] = (LIS_SCALAR *)malloc(bs*len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jj = iz[j];
			U->index[k][j] = jj;
/*			U->value[k][j] = z[jj];*/
			memcpy(&U->value[k][bs*j],&z[bs*jj],bs*sizeof(LIS_SCALAR));
/*			for(i=0;i<bnr;i++)
			{
				lis_printf(comm,"U:k=%d j=%d",k,jj);
				for(ii=0;ii<bnr;ii++)
				{
					lis_printf(comm," %e",z[bs*jj+ii*bnr+i]);
				}
				lis_printf(comm,"\n");
			}*/
		}
		for(j=len;j<cz;j++) zc[iz[j]] = 0;
		cz = nnz;
/*		lis_printf(comm,"k=%d U:lfil=%d nnz=%d len=%d\n",k,lfil,nnz,len);*/

		/* drop L */
		nnz = 0;
		for(j=0;j<cw;j++)
		{
			jj = iw[j];
			switch(bnr)
			{
			case 1:
				t       = fabs(w[jj]);
				break;
			case 2:
				t       = fabs(w[4*jj+0]) + fabs(w[4*jj+1]);
				t      += fabs(w[4*jj+2]) + fabs(w[4*jj+3]);
/*				t       = sqrt(t);*/
				break;
			case 3:
				t       = fabs(w[9*jj+0]) + fabs(w[9*jj+1]) + fabs(w[9*jj+2]);
				t      += fabs(w[9*jj+3]) + fabs(w[9*jj+4]) + fabs(w[9*jj+5]);
				t      += fabs(w[9*jj+6]) + fabs(w[9*jj+7]) + fabs(w[9*jj+8]);
/*				t       = sqrt(t);*/
				break;
			}
/*			lis_printf(comm,"k=%d L:jj=%d t=%e tol=%e\n",k,jj,t,toldd);*/
			if( t>toldd )
			{
				iw[nnz++]    = jj;
				tmp[jj] = t;
			}
			else
			{
				wc[jj] = 0;
				tmp[jj] = 0.0;
			}
		}
		len = _min(lfil,nnz);
		lis_sort_di(0,nnz-1,tmp,iw);
		lis_sort_i(0,len-1,iw);
		L->nnz[k] = len;
		if( len>0 )
		{
			L->index[k] = (LIS_INT *)malloc(len*sizeof(LIS_INT));
			L->value[k] = (LIS_SCALAR *)malloc(bs*len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jj = iw[j];
			L->index[k][j] = jj;

#if 1
			for(i=0;i<bnr;i++)
			{
				L->value[k][bs*j+i] = -w[bs*jj+i]*lu[0];
				for(ii=1;ii<bnr;ii++)
				{
					a = -w[bs*jj+ii*bnr+i];
					for(kk=0;kk<ii-1;kk++)
					{
						a -= L->value[k][bs*j+kk*bnr+i] * lu[ii*bnr+kk];
					}
					L->value[k][bs*j+ii*bnr+i] = a * lu[ii*bnr+ii];
				}
			}
			for(i=0;i<bnr;i++)
			{
				for(ii=bnr-1;ii>=0;ii--)
				{
					a = L->value[k][bs*j+ii*bnr+i];
					for(kk=ii+1;kk<bnr;kk++)
					{
						a -= L->value[k][bs*j+kk*bnr+i] * lu[ii*bnr+kk];
					}
					L->value[k][bs*j+ii*bnr+i] = a;
				}
			}
#endif
#if 0
			switch(bnr)
			{
			case 1:
				L->value[k][j] = w[jj]*tt[0];
				break;
			case 2:
				L->value[k][4*j+0] = w[4*jj+0]*tt[0] + w[4*jj+2]*tt[1];
				L->value[k][4*j+1] = w[4*jj+1]*tt[0] + w[4*jj+3]*tt[1];
				L->value[k][4*j+2] = w[4*jj+0]*tt[2] + w[4*jj+2]*tt[3];
				L->value[k][4*j+3] = w[4*jj+1]*tt[2] + w[4*jj+3]*tt[3];
				break;
			case 3:
				L->value[k][9*j+0] = w[9*jj+0]*tt[0] + w[9*jj+3]*tt[1] + w[9*jj+6]*tt[2];
				L->value[k][9*j+1] = w[9*jj+1]*tt[0] + w[9*jj+4]*tt[1] + w[9*jj+7]*tt[2];
				L->value[k][9*j+2] = w[9*jj+2]*tt[0] + w[9*jj+5]*tt[1] + w[9*jj+8]*tt[2];
				L->value[k][9*j+3] = w[9*jj+0]*tt[3] + w[9*jj+3]*tt[4] + w[9*jj+6]*tt[5];
				L->value[k][9*j+4] = w[9*jj+1]*tt[3] + w[9*jj+4]*tt[4] + w[9*jj+7]*tt[5];
				L->value[k][9*j+5] = w[9*jj+2]*tt[3] + w[9*jj+5]*tt[4] + w[9*jj+8]*tt[5];
				L->value[k][9*j+6] = w[9*jj+0]*tt[6] + w[9*jj+3]*tt[7] + w[9*jj+6]*tt[8];
				L->value[k][9*j+7] = w[9*jj+1]*tt[6] + w[9*jj+4]*tt[7] + w[9*jj+7]*tt[8];
				L->value[k][9*j+8] = w[9*jj+2]*tt[6] + w[9*jj+5]*tt[7] + w[9*jj+8]*tt[8];
				break;
			}
/*			for(i=0;i<bnr;i++)
			{
				lis_printf(comm,"L:k=%d j=%d",k,jj);
				for(ii=0;ii<bnr;ii++)
				{
					lis_printf(comm," %e",L->value[k][bs*j+ii*bnr+i]);
				}
				lis_printf(comm,"\n");
			}*/
#endif
		}
		for(j=len;j<cw;j++) wc[iw[j]] = 0;
		cw = nnz;
/*		lis_printf(comm,"k=%d L:lfil=%d nnz=%d len=%d\n",k,lfil,nnz,len);*/


		/**/
		for(j=0;j<cz;j++) zc[iz[j]] = 0;
		for(j=0;j<cw;j++) wc[iw[j]] = 0;
		if(U->nnz[k]>0)
		{
			jj     = iz[0];
			zc[k]  = 0;
			zl[k]  = zl[jj];
			zl[jj] = k;
		}
		if(L->nnz[k]>0)
		{
			jj     = iw[0];
			wc[k]  = 0;
			wl[k]  = wl[jj];
			wl[jj] = k;
		}
	}

	precon->L  = L;
	precon->U  = U;
	precon->WD  = D;

	lis_free2(12,tmp,w,z,iw,wc,wl,iz,zc,zl,ptr,index,value);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}



#undef __FUNC__
#define __FUNC__ "lis_psolve_iluc"
LIS_INT lis_psolve_iluc(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
#ifdef _OPENMP
	LIS_INT i,j,jj,n;
	LIS_INT nprocs,my_rank,is,ie;
	LIS_SCALAR w;
	LIS_SCALAR *x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_QUAD w1;
		LIS_SCALAR *xl;
	#endif


	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	n = L->n;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			nprocs = omp_get_max_threads();
			#pragma omp parallel private(i,j,jj,is,ie,my_rank,w)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

				for(i=is; i<ie; i++)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj     = L->index[i][j];
						x[jj] -= L->value[i][j] * x[i];
					}
				}
				for(i=ie-1; i>=is; i--)
				{
					w = x[i];
					for(j=0;j<U->nnz[i];j++)
					{
						jj = U->index[i][j];
						w -= U->value[i][j] * x[jj];
					}
					x[i] = w * D->value[i];
				}
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,X);
			nprocs = omp_get_max_threads();
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,w1,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,w1,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				for(i=is; i<ie; i++)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj     = L->index[i][j];
/*						if( jj>=is && jj<ie )*/
						{
							#ifndef USE_SSE2
								LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
							#else
								LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
							#endif
						}
					}
				}
				for(i=ie-1; i>=is; i--)
				{
					w1.hi = x[i];
					w1.lo = xl[i];
					for(j=0;j<U->nnz[i];j++)
					{
						jj = U->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-U->value[i][j]);
						#else
							LIS_QUAD_FMAD_SSE2(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-U->value[i][j]);
						#endif
					}
					#ifndef USE_SSE2
						LIS_QUAD_MULD(x[i],xl[i],w1.hi,w1.lo,D->value[i]);
					#else
						LIS_QUAD_MULD_SSE2(x[i],xl[i],w1.hi,w1.lo,D->value[i]);
					#endif
				}
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;

#else
	LIS_INT i,j,jj,n;
	LIS_SCALAR w;
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_QUAD w1;
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
	n = L->n;
	b = B->value;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif

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
					x[jj] -= L->value[i][j] * x[i];
				}
/*				x[i] = x[i] * L->work[i];*/
			}
			for(i=n-1; i>=0; i--)
			{
				w = x[i];
				for(j=0;j<U->nnz[i];j++)
				{
					jj = U->index[i][j];
					w -= U->value[i][j] * x[jj];
				}
				x[i] = w * D->value[i];
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,X);
			for(i=0; i<n; i++)
			{
				for(j=0;j<L->nnz[i];j++)
				{
					jj     = L->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
					#else
						LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
					#endif
					/* x[jj] -= L->value[i][j] * x[i]; */
				}
			}
			for(i=n-1; i>=0; i--)
			{
				w1.hi = x[i];
				w1.lo = xl[i];
				/* t = x[i]; */
				for(j=0;j<U->nnz[i];j++)
				{
					jj = U->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-U->value[i][j]);
					#else
						LIS_QUAD_FMAD_SSE2(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-U->value[i][j]);
					#endif
					/* t -= U->value[i][j] * x[jj]; */
				}
				#ifndef USE_SSE2
					LIS_QUAD_MULD(x[i],xl[i],w1.hi,w1.lo,D->value[i]);
				#else
					LIS_QUAD_MULD_SSE2(x[i],xl[i],w1.hi,w1.lo,D->value[i]);
				#endif
				/* x[i] = t * D->value[i]; */
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;

#endif
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_iluc_bsr"
LIS_INT lis_psolve_iluc_bsr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_INT i,j,jj,nr,bnr,ii,bs;
	LIS_SCALAR w[9],t;
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

	lis_vector_copy(B,X);
	for(i=0; i<nr; i++)
	{
		for(j=0;j<L->nnz[i];j++)
		{
			jj     = L->index[i][j];
/*			memcpy(w,&x[bnr*i],bnr*sizeof(LIS_SCALAR));*/
			switch(bnr)
			{
			case 1:
				x[jj] -= L->value[i][j] * x[i];
				break;
			case 2:
				x[2*jj+0] -= L->value[i][4*j]   * x[2*i+0];
				x[2*jj+0] -= L->value[i][4*j+2] * x[2*i+1];
				x[2*jj+1] -= L->value[i][4*j+1] * x[2*i+0];
				x[2*jj+1] -= L->value[i][4*j+3] * x[2*i+1];
				break;
			case 3:
				x[3*jj]   -= L->value[i][9*j]   * x[3*i+0] + L->value[i][9*j+3] * x[3*i+1] + L->value[i][9*j+6] * x[3*i+2];
				x[3*jj+1] -= L->value[i][9*j+1] * x[3*i+0] + L->value[i][9*j+4] * x[3*i+1] + L->value[i][9*j+7] * x[3*i+2];
				x[3*jj+2] -= L->value[i][9*j+2] * x[3*i+0] + L->value[i][9*j+5] * x[3*i+1] + L->value[i][9*j+8] * x[3*i+2];
				break;
			}
		}
	}
	for(i=nr-1; i>=0; i--)
	{
/*		w = x[i];*/
		memcpy(w,&x[bnr*i],bnr*sizeof(LIS_SCALAR));
		for(j=0;j<U->nnz[i];j++)
		{
			jj = U->index[i][j];
			switch(bnr)
			{
			case 1:
				w[0] -= U->value[i][j] * x[jj];
				break;
			case 2:
				w[0] -= U->value[i][4*j]   * x[2*jj];
				w[0] -= U->value[i][4*j+2] * x[2*jj+1];
				w[1] -= U->value[i][4*j+1] * x[2*jj];
				w[1] -= U->value[i][4*j+3] * x[2*jj+1];
				break;
			case 3:
				w[0] -= U->value[i][9*j]   * x[3*jj] + U->value[i][9*j+3] * x[3*jj+1] + U->value[i][9*j+6] * x[3*jj+2];
				w[1] -= U->value[i][9*j+1] * x[3*jj] + U->value[i][9*j+4] * x[3*jj+1] + U->value[i][9*j+7] * x[3*jj+2];
				w[2] -= U->value[i][9*j+2] * x[3*jj] + U->value[i][9*j+5] * x[3*jj+1] + U->value[i][9*j+8] * x[3*jj+2];
				break;
			}
		}
#if 0
		switch(bnr)
		{
		case 1:
			x[i] = w[0] * D->value[i];
			break;
		case 2:
			x[2*i+0] = D->value[4*i]   * w[0] + D->value[4*i+2] * w[1];
			x[2*i+1] = D->value[4*i+1] * w[0] + D->value[4*i+3] * w[1];
			break;
		case 3:
			x[3*i+0] = D->value[9*i]   * w[0] + D->value[9*i+3] * w[1] + D->value[9*i+6] * w[2];
			x[3*i+1] = D->value[9*i+1] * w[0] + D->value[9*i+4] * w[1] + D->value[9*i+7] * w[2];
			x[3*i+2] = D->value[9*i+2] * w[0] + D->value[9*i+5] * w[1] + D->value[9*i+8] * w[2];
			break;
		}
#else
		for(ii=0;ii<bnr;ii++)
		{
			t = w[ii];
			for(jj=0;jj<ii;jj++)
			{
				t -= D->value[bs*i+jj*bnr+ii] * x[bnr*i+jj];
			}
			x[bnr*i+ii] = t;
		}
		for(ii=bnr-1;ii>=0;ii--)
		{
			t = x[bnr*i+ii];
			for(jj=ii+1;jj<bnr;jj++)
			{
				t -= D->value[bs*i+jj*bnr+ii] * x[bnr*i+jj];
			}
			x[bnr*i+ii] = t * D->value[bs*i+ii*bnr+ii];
		}
#endif
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_iluc"
LIS_INT lis_psolveh_iluc(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
#ifdef _OPENMP
	LIS_INT i,j,jj,n,is,ie,my_rank,nprocs;
	LIS_SCALAR w;
	LIS_SCALAR *x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_QUAD w1;
		LIS_SCALAR *xl;
	#endif

	/*
	 *  (LU)'x = b
	 *  U'x=b
	 *  L'x=b
	 */


	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	n = L->n;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			nprocs = omp_get_max_threads();
			#pragma omp parallel private(i,j,jj,is,ie,my_rank,w)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				for(i=is; i<ie; i++)
				{
					x[i] = x[i] * D->value[i];
					for(j=0;j<U->nnz[i];j++)
					{
						jj     = U->index[i][j];
						x[jj] -= conj(U->value[i][j]) * x[i];
					}
				}
				for(i=ie-2; i>=is; i--)
				{
					w = x[i];
					for(j=0;j<L->nnz[i];j++)
					{
						jj  = L->index[i][j];
						w  -= conj(L->value[i][j]) * x[jj];
					}
					x[i] = w;
				}
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,X);
			nprocs = omp_get_max_threads();
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,w1,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,w1,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				for(i=is; i<ie; i++)
				{
					#ifndef USE_SSE2
						LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],D->value[i]);
					#else
						LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],conj(D->value[i]));
					#endif
					/* x[i] = x[i] * conj(D->value[i]); */
					for(j=0;j<U->nnz[i];j++)
					{
						jj     = U->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(U->value[i][j]));
						#else
							LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(U->value[i][j]));
						#endif
						/* x[jj] -= conj(U->value[i][j]) * x[i]; */
					}
				}
				for(i=ie-2; i>=is; i--)
				{
					w1.hi = x[i];
					w1.lo = xl[i];
					/* t = x[i]; */
					for(j=0;j<L->nnz[i];j++)
					{
						jj  = L->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-conj(L->value[i][j]));
						#else
							LIS_QUAD_FMAD_SSE2(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-conj(L->value[i][j]));
						#endif
						/* t  -= conj(L->value[i][j]) * x[jj]; */
					}
					x[i]  = w1.hi;
					xl[i] = w1.lo;
					/* x[i] = t; */
				}
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT i,j,jj,n;
	LIS_SCALAR w;
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_QUAD w1;
		LIS_SCALAR *xl;
	#endif

	/*
	 *  (LU)'x = b
	 *  U'x=b
	 *  L'x=b
	 */


	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	n = L->n;
	b = B->value;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			for(i=0; i<n; i++)
			{
				x[i] = x[i] * D->value[i];
				for(j=0;j<U->nnz[i];j++)
				{
					jj     = U->index[i][j];
					x[jj] -= conj(U->value[i][j]) * x[i];
				}
			}
			for(i=n-2; i>=0; i--)
			{
				w = x[i];
				for(j=0;j<L->nnz[i];j++)
				{
					jj  = L->index[i][j];
					w  -= conj(L->value[i][j]) * x[jj];
				}
				x[i] = w;
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,X);
			for(i=0; i<n; i++)
			{
				#ifndef USE_SSE2
					LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],D->value[i]);
				#else
					LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],conj(D->value[i]));
				#endif
				/* x[i] = x[i] * conj(D->value[i]); */
				for(j=0;j<U->nnz[i];j++)
				{
					jj     = U->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(U->value[i][j]));
					#else
						LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(U->value[i][j]));
					#endif
					/* x[jj] -= conj(U->value[i][j]) * x[i]; */
				}
			}
			for(i=n-2; i>=0; i--)
			{
				w1.hi = x[i];
				w1.lo = xl[i];
				/* t = x[i]; */
				for(j=0;j<L->nnz[i];j++)
				{
					jj  = L->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-conj(L->value[i][j]));
					#else
						LIS_QUAD_FMAD_SSE2(w1.hi,w1.lo,w1.hi,w1.lo,x[jj],xl[jj],-conj(L->value[i][j]));
					#endif
					/* t  -= conj(L->value[i][j]) * x[jj]; */
				}
				x[i]  = w1.hi;
				xl[i] = w1.lo;
				/* x[i] = t; */
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

