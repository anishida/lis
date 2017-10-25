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
 * lis_precon_create
 * lis_psolve
 * lis_psolveh
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_precon_create_is"
LIS_INT lis_precon_create_is(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	err;
	LIS_INT	k,nsol;
	LIS_MATRIX A,B;

	LIS_DEBUG_FUNC_IN;

	k    = solver->options[LIS_OPTIONS_ISLEVEL];
	nsol = solver->options[LIS_OPTIONS_SOLVER];

	if( k!=0 && (nsol<LIS_SOLVER_JACOBI || nsol>LIS_SOLVER_SOR) )
	{
		lis_psolve_xxx[LIS_PRECON_TYPE_IS]  = lis_psolve_is;
		lis_psolveh_xxx[LIS_PRECON_TYPE_IS] = lis_psolveh_is;

		if( solver->A->matrix_type!=LIS_MATRIX_CSR )
		{
			A = solver->A;
			err = lis_matrix_duplicate(A,&B);
			if( err ) return err;
			lis_matrix_set_type(B,LIS_MATRIX_CSR);
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
		}
		err = lis_matrix_split(solver->A);
		if( err )
		{
			return err;
		}
		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}
	else
	{
		lis_psolve_xxx[LIS_PRECON_TYPE_IS]  = lis_psolve_none;
		lis_psolveh_xxx[LIS_PRECON_TYPE_IS] = lis_psolveh_none;
	}

	switch( solver->A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_matrix_split(solver->A);
		if( err )
		{
			return err;
		}
		err = lis_precon_create_is_csr(solver,precon);
		break;
	default:
		A = solver->A;
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,LIS_MATRIX_CSR);
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
		err = lis_matrix_split(solver->A);
		if( err )
		{
			return err;
		}
		err = lis_precon_create_is_csr(solver,precon);
		break;
	}
/*
	err = lis_matrix_diag_duplicate(precon->A->D,&precon->A->WD);
	if( err ) return err;
	lis_matrix_diag_copy(precon->A->D,precon->A->WD);
*/
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_create_is_csr"
LIS_INT lis_precon_create_is_csr(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT i,j,k,m;
	LIS_INT n,gn,ja,jj,jb,jcol,jpos,err;
	LIS_INT nnzl,nnzu,kl,ku;
	LIS_INT *iw,*iw2;
	LIS_INT *lptr,*lindex,*uptr,*uindex;
	LIS_INT n2;
	LIS_SCALAR val,t;
	LIS_SCALAR w;
	LIS_SCALAR *lvalue,*uvalue,*diag;
	LIS_MATRIX A,P;
	LIS_VECTOR b,pb;
	LIS_Comm comm;

	LIS_DEBUG_FUNC_IN;


	A    = solver->A;
	n    = A->n;
	gn   = A->gn;
	w    = solver->params[LIS_PARAMS_ALPHA-LIS_OPTIONS_LEN]; 
	m    = solver->options[LIS_OPTIONS_M] + 1; 
	comm = A->comm;
	b    = solver->b;

	err  = lis_matrix_create(comm,&P);
	if( err ) return err;
	err  = lis_matrix_set_size(P,n,0);
	if( err ) return err;
	err  = lis_matrix_diag_mallocM(A,&diag);
	if( err ) return err;
	err  = lis_vector_duplicate(b,&pb);
	if( err ) return err;
    
	lptr  = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_precon_create_is_csr::lptr");
	uptr  = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_precon_create_is_csr::uptr");


	iw = (LIS_INT *)lis_malloc( 2*gn*sizeof(LIS_INT),"lis_precon_create_is_csr::iw" );
	memset( iw,0,gn*sizeof(LIS_INT) );
	iw2   = iw + gn;
	/*
	 *  C <- (I+S)A
	 */
	for(i=0;i<n;i++)
	{
		k     = 0;
		nnzl  = 0;
		nnzu  = 0;
		n2    = _min(A->U->ptr[i]+m,A->U->ptr[i+1]);
		/*
		 *  I*A
		 */
		for(jb=A->L->ptr[i];jb<A->L->ptr[i+1];jb++)
		{
			jcol     = A->L->index[jb];
			iw2[k++] = jcol;
			iw[jcol] = k;
			nnzl++;
		}
		for(jb=A->U->ptr[i];jb<A->U->ptr[i+1];jb++)
		{
			jcol     = A->U->index[jb];
			iw2[k++] = jcol;
			iw[jcol] = k;
			nnzu++;
		}
		/*
		 *  S*A
		 */
		for(ja=A->U->ptr[i];ja<n2;ja++)
		{
			jj = A->U->index[ja];
			#if USE_MPI
				if( jj>=n ) break;
			#endif
			for(jb=A->L->ptr[jj];jb<A->L->ptr[jj+1];jb++)
			{
				jcol = A->L->index[jb];
				jpos = iw[jcol];
				if( jpos==0 )
				{
					iw2[k++] = jcol;
					iw[jcol] = k;
					if( jcol<i )
					{
						nnzl++;
					}
					else if( jcol>i )
					{
						nnzu++;
					}
				}
			}
			jpos = iw[jj];
			if( jpos==0 )
			{
				iw2[k++] = jj;
				iw[jj]   = k;
				nnzu++;
			}
			for(jb=A->U->ptr[jj];jb<A->U->ptr[jj+1];jb++)
			{
				jcol = A->U->index[jb];
				jpos = iw[jcol];
				if( jpos==0 )
				{
					iw2[k++] = jcol;
					iw[jcol] = k;
					nnzu++;
				}
			}
		}
		lptr[i+1] = nnzl;
		uptr[i+1] = nnzu;
		for(j=0;j<k;j++)
		{
			iw[iw2[j]] = 0;
		}
	}

	lptr[0] = 0;
	uptr[0] = 0;
	for(i=0;i<n;i++)
	{
		lptr[i+1] += lptr[i]; 
		uptr[i+1] += uptr[i]; 
	}
	nnzl = lptr[n];
	nnzu = uptr[n];

	lindex  = (LIS_INT *)lis_malloc(nnzl*sizeof(LIS_INT),"lis_precon_create_is_csr::lindex");
	lvalue  = (LIS_SCALAR *)lis_malloc(nnzl*sizeof(LIS_SCALAR),"lis_precon_create_is_csr::lvalue");
	uindex  = (LIS_INT *)lis_malloc(nnzu*sizeof(LIS_INT),"lis_precon_create_is_csr::uindex");
	uvalue  = (LIS_SCALAR *)lis_malloc(nnzu*sizeof(LIS_SCALAR),"lis_precon_create_is_csr::uvalue");

	memset( iw,0,gn*sizeof(LIS_INT) );
	/*
	 *  C <- (I+S)A
	 */
	for(i=0;i<n;i++) diag[i] = A->D->value[i];
	for(i=0;i<n;i++)
	{
		kl = lptr[i];
		ku = uptr[i];
		n2 = _min(A->U->ptr[i]+m,A->U->ptr[i+1]);
		/*
		 *  I*A
		 */
		val       = A->D->value[i];
		t         = val * b->value[i];
		for(jb=A->L->ptr[i];jb<A->L->ptr[i+1];jb++)
		{
			jcol         = A->L->index[jb];
			lindex[kl]   = jcol;
			lvalue[kl++] = val * A->L->value[jb];
			iw[jcol]     = kl;
		}
		for(jb=A->U->ptr[i];jb<A->U->ptr[i+1];jb++)
		{
			jcol         = A->U->index[jb];
			uindex[ku]   = jcol;
			uvalue[ku++] = val * A->U->value[jb];
			iw[jcol]     = ku;
		}
		/*
		 *  S*A
		 */
		for(ja=A->U->ptr[i];ja<n2;ja++)
		{
			jj = A->U->index[ja];
			#if USE_MPI
				if( jj>=n ) break;
			#endif
			val  = -w * A->U->value[ja];
			t   += val * b->value[jj];
			for(jb=A->L->ptr[jj];jb<A->L->ptr[jj+1];jb++)
			{
				jcol = A->L->index[jb];
				if( jcol<i )
				{
					jpos = iw[jcol];
					if( jpos==0 )
					{
						lindex[kl]   = jcol;
						lvalue[kl++] = val * A->L->value[jb];
						iw[jcol]     = kl;
					}
					else
					{
						lvalue[jpos-1] += val * A->L->value[jb];
					}
				}
				else if( jcol>i )
				{
					jpos = iw[jcol];
					if( jpos==0 )
					{
						uindex[ku]   = jcol;
						uvalue[ku++] = val * A->L->value[jb];
						iw[jcol]     = ku;
					}
					else
					{
						uvalue[jpos-1] += val * A->L->value[jb];
					}
				}
				else
				{
					diag[i] += val * A->L->value[jb];
				}
			}
			jpos = iw[jj];
			if( jpos==0 )
			{
				uindex[ku]   = jj;
				uvalue[ku++] = val * A->D->value[jj];
				iw[jj]       = ku;
			}
			else
			{
				uvalue[jpos-1] += val * A->D->value[jj];
			}
			for(jb=A->U->ptr[jj];jb<A->U->ptr[jj+1];jb++)
			{
				jcol = A->U->index[jb];
				jpos = iw[jcol];
				if( jpos==0 )
				{
					uindex[ku]   = jcol;
					uvalue[ku++] = val * A->U->value[jb];
					iw[jcol]     = ku;
				}
				else
				{
					uvalue[jpos-1] += val * A->U->value[jb];
				}
			}
		}
		for(j=lptr[i];j<kl;j++)
		{
			iw[lindex[j]] = 0;
		}
		for(j=uptr[i];j<ku;j++)
		{
			iw[uindex[j]] = 0;
		}
		pb->value[i] = t;
	}
	lis_matrix_setDLU_csr(lptr[n],uptr[n],diag,lptr,lindex,lvalue,uptr,uindex,uvalue,P);
	lis_matrix_merge_csr(P);
	lis_matrix_assemble(P);

	precon->A  = P;
	precon->Pb = pb;
	precon->is_copy = LIS_TRUE;

	lis_free(iw);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_is"
LIS_INT lis_psolve_is(LIS_SOLVER solver, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_MATRIX A;
	LIS_INT i,j,jj,n,m;
	LIS_SCALAR t;
	LIS_SCALAR w;
	LIS_SCALAR *y,*x;

	/*
	 *  y = M^{-1}x
	 *  M^{-1}  = (I+S)
	 */

	LIS_DEBUG_FUNC_IN;

	A      = solver->A;
	n      = A->n;
	w      = solver->params[LIS_PARAMS_ALPHA-LIS_OPTIONS_LEN]; 
	m      = solver->options[LIS_OPTIONS_M] + 1; 
	y      = Y->value;
	x      = X->value;

	#ifdef USE_MPI
		LIS_MATVEC_SENDRECV;
	#endif
	#ifdef _OPENMP
	#pragma omp parallel private(i,jj,j,t)
	#endif
	{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<n;i++)
		{
			t = 0.0;
			for(j=A->U->ptr[i];j<_min(A->U->ptr[i]+m,A->U->ptr[i+1]);j++)
			{
				jj = A->U->index[j];
				t += A->U->value[j] * x[jj];
			}
			y[i] = x[i] - w*t;
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_is"
LIS_INT lis_psolveh_is(LIS_SOLVER solver, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_MATRIX A;
	LIS_INT i,j,jj,n,m,np;
	LIS_SCALAR t;
	LIS_SCALAR w;
	LIS_SCALAR *y,*x;
	#ifdef _OPENMP
		LIS_INT k,nprocs;
		LIS_SCALAR *tmp;
	#endif

	/*
	 *  y = M^{-1}x
	 *  M^{-1}  = (I+S)
	 */

	LIS_DEBUG_FUNC_IN;

	A      = solver->A;
	n      = A->n;
	np     = A->np;
	w      = solver->params[LIS_PARAMS_ALPHA-LIS_OPTIONS_LEN]; 
	m      = solver->options[LIS_OPTIONS_M] + 1; 
	y      = Y->value;
	x      = X->value;

	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		tmp = (LIS_SCALAR *)lis_malloc( nprocs*np*sizeof(LIS_SCALAR),"lis_psolveh_is::tmp" );
		#pragma omp parallel private(i,j,t,jj,k)
		{
			k = omp_get_thread_num();
			#pragma omp for
			for(j=0;j<nprocs;j++)
			{
				memset( &tmp[j*np], 0, np*sizeof(LIS_SCALAR) );
			}
			#pragma omp for 
			for(i=0; i<n; i++)
			{
				t = x[i];
				for(j=A->U->ptr[i];j<_min(A->U->ptr[i]+m,A->U->ptr[i+1]);j++)
				{
					jj  = k*np+A->U->index[j];
					tmp[jj] += w*conj(A->U->value[j]) * t;
				}
			}
			#pragma omp for 
			for(i=0;i<np;i++)
			{
				t = 0.0;
				for(j=0;j<nprocs;j++)
				{
					t += tmp[j*np+i];
				}
				y[i] = x[i] - t;
			}
		}
		lis_free(tmp);
	#else
		for(i=0; i<np; i++)
		{
			y[i] = x[i];
		}
		for(i=0; i<n; i++)
		{
			t = x[i];
			for(j=A->U->ptr[i];j<_min(A->U->ptr[i]+m,A->U->ptr[i+1]);j++)
			{
				jj  = A->U->index[j];
				y[jj] -= w*conj(A->U->value[j]) * t;
			}
		}
	#endif
	#ifdef USE_MPI
		LIS_MATVEC_REDUCE;
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
