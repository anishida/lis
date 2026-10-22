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

#define EBABLE_BSR	0

static void
lis_precon_release_ilut_numeric_values(
	LIS_MATRIX_ILU A)
{
#ifndef USE_QUAD_PRECISION

	LIS_INT i;

	if( A==NULL ||
		A->value==NULL )
	{
		return;
	}

	for(i=0;i<A->n;i++)
	{
		if( A->nnz[i]>0 &&
			A->value[i]!=NULL )
		{
			lis_free(
				A->value[i]);

			A->value[i] = NULL;
		}
	}

#else

	(void)A;

#endif
}


static LIS_INT
lis_precon_refresh_half_ilut(
	LIS_SOLVER solver,
	LIS_PRECON precon)
{
	LIS_HALF_ILU half;

	LIS_REAL underflow;
	LIS_REAL probe_error;
	LIS_REAL probe_tolerance;

	LIS_INT accepted;
	LIS_INT err;
	LIS_INT mode;

	if( solver==NULL ||
		precon==NULL )
	{
		return LIS_SUCCESS;
	}

	if( precon->half!=NULL )
	{
		lis_half_ilu_destroy(
			precon->half);

		precon->half = NULL;
	}

	mode =
		solver->precon_half_mode;

	if( mode==LIS_PRECON_HALF_OFF )
	{
		return LIS_SUCCESS;
	}

	if( precon->L==NULL ||
		precon->U==NULL ||
		precon->D==NULL ||
		precon->L->bs!=1 ||
		precon->U->bs!=1 )
	{
		return LIS_SUCCESS;
	}

	half = NULL;
	underflow = 0.0;

	err =
		lis_half_ilu_compress(
			precon->L,
			precon->U,
			&half,
			&underflow);

	/*
	 * Half storage is optional.
	 * Any compression failure keeps the normal factors.
	 */
	if( err || half==NULL )
	{
		lis_half_ilu_destroy(half);
		return LIS_SUCCESS;
	}

	accepted = LIS_FALSE;
	probe_error = 0.0;
	probe_tolerance = 0.0;

	err =
		lis_half_ilu_assess(
			precon->L,
			precon->U,
			precon->D,
			half,
			mode,
			&accepted,
			&probe_error,
			&probe_tolerance);

	if( err || !accepted )
	{
		lis_half_ilu_destroy(half);
		return LIS_SUCCESS;
	}

	precon->half = half;

	/*
	 * Non-quad builds no longer need the normal numerical
	 * L/U values after accepted binary16 storage.
	 *
	 * The symbolic pattern remains alive and is borrowed
	 * by the half representation.
	 */
	lis_precon_release_ilut_numeric_values(
		precon->L);

	lis_precon_release_ilut_numeric_values(
		precon->U);

	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_precon_create_ilut"
LIS_INT lis_precon_create_ilut(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	err;
	LIS_MATRIX A,B;

	LIS_DEBUG_FUNC_IN;

	switch( solver->A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err =
			lis_precon_create_ilut_csr(
				solver,
				precon);

		lis_psolve_xxx[LIS_PRECON_TYPE_ILUT] =
			lis_psolve_ilut_csr;

		lis_psolveh_xxx[LIS_PRECON_TYPE_ILUT] =
			lis_psolveh_ilut_csr;

		break;

	case LIS_MATRIX_USER:
		if(
			solver->A->nprocs==1
			&&
			solver->A->n==solver->A->gn)
		{
			err =
				lis_precon_create_ilut_csr(
					solver,
					precon);

			lis_psolve_xxx[LIS_PRECON_TYPE_ILUT] =
				lis_psolve_ilut_csr;

			lis_psolveh_xxx[LIS_PRECON_TYPE_ILUT] =
				lis_psolveh_ilut_csr;

			break;
		}

		/*
		 * USER callbacks expose global column indices.  Direct
		 * row access currently requires serial numbering.
		 * Distributed USER matrices retain the existing explicit
		 * conversion fallback.
		 */
		/* fall through */

	default:
		A = solver->A;
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,LIS_MATRIX_CSR);
		err = lis_matrix_convert(A,B);
		if( err ) return err;
		solver->A = B;
		err = lis_precon_create_ilut_csr(solver,precon);
		lis_psolve_xxx[LIS_PRECON_TYPE_ILUT]  = lis_psolve_ilut_csr;
		lis_psolveh_xxx[LIS_PRECON_TYPE_ILUT]  = lis_psolveh_ilut_csr;
		lis_matrix_destroy(B);
		solver->A = A;
		break;
	}

	if( err ) return err;

	err = lis_precon_refresh_half_ilut(
		solver,
		precon);

	if( err ) return err;

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}



#undef __FUNC__
#define __FUNC__ "lis_precon_create_ilut_csr"
LIS_INT lis_precon_create_ilut_csr(LIS_SOLVER solver, LIS_PRECON precon)
{
#ifdef _OPENMP
	LIS_INT	err;
	LIS_INT	i,j,k,jj;
	LIS_INT	is,ie,my_rank,nprocs;
	LIS_INT	n,lfil,len,matrix_nnz,row_err,get_err;
	LIS_SCALAR t,tol,m;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;

	LIS_REAL tnorm,tolnorm,pivot_tol,row_scale,rt;
	LIS_SCALAR fact,lxu,*wn,*w;
	LIS_INT	lenu,lenl,col,jpos,jrow,upos;
	LIS_INT	*jbuf,*iw;
	LIS_PRECON_ROW_SOURCE *rows,*row;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	tol    = solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
	pivot_tol = solver->params[
	    LIS_PARAMS_ILU_PIVOT_TOL-LIS_OPTIONS_LEN];

	err = lis_precon_check_ilu_pivot_tol(
	    A,pivot_tol);
	if( err ) return err;
	m      = solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
	err = lis_precon_row_source_count_nnz(
		A,
		&matrix_nnz);
	if( err ) return err;
	lfil   = (LIS_INT)((double)matrix_nnz/(2.0*n))*m;
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

	w   = (LIS_SCALAR *)lis_malloc(nprocs*(n+1)*sizeof(LIS_SCALAR),"lis_precon_create_ilut_csr::w");
	if( w==NULL )
	{
		LIS_SETERR_MEM(nprocs*(n+1)*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	wn = (LIS_SCALAR *)lis_malloc(nprocs*n*sizeof(LIS_SCALAR),"lis_precon_create_ilut_csr::w");
	if( wn==NULL )
	{
		LIS_SETERR_MEM(nprocs*n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	jbuf   = (LIS_INT *)lis_malloc(nprocs*n*sizeof(LIS_INT),"lis_precon_create_ilut_csr::iw");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(nprocs*n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (LIS_INT *)lis_malloc(nprocs*n*sizeof(LIS_INT),"lis_precon_create_ilut_csr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(nprocs*n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}



	rows = (LIS_PRECON_ROW_SOURCE *)lis_malloc(
		nprocs*sizeof(LIS_PRECON_ROW_SOURCE),
		"lis_precon_create_ilut_csr::rows");

	if( rows==NULL )
	{
		LIS_SETERR_MEM(
			nprocs*sizeof(LIS_PRECON_ROW_SOURCE));

		lis_free2(
			4,
			w,
			iw,
			wn,
			jbuf);

		lis_matrix_ilu_destroy(L);
		lis_matrix_ilu_destroy(U);
		lis_vector_destroy(D);

		return LIS_OUT_OF_MEMORY;
	}

	for(j=0;j<nprocs;j++)
	{
		err = lis_precon_row_source_init(
			A,
			&rows[j]);

		if( err )
		{
			while(j>0)
			{
				j--;

				lis_precon_row_source_destroy(
					&rows[j]);
			}

			lis_free(rows);

			lis_free2(
				4,
				w,
				iw,
				wn,
				jbuf);

			lis_matrix_ilu_destroy(L);
			lis_matrix_ilu_destroy(U);
			lis_vector_destroy(D);

			return err;
		}
	}

	row_err = LIS_SUCCESS;

	#pragma omp parallel private(is,ie,my_rank,i,j,k,jj,tnorm,tolnorm,len,lenu,lenl,col,t,jpos,jrow,fact,lxu,upos,get_err,row,row_scale,rt)
	{
		my_rank  = omp_get_thread_num();
		row = &rows[my_rank];
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		for(i=is;i<ie;i++) iw[my_rank*n+i] = -1;

		for(i=is;i<ie;i++)
		{
			get_err = LIS_SUCCESS;

			if( A->matrix_type==LIS_MATRIX_USER )
			{
				/*
				 * USER row callbacks are serialized until the
				 * callback API has an explicit thread-safety
				 * capability.  Row factorization remains parallel.
				 */
				#pragma omp critical(lis_precon_user_row_access)
				{
					get_err =
						lis_precon_row_source_get(
							row,
							i);

					if(
						get_err
						&&
						row_err==LIS_SUCCESS)
					{
						row_err = get_err;
					}
				}
			}
			else
			{
				get_err =
					lis_precon_row_source_get(
						row,
						i);

				if( get_err )
				{
					#pragma omp critical(lis_precon_row_error)
					{
						if( row_err==LIS_SUCCESS )
						{
							row_err = get_err;
						}
					}
				}
			}

			if( get_err )
			{
				break;
			}

			tnorm = 0;
			row_scale = 0.0;
			k = 0;

			for(j=0;j<row->nnz;j++)
			{
				jj = row->index[j];
				rt = fabs(row->value[j]);

				if( rt>row_scale )
				{
					row_scale = rt;
				}

				if( jj<is || jj>=ie )
				{
					continue;
				}

				tnorm += rt;
				k++;
			}

			tnorm   = tnorm / (double)k;
			tolnorm = tol * tnorm;

			lenu = 0;
			lenl = 0;

			jbuf[my_rank*n+i] = i;
			w[my_rank*n+i] = 0;
			iw[my_rank*n+i] = i;

			for(j=0;j<row->nnz;j++)
			{
				col = row->index[j];

				if( col<is || col>=ie )
				{
					continue;
				}

				t = row->value[j];

				if( col < i )
				{
					jbuf[my_rank*n+lenl] = col;
					iw[my_rank*n+col] = lenl;
					w[my_rank*n+lenl] = t;
					lenl++;
				}
				else if( col == i )
				{
					w[my_rank*n+i] = t;
				}
				else
				{
					lenu++;
					jpos = i + lenu;
					jbuf[my_rank*n+jpos] = col;
					iw[my_rank*n+col] = jpos;
					w[my_rank*n+jpos] = t;
				}
			}

			j = -1;
			len = 0;

			while( ++j < lenl )
			{
				jrow = jbuf[my_rank*n+j];
				jpos = j;
				for(k=j+1;k<lenl;k++)
				{
					if( jbuf[my_rank*n+k]<jrow )
					{
						jrow = jbuf[my_rank*n+k];
						jpos = k;
					}
				}
				if( jpos!=j )
				{
					col = jbuf[my_rank*n+j];
					jbuf[my_rank*n+j] = jbuf[my_rank*n+jpos];
					jbuf[my_rank*n+jpos] = col;
					iw[my_rank*n+jrow] = j;
					iw[my_rank*n+col] = jpos;
					t = w[my_rank*n+j];
					w[my_rank*n+j] = w[my_rank*n+jpos];
					w[my_rank*n+jpos] = t;
				}
				fact = w[my_rank*n+j] * D->value[jrow];
				w[my_rank*n+j] = fact;
				iw[my_rank*n+jrow] = -1;

				for(k=0;k<U->nnz[jrow];k++)
				{
					col = U->index[jrow][k];
					jpos = iw[my_rank*n+col];
					lxu = -fact * U->value[jrow][k];

					if( fabs(lxu) < tolnorm && jpos==-1 ) continue;
					if( col >= i )
					{
						if( jpos == -1 )
						{
							lenu++;
							upos = i + lenu;
							jbuf[my_rank*n+upos] = col;
							iw[my_rank*n+col] = upos;
							w[my_rank*n+upos] = lxu;
						}
						else
						{
							w[my_rank*n+jpos] += lxu;
						}
					}
					else
					{
						if( jpos == -1 )
						{
							jbuf[my_rank*n+lenl] = col;
							iw[my_rank*n+col] = lenl;
							w[my_rank*n+lenl] = lxu;
							lenl++;
						}
						else
						{
							w[my_rank*n+jpos] += lxu;
						}
					}
				}
			}

			iw[my_rank*n+i] = -1;
			for(j=0;j<lenu;j++)
			{
				iw[ my_rank*n+jbuf[my_rank*n+i+j+1] ] = -1;
			}

			if( pivot_tol>0.0 )
			{
			    w[my_rank*n+i] =
			        lis_precon_regularize_ilu_pivot(
			            w[my_rank*n+i],
			            row_scale,
			            pivot_tol);
			}
			D->value[i] = 1.0 / w[my_rank*n+i];


			len = _min(lfil,lenl);
			for(j=0;j<lenl;j++)
			{
				wn[my_rank*n+j] = fabs(w[my_rank*n+j]);
				iw[my_rank*n+j] = j;
			}
			lis_sort_di(0,lenl-1,&wn[my_rank*n],&iw[my_rank*n]);
			lis_sort_i(0,len-1,&iw[my_rank*n]);
			
			L->nnz[i] = len;
			if( len>0 )
			{
				L->index[i] = (LIS_INT *)malloc(len*sizeof(LIS_INT));
				L->value[i] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
			}
			for(j=0;j<len;j++)
			{
				jpos = iw[my_rank*n+j];
				L->index[i][j] = jbuf[my_rank*n+jpos];
				L->value[i][j] = w[my_rank*n+jpos];
			}
			for(j=0;j<lenl;j++) iw[my_rank*n+j] = -1;

			len = _min(lfil,lenu);
			for(j=0;j<lenu;j++)
			{
				wn[my_rank*n+j] = fabs(w[my_rank*n+i+j+1]);
				iw[my_rank*n+j] = i+j+1;
			}
			lis_sort_di(0,lenu-1,&wn[my_rank*n],&iw[my_rank*n]);
			lis_sort_i(0,len-1,&iw[my_rank*n]);
			
			U->nnz[i] = len;
			if( len>0 )
			{
				U->index[i] = (LIS_INT *)malloc(len*sizeof(LIS_INT));
				U->value[i] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
			}
			for(j=0;j<len;j++)
			{
				jpos = iw[my_rank*n+j];
				U->index[i][j] = jbuf[my_rank*n+jpos];
				U->value[i][j] = w[my_rank*n+jpos];
			}
			for(j=0;j<lenu;j++) iw[my_rank*n+j] = -1;
		}
	}


	for(j=0;j<nprocs;j++)
	{
		lis_precon_row_source_destroy(
			&rows[j]);
	}

	lis_free(rows);

	if( row_err )
	{
		lis_free2(
			4,
			w,
			iw,
			wn,
			jbuf);

		lis_matrix_ilu_destroy(L);
		lis_matrix_ilu_destroy(U);
		lis_vector_destroy(D);

		return row_err;
	}

	precon->L  = L;
	precon->U  = U;
	precon->D  = D;

	lis_free2(4,w,iw,wn,jbuf);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_INT	err;
	LIS_INT	i,j,k;
	LIS_INT	n,lfil,len,matrix_nnz,get_err;
	LIS_SCALAR gamma,t,tol,m;
	LIS_MATRIX A;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;

	LIS_REAL tnorm,tolnorm,pivot_tol,row_scale,rt;
	LIS_SCALAR fact,lxu,*wn,*w;
	LIS_INT	lenu,lenl,col,jpos,jrow,upos;
	LIS_INT	*jbuf,*iw;
	LIS_PRECON_ROW_SOURCE row;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	tol    = solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
	pivot_tol = solver->params[
	    LIS_PARAMS_ILU_PIVOT_TOL-LIS_OPTIONS_LEN];

	err = lis_precon_check_ilu_pivot_tol(
	    A,pivot_tol);
	if( err ) return err;
	m      = solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
	gamma  = solver->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN];
	err = lis_precon_row_source_count_nnz(
		A,
		&matrix_nnz);
	if( err ) return err;
	lfil   = (LIS_INT)(((double)matrix_nnz/(2.0*n))*m);

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

	w   = (LIS_SCALAR *)lis_malloc((n+1)*sizeof(LIS_SCALAR),"lis_precon_create_ilut_csr::w");
	if( w==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	wn = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_ilut_csr::w");
	if( wn==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	jbuf   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_ilut_csr::iw");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_ilut_csr::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}



	err = lis_precon_row_source_init(
		A,
		&row);

	if( err )
	{
		lis_free2(
			4,
			w,
			iw,
			wn,
			jbuf);

		lis_matrix_ilu_destroy(L);
		lis_matrix_ilu_destroy(U);
		lis_vector_destroy(D);

		return err;
	}

	for(i=0;i<n;i++) iw[i] = -1;

	for(i=0;i<n;i++)
	{
		get_err =
			lis_precon_row_source_get(
				&row,
				i);

		if( get_err )
		{
			lis_precon_row_source_destroy(
				&row);

			lis_free2(
				4,
				w,
				iw,
				wn,
				jbuf);

			lis_matrix_ilu_destroy(L);
			lis_matrix_ilu_destroy(U);
			lis_vector_destroy(D);

			return get_err;
		}

		tnorm = 0;
		row_scale = 0.0;

		for(j=0;j<row.nnz;j++)
		{
			rt = fabs(row.value[j]);

			if( rt>row_scale )
			{
				row_scale = rt;
			}

			tnorm += rt;
		}

		tnorm   = tnorm / (double)row.nnz;
		tolnorm = tol * tnorm;

		lenu = 0;
		lenl = 0;

		jbuf[i] = i;
		w[i] = 0;
		iw[i] = i;

		for(j=0;j<row.nnz;j++)
		{
			col = row.index[j];

			#ifdef USE_MPI
				if( col>n-1 ) continue;
			#endif

			t = row.value[j];

			if( col < i )
			{
				jbuf[lenl] = col;
				iw[col] = lenl;
				w[lenl] = t;
				lenl++;
			}
			else if( col == i )
			{
				w[i] = t;
			}
			else
			{
				lenu++;
				jpos = i + lenu;
				jbuf[jpos] = col;
				iw[col] = jpos;
				w[jpos] = t;
			}
		}

		j = -1;
		len = 0;

		while( ++j < lenl )
		{
			jrow = jbuf[j];
			jpos = j;
			for(k=j+1;k<lenl;k++)
			{
				if( jbuf[k]<jrow )
				{
					jrow = jbuf[k];
					jpos = k;
				}
			}
			if( jpos!=j )
			{
				col = jbuf[j];
				jbuf[j] = jbuf[jpos];
				jbuf[jpos] = col;
				iw[jrow] = j;
				iw[col] = jpos;
				t = w[j];
				w[j] = w[jpos];
				w[jpos] = t;
			}
			fact = w[j] * D->value[jrow];
			w[j] = fact;
			iw[jrow] = -1;

			for(k=0;k<U->nnz[jrow];k++)
			{
				col = U->index[jrow][k];
				jpos = iw[col];
				lxu = -fact * U->value[jrow][k];

				if( fabs(lxu) < tolnorm && jpos==-1 ) continue;
				if( col >= i )
				{
					if( jpos == -1 )
					{
						lenu++;
						upos = i + lenu;
						jbuf[upos] = col;
						iw[col] = upos;
						w[upos] = lxu;
					}
					else
					{
						w[jpos] += lxu;
					}
				}
				else
				{
					if( jpos == -1 )
					{
						jbuf[lenl] = col;
						iw[col] = lenl;
						w[lenl] = lxu;
						lenl++;
					}
					else
					{
						w[jpos] += lxu;
					}
				}
			}
/*			for(kk=0;kk<bs;kk++)
			{
				w[bs*len+kk] = -buf_fact[kk];
			}
			jbuf[len] = jrow;
			len++;*/
		}

		iw[i] = -1;
		for(j=0;j<lenu;j++)
		{
			iw[ jbuf[i+j+1] ] = -1;
		}

		if( pivot_tol>0.0 )
		{
		    w[i] =
		        lis_precon_regularize_ilu_pivot(
		            w[i],
		            row_scale,
		            pivot_tol);
		}
		D->value[i] = 1.0 / w[i];


		len = _min(lfil,lenl);
		for(j=0;j<lenl;j++)
		{
			wn[j] = fabs(w[j]);
			iw[j] = j;
		}
		lis_sort_di(0,lenl-1,wn,iw);
		lis_sort_i(0,len-1,iw);
		
		L->nnz[i] = len;
		if( len>0 )
		{
			L->index[i] = (LIS_INT *)malloc(len*sizeof(LIS_INT));
			L->value[i] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jpos = iw[j];
			L->index[i][j] = jbuf[jpos];
			L->value[i][j] = w[jpos];
		}
		for(j=0;j<lenl;j++) iw[j] = -1;

		len = _min(lfil,lenu);
		for(j=0;j<lenu;j++)
		{
			wn[j] = fabs(w[i+j+1]);
			iw[j] = i+j+1;
		}
		lis_sort_di(0,lenu-1,wn,iw);
		lis_sort_i(0,len-1,iw);
		
		U->nnz[i] = len;
		if( len>0 )
		{
			U->index[i] = (LIS_INT *)malloc(len*sizeof(LIS_INT));
			U->value[i] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jpos = iw[j];
			U->index[i][j] = jbuf[jpos];
			U->value[i][j] = w[jpos];
		}
		for(j=0;j<lenu;j++) iw[j] = -1;
	}

	lis_precon_row_source_destroy(
		&row);

	precon->L  = L;
	precon->U  = U;
	precon->D  = D;

	lis_free2(4,w,iw,wn,jbuf);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_ilut_csr"
LIS_INT lis_psolve_ilut_csr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
#ifdef _OPENMP
	LIS_INT i,j,jj,n;
	LIS_INT nprocs,my_rank,is,ie;
	LIS_SCALAR *x;
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

#ifdef USE_QUAD_PRECISION
	if( precon->half!=NULL &&
		B->precision==LIS_PRECISION_DEFAULT )
#else
	if( precon->half!=NULL )
#endif
	{
		return lis_half_ilu_psolve(
			solver,
			B,
			X);
	}

	L = precon->L;
	U = precon->U;
	D = precon->D;
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
			nprocs = omp_get_max_threads();
			#pragma omp parallel private(i,j,jj,is,ie,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				for(i=is; i<ie; i++)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj    = L->index[i][j];
						x[i] -= L->value[i][j] * x[jj];
					}
				}
				for(i=ie-1; i>=is; i--)
				{
					for(j=0;j<U->nnz[i];j++)
					{
						jj    = U->index[i][j];
						x[i] -= U->value[i][j] * x[jj];
					}
					x[i] = D->value[i] * x[i];
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
				for(i=is; i<ie; i++)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj    = L->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-L->value[i][j]);
						#else
							LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-L->value[i][j]);
						#endif
/*						x[i] -= L->value[i][j] * x[jj];*/
					}
				}
				for(i=ie-1; i>=is; i--)
				{
					for(j=0;j<U->nnz[i];j++)
					{
						jj    = U->index[i][j];
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
/*					x[i] = D->value[i] * x[i];*/
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

#ifdef USE_QUAD_PRECISION
	if( precon->half!=NULL &&
		B->precision==LIS_PRECISION_DEFAULT )
#else
	if( precon->half!=NULL )
#endif
	{
		return lis_half_ilu_psolve(
			solver,
			B,
			X);
	}

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
					jj    = L->index[i][j];
					x[i] -= L->value[i][j] * x[jj];
				}
			}
			for(i=n-1; i>=0; i--)
			{
				for(j=0;j<U->nnz[i];j++)
				{
					jj    = U->index[i][j];
					x[i] -= U->value[i][j] * x[jj];
				}
				x[i] = D->value[i] * x[i];
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
#define __FUNC__ "lis_psolveh_ilut_csr"
LIS_INT lis_psolveh_ilut_csr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
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

#ifdef USE_QUAD_PRECISION
	if( precon->half!=NULL &&
		B->precision==LIS_PRECISION_DEFAULT )
#else
	if( precon->half!=NULL )
#endif
	{
		return lis_half_ilu_psolveh(
			solver,
			B,
			X);
	}

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
					x[i] = conj(D->value[i]*x[i]);
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
/*						x[jj] -= conj(U->value[i][j]) * x[i];*/
					}
				}
				for(i=ie-1;i>=is;i--)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj     = L->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(L->value[i][j]));
						#else
							LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-conj(L->value[i][j]));
						#endif
/*						x[jj] -= conj(L->value[i][j]) * x[i];*/
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

#ifdef USE_QUAD_PRECISION
	if( precon->half!=NULL &&
		B->precision==LIS_PRECISION_DEFAULT )
#else
	if( precon->half!=NULL )
#endif
	{
		return lis_half_ilu_psolveh(
			solver,
			B,
			X);
	}

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
				x[i] = D->value[i]*x[i];
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
/*				x[i] = D->value[i]*x[i];*/
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
						LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
					#else
						LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
					#endif
/*					x[jj] -= L->value[i][j] * x[i];*/
				}
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

