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
 * lis_precon_create
 * lis_psolve
 * lis_psolveh
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_precon_create_sainv"
LIS_INT lis_precon_create_sainv(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	err;
	LIS_MATRIX A,B;


	LIS_DEBUG_FUNC_IN;

	switch( solver->A->matrix_type )
	{
	case LIS_MATRIX_CSR:
		err = lis_precon_create_sainv_csr(solver,precon);
		break;
	default:
		A = solver->A;
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,LIS_MATRIX_CSR);
		err = lis_matrix_convert(A,B);
		if( err ) return err;
		solver->A = B;
		err = lis_precon_create_sainv_csr(solver,precon);
		lis_matrix_destroy(B);
		solver->A = A;
		break;
	}

	#ifndef USE_QUAD_PRECISION
		err = lis_vector_duplicate(solver->A,&precon->temp);
	#else
		if( solver->precision==LIS_PRECISION_DEFAULT )
		{
			err = lis_vector_duplicate(solver->A,&precon->temp);
		}
		else
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,solver->A,&precon->temp);
		}
	#endif
	if( err ) return err;

	precon->A       = solver->A;
	precon->is_copy = LIS_FALSE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/********************************************
 W  = I
 Z  = I
 for i=1,...,n
   l    = A * Z_i
   u    = W_i^T * A
   D_ii = u * Z_i
   for j>i, l_j!=0
     W_j = W_j - drop( (l_j/D_ii)*W_i, tol )
   for j>i, u_j!=0
     Z_j = Z_j - drop( (u_j/D_ii)*Z_i, tol )
 ********************************************/
#if 1
#undef __FUNC__
#define __FUNC__ "lis_precon_create_sainv_csr"
LIS_INT lis_precon_create_sainv_csr(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	err;
	LIS_INT	i,j,k,ii,jj,ik,jk;
	LIS_INT	n,annz,cl,cu;
	LIS_INT	*ww,*il,*iu;
	LIS_SCALAR t,dd;
	LIS_REAL tol,nrm;
	LIS_SCALAR *d,*l,*u;
	LIS_MATRIX A,B;
	LIS_MATRIX_ILU W,Z;
	LIS_VECTOR D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	tol    = solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
	annz   = A->n / 10;

	W      = NULL;
	ww     = NULL;
	d      = NULL;
	l      = NULL;
	u      = NULL;
	il     = NULL;
	iu     = NULL;

	err = lis_matrix_ilu_create(n,1,&W);
	if( err ) return err;
	err = lis_matrix_ilu_create(n,1,&Z);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(W);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(Z);
	if( err ) return err;
	err = lis_vector_duplicate(A,&D);
	if( err ) return err;
	d = D->value;
	err = lis_matrix_ilu_premalloc(annz,W);
	if( err ) return err;
	err = lis_matrix_ilu_premalloc(annz,Z);
	if( err ) return err;
	l   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_sainv_csr::l");
	if( l==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	u   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_sainv_csr::u");
	if( u==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	il   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_sainv_csr::il");
	if( il==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iu   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_sainv_csr::iu");
	if( iu==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	ww   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_sainv_csr::ww");
	if( ww==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	err = lis_matrix_duplicate(A,&B);
	if( err ) return err;
	err = lis_matrix_convert_csr2csc(A,B);
	if( err )
	{
		return err;
	}

	for(i=0;i<n;i++) ww[i] = 0;
	for(i=0;i<n;i++)
	{
		W->value[i][0] = 1.0;
		W->index[i][0] = i;
		W->nnz[i]      = 1;
		Z->value[i][0] = 1.0;
		Z->index[i][0] = i;
		Z->nnz[i]      = 1;
	}
	for(i=0; i<n; i++)
	{
		/* nrm_inf(A[i,:]) */
		nrm = 0.0;
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			nrm = _max(nrm,fabs(A->value[j]));
		}
		nrm = 1.0/nrm;

		/* l = AZ_i */
		cl = 0;
		memset(l,0,n*sizeof(LIS_SCALAR));
		for(k=0;k<Z->nnz[i];k++)
		{
			ii = Z->index[i][k];
			for(j=B->ptr[ii];j<B->ptr[ii+1];j++)
			{
				jj     = B->index[j];
				if( jj>i )
				{
					l[jj] += B->value[j]*Z->value[i][k];
					if( ww[jj]==0 )
					{
						ww[jj]   = 1;
						il[cl++] = jj;
					}
				}
			}
		}
		for(k=0;k<cl;k++) ww[il[k]] = 0;

		/* u = W_i'A */
		cu = 0;
		memset(u,0,n*sizeof(LIS_SCALAR));
		for(k=0;k<W->nnz[i];k++)
		{
			ii = W->index[i][k];
			for(j=A->ptr[ii];j<A->ptr[ii+1];j++)
			{
				jj     = A->index[j];
				#ifdef USE_MPI
					if( jj>n-1 ) continue;
				#endif
				u[jj] += A->value[j]*W->value[i][k];
				if( jj>i && ww[jj]==0 )
				{
					ww[jj]   = 1;
					iu[cu++] = jj;
				}
			}
		}
		for(k=0;k<cu;k++) ww[iu[k]] = 0;

		/* d_ii = uZ_i or W_i'l  */
		t = 0.0;
		for(k=0;k<Z->nnz[i];k++)
		{
			t += u[Z->index[i][k]]*Z->value[i][k];
		}
		d[i] = 1.0/t;

		/* for j>i, l_j!=0            */
		/* w_j = w_j - (l_j/d_ii)*w_i */
		for(jj=0;jj<cl;jj++)
		{
			j = il[jj];
			dd = l[j]*d[i];
			for(k=0;k<W->nnz[j];k++)
			{
				ww[W->index[j][k]] = k+1;
			}
			for(ik=0;ik<W->nnz[i];ik++)
			{
				jk = ww[W->index[i][ik]];
				if( jk!=0 )
				{
					t = dd*W->value[i][ik];
					if( fabs(t)*nrm > tol )
					{
						W->value[j][jk-1] -= t;
					}
				}
				else
				{
					t = dd*W->value[i][ik];
					if( fabs(t)*nrm > tol )
					{
						if( W->nnz[j] == W->nnz_ma[j] )
						{
							W->nnz_ma[j] += annz;
							err = lis_matrix_ilu_realloc(j,W->nnz_ma[j],W);
							if( err ) return err;
						}
						jk                = W->nnz[j];
						W->index[j][jk] = W->index[i][ik];
						W->value[j][jk] = -t;
						W->nnz[j]++;
					}
				}
			}
			for(k=0;k<W->nnz[j];k++)
			{
				ww[W->index[j][k]] = 0;
			}
		}

		/* for j>i, u_j!=0            */
		/* z_j = z_j - (u_j/d_ii)*z_i */
		for(jj=0;jj<cu;jj++)
		{
			j = iu[jj];
			dd = u[j]*d[i];
			for(k=0;k<Z->nnz[j];k++)
			{
				ww[Z->index[j][k]] = k+1;
			}
			for(ik=0;ik<Z->nnz[i];ik++)
			{
				jk = ww[Z->index[i][ik]];
				if( jk!=0 )
				{
					t = dd*Z->value[i][ik];
					if( fabs(t)*nrm > tol )
					{
						Z->value[j][jk-1] -= t;
					}
				}
				else
				{
					t = dd*Z->value[i][ik];
					if( fabs(t)*nrm > tol )
					{
						if( Z->nnz[j] == Z->nnz_ma[j] )
						{
							Z->nnz_ma[j] += annz;
							err = lis_matrix_ilu_realloc(j,Z->nnz_ma[j],Z);
							if( err ) return err;
						}
						jk                = Z->nnz[j];
						Z->index[j][jk] = Z->index[i][ik];
						Z->value[j][jk] = -t;
						Z->nnz[j]++;
					}
				}
			}
			for(k=0;k<Z->nnz[j];k++)
			{
				ww[Z->index[j][k]] = 0;
			}
		}
	}

	lis_matrix_destroy(B);
	lis_free2(5,l,u,ww,il,iu);


	precon->L  = W;
	precon->U  = Z;
	precon->D  = D;


	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#else
/********************************************
 for i=1,...,n
   W_i  = e_i
   Z_i  = e_i
   r    = e_i^T * A
   c    = A * e_i
   for j=1,...,i-1
     W_i = W_i - drop( (r*Z_j/D_jj)*W_j, tol )
   for j=1,...,i-1
     Z_i = Z_i - drop( (W_j^T*c/D_jj)*Z_j, tol )
   l    = A * Z_i
   D_ii = W_i^T * A * Z_i
 ********************************************/
#undef __FUNC__
#define __FUNC__ "lis_precon_create_sainv_csr"
LIS_INT lis_precon_create_sainv_csr(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	err;
	LIS_INT	i,j,k,ii,jj,len,lfil;
	LIS_INT	n,nnz,annz,cl,cu,cc,m;
	LIS_INT	*wu,*wl,*il,*iu,*ic,*pc;
	LIS_SCALAR t,v;
	LIS_REAL tol,tol_dd,nrm;
	LIS_SCALAR *d,*r,*c,*l,*u,*tmp;
	LIS_MATRIX A,B;
	LIS_MATRIX_ILU W,Z;
	LIS_VECTOR D;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nnz    = A->nnz;
	tol    = solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
	m      = solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
	annz   = 10+A->nnz / A->n;
	lfil   = (LIS_INT)((double)A->nnz/(2.0*n))*m;

	W      = NULL;
	Z      = NULL;
	wu     = NULL;
	wl     = NULL;
	d      = NULL;
	l      = NULL;
	u      = NULL;
	il     = NULL;
	iu     = NULL;

	err = lis_matrix_ilu_create(n,1,&W);
	if( err ) return err;
	err = lis_matrix_ilu_create(n,1,&Z);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(W);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(Z);
	if( err ) return err;
	err = lis_vector_duplicate(A,&D);
	if( err ) return err;
	d = D->value;

	tmp   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_sainv_csr::l");
	if( tmp==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	r   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_sainv_csr::l");
	if( r==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	c   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_sainv_csr::u");
	if( c==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	l   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_sainv_csr::l");
	if( l==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	u   = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_sainv_csr::u");
	if( u==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	il   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_sainv_csr::il");
	if( il==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	iu   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_sainv_csr::iu");
	if( iu==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	ic   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_sainv_csr::iu");
	if( ic==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	wu   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_sainv_csr::ww");
	if( wu==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	wl   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_sainv_csr::ww");
	if( wl==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	pc   = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_precon_create_sainv_csr::iu");
	if( pc==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	lis_matrix_sort_csr(A);
	err = lis_matrix_duplicate(A,&B);
	if( err ) return err;
	err = lis_matrix_convert_csr2csc(A,B);
	if( err ) return err;

	for(i=0;i<n;i++)
	{
		wu[i] = 0;
		wl[i] = 0;
		pc[i] = A->ptr[i];
	}
	for(i=0; i<n; i++)
	{
		/* nrm_inf(A[i,:]) */
		nrm = 0.0;
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			nrm = _max(nrm,fabs(A->value[j]));
		}
		tol_dd = nrm * tol;

		/* l = e_i  */
		/* u = e_i  */
		l[i]  = 1.0;
		u[i]  = 1.0;
		il[0] = i;
		iu[0] = i;
		cl    = 1;
		cu    = 1;
		wu[i] = 1;
		wl[i] = 1;
		cc    = 0;

		/* r = e_i^T*A */
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			jj    = A->index[j];
			r[jj] = A->value[j];
		}
		/* c = A_i = A*e_i */
		for(j=B->ptr[i];j<B->ptr[i+1];j++)
		{
			jj    = B->index[j];
			c[jj] = B->value[j];
		}

	    /* W_i = W_i - (r*Z_j/D_jj)*W_j */
		for(j=0;j<i;j++)
		{
			t = 0.0;
			for(k=0;k<Z->nnz[j];k++)
			{
				t += r[Z->index[j][k]]*Z->value[j][k];
			}
			t = t * d[j];
			if( fabs(t) > tol_dd )
			{
				for(k=0;k<W->nnz[j];k++)
				{
					v      = t * W->value[j][k];
					if( fabs(v) > tol_dd )
					{
						jj     = W->index[j][k];
						if( wl[jj]==1 )
						{
							l[jj] -= v;
						}
						else
						{
							l[jj]    = -v;
							il[cl++] = jj;
							wl[jj]   = 1;
						}
					}
				}
			}
		}

		/* Z_i = Z_i - (W_j^T*c/D_jj)*Z_j */
		for(j=0;j<i;j++)
		{
			t = 0.0;
			for(k=0;k<W->nnz[j];k++)
			{
				t += c[W->index[j][k]]*W->value[j][k];
			}
			t = t * d[j];
			if( fabs(t) > tol_dd )
			{
				for(k=0;k<Z->nnz[j];k++)
				{
					v      = t * Z->value[j][k];
					if( fabs(v) > tol_dd )
					{
						jj     = Z->index[j][k];
						if( wu[jj]==1 )
						{
							u[jj] -= v;
						}
						else
						{
							u[jj]    = -v;
							iu[cu++] = jj;
							wu[jj]   = 1;
						}
					}
				}
			}
		}
/*
		len = _min(lfil,cl);
		for(j=0;j<cl;j++) tmp[j] = fabs(l[il[j]]);
		lis_sort_di(0,cl-1,tmp,il);
		lis_sort_i(0,len-1,il);
		cl = len;
		*/
		/*
		k = cl;
		for(j=0;j<cl;j++)
		{
			if( fabs(l[il[j]])<= tol_dd )
			{
				wl[il[j]] = 0;
				il[j] = n;
				k--;
			}
		}
		lis_sort_i(0,cl-1,il);
		cl = k;
		

		k = cu;
		for(j=0;j<cu;j++)
		{
			if( fabs(u[iu[j]])<= tol_dd )
			{
				wu[iu[j]] = 0;
				iu[j] = n;
				k--;
			}
		}
		lis_sort_i(0,cu-1,iu);
		cu = k;
		*/

		W->nnz[i] = cl;
		if( cl > 0 )
		{
			W->index[i] = (LIS_INT *)malloc(cl*sizeof(LIS_INT));
			W->value[i] = (LIS_SCALAR *)malloc(cl*sizeof(LIS_SCALAR));
			memcpy(W->index[i],il,cl*sizeof(LIS_INT));
			for(j=0;j<cl;j++)
			{
				W->value[i][j] = l[il[j]];
			}
		}
		Z->nnz[i] = cu;
		if( cu > 0 )
		{
			Z->index[i] = (LIS_INT *)malloc(cu*sizeof(LIS_INT));
			Z->value[i] = (LIS_SCALAR *)malloc(cu*sizeof(LIS_SCALAR));
			memcpy(Z->index[i],iu,cu*sizeof(LIS_INT));
			for(j=0;j<cu;j++)
			{
				Z->value[i][j] = u[iu[j]];
			}
		}

		for(j=A->ptr[i];j<A->ptr[i+1];j++) r[A->index[j]] = 0.0;
		for(j=B->ptr[i];j<B->ptr[i+1];j++) c[B->index[j]] = 0.0;
		for(j=0;j<cl;j++)
		{
			wl[il[j]] = 0;
			l[il[j]] = 0.0;
		}
		for(j=0;j<cu;j++)
		{
			wu[iu[j]] = 0;
		}

		/* D_ii = W_i^T * A * Z_i */
		cl = 0;
		for(k=0;k<Z->nnz[i];k++)
		{
			ii = Z->index[i][k];
			for(j=B->ptr[ii];j<B->ptr[ii+1];j++)
			{
				jj     = B->index[j];
				if( wl[jj]==0 )
				{
					l[jj] = B->value[j]*Z->value[i][k];
					wl[jj]   = 1;
					il[cl++] = jj;
				}
				else
				{
					l[jj] += B->value[j]*Z->value[i][k];
				}
			}
		}
		t = 0.0;
		for(j=0;j<W->nnz[i];j++)
		{
			k  = W->index[i][j];
			t += W->value[i][j] * l[k];
		}
		d[i] = 1.0 / t;
		for(j=0;j<cl;j++) wl[il[j]] = 0;

	}

	lis_matrix_destroy(B);
	lis_free2(11,r,c,il,l,wl,iu,u,wu,ic,pc,tmp);


	precon->L  = W;
	precon->U  = Z;
	precon->D  = D;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#endif

#undef __FUNC__
#define __FUNC__ "lis_psolve_sainv"
LIS_INT lis_psolve_sainv(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_INT i,n;
	LIS_MATRIX A;
	LIS_MATRIX_ILU W,Z;
	LIS_VECTOR t,d;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;

	/*
	 *  x  = Mb
	 *  M  = ZD^{-1}W'
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	A = precon->A;
	W = precon->L;
	Z = precon->U;
	d = precon->D;
	t = precon->temp;
	n = precon->L->n;

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_matvech_ilu(A,W,B,X);
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<n;i++)
			{
				t->value[i] = X->value[i]*d->value[i];
			}
			lis_matvec_ilu(A,Z,t,X);
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_matvech_ilu(A,W,B,X);
			#ifdef _OPENMP
			#ifndef USE_SSE2
				#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel for private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			#endif
			for(i=0;i<n;i++)
			{
				#ifndef USE_SSE2
					LIS_QUAD_MULD(t->value[i],t->value_lo[i],X->value[i],X->value_lo[i],d->value[i]);
				#else
					LIS_QUAD_MULD_SSE2(t->value[i],t->value_lo[i],X->value[i],X->value_lo[i],d->value[i]);
				#endif
				/* t->value[i] = X->value[i]*d->value[i]; */
			}
			lis_matvec_ilu(A,Z,t,X);
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_sainv"
LIS_INT lis_psolveh_sainv(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_INT i,n;
	LIS_MATRIX A;
	LIS_MATRIX_ILU W,Z;
	LIS_VECTOR t,d;
	LIS_PRECON precon;
	LIS_QUAD_DECLAR;

	/*
	 *  x  = M'b
	 *  M' = WD^{-1}Z'
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	A = precon->A;
	W = precon->L;
	Z = precon->U;
	d = precon->D;
	t = precon->temp;
	n = precon->L->n;

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_matvech_ilu(A,Z,B,X);
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<n;i++)
			{
				t->value[i] = X->value[i]*conj(d->value[i]);
			}
			lis_matvec_ilu(A,W,t,X);
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_matvech_ilu(A,Z,B,X);
			#ifdef _OPENMP
			#ifndef USE_SSE2
				#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel for private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			#endif
			for(i=0;i<n;i++)
			{
				#ifndef USE_SSE2
			  		LIS_QUAD_MULD(t->value[i],t->value_lo[i],X->value[i],X->value_lo[i],conj(d->value[i]));
				#else
					LIS_QUAD_MULD_SSE2(t->value[i],t->value_lo[i],X->value[i],X->value_lo[i],conj(d->value[i]));
				#endif
				/* t->value[i] = X->value[i]*conj(d->value[i]); */
			}
			lis_matvec_ilu(A,W,t,X);
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
