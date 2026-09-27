/*
 * Native CalculiX sparse matrix backend for LIS.
 *
 * The numerical arrays ad/au and the CalculiX structure jq/irow remain
 * owned by the caller.  LIS allocates only an integer row-gather map.
 *
 * Expected CalculiX layout:
 *   ad[i]                  diagonal
 *   jq[j]..jq[j+1]-1       off-diagonal entries stored by column
 *   irow[k]                row of entry k
 *
 * nasym == 0:
 *   au[k]                  stored triangle, mirrored symmetrically
 *
 * nasym != 0:
 *   au[k]                  lower counterpart
 *   au[nnz_offdiag + k]    upper counterpart
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
#include <string.h>
#include "lislib.h"

#ifdef USE_CCX

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_ccx"
LIS_INT lis_matrix_set_ccx(LIS_INT nnz_offdiag,
	LIS_SCALAR *ad, LIS_SCALAR *au,
	LIS_INT *jq, LIS_INT *irow,
	LIS_INT nasym, LIS_INT index_base,
	LIS_MATRIX A)
{
	LIS_INT err, n, i, j, k, p, base;
	LIS_INT kb, ke;
	LIS_INT *next = NULL;
	LIS_INT *tptr = NULL;
	LIS_INT *tcol = NULL;
	LIS_INT *tindex = NULL;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_SET);
	if( err ) return err;

	if( nnz_offdiag < 0 || ad==NULL || jq==NULL ||
	    (nnz_offdiag>0 && (au==NULL || irow==NULL)) )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"invalid CalculiX matrix arrays\n");
		return LIS_ERR_ILL_ARG;
	}
	if( index_base!=0 && index_base!=1 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,
			"CalculiX index base must be 0 or 1 (got %D)\n",index_base);
		return LIS_ERR_ILL_ARG;
	}

#ifdef USE_MPI
	if( A->nprocs != 1 )
	{
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
	}
#endif

	n = A->n;
	base = index_base;

	tptr = (LIS_INT *)lis_calloc((n+1)*sizeof(LIS_INT),
		"lis_matrix_set_ccx::tptr");
	if( tptr==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	if( nnz_offdiag>0 )
	{
		tcol = (LIS_INT *)lis_malloc(nnz_offdiag*sizeof(LIS_INT),
			"lis_matrix_set_ccx::tcol");
		tindex = (LIS_INT *)lis_malloc(nnz_offdiag*sizeof(LIS_INT),
			"lis_matrix_set_ccx::tindex");
		if( tcol==NULL || tindex==NULL )
		{
			LIS_SETERR_MEM(2*nnz_offdiag*sizeof(LIS_INT));
			goto oom;
		}
	}

	if( n>0 )
	{
		next = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),
			"lis_matrix_set_ccx::next");
		if( next==NULL )
		{
			LIS_SETERR_MEM(n*sizeof(LIS_INT));
			goto oom;
		}
	}

	/* Count entries of the stored triangle by destination row. */
	for(j=0;j<n;j++)
	{
		kb = jq[j]   - base;
		ke = jq[j+1] - base;
		if( kb<0 || ke<kb || ke>nnz_offdiag )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"invalid CalculiX jq[]\n");
			goto illarg;
		}
		for(k=kb;k<ke;k++)
		{
			i = irow[k] - base;
			if( i<0 || i>=n )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"invalid CalculiX irow[]\n");
				goto illarg;
			}
			tptr[i+1]++;
		}
	}

	for(i=0;i<n;i++)
		tptr[i+1] += tptr[i];
	for(i=0;i<n;i++)
		next[i] = tptr[i];

	/* Build row -> (source column, original au index) gather map. */
	for(j=0;j<n;j++)
	{
		kb = jq[j]   - base;
		ke = jq[j+1] - base;
		for(k=kb;k<ke;k++)
		{
			i = irow[k] - base;
			p = next[i]++;
			tcol[p] = j;
			tindex[p] = k;
		}
	}
	if( next ) lis_free(next);

	A->ccx_nnz_offdiag = nnz_offdiag;
	A->ccx_nasym = (nasym!=0);
	A->ccx_index_base = index_base;
	A->ccx_jq = jq;
	A->ccx_irow = irow;
	A->ccx_ad = ad;
	A->ccx_au = au;
	A->ccx_tptr = tptr;
	A->ccx_tcol = tcol;
	A->ccx_tindex = tindex;

	A->is_copy = LIS_FALSE;
	A->is_sorted = LIS_TRUE;
	A->nnz = n + 2*nnz_offdiag;

	/* set_ccx() returns a fully assembled native matrix. */
	A->matrix_type = LIS_MATRIX_CCX;
	A->status = LIS_MATRIX_CCX;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;

illarg:
	if( next ) lis_free(next);
	if( tptr ) lis_free(tptr);
	if( tcol ) lis_free(tcol);
	if( tindex ) lis_free(tindex);
	return LIS_ERR_ILL_ARG;

oom:
	if( next ) lis_free(next);
	if( tptr ) lis_free(tptr);
	if( tcol ) lis_free(tcol);
	if( tindex ) lis_free(tindex);
	return LIS_OUT_OF_MEMORY;
}
#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal_ccx"
LIS_INT lis_matrix_get_diagonal_ccx(LIS_MATRIX A, LIS_SCALAR d[])
{
	LIS_INT i;
	LIS_INT n = A->n;

#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	for(i=0;i<n;i++)
		d[i] = A->ccx_ad[i];

	return LIS_SUCCESS;
}

/*
 * Build a conventional full CSR copy of the live CalculiX matrix.
 *
 * This path is intentionally NOT used for Krylov matrix-vector products.
 * It exists so LIS preconditioners which require CSR (notably ILUT/ILUC/ILU)
 * can build their factors while the solver itself keeps using the zero-copy
 * LIS_MATRIX_CCX operator.
 */
#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_ccx2csr"
LIS_INT lis_matrix_convert_ccx2csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT err;
	LIS_INT n, nnz_off, nnz_full, base;
	LIS_INT i, j, k, kb, ke, p;
	LIS_INT *ptr = NULL, *index = NULL, *next = NULL;
	LIS_SCALAR *value = NULL;

	LIS_DEBUG_FUNC_IN;

	if( Ain==NULL || Aout==NULL || Ain->matrix_type!=LIS_MATRIX_CCX )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"CCX->CSR requires LIS_MATRIX_CCX input\n");
		return LIS_ERR_ILL_ARG;
	}

	n       = Ain->n;
	nnz_off = Ain->ccx_nnz_offdiag;
	base    = Ain->ccx_index_base;
	nnz_full = n + 2*nnz_off;

	err = lis_matrix_malloc_csr(n,nnz_full,&ptr,&index,&value);
	if( err ) return err;

	next = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),
		"lis_matrix_convert_ccx2csr::next");
	if( next==NULL )
	{
		lis_free2(3,ptr,index,value);
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}

	/* One diagonal entry in every row. */
	for(i=0;i<=n;i++) ptr[i] = 0;
	for(i=0;i<n;i++) ptr[i+1] = 1;

	/*
	 * Each native CCX off-diagonal entry produces two CSR entries:
	 *   (i,j) = lower/stored value
	 *   (j,i) = mirrored value, or the nonsymmetric upper value.
	 */
	for(j=0;j<n;j++)
	{
		kb = Ain->ccx_jq[j]   - base;
		ke = Ain->ccx_jq[j+1] - base;
		if( kb<0 || ke<kb || ke>nnz_off )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"invalid CalculiX jq[] in CCX->CSR\n");
			goto illarg;
		}
		for(k=kb;k<ke;k++)
		{
			i = Ain->ccx_irow[k] - base;
			if( i<0 || i>=n )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"invalid CalculiX irow[] in CCX->CSR\n");
				goto illarg;
			}
			ptr[i+1]++;
			ptr[j+1]++;
		}
	}

	for(i=0;i<n;i++) ptr[i+1] += ptr[i];
	for(i=0;i<n;i++) next[i] = ptr[i];

	/* Diagonal. */
	for(i=0;i<n;i++)
	{
		p = next[i]++;
		index[p] = i;
		value[p] = Ain->ccx_ad[i];
	}

	/* Lower/stored and matching upper entries. */
	for(j=0;j<n;j++)
	{
		kb = Ain->ccx_jq[j]   - base;
		ke = Ain->ccx_jq[j+1] - base;
		for(k=kb;k<ke;k++)
		{
			i = Ain->ccx_irow[k] - base;

			p = next[i]++;
			index[p] = j;
			value[p] = Ain->ccx_au[k];

			p = next[j]++;
			index[p] = i;
			value[p] = Ain->ccx_nasym
				? Ain->ccx_au[nnz_off+k]
				: Ain->ccx_au[k];
		}
	}

	lis_free(next);
	next = NULL;

	err = lis_matrix_set_csr(nnz_full,ptr,index,value,Aout);
	if( err ) goto fail_after_set;
	err = lis_matrix_assemble(Aout);
	if( err ) return err;

	Aout->is_sorted = LIS_FALSE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;

illarg:
	lis_free(next);
	lis_free2(3,ptr,index,value);
	return LIS_ERR_ILL_ARG;

fail_after_set:
	/* set_csr() did not take ownership on error. */
	lis_free2(3,ptr,index,value);
	return err;
}

#endif /* USE_CCX */