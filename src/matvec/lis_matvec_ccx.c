/*
 * High-throughput CalculiX SpMV for LIS.
 *
 * The algorithm deliberately uses two gather passes instead of scatter
 * updates.  Therefore OpenMP needs no atomics and no per-thread y[] copy.
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
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lislib.h"

#ifdef USE_CCX

void lis_matvec_ccx(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k,p,n,base,nnz;
	LIS_SCALAR sum;
	LIS_INT *tptr, *tcol, *tindex, *jq, *irow;
	LIS_SCALAR *ad, *au;

	n = A->n;
	base = A->ccx_index_base;
	nnz = A->ccx_nnz_offdiag;
	tptr = A->ccx_tptr;
	tcol = A->ccx_tcol;
	tindex = A->ccx_tindex;
	jq = A->ccx_jq;
	irow = A->ccx_irow;
	ad = A->ccx_ad;
	au = A->ccx_au;

	/*
	 * Single-threaded fast path.
	 *
	 * With one OpenMP thread the two gather passes below are unnecessary.
	 * Walk native CalculiX CCS storage once: the lower contribution is a
	 * scatter to y[i], while the matching upper contribution is accumulated
	 * locally into y[j]. No atomics are required because this path is
	 * strictly single-threaded.
	 */
#ifdef _OPENMP
	if (omp_get_max_threads() == 1)
#endif
	{
		LIS_SCALAR xj, aij, aji;

		for (i=0; i<n; i++)
			y[i] = ad[i] * x[i];

		for (j=0; j<n; j++)
		{
			sum = y[j];
			xj = x[j];

			for (k=jq[j]-base; k<jq[j+1]-base; k++)
			{
				i = irow[k]-base;
				aij = au[k];

				/* Stored lower triangle: A(i,j) * x(j). */
				y[i] += aij * xj;

				/* Matching upper triangle: A(j,i) * x(i). */
				if (A->ccx_nasym)
					aji = au[nnz+k];
				else
					aji = aij;

				sum += aji * x[i];
			}
			y[j] = sum;
		}
		return;
	}

	/* Diagonal plus the stored (lower) triangle, gathered by row. */
#ifdef _OPENMP
#pragma omp parallel for private(i,p,k,j,sum) schedule(static)
#endif
	for(i=0;i<n;i++)
	{
		sum = ad[i] * x[i];
		for(p=tptr[i];p<tptr[i+1];p++)
		{
			k = tindex[p];
			j = tcol[p];
			sum += au[k] * x[j];
		}
		y[i] = sum;
	}

	/* Matching upper triangle, gathered in native CCX column order. */
#ifdef _OPENMP
#pragma omp parallel for private(j,k,i,sum) schedule(static)
#endif
	for(j=0;j<n;j++)
	{
		sum = (LIS_SCALAR)0.0;
		for(k=jq[j]-base;k<jq[j+1]-base;k++)
		{
			i = irow[k]-base;
			if( A->ccx_nasym )
				sum += au[nnz+k] * x[i];
			else
				sum += au[k] * x[i];
		}
		y[j] += sum;
	}
}

void lis_matvech_ccx(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k,p,n,base,nnz;
	LIS_SCALAR sum;
	LIS_INT *tptr, *tcol, *tindex, *jq, *irow;
	LIS_SCALAR *ad, *au;

	n = A->n;
	base = A->ccx_index_base;
	nnz = A->ccx_nnz_offdiag;
	tptr = A->ccx_tptr;
	tcol = A->ccx_tcol;
	tindex = A->ccx_tindex;
	jq = A->ccx_jq;
	irow = A->ccx_irow;
	ad = A->ccx_ad;
	au = A->ccx_au;

	/*
	 * Serial conjugate-transpose fast path. For real CalculiX builds conj()
	 * is a no-op, but keeping the general form preserves LIS semantics.
	 */
#ifdef _OPENMP
	if (omp_get_max_threads() == 1)
#endif
	{
		LIS_SCALAR xj, aij, aji;

		for (i=0; i<n; i++)
			y[i] = conj(ad[i]) * x[i];

		for (j=0; j<n; j++)
		{
			sum = y[j];
			xj = x[j];

			for (k=jq[j]-base; k<jq[j+1]-base; k++)
			{
				i = irow[k]-base;
				aij = au[k];

				if (A->ccx_nasym)
					aji = au[nnz+k];
				else
					aji = aij;

				/* A^H(i,j) = conj(A(j,i)). */
				y[i] += conj(aji) * xj;

				/* A^H(j,i) = conj(A(i,j)). */
				sum += conj(aij) * x[i];
			}
			y[j] = sum;
		}
		return;
	}

#ifdef _OPENMP
#pragma omp parallel for private(i,p,k,j,sum) schedule(static)
#endif
	for(i=0;i<n;i++)
	{
		sum = conj(ad[i]) * x[i];
		for(p=tptr[i];p<tptr[i+1];p++)
		{
			k = tindex[p];
			j = tcol[p];
			if( A->ccx_nasym )
				sum += conj(au[nnz+k]) * x[j];
			else
				sum += conj(au[k]) * x[j];
		}
		y[i] = sum;
	}

#ifdef _OPENMP
#pragma omp parallel for private(j,k,i,sum) schedule(static)
#endif
	for(j=0;j<n;j++)
	{
		sum = (LIS_SCALAR)0.0;
		for(k=jq[j]-base;k<jq[j+1]-base;k++)
		{
			i = irow[k]-base;
			sum += conj(au[k]) * x[i];
		}
		y[j] += sum;
	}
}

#endif /* USE_CCX */