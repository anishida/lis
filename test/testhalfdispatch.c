#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lislib.h"


static int
set_one(
	LIS_MATRIX_ILU A,
	LIS_INT row,
	LIS_INT column,
	LIS_SCALAR value)
{
	A->index[row] =
		(LIS_INT *)malloc(
			sizeof(LIS_INT));

	A->value[row] =
		(LIS_SCALAR *)malloc(
			sizeof(LIS_SCALAR));

	if( A->index[row]==NULL ||
		A->value[row]==NULL )
	{
		free(A->index[row]);
		free(A->value[row]);

		A->index[row] = NULL;
		A->value[row] = NULL;

		return 1;
	}

	A->nnz[row] = 1;
	A->index[row][0] = column;
	A->value[row][0] = value;

	return 0;
}


int
main(
	int argc,
	char **argv)
{
	const LIS_INT n = 8;

	struct LIS_SOLVER_STRUCT solver_store;
	struct LIS_PRECON_STRUCT precon_store;
	struct LIS_MATRIX_STRUCT matrix_store;

	LIS_SOLVER solver;
	LIS_PRECON precon;

	LIS_MATRIX_ILU L,U;
	LIS_HALF_ILU half;

	LIS_VECTOR D;
	LIS_VECTOR B;
	LIS_VECTOR Xhalf;
	LIS_VECTOR Xnormal;

	LIS_REAL underflow;
	LIS_REAL probe_error;
	LIS_REAL tolerance;

	LIS_INT accepted;
	LIS_INT i;
	LIS_INT err;
	LIS_INT rank;
	LIS_INT is;
	LIS_INT ie;
	LIS_INT nprocs;

	double nr;
	double nd;
	double d;

	int failed;

	memset(
		&solver_store,
		0,
		sizeof(solver_store));

	memset(
		&precon_store,
		0,
		sizeof(precon_store));

	memset(
		&matrix_store,
		0,
		sizeof(matrix_store));

	solver = &solver_store;
	precon = &precon_store;

	L = NULL;
	U = NULL;
	half = NULL;

	D = NULL;
	B = NULL;
	Xhalf = NULL;
	Xnormal = NULL;

	failed = 0;

	if( lis_initialize(&argc,&argv) )
	{
		return 1;
	}

	if( lis_matrix_ilu_create(n,1,&L) ||
		lis_matrix_ilu_setCR(L) ||
		lis_matrix_ilu_create(n,1,&U) ||
		lis_matrix_ilu_setCR(U) )
	{
		failed = 1;
		goto cleanup;
	}

#ifdef _OPENMP
	nprocs = omp_get_max_threads();
#else
	nprocs = 1;
#endif

	/*
	 * OpenMP ILU(k) is block local.  Build the synthetic
	 * factors with the same partition contract so that this
	 * test measures FP16 quantization, not cross-block races.
	 */
	for(rank=0;rank<nprocs;rank++)
	{
		LIS_GET_ISIE(
			rank,
			nprocs,
			n,
			is,
			ie);

		for(i=is+1;i<ie;i++)
		{
			if( set_one(
					L,
					i,
					i-1,
					(LIS_SCALAR)
					(0.10 +
					 0.013*(double)i)) )
			{
				failed = 1;
				goto cleanup;
			}
		}

		for(i=is;i<ie-1;i++)
		{
			if( set_one(
					U,
					i,
					i+1,
					(LIS_SCALAR)
					(-0.08 -
					 0.009*(double)i)) )
			{
				failed = 1;
				goto cleanup;
			}
		}
	}

	if( lis_vector_create(
			LIS_COMM_WORLD,
			&D) ||
		lis_vector_set_size(
			D,
			0,
			n) ||
		lis_vector_create(
			LIS_COMM_WORLD,
			&B) ||
		lis_vector_set_size(
			B,
			0,
			n) ||
		lis_vector_duplicate(
			B,
			&Xhalf) ||
		lis_vector_duplicate(
			B,
			&Xnormal) )
	{
		failed = 1;
		goto cleanup;
	}

	for(i=0;i<n;i++)
	{
		D->value[i] =
			(LIS_SCALAR)
			(1.0 /
			 (1.2 +
			  0.07*(double)i));

		B->value[i] =
			(LIS_SCALAR)
			(0.2 +
			 0.11*(double)i);
	}

	err = lis_half_ilu_compress(
		L,
		U,
		&half,
		&underflow);

	if( err || half==NULL )
	{
		fprintf(
			stderr,
			"compress failed rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	err = lis_half_ilu_assess(
		L,
		U,
		D,
		half,
		LIS_PRECON_HALF_ON,
		&accepted,
		&probe_error,
		&tolerance);

	if( err || !accepted )
	{
		fprintf(
			stderr,
			"half assess failed rc=%d "
			"accepted=%d error=%.17g\n",
			(int)err,
			(int)accepted,
			(double)probe_error);

		failed = 1;
		goto cleanup;
	}

	matrix_store.n = n;

	solver->A = &matrix_store;
	solver->precon = precon;

	precon->L = L;
	precon->U = U;
	precon->D = D;
	precon->half = half;

	err = lis_psolve_iluk_csr(
		solver,
		B,
		Xhalf);

	if( err )
	{
		fprintf(
			stderr,
			"half dispatch rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	precon->half = NULL;

	err = lis_psolve_iluk_csr(
		solver,
		B,
		Xnormal);

	precon->half = half;

	if( err )
	{
		fprintf(
			stderr,
			"normal dispatch rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	nr = 0.0;
	nd = 0.0;

	for(i=0;i<n;i++)
	{
		d =
			(double)Xhalf->value[i] -
			(double)Xnormal->value[i];

		nr = hypot(
			nr,
			(double)Xnormal->value[i]);

		nd = hypot(
			nd,
			d);
	}

	if( nr==0.0 )
	{
		if( nd!=0.0 )
		{
			failed = 1;
		}
	}
	else if( nd/nr>1.0e-3 )
	{
		fprintf(
			stderr,
			"dispatch difference too large %.17g\n",
			nd/nr);

		failed = 1;
	}

	if( !failed )
	{
		printf(
			"HALF_ILUK_C4A PASSED "
			"relerr=%.17g probe=%.17g "
			"underflow=%.17g\n",
			nr==0.0 ? nd : nd/nr,
			(double)probe_error,
			(double)underflow);
	}


cleanup:

	lis_half_ilu_destroy(half);

	if( Xnormal!=NULL )
	{
		lis_vector_destroy(Xnormal);
	}

	if( Xhalf!=NULL )
	{
		lis_vector_destroy(Xhalf);
	}

	if( B!=NULL )
	{
		lis_vector_destroy(B);
	}

	if( D!=NULL )
	{
		lis_vector_destroy(D);
	}

	if( L!=NULL )
	{
		lis_matrix_ilu_destroy(L);
	}

	if( U!=NULL )
	{
		lis_matrix_ilu_destroy(U);
	}

	lis_finalize();

	return failed ? 1 : 0;
}
