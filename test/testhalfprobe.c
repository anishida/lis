#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lislib.h"

extern LIS_INT lis_half_ilu_compress(
	LIS_MATRIX_ILU,
	LIS_MATRIX_ILU,
	LIS_HALF_ILU *,
	LIS_REAL *);

extern LIS_INT lis_half_ilu_destroy(
	LIS_HALF_ILU);

extern LIS_INT lis_half_ilu_assess(
	LIS_MATRIX_ILU,
	LIS_MATRIX_ILU,
	LIS_VECTOR,
	LIS_HALF_ILU,
	LIS_INT,
	LIS_INT *,
	LIS_REAL *,
	LIS_REAL *);


static int
set_one(
	LIS_MATRIX_ILU A,
	LIS_INT row,
	LIS_INT column,
	LIS_SCALAR value)
{
	A->index[row] =
		(LIS_INT *)malloc(sizeof(LIS_INT));

	A->value[row] =
		(LIS_SCALAR *)malloc(sizeof(LIS_SCALAR));

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

	LIS_MATRIX_ILU L,U;
	LIS_HALF_ILU half;
	LIS_VECTOR D;

	LIS_INT i;
	LIS_INT err;
	LIS_INT accepted;

	LIS_REAL underflow;
	LIS_REAL probe_error;
	LIS_REAL tolerance;

	int strict;
	int failed;

	L = NULL;
	U = NULL;
	half = NULL;
	D = NULL;
	failed = 0;

	strict =
		argc>1 &&
		strcmp(argv[1],"strict")==0;

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

	for(i=1;i<n;i++)
	{
		if( set_one(
				L,
				i,
				i-1,
				(LIS_SCALAR)
				(0.1 + 0.013*(double)i)) )
		{
			failed = 1;
			goto cleanup;
		}
	}

	for(i=0;i<n-1;i++)
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

	err = lis_vector_create(
		LIS_COMM_WORLD,
		&D);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_set_size(
		D,
		0,
		n);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	for(i=0;i<n;i++)
	{
		D->value[i] =
			(LIS_SCALAR)
			(1.0 /
			 (1.2 + 0.07*(double)i));
	}

	err = lis_half_ilu_compress(
		L,
		U,
		&half,
		&underflow);

	if( err )
	{
		fprintf(
			stderr,
			"compress rc=%d\n",
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

	if( err )
	{
		fprintf(
			stderr,
			"assess rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	if( strict )
	{
		if( tolerance>1.1e-12 )
		{
			fprintf(
				stderr,
				"strict tolerance not applied: %.17g\n",
				(double)tolerance);

			failed = 1;
			goto cleanup;
		}

		if( accepted )
		{
			fprintf(
				stderr,
				"strict probe unexpectedly accepted: "
				"error=%.17g tol=%.17g\n",
				(double)probe_error,
				(double)tolerance);

			failed = 1;
			goto cleanup;
		}
	}
	else
	{
		if( !accepted )
		{
			fprintf(
				stderr,
				"default probe rejected: "
				"error=%.17g tol=%.17g\n",
				(double)probe_error,
				(double)tolerance);

			failed = 1;
			goto cleanup;
		}

		if( !(probe_error>0.0) ||
			probe_error>tolerance )
		{
			fprintf(
				stderr,
				"unexpected probe error %.17g "
				"tol %.17g\n",
				(double)probe_error,
				(double)tolerance);

			failed = 1;
			goto cleanup;
		}

		/*
		 * NOPROBE must skip quality probing but retain
		 * all storage safety checks.
		 */
		accepted = LIS_FALSE;
		probe_error = -1.0;

		err = lis_half_ilu_assess(
			L,
			U,
			D,
			half,
			LIS_PRECON_HALF_NOPROBE,
			&accepted,
			&probe_error,
			&tolerance);

		if( err || !accepted )
		{
			fprintf(
				stderr,
				"NOPROBE failed rc=%d accepted=%d\n",
				(int)err,
				(int)accepted);

			failed = 1;
			goto cleanup;
		}
	}

	printf(
		"HALF_PROBE_C3 PASSED "
		"strict=%d error=%.17g tol=%.17g "
		"underflow=%.17g\n",
		strict,
		(double)probe_error,
		(double)tolerance,
		(double)underflow);


cleanup:

	lis_half_ilu_destroy(half);

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
