#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lislib.h"

extern LIS_INT lis_half_ilu_compress(
	LIS_MATRIX_ILU L,
	LIS_MATRIX_ILU U,
	LIS_HALF_ILU *half,
	LIS_REAL *underflow_ratio);

extern LIS_INT lis_half_ilu_expand(
	LIS_HALF_ILU half,
	LIS_MATRIX_ILU *L,
	LIS_MATRIX_ILU *U);

extern LIS_INT lis_half_ilu_destroy(
	LIS_HALF_ILU half);


static int
set_row(
	LIS_MATRIX_ILU A,
	LIS_INT row,
	LIS_INT nnz,
	const LIS_INT *index,
	const LIS_SCALAR *value)
{
	LIS_INT j;

	if( nnz==0 )
	{
		return 0;
	}

	A->index[row] = (LIS_INT *)malloc(
		nnz*sizeof(LIS_INT));

	A->value[row] = (LIS_SCALAR *)malloc(
		nnz*sizeof(LIS_SCALAR));

	if( A->index[row]==NULL ||
		A->value[row]==NULL )
	{
		free(A->index[row]);
		free(A->value[row]);

		A->index[row] = NULL;
		A->value[row] = NULL;

		return 1;
	}

	A->nnz[row] = nnz;

	for(j=0;j<nnz;j++)
	{
		A->index[row][j] = index[j];
		A->value[row][j] = value[j];
	}

	return 0;
}


static int
check_structure(
	LIS_MATRIX_ILU A,
	LIS_MATRIX_ILU B)
{
	LIS_INT i,j;

	if( A->n!=B->n )
	{
		return 1;
	}

	for(i=0;i<A->n;i++)
	{
		if( A->nnz[i]!=B->nnz[i] )
		{
			return 1;
		}

		for(j=0;j<A->nnz[i];j++)
		{
			if( A->index[i][j]!=
				B->index[i][j] )
			{
				return 1;
			}
		}
	}

	return 0;
}


int
main(
	int argc,
	char **argv)
{
	LIS_MATRIX_ILU L,U,L2,U2;
	LIS_HALF_ILU half;
	LIS_REAL underflow;
	LIS_INT idx[2];
	LIS_SCALAR val[2];
	LIS_INT err;
	double rel;
	int failed;

	L = NULL;
	U = NULL;
	L2 = NULL;
	U2 = NULL;
	half = NULL;
	failed = 0;

	if( lis_initialize(&argc,&argv) )
	{
		return 1;
	}

	if( lis_matrix_ilu_create(4,1,&L) ||
		lis_matrix_ilu_setCR(L) ||
		lis_matrix_ilu_create(4,1,&U) ||
		lis_matrix_ilu_setCR(U) )
	{
		failed = 1;
		goto cleanup;
	}

	/* Extreme minimum binary64 value. */
	idx[0] = 0;
	val[0] = (LIS_SCALAR)nextafter(0.0,1.0);

	if( set_row(L,1,1,idx,val) )
	{
		failed = 1;
		goto cleanup;
	}

	/* Exactly representable after row scaling. */
	idx[0] = 0;
	idx[1] = 1;
	val[0] = 3.0;
	val[1] = -1.5;

	if( set_row(L,2,2,idx,val) )
	{
		failed = 1;
		goto cleanup;
	}

	/* Deliberate FP16 underflow of the small row entry. */
	idx[0] = 0;
	idx[1] = 2;
	val[0] = 1.0;
	val[1] = (LIS_SCALAR)ldexp(1.0,-50);

	if( set_row(L,3,2,idx,val) )
	{
		failed = 1;
		goto cleanup;
	}

	/* Very large finite values. */
	idx[0] = 1;
	idx[1] = 3;
	val[0] = (LIS_SCALAR)ldexp(1.0,1023);
	val[1] = (LIS_SCALAR)-ldexp(1.0,1022);

	if( set_row(U,0,2,idx,val) )
	{
		failed = 1;
		goto cleanup;
	}

	/* Ordinary non-exact quantization case. */
	idx[0] = 2;
	val[0] = (LIS_SCALAR)(1.0/3.0);

	if( set_row(U,1,1,idx,val) )
	{
		failed = 1;
		goto cleanup;
	}

	idx[0] = 3;
	val[0] = 0.25;

	if( set_row(U,2,1,idx,val) )
	{
		failed = 1;
		goto cleanup;
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
			"compress failed rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	if( !(underflow>0.0 &&
		  underflow<1.0e-6) )
	{
		fprintf(
			stderr,
			"unexpected underflow ratio %.17g\n",
			(double)underflow);

		failed = 1;
		goto cleanup;
	}

	err = lis_half_ilu_expand(
		half,
		&L2,
		&U2);

	if( err )
	{
		fprintf(
			stderr,
			"expand failed rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	if( check_structure(L,L2) ||
		check_structure(U,U2) )
	{
		fprintf(stderr,"structure mismatch\n");
		failed = 1;
		goto cleanup;
	}

	if( L2->value[1][0]!=L->value[1][0] ||
		L2->value[2][0]!=L->value[2][0] ||
		L2->value[2][1]!=L->value[2][1] ||
		L2->value[3][0]!=1.0 ||
		L2->value[3][1]!=0.0 )
	{
		fprintf(
			stderr,
			"L roundtrip mismatch\n");

		failed = 1;
		goto cleanup;
	}

	if( U2->value[0][0]!=U->value[0][0] ||
		U2->value[0][1]!=U->value[0][1] ||
		U2->value[2][0]!=U->value[2][0] )
	{
		fprintf(
			stderr,
			"U exact roundtrip mismatch\n");

		failed = 1;
		goto cleanup;
	}

	rel = fabs(
		(double)U2->value[1][0] -
		(1.0/3.0)) / (1.0/3.0);

	if( rel>5.0e-4 )
	{
		fprintf(
			stderr,
			"quantization error too large %.17g\n",
			rel);

		failed = 1;
		goto cleanup;
	}

	/* NaN factors must always be rejected. */
	U->value[1][0] = (LIS_SCALAR)NAN;

	lis_half_ilu_destroy(half);
	half = NULL;

	err = lis_half_ilu_compress(
		L,
		U,
		&half,
		NULL);

	if( err!=LIS_ERR_ILL_ARG ||
		half!=NULL )
	{
		fprintf(
			stderr,
			"NaN safety rejection failed rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	printf(
		"HALF_STORAGE_C2 PASSED "
		"underflow=%.17g\n",
		(double)underflow);


cleanup:

	lis_half_ilu_destroy(half);

	if( L2!=NULL )
	{
		lis_matrix_ilu_destroy(L2);
	}

	if( U2!=NULL )
	{
		lis_matrix_ilu_destroy(U2);
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
