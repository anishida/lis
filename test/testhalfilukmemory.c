#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "lislib.h"


static LIS_INT
create_matrix(
	LIS_MATRIX *A)
{
	const LIS_INT n = 32;

	LIS_INT i;
	LIS_INT err;

	err =
		lis_matrix_create(
			LIS_COMM_WORLD,
			A);

	if( err ) return err;

	err =
		lis_matrix_set_size(
			*A,
			0,
			n);

	if( err ) return err;

	for(i=0;i<n;i++)
	{
		if( i>0 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i-1,
					(LIS_SCALAR)-0.71,
					*A);

			if( err ) return err;
		}

		err =
			lis_matrix_set_value(
				LIS_INS_VALUE,
				i,
				i,
				(LIS_SCALAR)
				(4.0 +
				 0.013*(double)i),
				*A);

		if( err ) return err;

		if( i<n-1 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i+1,
					(LIS_SCALAR)-0.63,
					*A);

			if( err ) return err;
		}
	}

	err =
		lis_matrix_set_type(
			*A,
			LIS_MATRIX_CSR);

	if( err ) return err;

	return lis_matrix_assemble(*A);
}


static void
change_values(
	LIS_MATRIX A)
{
	LIS_INT i,j,col;

	for(i=0;i<A->n;i++)
	{
		for(j=A->ptr[i];
		    j<A->ptr[i+1];
		    j++)
		{
			col = A->index[j];

			if( col==i )
			{
				A->value[j] =
					(LIS_SCALAR)
					(5.25 +
					 0.017*(double)i);
			}
			else
			{
				A->value[j] =
					(LIS_SCALAR)-0.29;
			}
		}
	}
}


static int
normal_values_present(
	LIS_PRECON precon)
{
	LIS_INT i;

	for(i=0;i<precon->L->n;i++)
	{
		if( precon->L->nnz[i]>0 &&
			precon->L->value[i]==NULL )
		{
			return 0;
		}

		if( precon->U->nnz[i]>0 &&
			precon->U->value[i]==NULL )
		{
			return 0;
		}
	}

	return 1;
}


static int
normal_values_released(
	LIS_PRECON precon)
{
	LIS_INT i;

	for(i=0;i<precon->L->n;i++)
	{
		if( precon->L->nnz[i]>0 &&
			precon->L->value[i]!=NULL )
		{
			return 0;
		}

		if( precon->U->nnz[i]>0 &&
			precon->U->value[i]!=NULL )
		{
			return 0;
		}
	}

	return 1;
}


int
main(
	int argc,
	char **argv)
{
	LIS_MATRIX A;
	LIS_SOLVER solver;
	LIS_PRECON precon;

	LIS_VECTOR B;
	LIS_VECTOR X1;
	LIS_VECTOR X2;

	LIS_INT i;
	LIS_INT err;

	double nr;
	double nd;
	double delta;

	int fallback;
	int failed;

	A = NULL;
	solver = NULL;
	precon = NULL;

	B = NULL;
	X1 = NULL;
	X2 = NULL;

	failed = 0;

	fallback =
		argc>1 &&
		strcmp(
			argv[1],
			"fallback")==0;

	if( lis_initialize(
			&argc,
			&argv) )
	{
		return 1;
	}

	err = create_matrix(&A);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err =
		lis_solver_create(
			&solver);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err =
		lis_solver_set_option(
			"-i gmres "
			"-p ilu "
			"-ilu_fill 0 "
			"-scale none "
			"-print none "
			"-precon_half 2",
			solver);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err =
		lis_solver_set_matrix(
			A,
			solver);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err =
		lis_precon_psd_create(
			solver,
			&precon);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	solver->precon = precon;

	err = lis_vector_duplicate(A,&B);
	if( err ) goto vector_fail;

	err = lis_vector_duplicate(A,&X1);
	if( err ) goto vector_fail;

	err = lis_vector_duplicate(A,&X2);
	if( err ) goto vector_fail;

	for(i=0;i<A->n;i++)
	{
		B->value[i] =
			(LIS_SCALAR)
			(0.31 +
			 0.027*(double)i);
	}

	err =
		lis_precon_psd_update(
			solver,
			precon);

	if( err )
	{
		fprintf(
			stderr,
			"first update rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	if( precon->half==NULL )
	{
		fprintf(
			stderr,
			"first half factors missing\n");

		failed = 1;
		goto cleanup;
	}

#ifdef USE_QUAD_PRECISION

	if( !normal_values_present(precon) )
	{
		fprintf(
			stderr,
			"quad build lost normal factors\n");

		failed = 1;
		goto cleanup;
	}

#else

	if( !normal_values_released(precon) )
	{
		fprintf(
			stderr,
			"normal values were not released\n");

		failed = 1;
		goto cleanup;
	}

#endif

	err =
		lis_psolve_iluk_csr(
			solver,
			B,
			X1);

	if( err )
	{
		fprintf(
			stderr,
			"first apply rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	change_values(A);

	if( fallback )
	{
		solver->precon_half_mode =
			LIS_PRECON_HALF_ON;
	}

	err =
		lis_precon_psd_update(
			solver,
			precon);

	if( err )
	{
		fprintf(
			stderr,
			"second update rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	if( fallback )
	{
		if( precon->half!=NULL )
		{
			fprintf(
				stderr,
				"strict refresh unexpectedly kept half\n");

			failed = 1;
			goto cleanup;
		}

		if( !normal_values_present(precon) )
		{
			fprintf(
				stderr,
				"fresh fallback factors missing\n");

			failed = 1;
			goto cleanup;
		}
	}
	else
	{
		if( precon->half==NULL )
		{
			fprintf(
				stderr,
				"refreshed half factors missing\n");

			failed = 1;
			goto cleanup;
		}

#ifdef USE_QUAD_PRECISION

		if( !normal_values_present(precon) )
		{
			fprintf(
				stderr,
				"quad refresh lost normal factors\n");

			failed = 1;
			goto cleanup;
		}

#else

		if( !normal_values_released(precon) )
		{
			fprintf(
				stderr,
				"refreshed normal values remain allocated\n");

			failed = 1;
			goto cleanup;
		}

#endif
	}

	err =
		lis_psolve_iluk_csr(
			solver,
			B,
			X2);

	if( err )
	{
		fprintf(
			stderr,
			"second apply rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	nr = 0.0;
	nd = 0.0;

	for(i=0;i<A->n;i++)
	{
		double d;

		nr =
			hypot(
				nr,
				(double)X2->value[i]);

		d =
			(double)X2->value[i] -
			(double)X1->value[i];

		nd = hypot(nd,d);
	}

	delta =
		nr==0.0 ?
		nd :
		nd/nr;

	if( delta<1.0e-3 )
	{
		fprintf(
			stderr,
			"updated factors had no visible effect %.17g\n",
			delta);

		failed = 1;
		goto cleanup;
	}

	printf(
		"HALF_ILUK_MEMORY_STAGE2 PASSED "
		"fallback=%d delta=%.17g half=%d",
		fallback,
		delta,
		precon->half!=NULL);

#ifdef USE_QUAD_PRECISION

	printf(
		" normal=retained quad=1\n");

#else

	printf(
		" normal=%s quad=0\n",
		precon->half!=NULL ?
			"released" :
			"fallback");

#endif

	goto cleanup;


vector_fail:

	fprintf(
		stderr,
		"vector allocation rc=%d\n",
		(int)err);

	failed = 1;


cleanup:

	if( X2!=NULL )
		lis_vector_destroy(X2);

	if( X1!=NULL )
		lis_vector_destroy(X1);

	if( B!=NULL )
		lis_vector_destroy(B);

	if( precon!=NULL )
	{
		solver->precon = NULL;
		lis_precon_destroy(precon);
	}

	if( solver!=NULL )
		lis_solver_destroy(solver);

	if( A!=NULL )
		lis_matrix_destroy(A);

	lis_finalize();

	return failed ? 1 : 0;
}
