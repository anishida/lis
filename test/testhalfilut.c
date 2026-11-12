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


static LIS_INT
create_matrix(
	LIS_MATRIX *A)
{
	const LIS_INT n = 64;

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
		if( i>1 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i-2,
					(LIS_SCALAR)-0.08,
					*A);

			if( err ) return err;
		}

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
				(4.5 +
				 0.01*(double)i),
				*A);

		if( err ) return err;

		if( i<n-1 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i+1,
					(LIS_SCALAR)-0.53,
					*A);

			if( err ) return err;
		}

		if( i<n-2 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i+2,
					(LIS_SCALAR)-0.06,
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


static double
relative_difference(
	LIS_VECTOR A,
	const LIS_SCALAR *ref,
	LIS_INT n)
{
	LIS_INT i;

	double nr;
	double nd;
	double d;

	nr = 0.0;
	nd = 0.0;

	for(i=0;i<n;i++)
	{
		nr =
			hypot(
				nr,
				(double)ref[i]);

		d =
			(double)A->value[i] -
			(double)ref[i];

		nd =
			hypot(
				nd,
				d);
	}

	if( nr==0.0 )
	{
		return nd;
	}

	return nd/nr;
}


static int
run_mode(
	LIS_MATRIX A,
	LIS_INT mode,
	int strict,
	LIS_SCALAR *ref,
	LIS_SCALAR *refh)
{
	LIS_SOLVER solver;
	LIS_PRECON precon;

	LIS_VECTOR B;
	LIS_VECTOR X;
	LIS_VECTOR XH;

	LIS_INT i;
	LIS_INT err;

	double rel;
	double relh;

	char options[256];

	int failed;

	solver = NULL;
	precon = NULL;

	B = NULL;
	X = NULL;
	XH = NULL;

	failed = 0;

	err =
		lis_solver_create(
			&solver);

	if( err )
	{
		return 1;
	}

	snprintf(
		options,
		sizeof(options),
		"-i gmres "
		"-p ilut "
		"-scale none "
		"-print none "
		"-precon_half %d",
		(int)mode);

	err =
		lis_solver_set_option(
			options,
			solver);

	if( err )
	{
		fprintf(
			stderr,
			"option parse mode=%d rc=%d\n",
			(int)mode,
			(int)err);

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
		lis_precon_create(
			solver,
			&precon);

	if( err )
	{
		fprintf(
			stderr,
			"ILUT create mode=%d rc=%d\n",
			(int)mode,
			(int)err);

		precon = NULL;
		failed = 1;
		goto cleanup;
	}

	solver->precon = precon;

	if( mode==LIS_PRECON_HALF_OFF ||
		strict )
	{
		if( precon->half!=NULL )
		{
			fprintf(
				stderr,
				"unexpected half mode=%d strict=%d\n",
				(int)mode,
				strict);

			failed = 1;
			goto cleanup;
		}

		if( !normal_values_present(precon) )
		{
			fprintf(
				stderr,
				"normal factors missing mode=%d strict=%d\n",
				(int)mode,
				strict);

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
				"half factors missing mode=%d\n",
				(int)mode);

			failed = 1;
			goto cleanup;
		}

#ifdef USE_QUAD_PRECISION

		if( !normal_values_present(precon) )
		{
			fprintf(
				stderr,
				"quad normal factors missing\n");

			failed = 1;
			goto cleanup;
		}

#else

		if( !normal_values_released(precon) )
		{
			fprintf(
				stderr,
				"normal factor values not released\n");

			failed = 1;
			goto cleanup;
		}

#endif
	}

	err = lis_vector_duplicate(A,&B);
	if( err ) goto vector_fail;

	err = lis_vector_duplicate(A,&X);
	if( err ) goto vector_fail;

	err = lis_vector_duplicate(A,&XH);
	if( err ) goto vector_fail;

	for(i=0;i<A->n;i++)
	{
		B->value[i] =
			(LIS_SCALAR)
			(0.2 +
			 0.013*(double)i);
	}

	err =
		lis_psolve_ilut_csr(
			solver,
			B,
			X);

	if( err )
	{
		fprintf(
			stderr,
			"psolve mode=%d rc=%d\n",
			(int)mode,
			(int)err);

		failed = 1;
		goto cleanup;
	}

	err =
		lis_psolveh_ilut_csr(
			solver,
			B,
			XH);

	if( err )
	{
		fprintf(
			stderr,
			"psolveh mode=%d rc=%d\n",
			(int)mode,
			(int)err);

		failed = 1;
		goto cleanup;
	}

	if( mode==LIS_PRECON_HALF_OFF )
	{
		for(i=0;i<A->n;i++)
		{
			ref[i] =
				X->value[i];

			refh[i] =
				XH->value[i];
		}

		rel = 0.0;
		relh = 0.0;
	}
	else
	{
		rel =
			relative_difference(
				X,
				ref,
				A->n);

		relh =
			relative_difference(
				XH,
				refh,
				A->n);

		if( strict )
		{
			if( rel>1.0e-14 ||
				relh>1.0e-14 )
			{
				fprintf(
					stderr,
					"fallback mismatch "
					"rel=%.17g relh=%.17g\n",
					rel,
					relh);

				failed = 1;
				goto cleanup;
			}
		}
		else if( rel>1.0e-3 ||
			relh>1.0e-3 )
		{
			fprintf(
				stderr,
				"half mismatch mode=%d "
				"rel=%.17g relh=%.17g\n",
				(int)mode,
				rel,
				relh);

			failed = 1;
			goto cleanup;
		}
	}

	printf(
		"HALF_ILUT_STAGE3 PASSED "
		"mode=%d strict=%d "
		"rel=%.17g relh=%.17g "
		"half=%d",
		(int)mode,
		strict,
		rel,
		relh,
		precon->half!=NULL);

#ifdef USE_QUAD_PRECISION

	printf(
		" normal=retained quad=1\n");

#else

	if( precon->half!=NULL )
	{
		printf(
			" normal=released quad=0\n");
	}
	else
	{
		printf(
			" normal=retained quad=0\n");
	}

#endif

	goto cleanup;


vector_fail:

	fprintf(
		stderr,
		"vector allocation mode=%d rc=%d\n",
		(int)mode,
		(int)err);

	failed = 1;


cleanup:

	if( XH!=NULL )
		lis_vector_destroy(XH);

	if( X!=NULL )
		lis_vector_destroy(X);

	if( B!=NULL )
		lis_vector_destroy(B);

	if( precon!=NULL )
	{
		solver->precon = NULL;
		lis_precon_destroy(precon);
	}

	if( solver!=NULL )
	{
		lis_solver_destroy(solver);
	}

	return failed;
}



/* Exercise the opt-in symmetric-graph RCM and row-equilibrated ILUT path.
 * The old ILUT path remains unchanged when neither option is set.
 */
static int
check_ilut_ordering(LIS_MATRIX A)
{
    const char *variants[] = {
        "-ilut_rcm 1 -ilut_row_scale 0",
        "-ilut_rcm 0 -ilut_row_scale 1",
        "-ilut_rcm 1 -ilut_row_scale 1"
    };
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL,x=NULL,h=NULL,r=NULL,e=NULL;
    LIS_INT i,j,err,seen[64];
    LIS_REAL rn=0.0,bn=0.0;
    int v,failed=0;
    char opts[320];

    if(A->n!=64 || A->nprocs!=1) return 0;
    if(lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&h) ||
       lis_vector_duplicate(A,&r) ||
       lis_vector_duplicate(A,&e))
    { failed=1; goto cleanup; }
    lis_vector_set_all((LIS_SCALAR)1.0,e);
    err=lis_matvec(A,e,b);
    if(err) { failed=1; goto cleanup; }
    lis_vector_nrm2(b,&bn);

    for(v=0;v<3;v++)
    {
        solver=NULL;pc=NULL;
        err=lis_solver_create(&solver);
        if(err) { failed=1; goto cleanup; }
        snprintf(opts,sizeof(opts),
                 "-i gmres -p ilut -scale none -print none "
                 "-ilut_drop 0 -ilut_rate 10 -ilu_pivot_tol 1e-12 %s",
                 variants[v]);
        err=lis_solver_set_option(opts,solver);
        if(err) { failed=1; goto cleanup; }
        err=lis_solver_set_matrix(A,solver);
        if(err) { failed=1; goto cleanup; }
        err=lis_precon_create(solver,&pc);
        if(err) { failed=1; goto cleanup; }
        solver->precon=pc;
        if(!pc->ilut_perm || !pc->ilut_row_scale)
        { failed=1; goto cleanup; }

        memset(seen,0,sizeof(seen));
        for(i=0;i<64;i++)
        {
            j=pc->ilut_perm[i];
            if(j<0 || j>=64 || seen[j]++ ||
               !isfinite((double)pc->ilut_row_scale[i]) ||
               pc->ilut_row_scale[i]<=0.0)
            { failed=1; goto cleanup; }
        }

        lis_vector_set_all((LIS_SCALAR)0.0,x);
        err=lis_psolve(solver,b,x);
        if(err) { failed=1; goto cleanup; }
        err=lis_matvec(A,x,r);
        if(err) { failed=1; goto cleanup; }
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) { failed=1; goto cleanup; }
        lis_vector_nrm2(r,&rn);
        if(!(rn/bn<0.10))
        {
            fprintf(stderr,"ILUT RCM variant %d poor residual %.3g\n",
                    v,(double)(rn/bn));
            failed=1; goto cleanup;
        }

        /* Check transpose mapping independently from the normal solve. */
        lis_vector_set_all((LIS_SCALAR)0.0,h);
        err=lis_psolveh(solver,b,h);
        if(err) { failed=1; goto cleanup; }
        err=lis_matvech(A,h,r);
        if(err) { failed=1; goto cleanup; }
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) { failed=1; goto cleanup; }
        lis_vector_nrm2(r,&rn);
        if(!(rn/bn<0.10))
        {
            fprintf(stderr,"ILUT RCM adjoint variant %d residual %.3g\n",
                    v,(double)(rn/bn));
            failed=1; goto cleanup;
        }
        solver->precon=NULL;
        lis_precon_destroy(pc);pc=NULL;
        lis_solver_destroy(solver);solver=NULL;
    }

    printf("ILUT_RCM_ROW_SCALE PASSED\n");
cleanup:
    if(pc)
    {
        if(solver) solver->precon=NULL;
        lis_precon_destroy(pc);
    }
    if(solver) lis_solver_destroy(solver);
    if(e) lis_vector_destroy(e);
    if(r) lis_vector_destroy(r);
    if(h) lis_vector_destroy(h);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    return failed;
}

int
main(
	int argc,
	char **argv)
{
	LIS_MATRIX A;

	LIS_SCALAR *ref;
	LIS_SCALAR *refh;

	LIS_INT err;

	int strict;
	int failed;

	A = NULL;
	ref = NULL;
	refh = NULL;

	failed = 0;

	strict =
		argc>1 &&
		strcmp(
			argv[1],
			"fallback")==0;

	err =
		lis_initialize(
			&argc,
			&argv);

	if( err )
	{
		return 1;
	}

	err = create_matrix(&A);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	ref =
		(LIS_SCALAR *)malloc(
			A->n*sizeof(LIS_SCALAR));

	refh =
		(LIS_SCALAR *)malloc(
			A->n*sizeof(LIS_SCALAR));

	if( ref==NULL ||
		refh==NULL )
	{
		failed = 1;
		goto cleanup;
	}

	failed |= check_ilut_ordering(A);

	failed |=
		run_mode(
			A,
			LIS_PRECON_HALF_OFF,
			0,
			ref,
			refh);

	if( strict )
	{
		failed |=
			run_mode(
				A,
				LIS_PRECON_HALF_ON,
				1,
				ref,
				refh);
	}
	else
	{
		failed |=
			run_mode(
				A,
				LIS_PRECON_HALF_ON,
				0,
				ref,
				refh);

		failed |=
			run_mode(
				A,
				LIS_PRECON_HALF_NOPROBE,
				0,
				ref,
				refh);
	}

	if( !failed )
	{
		printf(
			"HALF_ILUT_STAGE3 COMPLETE PASS\n");
	}


cleanup:

	free(ref);
	free(refh);

	if( A!=NULL )
		lis_matrix_destroy(A);

	lis_finalize();

	return failed ? 1 : 0;
}
