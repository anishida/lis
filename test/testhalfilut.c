#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <math.h>
#include <stdio.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
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
    LIS_INT i,j,err,seen[64],nnz_l=0,nnz_u=0;
    LIS_REAL rn=0.0,bn=0.0,min_piv=0.0,max_piv=0.0;
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
        err=lis_precon_ilut_factor_stats(pc,&nnz_l,&nnz_u,
                                         &min_piv,&max_piv);
        if(err || nnz_l<=0 || nnz_u<=0 ||
           !(min_piv>0.0 && max_piv>=min_piv))
        { failed=1; goto cleanup; }
        printf("ILUT_RCM_FACTOR_STATS variant=%d nnz_l=%d nnz_u=%d min_pivot=%.6g max_pivot=%.6g\n",
               v,(int)nnz_l,(int)nnz_u,(double)min_piv,(double)max_piv);
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


/*
 * Deterministic stabilized contact-KKT-like system:
 * [ K  B^T ][u] = [f]
 * [ B -0.1I][λ]   [g]
 *
 * K is an SPD spring chain, B enforces displacement differences.
 * The negative multiplier block makes the matrix indefinite.
 * This is an algebraic regression, NOT a full nonlinear contact FEA run.
 */
static int
check_ilut_contact_kkt(void)
{
    const LIS_INT nu=36,nlambda=12,n=48;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,u=NULL,r=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,i,q,iters=0,variant;
    LIS_REAL bn=0.0,rn=0.0;
    char opts[320];
    int failed=1;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,n);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<nu;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)4.0,A);
        if(err) goto done;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                                     (LIS_SCALAR)-0.8,A);
            if(err) goto done;
        }
        if(i+1<nu)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)-0.8,A);
            if(err) goto done;
        }
    }
    for(q=0;q<nlambda;q++)
    {
        LIS_INT l=nu+q,a=3*q,c=3*q+1;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,l,
                                 (LIS_SCALAR)-0.1,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,a,(LIS_SCALAR)1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,c,(LIS_SCALAR)-1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,a,l,(LIS_SCALAR)1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,c,l,(LIS_SCALAR)-1.0,A);
        if(err) goto done;
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    if(lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&u) ||
       lis_vector_duplicate(A,&r)) goto done;
    lis_vector_set_all((LIS_SCALAR)1.0,u);
    err=lis_matvec(A,u,b);
    if(err) goto done;
    err=lis_vector_nrm2(b,&bn);
    if(err || bn<=0.0) goto done;
    for(variant=0;variant<2;variant++)
    {
        solver=NULL;
        err=lis_solver_create(&solver);
        if(err) goto done;
        snprintf(opts,sizeof(opts),
            "-i fgmres -p ilut -scale none -print none "
            "-tol 1e-8 -maxiter 120 -restart 30 "
            "-ilut_drop 0.005 -ilut_rate 8 -ilu_pivot_tol 1e-9 %s",
            variant ? "-ilut_rcm 1 -ilut_row_scale 1" : "");
        err=lis_solver_set_option(opts,solver);
        if(err) goto done;
        lis_vector_set_all((LIS_SCALAR)0.0,x);
        err=lis_solve(A,b,x,solver);
        if(err) goto done;
        err=lis_matvec(A,x,r);
        if(err) goto done;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) goto done;
        err=lis_vector_nrm2(r,&rn);
        if(err) goto done;
        err=lis_solver_get_iter(solver,&iters);
        if(err || rn/bn>1e-7 || iters<1)
        {
            fprintf(stderr,"ILUT_CONTACT_KKT variant=%d iter=%d residual=%.3e\n",
                    variant,(int)iters,(double)(rn/bn));
            goto done;
        }
        printf("ILUT_CONTACT_KKT variant=%d iter=%d residual=%.6e\n",
               variant,(int)iters,(double)(rn/bn));
        lis_solver_destroy(solver); solver=NULL;
    }
    failed=0;
done:
    if(failed)
        fprintf(stderr,"ILUT_CONTACT_KKT_FAILED rc=%d\n",(int)err);
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(u) lis_vector_destroy(u);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    return failed;
}

/* Opt-in serial partitioned RAS ILUT with diagonal-approximate Schur
 * correction must remain separate from the legacy ILUT default.
 */
static int
check_ilut_ras_schur(LIS_MATRIX A)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL,x=NULL,e=NULL,r=NULL;
    LIS_INT err,overlap,schur;
    LIS_REAL bn=0.0,rn=0.0;
    int failed=0;
    char opts[300];

    if(A->nprocs!=1 || A->n!=64) return 0;
    if(lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&e) || lis_vector_duplicate(A,&r))
    { failed=1; goto cleanup; }
    lis_vector_set_all((LIS_SCALAR)1.0,e);
    err=lis_matvec(A,e,b);
    if(err) { failed=1; goto cleanup; }
    err=lis_vector_nrm2(b,&bn);
    if(err) { failed=1; goto cleanup; }

    for(overlap=0;overlap<=2;overlap++)
    for(schur=0;schur<=1;schur++)
    {
        solver=NULL;pc=NULL;
        err=lis_solver_create(&solver);
        if(err) { failed=1; goto cleanup; }
        snprintf(opts,sizeof(opts),
            "-i gmres -p ilut -scale none -print none "
            "-iluc_drop 0 -iluc_rate 8 -ilu_pivot_tol 1e-12 "
            "-ilut_ras_domains 4 -ilut_ras_overlap %d -ilut_schur %d",
            (int)overlap,(int)schur);
        err=lis_solver_set_option(opts,solver);
        if(err) { failed=1; goto cleanup; }
        err=lis_solver_set_matrix(A,solver);
        if(err) { failed=1; goto cleanup; }
        err=lis_precon_create(solver,&pc);
        if(err) { failed=1; goto cleanup; }
        solver->precon=pc;
        if(pc->ilut_ras==NULL)
        { failed=1; goto cleanup; }
        err=lis_psolve(solver,b,x);
        if(err) { failed=1; goto cleanup; }
        err=lis_matvec(A,x,r);
        if(err) { failed=1; goto cleanup; }
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) { failed=1; goto cleanup; }
        err=lis_vector_nrm2(r,&rn);
        if(err || !(rn/bn<0.80))
        {
            fprintf(stderr,"ILUT RAS overlap=%d schur=%d residual %.5e\n",
                    (int)overlap,(int)schur,(double)(rn/bn));
            failed=1; goto cleanup;
        }
        /* No mathematically incorrect adjoint implementation. */
        err=lis_psolveh(solver,b,x);
        if(err!=LIS_ERR_NOT_IMPLEMENTED)
        {
            fprintf(stderr,"ILUT RAS unexpectedly accepted adjoint\n");
            failed=1; goto cleanup;
        }
        solver->precon=NULL;
        lis_precon_destroy(pc);pc=NULL;
        lis_solver_destroy(solver);solver=NULL;
    }
    printf("ILUT_RAS_SCHUR PASSED\n");

cleanup:
    if(pc)
    {
        if(solver) solver->precon=NULL;
        lis_precon_destroy(pc);
    }
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(e) lis_vector_destroy(e);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    return failed;
}


/*
 * Exact-interface correctness witness, not a general performance claim.
 * A = 4I with symmetric -1 edges (2,3), (3,4), and (4,5).
 * Partition [0,4) / [4,8) or [0,4) / [4,9):
 * Gamma={3,4}, A_II=4I, even for unequal partition sizes.
 * exact S=[[3.75,-1],[-1,3.75]].
 * For b=A*(e3+e4), overlap-0 RAS alone has residual ~0.2319;
 * the additional Schur step should recover the solution to roundoff.
 */
static int
check_ilut_ras_exact_interface_one(LIS_INT n,LIS_INT coupled,LIS_INT local_ilutp)
{
    const LIS_INT ei[3]={2,3,4},ej[3]={3,4,5};
    LIS_MATRIX A=NULL;
    LIS_VECTOR exact=NULL,b=NULL,x=NULL,r=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT err=0,i,v,edge_count;
    LIS_REAL bn=0.0,rn=0.0,baseline=0.0,schur=0.0,xerr=0.0;
    char opts[256];
    int failed=1;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,n);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<n;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)4.0,A);
        if(err) goto done;
    }
    edge_count=coupled ? n-1 : 3;
    for(i=0;i<edge_count;i++)
    {
        LIS_INT u=coupled ? i : ei[i],v=coupled ? i+1 : ej[i];
        err=lis_matrix_set_value(LIS_INS_VALUE,u,v,(LIS_SCALAR)-1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,v,u,(LIS_SCALAR)-1.0,A);
        if(err) goto done;
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    if(lis_vector_duplicate(A,&exact) || lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&x) || lis_vector_duplicate(A,&r))
        goto done;
    err=lis_vector_set_all((LIS_SCALAR)0.0,exact);
    if(err) goto done;
    exact->value[3]=(LIS_SCALAR)1.0;
    exact->value[4]=(LIS_SCALAR)1.0;
    err=lis_matvec(A,exact,b);
    if(err) goto done;
    err=lis_vector_nrm2(b,&bn);
    if(err || !(bn>0.0)) goto done;

    for(v=0;v<2;v++)
    {
        err=lis_solver_create(&solver);
        if(err) goto done;
        snprintf(opts,sizeof(opts),
            "-i gmres -p ilut -scale none -print none "
            "-iluc_drop 0 -iluc_rate 8 -ilu_pivot_tol 1e-12 "
            "-ilut_ras_domains 2 -ilut_ras_overlap 0 -ilut_schur %d %s",
            (int)v,local_ilutp ?
            "-ilutp 1 -ilutp_threshold 0.5 -ilutp_full 1 "
            "-ilutp_full_limit 128" : "");
        err=lis_solver_set_option(opts,solver);
        if(err) goto done;
        err=lis_solver_set_matrix(A,solver);
        if(err) goto done;
        err=lis_precon_create(solver,&pc);
        if(err) { pc=NULL; goto done; }
        solver->precon=pc;
        err=lis_psolve(solver,b,x);
        if(err) goto done;
        err=lis_matvec(A,x,r);
        if(err) goto done;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) goto done;
        err=lis_vector_nrm2(r,&rn);
        if(err) goto done;
        if(v==0) baseline=rn/bn;
        else
        {
            schur=rn/bn;
            for(i=0;i<n;i++)
            {
                LIS_REAL delta=fabs(x->value[i]-exact->value[i]);
                if(delta>xerr) xerr=delta;
            }
        }
        printf("ILUT_RAS_INTERFACE n=%ld coupled=%d ilutp=%d "
               "schur=%d true_rel_res=%.9e\n",
               (long)n,(int)coupled,(int)local_ilutp,
               (int)v,(double)(rn/bn));
        solver->precon=NULL;
        lis_precon_destroy(pc);pc=NULL;
        lis_solver_destroy(solver);solver=NULL;
    }
    /* Exact diagonal interior or >=5x improvement with coupled interior. */
    if(!(baseline>0.20 && baseline<0.26 &&
         (coupled ? (schur>0.001 && schur<0.05 &&
                     schur<baseline/5.0) :
                    (schur<1e-10 && xerr<1e-10))))
    {
        fprintf(stderr,"ILUT_RAS_INTERFACE_FAILED n=%ld coupled=%d baseline=%.9e "
                "schur=%.9e xerr=%.9e\n",
                (long)n,(int)coupled,(double)baseline,(double)schur,(double)xerr);
        goto done;
    }
    printf("ILUT_RAS_INTERFACE n=%ld coupled=%d ilutp=%d PASSED\n",
           (long)n,(int)coupled,(int)local_ilutp);
    failed=0;
done:
    if(pc)
    {
        if(solver) solver->precon=NULL;
        lis_precon_destroy(pc);
    }
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(exact) lis_vector_destroy(exact);
    if(A) lis_matrix_destroy(A);
    return failed;
}

static int
check_ilut_ras_exact_interface(void)
{
#ifdef _OPENMP
    int original_threads;
#endif
    int failed=check_ilut_ras_exact_interface_one(8,0,0);
    failed |= check_ilut_ras_exact_interface_one(9,0,0);
    failed |= check_ilut_ras_exact_interface_one(8,1,0);
    failed |= check_ilut_ras_exact_interface_one(9,1,0);
#ifdef _OPENMP
    original_threads=omp_get_max_threads();
    /* Native ILUTP (including nested RAS factors) is single-thread only. */
    omp_set_num_threads(1);
#endif
    /* Schur on/off with genuinely coupled interior and local ILUTP. */
    failed |= check_ilut_ras_exact_interface_one(8,1,1);
    failed |= check_ilut_ras_exact_interface_one(9,1,1);
#ifdef _OPENMP
    omp_set_num_threads(original_threads);
#endif
    return failed;
}

/*
 * Integration regression for #225 RAS with #229 ILUTP.
 * Sixteen reversible 2x2 off-diagonal blocks have zero diagonals.
 * Eight independent four-row domains exercise up to eight OpenMP workers.
 * Dynamic column pivoting is required for exact local complete-fill LU.
 * Without propagating -ilutp from the outer solver, this deliberately
 * difficult local factorization cannot pass the true residual test.
 */
static int check_ilut_ras_local_ilutp(void)
{
#ifdef _OPENMP
    int original_threads=omp_get_max_threads();
    int original_dynamic=omp_get_dynamic();
#endif
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL,x=NULL,exact=NULL,r=NULL;
    LIS_INT err=0,i,variant,limit;
    int expected_threads;
    LIS_REAL bn=0.0,rn=0.0,exerr=0.0;
    int failed=1;
    char opts[350];

#ifdef _OPENMP
    /* Test serial and genuinely parallel independent ILUTP subdomains. */
    omp_set_dynamic(0);
    omp_set_num_threads(1);
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto cleanup;
    err=lis_matrix_set_size(A,0,32);
    if(err) goto cleanup;
    if(A->nprocs!=1) { failed=0; goto cleanup; }
    for(i=0;i<32;i+=2)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,(LIS_SCALAR)1.0,A);
        if(err) goto cleanup;
        err=lis_matrix_set_value(LIS_INS_VALUE,i+1,i,(LIS_SCALAR)1.0,A);
        if(err) goto cleanup;
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto cleanup;
    err=lis_matrix_assemble(A);
    if(err) goto cleanup;

    if(lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&exact) || lis_vector_duplicate(A,&r))
        goto cleanup;
    for(i=0;i<32;i++) exact->value[i]=(LIS_SCALAR)(1.0+0.125*i);
    err=lis_matvec(A,exact,b);
    if(err) goto cleanup;
    err=lis_vector_nrm2(b,&bn);
    if(err || !(bn>0.0)) goto cleanup;
#ifdef _OPENMP
    limit=4;
#else
    limit=1; /* Never report 2/4/8 threads when OpenMP is disabled. */
#endif
    for(variant=0;variant<=limit;variant++)
    {
#ifdef _OPENMP
        /* Positive variants exercise 1, 2, 4 and 8 OpenMP workers.
         * Eight domains of four rows give eight independent pivoted LU
         * factorizations with a zero diagonal in every 2x2 block. */
        expected_threads=variant ? (1 << (variant-1)) : 1;
        omp_set_num_threads(expected_threads);
#endif
#ifndef _OPENMP
        expected_threads=1;
#endif
        exerr=0.0;
        err=lis_solver_create(&solver);
        if(err) goto cleanup;
        snprintf(opts,sizeof(opts),
            "-i fgmres -p ilut -scale none -print none "
            "-ilut_drop 0 -ilut_rate 8 "
            "-ilut_ras_domains 8 -ilut_ras_overlap 0 -ilut_schur 0 "
            "-ilutp %d -ilutp_threshold 0.5 -ilutp_full 1 "
            "-ilutp_full_limit 128",
            (int)(variant!=0));
        err=lis_solver_set_option(opts,solver);
        if(err) goto cleanup;
        err=lis_solver_set_matrix(A,solver);
        if(err) goto cleanup;
        err=lis_precon_create(solver,&pc);
        if(variant==0)
        {
            /* Full fill may never be accepted with ILUTP disabled. */
            if(err!=LIS_ERR_ILL_ARG)
            {
                fprintf(stderr,"RAS+ILUTP silently ignored -ilutp_full\n");
                goto cleanup;
            }
            pc=NULL;
            lis_solver_destroy(solver);solver=NULL;
            continue;
        }
        if(err) { pc=NULL; goto cleanup; }
        solver->precon=pc;
        err=lis_psolve(solver,b,x);
        if(err) goto cleanup;
        err=lis_matvec(A,x,r);
        if(err) goto cleanup;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) goto cleanup;
        err=lis_vector_nrm2(r,&rn);
        if(err) goto cleanup;
        for(i=0;i<32;i++)
        {
            LIS_REAL v=fabs(x->value[i]-exact->value[i]);
            if(v>exerr) exerr=v;
        }
        if(!(rn/bn<1.0e-10 && exerr<1.0e-10))
        {
            fprintf(stderr,
                "RAS+ILUTP pivoted local solve residual=%.6e xerr=%.6e\n",
                (double)(rn/bn),(double)exerr);
            goto cleanup;
        }
        printf("ILUT_RAS_LOCAL_ILUTP threads=%d residual=%.9e xerr=%.9e PASSED\n",
               expected_threads,
               (double)(rn/bn),(double)exerr);
#ifdef _OPENMP
        if(omp_get_max_threads()!=expected_threads)
        {
            fprintf(stderr,"ILUT RAS altered outer OpenMP thread budget\n");
            goto cleanup;
        }
#endif
        solver->precon=NULL;
        lis_precon_destroy(pc);pc=NULL;
        lis_solver_destroy(solver);solver=NULL;
    }
    failed=0;
cleanup:
    if(pc)
    {
        if(solver) solver->precon=NULL;
        lis_precon_destroy(pc);
    }
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(exact) lis_vector_destroy(exact);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
#ifdef _OPENMP
    omp_set_num_threads(original_threads);
    omp_set_dynamic(original_dynamic);
#endif
    return failed;
}

/*
 * Real cross-domain coupling regression for parallel local ILUTP:
 * nonzero links connect the eight contiguous four-row RAS cores.
 * With overlap 1/2, local subdomains share imported RHS rows while
 * correction writes remain restricted to disjoint owner cores.
 * Test both complete-fill and dropped ILUTP with RCM/row scaling, and
 * demand bitwise deterministic RAS apply across 1/2/4/8 workers.
 */
static int check_ilut_ras_parallel_overlap(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL,x=NULL,exact=NULL,r=NULL;
    LIS_SCALAR reference[32];
    LIS_INT err=0,i,overlap,full,trial,repeat;
    LIS_REAL bn=0.0,rn=0.0;
    int failed=1,expected_threads,limit;
    char opts[384];
#ifdef _OPENMP
    int original_threads=omp_get_max_threads();
    int original_dynamic=omp_get_dynamic();
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    limit=3;
#else
    limit=0;
#endif

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto cleanup;
    err=lis_matrix_set_size(A,0,32);
    if(err) goto cleanup;
    if(A->nprocs!=1) { failed=0; goto cleanup; }
    for(i=0;i<32;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)0.2,A);
        if(err) goto cleanup;
        if((i%2)==0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)1.0,A);
            if(err) goto cleanup;
            err=lis_matrix_set_value(LIS_INS_VALUE,i+1,i,
                                     (LIS_SCALAR)1.0,A);
            if(err) goto cleanup;
        }
        /* Deliberately nonsymmetric couplings across owner boundaries. */
        if((i%4)==3 && i+1<32)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)0.04,A);
            if(err) goto cleanup;
            err=lis_matrix_set_value(LIS_INS_VALUE,i+1,i,
                                     (LIS_SCALAR)0.025,A);
            if(err) goto cleanup;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto cleanup;
    err=lis_matrix_assemble(A);
    if(err) goto cleanup;
    if(lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&exact) || lis_vector_duplicate(A,&r))
        goto cleanup;
    for(i=0;i<32;i++) exact->value[i]=(LIS_SCALAR)(0.7+0.03125*i);
    err=lis_matvec(A,exact,b);
    if(err) goto cleanup;
    err=lis_vector_nrm2(b,&bn);
    if(err || !(bn>0.0)) goto cleanup;

    for(overlap=1;overlap<=2;overlap++)
    for(full=0;full<=1;full++)
    for(trial=0;trial<=limit;trial++)
    {
        expected_threads=1<<trial;
#ifdef _OPENMP
        omp_set_num_threads(expected_threads);
#endif
        err=lis_solver_create(&solver);
        if(err) goto cleanup;
        snprintf(opts,sizeof(opts),
            "-i fgmres -p ilut -scale none -print none "
            "-ilut_drop 0.001 -ilut_rate 24 -ilu_pivot_tol 1e-9 "
            "-ilut_ras_domains 8 -ilut_ras_overlap %d -ilut_schur 0 "
            "-ilut_rcm 1 -ilut_row_scale 1 -ilutp 1 "
            "-ilutp_threshold 0.5 -ilutp_full %d -ilutp_full_limit 1024",
            (int)overlap,(int)full);
        err=lis_solver_set_option(opts,solver);
        if(err) goto cleanup;
        err=lis_solver_set_matrix(A,solver);
        if(err) goto cleanup;
        err=lis_precon_create(solver,&pc);
        if(err) { pc=NULL; goto cleanup; }
        solver->precon=pc;
        for(repeat=0;repeat<3;repeat++)
        {
            err=lis_psolve(solver,b,x);
            if(err) goto cleanup;
            err=lis_matvec(A,x,r);
            if(err) goto cleanup;
            err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
            if(err) goto cleanup;
            err=lis_vector_nrm2(r,&rn);
            if(err || !(rn>=0.0 && rn/bn<0.30))
            {
                fprintf(stderr,
                    "RAS_OVERLAP_BAD_RESIDUAL overlap=%d full=%d "
                    "threads=%d repeat=%d relres=%.8e\n",
                    (int)overlap,(int)full,expected_threads,
                    (int)repeat,(double)(rn/bn));
                goto cleanup;
            }
            if(trial==0 && repeat==0)
                memcpy(reference,x->value,sizeof(reference));
            else if(memcmp(reference,x->value,sizeof(reference))!=0)
            {
                fprintf(stderr,
                    "RAS_OVERLAP_NONDETERMINISTIC overlap=%d full=%d "
                    "threads=%d repeat=%d\n",
                    (int)overlap,(int)full,expected_threads,(int)repeat);
                goto cleanup;
            }
        }
#ifdef _OPENMP
        if(omp_get_max_threads()!=expected_threads)
        {
            fprintf(stderr,"RAS_OVERLAP leaked OpenMP thread setting\n");
            goto cleanup;
        }
#endif
        printf("ILUT_RAS_OVERLAP overlap=%d full=%d threads=%d "
               "residual=%.9e DETERMINISTIC PASS\n",
               (int)overlap,(int)full,expected_threads,(double)(rn/bn));
        solver->precon=NULL;
        lis_precon_destroy(pc);pc=NULL;
        lis_solver_destroy(solver);solver=NULL;
    }
    failed=0;

cleanup:
    if(pc)
    {
        if(solver) solver->precon=NULL;
        lis_precon_destroy(pc);
    }
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(exact) lis_vector_destroy(exact);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
#ifdef _OPENMP
    omp_set_num_threads(original_threads);
    omp_set_dynamic(original_dynamic);
#endif
    return failed;
}


/*
 * An RCO-to-CSR assembly owns an intermediate duplicated matrix. Force
 * OOM at every allocation in that conversion and verify no matrix leaks.
 */
static LIS_INT make_unassembled_rco_oom(LIS_MATRIX *out)
{
    LIS_INT i,err;
    *out=NULL;
    err=lis_matrix_create(LIS_COMM_WORLD,out);
    if(err)return err;
    err=lis_matrix_set_size(*out,0,8);
    if(err)return err;
    err=lis_matrix_set_type(*out,LIS_MATRIX_CSR);
    if(err)return err;
    for(i=0;i<8;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,*out);
        if(err)return err;
        if(i+1<8)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,*out);
            if(err)return err;
        }
    }
    return LIS_SUCCESS;
}

static int check_rco_assemble_oom_cleanup(void)
{
#ifndef LIS_TEST_OOM
    return 0;
#else
    LIS_MATRIX A=NULL;
    LIS_INT err;
    size_t baseline,calls,target;
    int failed=1;
    baseline=lis_test_oom_live();
    err=make_unassembled_rco_oom(&A);
    if(err)goto done;
    lis_test_oom_begin(0);
    err=lis_matrix_assemble(A);
    lis_test_oom_end();
    calls=lis_test_oom_calls();
    lis_matrix_destroy(A);A=NULL;
    if(err || calls<4 || lis_test_oom_live()!=baseline)
        goto done;
    for(target=1;target<=calls;target++)
    {
        err=make_unassembled_rco_oom(&A);
        if(err)goto done;
        lis_test_oom_begin(target);
        err=lis_matrix_assemble(A);
        lis_test_oom_end();
        lis_matrix_destroy(A);A=NULL;
        if(err!=LIS_OUT_OF_MEMORY || lis_test_oom_live()!=baseline)
        {
            fprintf(stderr,"RCO_ASSEMBLE_OOM_FAIL at %lu error=%ld live=%lu baseline=%lu\n",
                (unsigned long)target,(long)err,
                (unsigned long)lis_test_oom_live(),(unsigned long)baseline);
            goto done;
        }
    }
    printf("RCO_ASSEMBLE_OOM_CLEANUP checked=%lu PASS\n",(unsigned long)calls);
    failed=0;
done:
    lis_test_oom_end();
    if(A)lis_matrix_destroy(A);
    return failed;
#endif
}

/*
 * Regression for option-tokenizer OOM: lis_text2args must clean up a
 * partially constructed linked list and lis_solver_set_option must
 * propagate the error, not dereference a null list.
 */
static int check_ilutp_options_oom(void)
{
#ifndef LIS_TEST_OOM
    return 0;
#else
    LIS_SOLVER solver=NULL;
    LIS_INT err=0;
    size_t baseline,calls,fail_at,hit;
    int failed=1;
    const char *options=
        "-i fgmres -p ilut -ilutp 1 -ilutp_threshold 0.5 "
        "-ilut_ras_domains 8 -ilut_ras_overlap 1 -scale none";

    err=lis_solver_create(&solver);
    if(err) goto done;
    baseline=lis_test_oom_live();
    lis_test_oom_begin(0);
    err=lis_solver_set_option((char *)options,solver);
    lis_test_oom_end();
    calls=lis_test_oom_calls();
    if(err || calls<7 || lis_test_oom_live()!=baseline)
        goto done;

    for(fail_at=1;fail_at<=calls;fail_at++)
    {
        lis_test_oom_begin(fail_at);
        err=lis_solver_set_option((char *)options,solver);
        lis_test_oom_end();
        hit=lis_test_oom_calls();
        if(err!=LIS_OUT_OF_MEMORY || hit!=fail_at ||
           lis_test_oom_live()!=baseline)
        {
            fprintf(stderr,
                    "ILUT_OPTIONS_OOM_FAIL fail_at=%lu err=%ld count=%lu\n",
                    (unsigned long)fail_at,(long)err,(unsigned long)hit);
            goto done;
        }
    }
    printf("ILUT_OPTIONS_OOM_PASS checked=%lu\n",(unsigned long)calls);
    failed=0;
done:
    lis_test_oom_end();
    if(solver) lis_solver_destroy(solver);
    return failed;
#endif
}

/*
 * Test-only fault injection: allow one local ILUTP worker to see OOM
 * while the other OpenMP workers are still active. The creator must join
 * all workers and free every partially constructed domain.
 *
 * LIS_TEST_OOM is not compiled into production libraries.
 */
static int check_ilut_ras_parallel_oom(void)
{
#if !defined(LIS_TEST_OOM) || !defined(_OPENMP)
    return 0;
#else
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT err=0,i;
    size_t baseline,calls,targets[3],t,hit;
    int failed=1,old_threads=omp_get_max_threads();
    int old_dynamic=omp_get_dynamic();

    omp_set_dynamic(0);
    omp_set_num_threads(8);
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto cleanup;
    err=lis_matrix_set_size(A,0,32);
    if(err) goto cleanup;
    if(A->nprocs!=1) { failed=0; goto cleanup; }
    for(i=0;i<32;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,
                                 (LIS_SCALAR)4.0,A);
        if(err) goto cleanup;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                                     (LIS_SCALAR)-0.15,A);
            if(err) goto cleanup;
        }
        if(i+1<32)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)-0.1,A);
            if(err) goto cleanup;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto cleanup;
    err=lis_matrix_assemble(A);
    if(err) goto cleanup;
    err=lis_solver_create(&solver);
    if(err) goto cleanup;
    err=lis_solver_set_option(
        "-i fgmres -p ilut -scale none -print none "
        "-ilutp 1 -ilutp_full 1 -ilutp_full_limit 512 "
        "-ilut_ras_domains 8 -ilut_ras_overlap 1 -ilut_schur 0",
        solver);
    if(err) goto cleanup;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto cleanup;

    baseline=lis_test_oom_live();
    lis_test_oom_begin(0);
    err=lis_precon_create(solver,&pc);
    lis_test_oom_end();
    calls=lis_test_oom_calls();
    if(err) { pc=NULL; goto cleanup; }
    lis_precon_destroy(pc);pc=NULL;
    if(calls<20 || lis_test_oom_live()!=baseline)
    {
        fprintf(stderr,"PARALLEL_RAS_OOM baseline invalid calls=%lu\n",
                (unsigned long)calls);
        goto cleanup;
    }
    targets[0]=calls/3;
    targets[1]=(2*calls)/3;
    targets[2]=calls-1;
    for(t=0;t<3;t++)
    {
        lis_test_oom_begin(targets[t]);
        err=lis_precon_create(solver,&pc);
        lis_test_oom_end();
        hit=lis_test_oom_calls();
        /* The constructor frees the partially initialized PC on error. */
        if(err==LIS_SUCCESS)
        {
            if(pc) { lis_precon_destroy(pc);pc=NULL; }
            fprintf(stderr,"PARALLEL_RAS_OOM did not fail at %lu\n",
                    (unsigned long)targets[t]);
            goto cleanup;
        }
        pc=NULL;
        if(err!=LIS_OUT_OF_MEMORY || hit<targets[t] ||
           lis_test_oom_live()!=baseline)
        {
            fprintf(stderr,
                    "PARALLEL_RAS_OOM leak/status at %lu err=%ld hit=%lu\n",
                    (unsigned long)targets[t],(long)err,
                    (unsigned long)hit);
            goto cleanup;
        }
        printf("PARALLEL_RAS_OOM injected=%lu live=0 PASS\n",
               (unsigned long)targets[t]);
    }
    failed=0;
cleanup:
    lis_test_oom_end();
    if(pc) lis_precon_destroy(pc);
    if(solver) lis_solver_destroy(solver);
    if(A) lis_matrix_destroy(A);
    omp_set_num_threads(old_threads);
    omp_set_dynamic(old_dynamic);
    return failed;
#endif
}

/*
 * Deterministic stabilized contact-KKT-like system:
 * [ K  B^T ][u] = [f]
 * [ B -0.1I][λ]   [g]
 *
 * K is an SPD spring chain, B enforces displacement differences.
 * The negative multiplier block makes the matrix indefinite.
 * This is an algebraic regression, NOT a full nonlinear contact FEA run.
 */
static int
check_ilut_ras_contact_kkt(void)
{
    const LIS_INT nu=36,nlambda=12,n=48;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,u=NULL,r=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,i,q,iters=0,variant;
    LIS_REAL bn=0.0,rn=0.0;
    char opts[320];
    int failed=1;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,n);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<nu;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)4.0,A);
        if(err) goto done;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                                     (LIS_SCALAR)-0.8,A);
            if(err) goto done;
        }
        if(i+1<nu)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)-0.8,A);
            if(err) goto done;
        }
    }
    for(q=0;q<nlambda;q++)
    {
        LIS_INT l=nu+q,a=3*q,c=3*q+1;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,l,
                                 (LIS_SCALAR)-0.1,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,a,(LIS_SCALAR)1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,c,(LIS_SCALAR)-1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,a,l,(LIS_SCALAR)1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,c,l,(LIS_SCALAR)-1.0,A);
        if(err) goto done;
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    if(lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&u) ||
       lis_vector_duplicate(A,&r)) goto done;
    lis_vector_set_all((LIS_SCALAR)1.0,u);
    err=lis_matvec(A,u,b);
    if(err) goto done;
    err=lis_vector_nrm2(b,&bn);
    if(err || bn<=0.0) goto done;
    for(variant=0;variant<2;variant++)
    {
        solver=NULL;
        err=lis_solver_create(&solver);
        if(err) goto done;
        snprintf(opts,sizeof(opts),
            "-i fgmres -p ilut -scale none -print none "
            "-tol 1e-8 -maxiter 120 -restart 30 "
            "-iluc_drop 0.005 -iluc_rate 8 -ilu_pivot_tol 1e-9 %s",
            variant ? "-ilut_ras_domains 4 -ilut_ras_overlap 1 -ilut_schur 1" : "");
        err=lis_solver_set_option(opts,solver);
        if(err) goto done;
        lis_vector_set_all((LIS_SCALAR)0.0,x);
        err=lis_solve(A,b,x,solver);
        if(err) goto done;
        err=lis_matvec(A,x,r);
        if(err) goto done;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) goto done;
        err=lis_vector_nrm2(r,&rn);
        if(err) goto done;
        err=lis_solver_get_iter(solver,&iters);
        if(err || rn/bn>1e-7 || iters<1)
        {
            fprintf(stderr,"ILUT_RAS_CONTACT_KKT variant=%d iter=%d residual=%.3e\n",
                    variant,(int)iters,(double)(rn/bn));
            goto done;
        }
        printf("ILUT_RAS_CONTACT_KKT variant=%d iter=%d residual=%.6e\n",
               variant,(int)iters,(double)(rn/bn));
        lis_solver_destroy(solver); solver=NULL;
    }
    failed=0;
done:
    if(failed)
        fprintf(stderr,"ILUT_RAS_CONTACT_KKT_FAILED rc=%d\n",(int)err);
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(u) lis_vector_destroy(u);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    return failed;
}


/*
 * Negative regression: the Schur approximation MUST reject a weak
 * interior diagonal instead of quietly treating diag(A_II)^-1 as zero.
 * Row 1 is interior for the 2-domain partition of this 8x8 system.
 */
static int
check_ilut_ras_unsafe_interior(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT i,err=0;
    int failed=1;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,8);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<8;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,
             i==1 ? (LIS_SCALAR)1.0e-20 : (LIS_SCALAR)3.0,A);
        if(err) goto done;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                                     (LIS_SCALAR)-0.4,A);
            if(err) goto done;
        }
        if(i<7)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)-0.4,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_option(
        "-i fgmres -p ilut -scale none -print none "
        "-iluc_drop 0 -iluc_rate 8 -ilu_pivot_tol 1e-9 "
        "-ilut_ras_domains 2 -ilut_ras_overlap 1 -ilut_schur 1",
        solver);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;
    err=lis_precon_create(solver,&pc);
    if(err==LIS_ERR_ILL_ARG)
    {
        pc=NULL; /* failing constructor releases its partial preconditioner */
        printf("ILUT_RAS_WEAK_DIAG_REJECTED\n");
        failed=0;
        goto done;
    }
    fprintf(stderr,"Schur accepted weak interior pivot; rc=%ld\n",(long)err);
done:
    if(pc) lis_precon_destroy(pc);
    if(solver) lis_solver_destroy(solver);
    if(A) lis_matrix_destroy(A);
    return failed;
}


/*
 * High-unit-contrast 48x48 contact-KKT algebraic regression:
 *       [ 1e8 K   B^T ] [u] = [f]
 *       [ B     -1e-8 I] [l]   [g]
 *
 * Compare unscaled, row equilibration, RCM+row equilibration, and
 * a manually block-balanced congruence D A D, where
 * D_u=1e-4 and D_l=1e4. The latter is a diagnostic reference,
 * NOT an implemented automatic two-sided equilibrator.
 *
 * Always check both the original true residual and D(Ax-b); the
 * latter gives multiplier equations meaningful weight despite
 * displacement stiffness being orders of magnitude larger.
 * Only the block-balanced reference is a hard convergence gate:
 * the other variants intentionally record numerical failures.
 * This is not a production contact FEA performance benchmark.
 */
static int
check_ilut_kkt_unit_contrast(void)
{
    const LIS_INT nu=36,nlambda=12,n=48;
    const double stiffness=1.0e8;
    const double du=1.0e-4,dl=1.0e4;
    const char *labels[]={"baseline","row","rcm_row","block_balanced"};
    const char *flags[]={"","-ilut_row_scale 1",
                          "-ilut_rcm 1 -ilut_row_scale 1",
                          "-ilut_rcm 1"};
    LIS_MATRIX A=NULL,Ab=NULL,M;
    LIS_VECTOR truth=NULL,b=NULL,bb=NULL,x=NULL,xb=NULL,r=NULL,wr=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,pass,i,q,variant,iters=0;
    LIS_REAL bn=0.0,bbn=0.0,rn=0.0,wrn=0.0;
    char opts[320];
    int failed=1,balanced_ok=0;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,n);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    err=lis_matrix_create(LIS_COMM_WORLD,&Ab);
    if(err) goto done;
    err=lis_matrix_set_size(Ab,0,n);
    if(err) goto done;

    for(pass=0;pass<2;pass++)
    {
        M=pass ? Ab : A;
        for(i=0;i<nu;i++)
        {
            double factor=pass ? du*du : 1.0;
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i,
                                     (LIS_SCALAR)(4.0*stiffness*factor),M);
            if(err) goto done;
            if(i>0)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                         (LIS_SCALAR)(-0.8*stiffness*factor),M);
                if(err) goto done;
            }
            if(i+1<nu)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                         (LIS_SCALAR)(-0.8*stiffness*factor),M);
                if(err) goto done;
            }
        }
        for(q=0;q<nlambda;q++)
        {
            LIS_INT l=nu+q,a=3*q,c=3*q+1;
            double cross=pass ? du*dl : 1.0;
            double cc=pass ? dl*dl : 1.0;
            err=lis_matrix_set_value(LIS_INS_VALUE,l,l,
                                     (LIS_SCALAR)(-1e-8*cc),M);
            if(err) goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,l,a,
                                     (LIS_SCALAR)cross,M);
            if(err) goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,l,c,
                                     (LIS_SCALAR)(-cross),M);
            if(err) goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,a,l,
                                     (LIS_SCALAR)cross,M);
            if(err) goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,c,l,
                                     (LIS_SCALAR)(-cross),M);
            if(err) goto done;
        }
        err=lis_matrix_set_type(M,LIS_MATRIX_CSR);
        if(err) goto done;
        err=lis_matrix_assemble(M);
        if(err) goto done;
    }

    if(lis_vector_duplicate(A,&truth) ||
       lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&bb) ||
       lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&xb) ||
       lis_vector_duplicate(A,&r) ||
       lis_vector_duplicate(A,&wr)) goto done;

    for(i=0;i<n;i++)
    {
        double di=i<nu ? du : dl;
        double yi=(i%2 ? -1.0 : 1.0)*(1.0+0.01*i);
        truth->value[i]=(LIS_SCALAR)(di*yi);
    }
    err=lis_matvec(A,truth,b);
    if(err) goto done;
    for(i=0;i<n;i++)
        bb->value[i]=b->value[i]*(i<nu ? du : dl);
    if(lis_vector_nrm2(b,&bn) || lis_vector_nrm2(bb,&bbn) ||
       !(bn>0.0) || !(bbn>0.0)) goto done;

    for(variant=0;variant<4;variant++)
    {
        LIS_MATRIX target=variant==3 ? Ab : A;
        LIS_VECTOR rhs=variant==3 ? bb : b;
        LIS_VECTOR solution=variant==3 ? xb : x;
        LIS_INT nnzl=0,nnzu=0;
        LIS_REAL minpivot=0.0,maxpivot=0.0;
        int solve_rc;
        solver=NULL;
        err=lis_solver_create(&solver);
        if(err) goto done;
        snprintf(opts,sizeof(opts),
                 "-i fgmres -p ilut -scale none -print none "
                 "-tol 1e-9 -maxiter 250 -restart 40 "
                 "-ilut_drop 0.001 -ilut_rate 12 "
                 "-ilu_pivot_tol 1e-12 %s",flags[variant]);
        err=lis_solver_set_option(opts,solver);
        if(err) goto done;
        lis_vector_set_all((LIS_SCALAR)0.0,solution);
        solve_rc=lis_solve(target,rhs,solution,solver);
        if(solve_rc)
        {
            printf("ILUT_KKT_CONTRAST variant=%s solve_rc=%d\n",
                   labels[variant],(int)solve_rc);
            lis_solver_destroy(solver); solver=NULL;
            continue;
        }
        if(variant==3)
            for(i=0;i<n;i++)
                x->value[i]=xb->value[i]*(i<nu ? du : dl);
        err=lis_matvec(A,x,r);
        if(err) goto done;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) goto done;
        for(i=0;i<n;i++)
            wr->value[i]=r->value[i]*(i<nu ? du : dl);
        if(lis_vector_nrm2(r,&rn) || lis_vector_nrm2(wr,&wrn))
            goto done;
        err=lis_solver_get_iter(solver,&iters);
        if(err) goto done;
        if(solver->precon)
            (void)lis_precon_ilut_factor_stats(solver->precon,
                    &nnzl,&nnzu,&minpivot,&maxpivot);
        printf("ILUT_KKT_CONTRAST variant=%s iter=%d "
               "true_residual=%.6e balanced_residual=%.6e "
               "nnz_lu=%d min_pivot=%.6e max_pivot=%.6e\n",
               labels[variant],(int)iters,(double)(rn/bn),
               (double)(wrn/bbn),(int)(nnzl+nnzu),
               (double)minpivot,(double)maxpivot);
        if(variant==3 && isfinite((double)rn) &&
           isfinite((double)wrn) && rn/bn<1e-7 && wrn/bbn<1e-7 &&
           iters>=1)
            balanced_ok=1;
        lis_solver_destroy(solver); solver=NULL;
    }
    failed=!balanced_ok;
done:
    if(failed)
        fprintf(stderr,"ILUT_KKT_CONTRAST_FAILED rc=%d balanced_ok=%d\n",
                (int)err,balanced_ok);
    if(solver) lis_solver_destroy(solver);
    if(wr) lis_vector_destroy(wr);
    if(r) lis_vector_destroy(r);
    if(xb) lis_vector_destroy(xb);
    if(x) lis_vector_destroy(x);
    if(bb) lis_vector_destroy(bb);
    if(b) lis_vector_destroy(b);
    if(truth) lis_vector_destroy(truth);
    if(Ab) lis_matrix_destroy(Ab);
    if(A) lis_matrix_destroy(A);
    return failed;
}


/*
 * Manually-invoked external Matrix Market FEA review harness.
 *
 * ./test/testhalfilut --fea-matrix A.mtx BASE_RATE RCM_RATE [DROP]
 *
 * Both modes solve the *original* system A x = A x_true, not the reordered
 * or equilibrated surrogate. Adjust the rates until factor nnz match to
 * within 5%; the output states explicitly if they do not.
 *
 * Intentionally NOT part of make check: no private FEA data are shipped.
 */
typedef struct {
    LIS_INT iterations,nnz_l,nnz_u,error,threads;
    LIS_REAL residual,min_pivot,max_pivot;
    double setup_cpu_s,solve_cpu_s;
    double setup_wall_s,solve_wall_s;
    int target_met;
} LIS_ILUT_FEA_METRICS;

/* Optional explicit thread budgets for the two ILUT variants, without
 * changing any global solver defaults or the standalone CLI contract.
 * OpenMP thread counts used to build and apply one ILUT factor must match.
 */
static LIS_INT
fea_variant_threads(LIS_INT reordered,int *threads)
{
    const char *env=getenv(reordered ?
                         "LIS_FEA_RCM_THREADS":"LIS_FEA_BASELINE_THREADS");
    char *end=NULL;
    long requested;
#ifdef _OPENMP
    *threads=omp_get_max_threads();
#else
    *threads=1;
#endif
    if(!env || !*env) return LIS_SUCCESS;
    requested=strtol(env,&end,10);
    if(end==env || *end!='\0' || requested<1 || requested>64)
        return LIS_ERR_ILL_ARG;
#ifndef _OPENMP
    if(requested!=1) return LIS_ERR_NOT_IMPLEMENTED;
#endif
    *threads=(int)requested;
    return LIS_SUCCESS;
}

static int
run_ilut_fea_variant(LIS_MATRIX A,LIS_VECTOR b,LIS_INT reordered,
                     double rate,double drop,LIS_ILUT_FEA_METRICS *m)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR x=NULL,r=NULL;
    LIS_INT err=LIS_SUCCESS;
    LIS_REAL bn=0.0,rn=0.0;
    clock_t c0,c1,c2;
    double w0,w1,w2;
    char opts[360];
    int failed=1,thread_count=1;
#ifdef _OPENMP
    int original_threads=omp_get_max_threads();
#endif
    if(!m) return 1;
    memset(m,0,sizeof(*m));
    m->iterations=-1;
    m->residual=-1.0;
    if(!isfinite(rate) || !isfinite(drop) ||
       rate<=0.0 || drop<0.0 || drop>=1.0)
    { err=LIS_ERR_ILL_ARG; goto done; }
    err=fea_variant_threads(reordered,&thread_count);
    if(err) goto done;
    m->threads=thread_count;
#ifdef _OPENMP
    omp_set_num_threads(thread_count);
#endif
    snprintf(opts,sizeof(opts),
             "-i fgmres -p ilut -scale none -print none -tol 1e-8 "
             "-maxiter 500 -restart 40 -ilu_pivot_tol 1e-10 "
             "-ilut_rate %.12g -ilut_drop %.12g %s",
             rate,drop,
             reordered ? "-ilut_rcm 1 -ilut_row_scale 1" : "");
    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_option(opts,solver);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;
    err=lis_vector_duplicate(A,&x);
    if(err) goto done;
    err=lis_vector_duplicate(A,&r);
    if(err) goto done;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);
    if(err) goto done;

    c0=clock(); w0=lis_wtime();
    err=lis_precon_create(solver,&pc);
    c1=clock(); w1=lis_wtime();
    m->setup_cpu_s=(double)(c1-c0)/CLOCKS_PER_SEC;
    m->setup_wall_s=w1-w0;
    if(err) { pc=NULL; goto done; }
    solver->precon=pc;
    err=lis_precon_ilut_factor_stats(pc,&m->nnz_l,&m->nnz_u,
                                     &m->min_pivot,&m->max_pivot);
    if(err) goto done;

    err=lis_solve_kernel(A,b,x,solver,pc);
    c2=clock(); w2=lis_wtime();
    m->solve_cpu_s=(double)(c2-c1)/CLOCKS_PER_SEC;
    m->solve_wall_s=w2-w1;
    if(err) goto done;
    err=lis_solver_get_iter(solver,&m->iterations);
    if(err) goto done;
    err=lis_matvec(A,x,r);
    if(err) goto done;
    err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
    if(err) goto done;
    err=lis_vector_nrm2(r,&rn);
    if(err) goto done;
    err=lis_vector_nrm2(b,&bn);
    if(err || !(bn>0.0)) { err=LIS_ERR_ILL_ARG; goto done; }
    m->residual=rn/bn;
    if(!isfinite((double)m->residual))
    { err=LIS_ERR_ILL_ARG; goto done; }

    m->target_met=m->residual<=1e-8;
    /* Preserve the previous smoke-test acceptance band (1e-7),
     * but never confuse it with the requested 1e-8 solve target. */
    if(m->residual>1e-7) goto done;
    failed=0;
done:
    m->error=err;
    if(failed)
        fprintf(stderr,
                "ILUT_FEA_BENCH_FAILED mode=%s rc=%ld true_relres=%.9g\n",
                reordered?"rcm":"baseline",(long)err,(double)m->residual);
    if(pc)
    {
        if(solver) solver->precon=NULL;
        lis_precon_destroy(pc);
    }
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(x) lis_vector_destroy(x);
#ifdef _OPENMP
    omp_set_num_threads(original_threads);
#endif
    return failed;
}

static int
run_ilut_fea_matrix(const char *file,double baseline_rate,
                    double rcm_rate,double drop)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR exact=NULL,b=NULL;
    LIS_ILUT_FEA_METRICS m[2];
    LIS_INT i,err=0;
    double fill_gap=-1.0;
    int failed=1,failed_baseline,failed_rcm;
    if(!file || !isfinite(baseline_rate) || !isfinite(rcm_rate) ||
       baseline_rate<=0.0 || rcm_rate<=0.0 ||
       !isfinite(drop) || drop<0.0 || drop>=1.0)
        return 1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_input_matrix(A,(char *)file);
    if(err) goto done;
    if(A->nprocs!=1 || A->matrix_type!=LIS_MATRIX_CSR ||
       A->n!=A->gn || A->n<2)
    {
        err=LIS_ERR_NOT_IMPLEMENTED;
        goto done;
    }
    err=lis_vector_duplicate(A,&exact);
    if(err) goto done;
    err=lis_vector_duplicate(A,&b);
    if(err) goto done;
    /* Avoid an all-ones null mode in elasticity stiffness matrices. */
    for(i=0;i<A->n;i++)
        exact->value[i]=(LIS_SCALAR)(1.0+0.017*(double)((i*19)%29));
    err=lis_matvec(A,exact,b);
    if(err) goto done;

    /* Execute both variants even if the first fails or does not converge.
     * A negative numerical benchmark is still useful evidence. */
    failed_baseline=run_ilut_fea_variant(A,b,0,baseline_rate,drop,&m[0]);
    failed_rcm=run_ilut_fea_variant(A,b,1,rcm_rate,drop,&m[1]);
    if(m[0].nnz_l+m[0].nnz_u>0 && m[1].nnz_l+m[1].nnz_u>0)
        fill_gap=fabs((double)(m[0].nnz_l+m[0].nnz_u) -
                      (double)(m[1].nnz_l+m[1].nnz_u))/
                      (double)(m[0].nnz_l+m[0].nnz_u);
    printf("ILUT_FEA_MATRIX file=%s dof=%ld nnz_A=%ld drop=%.9g "
           "fill_gap=%.5f fair_fill=%s\n",
           file,(long)A->n,(long)A->nnz,drop,fill_gap,
           fill_gap<0.0 ? "unknown" : fill_gap<=0.05 ? "yes" : "no");
    /* Preserve legacy columns and append wall, threading, and status. */
    printf("variant,rate,nnz_L,nnz_U,nnz_LU,iterations,"
           "true_relres,min_pivot,max_pivot,setup_cpu_s,solve_cpu_s,"
           "setup_wall_s,solve_wall_s,total_wall_s,threads,"
           "target_1e-8_met,solver_error,status\n");
    for(i=0;i<2;i++)
    {
        const LIS_ILUT_FEA_METRICS *v=&m[i];
        const char *status=v->error ? "solver_error" :
                           v->target_met ? "converged" : "not_converged";
        printf("%s,%.9g,%ld,%ld,%ld,%ld,%.9g,%.9g,%.9g,"
               "%.6f,%.6f,%.6f,%.6f,%.6f,%ld,%s,%ld,%s\n",
               i==0?"baseline":"rcm_scaled",
               i==0?baseline_rate:rcm_rate,
               (long)v->nnz_l,(long)v->nnz_u,
               (long)(v->nnz_l+v->nnz_u),
               (long)v->iterations,(double)v->residual,
               (double)v->min_pivot,(double)v->max_pivot,
               v->setup_cpu_s,v->solve_cpu_s,
               v->setup_wall_s,v->solve_wall_s,
               v->setup_wall_s+v->solve_wall_s,
               (long)v->threads,v->target_met?"yes":"no",
               (long)v->error,status);
    }
    failed=failed_baseline || failed_rcm;
done:
    if(failed)
        fprintf(stderr,"ILUT_FEA_MATRIX_FAILED file=%s rc=%ld\n",
                file?file:"(null)",(long)err);
    if(b) lis_vector_destroy(b);
    if(exact) lis_vector_destroy(exact);
    if(A) lis_matrix_destroy(A);
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

	if(argc>=5 && argc<=6 &&
	   strcmp(argv[1],"--fea-matrix")==0)
	{
		char *end1=NULL,*end2=NULL,*end3=NULL;
		double baseline_rate=strtod(argv[3],&end1);
		double rcm_rate=strtod(argv[4],&end2);
		double drop=argc==6 ? strtod(argv[5],&end3) : 0.005;
		if(end1==argv[3] || *end1!='\0' ||
		   end2==argv[4] || *end2!='\0' ||
		   (argc==6 && (end3==argv[5] || *end3!='\0')))
			failed=1;
		else
			failed=run_ilut_fea_matrix(argv[2],baseline_rate,rcm_rate,drop);
		lis_finalize();
		return failed ? 1 : 0;
	}

    /* Run the matrix conversion OOM regression independently of ILUTP. */
    if(argc==2 && strcmp(argv[1],"--rco-assemble-oom")==0)
    {
#ifndef LIS_TEST_OOM
        fprintf(stderr,"RCO assembly OOM test requires --enable-oom-tests\n");
        failed=1;
#else
        failed=check_rco_assemble_oom_cleanup();
#endif
        lis_finalize();
        return failed ? 1 : 0;
    }

	/* Isolate parallel ILUTP from unrelated legacy half-ILUT tests:
	 * the latter were not originally designed for OMP_NUM_THREADS=8. */
	if(argc==2 && strcmp(argv[1],"--parallel-ras-ilutp")==0)
	{
		failed=check_ilut_ras_local_ilutp();
		failed |= check_ilut_ras_parallel_overlap();
		failed |= check_ilut_ras_parallel_oom();
		failed |= check_ilutp_options_oom();
		failed |= check_rco_assemble_oom_cleanup();
		lis_finalize();
		return failed ? 1 : 0;
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
	failed |= check_ilut_contact_kkt();
	failed |= check_ilut_kkt_unit_contrast();
	failed |= check_ilut_ras_schur(A);
	failed |= check_ilut_ras_exact_interface();
	failed |= check_ilut_ras_local_ilutp();
	failed |= check_ilut_ras_parallel_overlap();
	failed |= check_ilut_ras_unsafe_interior();
	failed |= check_ilut_ras_contact_kkt();

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
