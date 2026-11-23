#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include "lis.h"
#include "lis_solver.h"
#include "lis_precon.h"
#include <stdio.h>

typedef struct
{
	LIS_INT n;
	LIS_SCALAR diag;
	LIS_SCALAR offdiag;
	LIS_INT fail_next_matvec;
	LIS_INT fail_next_psolve;
} USER_TRIDIAG;

static LIS_INT shell_matvec(void *ctx, const LIS_SCALAR *x, LIS_SCALAR *y)
{
	USER_TRIDIAG *A = (USER_TRIDIAG *)ctx;
	LIS_INT i;

	if( A->fail_next_matvec )
	{
		A->fail_next_matvec = LIS_FALSE;
		return LIS_ERR_ILL_ARG;
	}

	if( A->n==1 )
	{
		y[0] = A->diag*x[0];
		return LIS_SUCCESS;
	}

	y[0] = A->diag*x[0] + A->offdiag*x[1];
	for(i=1;i<A->n-1;i++)
	{
		y[i] = A->offdiag*x[i-1] + A->diag*x[i] + A->offdiag*x[i+1];
	}
	y[A->n-1] = A->offdiag*x[A->n-2] + A->diag*x[A->n-1];
	return LIS_SUCCESS;
}

/* Real symmetric test operator: A^H == A. */
static LIS_INT shell_matvech(void *ctx, const LIS_SCALAR *x, LIS_SCALAR *y)
{
	return shell_matvec(ctx,x,y);
}

static LIS_INT userdiag_create(LIS_SOLVER solver, LIS_PRECON precon)
{
	(void)solver;
	(void)precon;
	return LIS_SUCCESS;
}

static LIS_INT userdiag_psolve(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
	USER_TRIDIAG *A;
	void *ctx = NULL;
	LIS_INT i, err;

	err = lis_matrix_get_user_data(solver->A,&ctx);
	if( err ) return err;
	A = (USER_TRIDIAG *)ctx;
	if( A->fail_next_psolve )
	{
		A->fail_next_psolve = LIS_FALSE;
		return LIS_ERR_ILL_ARG;
	}

	for(i=0;i<b->n;i++) x->value[i] = b->value[i] / A->diag;
	return LIS_SUCCESS;
}

static LIS_INT userdiag_psolveh(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
	return userdiag_psolve(solver,b,x);
}

/* Verify an error on the first nonzero-initial-guess SpMV cannot be
 * misreported as convergence, even when certification's later SpMV succeeds.
 * An independent successful solve also checks stale status is not retained. */
static int test_initial_matvec_error(LIS_MATRIX A, LIS_VECTOR b,
                                     LIS_VECTOR x, USER_TRIDIAG *ctx)
{
    const char *methods[] = {
        "cg", "bicg", "cgs", "bicgstab", "bicgsafe",
        "bicr", "bicrsafe", "gmres", "fgmres"
    };
    LIS_INT k,err,status;
    LIS_SOLVER s=NULL;
    char options[160];

    /* Rank-local failure injection is intentionally not an MPI test. */
    if(A->nprocs!=1) return 0;

    for(k=0;k<(LIS_INT)(sizeof(methods)/sizeof(methods[0]));k++)
    {
        err=lis_solver_create(&s);
        if(err) return 1;
        snprintf(options,sizeof(options),
                 "-i %s -p none -tol 1e-11 -maxiter 1000 -print none",
                 methods[k]);
        err=lis_solver_set_option(options,s);
        if(err) goto fail;

        s->options[LIS_OPTIONS_INITGUESS_ZEROS]=LIS_FALSE;
        err=lis_vector_set_all((LIS_SCALAR)0.5,x);
        if(err) goto fail;
        ctx->fail_next_matvec=LIS_TRUE;

        err=lis_solve(A,b,x,s);
        if(err) goto fail;
        err=lis_solver_get_status(s,&status);
        if(err || status!=LIS_ERR_ILL_ARG ||
           ctx->fail_next_matvec || s->iter!=0)
        {
            fprintf(stderr,"initial USER callback error lost: %s status=%lld iter=%lld\n",
                    methods[k],(long long)status,(long long)s->iter);
            goto fail;
        }

        /* Warm-start from the exact solution: a successful early exit
         * must clear the previous callback failure deterministically. */
        ctx->fail_next_matvec=LIS_FALSE;
        s->options[LIS_OPTIONS_INITGUESS_ZEROS]=LIS_FALSE;
        err=lis_vector_set_all((LIS_SCALAR)1.0,x);
        if(err) goto fail;
        err=lis_solve(A,b,x,s);
        if(err) goto fail;
        err=lis_solver_get_status(s,&status);
        if(err || status!=LIS_SUCCESS)
        {
            fprintf(stderr,"USER callback recovery failed: %s status=%lld\n",
                    methods[k],(long long)status);
            goto fail;
        }

        lis_solver_destroy(s);
        s=NULL;
    }
    printf("USER_INITIAL_RESIDUAL_ERROR PASS\n");
    return 0;
fail:
    ctx->fail_next_matvec=LIS_FALSE;
    if(s) lis_solver_destroy(s);
    return 1;
}

/* Directly exercise the shared initial-residual M!=NULL branch.
 * The ordinary default-precision Krylov paths pass M==NULL, so
 * injecting a psolve error in their iterations would not test this fix.
 * QMR's quad variant uses M!=NULL; this is a default-precision helper
 * test so it is also runnable in configurations with quad disabled. */
static int test_initial_psolve_error(LIS_MATRIX A, LIS_VECTOR b,
                                     LIS_VECTOR x, USER_TRIDIAG *ctx)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_VECTOR scratch=NULL,residual=NULL,rhs_copy=NULL;
    LIS_REAL bnrm2=0,rhs_diff=0;
    LIS_INT err,result;
    int failed=1;
    char options[]="-i cg -p userdiag -tol 1e-12 -print none";

    if(A->nprocs!=1) return 0;
    err=lis_solver_create(&solver);if(err)goto done;
    err=lis_solver_set_option(options,solver);if(err)goto done;
    solver->A=A;
    solver->b=b;
    solver->x=x;
    solver->precision=LIS_PRECISION_DEFAULT;
    solver->options[LIS_OPTIONS_INITGUESS_ZEROS]=LIS_TRUE;
    err=lis_precon_create(solver,&precon);if(err)goto done;
    solver->precon=precon;
    err=lis_vector_duplicate(A,&scratch);if(err)goto done;
    err=lis_vector_duplicate(A,&residual);if(err)goto done;
    err=lis_vector_duplicate(A,&rhs_copy);if(err)goto done;
    err=lis_vector_copy(b,rhs_copy);if(err)goto done;

    solver->retcode=LIS_SUCCESS;
    solver->iter=101;
    ctx->fail_next_psolve=LIS_TRUE;
    result=lis_solver_get_initial_residual(
        solver,precon,scratch,residual,&bnrm2);
    if(result!=LIS_FAILS || solver->retcode!=LIS_ERR_ILL_ARG ||
       solver->iter!=0 || ctx->fail_next_psolve)
    {
        fprintf(stderr,
                "USER initial psolve error was lost: result=%lld status=%lld iter=%lld\n",
                (long long)result,(long long)solver->retcode,
                (long long)solver->iter);
        goto done;
    }

    /* A second call must succeed with the same preconditioner, and
     * the input RHS must not be accidentally changed by the failure. */
    solver->retcode=LIS_SUCCESS;
    result=lis_solver_get_initial_residual(
        solver,precon,scratch,residual,&bnrm2);
    if(result!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS ||
       ctx->fail_next_psolve)
    {
        fprintf(stderr,"USER initial psolve recovery failed: %lld\n",
                (long long)result);
        goto done;
    }
    err=lis_vector_axpy((LIS_SCALAR)-1,b,rhs_copy);if(err)goto done;
    err=lis_vector_nrm2(rhs_copy,&rhs_diff);if(err)goto done;
    if(rhs_diff!=(LIS_REAL)0.0)
    {
        fprintf(stderr,"USER initial psolve changed RHS: %.17g\n",
                (double)rhs_diff);
        goto done;
    }
    printf("USER_INITIAL_PSOLVE_ERROR PASS\n");
    failed=0;

done:
    ctx->fail_next_psolve=LIS_FALSE;
    if(solver) solver->precon=NULL; /* caller owns this test preconditioner */
    if(precon) lis_precon_destroy(precon);
    if(rhs_copy) lis_vector_destroy(rhs_copy);
    if(residual) lis_vector_destroy(residual);
    if(scratch) lis_vector_destroy(scratch);
    if(solver) lis_solver_destroy(solver);
    return failed;
}

int main(int argc, char *argv[])
{
	LIS_MATRIX A;
	LIS_VECTOR b, x, exact;
	LIS_SOLVER solver;
	USER_TRIDIAG shell;
	LIS_INT i, err;
	LIS_INT local_n, global_n;
	LIS_REAL nrm2;
	const LIS_INT n = 100;

	lis_initialize(&argc,&argv);

	shell.n = n;
	shell.diag = 2.0;
	shell.offdiag = -1.0;
	shell.fail_next_matvec = LIS_FALSE;
	shell.fail_next_psolve = LIS_FALSE;

	err = lis_precon_register("userdiag",
	                          userdiag_create,
	                          userdiag_psolve,
	                          userdiag_psolveh);
	if( err )
	{
		fprintf(stderr,"lis_precon_register failed: %d\n",(int)err);
		lis_finalize();
		return 1;
	}

	lis_matrix_create(LIS_COMM_WORLD,&A);
	lis_matrix_set_size(A,0,n);
	lis_matrix_get_size(A,&local_n,&global_n);
	(void)global_n;
	shell.n = local_n;
	err = lis_matrix_set_user(A,&shell,shell_matvec,shell_matvech);
	if( err )
	{
		fprintf(stderr,"lis_matrix_set_user failed: %d\n",(int)err);
		lis_matrix_destroy(A);
		lis_finalize();
		return 1;
	}

	lis_vector_duplicate(A,&b);
	lis_vector_duplicate(A,&x);
	lis_vector_duplicate(A,&exact);
	for(i=0;i<exact->n;i++) exact->value[i] = 1.0;

	err = lis_matvec(A,exact,b);
	if( err )
	{
		fprintf(stderr,"lis_matvec failed: %d\n",(int)err);
		return 1;
	}
	lis_vector_set_all(0.0,x);

	lis_solver_create(&solver);
	err = lis_solver_set_option(
		"-i bicg -p userdiag -tol 1.0e-12 -maxiter 1000 -print none",
		solver);
	if( err )
	{
		fprintf(stderr,"lis_solver_set_option failed: %d\n",(int)err);
		return 1;
	}

	err = lis_solve(A,b,x,solver);
	if( err )
	{
		fprintf(stderr,"lis_solve failed: %d\n",(int)err);
		return 1;
	}

	lis_vector_axpy(-1.0,exact,x);
	lis_vector_nrm2(x,&nrm2);
	printf("LIS_MATRIX_USER solution error = %.15e\n",(double)nrm2);

	if( nrm2<=1.0e-8 && (test_initial_matvec_error(A,b,x,&shell) ||
	                        test_initial_psolve_error(A,b,x,&shell)) )
	{
		fprintf(stderr,"USER initial residual regression FAILED\n");
		return 1;
	}

	lis_solver_destroy(solver);
	lis_vector_destroy(exact);
	lis_vector_destroy(x);
	lis_vector_destroy(b);
	lis_matrix_destroy(A);
	lis_precon_register_free();
	lis_finalize();

	if( nrm2>1.0e-8 )
	{
		fprintf(stderr,"FAILED\n");
		return 1;
	}
	printf("PASSED\n");
	return 0;
}
