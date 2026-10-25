#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include <math.h>
#include <stdio.h>
#include "lis.h"
#include "lis_precon.h"

static LIS_INT
create_nonsymmetric_matrix(LIS_MATRIX *A, LIS_INT n)
{
    LIS_INT err;
    LIS_INT i,is,ie;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,n);
    if( err ) return err;

    err = lis_matrix_get_range(*A,&is,&ie);
    if( err ) return err;

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)4.0,*A);
        if( err ) return err;

        if( i>0 )
        {
            err = lis_matrix_set_value(LIS_INS_VALUE,i,i-1,(LIS_SCALAR)-0.5,*A);
            if( err ) return err;
        }

        if( i<n-1 )
        {
            err = lis_matrix_set_value(LIS_INS_VALUE,i,i+1,(LIS_SCALAR)1.0,*A);
            if( err ) return err;
        }
    }

    err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
    if( err ) return err;

    return lis_matrix_assemble(*A);
}

static int
check_solution(LIS_VECTOR x, LIS_REAL tol, const char *label)
{
    LIS_INT i;
    LIS_REAL e;

    for(i=0;i<x->n;i++)
    {
        e = (LIS_REAL)fabs(x->value[i]-(LIS_SCALAR)1.0);
        if( e>tol )
        {
            fprintf(stderr,"%s solution mismatch at %d: error=%e\n",
                    label,(int)i,(double)e);
            return 1;
        }
    }

    return 0;
}

static int
run_nested_success(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x)
{
    LIS_SOLVER solver = NULL;
    LIS_INT err;
    int failed = 0;

    err = lis_solver_create(&solver);
    if( err ) return 1;

    err = lis_solver_set_option(
        "-i fgmres -p hybrid "
        "-hybrid_i gmres -hybrid_p jacobi "
        "-hybrid_maxiter 100 -hybrid_restart 32 "
        "-hybrid_tol 1.0e-10 "
        "-restart 16 -maxiter 80 -tol 1.0e-9 "
        "-scale none -print none",
        solver);
    if( err )
    {
        fprintf(stderr,"nested option setup failed: %d\n",(int)err);
        failed = 1;
        goto cleanup;
    }

    lis_vector_set_all((LIS_SCALAR)0.0,x);
    err = lis_solve(A,b,x,solver);
    if( err || solver->retcode!=LIS_SUCCESS )
    {
        fprintf(stderr,"nested solve failed: err=%d retcode=%d\n",
                (int)err,(int)solver->retcode);
        failed = 1;
        goto cleanup;
    }

    failed |= check_solution(x,(LIS_REAL)1.0e-7,"FGMRES nested Hybrid");

cleanup:
    if( solver ) lis_solver_destroy(solver);
    return failed;
}

static int
run_nested_mfamg_success(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR exact,
    LIS_VECTOR x)
{
    LIS_SOLVER solver = NULL;
    LIS_VECTOR near_nullspace[1];
    LIS_INT err;
    int failed = 0;

    near_nullspace[0] = exact;

    err = lis_solver_create(&solver);
    if( err ) return 1;

    err = lis_solver_set_option(
        "-i fgmres -p hybrid "
        "-hybrid_i gmres -hybrid_p mfamg "
        "-hybrid_maxiter 100 -hybrid_restart 32 "
        "-hybrid_tol 1.0e-9 "
        "-restart 16 -maxiter 80 -tol 1.0e-8 "
        "-scale none -print none",
        solver);
    if( err )
    {
        fprintf(stderr,"nested MFAMG option setup failed: %d\n",(int)err);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_near_nullspace(solver,1,near_nullspace);
    if( err )
    {
        fprintf(stderr,"nested MFAMG near-nullspace setup failed: %d\n",(int)err);
        failed = 1;
        goto cleanup;
    }

    lis_vector_set_all((LIS_SCALAR)0.0,x);
    err = lis_solve(A,b,x,solver);
    if( err || solver->retcode!=LIS_SUCCESS )
    {
        fprintf(stderr,"nested MFAMG solve failed: err=%d retcode=%d\n",
                (int)err,(int)solver->retcode);
        failed = 1;
        goto cleanup;
    }

    failed |= check_solution(x,(LIS_REAL)1.0e-6,
                             "FGMRES Hybrid GMRES MFAMG");

cleanup:
    if( solver ) lis_solver_destroy(solver);
    return failed;
}

static int
run_nested_inner_gcrodr_mfamg_success(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR exact,
    LIS_VECTOR x)
{
    LIS_SOLVER solver = NULL;
    LIS_VECTOR near_nullspace[1];
    LIS_INT err;
    int failed = 0;

    near_nullspace[0] = exact;

    err = lis_solver_create(&solver);
    if( err ) return 1;

    /*
     * Exercise GCRO-DR as the inner Krylov solver:
     *
     *   FGMRES -> Hybrid -> GCRO-DR -> MFAMG
     *
     * The Hybrid-owned inner solver object persists across preconditioner
     * applications, so this also exercises GCRO-DR state lifetime inside
     * the nested composition.
     */
    err = lis_solver_set_option(
        "-i fgmres -p hybrid "
        "-hybrid_i gcrodr -hybrid_p mfamg "
        "-hybrid_maxiter 64 -hybrid_restart 16 "
        "-hybrid_tol 1.0e-9 "
        "-gcrodr_keep 0 "
        "-restart 16 -maxiter 80 -tol 1.0e-8 "
        "-scale none -print none",
        solver);
    if( err )
    {
        fprintf(stderr,
                "nested inner GCRO-DR option setup failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_near_nullspace(
        solver,1,near_nullspace);
    if( err )
    {
        fprintf(stderr,
                "nested inner GCRO-DR near-nullspace setup failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    lis_vector_set_all((LIS_SCALAR)0.0,x);

    err = lis_solve(A,b,x,solver);
    if( err || solver->retcode!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "nested inner GCRO-DR solve failed: "
                "err=%d retcode=%d\n",
                (int)err,
                (int)solver->retcode);
        failed = 1;
        goto cleanup;
    }

    failed |= check_solution(
        x,
        (LIS_REAL)1.0e-6,
        "FGMRES Hybrid GCRO-DR MFAMG");

cleanup:
    if( solver ) lis_solver_destroy(solver);
    return failed;
}


static int
run_outer_fgcrodr_hybrid_mfamg_success(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR exact,
    LIS_VECTOR x)
{
    LIS_SOLVER solver = NULL;
    LIS_VECTOR near_nullspace[1];
    LIS_INT err;
    int failed = 0;

    near_nullspace[0] = exact;

    err = lis_solver_create(&solver);
    if( err ) return 1;

    /*
     * Exercise FGCRO-DR as the outer flexible Krylov solver:
     *
     *   FGCRO-DR -> Hybrid -> GMRES -> MFAMG
     *
     * Run the same solver twice so the second solve also exercises the
     * persistent GCRO-DR recycle state together with Hybrid.
     */
    err = lis_solver_set_option(
        "-i fgcrodr -p hybrid "
        "-hybrid_i gmres -hybrid_p mfamg "
        "-hybrid_maxiter 64 -hybrid_restart 16 "
        "-hybrid_tol 1.0e-9 "
        "-fgcrodr_k 4 "
        "-fgcrodr_keep 0.90 "
        "-fgcrodr_reject 0.35 "
        "-fgcrodr_update true "
        "-restart 16 -maxiter 80 -tol 1.0e-8 "
        "-scale none -print none",
        solver);
    if( err )
    {
        fprintf(stderr,
                "outer FGCRO-DR option setup failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_near_nullspace(
        solver,1,near_nullspace);
    if( err )
    {
        fprintf(stderr,
                "outer FGCRO-DR near-nullspace setup failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    /*
     * First solve.
     */
    lis_vector_set_all((LIS_SCALAR)0.0,x);

    err = lis_solve(A,b,x,solver);
    if( err || solver->retcode!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "outer FGCRO-DR first solve failed: "
                "err=%d retcode=%d\n",
                (int)err,
                (int)solver->retcode);
        failed = 1;
        goto cleanup;
    }

    failed |= check_solution(
        x,
        (LIS_REAL)1.0e-6,
        "FGCRO-DR Hybrid GMRES MFAMG first");

    if( failed ) goto cleanup;

    /*
     * Second solve on the same solver object.
     * This exercises persistent GCRO-DR recycling.
     */
    lis_vector_set_all((LIS_SCALAR)0.0,x);

    err = lis_solve(A,b,x,solver);
    if( err || solver->retcode!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "outer FGCRO-DR recycled solve failed: "
                "err=%d retcode=%d\n",
                (int)err,
                (int)solver->retcode);
        failed = 1;
        goto cleanup;
    }

    failed |= check_solution(
        x,
        (LIS_REAL)1.0e-6,
        "FGCRO-DR Hybrid GMRES MFAMG recycled");

cleanup:
    if( solver ) lis_solver_destroy(solver);
    return failed;
}


static int
run_adjoint_rejection(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x)
{
    LIS_SOLVER solver = NULL;
    LIS_PRECON precon = NULL;
    LIS_INT err;
    int failed = 0;

    err = lis_solver_create(&solver);
    if( err ) return 1;

    err = lis_solver_set_option(
        "-i gmres -p hybrid "
        "-hybrid_i gmres -hybrid_p jacobi "
        "-hybrid_maxiter 32 -hybrid_restart 16 "
        "-hybrid_tol 1.0e-10 "
        "-scale none -print none",
        solver);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    solver->A = A;
    err = lis_precon_create(solver,&precon);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    solver->precon = precon;
    lis_vector_set_all((LIS_SCALAR)0.0,x);
    err = lis_psolveh(solver,b,x);

    if( err!=LIS_ERR_NOT_IMPLEMENTED )
    {
        fprintf(stderr,
                "Hybrid adjoint returned %d, expected %d\n",
                (int)err,(int)LIS_ERR_NOT_IMPLEMENTED);
        failed = 1;
    }

cleanup:
    if( precon )
    {
        lis_precon_destroy(precon);
        precon = NULL;
        if( solver ) solver->precon = NULL;
    }
    if( solver ) lis_solver_destroy(solver);
    return failed;
}

static int
run_recursive_rejection(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x)
{
    LIS_SOLVER solver = NULL;
    LIS_INT err;
    int failed = 0;

    err = lis_solver_create(&solver);
    if( err ) return 1;

    err = lis_solver_set_option(
        "-i gmres -p hybrid "
        "-hybrid_i gmres -hybrid_p hybrid "
        "-hybrid_maxiter 4 -hybrid_restart 4 "
        "-hybrid_tol 1.0e-6 "
        "-restart 4 -maxiter 10 -tol 1.0e-8 "
        "-scale none -print none",
        solver);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    lis_vector_set_all((LIS_SCALAR)0.0,x);
    err = lis_solve(A,b,x,solver);

    if( err!=LIS_ERR_ILL_ARG )
    {
        fprintf(stderr,
                "recursive Hybrid returned %d, expected %d\n",
                (int)err,(int)LIS_ERR_ILL_ARG);
        failed = 1;
    }

cleanup:
    if( solver ) lis_solver_destroy(solver);
    return failed;
}

int main(int argc, char **argv)
{
    LIS_MATRIX A = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR exact = NULL;
    LIS_INT err;
    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if( err ) return 1;

    err = create_nonsymmetric_matrix(&A,32);
    if( err ) { failed = 1; goto cleanup; }

    err = lis_vector_duplicate(A,&b);
    if( err ) { failed = 1; goto cleanup; }

    err = lis_vector_duplicate(A,&x);
    if( err ) { failed = 1; goto cleanup; }

    err = lis_vector_duplicate(A,&exact);
    if( err ) { failed = 1; goto cleanup; }

    lis_vector_set_all((LIS_SCALAR)1.0,exact);
    err = lis_matvec(A,exact,b);
    if( err ) { failed = 1; goto cleanup; }

    failed |= run_nested_success(A,b,x);
    if( failed ) goto cleanup;

    /* MFAMG must also be composable as the inner preconditioner. */
    failed |= run_nested_mfamg_success(A,b,exact,x);
    if( failed ) goto cleanup;

    /*
     * GCRO-DR must work both as an inner Hybrid KSP and as an outer
     * recycled Krylov solver using Hybrid as its preconditioner.
     */
    failed |= run_nested_inner_gcrodr_mfamg_success(
        A,b,exact,x);
    if( failed ) goto cleanup;

    failed |= run_outer_fgcrodr_hybrid_mfamg_success(
        A,b,exact,x);
    if( failed ) goto cleanup;

    /* The legacy adjoint path was mathematically wrong for nonsymmetric A. */
    failed |= run_adjoint_rejection(A,b,x);
    if( failed ) goto cleanup;

    /* Normal nested KSP/PC operation must remain unaffected. */
    failed |= run_nested_success(A,b,x);
    if( failed ) goto cleanup;

    failed |= run_recursive_rejection(A,b,x);
    if( failed ) goto cleanup;

    printf("LIS_NESTED_KSP_PC PASSED\n");

cleanup:
    if( exact ) lis_vector_destroy(exact);
    if( x ) lis_vector_destroy(x);
    if( b ) lis_vector_destroy(b);
    if( A ) lis_matrix_destroy(A);

    lis_finalize();
    return failed ? 1 : 0;
}
