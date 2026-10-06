#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_solver.h"

#include <stdio.h>


static LIS_INT
create_diag_matrix(
    LIS_INT n,
    LIS_SCALAR diag,
    LIS_MATRIX *A)
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
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            diag,
            *A);

        if( err ) return err;
    }

    err = lis_matrix_set_type(
        *A,
        LIS_MATRIX_CSR);

    if( err ) return err;

    return lis_matrix_assemble(*A);
}


static int
check_solve(
    const char *name,
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_SOLVER solver)
{
    LIS_INT err;

    lis_vector_set_all(0.0,x);

    err = lis_solve(
        A,
        b,
        x,
        solver);

    if( err )
    {
        fprintf(
            stderr,
            "%s: lis_solve API error=%d\n",
            name,
            (int)err);

        return 1;
    }

    if( solver->retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "%s: solver retcode=%d\n",
            name,
            (int)solver->retcode);

        return 1;
    }

    if( solver->work==NULL )
    {
        fprintf(
            stderr,
            "%s: workspace was destroyed after successful solve\n",
            name);

        return 1;
    }

    if( !lis_solver_work_compatible(solver) )
    {
        fprintf(
            stderr,
            "%s: surviving workspace is incompatible\n",
            name);

        return 1;
    }

    return 0;
}


int
main(int argc, char **argv)
{
    LIS_MATRIX A = NULL;
    LIS_MATRIX B = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_SOLVER solver = NULL;

    LIS_VECTOR *saved_work = NULL;
    LIS_VECTOR saved_work0 = NULL;
    LIS_VECTOR saved_work1 = NULL;

    LIS_INT err;
    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if( err ) return 1;


    /*
     * Same distribution, different matrix objects and values.
     */
    err = create_diag_matrix(
        32,
        (LIS_SCALAR)2.0,
        &A);

    if( err )
    {
        fprintf(stderr,
                "create A failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }


    err = create_diag_matrix(
        32,
        (LIS_SCALAR)3.0,
        &B);

    if( err )
    {
        fprintf(stderr,
                "create B failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }


    err = lis_vector_duplicate(A,&b);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_duplicate(A,&x);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_solver_create(&solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_GMRES;

    solver->options[LIS_OPTIONS_PRECON] =
        LIS_PRECON_TYPE_NONE;

    solver->options[LIS_OPTIONS_SCALE] =
        LIS_SCALE_NONE;

    solver->options[LIS_OPTIONS_RESTART] = 5;

    solver->options[LIS_OPTIONS_MAXITER] = 100;

    solver->options[LIS_OPTIONS_OUTPUT] = 0;


    /*
     * Solve #1.
     *
     * GMRES:
     * NWORK + (restart+1) = 4 + 6 = 10.
     *
     * Before Stage 5C solver->work was destroyed at the
     * successful end of lis_solve_kernel().
     */
    lis_vector_set_all(
        (LIS_SCALAR)2.0,
        b);

    if( check_solve(
            "solve1",
            A,
            b,
            x,
            solver) )
    {
        failed = 1;
        goto cleanup;
    }


    if( solver->worklen!=10 )
    {
        fprintf(
            stderr,
            "solve1: expected worklen=10, got %d\n",
            (int)solver->worklen);

        failed = 1;
        goto cleanup;
    }


    if( solver->work[0]->n!=6 )
    {
        fprintf(
            stderr,
            "solve1: expected GMRES work[0]->n=6, got %d\n",
            (int)solver->work[0]->n);

        failed = 1;
        goto cleanup;
    }


    saved_work  = solver->work;
    saved_work0 = solver->work[0];
    saved_work1 = solver->work[1];


    /*
     * Solve #2: exactly the same layout/options.
     * Workspace must be reusable.
     */
    if( check_solve(
            "solve2",
            A,
            b,
            x,
            solver) )
    {
        failed = 1;
        goto cleanup;
    }


    if( solver->work!=saved_work ||
        solver->work[0]!=saved_work0 ||
        solver->work[1]!=saved_work1 )
    {
        fprintf(
            stderr,
            "solve2: compatible workspace was not retained\n");

        failed = 1;
        goto cleanup;
    }


    /*
     * Solve #3: different matrix object and numerical values,
     * but identical distribution/layout.
     *
     * Matrix pointer identity must not force reallocation.
     */
    lis_vector_set_all(
        (LIS_SCALAR)3.0,
        b);

    if( check_solve(
            "solve3",
            B,
            b,
            x,
            solver) )
    {
        failed = 1;
        goto cleanup;
    }


    if( solver->work!=saved_work ||
        solver->work[0]!=saved_work0 ||
        solver->work[1]!=saved_work1 )
    {
        fprintf(
            stderr,
            "solve3: equivalent matrix layout did not reuse workspace\n");

        failed = 1;
        goto cleanup;
    }


    /*
     * Solve #4: restart changes from 5 to 6.
     *
     * Existing workspace is structurally incompatible and
     * must be replaced with:
     *
     *   worklen = 4 + (6+1) = 11
     *   work[0]->n = 7
     *
     * Do not compare pointer inequality here: an allocator is
     * allowed to return the same address after free+malloc.
     */
    solver->options[LIS_OPTIONS_RESTART] = 6;

    if( check_solve(
            "solve4",
            B,
            b,
            x,
            solver) )
    {
        failed = 1;
        goto cleanup;
    }


    if( solver->worklen!=11 )
    {
        fprintf(
            stderr,
            "solve4: expected worklen=11, got %d\n",
            (int)solver->worklen);

        failed = 1;
    }


    if( solver->work[0]->n!=7 )
    {
        fprintf(
            stderr,
            "solve4: expected GMRES work[0]->n=7, got %d\n",
            (int)solver->work[0]->n);

        failed = 1;
    }


cleanup:

    if( solver )
    {
        lis_solver_destroy(solver);
    }

    if( x )
    {
        lis_vector_destroy(x);
    }

    if( b )
    {
        lis_vector_destroy(b);
    }

    if( B )
    {
        lis_matrix_destroy(B);
    }

    if( A )
    {
        lis_matrix_destroy(A);
    }

    lis_finalize();


    if( failed )
    {
        fprintf(
            stderr,
            "LIS_WORK_REUSE_STAGE5C FAILED\n");

        return 1;
    }


    printf(
        "LIS_WORK_REUSE_STAGE5C PASSED\n");

    return 0;
}
