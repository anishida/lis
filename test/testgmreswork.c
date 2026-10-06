#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_solver.h"

#include <stdio.h>


static LIS_INT
create_matrix(LIS_MATRIX *A)
{
    LIS_INT err;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,2);
    if( err ) return err;

    err = lis_matrix_set_value(
        LIS_INS_VALUE,
        0,0,
        (LIS_SCALAR)2.0,
        *A);
    if( err ) return err;

    err = lis_matrix_set_value(
        LIS_INS_VALUE,
        1,1,
        (LIS_SCALAR)3.0,
        *A);
    if( err ) return err;

    err = lis_matrix_set_type(
        *A,
        LIS_MATRIX_CSR);
    if( err ) return err;

    return lis_matrix_assemble(*A);
}


static int
check_failed_work_state(
    const char *name,
    LIS_SOLVER solver,
    LIS_INT err)
{
    int failed = 0;

    printf(
        "%s returned=%d work=%p worklen=%d\n",
        name,
        (int)err,
        (void *)solver->work,
        (int)solver->worklen);

    if( err!=LIS_ERR_ILL_ARG )
    {
        fprintf(
            stderr,
            "%s returned %d, expected LIS_ERR_ILL_ARG=%d\n",
            name,
            (int)err,
            (int)LIS_ERR_ILL_ARG);

        failed = 1;
    }

    if( solver->work!=NULL )
    {
        fprintf(
            stderr,
            "%s published work after allocation failure\n",
            name);

        failed = 1;
    }

    if( solver->worklen!=0 )
    {
        fprintf(
            stderr,
            "%s worklen=%d, expected 0\n",
            name,
            (int)solver->worklen);

        failed = 1;
    }

    return failed;
}


int
main(int argc, char **argv)
{
    LIS_MATRIX A = NULL;
    LIS_SOLVER solver = NULL;
    LIS_INT err;
    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if( err ) return 1;

    err = create_matrix(&A);
    if( err )
    {
        fprintf(
            stderr,
            "create_matrix failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_solver_create(&solver);
    if( err )
    {
        fprintf(
            stderr,
            "lis_solver_create failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_matrix(A,solver);
    if( err )
    {
        fprintf(
            stderr,
            "lis_solver_set_matrix failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    /*
     * This test calls malloc_work directly on purpose.
     *
     * Normal solve parameter validation rejects negative restart
     * before malloc_work is entered.  Here restart=-2 is only a
     * deterministic, allocation-safe way to make the special
     * work[0] vector call:
     *
     *     lis_vector_set_size(work[0], -1, 0)
     *
     * which returns LIS_ERR_ILL_ARG.
     *
     * malloc_work must propagate that error and must not publish a
     * partially initialized solver->work array.
     */


    /* -------------------------------------------------------- */
    /* GMRES                                                    */
    /* NWORK=4, restart=-2 -> worklen=3, work[0] size=-1       */
    /* -------------------------------------------------------- */

    solver->options[LIS_OPTIONS_RESTART] = -2;

    err = lis_gmres_malloc_work(solver);

    if( check_failed_work_state(
            "GMRES malloc_work",
            solver,
            err) )
    {
        failed = 1;
    }

    /*
     * Current RED implementation publishes work even though
     * lis_vector_set_size() failed.  Clean it so the FGMRES part
     * can run safely.  After the fix this is a no-op.
     */
    if( solver->work )
    {
        lis_solver_work_destroy(solver);
    }


    /* -------------------------------------------------------- */
    /* FGMRES                                                   */
    /* NWORK=4, restart=-2 -> worklen=1, work[0] size=-1       */
    /* -------------------------------------------------------- */

    solver->options[LIS_OPTIONS_RESTART] = -2;

    err = lis_fgmres_malloc_work(solver);

    if( check_failed_work_state(
            "FGMRES malloc_work",
            solver,
            err) )
    {
        failed = 1;
    }

    /*
     * Same RED cleanup for the current implementation.
     */
    if( solver->work )
    {
        lis_solver_work_destroy(solver);
    }


cleanup:

    if( solver )
    {
        lis_solver_destroy(solver);
        solver = NULL;
    }

    if( A )
    {
        lis_matrix_destroy(A);
        A = NULL;
    }

    lis_finalize();

    if( failed )
    {
        fprintf(
            stderr,
            "LIS_GMRES_WORK_STAGE5A FAILED\n");

        return 1;
    }

    printf(
        "LIS_GMRES_WORK_STAGE5A PASSED\n");

    return 0;
}
