#include <stdio.h>
#include "lis.h"

static int create_count = 0;
static int destroy_count = 0;
static int psolve_count = 0;

static LIS_INT create_threads_seen = -1;
static LIS_INT psolve_threads_seen = -1;


static LIS_INT
thread_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    (void)solver;
    (void)precon;

    create_count++;

#ifdef _OPENMP
    create_threads_seen = (LIS_INT)omp_get_max_threads();
#else
    create_threads_seen = 1;
#endif

    return LIS_SUCCESS;
}


static LIS_INT
thread_psolve(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    (void)solver;

    psolve_count++;

#ifdef _OPENMP
    psolve_threads_seen = (LIS_INT)omp_get_max_threads();
#else
    psolve_threads_seen = 1;
#endif

    return lis_vector_copy(b,x);
}


static LIS_INT
thread_destroy(LIS_PRECON precon)
{
    (void)precon;

    destroy_count++;

    return LIS_SUCCESS;
}


static LIS_INT
create_matrix(LIS_MATRIX *A)
{
    LIS_INT err;
    LIS_INT i,is,ie;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,2);
    if( err ) return err;

    err = lis_matrix_get_range(*A,&is,&ie);
    if( err ) return err;

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,i,
            (LIS_SCALAR)(i==0 ? 2.0 : 3.0),
            *A);

        if( err ) return err;
    }

    err = lis_matrix_set_type(
        *A,
        LIS_MATRIX_CSR);

    if( err ) return err;

    return lis_matrix_assemble(*A);
}


static LIS_INT
create_solver(LIS_SOLVER *solver)
{
    LIS_INT err;

    err = lis_solver_create(solver);
    if( err ) return err;

    return lis_solver_set_option(
        "-i cg -p thrs9 "
        "-scale none -print none "
        "-conv_cond nrm2_b -maxiter 1 "
        "-maxiter_noimp 0 -tol 0.25",
        *solver);
}


static void
reset_observed(void)
{
    create_threads_seen = -1;
    psolve_threads_seen = -1;
    psolve_count = 0;
}


static LIS_INT
runtime_threads(void)
{
#ifdef _OPENMP
    return (LIS_INT)omp_get_max_threads();
#else
    return 1;
#endif
}


static int
check_runtime(
    const char *name,
    LIS_INT expected)
{
    LIS_INT got = runtime_threads();

    if( got!=expected )
    {
        fprintf(
            stderr,
            "%s runtime=%d expected=%d\n",
            name,
            (int)got,
            (int)expected);
        return 1;
    }

    return 0;
}


static int
check_seen(
    const char *name,
    LIS_INT expected,
    LIS_INT require_create)
{
    if( require_create &&
        create_threads_seen!=expected )
    {
        fprintf(
            stderr,
            "%s create_threads=%d expected=%d\n",
            name,
            (int)create_threads_seen,
            (int)expected);
        return 1;
    }

    if( psolve_count<=0 ||
        psolve_threads_seen!=expected )
    {
        fprintf(
            stderr,
            "%s psolve_count=%d threads=%d expected=%d\n",
            name,
            psolve_count,
            (int)psolve_threads_seen,
            (int)expected);
        return 1;
    }

    return 0;
}


int main(int argc, char **argv)
{
    LIS_MATRIX A = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT value;
    LIS_INT runtime_default;
    LIS_INT second_threads;

    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if( err ) return 1;

    runtime_default = runtime_threads();

    printf(
        "THREAD_RUNTIME_DEFAULT=%d\n",
        (int)runtime_default);

    err = lis_precon_register_ex(
        "thrs9",
        thread_create,
        thread_psolve,
        thread_psolve,
        thread_destroy);

    if( err )
    {
        fprintf(
            stderr,
            "thread preconditioner registration failed rc=%d\n",
            (int)err);
        failed = 1;
        goto cleanup;
    }

    err = create_matrix(&A);
    if( err )
    {
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

    lis_vector_set_all(
        (LIS_SCALAR)2.0,
        b);


    /* ======================================================== */
    /* Public API defaults and validation.                      */
    /* ======================================================== */

    err = create_solver(&solver);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    value = -1;

    err = lis_solver_get_num_threads(
        solver,&value);

    if( err!=LIS_SUCCESS || value!=0 )
    {
        fprintf(
            stderr,
            "DEFAULT getter rc=%d value=%d expected=0\n",
            (int)err,
            (int)value);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_num_threads(
        solver,-1);

    if( err!=LIS_ERR_ILL_ARG )
    {
        fprintf(
            stderr,
            "NEGATIVE rc=%d expected=%d\n",
            (int)err,
            (int)LIS_ERR_ILL_ARG);
        failed = 1;
        goto cleanup;
    }

#ifdef _OPENMP
    err = lis_solver_set_num_threads(
        solver,runtime_default+1);

    if( err!=LIS_ERR_ILL_ARG )
    {
        fprintf(
            stderr,
            "OVER_CAPACITY rc=%d expected=%d "
            "capacity=%d requested=%d\n",
            (int)err,
            (int)LIS_ERR_ILL_ARG,
            (int)runtime_default,
            (int)(runtime_default+1));
        failed = 1;
        goto cleanup;
    }
#else
    err = lis_solver_set_num_threads(
        solver,2);

    if( err!=LIS_ERR_NOT_IMPLEMENTED )
    {
        fprintf(
            stderr,
            "NO_OPENMP rc=%d expected=%d\n",
            (int)err,
            (int)LIS_ERR_NOT_IMPLEMENTED);
        failed = 1;
        goto cleanup;
    }
#endif

    printf(
        "THREAD_API_VALIDATION PASSED default=0 capacity=%d\n",
        (int)runtime_default);


#ifdef _OPENMP
    /*
     * num_threads == 0 inherits the current OpenMP setting.
     * A runtime setting above the capacity established by
     * lis_initialize() must therefore be rejected before any
     * solver work starts.
     */
    omp_set_num_threads(runtime_default+1);

    if( runtime_threads()!=runtime_default+1 )
    {
        fprintf(
            stderr,
            "RUNTIME_OVER_CAPACITY setup=%d expected=%d\n",
            (int)runtime_threads(),
            (int)(runtime_default+1));

        omp_set_num_threads(runtime_default);

        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_num_threads(
        solver,0);

    if( err )
    {
        omp_set_num_threads(runtime_default);

        failed = 1;
        goto cleanup;
    }

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    err = lis_solve(
        A,b,x,solver);

    omp_set_num_threads(runtime_default);

    if( err!=LIS_ERR_ILL_ARG )
    {
        fprintf(
            stderr,
            "RUNTIME_OVER_CAPACITY rc=%d expected=%d\n",
            (int)err,
            (int)LIS_ERR_ILL_ARG);

        failed = 1;
        goto cleanup;
    }

    if( check_runtime(
            "RUNTIME_OVER_CAPACITY_RESTORE",
            runtime_default) )
    {
        failed = 1;
        goto cleanup;
    }

    printf(
        "RUNTIME_OVER_CAPACITY PASSED "
        "requested=%d capacity=%d restored=%d\n",
        (int)(runtime_default+1),
        (int)runtime_default,
        (int)runtime_threads());
#endif


    /* ======================================================== */
    /* Legacy solve: explicit one-thread override.              */
    /* ======================================================== */

    err = lis_solver_set_num_threads(
        solver,1);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    reset_observed();

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    err = lis_solve(
        A,b,x,solver);

    if( err!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "LEGACY_ONE failed api=%d status=%d\n",
            (int)err,
            (int)solver->retcode);
        failed = 1;
        goto cleanup;
    }

    if( check_seen(
            "LEGACY_ONE",
            1,
            LIS_TRUE) )
    {
        failed = 1;
        goto cleanup;
    }

    if( check_runtime(
            "LEGACY_ONE_RESTORE",
            runtime_default) )
    {
        failed = 1;
        goto cleanup;
    }

    printf(
        "LEGACY_ONE PASSED create_threads=%d "
        "psolve_threads=%d restored=%d\n",
        (int)create_threads_seen,
        (int)psolve_threads_seen,
        (int)runtime_threads());


    /* ======================================================== */
    /* Legacy solve: zero means inherit current LIS runtime.    */
    /* ======================================================== */

    err = lis_solver_set_num_threads(
        solver,0);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    reset_observed();

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    err = lis_solve(
        A,b,x,solver);

    if( err!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "LEGACY_INHERIT failed api=%d status=%d\n",
            (int)err,
            (int)solver->retcode);
        failed = 1;
        goto cleanup;
    }

    if( check_seen(
            "LEGACY_INHERIT",
            runtime_default,
            LIS_TRUE) )
    {
        failed = 1;
        goto cleanup;
    }

    if( check_runtime(
            "LEGACY_INHERIT_RESTORE",
            runtime_default) )
    {
        failed = 1;
        goto cleanup;
    }

    printf(
        "LEGACY_INHERIT PASSED threads=%d\n",
        (int)runtime_default);


    /* ======================================================== */
    /* Legacy error path must also restore runtime.             */
    /* ======================================================== */

    if( runtime_default>1 )
    {
        err = lis_solver_set_num_threads(
            solver,1);

        if( err )
        {
            failed = 1;
            goto cleanup;
        }

        solver->options[LIS_OPTIONS_PRECON] = -1;

        err = lis_solve(
            A,b,x,solver);

        if( err!=LIS_ERR_ILL_ARG )
        {
            fprintf(
                stderr,
                "LEGACY_ERROR rc=%d expected=%d\n",
                (int)err,
                (int)LIS_ERR_ILL_ARG);
            failed = 1;
            goto cleanup;
        }

        if( check_runtime(
                "LEGACY_ERROR_RESTORE",
                runtime_default) )
        {
            failed = 1;
            goto cleanup;
        }

        printf(
            "LEGACY_ERROR_RESTORE PASSED restored=%d\n",
            (int)runtime_threads());

        lis_solver_destroy(solver);
        solver = NULL;

        err = create_solver(&solver);
        if( err )
        {
            failed = 1;
            goto cleanup;
        }
    }


    /* ======================================================== */
    /* Sequence solve: thread changes must not rebuild precon.  */
    /* ======================================================== */

    err = lis_solver_set_num_threads(
        solver,1);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    create_count = 0;
    destroy_count = 0;
    reset_observed();

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    err = lis_solve_sequence(
        A,b,x,solver);

    if( err!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "SEQUENCE_ONE failed api=%d status=%d\n",
            (int)err,
            (int)solver->retcode);
        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        check_seen(
            "SEQUENCE_ONE",
            1,
            LIS_TRUE) ||
        check_runtime(
            "SEQUENCE_ONE_RESTORE",
            runtime_default) )
    {
        fprintf(
            stderr,
            "SEQUENCE_ONE lifecycle create=%d destroy=%d\n",
            create_count,
            destroy_count);
        failed = 1;
        goto cleanup;
    }

    printf(
        "SEQUENCE_ONE PASSED create=%d "
        "threads=%d restored=%d\n",
        create_count,
        (int)psolve_threads_seen,
        (int)runtime_threads());


#ifdef _OPENMP
    if( runtime_default>=2 )
    {
        second_threads = 2;
    }
    else
    {
        second_threads = 1;
    }
#else
    second_threads = 1;
#endif

    err = lis_solver_set_num_threads(
        solver,
        second_threads);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    reset_observed();

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    err = lis_solve_sequence(
        A,b,x,solver);

    if( err!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "SEQUENCE_CHANGE failed api=%d status=%d\n",
            (int)err,
            (int)solver->retcode);
        failed = 1;
        goto cleanup;
    }

    /*
     * Changing only the thread override must not invalidate the
     * persistent preconditioner.
     */
    if( create_count!=1 )
    {
        fprintf(
            stderr,
            "SEQUENCE_CHANGE rebuilt preconditioner create=%d\n",
            create_count);
        failed = 1;
        goto cleanup;
    }

    if( check_seen(
            "SEQUENCE_CHANGE",
            second_threads,
            LIS_FALSE) )
    {
        failed = 1;
        goto cleanup;
    }

    if( check_runtime(
            "SEQUENCE_CHANGE_RESTORE",
            runtime_default) )
    {
        failed = 1;
        goto cleanup;
    }

    printf(
        "SEQUENCE_CHANGE PASSED threads=%d "
        "create=%d destroy=%d restored=%d\n",
        (int)second_threads,
        create_count,
        destroy_count,
        (int)runtime_threads());


    /* ======================================================== */
    /* Reset must preserve thread configuration itself.         */
    /* ======================================================== */

    err = lis_solver_reset_sequence(
        solver);

    if( err )
    {
        fprintf(
            stderr,
            "SEQUENCE_RESET failed rc=%d\n",
            (int)err);
        failed = 1;
        goto cleanup;
    }

    value = -1;

    err = lis_solver_get_num_threads(
        solver,&value);

    if( err ||
        value!=second_threads )
    {
        fprintf(
            stderr,
            "SEQUENCE_RESET thread setting=%d expected=%d\n",
            (int)value,
            (int)second_threads);
        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        destroy_count!=1 )
    {
        fprintf(
            stderr,
            "SEQUENCE_RESET lifecycle create=%d destroy=%d\n",
            create_count,
            destroy_count);
        failed = 1;
        goto cleanup;
    }

    printf(
        "SEQUENCE_RESET PASSED setting=%d "
        "create=%d destroy=%d\n",
        (int)value,
        create_count,
        destroy_count);


    printf(
        "LIS_SOLVER_THREADS_S9 PASSED\n");


cleanup:

    if( solver!=NULL )
    {
        lis_solver_destroy(solver);
        solver = NULL;
    }

    if( x!=NULL )
    {
        lis_vector_destroy(x);
        x = NULL;
    }

    if( b!=NULL )
    {
        lis_vector_destroy(b);
        b = NULL;
    }

    if( A!=NULL )
    {
        lis_matrix_destroy(A);
        A = NULL;
    }

    lis_finalize();

    return failed ? 1 : 0;
}
