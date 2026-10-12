#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_solver.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 64

typedef struct
{
    LIS_INT api;
    LIS_INT status;
    LIS_INT iter;
    LIS_REAL reported;
    LIS_REAL true_residual;

    LIS_INT recycle_dim;
    LIS_INT recycle_action;

    LIS_REAL recycle_quality;
    LIS_REAL recycle_quality_min;

    LIS_REAL recycle_orth_error;
    LIS_REAL recycle_relation_error;
} RESULT;


static LIS_INT
build_matrix(LIS_MATRIX *Aout)
{
    LIS_MATRIX A = NULL;
    LIS_INT i;
    LIS_INT err;

    err = lis_matrix_create(LIS_COMM_WORLD,&A);
    if( err ) return err;

    err = lis_matrix_set_size(A,0,N);
    if( err )
    {
        lis_matrix_destroy(A);
        return err;
    }

    /*
     * Mildly nonsymmetric tridiagonal system.
     */
    for(i=0;i<N;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,i,
            (LIS_SCALAR)4.0,
            A);
        if( err ) goto fail;

        if( i>0 )
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,
                i,i-1,
                (LIS_SCALAR)-1.0,
                A);
            if( err ) goto fail;
        }

        if( i<N-1 )
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,
                i,i+1,
                (LIS_SCALAR)-0.5,
                A);
            if( err ) goto fail;
        }
    }

    err = lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if( err ) goto fail;

    err = lis_matrix_assemble(A);
    if( err ) goto fail;

    *Aout = A;
    return LIS_SUCCESS;

fail:
    lis_matrix_destroy(A);
    return err;
}


static void
set_exact_solution(LIS_VECTOR u)
{
    LIS_INT i;

    for(i=0;i<u->n;i++)
    {
        LIS_REAL g = (LIS_REAL)(u->is+i+1);

        u->value[i] =
            (LIS_SCALAR)(
                1.0
                + 0.15*sin(0.071*g)
                + 0.04*cos(0.037*g));
    }
}


static LIS_REAL
true_relative_residual(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_VECTOR r)
{
    LIS_REAL rnrm = 0.0;
    LIS_REAL bnrm = 0.0;

    lis_matvec(A,x,r);
    lis_vector_axpy((LIS_SCALAR)-1.0,b,r);

    lis_vector_nrm2(r,&rnrm);
    lis_vector_nrm2(b,&bnrm);

    if( bnrm!=0.0 )
        return rnrm/bnrm;

    return rnrm;
}


static LIS_INT
run_solver(
    const char *method,
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_VECTOR r,
    RESULT *result)
{
    LIS_SOLVER solver = NULL;
    LIS_INT err;
    char option[512];

    memset(result,0,sizeof(*result));

    result->api = -999;
    result->status = -999;
    result->iter = -999;
    result->reported = -1.0;
    result->true_residual = -1.0;

    err = lis_solver_create(&solver);
    if( err ) return err;

    snprintf(
        option,
        sizeof(option),
        "-i %s -p jacobi "
        "-restart 10 -maxiter 500 "
        "-tol 1e-12 -print none %s",
        method,
        strcmp(method,"fgcrodr")==0
            ? "-fgcrodr_k 0"
            : "");

    err = lis_solver_set_option(option,solver);
    if( err )
    {
        lis_solver_destroy(solver);
        return err;
    }

    lis_vector_set_all((LIS_SCALAR)0.0,x);

    err = lis_solve(A,b,x,solver);

    result->api = err;
    result->status = solver->retcode;
    result->iter = solver->iter;
    result->reported = solver->resid;
    result->true_residual =
        true_relative_residual(A,b,x,r);

    lis_solver_destroy(solver);

    return LIS_SUCCESS;
}



static LIS_INT
run_fgcrodr_core(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_VECTOR r,
    LIS_INT recycle_dim,
    RESULT *result)
{
    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT qerr;

    char option[512];


    memset(
        result,
        0,
        sizeof(*result));


    result->api = -999;
    result->status = -999;
    result->iter = -999;
    result->reported = -1.0;
    result->true_residual = -1.0;


    err = lis_solver_create(
        &solver);

    if( err )
        return err;


    snprintf(
        option,
        sizeof(option),
        "-i fgcrodr "
        "-p jacobi "
        "-restart 10 "
        "-maxiter 500 "
        "-tol 1e-12 "
        "-print none "
        "-fgcrodr_k %d",
        (int)recycle_dim);


    err = lis_solver_set_option(
        option,
        solver);

    if( err )
    {
        lis_solver_destroy(solver);
        return err;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    err = lis_solve(
        A,
        b,
        x,
        solver);


    result->api = err;
    result->status = solver->retcode;
    result->iter = solver->iter;
    result->reported = solver->resid;

    result->true_residual =
        true_relative_residual(
            A,
            b,
            x,
            r);


    qerr =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &result->recycle_dim,
            &result->recycle_action,
            &result->recycle_quality,
            &result->recycle_quality_min);

    if( qerr )
    {
        lis_solver_destroy(solver);
        return qerr;
    }


    qerr =
        lis_fgcrodr_validate_recycle(
            solver,
            &result->recycle_orth_error,
            &result->recycle_relation_error);

    if( qerr )
    {
        lis_solver_destroy(solver);
        return qerr;
    }


    lis_solver_destroy(solver);

    return LIS_SUCCESS;
}


static int
check_text_registration(void)
{
    LIS_SOLVER solver = NULL;
    LIS_INT err;
    int failed = 0;

    if( LIS_SOLVER_FGCRODR!=27 )
    {
        fprintf(
            stderr,
            "LIS_SOLVER_FGCRODR=%d, expected 27\n",
            (int)LIS_SOLVER_FGCRODR);

        return 1;
    }

    if( LIS_SOLVER_LEN!=27 )
    {
        fprintf(
            stderr,
            "LIS_SOLVER_LEN=%d, expected 27\n",
            (int)LIS_SOLVER_LEN);

        return 1;
    }

    err = lis_solver_create(&solver);
    if( err )
        return 1;

    err = lis_solver_set_option(
        "-i fgcrodr -print none",
        solver);

    if( err )
    {
        fprintf(
            stderr,
            "textual FGCRO-DR registration failed: %d\n",
            (int)err);

        failed = 1;
    }
    else if(
        solver->options[LIS_OPTIONS_SOLVER]
        != LIS_SOLVER_FGCRODR )
    {
        fprintf(
            stderr,
            "fgcrodr selected solver %d, expected %d\n",
            (int)solver->options[LIS_OPTIONS_SOLVER],
            (int)LIS_SOLVER_FGCRODR);

        failed = 1;
    }

    lis_solver_destroy(solver);

    return failed;
}



static int
check_state_configuration(void)
{
    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT k;
    LIS_INT update;
    LIS_INT dim;
    LIS_INT action;

    LIS_REAL keep;
    LIS_REAL reject;
    LIS_REAL quality;
    LIS_REAL qmin;

    int failed = 0;


    err = lis_solver_create(&solver);
    if( err )
        return 1;


    /* -------------------------------------------------------- */
    /* Defaults                                                  */
    /* -------------------------------------------------------- */

    err = lis_solver_get_fgcrodr_recycle_config(
        solver,
        &k,
        &keep,
        &reject,
        &update);

    if( err )
    {
        fprintf(
            stderr,
            "default FGCRO-DR config query failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    printf(
        "FGCRO-DR defaults k=%d keep=%.2f reject=%.2f update=%d\n",
        (int)k,
        (double)keep,
        (double)reject,
        (int)update);

    if( k!=8 ||
        fabs((double)keep-0.90)>1.0e-14 ||
        fabs((double)reject-0.10)>1.0e-14 ||
        update!=LIS_TRUE )
    {
        fprintf(
            stderr,
            "unexpected FGCRO-DR defaults\n");

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Textual configuration                                     */
    /* -------------------------------------------------------- */

    err = lis_solver_set_option(
        "-i fgcrodr "
        "-fgcrodr_k 5 "
        "-fgcrodr_keep 0.80 "
        "-fgcrodr_reject 0.20 "
        "-fgcrodr_update false",
        solver);

    if( err )
    {
        fprintf(
            stderr,
            "valid FGCRO-DR textual options failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_fgcrodr_check_params(solver);

    if( err )
    {
        fprintf(
            stderr,
            "valid FGCRO-DR configuration rejected: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_solver_get_fgcrodr_recycle_config(
        solver,
        &k,
        &keep,
        &reject,
        &update);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    if( k!=5 ||
        fabs((double)keep-0.80)>1.0e-14 ||
        fabs((double)reject-0.20)>1.0e-14 ||
        update!=LIS_FALSE )
    {
        fprintf(
            stderr,
            "textual FGCRO-DR configuration mismatch\n");

        failed = 1;
        goto cleanup;
    }


    /* No recycle vectors exist before the numerical core builds them. */

    err = lis_solver_get_fgcrodr_recycle(
        solver,
        &dim,
        &action,
        &quality,
        &qmin);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    if( dim!=0 ||
        action!=LIS_FGCRODR_RECYCLE_NONE )
    {
        fprintf(
            stderr,
            "unexpected recycle state before first recycle build\n");

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Sequence reset clears vectors, not user configuration.    */
    /* -------------------------------------------------------- */

    err = lis_solver_reset_sequence(solver);

    if( err )
    {
        fprintf(
            stderr,
            "sequence reset failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_solver_get_fgcrodr_recycle_config(
        solver,
        &k,
        &keep,
        &reject,
        &update);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    if( k!=5 ||
        fabs((double)keep-0.80)>1.0e-14 ||
        fabs((double)reject-0.20)>1.0e-14 ||
        update!=LIS_FALSE )
    {
        fprintf(
            stderr,
            "sequence reset destroyed FGCRO-DR configuration\n");

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Atomic public configuration API                           */
    /* -------------------------------------------------------- */

    err = lis_solver_set_fgcrodr_recycle(
        solver,
        6,
        (LIS_REAL)0.75,
        (LIS_REAL)0.25);

    if( err )
    {
        fprintf(
            stderr,
            "public recycle configuration API failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_fgcrodr_recycle_update(
        solver,
        LIS_TRUE);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_get_fgcrodr_recycle_config(
        solver,
        &k,
        &keep,
        &reject,
        &update);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    if( k!=6 ||
        fabs((double)keep-0.75)>1.0e-14 ||
        fabs((double)reject-0.25)>1.0e-14 ||
        update!=LIS_TRUE )
    {
        fprintf(
            stderr,
            "public recycle configuration API mismatch\n");

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Validation                                                */
    /* -------------------------------------------------------- */

    err = lis_solver_set_option(
        "-fgcrodr_keep 1.25",
        solver);

    if( err==LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "invalid -fgcrodr_keep accepted\n");

        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_option(
        "-fgcrodr_reject nan",
        solver);

    if( err==LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "NaN -fgcrodr_reject accepted\n");

        failed = 1;
        goto cleanup;
    }


    /*
     * Individual threshold options are order-independent.
     * Cross-threshold validation is performed by check_params().
     */
    err = lis_solver_set_option(
        "-fgcrodr_keep 0.10 "
        "-fgcrodr_reject 0.90",
        solver);

    if( err )
    {
        fprintf(
            stderr,
            "order-independent threshold assignment failed\n");

        failed = 1;
        goto cleanup;
    }

    err = lis_fgcrodr_check_params(solver);

    if( err==LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "reject > keep was not rejected\n");

        failed = 1;
        goto cleanup;
    }


    /* Restore valid configuration. */

    err = lis_solver_set_option(
        "-fgcrodr_k 0 "
        "-fgcrodr_keep 0.90 "
        "-fgcrodr_reject 0.10 "
        "-fgcrodr_update true",
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_fgcrodr_check_params(solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    printf(
        "FGCRO-DR state/config lifecycle: PASS\n");


cleanup:

    if( solver )
        lis_solver_destroy(solver);

    return failed;
}



static void
set_second_exact_solution(
    LIS_VECTOR u)
{
    LIS_INT i;


    for(i=0;i<u->n;i++)
    {
        LIS_REAL g =
            (LIS_REAL)(u->is+i+1);


        /*
         * Related, but deliberately non-identical, second system.
         *
         * The matrix remains unchanged.  Only the exact solution and
         * therefore the right-hand side change.
         */
        u->value[i] =
            (LIS_SCALAR)(
                1.0
                + 0.15*sin(0.071*g + 0.08)
                + 0.04*cos(0.037*g)
                + 0.01*sin(0.19*g));
    }
}


static int
check_sequence_recycle(
    LIS_MATRIX A)
{
    LIS_SOLVER solver = NULL;

    LIS_VECTOR u1 = NULL;
    LIS_VECTOR u2 = NULL;

    LIS_VECTOR b1 = NULL;
    LIS_VECTOR b2 = NULL;

    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;
    LIS_VECTOR diff = NULL;

    LIS_INT err;

    LIS_INT api1;
    LIS_INT api2;
    LIS_INT api3;

    LIS_INT status1;
    LIS_INT status2;
    LIS_INT status3;

    LIS_INT iter1;
    LIS_INT iter2;
    LIS_INT iter3;

    LIS_INT dim1;
    LIS_INT dim2;
    LIS_INT dim3;

    LIS_INT action1;
    LIS_INT action2;
    LIS_INT action3;

    LIS_REAL quality1;
    LIS_REAL quality2;

    LIS_REAL qmin1;
    LIS_REAL qmin2;

    LIS_REAL orth1;
    LIS_REAL orth2;
    LIS_REAL orth3;

    LIS_REAL relation1;
    LIS_REAL relation2;
    LIS_REAL relation3;

    LIS_REAL true1;
    LIS_REAL true2;
    LIS_REAL true3;

    LIS_UNSIGNED_INT builds1;
    LIS_UNSIGNED_INT builds2;
    LIS_UNSIGNED_INT builds3;

    LIS_UNSIGNED_INT updates1;
    LIS_UNSIGNED_INT updates2;
    LIS_UNSIGNED_INT updates3;

    LIS_UNSIGNED_INT uses1;
    LIS_UNSIGNED_INT uses2;
    LIS_UNSIGNED_INT uses3;

    LIS_UNSIGNED_INT rejects1;
    LIS_UNSIGNED_INT rejects2;
    LIS_UNSIGNED_INT rejects3;

    LIS_REAL diff_norm;
    LIS_REAL exact_norm;
    LIS_REAL solution_error;

    int failed = 0;


    err =
        lis_vector_duplicate(
            A,
            &u1);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &u2);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b1);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b2);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &x);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &r);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &diff);

    if( err )
        goto api_fail;


    set_exact_solution(
        u1);

    set_second_exact_solution(
        u2);


    err =
        lis_matvec(
            A,
            u1,
            b1);

    if( err )
        goto api_fail;


    err =
        lis_matvec(
            A,
            u2,
            b2);

    if( err )
        goto api_fail;


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p jacobi "
            "-restart 10 "
            "-maxiter 500 "
            "-tol 1e-12 "
            "-print none "
            "-fgcrodr_k 5",
            solver);

    if( err )
        goto api_fail;


    /* ======================================================== */
    /* Solve 1: no recycle on entry -> BUILD                    */
    /* ======================================================== */

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    api1 =
        lis_solve_sequence(
            A,
            b1,
            x,
            solver);

    status1 =
        solver->retcode;

    iter1 =
        solver->iter;

    true1 =
        true_relative_residual(
            A,
            b1,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim1,
            &action1,
            &quality1,
            &qmin1);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_validate_recycle(
            solver,
            &orth1,
            &relation1);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds1,
            &updates1,
            &uses1,
            &rejects1);

    if( err )
        goto api_fail;


    /* ======================================================== */
    /* Solve 2: same A generation -> active KEEP                */
    /*                                                          */
    /* Start again from zero so a nontrivial Krylov solve is     */
    /* required and reuse cannot be hidden by an already exact x.*/
    /* ======================================================== */

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    api2 =
        lis_solve_sequence(
            A,
            b2,
            x,
            solver);

    status2 =
        solver->retcode;

    iter2 =
        solver->iter;

    true2 =
        true_relative_residual(
            A,
            b2,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim2,
            &action2,
            &quality2,
            &qmin2);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_validate_recycle(
            solver,
            &orth2,
            &relation2);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds2,
            &updates2,
            &uses2,
            &rejects2);

    if( err )
        goto api_fail;


    /*
     * Independent solution check against the known exact second
     * solution, not only against the solver's reported residual.
     */
    lis_vector_copy(
        x,
        diff);

    lis_vector_axpy(
        (LIS_SCALAR)-1.0,
        u2,
        diff);

    lis_vector_nrm2(
        diff,
        &diff_norm);

    lis_vector_nrm2(
        u2,
        &exact_norm);

    solution_error =
        exact_norm>0.0
        ? diff_norm/exact_norm
        : diff_norm;


    /* ======================================================== */
    /* Solve 3: reuse remains active, harmonic UPDATE disabled  */
    /* ======================================================== */

    err =
        lis_solver_set_fgcrodr_recycle_update(
            solver,
            LIS_FALSE);

    if( err )
        goto api_fail;


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    api3 =
        lis_solve_sequence(
            A,
            b1,
            x,
            solver);

    status3 =
        solver->retcode;

    iter3 =
        solver->iter;

    true3 =
        true_relative_residual(
            A,
            b1,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim3,
            &action3,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_validate_recycle(
            solver,
            &orth3,
            &relation3);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds3,
            &updates3,
            &uses3,
            &rejects3);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR sequence solve1 "
        "api=%d status=%d iter=%d "
        "true=%.17e dim=%d action=%d "
        "builds=%llu updates=%llu uses=%llu "
        "orth=%.3e AZ-C=%.3e\n",
        (int)api1,
        (int)status1,
        (int)iter1,
        (double)true1,
        (int)dim1,
        (int)action1,
        (unsigned long long)builds1,
        (unsigned long long)updates1,
        (unsigned long long)uses1,
        (double)orth1,
        (double)relation1);


    printf(
        "FGCRO-DR sequence solve2 "
        "api=%d status=%d iter=%d "
        "true=%.17e dim=%d action=%d "
        "quality=%.3f solution_error=%.3e "
        "builds=%llu updates=%llu uses=%llu "
        "orth=%.3e AZ-C=%.3e\n",
        (int)api2,
        (int)status2,
        (int)iter2,
        (double)true2,
        (int)dim2,
        (int)action2,
        (double)quality2,
        (double)solution_error,
        (unsigned long long)builds2,
        (unsigned long long)updates2,
        (unsigned long long)uses2,
        (double)orth2,
        (double)relation2);


    printf(
        "FGCRO-DR sequence solve3 update=off "
        "api=%d status=%d iter=%d "
        "true=%.17e dim=%d action=%d "
        "builds=%llu updates=%llu uses=%llu "
        "orth=%.3e AZ-C=%.3e\n",
        (int)api3,
        (int)status3,
        (int)iter3,
        (double)true3,
        (int)dim3,
        (int)action3,
        (unsigned long long)builds3,
        (unsigned long long)updates3,
        (unsigned long long)uses3,
        (double)orth3,
        (double)relation3);


    /* ======================================================== */
    /* Stage 2C contract                                       */
    /* ======================================================== */

    if(
        api1!=LIS_SUCCESS ||
        status1!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "FGCRO-DR sequence solve1 failed: "
            "api=%d status=%d\n",
            (int)api1,
            (int)status1);

        failed = 1;
    }


    if(
        dim1!=5 ||
        action1!=LIS_FGCRODR_RECYCLE_BUILD )
    {
        fprintf(
            stderr,
            "FGCRO-DR sequence BUILD mismatch: "
            "dim=%d action=%d\n",
            (int)dim1,
            (int)action1);

        failed = 1;
    }


    if(
        builds1!=1 ||
        updates1!=0 ||
        uses1!=0 )
    {
        fprintf(
            stderr,
            "FGCRO-DR initial harmonic statistics mismatch: "
            "builds=%llu updates=%llu uses=%llu\n",
            (unsigned long long)builds1,
            (unsigned long long)updates1,
            (unsigned long long)uses1);

        failed = 1;
    }


    if( true1>1.0e-10 )
    {
        fprintf(
            stderr,
            "FGCRO-DR sequence solve1 true residual "
            "too large: %.17e\n",
            (double)true1);

        failed = 1;
    }


    if(
        api2!=LIS_SUCCESS ||
        status2!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "FGCRO-DR sequence solve2 failed: "
            "api=%d status=%d\n",
            (int)api2,
            (int)status2);

        failed = 1;
    }


    if( dim2!=5 )
    {
        fprintf(
            stderr,
            "FGCRO-DR sequence KEEP dimension mismatch: %d\n",
            (int)dim2);

        failed = 1;
    }


    if(
        action2
        != LIS_FGCRODR_RECYCLE_KEEP )
    {
        fprintf(
            stderr,
            "FGCRO-DR sequence expected KEEP, got action=%d\n",
            (int)action2);

        failed = 1;
    }


    if(
        builds2!=builds1 ||
        updates2<=updates1 ||
        uses2!=uses1+1 )
    {
        fprintf(
            stderr,
            "FGCRO-DR harmonic UPDATE statistics mismatch: "
            "before=%llu after=%llu uses=%llu->%llu\n",
            (unsigned long long)updates1,
            (unsigned long long)updates2,
            (unsigned long long)uses1,
            (unsigned long long)uses2);

        failed = 1;
    }


    if( true2>1.0e-10 )
    {
        fprintf(
            stderr,
            "FGCRO-DR sequence solve2 true residual "
            "too large: %.17e\n",
            (double)true2);

        failed = 1;
    }


    if( solution_error>1.0e-9 )
    {
        fprintf(
            stderr,
            "FGCRO-DR sequence solution error "
            "too large: %.17e\n",
            (double)solution_error);

        failed = 1;
    }


    if(
        orth2>1.0e-12 ||
        relation2>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR sequence recycle invariant failure: "
            "orth=%.17e AZ-C=%.17e\n",
            (double)orth2,
            (double)relation2);

        failed = 1;
    }


    if(
        api3!=LIS_SUCCESS ||
        status3!=LIS_SUCCESS ||
        dim3!=5 ||
        action3!=LIS_FGCRODR_RECYCLE_KEEP )
    {
        fprintf(
            stderr,
            "FGCRO-DR update=off solve failure: "
            "api=%d status=%d dim=%d action=%d\n",
            (int)api3,
            (int)status3,
            (int)dim3,
            (int)action3);

        failed = 1;
    }


    if(
        updates3!=updates2 ||
        builds3!=builds2 ||
        uses3!=uses2+1 )
    {
        fprintf(
            stderr,
            "FGCRO-DR update=off statistics mismatch: "
            "builds=%llu->%llu updates=%llu->%llu "
            "uses=%llu->%llu\n",
            (unsigned long long)builds2,
            (unsigned long long)builds3,
            (unsigned long long)updates2,
            (unsigned long long)updates3,
            (unsigned long long)uses2,
            (unsigned long long)uses3);

        failed = 1;
    }


    if(
        true3>1.0e-10 ||
        orth3>1.0e-12 ||
        relation3>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR update=off invariant failure: "
            "true=%.17e orth=%.17e AZ-C=%.17e\n",
            (double)true3,
            (double)orth3,
            (double)relation3);

        failed = 1;
    }


    if( failed )
    {
        fprintf(
            stderr,
            "FGCRO-DR SEQUENCE REUSE FAILED\n");
    }
    else
    {
        printf(
            "FGCRO-DR SEQUENCE REUSE PASSED\n");
    }


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( diff )
        lis_vector_destroy(
            diff);

    if( r )
        lis_vector_destroy(
            r);

    if( x )
        lis_vector_destroy(
            x);

    if( b2 )
        lis_vector_destroy(
            b2);

    if( b1 )
        lis_vector_destroy(
            b1);

    if( u2 )
        lis_vector_destroy(
            u2);

    if( u1 )
        lis_vector_destroy(
            u1);


    return failed;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR sequence API failure: %d\n",
        (int)err);

    failed = 1;

    goto cleanup;
}


int
main(int argc, char **argv)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR u = NULL;
    LIS_VECTOR b = NULL;

    LIS_VECTOR xf = NULL;
    LIS_VECTOR xg = NULL;
    LIS_VECTOR xk = NULL;
    LIS_VECTOR work = NULL;
    LIS_VECTOR diff = NULL;

    RESULT fgmres;
    RESULT fgcrodr;
    RESULT fgcrodr_core;

    LIS_REAL diff_norm = 0.0;
    LIS_REAL x_norm = 0.0;
    LIS_REAL relative_difference;

    LIS_REAL reported_difference;
    LIS_REAL true_difference;

    LIS_INT err;
    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if( err )
        return (int)err;

    if( check_text_registration() )
    {
        failed = 1;
        goto cleanup;
    }

    if( check_state_configuration() )
    {
        failed = 1;
        goto cleanup;
    }

    err = build_matrix(&A);
    if( err )
    {
        fprintf(
            stderr,
            "build_matrix failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    if( (err=lis_vector_duplicate(A,&u)) ) goto api_fail;
    if( (err=lis_vector_duplicate(A,&b)) ) goto api_fail;
    if( (err=lis_vector_duplicate(A,&xf)) ) goto api_fail;
    if( (err=lis_vector_duplicate(A,&xg)) ) goto api_fail;
    if( (err=lis_vector_duplicate(A,&xk)) ) goto api_fail;
    if( (err=lis_vector_duplicate(A,&work)) ) goto api_fail;
    if( (err=lis_vector_duplicate(A,&diff)) ) goto api_fail;

    if( check_sequence_recycle(A) )
    {
        failed = 1;
    }


    set_exact_solution(u);

    err = lis_matvec(A,u,b);
    if( err ) goto api_fail;

    err = run_solver(
        "fgmres",
        A,b,xf,work,
        &fgmres);

    if( err ) goto api_fail;

    err = run_solver(
        "fgcrodr",
        A,b,xg,work,
        &fgcrodr);

    if( err ) goto api_fail;

    err = run_fgcrodr_core(
        A,
        b,
        xk,
        work,
        5,
        &fgcrodr_core);

    if( err ) goto api_fail;


    printf(
        "FGMRES   api=%d status=%d iter=%d "
        "reported=%.17e true=%.17e\n",
        (int)fgmres.api,
        (int)fgmres.status,
        (int)fgmres.iter,
        (double)fgmres.reported,
        (double)fgmres.true_residual);

    printf(
        "FGCRO-DR api=%d status=%d iter=%d "
        "reported=%.17e true=%.17e\n",
        (int)fgcrodr.api,
        (int)fgcrodr.status,
        (int)fgcrodr.iter,
        (double)fgcrodr.reported,
        (double)fgcrodr.true_residual);


    printf(
        "FGCRO-DR core(k=5) api=%d status=%d iter=%d "
        "reported=%.17e true=%.17e "
        "recycle_dim=%d action=%d quality=%.3f "
        "orth=%.3e AZ-C=%.3e\n",
        (int)fgcrodr_core.api,
        (int)fgcrodr_core.status,
        (int)fgcrodr_core.iter,
        (double)fgcrodr_core.reported,
        (double)fgcrodr_core.true_residual,
        (int)fgcrodr_core.recycle_dim,
        (int)fgcrodr_core.recycle_action,
        (double)fgcrodr_core.recycle_quality,
        (double)fgcrodr_core.recycle_orth_error,
        (double)fgcrodr_core.recycle_relation_error);

    if( fgmres.api!=fgcrodr.api )
    {
        fprintf(stderr,"API return mismatch\n");
        failed = 1;
    }

    if( fgmres.status!=fgcrodr.status )
    {
        fprintf(stderr,"solver status mismatch\n");
        failed = 1;
    }

    if( fgmres.iter!=fgcrodr.iter )
    {
        fprintf(
            stderr,
            "iteration mismatch: %d vs %d\n",
            (int)fgmres.iter,
            (int)fgcrodr.iter);

        failed = 1;
    }

    if( fgmres.status!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "FGMRES reference did not converge: %d\n",
            (int)fgmres.status);

        failed = 1;
    }

    if( fgcrodr.status!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "FGCRO-DR reference path did not converge: %d\n",
            (int)fgcrodr.status);

        failed = 1;
    }

    reported_difference =
        fabs(
            (double)fgmres.reported
            - (double)fgcrodr.reported);

    true_difference =
        fabs(
            (double)fgmres.true_residual
            - (double)fgcrodr.true_residual);

    if( reported_difference>1.0e-14 )
    {
        fprintf(
            stderr,
            "reported residual mismatch: %.17e\n",
            (double)reported_difference);

        failed = 1;
    }

    if( true_difference>1.0e-14 )
    {
        fprintf(
            stderr,
            "true residual mismatch: %.17e\n",
            (double)true_difference);

        failed = 1;
    }

    if( fgmres.true_residual>1.0e-10 ||
        fgcrodr.true_residual>1.0e-10 )
    {
        fprintf(
            stderr,
            "true residual too large\n");

        failed = 1;
    }



    /*
     * Stage 2B contract:
     *
     * The first positive-k solve must construct a persistent
     * recycle space from the actual flexible Z directions.
     */

    if( fgcrodr_core.recycle_dim!=5 )
    {
        fprintf(
            stderr,
            "FGCRO-DR recycle BUILD dimension mismatch: %d\n",
            (int)fgcrodr_core.recycle_dim);

        failed = 1;
    }

    if(
        fgcrodr_core.recycle_action
        != LIS_FGCRODR_RECYCLE_BUILD )
    {
        fprintf(
            stderr,
            "FGCRO-DR recycle BUILD action mismatch: %d\n",
            (int)fgcrodr_core.recycle_action);

        failed = 1;
    }

    if(
        fabs(
            (double)fgcrodr_core.recycle_quality
            - 1.0)
        > 1.0e-14 )
    {
        fprintf(
            stderr,
            "FGCRO-DR initial recycle quality mismatch\n");

        failed = 1;
    }



    /*
     * Stage 2B recycle invariants.
     */

    if(
        fgcrodr_core.recycle_orth_error
        > 1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR C orthogonality error too large: %.17e\n",
            (double)fgcrodr_core.recycle_orth_error);

        failed = 1;
    }


    if(
        fgcrodr_core.recycle_relation_error
        > 1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR A*Z-C error too large: %.17e\n",
            (double)fgcrodr_core.recycle_relation_error);

        failed = 1;
    }


    /*
     * Stage 2A contract.
     *
     * With no active recycle projection, positive-k FGCRO-DR uses
     * its own flexible Arnoldi core and must remain numerically
     * equivalent to FGMRES for the same fixed preconditioner.
     */

    if( fgcrodr_core.api!=fgmres.api )
    {
        fprintf(
            stderr,
            "FGCRO-DR core API return mismatch\n");

        failed = 1;
    }

    if( fgcrodr_core.status!=fgmres.status )
    {
        fprintf(
            stderr,
            "FGCRO-DR core status mismatch\n");

        failed = 1;
    }

    if( fgcrodr_core.iter!=fgmres.iter )
    {
        fprintf(
            stderr,
            "FGCRO-DR core iteration mismatch: %d vs %d\n",
            (int)fgcrodr_core.iter,
            (int)fgmres.iter);

        failed = 1;
    }

    if(
        fabs(
            (double)fgcrodr_core.reported
            - (double)fgmres.reported)
        > 1.0e-14 )
    {
        fprintf(
            stderr,
            "FGCRO-DR core reported residual mismatch\n");

        failed = 1;
    }

    if(
        fabs(
            (double)fgcrodr_core.true_residual
            - (double)fgmres.true_residual)
        > 1.0e-14 )
    {
        fprintf(
            stderr,
            "FGCRO-DR core true residual mismatch\n");

        failed = 1;
    }


    /*
     * diff = xf - xk
     */
    lis_vector_copy(
        xf,
        diff);

    lis_vector_axpy(
        (LIS_SCALAR)-1.0,
        xk,
        diff);

    lis_vector_nrm2(
        diff,
        &diff_norm);

    lis_vector_nrm2(
        xf,
        &x_norm);

    relative_difference =
        (x_norm>0.0)
        ? diff_norm/x_norm
        : diff_norm;

    printf(
        "FGCRO-DR core solution relative difference=%.17e\n",
        (double)relative_difference);

    if( relative_difference>1.0e-13 )
    {
        fprintf(
            stderr,
            "FGCRO-DR core solution mismatch\n");

        failed = 1;
    }


    /*
     * diff = xf - xg
     */
    lis_vector_copy(xf,diff);
    lis_vector_axpy(
        (LIS_SCALAR)-1.0,
        xg,
        diff);

    lis_vector_nrm2(diff,&diff_norm);
    lis_vector_nrm2(xf,&x_norm);

    relative_difference =
        (x_norm>0.0)
        ? diff_norm/x_norm
        : diff_norm;

    printf(
        "solution relative difference=%.17e\n",
        (double)relative_difference);

    if( relative_difference>1.0e-13 )
    {
        fprintf(
            stderr,
            "solution mismatch\n");

        failed = 1;
    }

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "LIS API failed: %d\n",
        (int)err);

    failed = 1;


cleanup:

    if( diff ) lis_vector_destroy(diff);
    if( work ) lis_vector_destroy(work);
    if( xk ) lis_vector_destroy(xk);
    if( xg ) lis_vector_destroy(xg);
    if( xf ) lis_vector_destroy(xf);
    if( b ) lis_vector_destroy(b);
    if( u ) lis_vector_destroy(u);

    if( A ) lis_matrix_destroy(A);

    lis_finalize();

    if( failed )
    {
        fprintf(
            stderr,
            "LIS_FGCRODR_CORE FAILED\n");

        return 1;
    }

    printf(
        "LIS_FGCRODR_CORE PASSED\n");

    return 0;
}
