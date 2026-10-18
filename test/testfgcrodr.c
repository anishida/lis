#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include <math.h>
#include "lis.h"
#include "lis_solver.h"
#include "lis_precon.h"

#include <stdio.h>
#include <string.h>

#ifdef _COMPLEX
int
main(int argc, char **argv)
{
    (void)argc;
    (void)argv;

    printf("SKIP: FGCRO-DR recycling is not implemented for complex scalar builds\n");

    return 77;
}
#else

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



/*
 * Stage 2F-C deterministic recovery preconditioner.
 *
 * The first psolve intentionally distorts one component so a
 * one-iteration FGCRO-DR attempt cannot satisfy tol=0.2.
 *
 * All later psolve calls are identity.  Therefore Attempt 1 uses the
 * SAME persistent preconditioner object but receives an undistorted
 * flexible direction and converges.
 */
static LIS_INT fgcrodr_recovery_create_count = 0;
static LIS_INT fgcrodr_recovery_destroy_count = 0;
static LIS_INT fgcrodr_recovery_psolve_count = 0;


static LIS_INT
fgcrodr_recovery_create(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    (void)solver;
    (void)precon;

    fgcrodr_recovery_create_count++;

    return LIS_SUCCESS;
}


static LIS_INT
fgcrodr_recovery_psolve(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    LIS_INT err;

    (void)solver;


    err =
        lis_vector_copy(
            b,
            x);

    if( err )
        return err;


    fgcrodr_recovery_psolve_count++;


    /*
     * The Stage 2F-C recovery matrix is a serial 2x2 test problem.
     *
     * First flexible direction:
     *
     *     z = [v0, 0.01*v1]
     *
     * gives a one-step residual reduction of about 0.696, which is
     * deliberately worse than tol=0.2.
     *
     * Subsequent calls are identity.
     */
    if(
        fgcrodr_recovery_psolve_count==1 &&
        x->n>=2 )
    {
        x->value[1] *=
            (LIS_SCALAR)0.01;
    }


    return LIS_SUCCESS;
}


static LIS_INT
fgcrodr_recovery_destroy(
    LIS_PRECON precon)
{
    (void)precon;

    fgcrodr_recovery_destroy_count++;

    return LIS_SUCCESS;
}


static void
fgcrodr_recovery_reset_counts(void)
{
    fgcrodr_recovery_create_count = 0;
    fgcrodr_recovery_destroy_count = 0;
    fgcrodr_recovery_psolve_count = 0;
}



/*
 * Stage 2F-D deterministic preconditioner-rebuild test.
 *
 * The first persistent preconditioner instance intentionally distorts
 * every flexible direction:
 *
 *     z = [v0, 0.01*v1]
 *
 * On the 2x2 recovery system a single FGCRO-DR iteration then reduces
 * the current residual by about 0.696, which remains above tol=0.2.
 *
 * Therefore:
 *
 *     Attempt 0 -> MAXITER
 *     Attempt 1 -> MAXITER
 *
 * Sequence recovery then destroys the first persistent preconditioner.
 *
 * The second preconditioner instance is identity:
 *
 *     Attempt 2 -> SUCCESS
 *
 * This isolates the second recovery tier:
 *
 *     recycle reset
 *     persistent preconditioner rebuild
 *     fresh FGCRO recycle BUILD
 */
static LIS_INT fgcrodr_rebuild_create_count = 0;
static LIS_INT fgcrodr_rebuild_destroy_count = 0;
static LIS_INT fgcrodr_rebuild_psolve_count = 0;


static LIS_INT
fgcrodr_rebuild_create(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    (void)solver;
    (void)precon;

    fgcrodr_rebuild_create_count++;

    return LIS_SUCCESS;
}


static LIS_INT
fgcrodr_rebuild_psolve(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    LIS_INT err;

    (void)solver;


    err =
        lis_vector_copy(
            b,
            x);

    if( err )
        return err;


    fgcrodr_rebuild_psolve_count++;


    /*
     * Only the FIRST persistent preconditioner instance is poor.
     *
     * Generic recovery destroys it before Attempt 2.  Once create
     * count reaches two, this callback becomes identity.
     */
    if(
        fgcrodr_rebuild_create_count==1 &&
        x->n>=2 )
    {
        x->value[1] *=
            (LIS_SCALAR)0.01;
    }


    return LIS_SUCCESS;
}


static LIS_INT
fgcrodr_rebuild_destroy(
    LIS_PRECON precon)
{
    (void)precon;

    fgcrodr_rebuild_destroy_count++;

    return LIS_SUCCESS;
}


static void
fgcrodr_rebuild_reset_counts(void)
{
    fgcrodr_rebuild_create_count = 0;
    fgcrodr_rebuild_destroy_count = 0;
    fgcrodr_rebuild_psolve_count = 0;
}


static LIS_INT
build_recovery_matrix(
    LIS_MATRIX *Aout)
{
    LIS_MATRIX A = NULL;

    LIS_INT err;


    if( Aout==NULL )
        return LIS_ERR_ILL_ARG;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        return err;


    err =
        lis_matrix_set_size(
            A,
            0,
            2);

    if( err )
        goto fail;


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            0,
            0,
            (LIS_SCALAR)2.0,
            A);

    if( err )
        goto fail;


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            1,
            1,
            (LIS_SCALAR)3.0,
            A);

    if( err )
        goto fail;


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if( err )
        goto fail;


    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto fail;


    *Aout =
        A;


    return LIS_SUCCESS;


fail:

    lis_matrix_destroy(
        A);

    return err;
}



/*
 * Change CSR numerical values without changing ptr/index.
 *
 * mode 1:
 *     moderate nonuniform perturbation.  The operator image changes,
 *     but the old flexible recycle space remains usable.
 *
 * mode 2:
 *     strongly different, still diagonally dominant operator used to
 *     exercise the REJECT policy.
 */
static LIS_INT
perturb_matrix_values(
    LIS_MATRIX A,
    LIS_INT mode)
{
    LIS_INT i;
    LIS_INT p;
    LIS_INT col;

    LIS_SCALAR factor;


    if(
        A==NULL ||
        A->matrix_type!=LIS_MATRIX_CSR ||
        A->ptr==NULL ||
        A->index==NULL ||
        A->value==NULL )
    {
        return LIS_ERR_ILL_ARG;
    }


    for(i=0;i<A->n;i++)
    {
        for(
            p=A->ptr[i];
            p<A->ptr[i+1];
            p++ )
        {
            col =
                A->index[p];


            if( mode==1 )
            {
                /*
                 * Preserve a well-behaved matrix, but rotate the
                 * operator images through nonuniform row/entry
                 * perturbations.
                 */
                if( col==i )
                {
                    factor =
                        (i%2)
                        ? (LIS_SCALAR)1.35
                        : (LIS_SCALAR)0.75;
                }
                else if( col<i )
                {
                    factor =
                        (i%3)
                        ? (LIS_SCALAR)0.55
                        : (LIS_SCALAR)1.45;
                }
                else
                {
                    factor =
                        (i%4)
                        ? (LIS_SCALAR)1.60
                        : (LIS_SCALAR)0.50;
                }


                A->value[p] *=
                    factor;
            }
            else if( mode==2 )
            {
                /*
                 * Strongly different but diagonally dominant
                 * tridiagonal operator.
                 */
                if( col==i )
                {
                    A->value[p] =
                        (LIS_SCALAR)(
                            7.0
                            +0.4*(i%5));
                }
                else if( col<i )
                {
                    A->value[p] =
                        (i%2)
                        ? (LIS_SCALAR)-2.05
                        : (LIS_SCALAR)-0.05;
                }
                else
                {
                    A->value[p] =
                        (i%2)
                        ? (LIS_SCALAR)-0.03
                        : (LIS_SCALAR)-1.80;
                }
            }
            else
            {
                return LIS_ERR_ILL_ARG;
            }
        }
    }


    return
        lis_matrix_mark_values_changed(
            A);
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
    LIS_INT api4;

    LIS_INT status1;
    LIS_INT status2;
    LIS_INT status3;
    LIS_INT status4;

    LIS_INT iter1;
    LIS_INT iter2;
    LIS_INT iter3;
    LIS_INT iter4;

    LIS_INT dim1;
    LIS_INT dim2;
    LIS_INT dim3;
    LIS_INT dim4;

    LIS_INT action1;
    LIS_INT action2;
    LIS_INT action3;
    LIS_INT action4;

    LIS_REAL quality1;
    LIS_REAL quality2;
    LIS_REAL quality4;

    LIS_REAL qmin1;
    LIS_REAL qmin2;
    LIS_REAL qmin4;

    LIS_REAL orth1;
    LIS_REAL orth2;
    LIS_REAL orth3;
    LIS_REAL orth4;

    LIS_REAL relation1;
    LIS_REAL relation2;
    LIS_REAL relation3;
    LIS_REAL relation4;

    LIS_REAL true1;
    LIS_REAL true2;
    LIS_REAL true3;
    LIS_REAL true4;

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

    LIS_UNSIGNED_INT cross4;
    LIS_UNSIGNED_INT keeps4;
    LIS_UNSIGNED_INT transports4;
    LIS_UNSIGNED_INT policy_rejects4;

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


    /* ======================================================== */
    /* Solve 4: value-generation change, identical numeric A    */
    /*                                                          */
    /* This isolates generation-aware cross-system KEEP.        */
    /* ======================================================== */

    err =
        lis_solver_set_fgcrodr_recycle_update(
            solver,
            LIS_TRUE);

    if( err )
        goto api_fail;


    err =
        lis_matrix_mark_values_changed(
            A);

    if( err )
        goto api_fail;


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    api4 =
        lis_solve_sequence(
            A,
            b2,
            x,
            solver);

    status4 =
        solver->retcode;

    iter4 =
        solver->iter;

    true4 =
        true_relative_residual(
            A,
            b2,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim4,
            &action4,
            &quality4,
            &qmin4);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_validate_recycle(
            solver,
            &orth4,
            &relation4);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_policy_stats(
            solver,
            &cross4,
            &keeps4,
            &transports4,
            &policy_rejects4);

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


    printf(
        "FGCRO-DR sequence solve4 value-generation "
        "api=%d status=%d iter=%d "
        "true=%.17e dim=%d action=%d "
        "quality=%.6f qmin=%.6f "
        "cross=%llu keeps=%llu transports=%llu rejects=%llu "
        "orth=%.3e AZ-C=%.3e\n",
        (int)api4,
        (int)status4,
        (int)iter4,
        (double)true4,
        (int)dim4,
        (int)action4,
        (double)quality4,
        (double)qmin4,
        (unsigned long long)cross4,
        (unsigned long long)keeps4,
        (unsigned long long)transports4,
        (unsigned long long)policy_rejects4,
        (double)orth4,
        (double)relation4);


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


    if(
        api4!=LIS_SUCCESS ||
        status4!=LIS_SUCCESS ||
        dim4!=5 ||
        action4!=LIS_FGCRODR_RECYCLE_KEEP )
    {
        fprintf(
            stderr,
            "FGCRO-DR generation KEEP failure: "
            "api=%d status=%d dim=%d action=%d\n",
            (int)api4,
            (int)status4,
            (int)dim4,
            (int)action4);

        failed = 1;
    }


    if(
        quality4<0.999999 ||
        qmin4<0.999999 )
    {
        fprintf(
            stderr,
            "FGCRO-DR generation KEEP quality failure: "
            "quality=%.17e qmin=%.17e\n",
            (double)quality4,
            (double)qmin4);

        failed = 1;
    }


    if(
        cross4!=1 ||
        transports4!=0 ||
        policy_rejects4!=0 )
    {
        fprintf(
            stderr,
            "FGCRO-DR generation policy counters mismatch: "
            "cross=%llu transports=%llu rejects=%llu\n",
            (unsigned long long)cross4,
            (unsigned long long)transports4,
            (unsigned long long)policy_rejects4);

        failed = 1;
    }


    if(
        true4>1.0e-10 ||
        orth4>1.0e-12 ||
        relation4>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR generation KEEP invariant failure: "
            "true=%.17e orth=%.17e AZ-C=%.17e\n",
            (double)true4,
            (double)orth4,
            (double)relation4);

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




static LIS_INT
check_generation_transport_reject(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR u = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT failed;

    LIS_INT api0;
    LIS_INT api1;
    LIS_INT api2;

    LIS_INT status0;
    LIS_INT status1;
    LIS_INT status2;

    LIS_INT dim0;
    LIS_INT dim1;
    LIS_INT dim2;

    LIS_INT action0;
    LIS_INT action1;
    LIS_INT action2;

    LIS_REAL quality1;
    LIS_REAL qmin1;

    LIS_REAL true0;
    LIS_REAL true1;
    LIS_REAL true2;

    LIS_REAL orth1;
    LIS_REAL relation1;

    LIS_REAL orth2;
    LIS_REAL relation2;

    LIS_UNSIGNED_INT cross1;
    LIS_UNSIGNED_INT keep1;
    LIS_UNSIGNED_INT transport1;
    LIS_UNSIGNED_INT reject1;

    LIS_UNSIGNED_INT cross2;
    LIS_UNSIGNED_INT keep2;
    LIS_UNSIGNED_INT transport2;
    LIS_UNSIGNED_INT reject2;

    LIS_UNSIGNED_INT builds0;
    LIS_UNSIGNED_INT builds2;

    failed = 0;


    err =
        build_matrix(
            &A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &u);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

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


    /*
     * Use the existing deterministic second exact solution.
     */
    set_second_exact_solution(
        u);


    lis_matvec(
        A,
        u,
        b);


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


    /*
     * TRANSPORT test:
     *
     * reject=0 guarantees a nondegenerate changed space is accepted;
     * keep≈1 ensures a genuinely changed operator does not classify as
     * KEEP.
     */
    err =
        lis_solver_set_fgcrodr_recycle(
            solver,
            5,
            (LIS_REAL)0.999999,
            (LIS_REAL)0.0);

    if( err )
        goto api_fail;


    /* -------------------------------------------------------- */
    /* Solve 0: establish the initial harmonic recycle space    */
    /* -------------------------------------------------------- */

    lis_vector_set_all(
        0.0,
        x);


    api0 =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    status0 =
        solver->retcode;


    true0 =
        true_relative_residual(
            A,
            b,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim0,
            &action0,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds0,
            NULL,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    /* -------------------------------------------------------- */
    /* Solve 1: real value change -> TRANSPORT                  */
    /* -------------------------------------------------------- */

    err =
        perturb_matrix_values(
            A,
            1);

    if( err )
        goto api_fail;


    /*
     * Preserve a known exact solution for the changed operator.
     */
    lis_matvec(
        A,
        u,
        b);


    lis_vector_set_all(
        0.0,
        x);


    api1 =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    status1 =
        solver->retcode;


    true1 =
        true_relative_residual(
            A,
            b,
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
        lis_fgcrodr_get_policy_stats(
            solver,
            &cross1,
            &keep1,
            &transport1,
            &reject1);

    if( err )
        goto api_fail;


    /* -------------------------------------------------------- */
    /* Solve 2: strong value change -> REJECT + fresh BUILD     */
    /* -------------------------------------------------------- */

    err =
        lis_solver_set_fgcrodr_recycle(
            solver,
            5,
            (LIS_REAL)1.0,
            (LIS_REAL)0.999999);

    if( err )
        goto api_fail;


    err =
        perturb_matrix_values(
            A,
            2);

    if( err )
        goto api_fail;


    lis_matvec(
        A,
        u,
        b);


    lis_vector_set_all(
        0.0,
        x);


    api2 =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    status2 =
        solver->retcode;


    true2 =
        true_relative_residual(
            A,
            b,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim2,
            &action2,
            NULL,
            NULL);

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
        lis_fgcrodr_get_policy_stats(
            solver,
            &cross2,
            &keep2,
            &transport2,
            &reject2);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds2,
            NULL,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR policy base "
        "api=%d status=%d dim=%d action=%d "
        "true=%.3e builds=%llu\n",
        (int)api0,
        (int)status0,
        (int)dim0,
        (int)action0,
        (double)true0,
        (unsigned long long)builds0);


    printf(
        "FGCRO-DR policy transport "
        "api=%d status=%d dim=%d action=%d "
        "quality=%.6f qmin=%.6f "
        "cross=%llu keep=%llu transport=%llu reject=%llu "
        "true=%.3e orth=%.3e AZ-C=%.3e\n",
        (int)api1,
        (int)status1,
        (int)dim1,
        (int)action1,
        (double)quality1,
        (double)qmin1,
        (unsigned long long)cross1,
        (unsigned long long)keep1,
        (unsigned long long)transport1,
        (unsigned long long)reject1,
        (double)true1,
        (double)orth1,
        (double)relation1);


    printf(
        "FGCRO-DR policy reject "
        "api=%d status=%d dim=%d action=%d "
        "cross=%llu keep=%llu transport=%llu reject=%llu "
        "builds=%llu true=%.3e orth=%.3e AZ-C=%.3e\n",
        (int)api2,
        (int)status2,
        (int)dim2,
        (int)action2,
        (unsigned long long)cross2,
        (unsigned long long)keep2,
        (unsigned long long)transport2,
        (unsigned long long)reject2,
        (unsigned long long)builds2,
        (double)true2,
        (double)orth2,
        (double)relation2);


    if(
        api0!=LIS_SUCCESS ||
        status0!=LIS_SUCCESS ||
        dim0!=5 ||
        action0!=LIS_FGCRODR_RECYCLE_BUILD ||
        true0>1.0e-10 )
    {
        fprintf(
            stderr,
            "FGCRO-DR policy base failure\n");

        failed = 1;
    }


    if(
        api1!=LIS_SUCCESS ||
        status1!=LIS_SUCCESS ||
        dim1!=5 ||
        action1!=
            LIS_FGCRODR_RECYCLE_TRANSPORT )
    {
        fprintf(
            stderr,
            "FGCRO-DR TRANSPORT action failure: "
            "api=%d status=%d dim=%d action=%d "
            "quality=%.17e qmin=%.17e\n",
            (int)api1,
            (int)status1,
            (int)dim1,
            (int)action1,
            (double)quality1,
            (double)qmin1);

        failed = 1;
    }


    if(
        qmin1>=0.999999 ||
        cross1!=1 ||
        transport1!=1 ||
        reject1!=0 )
    {
        fprintf(
            stderr,
            "FGCRO-DR TRANSPORT policy failure: "
            "qmin=%.17e cross=%llu transport=%llu reject=%llu\n",
            (double)qmin1,
            (unsigned long long)cross1,
            (unsigned long long)transport1,
            (unsigned long long)reject1);

        failed = 1;
    }


    if(
        true1>1.0e-10 ||
        orth1>1.0e-12 ||
        relation1>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR TRANSPORT invariant failure: "
            "true=%.17e orth=%.17e AZ-C=%.17e\n",
            (double)true1,
            (double)orth1,
            (double)relation1);

        failed = 1;
    }


    /*
     * REJECT must remain externally visible even though the current
     * solve is allowed to build a fresh recycle space afterward.
     */
    if(
        api2!=LIS_SUCCESS ||
        status2!=LIS_SUCCESS ||
        action2!=
            LIS_FGCRODR_RECYCLE_REJECT )
    {
        fprintf(
            stderr,
            "FGCRO-DR REJECT action failure: "
            "api=%d status=%d dim=%d action=%d\n",
            (int)api2,
            (int)status2,
            (int)dim2,
            (int)action2);

        failed = 1;
    }


    if(
        cross2!=cross1 ||
        transport2!=transport1 ||
        reject2!=reject1+1 )
    {
        fprintf(
            stderr,
            "FGCRO-DR REJECT counters failure: "
            "cross=%llu->%llu transport=%llu->%llu "
            "reject=%llu->%llu\n",
            (unsigned long long)cross1,
            (unsigned long long)cross2,
            (unsigned long long)transport1,
            (unsigned long long)transport2,
            (unsigned long long)reject1,
            (unsigned long long)reject2);

        failed = 1;
    }


    /*
     * With harmonic updating enabled, rejection of the stale basis
     * should be followed by construction of a fresh basis.
     */
    if(
        dim2!=5 ||
        builds2<=builds0 ||
        true2>1.0e-10 ||
        orth2>1.0e-12 ||
        relation2>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR REJECT rebuild failure: "
            "dim=%d builds=%llu->%llu "
            "true=%.17e orth=%.17e AZ-C=%.17e\n",
            (int)dim2,
            (unsigned long long)builds0,
            (unsigned long long)builds2,
            (double)true2,
            (double)orth2,
            (double)relation2);

        failed = 1;
    }


    if( failed )
    {
        fprintf(
            stderr,
            "FGCRO-DR GENERATION POLICY FAILED\n");

        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    printf(
        "FGCRO-DR GENERATION POLICY PASSED\n");

    err =
        LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR generation policy API failure: %d\n",
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( r )
        lis_vector_destroy(
            r);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( u )
        lis_vector_destroy(
            u);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}




static LIS_INT
check_sequence_reset_lifecycle(void)
{
    LIS_MATRIX A = NULL;
    LIS_MATRIX stored_matrix;

    LIS_VECTOR u = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    LIS_SEQUENCE_DIAGNOSTICS diag;

    LIS_INT err;
    LIS_INT failed;

    LIS_INT valid;

    LIS_INT api1;
    LIS_INT api2;

    LIS_INT status1;
    LIS_INT status2;

    LIS_INT dim1;
    LIS_INT dim0;
    LIS_INT dim2;

    LIS_INT action1;
    LIS_INT action0;
    LIS_INT action2;

    LIS_REAL true1;
    LIS_REAL true2;

    LIS_REAL orth2;
    LIS_REAL relation2;

    LIS_UNSIGNED_INT sg;
    LIS_UNSIGNED_INT vg;

    LIS_UNSIGNED_INT builds1;
    LIS_UNSIGNED_INT builds0;
    LIS_UNSIGNED_INT builds2;

    LIS_UNSIGNED_INT updates1;
    LIS_UNSIGNED_INT updates0;

    LIS_UNSIGNED_INT uses1;
    LIS_UNSIGNED_INT uses0;

    LIS_UNSIGNED_INT rejects1;
    LIS_UNSIGNED_INT rejects0;


    failed = 0;


    err =
        build_matrix(
            &A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &u);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

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


    set_second_exact_solution(
        u);


    lis_matvec(
        A,
        u,
        b);


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


    /* -------------------------------------------------------- */
    /* Solve 1: establish persistent FGCRO recycle state        */
    /* -------------------------------------------------------- */

    lis_vector_set_all(
        0.0,
        x);


    api1 =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    status1 =
        solver->retcode;


    true1 =
        true_relative_residual(
            A,
            b,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim1,
            &action1,
            NULL,
            NULL);

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


    err =
        lis_fgcrodr_get_lifecycle_state(
            solver,
            &valid,
            &stored_matrix,
            &sg,
            &vg);

    if( err )
        goto api_fail;


    if(
        api1!=LIS_SUCCESS ||
        status1!=LIS_SUCCESS ||
        dim1!=5 ||
        action1!=LIS_FGCRODR_RECYCLE_BUILD ||
        !valid ||
        stored_matrix!=A ||
        true1>1.0e-10 )
    {
        fprintf(
            stderr,
            "FGCRO-DR RESET initial BUILD failure: "
            "api=%d status=%d dim=%d action=%d "
            "valid=%d matrix=%p/%p true=%.17e\n",
            (int)api1,
            (int)status1,
            (int)dim1,
            (int)action1,
            (int)valid,
            (void *)stored_matrix,
            (void *)A,
            (double)true1);

        failed = 1;
    }


    /* -------------------------------------------------------- */
    /* Explicit sequence reset                                  */
    /* -------------------------------------------------------- */

    err =
        lis_solver_reset_sequence(
            solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim0,
            &action0,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds0,
            &updates0,
            &uses0,
            &rejects0);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_lifecycle_state(
            solver,
            &valid,
            &stored_matrix,
            &sg,
            &vg);

    if( err )
        goto api_fail;


    err =
        lis_solver_get_sequence_diagnostics(
            solver,
            &diag);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR reset cleared "
        "dim=%d action=%d valid=%d matrix=%p "
        "sg=%llu vg=%llu "
        "builds=%llu updates=%llu uses=%llu rejects=%llu "
        "diag=%d/%d/%d/%d\n",
        (int)dim0,
        (int)action0,
        (int)valid,
        (void *)stored_matrix,
        (unsigned long long)sg,
        (unsigned long long)vg,
        (unsigned long long)builds0,
        (unsigned long long)updates0,
        (unsigned long long)uses0,
        (unsigned long long)rejects0,
        (int)diag.attempts,
        (int)diag.retries,
        (int)diag.recycle_resets,
        (int)diag.precon_rebuilds);


    if(
        dim0!=0 ||
        action0!=LIS_FGCRODR_RECYCLE_NONE ||
        valid ||
        stored_matrix!=NULL ||
        sg!=0 ||
        vg!=0 )
    {
        fprintf(
            stderr,
            "FGCRO-DR RESET metadata failure: "
            "dim=%d action=%d valid=%d matrix=%p "
            "sg=%llu vg=%llu\n",
            (int)dim0,
            (int)action0,
            (int)valid,
            (void *)stored_matrix,
            (unsigned long long)sg,
            (unsigned long long)vg);

        failed = 1;
    }


    /*
     * Clearing recycle vectors must not erase lifetime FGCRO
     * statistics/configuration.
     */
    if(
        builds0!=builds1 ||
        updates0!=updates1 ||
        uses0!=uses1 ||
        rejects0!=rejects1 )
    {
        fprintf(
            stderr,
            "FGCRO-DR RESET lifetime counter failure\n");

        failed = 1;
    }


    /*
     * reset_sequence() clears per-call diagnostics, while lifetime
     * sequence totals intentionally survive.
     */
    if(
        diag.attempts!=0 ||
        diag.retries!=0 ||
        diag.recycle_resets!=0 ||
        diag.precon_rebuilds!=0 )
    {
        fprintf(
            stderr,
            "FGCRO-DR RESET sequence diagnostic failure: "
            "%d/%d/%d/%d\n",
            (int)diag.attempts,
            (int)diag.retries,
            (int)diag.recycle_resets,
            (int)diag.precon_rebuilds);

        failed = 1;
    }


    /* -------------------------------------------------------- */
    /* Solve 2: reset must cause a completely fresh BUILD       */
    /* -------------------------------------------------------- */

    lis_vector_set_all(
        0.0,
        x);


    api2 =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    status2 =
        solver->retcode;


    true2 =
        true_relative_residual(
            A,
            b,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim2,
            &action2,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds2,
            NULL,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_validate_recycle(
            solver,
            &orth2,
            &relation2);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR reset rebuild "
        "api=%d status=%d dim=%d action=%d "
        "builds=%llu->%llu true=%.3e "
        "orth=%.3e AZ-C=%.3e\n",
        (int)api2,
        (int)status2,
        (int)dim2,
        (int)action2,
        (unsigned long long)builds1,
        (unsigned long long)builds2,
        (double)true2,
        (double)orth2,
        (double)relation2);


    if(
        api2!=LIS_SUCCESS ||
        status2!=LIS_SUCCESS ||
        dim2!=5 ||
        action2!=LIS_FGCRODR_RECYCLE_BUILD ||
        builds2!=builds1+1 ||
        true2>1.0e-10 ||
        orth2>1.0e-12 ||
        relation2>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR RESET rebuild failure\n");

        failed = 1;
    }


    if( failed )
    {
        fprintf(
            stderr,
            "FGCRO-DR SEQUENCE RESET LIFECYCLE FAILED\n");

        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    printf(
        "FGCRO-DR SEQUENCE RESET LIFECYCLE PASSED\n");

    err =
        LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR sequence reset API failure: %d\n",
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( r )
        lis_vector_destroy(
            r);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( u )
        lis_vector_destroy(
            u);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}




static LIS_INT
check_sequence_recovery_lifecycle(void)
{
    LIS_MATRIX A = NULL;
    LIS_MATRIX stored_matrix = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    LIS_SEQUENCE_DIAGNOSTICS diag;

    LIS_INT err;
    LIS_INT failed;

    LIS_INT api;
    LIS_INT status;

    LIS_INT dim;
    LIS_INT action;
    LIS_INT valid;

    LIS_REAL true_residual;
    LIS_REAL orth;
    LIS_REAL relation;

    LIS_UNSIGNED_INT sg;
    LIS_UNSIGNED_INT vg;

    LIS_UNSIGNED_INT builds;
    LIS_UNSIGNED_INT updates;
    LIS_UNSIGNED_INT uses;
    LIS_UNSIGNED_INT rejects;


    failed =
        0;


    err =
        build_recovery_matrix(
            &A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

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


    lis_vector_set_all(
        (LIS_SCALAR)2.0,
        b);


    fgcrodr_recovery_reset_counts();


    err =
        lis_precon_register_ex(
            "fgc2fc",
            fgcrodr_recovery_create,
            fgcrodr_recovery_psolve,
            fgcrodr_recovery_psolve,
            fgcrodr_recovery_destroy);

    if( err )
        goto api_fail;


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    /*
     * Use the deterministic Stage 2F-C flexible preconditioner.
     *
     * The persistent preconditioner object must survive Attempt 0 and
     * Attempt 1.  Only its first psolve call is intentionally
     * distorted.
     */
    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p fgc2fc "
            "-scale none "
            "-print none "
            "-restart 2 "
            "-maxiter 20 "
            "-tol 1e-12 "
            "-fgcrodr_k 1",
            solver);

    if( err )
        goto api_fail;


    /*
     * Attempt 0:
     *
     *     distorted flexible direction
     *     one-step normalized residual ~= 0.696 > 0.2
     *     -> MAXITER
     *
     * Recovery:
     *
     *     clear recycle only
     *     retain persistent preconditioner
     *
     * Attempt 1:
     *
     *     same preconditioner object
     *     identity flexible direction
     *     improved warm start
     *     one-step normalized residual << 0.2
     *     -> SUCCESS
     */
    solver->options[LIS_OPTIONS_MAXITER] =
        1;

    solver->params[
        LIS_PARAMS_RESID
        -LIS_OPTIONS_LEN] =
            (LIS_SCALAR)2.0e-1;


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    api =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    status =
        solver->retcode;


    true_residual =
        true_relative_residual(
            A,
            b,
            x,
            r);


    err =
        lis_solver_get_sequence_diagnostics(
            solver,
            &diag);

    if( err )
        goto api_fail;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds,
            &updates,
            &uses,
            &rejects);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_lifecycle_state(
            solver,
            &valid,
            &stored_matrix,
            &sg,
            &vg);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_validate_recycle(
            solver,
            &orth,
            &relation);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR recovery "
        "api=%d status=%d "
        "attempts=%d retries=%d resets=%d rebuilds=%d "
        "first=%d final=%d warm=%d/%d "
        "dim=%d action=%d "
        "builds=%llu updates=%llu uses=%llu rejects=%llu "
        "precon=%d/%d psolve=%d "
        "true=%.6e orth=%.3e AZ-C=%.3e\n",
        (int)api,
        (int)status,
        (int)diag.attempts,
        (int)diag.retries,
        (int)diag.recycle_resets,
        (int)diag.precon_rebuilds,
        (int)diag.first_retcode,
        (int)diag.final_retcode,
        (int)diag.initial_warm_start_kept,
        (int)diag.last_warm_start_kept,
        (int)dim,
        (int)action,
        (unsigned long long)builds,
        (unsigned long long)updates,
        (unsigned long long)uses,
        (unsigned long long)rejects,
        (int)fgcrodr_recovery_create_count,
        (int)fgcrodr_recovery_destroy_count,
        (int)fgcrodr_recovery_psolve_count,
        (double)true_residual,
        (double)orth,
        (double)relation);


    /*
     * Generic recovery contract.
     */
    if(
        api!=LIS_SUCCESS ||
        status!=LIS_SUCCESS ||
        diag.attempts!=2 ||
        diag.retries!=1 ||
        diag.recycle_resets!=1 ||
        diag.precon_rebuilds!=0 ||
        diag.first_retcode!=LIS_MAXITER ||
        diag.final_retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery diagnostics failure: "
            "api=%d status=%d "
            "attempts=%d retries=%d resets=%d rebuilds=%d "
            "first=%d final=%d\n",
            (int)api,
            (int)status,
            (int)diag.attempts,
            (int)diag.retries,
            (int)diag.recycle_resets,
            (int)diag.precon_rebuilds,
            (int)diag.first_retcode,
            (int)diag.final_retcode);

        failed =
            1;
    }


    /*
     * Attempt 0 starts from zero.
     * Attempt 1 must preserve the useful partial solution produced by
     * the MAXITER attempt.
     */
    if(
        diag.initial_warm_start_kept!=LIS_FALSE ||
        diag.last_warm_start_kept!=LIS_TRUE )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery warm-start failure: "
            "initial=%d last=%d\n",
            (int)diag.initial_warm_start_kept,
            (int)diag.last_warm_start_kept);

        failed =
            1;
    }


    /*
     * Attempt 1 must retain the same persistent preconditioner.
     * Rebuild is reserved for Attempt 2.
     */
    if(
        fgcrodr_recovery_create_count!=1 ||
        fgcrodr_recovery_destroy_count!=0 ||
        fgcrodr_recovery_psolve_count<2 )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery preconditioner lifecycle failure: "
            "create=%d destroy=%d psolve=%d\n",
            (int)fgcrodr_recovery_create_count,
            (int)fgcrodr_recovery_destroy_count,
            (int)fgcrodr_recovery_psolve_count);

        failed =
            1;
    }


    /*
     * Both attempts reach the harmonic recycle BUILD point:
     *
     *     attempt 0 -> BUILD #1 -> MAXITER
     *     recovery  -> clear recycle
     *     attempt 1 -> BUILD #2 -> SUCCESS
     *
     * If recovery accidentally reused the first basis, builds would
     * not have this clean 2/0/0 lifecycle.
     */
    if(
        builds!=2 ||
        updates!=0 ||
        uses!=0 ||
        rejects!=0 ||
        dim!=1 ||
        action!=LIS_FGCRODR_RECYCLE_BUILD )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery recycle lifecycle failure: "
            "builds=%llu updates=%llu uses=%llu rejects=%llu "
            "dim=%d action=%d\n",
            (unsigned long long)builds,
            (unsigned long long)updates,
            (unsigned long long)uses,
            (unsigned long long)rejects,
            (int)dim,
            (int)action);

        failed =
            1;
    }


    /*
     * The successful retry owns a fresh, generation-compatible
     * recycle relation.
     */
    if(
        !valid ||
        stored_matrix!=A ||
        true_residual>2.0e-1 ||
        orth>1.0e-12 ||
        relation>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery final-state failure: "
            "valid=%d matrix=%p/%p "
            "true=%.17e orth=%.17e AZ-C=%.17e "
            "sg=%llu vg=%llu\n",
            (int)valid,
            (void *)stored_matrix,
            (void *)A,
            (double)true_residual,
            (double)orth,
            (double)relation,
            (unsigned long long)sg,
            (unsigned long long)vg);

        failed =
            1;
    }


    /*
     * Fresh solver: lifetime sequence diagnostics must reflect exactly
     * this one call and its one recovery retry.
     */
    if(
        diag.total_solves!=1 ||
        diag.total_attempts!=2 ||
        diag.total_retries!=1 ||
        diag.total_recycle_resets!=1 ||
        diag.total_precon_rebuilds!=0 )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery lifetime diagnostics failure: "
            "%llu/%llu/%llu/%llu/%llu\n",
            (unsigned long long)diag.total_solves,
            (unsigned long long)diag.total_attempts,
            (unsigned long long)diag.total_retries,
            (unsigned long long)diag.total_recycle_resets,
            (unsigned long long)diag.total_precon_rebuilds);

        failed =
            1;
    }


    if( failed )
    {
        fprintf(
            stderr,
            "FGCRO-DR SEQUENCE RECOVERY LIFECYCLE FAILED\n");

        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    printf(
        "FGCRO-DR SEQUENCE RECOVERY LIFECYCLE PASSED\n");


    err =
        LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR sequence recovery API failure: %d\n",
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( r )
        lis_vector_destroy(
            r);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}




static LIS_INT
check_sequence_recovery_rebuild_lifecycle(void)
{
    LIS_MATRIX A = NULL;
    LIS_MATRIX stored_matrix = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    LIS_SEQUENCE_DIAGNOSTICS diag;

    LIS_INT err;
    LIS_INT failed;

    LIS_INT api;
    LIS_INT status;

    LIS_INT dim;
    LIS_INT action;
    LIS_INT valid;

    LIS_REAL true_residual;
    LIS_REAL orth;
    LIS_REAL relation;

    LIS_UNSIGNED_INT sg;
    LIS_UNSIGNED_INT vg;

    LIS_UNSIGNED_INT builds;
    LIS_UNSIGNED_INT updates;
    LIS_UNSIGNED_INT uses;
    LIS_UNSIGNED_INT rejects;


    failed =
        0;


    err =
        build_recovery_matrix(
            &A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

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


    lis_vector_set_all(
        (LIS_SCALAR)2.0,
        b);


    fgcrodr_rebuild_reset_counts();


    err =
        lis_precon_register_ex(
            "fgc2fd",
            fgcrodr_rebuild_create,
            fgcrodr_rebuild_psolve,
            fgcrodr_rebuild_psolve,
            fgcrodr_rebuild_destroy);

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
            "-p fgc2fd "
            "-scale none "
            "-print none "
            "-restart 2 "
            "-maxiter 20 "
            "-tol 1e-12 "
            "-fgcrodr_k 1",
            solver);

    if( err )
        goto api_fail;


    /*
     * One iteration per sequence attempt.
     *
     * Inherited first preconditioner instance:
     *
     *     Attempt 0 remains above the 0.3 solver threshold
     *     Attempt 1 remains above the 0.3 solver threshold
     *
     * This forces the second recovery tier.
     *
     * Rebuilt preconditioner:
     *
     *     Attempt 2 uses identity and crosses tol=0.3.
     */
    solver->options[LIS_OPTIONS_MAXITER] =
        1;

    /*
     * The final identity-preconditioned Attempt 2 reaches a physical
     * relative residual near 9.51e-2.  On this test system
     * ||b||_2=sqrt(8), corresponding to an absolute residual near
     * 2.69e-1.  A threshold of 3.0e-1 cleanly separates Attempt 2
     * from the two deliberately poor earlier attempts.
     */
    solver->params[
        LIS_PARAMS_RESID
        -LIS_OPTIONS_LEN] =
            (LIS_SCALAR)3.0e-1;

    /*
     * Prime the persistent sequence preconditioner with an exact
     * zero-RHS solve.  The recovery call below must therefore start
     * with preconditioner state inherited from an earlier sequence
     * solve; rebuilding a preconditioner created during the same call
     * is intentionally not a recovery action.
     */
    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    api =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    if(
        api!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS ||
        solver->sequence_precon==NULL ||
        fgcrodr_rebuild_create_count!=1 ||
        fgcrodr_rebuild_destroy_count!=0 )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery rebuild prime failure: "
            "api=%d status=%d precon=%p create=%d destroy=%d\n",
            (int)api,
            (int)solver->retcode,
            (void *)solver->sequence_precon,
            (int)fgcrodr_rebuild_create_count,
            (int)fgcrodr_rebuild_destroy_count);

        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }

    lis_vector_set_all(
        (LIS_SCALAR)2.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    api =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    status =
        solver->retcode;


    true_residual =
        true_relative_residual(
            A,
            b,
            x,
            r);


    err =
        lis_solver_get_sequence_diagnostics(
            solver,
            &diag);

    if( err )
        goto api_fail;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds,
            &updates,
            &uses,
            &rejects);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_lifecycle_state(
            solver,
            &valid,
            &stored_matrix,
            &sg,
            &vg);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_validate_recycle(
            solver,
            &orth,
            &relation);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR recovery rebuild "
        "api=%d status=%d "
        "attempts=%d retries=%d resets=%d rebuilds=%d "
        "first=%d final=%d warm=%d/%d "
        "dim=%d action=%d "
        "builds=%llu updates=%llu uses=%llu rejects=%llu "
        "precon=%d/%d psolve=%d "
        "reported=%.6e true=%.6e "
        "orth=%.3e AZ-C=%.3e\n",
        (int)api,
        (int)status,
        (int)diag.attempts,
        (int)diag.retries,
        (int)diag.recycle_resets,
        (int)diag.precon_rebuilds,
        (int)diag.first_retcode,
        (int)diag.final_retcode,
        (int)diag.initial_warm_start_kept,
        (int)diag.last_warm_start_kept,
        (int)dim,
        (int)action,
        (unsigned long long)builds,
        (unsigned long long)updates,
        (unsigned long long)uses,
        (unsigned long long)rejects,
        (int)fgcrodr_rebuild_create_count,
        (int)fgcrodr_rebuild_destroy_count,
        (int)fgcrodr_rebuild_psolve_count,
        (double)solver->resid,
        (double)true_residual,
        (double)orth,
        (double)relation);


    /*
     * Exact generic second-tier recovery contract.
     */
    if(
        api!=LIS_SUCCESS ||
        status!=LIS_SUCCESS ||
        diag.attempts!=3 ||
        diag.retries!=2 ||
        diag.recycle_resets!=2 ||
        diag.precon_rebuilds!=1 ||
        diag.first_retcode!=LIS_MAXITER ||
        diag.final_retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery rebuild diagnostics failure: "
            "api=%d status=%d "
            "attempts=%d retries=%d resets=%d rebuilds=%d "
            "first=%d final=%d\n",
            (int)api,
            (int)status,
            (int)diag.attempts,
            (int)diag.retries,
            (int)diag.recycle_resets,
            (int)diag.precon_rebuilds,
            (int)diag.first_retcode,
            (int)diag.final_retcode);

        failed =
            1;
    }


    /*
     * The caller supplied x=0, while both retry attempts operate on
     * progressively improved warm starts.
     */
    if(
        diag.initial_warm_start_kept!=LIS_FALSE ||
        diag.last_warm_start_kept!=LIS_TRUE )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery rebuild warm-start failure: "
            "initial=%d last=%d\n",
            (int)diag.initial_warm_start_kept,
            (int)diag.last_warm_start_kept);

        failed =
            1;
    }


    /*
     * Persistent preconditioner lifecycle:
     *
     * create #1 : zero-RHS priming solve
     * reuse      : recovery Attempt 0
     * reuse      : recovery Attempt 1
     * destroy #1: before recovery Attempt 2
     * create #2 : recovery Attempt 2
     */
    if(
        fgcrodr_rebuild_create_count!=2 ||
        fgcrodr_rebuild_destroy_count!=1 ||
        fgcrodr_rebuild_psolve_count!=3 )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery rebuild preconditioner failure: "
            "create=%d destroy=%d psolve=%d\n",
            (int)fgcrodr_rebuild_create_count,
            (int)fgcrodr_rebuild_destroy_count,
            (int)fgcrodr_rebuild_psolve_count);

        failed =
            1;
    }


    /*
     * Every attempt reaches the harmonic recycle point.
     *
     * Attempt 0 BUILD is discarded.
     * Attempt 1 BUILD is discarded.
     * Attempt 2 creates the final fresh BUILD.
     */
    if(
        builds!=3 ||
        updates!=0 ||
        uses!=0 ||
        rejects!=0 ||
        dim!=1 ||
        action!=LIS_FGCRODR_RECYCLE_BUILD )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery rebuild recycle failure: "
            "builds=%llu updates=%llu uses=%llu rejects=%llu "
            "dim=%d action=%d\n",
            (unsigned long long)builds,
            (unsigned long long)updates,
            (unsigned long long)uses,
            (unsigned long long)rejects,
            (int)dim,
            (int)action);

        failed =
            1;
    }


    /*
     * Final recycle relation belongs to Attempt 2 and the rebuilt
     * preconditioner.
     */
    if(
        !valid ||
        stored_matrix!=A ||
        true_residual>2.0e-1 ||
        orth>1.0e-12 ||
        relation>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery rebuild final-state failure: "
            "valid=%d matrix=%p/%p "
            "true=%.17e orth=%.17e AZ-C=%.17e "
            "sg=%llu vg=%llu\n",
            (int)valid,
            (void *)stored_matrix,
            (void *)A,
            (double)true_residual,
            (double)orth,
            (double)relation,
            (unsigned long long)sg,
            (unsigned long long)vg);

        failed =
            1;
    }


    /*
     * Lifetime diagnostics include the one-attempt zero-RHS priming
     * solve plus this three-attempt recovery solve.
     */
    if(
        diag.total_solves!=2 ||
        diag.total_attempts!=4 ||
        diag.total_retries!=2 ||
        diag.total_recycle_resets!=2 ||
        diag.total_precon_rebuilds!=1 )
    {
        fprintf(
            stderr,
            "FGCRO-DR recovery rebuild lifetime failure: "
            "%llu/%llu/%llu/%llu/%llu\n",
            (unsigned long long)diag.total_solves,
            (unsigned long long)diag.total_attempts,
            (unsigned long long)diag.total_retries,
            (unsigned long long)diag.total_recycle_resets,
            (unsigned long long)diag.total_precon_rebuilds);

        failed =
            1;
    }


    if( failed )
    {
        fprintf(
            stderr,
            "FGCRO-DR SEQUENCE RECOVERY REBUILD FAILED\n");

        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    printf(
        "FGCRO-DR SEQUENCE RECOVERY REBUILD PASSED\n");


    err =
        LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR sequence recovery rebuild API failure: %d\n",
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( r )
        lis_vector_destroy(
            r);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}




static LIS_INT
check_correctness_edge_cases(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR u = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT failed;

    LIS_INT api;
    LIS_INT status;
    LIS_INT dim;
    LIS_INT action;

    LIS_INT dim_before;
    LIS_INT action_before;

    LIS_INT iter;

    LIS_REAL true_residual;
    LIS_REAL diff_norm;

    LIS_UNSIGNED_INT builds_before;
    LIS_UNSIGNED_INT builds_after;

    LIS_UNSIGNED_INT updates_before;
    LIS_UNSIGNED_INT updates_after;

    LIS_UNSIGNED_INT uses_before;
    LIS_UNSIGNED_INT uses_after;

    LIS_UNSIGNED_INT rejects_before;
    LIS_UNSIGNED_INT rejects_after;


    failed =
        0;


    err =
        build_recovery_matrix(
            &A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

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
            &u);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &r);

    if( err )
        goto api_fail;


    /*
     * ========================================================
     * Case 1: zero RHS, zero initial guess
     *
     * Must exit before Arnoldi and before any recycle BUILD.
     *
     * LIS reports this initially converged public solve as iter=1.
     * ========================================================
     */

    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 20 "
            "-tol 1e-12 "
            "-fgcrodr_k 2",
            solver);

    if( err )
        goto api_fail;


    lis_vector_set_all(
        0.0,
        b);

    lis_vector_set_all(
        0.0,
        x);


    api =
        lis_solve(
            A,
            b,
            x,
            solver);

    status =
        solver->retcode;

    iter =
        solver->iter;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR edge zero-rhs "
        "api=%d status=%d iter=%d "
        "resid=%.3e dim=%d action=%d\n",
        (int)api,
        (int)status,
        (int)iter,
        (double)solver->resid,
        (int)dim,
        (int)action);


    if(
        api!=LIS_SUCCESS ||
        status!=LIS_SUCCESS ||
        iter!=1 ||
        dim!=0 ||
        action!=LIS_FGCRODR_RECYCLE_NONE ||
        solver->resid!=solver->resid )
    {
        fprintf(
            stderr,
            "FGCRO-DR EDGE zero RHS failure\n");

        failed =
            1;
    }


    lis_solver_destroy(
        solver);

    solver =
        NULL;


    /*
     * ========================================================
     * Case 2: existing recycle + exact initial solution
     *
     * The exact x must return without consuming, replacing or
     * updating the already valid recycle space.  The persistent
     * basis remains valid, while the current-solve recycle action
     * is NONE.
     * ========================================================
     */

    lis_vector_set_all(
        0.0,
        u);

    u->value[0] =
        (LIS_SCALAR)1.0;

    u->value[1] =
        (LIS_SCALAR)(2.0/3.0);


    lis_matvec(
        A,
        u,
        b);


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 100 "
            "-tol 1e-12 "
            "-fgcrodr_k 1",
            solver);

    if( err )
        goto api_fail;


    lis_vector_set_all(
        0.0,
        x);


    api =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    if(
        api!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "FGCRO-DR EDGE setup solve failure\n");

        failed =
            1;

        goto edge_cleanup;
    }


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim_before,
            &action_before,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds_before,
            &updates_before,
            &uses_before,
            &rejects_before);

    if( err )
        goto api_fail;


    lis_vector_copy(
        u,
        x);


    api =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    status =
        solver->retcode;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds_after,
            &updates_after,
            &uses_after,
            &rejects_after);

    if( err )
        goto api_fail;


    lis_vector_axpy(
        -1.0,
        u,
        x);

    lis_vector_nrm2(
        x,
        &diff_norm);


    printf(
        "FGCRO-DR edge exact-x "
        "api=%d status=%d iter=%d "
        "dim=%d/%d action=%d/%d "
        "builds=%llu/%llu updates=%llu/%llu "
        "uses=%llu/%llu rejects=%llu/%llu "
        "error=%.3e\n",
        (int)api,
        (int)status,
        (int)solver->iter,
        (int)dim_before,
        (int)dim,
        (int)action_before,
        (int)action,
        (unsigned long long)builds_before,
        (unsigned long long)builds_after,
        (unsigned long long)updates_before,
        (unsigned long long)updates_after,
        (unsigned long long)uses_before,
        (unsigned long long)uses_after,
        (unsigned long long)rejects_before,
        (unsigned long long)rejects_after,
        (double)diff_norm);


    if(
        api!=LIS_SUCCESS ||
        status!=LIS_SUCCESS ||
        solver->iter!=1 ||
        dim_before!=1 ||
        dim!=dim_before ||
        action!=LIS_FGCRODR_RECYCLE_NONE ||
        builds_after!=builds_before ||
        updates_after!=updates_before ||
        uses_after!=uses_before ||
        rejects_after!=rejects_before ||
        diff_norm>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR EDGE exact initial x failure\n");

        failed =
            1;
    }


edge_cleanup:

    if( solver )
    {
        lis_solver_destroy(
            solver);

        solver =
            NULL;
    }


    /*
     * ========================================================
     * Case 3: one-step lucky/happy Arnoldi breakdown
     *
     * For A=diag(2,3), b=[2,0], x0=0 the first Krylov vector is
     * an exact eigenvector.  h(2,1)=0 is therefore a happy
     * breakdown and must produce SUCCESS, not LIS_BREAKDOWN.
     * ========================================================
     */

    lis_vector_set_all(
        0.0,
        b);

    b->value[0] =
        (LIS_SCALAR)2.0;


    lis_vector_set_all(
        0.0,
        x);


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 20 "
            "-tol 1e-12 "
            "-fgcrodr_k 1",
            solver);

    if( err )
        goto api_fail;


    api =
        lis_solve(
            A,
            b,
            x,
            solver);

    status =
        solver->retcode;


    true_residual =
        true_relative_residual(
            A,
            b,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR edge happy-breakdown "
        "api=%d status=%d iter=%d "
        "true=%.3e dim=%d action=%d\n",
        (int)api,
        (int)status,
        (int)solver->iter,
        (double)true_residual,
        (int)dim,
        (int)action);


    if(
        api!=LIS_SUCCESS ||
        status!=LIS_SUCCESS ||
        solver->iter!=1 ||
        true_residual>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR EDGE happy breakdown failure\n");

        failed =
            1;
    }


    lis_solver_destroy(
        solver);

    solver =
        NULL;


    /*
     * ========================================================
     * Case 4: requested k larger than restart
     *
     * Effective recycle dimension must never exceed restart-1.
     * ========================================================
     */

    lis_vector_set_all(
        (LIS_SCALAR)2.0,
        b);

    lis_vector_set_all(
        0.0,
        x);


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 3 "
            "-maxiter 100 "
            "-tol 1e-12 "
            "-fgcrodr_k 20",
            solver);

    if( err )
        goto api_fail;


    api =
        lis_solve(
            A,
            b,
            x,
            solver);

    status =
        solver->retcode;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    true_residual =
        true_relative_residual(
            A,
            b,
            x,
            r);


    printf(
        "FGCRO-DR edge k>=restart "
        "api=%d status=%d iter=%d "
        "dim=%d restart=3 true=%.3e action=%d\n",
        (int)api,
        (int)status,
        (int)solver->iter,
        (int)dim,
        (double)true_residual,
        (int)action);


    if(
        api!=LIS_SUCCESS ||
        status!=LIS_SUCCESS ||
        dim>2 ||
        true_residual>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR EDGE k>=restart failure\n");

        failed =
            1;
    }


    if( failed )
    {
        fprintf(
            stderr,
            "FGCRO-DR CORRECTNESS EDGE CASES FAILED\n");

        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    printf(
        "FGCRO-DR CORRECTNESS EDGE CASES PASSED\n");


    err =
        LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR correctness edge API failure: %d\n",
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( r )
        lis_vector_destroy(
            r);

    if( u )
        lis_vector_destroy(
            u);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}




static LIS_INT
check_runtime_configuration_transitions(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR u = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR xref = NULL;
    LIS_VECTOR work = NULL;

    LIS_SOLVER solver = NULL;
    LIS_SOLVER fgmres = NULL;

    LIS_INT err;
    LIS_INT failed;

    LIS_INT api;
    LIS_INT fgmres_api;

    LIS_INT dim;
    LIS_INT action;

    LIS_INT target_dim;

    LIS_INT iter_k0;
    LIS_INT iter_fgmres;

    LIS_REAL keep_quality;
    LIS_REAL reject_quality;

    LIS_REAL true_residual;
    LIS_REAL true_fgmres;

    LIS_REAL diff_norm;
    LIS_REAL ref_norm;
    LIS_REAL relative_difference;

    LIS_UNSIGNED_INT builds0;
    LIS_UNSIGNED_INT builds1;
    LIS_UNSIGNED_INT builds2;
    LIS_UNSIGNED_INT builds3;


    failed = 0;


    err =
        build_matrix(
            &A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &u);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

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
            &xref);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &work);

    if( err )
        goto api_fail;


    set_second_exact_solution(
        u);


    lis_matvec(
        A,
        u,
        b);


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 10 "
            "-maxiter 500 "
            "-tol 1e-12 "
            "-fgcrodr_k 5",
            solver);

    if( err )
        goto api_fail;


    /*
     * ========================================================
     * k = 5
     * ========================================================
     */

    lis_vector_set_all(
        0.0,
        x);


    api =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);


    true_residual =
        true_relative_residual(
            A,
            b,
            x,
            work);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds0,
            NULL,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR config k=5 "
        "api=%d status=%d iter=%d "
        "dim=%d action=%d builds=%llu true=%.3e\n",
        (int)api,
        (int)solver->retcode,
        (int)solver->iter,
        (int)dim,
        (int)action,
        (unsigned long long)builds0,
        (double)true_residual);


    if(
        api!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS ||
        dim!=5 ||
        action!=LIS_FGCRODR_RECYCLE_BUILD ||
        builds0!=1 ||
        true_residual>1.0e-10 )
    {
        fprintf(
            stderr,
            "FGCRO-DR CONFIG k=5 failure\n");

        failed = 1;
    }


    /*
     * ========================================================
     * k = 5 -> 2
     *
     * Changing the requested recycle dimension invalidates the
     * existing persistent basis immediately.
     * ========================================================
     */

    err =
        lis_solver_set_fgcrodr_recycle(
            solver,
            2,
            0.90,
            0.10);

    if( err )
        goto api_fail;


    err =
        lis_solver_get_fgcrodr_recycle_config(
            solver,
            &target_dim,
            &keep_quality,
            &reject_quality,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR config 5->2 clear "
        "target=%d dim=%d action=%d\n",
        (int)target_dim,
        (int)dim,
        (int)action);


    if(
        target_dim!=2 ||
        dim!=0 ||
        action!=LIS_FGCRODR_RECYCLE_NONE )
    {
        fprintf(
            stderr,
            "FGCRO-DR CONFIG 5->2 clear failure\n");

        failed = 1;
    }


    lis_vector_set_all(
        0.0,
        x);


    api =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);


    true_residual =
        true_relative_residual(
            A,
            b,
            x,
            work);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds1,
            NULL,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR config k=2 "
        "api=%d status=%d iter=%d "
        "dim=%d action=%d builds=%llu true=%.3e\n",
        (int)api,
        (int)solver->retcode,
        (int)solver->iter,
        (int)dim,
        (int)action,
        (unsigned long long)builds1,
        (double)true_residual);


    if(
        api!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS ||
        dim!=2 ||
        action!=LIS_FGCRODR_RECYCLE_BUILD ||
        builds1!=builds0+1 ||
        true_residual>1.0e-10 )
    {
        fprintf(
            stderr,
            "FGCRO-DR CONFIG k=2 rebuild failure\n");

        failed = 1;
    }


    /*
     * ========================================================
     * k = 2 -> 0
     *
     * k=0 is required to be the exact FGMRES reference path.
     * No FGCRO recycle state may survive or be rebuilt.
     * ========================================================
     */

    err =
        lis_solver_set_fgcrodr_recycle(
            solver,
            0,
            0.90,
            0.10);

    if( err )
        goto api_fail;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    if(
        dim!=0 ||
        action!=LIS_FGCRODR_RECYCLE_NONE )
    {
        fprintf(
            stderr,
            "FGCRO-DR CONFIG 2->0 clear failure: "
            "dim=%d action=%d\n",
            (int)dim,
            (int)action);

        failed = 1;
    }


    /*
     * Independent FGMRES reference.
     */
    err =
        lis_solver_create(
            &fgmres);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgmres "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 10 "
            "-maxiter 500 "
            "-tol 1e-12",
            fgmres);

    if( err )
        goto api_fail;


    lis_vector_set_all(
        0.0,
        xref);


    fgmres_api =
        lis_solve(
            A,
            b,
            xref,
            fgmres);


    iter_fgmres =
        fgmres->iter;


    true_fgmres =
        true_relative_residual(
            A,
            b,
            xref,
            work);


    /*
     * FGCRO-DR with k=0.
     */
    lis_vector_set_all(
        0.0,
        x);


    api =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);


    iter_k0 =
        solver->iter;


    true_residual =
        true_relative_residual(
            A,
            b,
            x,
            work);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds2,
            NULL,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    /*
     * Relative difference between the k=0 result and direct FGMRES.
     */
    lis_vector_copy(
        x,
        work);

    lis_vector_axpy(
        -1.0,
        xref,
        work);

    lis_vector_nrm2(
        work,
        &diff_norm);

    lis_vector_nrm2(
        xref,
        &ref_norm);


    relative_difference =
        ref_norm>0.0
        ? diff_norm/ref_norm
        : diff_norm;


    printf(
        "FGCRO-DR config k=0/FGMRES "
        "api=%d/%d status=%d/%d "
        "iter=%d/%d "
        "dim=%d action=%d builds=%llu "
        "true=%.3e/%.3e diff=%.3e\n",
        (int)api,
        (int)fgmres_api,
        (int)solver->retcode,
        (int)fgmres->retcode,
        (int)iter_k0,
        (int)iter_fgmres,
        (int)dim,
        (int)action,
        (unsigned long long)builds2,
        (double)true_residual,
        (double)true_fgmres,
        (double)relative_difference);


    if(
        api!=LIS_SUCCESS ||
        fgmres_api!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS ||
        fgmres->retcode!=LIS_SUCCESS ||
        dim!=0 ||
        action!=LIS_FGCRODR_RECYCLE_NONE ||
        builds2!=builds1 ||
        iter_k0!=iter_fgmres ||
        true_residual>1.0e-10 ||
        true_fgmres>1.0e-10 ||
        relative_difference>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR CONFIG k=0 FGMRES equivalence failure\n");

        failed = 1;
    }


    lis_solver_destroy(
        fgmres);

    fgmres = NULL;


    /*
     * ========================================================
     * k = 0 -> 3
     *
     * Returning from the reference path must start from an empty
     * recycle state and perform a new BUILD.
     * ========================================================
     */

    err =
        lis_solver_set_fgcrodr_recycle(
            solver,
            3,
            0.90,
            0.10);

    if( err )
        goto api_fail;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR config 0->3 clear "
        "dim=%d action=%d\n",
        (int)dim,
        (int)action);


    if(
        dim!=0 ||
        action!=LIS_FGCRODR_RECYCLE_NONE )
    {
        fprintf(
            stderr,
            "FGCRO-DR CONFIG 0->3 pre-build failure\n");

        failed = 1;
    }


    lis_vector_set_all(
        0.0,
        x);


    api =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);


    true_residual =
        true_relative_residual(
            A,
            b,
            x,
            work);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            solver,
            &builds3,
            NULL,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR config k=3 "
        "api=%d status=%d iter=%d "
        "dim=%d action=%d builds=%llu true=%.3e\n",
        (int)api,
        (int)solver->retcode,
        (int)solver->iter,
        (int)dim,
        (int)action,
        (unsigned long long)builds3,
        (double)true_residual);


    if(
        api!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS ||
        dim!=3 ||
        action!=LIS_FGCRODR_RECYCLE_BUILD ||
        builds3!=builds2+1 ||
        true_residual>1.0e-10 )
    {
        fprintf(
            stderr,
            "FGCRO-DR CONFIG k=3 rebuild failure\n");

        failed = 1;
    }


    if( failed )
    {
        fprintf(
            stderr,
            "FGCRO-DR RUNTIME CONFIGURATION TRANSITIONS FAILED\n");

        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    printf(
        "FGCRO-DR RUNTIME CONFIGURATION TRANSITIONS PASSED\n");


    err =
        LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR runtime configuration API failure: %d\n",
        (int)err);


cleanup:

    if( fgmres )
        lis_solver_destroy(
            fgmres);

    if( solver )
        lis_solver_destroy(
            solver);

    if( work )
        lis_vector_destroy(
            work);

    if( xref )
        lis_vector_destroy(
            xref);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( u )
        lis_vector_destroy(
            u);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}




static LIS_INT
build_scalar_matrix(
    LIS_MATRIX *Aout)
{
    LIS_MATRIX A = NULL;

    LIS_INT err;


    if( Aout==NULL )
        return LIS_ERR_ILL_ARG;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        return err;


    err =
        lis_matrix_set_size(
            A,
            0,
            1);

    if( err )
        goto fail;


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            0,
            0,
            (LIS_SCALAR)2.0,
            A);

    if( err )
        goto fail;


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if( err )
        goto fail;


    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto fail;


    *Aout =
        A;


    return LIS_SUCCESS;


fail:

    lis_matrix_destroy(
        A);

    return err;
}


#ifndef _COMPLEX

static LIS_INT
run_scaled_givens_case(
    LIS_SCALAR scale,
    const char *name)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT api;

    LIS_REAL true_residual;
    LIS_REAL xabs;

    LIS_INT finite;
    LIS_INT failed;


    failed =
        0;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_size(
            A,
            0,
            1);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            0,
            0,
            scale,
            A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if( err )
        goto api_fail;


    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &x);

    if( err )
        goto api_fail;


    b->value[0] =
        (LIS_SCALAR)1.0;

    x->value[0] =
        (LIS_SCALAR)0.0;


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 20 "
            "-tol 1e-12 "
            "-fgcrodr_k 1",
            solver);

    if( err )
        goto api_fail;


    api =
        lis_solve(
            A,
            b,
            x,
            solver);


    true_residual =
        fabs(
            (LIS_REAL)(
                1.0
                - scale*x->value[0]));

    xabs =
        fabs(
            (LIS_REAL)x->value[0]);


    finite =
        (
            xabs==xabs &&
            xabs<=LIS_SCALAR_MAX
        )
        ? 1
        : 0;


    printf(
        "FGCRO-DR scaled %-9s "
        "scale=% .3e api=%d status=%d iter=%d "
        "x=% .17e true=% .17e finite=%d\n",
        name,
        (double)scale,
        (int)api,
        (int)solver->retcode,
        (int)solver->iter,
        (double)x->value[0],
        (double)true_residual,
        (int)finite);


    if(
        api!=LIS_SUCCESS ||
        solver->retcode!=LIS_SUCCESS ||
        solver->iter!=1 ||
        true_residual>1.0e-12 ||
        !finite )
    {
        failed =
            1;
    }


    err =
        failed
        ? LIS_ERR_ILL_ARG
        : LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR scaled %s API failure: %d\n",
        name,
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}


static LIS_INT
fgcrodr_backsub_create(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    (void)solver;
    (void)precon;


    return LIS_SUCCESS;
}


static LIS_INT
fgcrodr_backsub_psolve(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    LIS_INT err;


    (void)solver;


    err =
        lis_vector_copy(
            b,
            x);

    if( err )
        return err;


    /*
     * Keep the physical flexible direction finite while making
     * the reduced triangular coefficient unrepresentable.
     *
     * For A = [1]:
     *
     *     z = 1e-310 v
     *     H = 1e-310
     *     y = 1 / H -> overflow
     *
     * The solver must detect the nonfinite reduced coefficient
     * before applying x += Z y.
     */
    lis_vector_scale(
        (LIS_SCALAR)1.0e-310,
        x);


    return LIS_SUCCESS;
}


static LIS_INT
fgcrodr_backsub_destroy(
    LIS_PRECON precon)
{
    (void)precon;


    return LIS_SUCCESS;
}


static LIS_INT
check_backsub_nonfinite_guard(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT api;

    LIS_REAL xabs;
    LIS_REAL true_residual;

    LIS_INT finite;
    LIS_INT failed;


    failed =
        0;


    err =
        lis_precon_register_ex(
            "fgc2ic",
            fgcrodr_backsub_create,
            fgcrodr_backsub_psolve,
            fgcrodr_backsub_psolve,
            fgcrodr_backsub_destroy);

    if( err )
        goto api_fail;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_size(
            A,
            0,
            1);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            0,
            0,
            (LIS_SCALAR)1.0,
            A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if( err )
        goto api_fail;


    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &x);

    if( err )
        goto api_fail;


    b->value[0] =
        (LIS_SCALAR)1.0;

    x->value[0] =
        (LIS_SCALAR)0.0;


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p fgc2ic "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 20 "
            "-tol 1e-12 "
            "-fgcrodr_k 1",
            solver);

    if( err )
        goto api_fail;


    api =
        lis_solve(
            A,
            b,
            x,
            solver);


    xabs =
        fabs(
            (LIS_REAL)x->value[0]);


    finite =
        (
            xabs==xabs &&
            xabs<=LIS_SCALAR_MAX
        )
        ? 1
        : 0;


    if( finite )
    {
        true_residual =
            fabs(
                1.0
                -
                (LIS_REAL)x->value[0]);
    }
    else
    {
        true_residual =
            LIS_SCALAR_MAX;
    }


    printf(
        "FGCRO-DR backsub nonfinite "
        "api=%d status=%d iter=%d "
        "x=% .17e finite=%d true=% .17e\n",
        (int)api,
        (int)solver->retcode,
        (int)solver->iter,
        (double)x->value[0],
        (int)finite,
        (double)true_residual);


    /*
     * lis_solve() reports API execution success while the solver
     * status carries the numerical breakdown.
     *
     * Most importantly, x must remain untouched and finite.
     */
    if(
        api!=LIS_SUCCESS ||
        solver->retcode!=LIS_BREAKDOWN ||
        !finite ||
        xabs!=0.0 ||
        true_residual!=1.0 )
    {
        fprintf(
            stderr,
            "FGCRO-DR BACKSUB NONFINITE GUARD FAILED\n");

        failed =
            1;
    }


    if( !failed )
    {
        printf(
            "FGCRO-DR BACKSUB NONFINITE GUARD PASSED\n");
    }


    err =
        failed
        ? LIS_ERR_ILL_ARG
        : LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR backsub nonfinite API failure: %d\n",
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}


static LIS_INT
fgcrodr_backsub_earlier_psolve_count = 0;


static LIS_INT
fgcrodr_backsub_earlier_create(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    (void)solver;
    (void)precon;


    fgcrodr_backsub_earlier_psolve_count =
        0;


    return LIS_SUCCESS;
}


static LIS_INT
fgcrodr_backsub_earlier_psolve(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    LIS_INT err;


    (void)solver;


    fgcrodr_backsub_earlier_psolve_count++;


    if(
        fgcrodr_backsub_earlier_psolve_count
        ==
        1 )
    {
        /*
         * Initial Arnoldi vector is e0.
         *
         * Use
         *
         *     z0 = 1e-308 * e0.
         *
         * With the matrix used by the regression:
         *
         *     A*z0 ~= e1.
         *
         * The first rotated diagonal is therefore O(1), so this
         * case must continue into a second Arnoldi iteration.
         */
        err =
            lis_vector_set_all(
                (LIS_SCALAR)0.0,
                x);

        if( err )
            return err;


        if( x->n<2 )
            return LIS_ERR_ILL_ARG;


        x->value[0] =
            (LIS_SCALAR)1.0e-308;


        return LIS_SUCCESS;
    }


    /*
     * At the second Arnoldi iteration v1=e1.
     *
     * Identity gives z1=e1.  The reduced triangular system then
     * has a finite final coefficient y1 but an unrepresentable
     * earlier coefficient y0.
     */
    return
        lis_vector_copy(
            b,
            x);
}


static LIS_INT
fgcrodr_backsub_earlier_destroy(
    LIS_PRECON precon)
{
    (void)precon;


    return LIS_SUCCESS;
}


static LIS_INT
check_backsub_earlier_coefficient_guard(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT api;

    LIS_REAL x0abs;
    LIS_REAL x1abs;

    LIS_INT finite;
    LIS_INT failed;


    failed =
        0;


    err =
        lis_precon_register_ex(
            "fgc2ic4",
            fgcrodr_backsub_earlier_create,
            fgcrodr_backsub_earlier_psolve,
            fgcrodr_backsub_earlier_psolve,
            fgcrodr_backsub_earlier_destroy);

    if( err )
        goto api_fail;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_size(
            A,
            0,
            2);

    if( err )
        goto api_fail;


    /*
     * A =
     *
     *     [ 0      0.5   ]
     *     [ 1e308  1e308 ]
     *
     * The matrix is nonsingular:
     *
     *     det(A) = -0.5e308.
     *
     * For b=[1,0], the exact physical solution is [-2,2].
     *
     * The flexible directions are deliberately chosen so the
     * reduced solve instead reaches
     *
     *     R00 ~= 1
     *     R01 ~= 1e308
     *     R11 ~= -0.5
     *
     * with finite y1 but overflowing earlier y0.
     */

    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            0,
            1,
            (LIS_SCALAR)0.5,
            A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            1,
            0,
            (LIS_SCALAR)1.0e308,
            A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            1,
            1,
            (LIS_SCALAR)1.0e308,
            A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if( err )
        goto api_fail;


    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &x);

    if( err )
        goto api_fail;


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    b->value[0] =
        (LIS_SCALAR)1.0;


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p fgc2ic4 "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 20 "
            "-tol 1e-12 "
            "-fgcrodr_k 1",
            solver);

    if( err )
        goto api_fail;


    api =
        lis_solve(
            A,
            b,
            x,
            solver);


    x0abs =
        fabs(
            (LIS_REAL)x->value[0]);

    x1abs =
        fabs(
            (LIS_REAL)x->value[1]);


    finite =
        (
            x0abs==x0abs &&
            x1abs==x1abs &&
            x0abs<=LIS_SCALAR_MAX &&
            x1abs<=LIS_SCALAR_MAX
        )
        ? 1
        : 0;


    printf(
        "FGCRO-DR backsub earlier "
        "api=%d status=%d iter=%d psolve=%d "
        "x=[% .17e,% .17e] finite=%d\n",
        (int)api,
        (int)solver->retcode,
        (int)solver->iter,
        (int)fgcrodr_backsub_earlier_psolve_count,
        (double)x->value[0],
        (double)x->value[1],
        (int)finite);


    /*
     * iter=2 and psolve=2 prove that the first reduced
     * coefficient was accepted and the failure occurred only
     * while back-substituting the earlier coefficient.
     *
     * x must remain untouched.
     */
    if(
        api!=LIS_SUCCESS ||
        solver->retcode!=LIS_BREAKDOWN ||
        solver->iter!=2 ||
        fgcrodr_backsub_earlier_psolve_count!=2 ||
        !finite ||
        x0abs!=0.0 ||
        x1abs!=0.0 )
    {
        fprintf(
            stderr,
            "FGCRO-DR BACKSUB EARLIER-COEFFICIENT GUARD FAILED\n");

        failed =
            1;
    }


    if( !failed )
    {
        printf(
            "FGCRO-DR BACKSUB EARLIER-COEFFICIENT GUARD PASSED\n");
    }


    err =
        failed
        ? LIS_ERR_ILL_ARG
        : LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR backsub earlier API failure: %d\n",
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}


static LIS_INT
check_scaled_givens_robustness(void)
{
    LIS_INT err;

    LIS_INT failed;


    failed =
        0;


    /*
     * Direct evaluation of
     *
     *     sqrt(a*a + b*b)
     *
     * in the Givens construction can overflow or underflow even
     * when the rotation and the linear-system solution are finite.
     *
     * Both systems below are nonsingular:
     *
     *     1e+200 * x = 1  ->  x = 1e-200
     *     1e-200 * x = 1  ->  x = 1e+200
     */

    err =
        run_scaled_givens_case(
            (LIS_SCALAR)1.0e200,
            "overflow");

    if( err!=LIS_SUCCESS )
    {
        failed =
            1;
    }


    err =
        run_scaled_givens_case(
            (LIS_SCALAR)1.0e-200,
            "underflow");

    if( err!=LIS_SUCCESS )
    {
        failed =
            1;
    }


    if( failed )
    {
        fprintf(
            stderr,
            "FGCRO-DR SCALED GIVENS ROBUSTNESS FAILED\n");

        return LIS_ERR_ILL_ARG;
    }


    printf(
        "FGCRO-DR SCALED GIVENS ROBUSTNESS PASSED\n");


    return LIS_SUCCESS;
}

#endif


static LIS_INT
check_small_dimension_and_restart_guard(void)
{
    LIS_MATRIX A = NULL;
    LIS_MATRIX A1 = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR xref = NULL;
    LIS_VECTOR work = NULL;

    LIS_VECTOR b1 = NULL;
    LIS_VECTOR x1 = NULL;
    LIS_VECTOR r1 = NULL;

    LIS_SOLVER bad = NULL;
    LIS_SOLVER k0 = NULL;
    LIS_SOLVER fgmres = NULL;
    LIS_SOLVER small = NULL;

    LIS_INT err;
    LIS_INT failed;

    LIS_INT option_rc;
    LIS_INT api;

    LIS_INT api_k0;
    LIS_INT api_fgmres;

    LIS_INT dim;
    LIS_INT action;

    LIS_REAL true_k0;
    LIS_REAL true_fgmres;
    LIS_REAL true_small;

    LIS_REAL diff_norm;
    LIS_REAL ref_norm;
    LIS_REAL relative_difference;

    LIS_REAL orth;
    LIS_REAL relation;

    LIS_UNSIGNED_INT builds;


    failed =
        0;


    /*
     * ========================================================
     * Case 1:
     *
     * Recycling requires restart >= 2.
     * ========================================================
     */

    err =
        build_recovery_matrix(
            &A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

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
            &xref);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &work);

    if( err )
        goto api_fail;


    lis_vector_set_all(
        (LIS_SCALAR)2.0,
        b);

    lis_vector_set_all(
        0.0,
        x);


    err =
        lis_solver_create(
            &bad);

    if( err )
        goto api_fail;


    option_rc =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 1 "
            "-maxiter 100 "
            "-tol 1e-12 "
            "-fgcrodr_k 1",
            bad);


    if( option_rc==LIS_SUCCESS )
    {
        api =
            lis_solve(
                A,
                b,
                x,
                bad);
    }
    else
    {
        api =
            option_rc;
    }


    err =
        lis_solver_get_fgcrodr_recycle(
            bad,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR edge restart=1,k=1 "
        "option=%d api=%d status=%d "
        "dim=%d action=%d\n",
        (int)option_rc,
        (int)api,
        (int)bad->retcode,
        (int)dim,
        (int)action);


    if(
        api!=LIS_ERR_ILL_ARG ||
        dim!=0 ||
        action!=LIS_FGCRODR_RECYCLE_NONE )
    {
        fprintf(
            stderr,
            "FGCRO-DR EDGE restart=1 recycling guard failure\n");

        failed =
            1;
    }


    lis_solver_destroy(
        bad);

    bad =
        NULL;


    /*
     * ========================================================
     * Case 2:
     *
     * restart=1 itself is valid for the k=0 reference path.
     * FGCRO-DR(k=0) must remain exactly FGMRES.
     * ========================================================
     */

    err =
        lis_solver_create(
            &k0);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 1 "
            "-maxiter 500 "
            "-tol 1e-12 "
            "-fgcrodr_k 0",
            k0);

    if( err )
        goto api_fail;


    err =
        lis_solver_create(
            &fgmres);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgmres "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 1 "
            "-maxiter 500 "
            "-tol 1e-12",
            fgmres);

    if( err )
        goto api_fail;


    lis_vector_set_all(
        0.0,
        x);

    lis_vector_set_all(
        0.0,
        xref);


    api_k0 =
        lis_solve(
            A,
            b,
            x,
            k0);


    api_fgmres =
        lis_solve(
            A,
            b,
            xref,
            fgmres);


    true_k0 =
        true_relative_residual(
            A,
            b,
            x,
            work);


    true_fgmres =
        true_relative_residual(
            A,
            b,
            xref,
            work);


    lis_vector_copy(
        x,
        work);

    lis_vector_axpy(
        -1.0,
        xref,
        work);

    lis_vector_nrm2(
        work,
        &diff_norm);

    lis_vector_nrm2(
        xref,
        &ref_norm);


    relative_difference =
        ref_norm>0.0
        ? diff_norm/ref_norm
        : diff_norm;


    err =
        lis_solver_get_fgcrodr_recycle(
            k0,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR edge restart=1,k=0/FGMRES "
        "api=%d/%d status=%d/%d "
        "iter=%d/%d dim=%d action=%d "
        "true=%.3e/%.3e diff=%.3e\n",
        (int)api_k0,
        (int)api_fgmres,
        (int)k0->retcode,
        (int)fgmres->retcode,
        (int)k0->iter,
        (int)fgmres->iter,
        (int)dim,
        (int)action,
        (double)true_k0,
        (double)true_fgmres,
        (double)relative_difference);


    if(
        api_k0!=LIS_SUCCESS ||
        api_fgmres!=LIS_SUCCESS ||
        k0->retcode!=LIS_SUCCESS ||
        fgmres->retcode!=LIS_SUCCESS ||
        k0->iter!=fgmres->iter ||
        dim!=0 ||
        action!=LIS_FGCRODR_RECYCLE_NONE ||
        true_k0>1.0e-10 ||
        true_fgmres>1.0e-10 ||
        relative_difference>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR EDGE restart=1 k=0 equivalence failure\n");

        failed =
            1;
    }


    lis_solver_destroy(
        fgmres);

    fgmres =
        NULL;


    lis_solver_destroy(
        k0);

    k0 =
        NULL;


    /*
     * ========================================================
     * Case 3:
     *
     * One-dimensional system with an oversized requested recycle
     * dimension.  Arnoldi terminates after one exact direction.
     *
     * The recycle implementation must accept this short Krylov
     * basis and reduce the effective dimension to one.
     * ========================================================
     */

    err =
        build_scalar_matrix(
            &A1);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A1,
            &b1);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A1,
            &x1);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A1,
            &r1);

    if( err )
        goto api_fail;


    b1->value[0] =
        (LIS_SCALAR)2.0;

    x1->value[0] =
        (LIS_SCALAR)0.0;


    err =
        lis_solver_create(
            &small);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 20 "
            "-tol 1e-12 "
            "-fgcrodr_k 8",
            small);

    if( err )
        goto api_fail;


    api =
        lis_solve(
            A1,
            b1,
            x1,
            small);


    true_small =
        true_relative_residual(
            A1,
            b1,
            x1,
            r1);


    err =
        lis_solver_get_fgcrodr_recycle(
            small,
            &dim,
            &action,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_get_stats(
            small,
            &builds,
            NULL,
            NULL,
            NULL);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_validate_recycle(
            small,
            &orth,
            &relation);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR edge 1x1 "
        "api=%d status=%d iter=%d "
        "dim=%d action=%d builds=%llu "
        "true=%.3e orth=%.3e AZ-C=%.3e\n",
        (int)api,
        (int)small->retcode,
        (int)small->iter,
        (int)dim,
        (int)action,
        (unsigned long long)builds,
        (double)true_small,
        (double)orth,
        (double)relation);


    if(
        api!=LIS_SUCCESS ||
        small->retcode!=LIS_SUCCESS ||
        small->iter!=1 ||
        dim!=1 ||
        action!=LIS_FGCRODR_RECYCLE_BUILD ||
        builds!=1 ||
        true_small>1.0e-12 ||
        orth>1.0e-12 ||
        relation>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR EDGE 1x1 short-basis failure\n");

        failed =
            1;
    }


    if( failed )
    {
        fprintf(
            stderr,
            "FGCRO-DR SMALL-DIMENSION/RESTART GUARDS FAILED\n");

        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    printf(
        "FGCRO-DR SMALL-DIMENSION/RESTART GUARDS PASSED\n");


    err =
        LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR small-dimension/restart API failure: %d\n",
        (int)err);


cleanup:

    if( small )
        lis_solver_destroy(
            small);

    if( fgmres )
        lis_solver_destroy(
            fgmres);

    if( k0 )
        lis_solver_destroy(
            k0);

    if( bad )
        lis_solver_destroy(
            bad);

    if( r1 )
        lis_vector_destroy(
            r1);

    if( x1 )
        lis_vector_destroy(
            x1);

    if( b1 )
        lis_vector_destroy(
            b1);

    if( A1 )
        lis_matrix_destroy(
            A1);

    if( work )
        lis_vector_destroy(
            work);

    if( xref )
        lis_vector_destroy(
            xref);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}


#ifndef _COMPLEX

static LIS_INT
check_active_recycle_noimp_progress(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT api1;
    LIS_INT api2;

    LIS_INT dim1;
    LIS_INT dim2;

    LIS_INT action1;
    LIS_INT action2;

    LIS_INT valid1;
    LIS_INT valid2;

    LIS_MATRIX stored1;
    LIS_MATRIX stored2;

    LIS_UNSIGNED_INT sg1;
    LIS_UNSIGNED_INT vg1;
    LIS_UNSIGNED_INT sg2;
    LIS_UNSIGNED_INT vg2;

    LIS_UNSIGNED_INT builds1;
    LIS_UNSIGNED_INT updates1;
    LIS_UNSIGNED_INT uses1;
    LIS_UNSIGNED_INT rejects1;

    LIS_UNSIGNED_INT builds2;
    LIS_UNSIGNED_INT updates2;
    LIS_UNSIGNED_INT uses2;
    LIS_UNSIGNED_INT rejects2;

    LIS_REAL true1;
    LIS_REAL true2;
    LIS_REAL xnorm2;

    LIS_REAL orth1;
    LIS_REAL relation1;
    LIS_REAL orth2;
    LIS_REAL relation2;

    LIS_INT failed;


    failed =
        0;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_size(
            A,
            0,
            5);

    if( err )
        goto api_fail;


#define SET_A(row,col,val)                                      \
    do                                                          \
    {                                                           \
        err =                                                   \
            lis_matrix_set_value(                               \
                LIS_INS_VALUE,                                  \
                (row),                                          \
                (col),                                          \
                (LIS_SCALAR)(val),                              \
                A);                                             \
                                                                \
        if( err )                                               \
            goto api_fail;                                      \
    }                                                           \
    while(0)


    /*
     * Block 0:
     *
     *     [2]
     *
     * Solve 1 uses b=e0 and creates an exact one-dimensional
     * recycle space.
     */
    SET_A(0,0,2.0);


    /*
     * Independent Stage-2J no-improvement block on indices 1..4:
     *
     *     [ 1  1  1  0 ]
     *     [ 1  1  1  1 ]
     *     [ 0  1  1  0 ]
     *     [ 0  0  1  0 ]
     *
     * Solve 2 uses b=e1.  The existing recycle space is valid and
     * active, but orthogonal to this block.
     */
    SET_A(1,1,1.0);
    SET_A(2,1,1.0);

    SET_A(1,2,1.0);
    SET_A(2,2,1.0);
    SET_A(3,2,1.0);

    SET_A(1,3,1.0);
    SET_A(2,3,1.0);
    SET_A(3,3,1.0);
    SET_A(4,3,1.0);

    SET_A(2,4,1.0);


#undef SET_A


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if( err )
        goto api_fail;


    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

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
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 4 "
            "-maxiter_noimp 0 "
            "-tol 1e-14 "
            "-fgcrodr_k 1 "
            "-fgcrodr_update 1",
            solver);

    if( err )
        goto api_fail;


    /*
     * --------------------------------------------------------
     * Solve 1: create the recycle space.
     * --------------------------------------------------------
     */

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    b->value[0] =
        (LIS_SCALAR)1.0;


    api1 =
        lis_solve(
            A,
            b,
            x,
            solver);


    if( api1!=LIS_SUCCESS )
    {
        err =
            api1;

        goto api_fail;
    }


    err =
        lis_matvec(
            A,
            x,
            r);

    if( err )
        goto api_fail;


    err =
        lis_vector_xpay(
            b,
            (LIS_SCALAR)-1.0,
            r);

    if( err )
        goto api_fail;


    err =
        lis_vector_nrm2(
            r,
            &true1);

    if( err )
        goto api_fail;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim1,
            &action1,
            NULL,
            NULL);

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


    err =
        lis_fgcrodr_get_lifecycle_state(
            solver,
            &valid1,
            &stored1,
            &sg1,
            &vg1);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_validate_recycle(
            solver,
            &orth1,
            &relation1);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR active recycle build "
        "api=%d status=%d iter=%d true=% .17e "
        "dim=%d action=%d valid=%d "
        "builds=%llu updates=%llu uses=%llu rejects=%llu "
        "orth=% .3e AZ-C=% .3e\n",
        (int)api1,
        (int)solver->retcode,
        (int)solver->iter,
        (double)true1,
        (int)dim1,
        (int)action1,
        (int)valid1,
        (unsigned long long)builds1,
        (unsigned long long)updates1,
        (unsigned long long)uses1,
        (unsigned long long)rejects1,
        (double)orth1,
        (double)relation1);


    if(
        solver->retcode!=LIS_SUCCESS ||
        solver->iter!=1 ||
        true1>1.0e-13 ||
        dim1!=1 ||
        action1!=LIS_FGCRODR_RECYCLE_BUILD ||
        !valid1 ||
        stored1!=A ||
        builds1!=1 ||
        updates1!=0 ||
        uses1!=0 ||
        rejects1!=0 ||
        orth1>1.0e-12 ||
        relation1>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR ACTIVE RECYCLE BUILD FAILED\n");

        failed =
            1;

        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    /*
     * --------------------------------------------------------
     * Solve 2: active recycle + deterministic no-improvement.
     * --------------------------------------------------------
     */

    solver->options[
        LIS_OPTIONS_MAXITER] =
            3;

    solver->options[
        LIS_OPTIONS_MAXITER_NO_IMP] =
            1;

    solver->params[
        LIS_PARAMS_RESID
        -LIS_OPTIONS_LEN] =
            (LIS_SCALAR)1.0e-14;


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    b->value[1] =
        (LIS_SCALAR)1.0;


    api2 =
        lis_solve(
            A,
            b,
            x,
            solver);


    if( api2!=LIS_SUCCESS )
    {
        err =
            api2;

        goto api_fail;
    }


    err =
        lis_matvec(
            A,
            x,
            r);

    if( err )
        goto api_fail;


    err =
        lis_vector_xpay(
            b,
            (LIS_SCALAR)-1.0,
            r);

    if( err )
        goto api_fail;


    err =
        lis_vector_nrm2(
            r,
            &true2);

    if( err )
        goto api_fail;


    err =
        lis_vector_nrm2(
            x,
            &xnorm2);

    if( err )
        goto api_fail;


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim2,
            &action2,
            NULL,
            NULL);

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


    err =
        lis_fgcrodr_get_lifecycle_state(
            solver,
            &valid2,
            &stored2,
            &sg2,
            &vg2);

    if( err )
        goto api_fail;


    err =
        lis_fgcrodr_validate_recycle(
            solver,
            &orth2,
            &relation2);

    if( err )
        goto api_fail;


    printf(
        "FGCRO-DR active recycle noimp "
        "api=%d status=%d iter=%d "
        "reported=% .17e true=% .17e xnorm=% .17e "
        "x=[% .17e,% .17e,% .17e,% .17e,% .17e] "
        "dim=%d action=%d valid=%d "
        "builds=%llu updates=%llu uses=%llu rejects=%llu "
        "orth=% .3e AZ-C=% .3e\n",
        (int)api2,
        (int)solver->retcode,
        (int)solver->iter,
        (double)solver->resid,
        (double)true2,
        (double)xnorm2,
        (double)x->value[0],
        (double)x->value[1],
        (double)x->value[2],
        (double)x->value[3],
        (double)x->value[4],
        (int)dim2,
        (int)action2,
        (int)valid2,
        (unsigned long long)builds2,
        (unsigned long long)updates2,
        (unsigned long long)uses2,
        (unsigned long long)rejects2,
        (double)orth2,
        (double)relation2);


    if(
        solver->retcode!=LIS_MAXITER ||
        solver->iter!=3 ||
        fabs(
            solver->resid
            -0.7071067811865476)
            >1.0e-12 ||
        fabs(
            true2
            -0.7071067811865476)
            >1.0e-12 ||
        fabs(
            xnorm2
            -0.5)
            >1.0e-12 ||
        fabs(
            (LIS_REAL)x->value[0])
            >1.0e-12 ||
        fabs(
            (LIS_REAL)x->value[1]
            -0.5)
            >1.0e-12 ||
        fabs(
            (LIS_REAL)x->value[2])
            >1.0e-12 ||
        fabs(
            (LIS_REAL)x->value[3])
            >1.0e-12 ||
        fabs(
            (LIS_REAL)x->value[4])
            >1.0e-12 ||
        dim2!=1 ||
        action2!=LIS_FGCRODR_RECYCLE_KEEP ||
        !valid2 ||
        stored2!=A ||
        sg2!=sg1 ||
        vg2!=vg1 ||
        builds2!=builds1 ||
        updates2!=updates1+1 ||
        uses2!=uses1+1 ||
        rejects2!=rejects1 ||
        orth2>1.0e-12 ||
        relation2>1.0e-12 )
    {
        fprintf(
            stderr,
            "FGCRO-DR ACTIVE RECYCLE NO-IMPROVEMENT FAILED\n");

        failed =
            1;
    }
    else
    {
        printf(
            "FGCRO-DR ACTIVE RECYCLE NO-IMPROVEMENT PASSED\n");
    }


    err =
        failed
        ? LIS_ERR_ILL_ARG
        : LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR active recycle no-improvement API failure: %d\n",
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( r )
        lis_vector_destroy(
            r);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}


static LIS_INT
check_reliable_candidate_active_recycle(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err = LIS_SUCCESS;

    LIS_INT api1;
    LIS_INT status1;
    LIS_INT iter1;

    LIS_INT api2;
    LIS_INT status2;
    LIS_INT iter2;

    LIS_INT dim1;
    LIS_INT action1;

    LIS_INT dim2;
    LIS_INT action2;

    LIS_UNSIGNED_INT builds1;
    LIS_UNSIGNED_INT updates1;
    LIS_UNSIGNED_INT uses1;
    LIS_UNSIGNED_INT rejects1;

    LIS_UNSIGNED_INT builds2;
    LIS_UNSIGNED_INT updates2;
    LIS_UNSIGNED_INT uses2;
    LIS_UNSIGNED_INT rejects2;

    LIS_REAL true1;
    LIS_REAL true2;

    LIS_REAL reported2;

    LIS_INT i;

    int failed = 0;

    const LIS_INT row[] = {
        0, 0, 0,
        1,
        2,
        3, 3, 3
    };

    const LIS_INT col[] = {
        0, 2, 3,
        1,
        2,
        1, 2, 3
    };

    const LIS_SCALAR val[] = {
         9.6635805078063965e+11,
         1.0655399464330062e-07,
        -2.9511571135475290e-12,

         8.0844035225373991e+03,

         2.0245180553697102e-05,

         1.1098830243238601e-03,
        -1.9650356386922263e+02,
         1.3661903277464329e+10
    };

    const LIS_SCALAR seed_rhs[] = {
         1.0,
        -1.0,
        -1.0,
        -1.0
    };

    const LIS_SCALAR active_rhs[] = {
        -1.0,
        -1.0,
        -1.0,
        -1.0
    };


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_size(
            A,
            0,
            4);

    if( err )
        goto api_fail;


    for(i=0;i<8;i++)
    {
        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                row[i],
                col[i],
                val[i],
                A);

        if( err )
            goto api_fail;
    }


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if( err )
        goto api_fail;


    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

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
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgcrodr "
            "-fgcrodr_k 1 "
            "-p jacobi "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 32 "
            "-maxiter_noimp 0 "
            "-tol 1e-10",
            solver);

    if( err )
        goto api_fail;


    /*
     * Solve 1 builds the persistent recycle space.
     */
    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    for(i=0;i<4;i++)
        b->value[i] = seed_rhs[i];


    api1 =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    status1 =
        solver->retcode;

    iter1 =
        solver->iter;

    true1 =
        true_relative_residual(
            A,
            b,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim1,
            &action1,
            NULL,
            NULL);

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


    /*
     * Solve 2 starts from zero and consumes the persistent
     * recycle space.  The projected residual is known to produce
     * a premature convergence candidate on the uncorrected core.
     */
    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    for(i=0;i<4;i++)
        b->value[i] = active_rhs[i];


    api2 =
        lis_solve_sequence(
            A,
            b,
            x,
            solver);

    status2 =
        solver->retcode;

    iter2 =
        solver->iter;

    reported2 =
        solver->resid;

    true2 =
        true_relative_residual(
            A,
            b,
            x,
            r);


    err =
        lis_solver_get_fgcrodr_recycle(
            solver,
            &dim2,
            &action2,
            NULL,
            NULL);

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


    printf(
        "FGCRO-DR reliable candidate "
        "seed_api=%d seed_status=%d seed_iter=%d "
        "seed_true=%.17e seed_dim=%d seed_action=%d "
        "api=%d status=%d iter=%d "
        "reported=%.17e true=%.17e "
        "dim=%d action=%d "
        "builds=%llu/%llu "
        "updates=%llu/%llu "
        "uses=%llu/%llu "
        "rejects=%llu/%llu\n",
        (int)api1,
        (int)status1,
        (int)iter1,
        (double)true1,
        (int)dim1,
        (int)action1,
        (int)api2,
        (int)status2,
        (int)iter2,
        (double)reported2,
        (double)true2,
        (int)dim2,
        (int)action2,
        (unsigned long long)builds1,
        (unsigned long long)builds2,
        (unsigned long long)updates1,
        (unsigned long long)updates2,
        (unsigned long long)uses1,
        (unsigned long long)uses2,
        (unsigned long long)rejects1,
        (unsigned long long)rejects2);


    if(
        api1!=LIS_SUCCESS ||
        status1!=LIS_SUCCESS ||
        true1>1.0e-10 ||
        dim1!=1 ||
        action1!=LIS_FGCRODR_RECYCLE_BUILD ||
        builds1!=1 ||
        uses1!=0 )
    {
        fprintf(
            stderr,
            "FGCRO-DR RELIABLE CANDIDATE SEED FAILED\n");

        failed = 1;
    }


    if(
        api2!=LIS_SUCCESS ||
        status2!=LIS_SUCCESS ||
        iter2<=1 ||
        reported2!=reported2 ||
        true2!=true2 ||
        reported2>1.0e-10 ||
        true2>1.0e-10 ||
        dim2!=1 ||
        action2!=LIS_FGCRODR_RECYCLE_KEEP ||
        builds2!=builds1 ||
        uses2<uses1+1 )
    {
        fprintf(
            stderr,
            "FGCRO-DR RELIABLE CANDIDATE FAILED\n");

        failed = 1;
    }


    if( !failed )
    {
        printf(
            "FGCRO-DR RELIABLE CANDIDATE PASSED\n");
    }


    err =
        failed
        ? LIS_ERR_ILL_ARG
        : LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGCRO-DR reliable candidate API failure: %d\n",
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( r )
        lis_vector_destroy(
            r);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}


#endif


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

    if(
        check_generation_transport_reject()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }


    if(
        check_sequence_reset_lifecycle()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }


    if(
        check_sequence_recovery_lifecycle()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }


    if(
        check_sequence_recovery_rebuild_lifecycle()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }


    if(
        check_correctness_edge_cases()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }


    if(
        check_runtime_configuration_transitions()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }


    if(
        check_small_dimension_and_restart_guard()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }


#ifndef _COMPLEX
    if(
        check_scaled_givens_robustness()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }
#endif


#ifndef _COMPLEX
    if(
        check_backsub_nonfinite_guard()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }


    if(
        check_backsub_earlier_coefficient_guard()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }
#endif


#ifndef _COMPLEX
    if(
        check_reliable_candidate_active_recycle()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }


    if(
        check_active_recycle_noimp_progress()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }
#endif


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
#endif /* !_COMPLEX */
