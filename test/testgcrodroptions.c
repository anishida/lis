#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <stdio.h>
#include <math.h>
#include "lis.h"
#include "lis_solver.h"

#ifdef _COMPLEX
int main(int argc, char **argv)
{
    (void)argc;
    (void)argv;
    printf("SKIP: GCRO-DR is not implemented for complex scalar builds\n");
    return 77;
}
#else

static int close_real(LIS_REAL a, LIS_REAL b)
{
    LIS_REAL d = a-b;
    if( d<0 ) d=-d;
    return d < (LIS_REAL)1.0e-12;
}

static int check_config(LIS_SOLVER solver, LIS_INT k,
                        LIS_REAL keep, LIS_REAL reject,
                        const char *label)
{
    if( solver->recycle_target_dim!=k ||
        !close_real(solver->recycle_keep_threshold,keep) ||
        !close_real(solver->recycle_reject_threshold,reject) )
    {
        fprintf(stderr,
                "%s: got k=%d keep=%.17g reject=%.17g\n",
                label,
                (int)solver->recycle_target_dim,
                (double)solver->recycle_keep_threshold,
                (double)solver->recycle_reject_threshold);
        return 1;
    }
    return 0;
}

int main(int argc, char **argv)
{
    LIS_SOLVER solver = NULL;
    LIS_INT err;

    err = lis_initialize(&argc,&argv);
    if( err ) return (int)err;

    err = lis_solver_create(&solver);
    if( err ) goto fail;

    err = lis_solver_set_option(
        "-i gcrodr -restart 30 "
        "-gcrodr_k 6 -gcrodr_keep 0.90 -gcrodr_reject 0.35",
        solver);
    if( err )
    {
        fprintf(stderr,"basic GCRO-DR textual options failed: %d\n",(int)err);
        goto fail;
    }
    if( check_config(solver,6,(LIS_REAL)0.90,(LIS_REAL)0.35,"basic") )
        goto fail_test;

    /*
     * Verify that threshold option ordering is independent.
     * This is important because textual/command-line options are not
     * required to present KEEP before REJECT.
     */
    err = lis_solver_set_option(
        "-gcrodr_reject 0.95 -gcrodr_k 12 -gcrodr_keep 0.99",
        solver);
    if( err )
    {
        fprintf(stderr,"order-independent GCRO-DR options failed: %d\n",(int)err);
        goto fail;
    }
    if( check_config(solver,12,(LIS_REAL)0.99,(LIS_REAL)0.95,"ordered") )
        goto fail_test;

    err = lis_gcrodr_check_params(solver);
    if( err )
    {
        fprintf(stderr,"valid GCRO-DR config rejected: %d\n",(int)err);
        goto fail;
    }

    /* Range validation must reject malformed solver-specific values. */
    err = lis_solver_set_option("-gcrodr_keep 1.25",solver);
    if( err==LIS_SUCCESS )
    {
        fprintf(stderr,"invalid -gcrodr_keep was accepted\n");
        goto fail_test;
    }

    err = lis_solver_set_option("-gcrodr_keep nan",solver);
    if( err==LIS_SUCCESS )
    {
        fprintf(stderr,"NaN -gcrodr_keep was accepted\n");
        goto fail_test;
    }

    err = lis_solver_set_option("-gcrodr_reject nan",solver);
    if( err==LIS_SUCCESS )
    {
        fprintf(stderr,"NaN -gcrodr_reject was accepted\n");
        goto fail_test;
    }

    /* Restore a valid value. */
    err = lis_solver_set_option("-gcrodr_keep 0.99",solver);
    if( err ) goto fail;

    /*
     * Cross-threshold validation is intentionally deferred until the
     * solver parameter check so option order remains independent.
     */
    err = lis_solver_set_option(
        "-gcrodr_keep 0.90 -gcrodr_reject 0.95",
        solver);
    if( err ) goto fail;

    err = lis_gcrodr_check_params(solver);
    if( err==LIS_SUCCESS )
    {
        fprintf(stderr,"reject > keep was not rejected by check_params\n");
        goto fail_test;
    }

    err = lis_solver_set_option(
        "-gcrodr_keep 0.90 -gcrodr_reject 0.35 -gcrodr_k 6",
        solver);
    if( err ) goto fail;
    err = lis_gcrodr_check_params(solver);
    if( err ) goto fail;

    printf("GCRO-DR textual option test: PASS\n");
    printf("k=%d keep=%.2f reject=%.2f\n",
           (int)solver->recycle_target_dim,
           (double)solver->recycle_keep_threshold,
           (double)solver->recycle_reject_threshold);

    lis_solver_destroy(solver);
    lis_finalize();
    return 0;

fail_test:
    err = 1;
fail:
    if( solver ) lis_solver_destroy(solver);
    lis_finalize();
    return (int)(err ? err : 1);
}
#endif /* !_COMPLEX */
