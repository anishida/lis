#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_solver.h"

#include <stdio.h>


/*
 * Internal function implemented in Stage 5B.
 */
extern LIS_INT lis_vector_duplicateex(
    LIS_INT precision,
    void *vin,
    LIS_VECTOR *vout);



static LIS_INT
create_matrix(LIS_INT n, LIS_MATRIX *A)
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
            (LIS_SCALAR)2.0,
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
check_vector_contract(
    LIS_MATRIX A,
    LIS_MATRIX B)
{
    LIS_VECTOR v = NULL;
    LIS_VECTOR q = NULL;
    LIS_INT err;
    int failed = 0;

    err = lis_vector_duplicate(A,&v);
    if( err )
    {
        fprintf(stderr,
                "lis_vector_duplicate failed: %d\n",
                (int)err);

        return 1;
    }


    /*
     * Same matrix layout.
     */
    if( !lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "same layout reported incompatible\n");

        failed = 1;
    }


    /*
     * Different matrix object, same distribution.
     */
    if( !lis_solver_work_vector_compatible(
            v,
            B,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "equivalent matrix layout reported incompatible\n");

        failed = 1;
    }


    /*
     * Precision mismatch.
     */
    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_QUAD) )
    {
        fprintf(stderr,
                "precision mismatch reported compatible\n");

        failed = 1;
    }


    err = lis_vector_duplicateex(
        LIS_PRECISION_QUAD,
        A,
        &q);

    if( err )
    {
        fprintf(stderr,
                "lis_vector_duplicateex failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    if( !lis_solver_work_vector_compatible(
            q,
            A,
            LIS_PRECISION_QUAD) )
    {
        fprintf(stderr,
                "quad vector reported incompatible\n");

        failed = 1;
    }


    /*
     * Structural mismatch tests.
     * Every field is restored before destruction.
     */

    v->n--;

    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "n mismatch reported compatible\n");

        failed = 1;
    }

    v->n++;


    v->np--;

    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "np mismatch reported compatible\n");

        failed = 1;
    }

    v->np++;


    v->pad++;

    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "pad mismatch reported compatible\n");

        failed = 1;
    }

    v->pad--;


    v->is++;

    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "range start mismatch reported compatible\n");

        failed = 1;
    }

    v->is--;


    v->ie++;

    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "range end mismatch reported compatible\n");

        failed = 1;
    }

    v->ie--;


    v->gn--;

    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "global size mismatch reported compatible\n");

        failed = 1;
    }

    v->gn++;


    v->origin = !v->origin;

    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "origin mismatch reported compatible\n");

        failed = 1;
    }

    v->origin = !v->origin;


    v->my_rank++;

    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "rank mismatch reported compatible\n");

        failed = 1;
    }

    v->my_rank--;


    v->nprocs++;

    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "process count mismatch reported compatible\n");

        failed = 1;
    }

    v->nprocs--;


#ifdef USE_MPI

    if( v->ranges!=NULL && A->ranges!=NULL )
    {
        v->ranges[0]++;

        if( lis_solver_work_vector_compatible(
                v,
                A,
                LIS_PRECISION_DEFAULT) )
        {
            fprintf(stderr,
                    "MPI ranges mismatch reported compatible\n");

            failed = 1;
        }

        v->ranges[0]--;
    }

#else

    v->comm++;

    if( lis_solver_work_vector_compatible(
            v,
            A,
            LIS_PRECISION_DEFAULT) )
    {
        fprintf(stderr,
                "communicator mismatch reported compatible\n");

        failed = 1;
    }

    v->comm--;

#endif


cleanup:

    if( q )
    {
        lis_vector_destroy(q);
    }

    if( v )
    {
        lis_vector_destroy(v);
    }

    return failed;
}


static int
check_workspace_contract(
    LIS_MATRIX A)
{
    LIS_SOLVER solver = NULL;
    LIS_INT err;
    LIS_INT saved;
    int failed = 0;

    err = lis_solver_create(&solver);
    if( err )
    {
        fprintf(stderr,
                "lis_solver_create failed: %d\n",
                (int)err);

        return 1;
    }

    err = lis_solver_set_matrix(A,solver);
    if( err )
    {
        fprintf(stderr,
                "lis_solver_set_matrix failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    solver->precision = LIS_PRECISION_DEFAULT;


    /* ====================================================== */
    /* CG: fixed workspace                                    */
    /* ====================================================== */

    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_CG;

    err = lis_cg_malloc_work(solver);

    if( err )
    {
        fprintf(stderr,
                "CG malloc_work failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    if( !lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "CG workspace reported incompatible\n");

        failed = 1;
    }

    /*
     * COCG also requires four ordinary matrix-layout vectors.
     * The workspace contract is structural, not tied to the
     * solver that originally allocated it.
     */
    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_COCG;

    if( !lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "CG workspace not reusable by compatible COCG layout\n");

        failed = 1;
    }

    /*
     * BiCG requires six vectors, so the same work array cannot
     * be compatible.
     */
    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_BICG;

    if( lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "wrong worklen reported compatible for BiCG\n");

        failed = 1;
    }

    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_CG;

    solver->work[0]->n--;

    if( lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "corrupted CG work vector reported compatible\n");

        failed = 1;
    }

    solver->work[0]->n++;

    /*
     * Actual allocated precision must match the current solver
     * precision requirement.
     */
    solver->precision = LIS_PRECISION_QUAD;

    if( lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "workspace precision mismatch reported compatible\n");

        failed = 1;
    }

    solver->precision = LIS_PRECISION_DEFAULT;

    lis_solver_work_destroy(solver);


    /* ====================================================== */
    /* GMRES                                                  */
    /* ====================================================== */

    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_GMRES;

    solver->options[LIS_OPTIONS_RESTART] = 5;

    err = lis_gmres_malloc_work(solver);

    if( err )
    {
        fprintf(stderr,
                "GMRES malloc_work failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    if( !lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "GMRES workspace reported incompatible\n");

        failed = 1;
    }

    saved = solver->options[LIS_OPTIONS_RESTART];

    solver->options[LIS_OPTIONS_RESTART] = 6;

    if( lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "GMRES restart change reported compatible\n");

        failed = 1;
    }

    solver->options[LIS_OPTIONS_RESTART] = saved;


    /*
     * work[0] is the special restart+1 vector.
     */
    solver->work[0]->n--;

    if( lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "GMRES special work[0] mismatch reported compatible\n");

        failed = 1;
    }

    solver->work[0]->n++;

    lis_solver_work_destroy(solver);


    /* ====================================================== */
    /* FGMRES                                                 */
    /* ====================================================== */

    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_FGMRES;

    solver->options[LIS_OPTIONS_RESTART] = 5;

    err = lis_fgmres_malloc_work(solver);

    if( err )
    {
        fprintf(stderr,
                "FGMRES malloc_work failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    if( !lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "FGMRES workspace reported incompatible\n");

        failed = 1;
    }

    solver->options[LIS_OPTIONS_RESTART] = 6;

    if( lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "FGMRES restart change reported compatible\n");

        failed = 1;
    }

    solver->options[LIS_OPTIONS_RESTART] = 5;

    lis_solver_work_destroy(solver);


    /* ====================================================== */
    /* BiCGSTAB(l)                                            */
    /* ====================================================== */

    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_BICGSTABL;

    solver->options[LIS_OPTIONS_ELL] = 2;

    err = lis_bicgstabl_malloc_work(solver);

    if( err )
    {
        fprintf(stderr,
                "BiCGSTAB(l) malloc_work failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    if( !lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "BiCGSTAB(l) workspace reported incompatible\n");

        failed = 1;
    }

    solver->options[LIS_OPTIONS_ELL] = 3;

    if( lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "BiCGSTAB(l) ell change reported compatible\n");

        failed = 1;
    }

    solver->options[LIS_OPTIONS_ELL] = 2;

    lis_solver_work_destroy(solver);


    /* ====================================================== */
    /* Orthomin                                               */
    /* ====================================================== */

    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_ORTHOMIN;

    solver->options[LIS_OPTIONS_RESTART] = 5;

    err = lis_orthomin_malloc_work(solver);

    if( err )
    {
        fprintf(stderr,
                "Orthomin malloc_work failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    if( !lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "Orthomin workspace reported incompatible\n");

        failed = 1;
    }

    solver->options[LIS_OPTIONS_RESTART] = 6;

    if( lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "Orthomin restart change reported compatible\n");

        failed = 1;
    }

    solver->options[LIS_OPTIONS_RESTART] = 5;

    lis_solver_work_destroy(solver);


    /* ====================================================== */
    /* IDR(s)                                                 */
    /* ====================================================== */

    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_IDRS;

    solver->options[LIS_OPTIONS_IDRS_RESTART] = 2;

    err = lis_idrs_malloc_work(solver);

    if( err )
    {
        fprintf(stderr,
                "IDR(s) malloc_work failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    if( !lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "IDR(s) workspace reported incompatible\n");

        failed = 1;
    }

    solver->options[LIS_OPTIONS_IDRS_RESTART] = 3;

    if( lis_solver_work_compatible(solver) )
    {
        fprintf(stderr,
                "IDR(s) dimension change reported compatible\n");

        failed = 1;
    }

    solver->options[LIS_OPTIONS_IDRS_RESTART] = 2;

    lis_solver_work_destroy(solver);


cleanup:

    if( solver )
    {
        lis_solver_destroy(solver);
    }

    return failed;
}


int
main(int argc, char **argv)
{
    LIS_MATRIX A = NULL;
    LIS_MATRIX B = NULL;
    LIS_INT err;
    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if( err ) return 1;


    err = create_matrix(8,&A);

    if( err )
    {
        fprintf(stderr,
                "create_matrix(A) failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }


    err = create_matrix(8,&B);

    if( err )
    {
        fprintf(stderr,
                "create_matrix(B) failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }


    if( check_vector_contract(A,B) )
    {
        failed = 1;
    }


    if( check_workspace_contract(A) )
    {
        failed = 1;
    }


cleanup:

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
        fprintf(stderr,
                "LIS_WORK_COMPAT_STAGE5B FAILED\n");

        return 1;
    }


    printf(
        "LIS_WORK_COMPAT_STAGE5B PASSED\n");

    return 0;
}
