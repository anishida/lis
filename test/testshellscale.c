#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include "lis.h"
#include <stdio.h>


typedef struct
{
    LIS_INT n;
    LIS_INT is;

    LIS_SCALAR even_diag;
    LIS_SCALAR odd_diag;

    LIS_INT capture_input;
    LIS_INT captured_input;
    LIS_SCALAR captured_x0;
} USER_DIAG;


static LIS_SCALAR
diag_value(
    USER_DIAG *A,
    LIS_INT local_i)
{
    LIS_INT global_i = A->is + local_i;

    /*
     * Alternating positive diagonal:
     *
     *   4, 9, 4, 9, ...
     *
     * This deliberately gives different symmetric
     * scaling factors:
     *
     *   1/sqrt(4) = 1/2
     *   1/sqrt(9) = 1/3
     */
    if( (global_i % 2)==0 )
        return A->even_diag;

    return A->odd_diag;
}


static LIS_INT
user_matvec(
    void *ctx,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
    USER_DIAG *A = (USER_DIAG *)ctx;
    LIS_INT i;

    if( A->capture_input &&
        !A->captured_input &&
        A->n>0 )
    {
        A->captured_x0 = x[0];
        A->captured_input = 1;
    }

    for(i=0;i<A->n;i++)
        y[i] = diag_value(A,i)*x[i];

    return LIS_SUCCESS;
}


static LIS_INT
user_matvech(
    void *ctx,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
#ifdef _COMPLEX
    USER_DIAG *A = (USER_DIAG *)ctx;
    LIS_INT i;

    for(i=0;i<A->n;i++)
        y[i] = conj(diag_value(A,i))*x[i];

    return LIS_SUCCESS;
#else
    return user_matvec(ctx,x,y);
#endif
}


static LIS_INT
user_get_diagonal(
    void *ctx,
    LIS_SCALAR *d)
{
    USER_DIAG *A = (USER_DIAG *)ctx;
    LIS_INT i;

    for(i=0;i<A->n;i++)
        d[i] = diag_value(A,i);

    return LIS_SUCCESS;
}


static int
scalar_close(
    LIS_SCALAR got,
    LIS_SCALAR expected)
{
    const LIS_REAL tol = (LIS_REAL)1.0e-9;

#ifdef _COMPLEX
    LIS_REAL dr = (LIS_REAL)creal(got-expected);
    LIS_REAL di = (LIS_REAL)cimag(got-expected);

    if( dr<0 ) dr = -dr;
    if( di<0 ) di = -di;

    return dr<=tol && di<=tol;
#else
    LIS_REAL diff = (LIS_REAL)(got-expected);

    if( diff<0 ) diff = -diff;

    return diff<=tol;
#endif
}


static void
print_value(
    const char *name,
    LIS_INT global_i,
    LIS_SCALAR got,
    LIS_SCALAR expected)
{
#ifdef _COMPLEX
    fprintf(
        stderr,
        "%s mismatch at global index %d: "
        "got=(%.15e,%.15e) expected=(%.15e,%.15e)\n",
        name,
        (int)global_i,
        (double)creal(got),
        (double)cimag(got),
        (double)creal(expected),
        (double)cimag(expected));
#else
    fprintf(
        stderr,
        "%s mismatch at global index %d: "
        "got=%.15e expected=%.15e\n",
        name,
        (int)global_i,
        (double)got,
        (double)expected);
#endif
}


static int
check_constant_vector(
    const char *name,
    LIS_VECTOR v,
    LIS_INT is,
    LIS_INT n,
    LIS_SCALAR expected)
{
    LIS_INT i;
    int failed = 0;

    for(i=0;i<n;i++)
    {
        if( !scalar_close(v->value[i],expected) )
        {
            print_value(
                name,
                is+i,
                v->value[i],
                expected);

            failed = 1;
        }
    }

    return failed;
}


static int
check_scale_vector(
    USER_DIAG *A,
    LIS_VECTOR d,
    LIS_SCALAR even_scale,
    LIS_SCALAR odd_scale)
{
    LIS_INT i;
    int failed = 0;

    if( d==NULL )
    {
        fprintf(
            stderr,
            "solver scaling vector is NULL\n");

        return 1;
    }

    for(i=0;i<A->n;i++)
    {
        LIS_INT global_i = A->is+i;

        LIS_SCALAR expected =
            ((global_i % 2)==0)
            ? even_scale
            : odd_scale;

        if( !scalar_close(
                d->value[i],
                expected) )
        {
            print_value(
                "solver scaling vector",
                global_i,
                d->value[i],
                expected);

            failed = 1;
        }
    }

    return failed;
}


static int
check_rhs(
    USER_DIAG *A,
    LIS_VECTOR b,
    LIS_SCALAR factor)
{
    LIS_INT i;
    int failed = 0;

    for(i=0;i<A->n;i++)
    {
        LIS_SCALAR expected =
            factor*diag_value(A,i);

        if( !scalar_close(
                b->value[i],
                expected) )
        {
            print_value(
                "physical RHS",
                A->is+i,
                b->value[i],
                expected);

            failed = 1;
        }
    }

    return failed;
}


/* STAGE6_REPEATED_DIAG */



int
main(
    int argc,
    char **argv)
{
    LIS_MATRIX A = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_SOLVER solver = NULL;

    USER_DIAG ctx;

    LIS_INT err;
    LIS_INT i,is,ie;
    LIS_INT global_n = 4;

    int failed = 0;


    err = lis_initialize(&argc,&argv);
    if( err )
        return 1;


    /*
     * Matrix:
     *
     *       [ 4  0  0  0 ]
     *       [ 0  9  0  0 ]
     *   A = [ 0  0  4  0 ]
     *       [ 0  0  0  9 ]
     *
     * with the same pattern extended for MPI if needed.
     *
     * Choose
     *
     *       b_i = 2*A_ii
     *
     * therefore the exact solution is
     *
     *       x_i = 2
     *
     * for every component.
     */

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);

    if( err )
    {
        fprintf(stderr,
                "lis_matrix_create failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    if( A->nprocs>global_n )
        global_n = A->nprocs;


    err = lis_matrix_set_size(
        A,
        0,
        global_n);

    if( err )
    {
        fprintf(stderr,
                "lis_matrix_set_size failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    err = lis_matrix_get_range(
        A,
        &is,
        &ie);

    if( err )
    {
        fprintf(stderr,
                "lis_matrix_get_range failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    ctx.n  = ie-is;
    ctx.is = is;

    ctx.even_diag = (LIS_SCALAR)4.0;
    ctx.odd_diag  = (LIS_SCALAR)9.0;

    ctx.capture_input = 0;
    ctx.captured_input = 0;
    ctx.captured_x0 = (LIS_SCALAR)0.0;


    err = lis_matrix_set_user(
        A,
        &ctx,
        user_matvec,
        user_matvech);

    if( err )
    {
        fprintf(stderr,
                "lis_matrix_set_user failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    err = lis_matrix_set_user_diagonal(
        A,
        user_get_diagonal);

    if( err )
    {
        fprintf(stderr,
                "lis_matrix_set_user_diagonal failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_duplicate(
        A,
        &b);

    if( err )
    {
        fprintf(stderr,
                "duplicate b failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_duplicate(
        A,
        &x);

    if( err )
    {
        fprintf(stderr,
                "duplicate x failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    for(i=0;i<ctx.n;i++)
    {
        b->value[i] =
            (LIS_SCALAR)2.0*
            diag_value(&ctx,i);

        x->value[i] =
            (LIS_SCALAR)0.0;
    }


    err = lis_solver_create(
        &solver);

    if( err )
    {
        fprintf(stderr,
                "lis_solver_create failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    err = lis_solver_set_option(
        "-i cg "
        "-p none "
        "-scale symm_diag "
        "-initx_zeros false "
        "-tol 1.0e-12 "
        "-maxiter 100 "
        "-print none",
        solver);

    if( err )
    {
        fprintf(stderr,
                "lis_solver_set_option failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    /*
     * Solve #1:
     *
     *   diag(A) = 4, 9, 4, 9, ...
     *   D       = 1/2, 1/3, 1/2, 1/3, ...
     *   b       = 2*diag(A)
     *   x0      = 0
     *   x       = 2
     */
    err = lis_solve(
        A,
        b,
        x,
        solver);

    if( err )
    {
        fprintf(
            stderr,
            "solve #1 API error: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    if( solver->retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "solve #1 retcode: %d\n",
            (int)solver->retcode);

        failed = 1;
        goto cleanup;
    }

    failed |= check_constant_vector(
        "solve #1 solution",
        x,
        is,
        ctx.n,
        (LIS_SCALAR)2.0);

    failed |= check_scale_vector(
        &ctx,
        solver->d,
        (LIS_SCALAR)0.5,
        (LIS_SCALAR)(1.0/3.0));

    failed |= check_rhs(
        &ctx,
        b,
        (LIS_SCALAR)2.0);

    if( A->is_scaled )
    {
        fprintf(
            stderr,
            "solve #1 modified A->is_scaled\n");

        failed = 1;
    }

    if( solver->A!=A || solver->b!=b )
    {
        fprintf(
            stderr,
            "solve #1 did not restore original A/b\n");

        failed = 1;
    }


    /*
     * Solve #2:
     *
     * Change values of the SAME USER matrix.
     *
     *   diag(A) = 16, 25, 16, 25, ...
     *   D       = 1/4, 1/5, 1/4, 1/5, ...
     *   b       = 3*diag(A)
     *   x0      = 0
     *   x       = 3
     *
     * Checking solver->d directly proves that scaling was
     * recomputed instead of reusing stale values from solve #1.
     */
    ctx.even_diag = (LIS_SCALAR)16.0;
    ctx.odd_diag  = (LIS_SCALAR)25.0;

    for(i=0;i<ctx.n;i++)
    {
        b->value[i] =
            (LIS_SCALAR)3.0*
            diag_value(&ctx,i);

        x->value[i] =
            (LIS_SCALAR)0.0;
    }

    err = lis_solve(
        A,
        b,
        x,
        solver);

    if( err )
    {
        fprintf(
            stderr,
            "solve #2 API error: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    if( solver->retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "solve #2 retcode: %d\n",
            (int)solver->retcode);

        failed = 1;
        goto cleanup;
    }

    failed |= check_constant_vector(
        "solve #2 solution",
        x,
        is,
        ctx.n,
        (LIS_SCALAR)3.0);

    failed |= check_scale_vector(
        &ctx,
        solver->d,
        (LIS_SCALAR)0.25,
        (LIS_SCALAR)0.20);

    failed |= check_rhs(
        &ctx,
        b,
        (LIS_SCALAR)3.0);

    if( A->is_scaled )
    {
        fprintf(
            stderr,
            "solve #2 modified A->is_scaled\n");

        failed = 1;
    }

    if( solver->A!=A || solver->b!=b )
    {
        fprintf(
            stderr,
            "solve #2 did not restore original A/b\n");

        failed = 1;
    }


    /*
     * Solve #3: physical warm start.
     *
     * Exact physical solution:
     *
     *     x = 5
     *
     * Caller gives:
     *
     *     x0 = 7
     *
     * Scaled unknown:
     *
     *     x = D*y
     *     y0 = D^-1*x0
     *
     * D*A*D therefore passes D*y0 = x0 to the underlying
     * application matvec.  Capture the FIRST physical USER
     * matvec input and require it to be 7.
     */
    for(i=0;i<ctx.n;i++)
    {
        b->value[i] =
            (LIS_SCALAR)5.0*
            diag_value(&ctx,i);

        x->value[i] =
            (LIS_SCALAR)7.0;
    }

    ctx.capture_input = 1;
    ctx.captured_input = 0;
    ctx.captured_x0 = (LIS_SCALAR)0.0;

    err = lis_solve(
        A,
        b,
        x,
        solver);

    ctx.capture_input = 0;

    if( err )
    {
        fprintf(
            stderr,
            "solve #3 API error: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    if( solver->retcode!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "solve #3 retcode: %d\n",
            (int)solver->retcode);

        failed = 1;
        goto cleanup;
    }

    if( !ctx.captured_input )
    {
        fprintf(
            stderr,
            "solve #3 did not capture USER matvec input\n");

        failed = 1;
    }
    else if( !scalar_close(
                ctx.captured_x0,
                (LIS_SCALAR)7.0) )
    {
        print_value(
            "warm-start physical USER input",
            is,
            ctx.captured_x0,
            (LIS_SCALAR)7.0);

        failed = 1;
    }

    failed |= check_constant_vector(
        "solve #3 solution",
        x,
        is,
        ctx.n,
        (LIS_SCALAR)5.0);

    failed |= check_scale_vector(
        &ctx,
        solver->d,
        (LIS_SCALAR)0.25,
        (LIS_SCALAR)0.20);

    failed |= check_rhs(
        &ctx,
        b,
        (LIS_SCALAR)5.0);

    if( A->is_scaled )
    {
        fprintf(
            stderr,
            "solve #3 modified A->is_scaled\n");

        failed = 1;
    }

    if( solver->A!=A || solver->b!=b )
    {
        fprintf(
            stderr,
            "solve #3 did not restore original A/b\n");

        failed = 1;
    }


    if( !failed )
    {
        printf(
            "LIS_SHELL_SCALE_STAGE6A PASSED\n");
    }


cleanup:

    if( solver )
        lis_solver_destroy(solver);

    if( x )
        lis_vector_destroy(x);

    if( b )
        lis_vector_destroy(b);

    if( A )
        lis_matrix_destroy(A);

    lis_finalize();

    return failed ? 1 : 0;
}
