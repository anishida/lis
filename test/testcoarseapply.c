#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#include <stdio.h>

#include "lislib.h"

#define TEST_N 32
#define TEST_EPS 1.0e-12


static LIS_REAL
test_scalar_real(
    LIS_SCALAR value)
{
#ifdef _COMPLEX
    return creal(value);
#else
    return value;
#endif
}


static LIS_INT
build_matrix(
    LIS_MATRIX *Aout)
{
    LIS_MATRIX A = NULL;
    LIS_INT i,err;
    LIS_SCALAR d;

    err = lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) return err;

    err = lis_matrix_set_size(A,0,TEST_N);
    if(err) goto fail;

    for(i=0;i<TEST_N;i++)
    {
        d = (i==0)
            ? TEST_EPS
            : 1.0 + 0.05*(LIS_REAL)i;

        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                d,
                A);

        if(err) goto fail;
    }

    err = lis_matrix_assemble(A);
    if(err) goto fail;

    *Aout = A;
    return LIS_SUCCESS;

fail:
    lis_matrix_destroy(A);
    return err;
}


static LIS_REAL
true_residual(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_VECTOR w)
{
    LIS_REAL rn = 0.0;
    LIS_REAL bn = 0.0;

    lis_matvec(A,x,w);
    lis_vector_axpy(-1.0,b,w);
    lis_vector_nrm2(w,&rn);
    lis_vector_nrm2(b,&bn);

    if(bn!=0.0)
        return rn/bn;

    return rn;
}


static int
check_apply(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR z = NULL;
    LIS_VECTOR Az = NULL;
    LIS_VECTOR out = NULL;
    LIS_VECTOR diff = NULL;
    LIS_VECTOR basis[1];

    LIS_SOLVER solver = NULL;
    LIS_PRECON precon = NULL;

    LIS_REAL error = -1.0;
    LIS_INT err;
    int failed = 0;

    err = build_matrix(&A);
    if(err) return 1;

    if(
        lis_vector_duplicate(A,&z)
        ||
        lis_vector_duplicate(A,&Az)
        ||
        lis_vector_duplicate(A,&out)
        ||
        lis_vector_duplicate(A,&diff))
    {
        failed = 1;
        goto cleanup;
    }

    lis_vector_set_all(0.0,z);
    z->value[0] = 1.0;

    err = lis_matvec(A,z,Az);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_create(&solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    err =
        lis_solver_set_option(
            "-i fgmres -p none "
            "-scale none -print none "
            "-restart 8 "
            "-maxiter 100 "
            "-maxiter_noimp 0 "
            "-tol 1.0e-12",
            solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    basis[0] = z;

    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    solver->A = A;
    solver->precision = LIS_PRECISION_DOUBLE;

    err =
        lis_precon_create(
            solver,
            &precon);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    solver->precon = precon;

    err =
        lis_solver_near_nullspace_coarse_setup(
            solver,
            A,
            LIS_SCALE_NONE);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    if(
        solver->near_nullspace_coarse==NULL
        ||
        !solver->near_nullspace_coarse->ready
        ||
        solver->near_nullspace_coarse->work==NULL
        ||
        solver->near_nullspace_coarse->rhs==NULL
        ||
        solver->near_nullspace_coarse->coeff==NULL)
    {
        printf("APPLY_STATE invalid\n");
        failed = 1;
        goto cleanup;
    }

    err =
        lis_psolve(
            solver,
            Az,
            out);

    if(err)
    {
        printf(
            "APPLY rc=%d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    lis_vector_copy(out,diff);
    lis_vector_axpy(-1.0,z,diff);
    lis_vector_nrm2(diff,&error);

    printf(
        "APPLY_BAZ "
        "error=% .12e "
        "E=% .12e\n",
        (double)error,
        (double)test_scalar_real(
            solver->near_nullspace_coarse->E[0]));

    if(
        !isfinite((double)error)
        ||
        error>1.0e-8)
    {
        failed = 1;
    }

cleanup:

    if(solver)
        solver->precon = NULL;

    if(precon)
        lis_precon_destroy(precon);

    if(solver)
        lis_solver_destroy(solver);

    if(diff)
        lis_vector_destroy(diff);

    if(out)
        lis_vector_destroy(out);

    if(Az)
        lis_vector_destroy(Az);

    if(z)
        lis_vector_destroy(z);

    if(A)
        lis_matrix_destroy(A);

    return failed;
}


static int
check_fgmres(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR z = NULL;
    LIS_VECTOR w = NULL;
    LIS_VECTOR basis[1];

    LIS_SOLVER solver = NULL;

    LIS_INT i,err,api;
    LIS_INT status = -999;
    LIS_INT iter = -999;

    LIS_REAL reported = -1.0;
    LIS_REAL tr = -1.0;
    LIS_REAL mode_error = -1.0;

    int failed = 0;

    err = build_matrix(&A);
    if(err) return 1;

    if(
        lis_vector_duplicate(A,&b)
        ||
        lis_vector_duplicate(A,&x)
        ||
        lis_vector_duplicate(A,&z)
        ||
        lis_vector_duplicate(A,&w))
    {
        failed = 1;
        goto cleanup;
    }

    lis_vector_set_all(0.0,z);
    z->value[0] = 1.0;

    b->value[0] = 1.0;

    for(i=1;i<TEST_N;i++)
    {
        b->value[i] =
            sin(0.31*(LIS_REAL)(i+1))
            +
            0.20*cos(0.17*(LIS_REAL)(i+1));
    }

    lis_vector_set_all(0.0,x);

    err = lis_solver_create(&solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    err =
        lis_solver_set_option(
            "-i fgmres -p none "
            "-scale none -print none "
            "-restart 8 "
            "-maxiter 200 "
            "-maxiter_noimp 0 "
            "-tol 1.0e-10",
            solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    basis[0] = z;

    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    api =
        lis_solve(
            A,
            b,
            x,
            solver);

    lis_solver_get_status(
        solver,
        &status);

    lis_solver_get_iter(
        solver,
        &iter);

    lis_solver_get_residualnorm(
        solver,
        &reported);

    tr =
        true_residual(
            A,
            b,
            x,
            w);

    mode_error =
        fabs(
            test_scalar_real(
                TEST_EPS*x->value[0]
                -
                1.0));

    printf(
        "FGMRES_TWO_LEVEL "
        "api=%d status=%d iter=%d "
        "reported=% .12e "
        "true=% .12e "
        "mode_error=% .12e "
        "ready=%d\n",
        (int)api,
        (int)status,
        (int)iter,
        (double)reported,
        (double)tr,
        (double)mode_error,
        solver->near_nullspace_coarse
            ? (int)solver->near_nullspace_coarse->ready
            : 0);

    if(api!=LIS_SUCCESS)
        failed = 1;

    if(status!=LIS_SUCCESS)
        failed = 1;

    if(
        !isfinite((double)reported)
        ||
        !isfinite((double)tr)
        ||
        !isfinite((double)mode_error))
    {
        failed = 1;
    }

    if(tr>1.0e-8)
        failed = 1;

    if(mode_error>1.0e-8)
        failed = 1;

    if(
        solver->near_nullspace_coarse==NULL
        ||
        !solver->near_nullspace_coarse->ready)
    {
        failed = 1;
    }

cleanup:

    if(solver)
        solver->precon = NULL;

    if(solver)
        lis_solver_destroy(solver);

    if(w)
        lis_vector_destroy(w);

    if(z)
        lis_vector_destroy(z);

    if(x)
        lis_vector_destroy(x);

    if(b)
        lis_vector_destroy(b);

    if(A)
        lis_matrix_destroy(A);

    return failed;
}



/* ============================================================ */
/* Stage 7C hardening tests                                     */
/* ============================================================ */

static int
test_scalar_near(
    LIS_SCALAR actual,
    LIS_SCALAR expected,
    LIS_REAL tol)
{
#ifdef _COMPLEX
    return
        cabs(actual-expected)
        <=
        tol;
#else
    return
        fabs(actual-expected)
        <=
        tol;
#endif
}


static LIS_INT
build_matrix_first(
    LIS_SCALAR first,
    LIS_MATRIX *Aout)
{
    LIS_MATRIX A = NULL;

    LIS_INT i,is,ie,err;

    LIS_SCALAR d;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if(err)
        return err;


    err =
        lis_matrix_set_size(
            A,
            0,
            TEST_N);

    if(err)
        goto fail;


    err =
        lis_matrix_get_range(
            A,
            &is,
            &ie);

    if(err)
        goto fail;


    for(i=is;i<ie;i++)
    {
        d =
            (i==0)
            ?
            first
            :
            (LIS_SCALAR)(
                1.0
                +
                0.05*(LIS_REAL)i);


        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                d,
                A);

        if(err)
            goto fail;
    }


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if(err)
        goto fail;


    err =
        lis_matrix_assemble(
            A);

    if(err)
        goto fail;


    *Aout =
        A;


    return LIS_SUCCESS;


fail:

    if(A)
        lis_matrix_destroy(A);


    return err;
}


static LIS_INT
build_pivot_matrix(
    LIS_MATRIX *Aout)
{
    LIS_MATRIX A = NULL;

    LIS_INT i,is,ie,err;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if(err)
        return err;


    err =
        lis_matrix_set_size(
            A,
            0,
            TEST_N);

    if(err)
        goto fail;


    err =
        lis_matrix_get_range(
            A,
            &is,
            &ie);

    if(err)
        goto fail;


    for(i=is;i<ie;i++)
    {
        if(i==0)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    0,
                    1,
                    (LIS_SCALAR)1.0,
                    A);

            if(err)
                goto fail;
        }
        else if(i==1)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    1,
                    0,
                    (LIS_SCALAR)1.0,
                    A);

            if(err)
                goto fail;


            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    1,
                    1,
                    (LIS_SCALAR)1.0,
                    A);

            if(err)
                goto fail;
        }
        else
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i,
                    (LIS_SCALAR)(i+3),
                    A);

            if(err)
                goto fail;
        }
    }


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if(err)
        goto fail;


    err =
        lis_matrix_assemble(
            A);

    if(err)
        goto fail;


    *Aout =
        A;


    return LIS_SUCCESS;


fail:

    if(A)
        lis_matrix_destroy(A);


    return err;
}


static int
check_pivot_coarse(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR z0 = NULL;
    LIS_VECTOR z1 = NULL;

    LIS_VECTOR basis[2];

    LIS_SOLVER solver = NULL;

    LIS_NEAR_NULLSPACE_COARSE coarse = NULL;

    LIS_SCALAR rhs[2];
    LIS_SCALAR coeff[2];

    LIS_INT err = -999;
    LIS_INT solve_rc = -999;

    int failed = 0;


    err =
        build_pivot_matrix(
            &A);

    if(err)
        return 1;


    if(
        lis_vector_duplicate(A,&z0)
        ||
        lis_vector_duplicate(A,&z1))
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        z0);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        z1);


    if(z0->n>0)
        z0->value[0] =
            (LIS_SCALAR)1.0;

    if(z1->n>1)
        z1->value[1] =
            (LIS_SCALAR)1.0;


    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    basis[0] = z0;
    basis[1] = z1;


    err =
        lis_solver_set_near_nullspace(
            solver,
            2,
            basis);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    solver->precision =
        LIS_PRECISION_DOUBLE;


    err =
        lis_solver_near_nullspace_coarse_setup(
            solver,
            A,
            LIS_SCALE_NONE);

    if(err)
    {
        failed = 1;
        goto report;
    }


    coarse =
        solver->near_nullspace_coarse;


    if(
        coarse==NULL
        ||
        !coarse->ready
        ||
        coarse->dim!=2)
    {
        failed = 1;
        goto report;
    }


    /*
     * E =
     *
     *     [ 0  1 ]
     *     [ 1  1 ]
     *
     * Hence the first LU pivot is zero before pivoting.
     * Solving
     *
     *     E [2,3]^T = [3,5]^T
     *
     * therefore exercises the row-swap path.
     */
    rhs[0] =
        (LIS_SCALAR)3.0;

    rhs[1] =
        (LIS_SCALAR)5.0;


    solve_rc =
        lis_solver_near_nullspace_coarse_solve(
            solver,
            rhs,
            coeff);


    if(
        solve_rc
        ||
        !test_scalar_near(
            coarse->E[0],
            (LIS_SCALAR)0.0,
            1.0e-12)
        ||
        !test_scalar_near(
            coarse->E[1],
            (LIS_SCALAR)1.0,
            1.0e-12)
        ||
        !test_scalar_near(
            coarse->E[2],
            (LIS_SCALAR)1.0,
            1.0e-12)
        ||
        !test_scalar_near(
            coarse->E[3],
            (LIS_SCALAR)1.0,
            1.0e-12)
        ||
        !test_scalar_near(
            coeff[0],
            (LIS_SCALAR)2.0,
            1.0e-11)
        ||
        !test_scalar_near(
            coeff[1],
            (LIS_SCALAR)3.0,
            1.0e-11))
    {
        failed = 1;
    }


report:

    printf(
        "PIVOT_COARSE "
        "setup=%d solve=%d "
        "E00=% .12e E01=% .12e "
        "E10=% .12e E11=% .12e "
        "c0=% .12e c1=% .12e\n",
        (int)err,
        (int)solve_rc,
        coarse
            ?
            (double)test_scalar_real(coarse->E[0])
            :
            -999.0,
        coarse
            ?
            (double)test_scalar_real(coarse->E[1])
            :
            -999.0,
        coarse
            ?
            (double)test_scalar_real(coarse->E[2])
            :
            -999.0,
        coarse
            ?
            (double)test_scalar_real(coarse->E[3])
            :
            -999.0,
        coarse
            ?
            (double)test_scalar_real(coeff[0])
            :
            -999.0,
        coarse
            ?
            (double)test_scalar_real(coeff[1])
            :
            -999.0);


cleanup:

    if(solver)
        lis_solver_destroy(
            solver);

    if(z1)
        lis_vector_destroy(
            z1);

    if(z0)
        lis_vector_destroy(
            z0);

    if(A)
        lis_matrix_destroy(
            A);


    return failed;
}


static int
check_changing_matrix_rebuild(void)
{
    const LIS_SCALAR first0 =
        (LIS_SCALAR)TEST_EPS;

    const LIS_SCALAR first1 =
        (LIS_SCALAR)(7.0*TEST_EPS);


    LIS_MATRIX A0 = NULL;
    LIS_MATRIX A1 = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR z = NULL;
    LIS_VECTOR exact = NULL;
    LIS_VECTOR w = NULL;

    LIS_VECTOR basis[1];

    LIS_SOLVER solver = NULL;

    LIS_SCALAR E0 =
        (LIS_SCALAR)0.0;

    LIS_SCALAR E1 =
        (LIS_SCALAR)0.0;

    LIS_REAL tr0 = -1.0;
    LIS_REAL tr1 = -1.0;

    LIS_INT api0 = -999;
    LIS_INT api1 = -999;
    LIS_INT err;

    int failed = 0;


    err =
        build_matrix_first(
            first0,
            &A0);

    if(err)
        return 1;


    err =
        build_matrix_first(
            first1,
            &A1);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    if(
        lis_vector_duplicate(A0,&b)
        ||
        lis_vector_duplicate(A0,&x)
        ||
        lis_vector_duplicate(A0,&z)
        ||
        lis_vector_duplicate(A0,&exact)
        ||
        lis_vector_duplicate(A0,&w))
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        z);

    z->value[0] =
        (LIS_SCALAR)1.0;


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        exact);


    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err =
        lis_solver_set_option(
            "-i fgmres -p none "
            "-scale none -print none "
            "-restart 8 "
            "-maxiter 200 "
            "-maxiter_noimp 0 "
            "-tol 1.0e-10",
            solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err =
        lis_matvec(
            A0,
            exact,
            b);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    api0 =
        lis_solve(
            A0,
            b,
            x,
            solver);


    if(
        api0==LIS_SUCCESS
        &&
        solver->near_nullspace_coarse
        &&
        solver->near_nullspace_coarse->ready)
    {
        E0 =
            solver
            ->
            near_nullspace_coarse
            ->
            E[0];

        tr0 =
            true_residual(
                A0,
                b,
                x,
                w);
    }
    else
    {
        failed = 1;
    }


    err =
        lis_matvec(
            A1,
            exact,
            b);

    if(err)
    {
        failed = 1;
        goto report;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    api1 =
        lis_solve(
            A1,
            b,
            x,
            solver);


    if(
        api1==LIS_SUCCESS
        &&
        solver->near_nullspace_coarse
        &&
        solver->near_nullspace_coarse->ready)
    {
        E1 =
            solver
            ->
            near_nullspace_coarse
            ->
            E[0];

        tr1 =
            true_residual(
                A1,
                b,
                x,
                w);
    }
    else
    {
        failed = 1;
    }


    if(
        !test_scalar_near(
            E0,
            first0,
            1.0e-18)
        ||
        !test_scalar_near(
            E1,
            first1,
            1.0e-18)
        ||
        !isfinite((double)tr0)
        ||
        !isfinite((double)tr1)
        ||
        tr0>1.0e-8
        ||
        tr1>1.0e-8)
    {
        failed = 1;
    }


report:

    printf(
        "REBUILD_A "
        "api0=%d api1=%d "
        "E0=% .12e E1=% .12e "
        "true0=% .12e true1=% .12e "
        "ready=%d\n",
        (int)api0,
        (int)api1,
        (double)test_scalar_real(E0),
        (double)test_scalar_real(E1),
        (double)tr0,
        (double)tr1,
        solver
            &&
        solver->near_nullspace_coarse
            ?
            (int)
                solver
                ->
                near_nullspace_coarse
                ->
                ready
            :
            0);


cleanup:

    if(solver)
        solver->precon = NULL;

    if(solver)
        lis_solver_destroy(
            solver);

    if(w)
        lis_vector_destroy(
            w);

    if(exact)
        lis_vector_destroy(
            exact);

    if(z)
        lis_vector_destroy(
            z);

    if(x)
        lis_vector_destroy(
            x);

    if(b)
        lis_vector_destroy(
            b);

    if(A1)
        lis_matrix_destroy(
            A1);

    if(A0)
        lis_matrix_destroy(
            A0);


    return failed;
}


/* ------------------------------------------------------------ */
/* USER preconditioner below the two-level wrapper               */
/* ------------------------------------------------------------ */

static int
stage7c_user_magic =
    0x37435550;

static int
stage7c_user_create_count =
    0;

static int
stage7c_user_psolve_count =
    0;


static LIS_INT
stage7c_user_create(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    (void)solver;

    stage7c_user_create_count++;


    return
        lis_precon_set_user_data(
            precon,
            &stage7c_user_magic);
}


static LIS_INT
stage7c_user_psolve(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    void *user_data = NULL;

    LIS_INT err;


    err =
        lis_precon_get_user_data(
            solver->precon,
            &user_data);

    if(err)
        return err;


    if(
        user_data
        !=
        &stage7c_user_magic)
    {
        return LIS_FAILS;
    }


    stage7c_user_psolve_count++;


    return
        lis_vector_copy(
            b,
            x);
}


static LIS_INT
stage7c_user_psolveh(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    return
        stage7c_user_psolve(
            solver,
            b,
            x);
}


static int
check_user_wrapper(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR z = NULL;
    LIS_VECTOR exact = NULL;
    LIS_VECTOR w = NULL;

    LIS_VECTOR basis[1];

    LIS_SOLVER solver = NULL;

    LIS_REAL tr = -1.0;

    LIS_INT api = -999;
    LIS_INT err;

    int registered = 0;
    int failed = 0;


    stage7c_user_create_count =
        0;

    stage7c_user_psolve_count =
        0;


    err =
        lis_precon_register(
            "c7u",
            stage7c_user_create,
            stage7c_user_psolve,
            stage7c_user_psolveh);

    if(err)
        return 1;


    registered =
        1;


    err =
        build_matrix(
            &A);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    if(
        lis_vector_duplicate(A,&b)
        ||
        lis_vector_duplicate(A,&x)
        ||
        lis_vector_duplicate(A,&z)
        ||
        lis_vector_duplicate(A,&exact)
        ||
        lis_vector_duplicate(A,&w))
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        z);

    z->value[0] =
        (LIS_SCALAR)1.0;


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        exact);


    err =
        lis_matvec(
            A,
            exact,
            b);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err =
        lis_solver_set_option(
            "-i fgmres -p c7u "
            "-scale none -print none "
            "-restart 8 "
            "-maxiter 200 "
            "-maxiter_noimp 0 "
            "-tol 1.0e-10",
            solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    api =
        lis_solve(
            A,
            b,
            x,
            solver);


    if(api==LIS_SUCCESS)
    {
        tr =
            true_residual(
                A,
                b,
                x,
                w);
    }


    printf(
        "USER_WRAPPER "
        "api=%d create=%d psolve=%d "
        "true=% .12e ready=%d\n",
        (int)api,
        stage7c_user_create_count,
        stage7c_user_psolve_count,
        (double)tr,
        solver
            &&
        solver->near_nullspace_coarse
            ?
            (int)
                solver
                ->
                near_nullspace_coarse
                ->
                ready
            :
            0);


    if(
        api!=LIS_SUCCESS
        ||
        stage7c_user_create_count!=1
        ||
        stage7c_user_psolve_count<1
        ||
        !isfinite((double)tr)
        ||
        tr>1.0e-8
        ||
        solver->near_nullspace_coarse==NULL
        ||
        !solver->near_nullspace_coarse->ready)
    {
        failed = 1;
    }


cleanup:

    if(solver)
        solver->precon = NULL;

    if(solver)
        lis_solver_destroy(
            solver);

    if(w)
        lis_vector_destroy(
            w);

    if(exact)
        lis_vector_destroy(
            exact);

    if(z)
        lis_vector_destroy(
            z);

    if(x)
        lis_vector_destroy(
            x);

    if(b)
        lis_vector_destroy(
            b);

    if(A)
        lis_matrix_destroy(
            A);


    if(registered)
        lis_precon_register_free();


    return failed;
}


static int
check_solver_gate(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR z = NULL;

    LIS_VECTOR basis[1];

    LIS_SOLVER solver = NULL;

    LIS_INT cg_rc = -999;
    LIS_INT invalid_rc = -999;
    LIS_INT err;

    int failed = 0;


    err =
        build_matrix(
            &A);

    if(err)
        return 1;


    if(
        lis_vector_duplicate(A,&b)
        ||
        lis_vector_duplicate(A,&x)
        ||
        lis_vector_duplicate(A,&z))
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        z);

    z->value[0] =
        (LIS_SCALAR)1.0;


    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err =
        lis_solver_set_option(
            "-i cg -p none "
            "-scale none -print none "
            "-maxiter 20 "
            "-tol 1.0e-10",
            solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    /*
     * Valid but currently unsupported Stage 7C solver.
     */
    cg_rc =
        lis_solve(
            A,
            b,
            x,
            solver);


    if(
        cg_rc!=LIS_ERR_NOT_IMPLEMENTED
        ||
        solver->near_nullspace_coarse!=NULL)
    {
        failed = 1;
    }


    /*
     * Invalid generic option must retain the normal LIS error.
     */
    solver->options[LIS_OPTIONS_SOLVER] =
        LIS_SOLVER_LEN + 1;


    invalid_rc =
        lis_solve(
            A,
            b,
            x,
            solver);


    if(
        invalid_rc!=LIS_ERR_ILL_ARG
        ||
        solver->near_nullspace_coarse!=NULL)
    {
        failed = 1;
    }


    printf(
        "SOLVER_GATE "
        "cg=%d expected_cg=%d "
        "invalid=%d expected_invalid=%d "
        "coarse=%p\n",
        (int)cg_rc,
        (int)LIS_ERR_NOT_IMPLEMENTED,
        (int)invalid_rc,
        (int)LIS_ERR_ILL_ARG,
        (void *)solver->near_nullspace_coarse);


cleanup:

    if(solver)
        solver->precon = NULL;

    if(solver)
        lis_solver_destroy(
            solver);

    if(z)
        lis_vector_destroy(
            z);

    if(x)
        lis_vector_destroy(
            x);

    if(b)
        lis_vector_destroy(
            b);

    if(A)
        lis_matrix_destroy(
            A);


    return failed;
}


int
main(
    int argc,
    char **argv)
{
    LIS_INT err;
    int failed = 0;

    err =
        lis_initialize(
            &argc,
            &argv);

    if(err)
        return 1;

    failed |=
        check_apply();

    failed |=
        check_fgmres();

    failed |=
        check_pivot_coarse();

    failed |=
        check_changing_matrix_rebuild();

    failed |=
        check_solver_gate();

    failed |=
        check_user_wrapper();

    if(failed)
    {
        printf(
            "LIS_COARSE_APPLY_STAGE7C FAILED\n");

        lis_finalize();

        return 1;
    }

    printf(
        "LIS_COARSE_APPLY_STAGE7C PASSED\n");

    lis_finalize();

    return 0;
}
