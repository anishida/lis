#include "lis.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static LIS_INT
make_system(
    LIS_MATRIX *A,
    LIS_VECTOR *b,
    LIS_VECTOR *q)
{
    LIS_INT err;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        A);

    if(err)
    {
        return err;
    }

    err = lis_matrix_set_size(
        *A,
        0,
        3);

    if(err)
    {
        return err;
    }

    lis_matrix_set_value(
        LIS_INS_VALUE,
        0,0, 1.0,
        *A);

    lis_matrix_set_value(
        LIS_INS_VALUE,
        0,1,-1.0,
        *A);

    lis_matrix_set_value(
        LIS_INS_VALUE,
        1,0,-1.0,
        *A);

    lis_matrix_set_value(
        LIS_INS_VALUE,
        1,1, 2.0,
        *A);

    lis_matrix_set_value(
        LIS_INS_VALUE,
        1,2,-1.0,
        *A);

    lis_matrix_set_value(
        LIS_INS_VALUE,
        2,1,-1.0,
        *A);

    lis_matrix_set_value(
        LIS_INS_VALUE,
        2,2, 1.0,
        *A);

    err = lis_matrix_assemble(
        *A);

    if(err)
    {
        return err;
    }

    err = lis_vector_duplicate(
        *A,
        b);

    if(err)
    {
        return err;
    }

    err = lis_vector_duplicate(
        *A,
        q);

    if(err)
    {
        return err;
    }

    /*
     * Compatible RHS:
     *
     * b = [1,-2,1]^T
     *
     * ker(A) = span([1,1,1]^T)
     *
     * gauge-fixed solution:
     *
     * x = [1/3,-2/3,1/3]^T
     */
    lis_vector_set_all(
        0.0,
        *b);

    lis_vector_set_value(
        LIS_INS_VALUE,
        0,
        1.0,
        *b);

    lis_vector_set_value(
        LIS_INS_VALUE,
        1,
        -2.0,
        *b);

    lis_vector_set_value(
        LIS_INS_VALUE,
        2,
        1.0,
        *b);

    lis_vector_set_all(
        1.0,
        *q);

    return LIS_SUCCESS;
}


static int
run_supported(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR q,
    const char *name)
{
    LIS_SOLVER solver = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR basis[1];

    LIS_INT err;
    LIS_INT api;

    double x0,x1,x2;
    double r0,r1,r2;
    double rel;
    double gauge;
    double solerr;

    char opts[512];

    int pass = 1;

    err = lis_vector_duplicate(
        A,
        &x);

    if(err)
    {
        return 0;
    }

    lis_vector_set_all(
        0.0,
        x);

    err = lis_solver_create(
        &solver);

    if(err)
    {
        lis_vector_destroy(x);
        return 0;
    }

    if(strcmp(name,"fgcrodr")==0)
    {
        snprintf(
            opts,
            sizeof(opts),
            "-i fgcrodr "
            "-p jacobi "
            "-scale none "
            "-tol 1.0e-12 "
            "-maxiter 100 "
            "-restart 6 "
            "-fgcrodr_k 1");
    }
    else
    {
        snprintf(
            opts,
            sizeof(opts),
            "-i fgmres "
            "-p jacobi "
            "-scale none "
            "-tol 1.0e-12 "
            "-maxiter 100 "
            "-restart 6");
    }

    err = lis_solver_set_option(
        opts,
        solver);

    if(err)
    {
        pass = 0;
    }

    basis[0] = q;

    err = lis_solver_set_nullspace(
        solver,
        1,
        basis);

    if(err)
    {
        pass = 0;
    }

    api = lis_solve(
        A,
        b,
        x,
        solver);

    x0 = (double)x->value[0];
    x1 = (double)x->value[1];
    x2 = (double)x->value[2];

    r0 =
        x0
        - x1
        - 1.0;

    r1 =
        -x0
        + 2.0*x1
        - x2
        + 2.0;

    r2 =
        -x1
        + x2
        - 1.0;

    rel =
        sqrt(
            r0*r0
            +
            r1*r1
            +
            r2*r2)
        /
        sqrt(6.0);

    gauge =
        fabs(
            x0
            +
            x1
            +
            x2);

    solerr =
        sqrt(
            (x0-1.0/3.0)*(x0-1.0/3.0)
            +
            (x1+2.0/3.0)*(x1+2.0/3.0)
            +
            (x2-1.0/3.0)*(x2-1.0/3.0));

    printf(
        "NULLSPACE_SOLVE solver=%s "
        "api=%d status=%d iter=%d "
        "x=[% .17e,% .17e,% .17e] "
        "relres=% .17e "
        "gauge=% .17e "
        "solerr=% .17e\n",
        name,
        (int)api,
        (int)solver->retcode,
        (int)solver->iter,
        x0,x1,x2,
        rel,
        gauge,
        solerr);

    if(
        api!=LIS_SUCCESS
        ||
        solver->retcode!=LIS_SUCCESS
        ||
        !isfinite(rel)
        ||
        !isfinite(gauge)
        ||
        !isfinite(solerr)
        ||
        rel>1.0e-9
        ||
        gauge>1.0e-9
        ||
        solerr>1.0e-8)
    {
        pass = 0;
    }

    lis_solver_destroy(
        solver);

    lis_vector_destroy(
        x);

    return pass;
}


static int
run_reject(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR q,
    const char *opts,
    const char *name)
{
    LIS_SOLVER solver = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR basis[1];

    LIS_INT err;
    LIS_INT api;
    LIS_INT status;

    int pass = 1;

    err = lis_vector_duplicate(
        A,
        &x);

    if(err)
    {
        return 0;
    }

    lis_vector_set_all(
        0.0,
        x);

    err = lis_solver_create(
        &solver);

    if(err)
    {
        lis_vector_destroy(x);
        return 0;
    }

    err = lis_solver_set_option(
        opts,
        solver);

    if(err)
    {
        pass = 0;
    }

    basis[0] = q;

    err = lis_solver_set_nullspace(
        solver,
        1,
        basis);

    if(err)
    {
        pass = 0;
    }

    api = lis_solve(
        A,
        b,
        x,
        solver);

    status = solver->retcode;

    printf(
        "NULLSPACE_REJECT case=%s "
        "api=%d status=%d expected=%d\n",
        name,
        (int)api,
        (int)status,
        (int)LIS_ERR_NOT_IMPLEMENTED);

    if(
        api!=LIS_ERR_NOT_IMPLEMENTED
        &&
        status!=LIS_ERR_NOT_IMPLEMENTED)
    {
        pass = 0;
    }

    lis_solver_destroy(
        solver);

    lis_vector_destroy(
        x);

    return pass;
}



static int
run_incompatible_rhs(
    LIS_MATRIX A,
    LIS_VECTOR q)
{
    LIS_SOLVER solver = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR basis[1];

    LIS_INT err;
    LIS_INT api;
    LIS_INT status;

    double x0,x1,x2;
    double r0,r1,r2;
    double rel;

    int pass = 1;

    err = lis_vector_duplicate(
        A,
        &b);

    if(err)
    {
        return 0;
    }

    err = lis_vector_duplicate(
        A,
        &x);

    if(err)
    {
        lis_vector_destroy(b);
        return 0;
    }

    /*
     * Incompatible RHS for ker(A)=span([1,1,1]^T):
     *
     *     b = [1,0,0]^T
     *
     * Right-nullspace deflation must not project b and must not
     * report convergence for the physical residual.
     */
    lis_vector_set_all(
        0.0,
        b);

    lis_vector_set_value(
        LIS_INS_VALUE,
        0,
        1.0,
        b);

    lis_vector_set_all(
        0.0,
        x);

    err = lis_solver_create(
        &solver);

    if(err)
    {
        lis_vector_destroy(x);
        lis_vector_destroy(b);
        return 0;
    }

    err = lis_solver_set_option(
        "-i fgmres "
        "-p jacobi "
        "-scale none "
        "-tol 1.0e-12 "
        "-maxiter 20 "
        "-restart 3",
        solver);

    if(err)
    {
        pass = 0;
    }

    basis[0] = q;

    err = lis_solver_set_nullspace(
        solver,
        1,
        basis);

    if(err)
    {
        pass = 0;
    }

    api = lis_solve(
        A,
        b,
        x,
        solver);

    status = solver->retcode;

    x0 = (double)x->value[0];
    x1 = (double)x->value[1];
    x2 = (double)x->value[2];

    r0 =
        x0
        - x1
        - 1.0;

    r1 =
        -x0
        + 2.0*x1
        - x2;

    r2 =
        -x1
        + x2;

    rel =
        sqrt(
            r0*r0
            +
            r1*r1
            +
            r2*r2);

    printf(
        "NULLSPACE_INCOMPATIBLE "
        "api=%d status=%d iter=%d "
        "relres=% .17e\n",
        (int)api,
        (int)status,
        (int)solver->iter,
        rel);

    if(
        status==LIS_SUCCESS
        ||
        !isfinite(rel)
        ||
        rel<=1.0e-8)
    {
        pass = 0;
    }

    lis_solver_destroy(
        solver);

    lis_vector_destroy(
        x);

    lis_vector_destroy(
        b);

    return pass;
}

int
main(
    int argc,
    char *argv[])
{
    LIS_MATRIX A = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR q = NULL;

    LIS_INT err;

    int fgmres_ok;
    int fgcrodr_ok;
    int solver_reject_ok;
    int scale_reject_ok;
    int incompatible_ok;

    int pass = 1;

    err = lis_initialize(
        &argc,
        &argv);

    if(err)
    {
        printf(
            "INIT_FAIL err=%d\n",
            (int)err);

        return 2;
    }

    err = make_system(
        &A,
        &b,
        &q);

    if(err)
    {
        printf(
            "SYSTEM_BUILD_FAIL err=%d\n",
            (int)err);

        lis_finalize();

        return 3;
    }

    fgmres_ok =
        run_supported(
            A,
            b,
            q,
            "fgmres");

    fgcrodr_ok =
        run_supported(
            A,
            b,
            q,
            "fgcrodr");

    solver_reject_ok =
        run_reject(
            A,
            b,
            q,
            "-i cg "
            "-p jacobi "
            "-scale none "
            "-tol 1.0e-12 "
            "-maxiter 100",
            "unsupported_solver");

    scale_reject_ok =
        run_reject(
            A,
            b,
            q,
            "-i fgmres "
            "-p jacobi "
            "-scale symm_diag "
            "-tol 1.0e-12 "
            "-maxiter 100 "
            "-restart 6",
            "unsupported_scaling");


    incompatible_ok =
        run_incompatible_rhs(
            A,
            q);

    printf(
        "NULLSPACE_SUMMARY "
        "fgmres=%d "
        "fgcrodr=%d "
        "solver_reject=%d "
        "scale_reject=%d "
        "incompatible=%d\n",
        fgmres_ok,
        fgcrodr_ok,
        solver_reject_ok,
        scale_reject_ok,
        incompatible_ok);

    if(
        !fgmres_ok
        ||
        !fgcrodr_ok
        ||
        !solver_reject_ok
        ||
        !scale_reject_ok
        ||
        !incompatible_ok)
    {
        pass = 0;
    }

    lis_vector_destroy(
        q);

    lis_vector_destroy(
        b);

    lis_matrix_destroy(
        A);

    lis_finalize();

    if(pass)
    {
        printf(
            "LIS_EXACT_NULLSPACE_DEFLATION_REGRESSION PASSED\n");

        return 0;
    }

    printf(
        "LIS_EXACT_NULLSPACE_DEFLATION_REGRESSION FAILED\n");

    return 1;
}
