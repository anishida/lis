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
    else if(strcmp(name,"gcrodr")==0)
    {
        snprintf(
            opts,
            sizeof(opts),
            "-i gcrodr "
            "-p jacobi "
            "-scale none "
            "-tol 1.0e-12 "
            "-maxiter 100 "
            "-restart 6 "
            "-gcrodr_keep 1");
    }
    else if(strcmp(name,"gmres")==0)
    {
        snprintf(
            opts,
            sizeof(opts),
            "-i gmres "
            "-p jacobi "
            "-scale none "
            "-tol 1.0e-12 "
            "-maxiter 100 "
            "-restart 6");
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
    char mutable_opts[512];
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

    snprintf(
        mutable_opts,
        sizeof(mutable_opts),
        "%s",
        opts);


    err = lis_solver_set_option(
        mutable_opts,
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
    LIS_VECTOR q,
    const char *name)
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

    char opts[512];

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

    if(strcmp(name,"fgcrodr")==0)
    {
        snprintf(
            opts,
            sizeof(opts),
            "-i fgcrodr "
            "-p jacobi "
            "-scale none "
            "-tol 1.0e-12 "
            "-maxiter 40 "
            "-restart 6 "
            "-fgcrodr_k 1");
    }
    else if(strcmp(name,"gcrodr")==0)
    {
        snprintf(
            opts,
            sizeof(opts),
            "-i gcrodr "
            "-p jacobi "
            "-scale none "
            "-tol 1.0e-12 "
            "-maxiter 40 "
            "-restart 6 "
            "-gcrodr_keep 1");
    }
    else if(strcmp(name,"gmres")==0)
    {
        snprintf(
            opts,
            sizeof(opts),
            "-i gmres "
            "-p jacobi "
            "-scale none "
            "-tol 1.0e-12 "
            "-maxiter 40 "
            "-restart 6");
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
            "-maxiter 40 "
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
        "solver=%s "
        "api=%d status=%d iter=%d "
        "relres=% .17e\n",
        name,
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

static int
run_gcrodr_repeated(
    LIS_MATRIX A,
    LIS_VECTOR q)
{
    LIS_SOLVER solver = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR basis[1];

    LIS_INT err;
    LIS_INT api1;
    LIS_INT api2;
    LIS_INT status1;
    LIS_INT status2;

    double x0,x1,x2;
    double r0,r1,r2;

    double rel1;
    double rel2;
    double gauge1;
    double gauge2;

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

    err = lis_solver_create(
        &solver);

    if(err)
    {
        lis_vector_destroy(x);
        lis_vector_destroy(b);
        return 0;
    }

    err = lis_solver_set_option(
        "-i gcrodr "
        "-p jacobi "
        "-scale none "
        "-tol 1.0e-12 "
        "-maxiter 40 "
        "-restart 6 "
        "-gcrodr_keep 1",
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

    lis_vector_set_all(
        0.0,
        b);

    lis_vector_set_value(
        LIS_INS_VALUE,
        0,
        1.0,
        b);

    lis_vector_set_value(
        LIS_INS_VALUE,
        1,
        -2.0,
        b);

    lis_vector_set_value(
        LIS_INS_VALUE,
        2,
        1.0,
        b);

    lis_vector_set_all(
        0.0,
        x);

    api1 = lis_solve(
        A,
        b,
        x,
        solver);

    status1 = solver->retcode;

    x0 = (double)x->value[0];
    x1 = (double)x->value[1];
    x2 = (double)x->value[2];

    r0 = x0-x1-1.0;
    r1 = -x0+2.0*x1-x2+2.0;
    r2 = -x1+x2-1.0;

    rel1 =
        sqrt(
            r0*r0
            +
            r1*r1
            +
            r2*r2)
        /
        sqrt(6.0);

    gauge1 =
        fabs(
            x0
            +
            x1
            +
            x2);

    if(
        api1!=LIS_SUCCESS
        ||
        status1!=LIS_SUCCESS
        ||
        !isfinite(rel1)
        ||
        !isfinite(gauge1)
        ||
        rel1>1.0e-9
        ||
        gauge1>1.0e-9)
    {
        pass = 0;
    }

    lis_vector_set_all(
        0.0,
        b);

    lis_vector_set_value(
        LIS_INS_VALUE,
        0,
        -1.0,
        b);

    lis_vector_set_value(
        LIS_INS_VALUE,
        1,
        2.0,
        b);

    lis_vector_set_value(
        LIS_INS_VALUE,
        2,
        -1.0,
        b);

    lis_vector_set_all(
        0.0,
        x);

    api2 = lis_solve(
        A,
        b,
        x,
        solver);

    status2 = solver->retcode;

    x0 = (double)x->value[0];
    x1 = (double)x->value[1];
    x2 = (double)x->value[2];

    r0 = x0-x1+1.0;
    r1 = -x0+2.0*x1-x2-2.0;
    r2 = -x1+x2+1.0;

    rel2 =
        sqrt(
            r0*r0
            +
            r1*r1
            +
            r2*r2)
        /
        sqrt(6.0);

    gauge2 =
        fabs(
            x0
            +
            x1
            +
            x2);

    if(
        api2!=LIS_SUCCESS
        ||
        status2!=LIS_SUCCESS
        ||
        !isfinite(rel2)
        ||
        !isfinite(gauge2)
        ||
        rel2>1.0e-9
        ||
        gauge2>1.0e-9)
    {
        pass = 0;
    }

    printf(
        "NULLSPACE_GCRODR_REPEAT "
        "api1=%d status1=%d "
        "rel1=% .17e gauge1=% .17e "
        "api2=%d status2=%d "
        "rel2=% .17e gauge2=% .17e "
        "pass=%d\n",
        (int)api1,
        (int)status1,
        rel1,
        gauge1,
        (int)api2,
        (int)status2,
        rel2,
        gauge2,
        pass);

    lis_solver_destroy(
        solver);

    lis_vector_destroy(
        x);

    lis_vector_destroy(
        b);

    return pass;
}


static LIS_INT
make_vector3(
    LIS_MATRIX A,
    LIS_SCALAR value0,
    LIS_SCALAR value1,
    LIS_SCALAR value2,
    LIS_VECTOR *vector)
{
    LIS_INT err;

    *vector = NULL;

    err =
        lis_vector_duplicate(
            A,
            vector);

    if(err)
    {
        return err;
    }

    err =
        lis_vector_set_all(
            0.0,
            *vector);

    if(err)
    {
        goto fail;
    }

    err =
        lis_vector_set_value(
            LIS_INS_VALUE,
            0,
            value0,
            *vector);

    if(err)
    {
        goto fail;
    }

    err =
        lis_vector_set_value(
            LIS_INS_VALUE,
            1,
            value1,
            *vector);

    if(err)
    {
        goto fail;
    }

    err =
        lis_vector_set_value(
            LIS_INS_VALUE,
            2,
            value2,
            *vector);

    if(err)
    {
        goto fail;
    }

    return LIS_SUCCESS;

fail:
    lis_vector_destroy(
        *vector);

    *vector = NULL;

    return err;
}


static int
run_left_nullspace_api(
    LIS_MATRIX A)
{
    LIS_SOLVER solver = NULL;

    LIS_VECTOR v0 = NULL;
    LIS_VECTOR v1 = NULL;
    LIS_VECTOR replacement = NULL;
    LIS_VECTOR dependent0 = NULL;
    LIS_VECTOR dependent1 = NULL;

    LIS_VECTOR basis[2];
    LIS_VECTOR one[1];
    LIS_VECTOR dependent[2];

    LIS_SCALAR before[3];
    LIS_SCALAR after[3];
    LIS_SCALAR dot = 0.0;

    LIS_REAL nrm0 = 0.0;
    LIS_REAL nrm1 = 0.0;
    LIS_REAL replacement_nrm = 0.0;
    LIS_REAL dot_abs = 0.0;
    LIS_REAL ownership_delta = 0.0;

    LIS_INT err;
    LIS_INT api_error = LIS_SUCCESS;
    LIS_INT dependent_api = LIS_SUCCESS;

    int i;
    int ownership_ok = 0;
    int orthonormal_ok = 0;
    int replacement_ok = 0;
    int clear_ok = 0;
    int zero_clear_ok = 0;
    int dependent_ok = 0;
    int right_untouched_ok = 0;

    int pass = 1;

    const LIS_REAL tol =
        (LIS_REAL)1.0e-10;

    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            1.0,
            0.0,
            &v0);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            0.0,
            1.0,
            &v1);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    basis[0] = v0;
    basis[1] = v1;

    err =
        lis_solver_set_left_nullspace(
            solver,
            2,
            basis);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    ownership_ok =
        solver->left_nullspace!=NULL
        &&
        solver->left_nullspace_dim==2
        &&
        solver->left_nullspace[0]!=v0
        &&
        solver->left_nullspace[1]!=v1;

    if(ownership_ok)
    {
        err =
            lis_vector_get_values(
                solver->left_nullspace[0],
                0,
                3,
                before);

        if(err)
        {
            api_error = err;
            pass = 0;
            goto cleanup;
        }

        err =
            lis_vector_set_all(
                9.0,
                v0);

        if(err)
        {
            api_error = err;
            pass = 0;
            goto cleanup;
        }

        err =
            lis_vector_get_values(
                solver->left_nullspace[0],
                0,
                3,
                after);

        if(err)
        {
            api_error = err;
            pass = 0;
            goto cleanup;
        }

        for(i=0;i<3;i++)
        {
            LIS_REAL delta =
                (LIS_REAL)fabs(
                    after[i]-before[i]);

            if(delta>ownership_delta)
            {
                ownership_delta = delta;
            }
        }

        if(ownership_delta>tol)
        {
            ownership_ok = 0;
        }
    }

    err =
        lis_vector_nrm2(
            solver->left_nullspace[0],
            &nrm0);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        lis_vector_nrm2(
            solver->left_nullspace[1],
            &nrm1);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        lis_vector_dot(
            solver->left_nullspace[0],
            solver->left_nullspace[1],
            &dot);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    dot_abs =
        (LIS_REAL)fabs(dot);

    orthonormal_ok =
        ((nrm0>=1.0 ? nrm0-1.0 : 1.0-nrm0)<=tol)
        &&
        ((nrm1>=1.0 ? nrm1-1.0 : 1.0-nrm1)<=tol)
        &&
        dot_abs<=tol;

    right_untouched_ok =
        solver->nullspace==NULL
        &&
        solver->nullspace_dim==0;

    err =
        make_vector3(
            A,
            1.0,
            -1.0,
            0.0,
            &replacement);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    one[0] = replacement;

    err =
        lis_solver_set_left_nullspace(
            solver,
            1,
            one);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        lis_vector_nrm2(
            solver->left_nullspace[0],
            &replacement_nrm);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    replacement_ok =
        solver->left_nullspace!=NULL
        &&
        solver->left_nullspace_dim==1
        &&
        solver->left_nullspace[0]!=replacement
        &&
        ((replacement_nrm>=1.0
            ? replacement_nrm-1.0
            : 1.0-replacement_nrm)<=tol);

    err =
        lis_solver_clear_left_nullspace(
            solver);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    clear_ok =
        solver->left_nullspace==NULL
        &&
        solver->left_nullspace_dim==0;

    err =
        lis_solver_set_left_nullspace(
            solver,
            1,
            one);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        lis_solver_set_left_nullspace(
            solver,
            0,
            NULL);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    zero_clear_ok =
        solver->left_nullspace==NULL
        &&
        solver->left_nullspace_dim==0;

    err =
        make_vector3(
            A,
            1.0,
            2.0,
            3.0,
            &dependent0);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            2.0,
            4.0,
            6.0,
            &dependent1);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    dependent[0] = dependent0;
    dependent[1] = dependent1;

    dependent_api =
        lis_solver_set_left_nullspace(
            solver,
            2,
            dependent);

    dependent_ok =
        dependent_api==LIS_ERR_ILL_ARG
        &&
        solver->left_nullspace==NULL
        &&
        solver->left_nullspace_dim==0;

    right_untouched_ok =
        right_untouched_ok
        &&
        solver->nullspace==NULL
        &&
        solver->nullspace_dim==0;

    if(
        !ownership_ok
        ||
        !orthonormal_ok
        ||
        !replacement_ok
        ||
        !clear_ok
        ||
        !zero_clear_ok
        ||
        !dependent_ok
        ||
        !right_untouched_ok)
    {
        pass = 0;
    }

cleanup:

    printf(
        "LEFT_NULLSPACE_API "
        "ownership=%d "
        "orthonormal=%d "
        "replacement=%d "
        "clear=%d "
        "zero_clear=%d "
        "dependent=%d "
        "right_untouched=%d "
        "api_error=%d "
        "dependent_api=%d "
        "owner_delta=% .3e "
        "nrm0=% .16e "
        "nrm1=% .16e "
        "dot=% .3e "
        "replacement_nrm=% .16e\n",
        ownership_ok,
        orthonormal_ok,
        replacement_ok,
        clear_ok,
        zero_clear_ok,
        dependent_ok,
        right_untouched_ok,
        (int)api_error,
        (int)dependent_api,
        (double)ownership_delta,
        (double)nrm0,
        (double)nrm1,
        (double)dot_abs,
        (double)replacement_nrm);

    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(dependent1)
    {
        lis_vector_destroy(
            dependent1);
    }

    if(dependent0)
    {
        lis_vector_destroy(
            dependent0);
    }

    if(replacement)
    {
        lis_vector_destroy(
            replacement);
    }

    if(v1)
    {
        lis_vector_destroy(
            v1);
    }

    if(v0)
    {
        lis_vector_destroy(
            v0);
    }

    return pass;
}

static int
run_rhs_compatibility_nonsymmetric(
    LIS_MATRIX A)
{
    LIS_SOLVER solver = NULL;

    LIS_VECTOR right_vector = NULL;
    LIS_VECTOR left_vector = NULL;
    LIS_VECTOR compatible_rhs = NULL;
    LIS_VECTOR incompatible_rhs = NULL;
    LIS_VECTOR zero_rhs = NULL;

    LIS_VECTOR right_basis[1];
    LIS_VECTOR left_basis[1];

    LIS_SCALAR right_coefficient = 0.0;
    LIS_SCALAR rhs_before[3];
    LIS_SCALAR rhs_after[3];

    LIS_REAL compatible_defect = 0.0;
    LIS_REAL compatible_relative = 0.0;
    LIS_REAL incompatible_defect = 0.0;
    LIS_REAL incompatible_relative = 0.0;
    LIS_REAL zero_defect = 0.0;
    LIS_REAL zero_relative = 0.0;
    LIS_REAL right_coefficient_abs = 0.0;
    LIS_REAL rhs_delta = 0.0;

    LIS_REAL expected_defect;
    LIS_REAL expected_relative;

    LIS_INT err;
    LIS_INT api_error = LIS_SUCCESS;
    LIS_INT missing_left_api = LIS_SUCCESS;

    int i;
    int missing_left_ok = 0;
    int compatible_ok = 0;
    int right_distinction_ok = 0;
    int incompatible_ok = 0;
    int zero_ok = 0;
    int rhs_unchanged_ok = 0;
    int state_ok = 0;

    int pass = 1;

    const LIS_REAL tol =
        (LIS_REAL)1.0e-10;

    /*
     * Conceptual nonsymmetric operator:
     *
     *       [ 1 0 0 ]
     *   A = [ 1 0 0 ]
     *       [ 0 0 1 ]
     *
     * ker(A)   = span([0,1,0]^T)
     * ker(A^H) = span([1,-1,0]^T).
     *
     * The compatibility query itself consumes the registered
     * left basis and RHS only, so the existing matrix A is used
     * here only as the vector-layout template.
     */

    expected_defect =
        (LIS_REAL)(1.0/sqrt(2.0));

    expected_relative =
        (LIS_REAL)(1.0/sqrt(10.0));

    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            0.0,
            1.0,
            0.0,
            &right_vector);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            -1.0,
            0.0,
            &left_vector);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            1.0,
            2.0,
            &compatible_rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            0.0,
            2.0,
            &incompatible_rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            0.0,
            0.0,
            0.0,
            &zero_rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    right_basis[0] = right_vector;
    left_basis[0] = left_vector;

    err =
        lis_solver_set_nullspace(
            solver,
            1,
            right_basis);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    /*
     * A registered right nullspace must not be used as an
     * implicit compatibility basis.
     */
    missing_left_api =
        lis_solver_check_rhs_compatibility(
            solver,
            compatible_rhs,
            &compatible_defect,
            &compatible_relative);

    missing_left_ok =
        missing_left_api==LIS_ERR_ILL_ARG;

    err =
        lis_solver_set_left_nullspace(
            solver,
            1,
            left_basis);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        lis_vector_dot(
            solver->nullspace[0],
            compatible_rhs,
            &right_coefficient);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    right_coefficient_abs =
        (LIS_REAL)fabs(right_coefficient);

    right_distinction_ok =
        fabs(right_coefficient_abs-1.0)<=tol;

    compatible_defect = -1.0;
    compatible_relative = -1.0;

    err =
        lis_solver_check_rhs_compatibility(
            solver,
            compatible_rhs,
            &compatible_defect,
            &compatible_relative);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    compatible_ok =
        compatible_defect<=tol
        &&
        compatible_relative<=tol;

    err =
        lis_vector_get_values(
            incompatible_rhs,
            0,
            3,
            rhs_before);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        lis_solver_check_rhs_compatibility(
            solver,
            incompatible_rhs,
            &incompatible_defect,
            &incompatible_relative);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    incompatible_ok =
        fabs(incompatible_defect-expected_defect)<=tol
        &&
        fabs(incompatible_relative-expected_relative)<=tol;

    err =
        lis_vector_get_values(
            incompatible_rhs,
            0,
            3,
            rhs_after);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    for(i=0;i<3;i++)
    {
        LIS_REAL delta =
            (LIS_REAL)fabs(
                rhs_after[i]-rhs_before[i]);

        if(delta>rhs_delta)
        {
            rhs_delta = delta;
        }
    }

    rhs_unchanged_ok =
        rhs_delta<=tol;

    err =
        lis_solver_check_rhs_compatibility(
            solver,
            zero_rhs,
            &zero_defect,
            &zero_relative);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    zero_ok =
        zero_defect<=tol
        &&
        zero_relative<=tol;

    state_ok =
        solver->nullspace!=NULL
        &&
        solver->nullspace_dim==1
        &&
        solver->left_nullspace!=NULL
        &&
        solver->left_nullspace_dim==1;

    if(
        !missing_left_ok
        ||
        !compatible_ok
        ||
        !right_distinction_ok
        ||
        !incompatible_ok
        ||
        !zero_ok
        ||
        !rhs_unchanged_ok
        ||
        !state_ok)
    {
        pass = 0;
    }

cleanup:

    printf(
        "RHS_COMPATIBILITY_NONSYMMETRIC "
        "missing_left=%d "
        "compatible=%d "
        "right_distinction=%d "
        "incompatible=%d "
        "zero=%d "
        "rhs_unchanged=%d "
        "state=%d "
        "api_error=%d "
        "missing_left_api=%d "
        "compatible_defect=% .16e "
        "right_coeff=% .16e "
        "incompatible_defect=% .16e "
        "incompatible_relative=% .16e "
        "rhs_delta=% .3e\n",
        missing_left_ok,
        compatible_ok,
        right_distinction_ok,
        incompatible_ok,
        zero_ok,
        rhs_unchanged_ok,
        state_ok,
        (int)api_error,
        (int)missing_left_api,
        (double)compatible_defect,
        (double)right_coefficient_abs,
        (double)incompatible_defect,
        (double)incompatible_relative,
        (double)rhs_delta);

    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(zero_rhs)
    {
        lis_vector_destroy(
            zero_rhs);
    }

    if(incompatible_rhs)
    {
        lis_vector_destroy(
            incompatible_rhs);
    }

    if(compatible_rhs)
    {
        lis_vector_destroy(
            compatible_rhs);
    }

    if(left_vector)
    {
        lis_vector_destroy(
            left_vector);
    }

    if(right_vector)
    {
        lis_vector_destroy(
            right_vector);
    }

    return pass;
}

static int
run_rhs_compatibility_multivector(
    LIS_MATRIX A)
{
    LIS_SOLVER solver = NULL;

    LIS_VECTOR w0 = NULL;
    LIS_VECTOR w1 = NULL;
    LIS_VECTOR rhs = NULL;

    LIS_VECTOR basis[2];

    LIS_REAL defect = 0.0;
    LIS_REAL relative_defect = 0.0;

    LIS_INT err;
    LIS_INT api_error = LIS_SUCCESS;

    int dimension_ok = 0;
    int value_ok = 0;
    int pass = 1;

    const LIS_REAL tol =
        (LIS_REAL)1.0e-10;

    const LIS_REAL expected_defect =
        (LIS_REAL)5.0;

    const LIS_REAL expected_relative =
        (LIS_REAL)(5.0/13.0);

    /*
     * Orthonormal two-vector left basis:
     *
     *   w0 = [1,0,0]^T
     *   w1 = [0,1,0]^T
     *
     * For rhs = [3,4,12]^T:
     *
     *   ||W^H rhs||_2 = sqrt(3^2 + 4^2) = 5
     *   ||rhs||_2     = 13
     *   relative      = 5/13.
     */

    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            0.0,
            0.0,
            &w0);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            0.0,
            1.0,
            0.0,
            &w1);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            3.0,
            4.0,
            12.0,
            &rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    basis[0] = w0;
    basis[1] = w1;

    err =
        lis_solver_set_left_nullspace(
            solver,
            2,
            basis);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    dimension_ok =
        solver->left_nullspace!=NULL
        &&
        solver->left_nullspace_dim==2;

    err =
        lis_solver_check_rhs_compatibility(
            solver,
            rhs,
            &defect,
            &relative_defect);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    value_ok =
        fabs(defect-expected_defect)<=tol
        &&
        fabs(relative_defect-expected_relative)<=tol;

    if(
        !dimension_ok
        ||
        !value_ok)
    {
        pass = 0;
    }

cleanup:

    printf(
        "RHS_COMPATIBILITY_MULTIVECTOR "
        "dimension=%d "
        "value=%d "
        "api_error=%d "
        "defect=% .16e "
        "relative=% .16e\n",
        dimension_ok,
        value_ok,
        (int)api_error,
        (double)defect,
        (double)relative_defect);

    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(rhs)
    {
        lis_vector_destroy(
            rhs);
    }

    if(w1)
    {
        lis_vector_destroy(
            w1);
    }

    if(w0)
    {
        lis_vector_destroy(
            w0);
    }

    return pass;
}

static int
run_rhs_projection_nonsymmetric(
    LIS_MATRIX A)
{
    LIS_SOLVER solver = NULL;

    LIS_VECTOR left_vector = NULL;
    LIS_VECTOR rhs = NULL;
    LIS_VECTOR projected_rhs = NULL;

    LIS_VECTOR left_basis[1];

    LIS_SCALAR rhs_before[3];
    LIS_SCALAR rhs_after[3];

    LIS_REAL pre_defect = 0.0;
    LIS_REAL pre_relative = 0.0;
    LIS_REAL post_defect = 0.0;
    LIS_REAL post_relative = 0.0;

    LIS_REAL rhs_delta = 0.0;
    LIS_REAL projected_delta = 0.0;

    LIS_INT err;
    LIS_INT api_error = LIS_SUCCESS;

    int i;
    int pre_incompatible_ok = 0;
    int projected_values_ok = 0;
    int post_compatible_ok = 0;
    int rhs_unchanged_ok = 0;
    int state_ok = 0;

    int pass = 1;

    const LIS_REAL tol =
        (LIS_REAL)1.0e-10;

    const LIS_REAL expected_pre_defect =
        (LIS_REAL)(1.0/sqrt(2.0));

    const LIS_REAL expected_pre_relative =
        (LIS_REAL)(1.0/sqrt(10.0));

    /*
     * Conceptual nonsymmetric operator:
     *
     *       [ 1 0 0 ]
     *   A = [ 1 0 0 ]
     *       [ 0 0 1 ]
     *
     * Left nullspace:
     *
     *   W = span([1,-1,0]^T).
     *
     * After normalization:
     *
     *   w = [1,-1,0]^T / sqrt(2).
     *
     * For
     *
     *   b = [1,0,2]^T
     *
     * the incompatible component is
     *
     *   w^H b = 1/sqrt(2),
     *
     * and therefore
     *
     *   b_p = b - w(w^H b)
     *       = [0.5,0.5,2]^T.
     *
     * The matrix A is used only as the vector-layout
     * template; the projection consumes the registered
     * left-nullspace basis and RHS.
     */

    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            -1.0,
            0.0,
            &left_vector);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            0.0,
            2.0,
            &rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            7.0,
            8.0,
            9.0,
            &projected_rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    for(i=0;i<3;i++)
    {
        rhs_before[i] = rhs->value[i];
    }

    left_basis[0] = left_vector;

    err =
        lis_solver_set_left_nullspace(
            solver,
            1,
            left_basis);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        lis_solver_check_rhs_compatibility(
            solver,
            rhs,
            &pre_defect,
            &pre_relative);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    pre_incompatible_ok =
        fabs(pre_defect-expected_pre_defect)<=tol
        &&
        fabs(pre_relative-expected_pre_relative)<=tol;

    err =
        lis_solver_project_rhs_compatible(
            solver,
            rhs,
            projected_rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    for(i=0;i<3;i++)
    {
        rhs_after[i] = rhs->value[i];

        rhs_delta +=
            (LIS_REAL)fabs(
                rhs_after[i]-rhs_before[i]);
    }

    rhs_unchanged_ok =
        rhs_delta<=tol;

    projected_delta =
        (LIS_REAL)fabs(
            projected_rhs->value[0]
            -(LIS_SCALAR)0.5)
        +
        (LIS_REAL)fabs(
            projected_rhs->value[1]
            -(LIS_SCALAR)0.5)
        +
        (LIS_REAL)fabs(
            projected_rhs->value[2]
            -(LIS_SCALAR)2.0);

    projected_values_ok =
        projected_delta<=tol;

    err =
        lis_solver_check_rhs_compatibility(
            solver,
            projected_rhs,
            &post_defect,
            &post_relative);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    post_compatible_ok =
        post_defect<=tol
        &&
        post_relative<=tol;

    state_ok =
        solver->left_nullspace!=NULL
        &&
        solver->left_nullspace_dim==1;

    if(
        !pre_incompatible_ok
        ||
        !projected_values_ok
        ||
        !post_compatible_ok
        ||
        !rhs_unchanged_ok
        ||
        !state_ok)
    {
        pass = 0;
    }

cleanup:

    printf(
        "RHS_PROJECTION_NONSYMMETRIC "
        "pre_incompatible=%d "
        "values=%d "
        "post_compatible=%d "
        "rhs_unchanged=%d "
        "state=%d "
        "api_error=%d "
        "pre_defect=% .16e "
        "pre_relative=% .16e "
        "post_defect=% .16e "
        "post_relative=% .16e "
        "rhs_delta=% .3e "
        "projected_delta=% .3e "
        "p0=% .16e "
        "p1=% .16e "
        "p2=% .16e\n",
        pre_incompatible_ok,
        projected_values_ok,
        post_compatible_ok,
        rhs_unchanged_ok,
        state_ok,
        (int)api_error,
        (double)pre_defect,
        (double)pre_relative,
        (double)post_defect,
        (double)post_relative,
        (double)rhs_delta,
        (double)projected_delta,
        projected_rhs
            ? (double)projected_rhs->value[0]
            : 0.0,
        projected_rhs
            ? (double)projected_rhs->value[1]
            : 0.0,
        projected_rhs
            ? (double)projected_rhs->value[2]
            : 0.0);

    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(projected_rhs)
    {
        lis_vector_destroy(
            projected_rhs);
    }

    if(rhs)
    {
        lis_vector_destroy(
            rhs);
    }

    if(left_vector)
    {
        lis_vector_destroy(
            left_vector);
    }

    return pass;
}


static int
run_rhs_projection_multivector(
    LIS_MATRIX A)
{
    LIS_SOLVER solver = NULL;

    LIS_VECTOR w0 = NULL;
    LIS_VECTOR w1 = NULL;
    LIS_VECTOR rhs = NULL;
    LIS_VECTOR projected_rhs = NULL;

    LIS_VECTOR basis[2];

    LIS_SCALAR rhs_before[3];

    LIS_REAL pre_defect = 0.0;
    LIS_REAL pre_relative = 0.0;
    LIS_REAL post_defect = 0.0;
    LIS_REAL post_relative = 0.0;

    LIS_REAL rhs_delta = 0.0;
    LIS_REAL projected_delta = 0.0;

    LIS_INT err;
    LIS_INT api_error = LIS_SUCCESS;

    int i;
    int dimension_ok = 0;
    int pre_ok = 0;
    int values_ok = 0;
    int post_ok = 0;
    int rhs_unchanged_ok = 0;

    int pass = 1;

    const LIS_REAL tol =
        (LIS_REAL)1.0e-10;

    const LIS_REAL expected_defect =
        (LIS_REAL)5.0;

    const LIS_REAL expected_relative =
        (LIS_REAL)(5.0/13.0);

    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            0.0,
            0.0,
            &w0);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            0.0,
            1.0,
            0.0,
            &w1);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            3.0,
            4.0,
            12.0,
            &rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            7.0,
            8.0,
            9.0,
            &projected_rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    for(i=0;i<3;i++)
    {
        rhs_before[i] =
            rhs->value[i];
    }

    basis[0] = w0;
    basis[1] = w1;

    err =
        lis_solver_set_left_nullspace(
            solver,
            2,
            basis);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    dimension_ok =
        solver->left_nullspace_dim==2;

    err =
        lis_solver_check_rhs_compatibility(
            solver,
            rhs,
            &pre_defect,
            &pre_relative);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    pre_ok =
        fabs(pre_defect-expected_defect)<=tol
        &&
        fabs(pre_relative-expected_relative)<=tol;

    err =
        lis_solver_project_rhs_compatible(
            solver,
            rhs,
            projected_rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    for(i=0;i<3;i++)
    {
        rhs_delta +=
            (LIS_REAL)fabs(
                rhs->value[i]
                -rhs_before[i]);
    }

    rhs_unchanged_ok =
        rhs_delta<=tol;

    projected_delta =
        (LIS_REAL)fabs(
            projected_rhs->value[0])
        +
        (LIS_REAL)fabs(
            projected_rhs->value[1])
        +
        (LIS_REAL)fabs(
            projected_rhs->value[2]
            -(LIS_SCALAR)12.0);

    values_ok =
        projected_delta<=tol;

    err =
        lis_solver_check_rhs_compatibility(
            solver,
            projected_rhs,
            &post_defect,
            &post_relative);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    post_ok =
        post_defect<=tol
        &&
        post_relative<=tol;

    if(
        !dimension_ok
        ||
        !pre_ok
        ||
        !values_ok
        ||
        !post_ok
        ||
        !rhs_unchanged_ok)
    {
        pass = 0;
    }

cleanup:

    printf(
        "RHS_PROJECTION_MULTIVECTOR "
        "dimension=%d "
        "pre=%d "
        "values=%d "
        "post=%d "
        "rhs_unchanged=%d "
        "api_error=%d "
        "pre_defect=% .16e "
        "pre_relative=% .16e "
        "post_defect=% .16e "
        "post_relative=% .16e "
        "rhs_delta=% .3e "
        "projected_delta=% .3e "
        "p0=% .16e "
        "p1=% .16e "
        "p2=% .16e\n",
        dimension_ok,
        pre_ok,
        values_ok,
        post_ok,
        rhs_unchanged_ok,
        (int)api_error,
        (double)pre_defect,
        (double)pre_relative,
        (double)post_defect,
        (double)post_relative,
        (double)rhs_delta,
        (double)projected_delta,
        projected_rhs
            ? (double)projected_rhs->value[0]
            : 0.0,
        projected_rhs
            ? (double)projected_rhs->value[1]
            : 0.0,
        projected_rhs
            ? (double)projected_rhs->value[2]
            : 0.0);

    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(projected_rhs)
    {
        lis_vector_destroy(
            projected_rhs);
    }

    if(rhs)
    {
        lis_vector_destroy(
            rhs);
    }

    if(w1)
    {
        lis_vector_destroy(
            w1);
    }

    if(w0)
    {
        lis_vector_destroy(
            w0);
    }

    return pass;
}


static int
run_rhs_projection_alias_rejection(
    LIS_MATRIX A)
{
    LIS_SOLVER solver = NULL;

    LIS_VECTOR w = NULL;
    LIS_VECTOR rhs = NULL;

    LIS_VECTOR basis[1];

    LIS_SCALAR rhs_before[3];
    LIS_SCALAR basis_before[3];

    LIS_REAL rhs_delta = 0.0;
    LIS_REAL basis_delta = 0.0;

    LIS_INT err;
    LIS_INT same_api = LIS_SUCCESS;
    LIS_INT basis_api = LIS_SUCCESS;
    LIS_INT api_error = LIS_SUCCESS;

    int i;
    int same_object_ok = 0;
    int basis_alias_ok = 0;
    int rhs_unchanged_ok = 0;
    int basis_unchanged_ok = 0;
    int state_ok = 0;

    int pass = 1;

    const LIS_REAL tol =
        (LIS_REAL)1.0e-10;

    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            -1.0,
            0.0,
            &w);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    err =
        make_vector3(
            A,
            1.0,
            0.0,
            2.0,
            &rhs);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    basis[0] = w;

    err =
        lis_solver_set_left_nullspace(
            solver,
            1,
            basis);

    if(err)
    {
        api_error = err;
        pass = 0;
        goto cleanup;
    }

    for(i=0;i<3;i++)
    {
        rhs_before[i] =
            rhs->value[i];

        basis_before[i] =
            solver->left_nullspace[0]->value[i];
    }

    same_api =
        lis_solver_project_rhs_compatible(
            solver,
            rhs,
            rhs);

    same_object_ok =
        same_api==LIS_ERR_ILL_ARG;

    basis_api =
        lis_solver_project_rhs_compatible(
            solver,
            rhs,
            solver->left_nullspace[0]);

    basis_alias_ok =
        basis_api==LIS_ERR_ILL_ARG;

    for(i=0;i<3;i++)
    {
        rhs_delta +=
            (LIS_REAL)fabs(
                rhs->value[i]
                -rhs_before[i]);

        basis_delta +=
            (LIS_REAL)fabs(
                solver->left_nullspace[0]->value[i]
                -basis_before[i]);
    }

    rhs_unchanged_ok =
        rhs_delta<=tol;

    basis_unchanged_ok =
        basis_delta<=tol;

    state_ok =
        solver->left_nullspace!=NULL
        &&
        solver->left_nullspace_dim==1;

    if(
        !same_object_ok
        ||
        !basis_alias_ok
        ||
        !rhs_unchanged_ok
        ||
        !basis_unchanged_ok
        ||
        !state_ok)
    {
        pass = 0;
    }

cleanup:

    printf(
        "RHS_PROJECTION_ALIAS "
        "same_object=%d "
        "basis_alias=%d "
        "rhs_unchanged=%d "
        "basis_unchanged=%d "
        "state=%d "
        "api_error=%d "
        "same_api=%d "
        "basis_api=%d "
        "rhs_delta=% .3e "
        "basis_delta=% .3e\n",
        same_object_ok,
        basis_alias_ok,
        rhs_unchanged_ok,
        basis_unchanged_ok,
        state_ok,
        (int)api_error,
        (int)same_api,
        (int)basis_api,
        (double)rhs_delta,
        (double)basis_delta);

    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(rhs)
    {
        lis_vector_destroy(
            rhs);
    }

    if(w)
    {
        lis_vector_destroy(
            w);
    }

    return pass;
}


static LIS_INT
configure_rhs_policy_solver(
    LIS_VECTOR null_vector,
    LIS_SOLVER *solver)
{
    LIS_VECTOR basis[1];
    LIS_INT err;

    char options[] =
        "-i gmres "
        "-p jacobi "
        "-scale none "
        "-tol 1.0e-12 "
        "-maxiter 100 "
        "-restart 6";

    *solver = NULL;

    err =
        lis_solver_create(
            solver);

    if(err)
    {
        return err;
    }

    err =
        lis_solver_set_option(
            options,
            *solver);

    if(err)
    {
        goto fail;
    }

    basis[0] = null_vector;

    err =
        lis_solver_set_nullspace(
            *solver,
            1,
            basis);

    if(err)
    {
        goto fail;
    }

    err =
        lis_solver_set_left_nullspace(
            *solver,
            1,
            basis);

    if(err)
    {
        goto fail;
    }

    return LIS_SUCCESS;

fail:

    lis_solver_destroy(
        *solver);

    *solver = NULL;

    return err;
}


static int
run_rhs_policy_solve(
    LIS_MATRIX A,
    LIS_VECTOR compatible_rhs,
    LIS_VECTOR null_vector)
{
    LIS_SOLVER solver = NULL;
    LIS_SOLVER sequence_solver = NULL;

    LIS_VECTOR incompatible_rhs = NULL;
    LIS_VECTOR projected_rhs = NULL;
    LIS_VECTOR residual = NULL;
    LIS_VECTOR x = NULL;

    LIS_INT err;
    LIS_INT api;
    LIS_INT policy = -1;

    LIS_REAL tolerance = 0.0;
    LIS_REAL project_residual = 0.0;
    LIS_REAL sequence_residual = 0.0;

    LIS_SCALAR rhs_before[3];

    int default_ok = 0;
    int invalid_tolerance_ok = 0;
    int assume_ok = 0;
    int reject_ok = 0;
    int project_ok = 0;
    int sequence_reject_ok = 0;
    int sequence_project_ok = 0;
    int rhs_unchanged_ok = 0;

    int pass = 1;
    int i;

    err =
        make_vector3(
            A,
            1.0,
            -2.0,
            2.0,
            &incompatible_rhs);

    if(err)
    {
        return 0;
    }

    err =
        lis_vector_duplicate(
            A,
            &projected_rhs);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    err =
        lis_vector_duplicate(
            A,
            &residual);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    err =
        lis_vector_duplicate(
            A,
            &x);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    for(i=0;i<3;i++)
    {
        rhs_before[i] =
            incompatible_rhs->value[i];
    }

    /*
     * Ordinary lis_solve() policy path.
     */
    err =
        configure_rhs_policy_solver(
            null_vector,
            &solver);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    err =
        lis_solver_get_rhs_compatibility_policy(
            solver,
            &policy);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    err =
        lis_solver_get_rhs_compatibility_tolerance(
            solver,
            &tolerance);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    default_ok =
        policy==
            LIS_RHS_POLICY_ASSUME_COMPATIBLE
        &&
        fabs(
            tolerance-
            (LIS_REAL)1.0e-12)
            <=
            (LIS_REAL)1.0e-18;

    err =
        lis_solver_set_rhs_compatibility_tolerance(
            solver,
            (LIS_REAL)2.0e-12);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    api =
        lis_solver_set_rhs_compatibility_tolerance(
            solver,
            (LIS_REAL)INFINITY);

    err =
        lis_solver_get_rhs_compatibility_tolerance(
            solver,
            &tolerance);

    invalid_tolerance_ok =
        api==LIS_ERR_ILL_ARG
        &&
        err==LIS_SUCCESS
        &&
        fabs(
            tolerance-
            (LIS_REAL)2.0e-12)
            <=
            (LIS_REAL)1.0e-18;

    lis_vector_set_all(
        0.0,
        x);

    api =
        lis_solve(
            A,
            compatible_rhs,
            x,
            solver);

    assume_ok =
        api==LIS_SUCCESS
        &&
        solver->retcode==LIS_SUCCESS;

    err =
        lis_solver_set_rhs_compatibility_policy(
            solver,
            LIS_RHS_POLICY_REJECT);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    err =
        lis_solver_set_rhs_compatibility_tolerance(
            solver,
            (LIS_REAL)1.0e-12);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    lis_vector_set_all(
        0.0,
        x);

    api =
        lis_solve(
            A,
            incompatible_rhs,
            x,
            solver);

    reject_ok =
        api==LIS_ERR_ILL_ARG
        &&
        solver->retcode==
            LIS_ERR_ILL_ARG;

    err =
        lis_solver_set_rhs_compatibility_policy(
            solver,
            LIS_RHS_POLICY_PROJECT);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    err =
        lis_solver_project_rhs_compatible(
            solver,
            incompatible_rhs,
            projected_rhs);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    lis_vector_set_all(
        0.0,
        x);

    api =
        lis_solve(
            A,
            incompatible_rhs,
            x,
            solver);

    if(
        api==LIS_SUCCESS
        &&
        solver->retcode==LIS_SUCCESS)
    {
        err =
            lis_matvec(
                A,
                x,
                residual);

        if(!err)
        {
            err =
                lis_vector_axpy(
                    -1.0,
                    projected_rhs,
                    residual);
        }

        if(!err)
        {
            err =
                lis_vector_nrm2(
                    residual,
                    &project_residual);
        }

        if(!err)
        {
            project_ok =
                project_residual
                <
                (LIS_REAL)1.0e-10;
        }
    }

    /*
     * Persistent sequence path.
     *
     * REJECT must stop before an attempt is made.
     * PROJECT must use the compatible effective RHS.
     */
    err =
        configure_rhs_policy_solver(
            null_vector,
            &sequence_solver);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    err =
        lis_solver_set_rhs_compatibility_policy(
            sequence_solver,
            LIS_RHS_POLICY_REJECT);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    lis_vector_set_all(
        0.0,
        x);

    api =
        lis_solve_sequence(
            A,
            incompatible_rhs,
            x,
            sequence_solver);

    sequence_reject_ok =
        api==LIS_ERR_ILL_ARG
        &&
        sequence_solver->retcode==
            LIS_ERR_ILL_ARG
        &&
        sequence_solver->
            sequence_diag_attempts==0;

    err =
        lis_solver_set_rhs_compatibility_policy(
            sequence_solver,
            LIS_RHS_POLICY_PROJECT);

    if(err)
    {
        pass = 0;
        goto cleanup;
    }

    lis_vector_set_all(
        0.0,
        x);

    api =
        lis_solve_sequence(
            A,
            incompatible_rhs,
            x,
            sequence_solver);

    if(
        api==LIS_SUCCESS
        &&
        sequence_solver->retcode==
            LIS_SUCCESS
        &&
        sequence_solver->
            sequence_diag_attempts>=1)
    {
        err =
            lis_matvec(
                A,
                x,
                residual);

        if(!err)
        {
            err =
                lis_vector_axpy(
                    -1.0,
                    projected_rhs,
                    residual);
        }

        if(!err)
        {
            err =
                lis_vector_nrm2(
                    residual,
                    &sequence_residual);
        }

        if(!err)
        {
            sequence_project_ok =
                sequence_residual
                <
                (LIS_REAL)1.0e-10;
        }
    }

    rhs_unchanged_ok = 1;

    for(i=0;i<3;i++)
    {
        if(
            incompatible_rhs->value[i]
            !=
            rhs_before[i])
        {
            rhs_unchanged_ok = 0;
        }
    }

    printf(
        "RHS_POLICY_SOLVE "
        "default=%d "
        "inf_tolerance=%d "
        "assume=%d "
        "reject=%d "
        "project=%d "
        "sequence_reject=%d "
        "sequence_project=%d "
        "rhs_unchanged=%d "
        "project_residual=% .16e "
        "sequence_residual=% .16e\n",
        default_ok,
        invalid_tolerance_ok,
        assume_ok,
        reject_ok,
        project_ok,
        sequence_reject_ok,
        sequence_project_ok,
        rhs_unchanged_ok,
        (double)project_residual,
        (double)sequence_residual);

    if(
        !default_ok
        ||
        !invalid_tolerance_ok
        ||
        !assume_ok
        ||
        !reject_ok
        ||
        !project_ok
        ||
        !sequence_reject_ok
        ||
        !sequence_project_ok
        ||
        !rhs_unchanged_ok)
    {
        pass = 0;
    }

cleanup:

    if(sequence_solver)
    {
        lis_solver_destroy(
            sequence_solver);
    }

    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(x)
    {
        lis_vector_destroy(
            x);
    }

    if(residual)
    {
        lis_vector_destroy(
            residual);
    }

    if(projected_rhs)
    {
        lis_vector_destroy(
            projected_rhs);
    }

    if(incompatible_rhs)
    {
        lis_vector_destroy(
            incompatible_rhs);
    }

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

    int gmres_ok;
    int fgmres_ok;
    int gcrodr_ok;
    int fgcrodr_ok;

    int solver_reject_ok;
    int scale_reject_ok;

    int incompatible_gmres_ok;
    int incompatible_fgmres_ok;
    int incompatible_gcrodr_ok;
    int incompatible_fgcrodr_ok;

    int gcrodr_repeat_ok;
    int left_api_ok;
    int rhs_compat_ok;
    int rhs_multi_ok;
    int rhs_project_ok;
    int rhs_project_multi_ok;
    int rhs_project_alias_ok;
    int rhs_policy_ok;

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

    gmres_ok =
        run_supported(
            A,
            b,
            q,
            "gmres");

    fgmres_ok =
        run_supported(
            A,
            b,
            q,
            "fgmres");

    gcrodr_ok =
        run_supported(
            A,
            b,
            q,
            "gcrodr");

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


    incompatible_gmres_ok =
        run_incompatible_rhs(
            A,
            q,
            "gmres");

    incompatible_fgmres_ok =
        run_incompatible_rhs(
            A,
            q,
            "fgmres");

    incompatible_gcrodr_ok =
        run_incompatible_rhs(
            A,
            q,
            "gcrodr");

    incompatible_fgcrodr_ok =
        run_incompatible_rhs(
            A,
            q,
            "fgcrodr");

    gcrodr_repeat_ok =
        run_gcrodr_repeated(
            A,
            q);

    left_api_ok =
        run_left_nullspace_api(
            A);

    rhs_compat_ok =
        run_rhs_compatibility_nonsymmetric(
            A);

    rhs_multi_ok =
        run_rhs_compatibility_multivector(
            A);

    rhs_project_ok =
        run_rhs_projection_nonsymmetric(
            A);

    rhs_project_multi_ok =
        run_rhs_projection_multivector(
            A);

    rhs_project_alias_ok =
        run_rhs_projection_alias_rejection(
            A);

    rhs_policy_ok =
        run_rhs_policy_solve(
            A,
            b,
            q);

    printf(
        "NULLSPACE_SUMMARY "
        "gmres=%d "
        "fgmres=%d "
        "gcrodr=%d "
        "fgcrodr=%d "
        "solver_reject=%d "
        "scale_reject=%d "
        "incompatible_gmres=%d "
        "incompatible_fgmres=%d "
        "incompatible_gcrodr=%d "
        "incompatible_fgcrodr=%d "
        "gcrodr_repeat=%d "
        "left_api=%d "
        "rhs_compat=%d "
        "rhs_multi=%d "
        "rhs_project=%d "
        "rhs_project_multi=%d "
        "rhs_project_alias=%d "
        "rhs_policy=%d\n",
        gmres_ok,
        fgmres_ok,
        gcrodr_ok,
        fgcrodr_ok,
        solver_reject_ok,
        scale_reject_ok,
        incompatible_gmres_ok,
        incompatible_fgmres_ok,
        incompatible_gcrodr_ok,
        incompatible_fgcrodr_ok,
        gcrodr_repeat_ok,
        left_api_ok,
        rhs_compat_ok,
        rhs_multi_ok,
        rhs_project_ok,
        rhs_project_multi_ok,
        rhs_project_alias_ok,
        rhs_policy_ok);

    if(
        !gmres_ok
        ||
        !fgmres_ok
        ||
        !gcrodr_ok
        ||
        !fgcrodr_ok
        ||
        !solver_reject_ok
        ||
        !scale_reject_ok
        ||
        !incompatible_gmres_ok
        ||
        !incompatible_fgmres_ok
        ||
        !incompatible_gcrodr_ok
        ||
        !incompatible_fgcrodr_ok
        ||
        !gcrodr_repeat_ok
        ||
        !left_api_ok
        ||
        !rhs_compat_ok
        ||
        !rhs_multi_ok
        ||
        !rhs_project_ok
        ||
        !rhs_project_multi_ok
        ||
        !rhs_project_alias_ok
        ||
        !rhs_policy_ok)
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
