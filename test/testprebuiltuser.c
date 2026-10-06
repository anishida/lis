#include "lislib.h"
#include <math.h>
#include <stdio.h>

typedef struct
{
    LIS_INT n;
    LIS_SCALAR diag;
    LIS_SCALAR offdiag;
    unsigned long matvec_calls;
} USER_TRIDIAG;

static LIS_INT shell_matvec(void *ctx, const LIS_SCALAR *x, LIS_SCALAR *y)
{
    USER_TRIDIAG *A = (USER_TRIDIAG *)ctx;
    LIS_INT i;
    A->matvec_calls++;
    if( A->n==1 ) { y[0] = A->diag*x[0]; return LIS_SUCCESS; }
    y[0] = A->diag*x[0] + A->offdiag*x[1];
    for(i=1;i<A->n-1;i++) y[i] = A->offdiag*x[i-1] + A->diag*x[i] + A->offdiag*x[i+1];
    y[A->n-1] = A->offdiag*x[A->n-2] + A->diag*x[A->n-1];
    return LIS_SUCCESS;
}

static LIS_INT shell_matvech(void *ctx, const LIS_SCALAR *x, LIS_SCALAR *y)
{
    return shell_matvec(ctx,x,y);
}

int main(int argc, char **argv)
{
    const LIS_INT n = 64;
    USER_TRIDIAG shell;
    LIS_MATRIX C = NULL, U = NULL;
    LIS_VECTOR b = NULL, x = NULL, exact = NULL;
    LIS_SOLVER plain = NULL, solver = NULL;
    LIS_PRECON precon = NULL;
    LIS_INT i, err, status;
    LIS_REAL nrm2;
    unsigned long calls_after_first;

    err = lis_initialize(&argc,&argv);
    if( err ) return 1;

    shell.n = n; shell.diag = 4.0; shell.offdiag = -1.0; shell.matvec_calls = 0;

    err = lis_matrix_create(LIS_COMM_WORLD,&C); if( err ) goto api_fail;
    err = lis_matrix_set_size(C,0,n); if( err ) goto api_fail;
    err = lis_matrix_set_type(C,LIS_MATRIX_CSR); if( err ) goto api_fail;
    for(i=0;i<n;i++)
    {
        err = lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,C); if( err ) goto api_fail;
        if(i>0) { err = lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,C); if( err ) goto api_fail; }
        if(i<n-1) { err = lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,C); if( err ) goto api_fail; }
    }
    err = lis_matrix_assemble(C); if( err ) goto api_fail;

    err = lis_matrix_create(LIS_COMM_WORLD,&U); if( err ) goto api_fail;
    err = lis_matrix_set_size(U,0,n); if( err ) goto api_fail;
    err = lis_matrix_set_user(U,&shell,shell_matvec,shell_matvech); if( err ) goto api_fail;

    err = lis_vector_duplicate(U,&b); if( err ) goto api_fail;
    err = lis_vector_duplicate(U,&x); if( err ) goto api_fail;
    err = lis_vector_duplicate(U,&exact); if( err ) goto api_fail;
    for(i=0;i<n;i++) exact->value[i] = 1.0;
    err = lis_matvec(U,exact,b); if( err ) goto api_fail;

    err = lis_solver_create(&plain); if( err ) goto api_fail;
    err = lis_solver_set_option("-i gmres -p ilu -fill 0 -scale none -tol 1.0e-12 -maxiter 200 -maxiter_noimp 0 -print none",plain); if( err ) goto api_fail;
    lis_vector_set_all(0.0,x);
    err = lis_solve(U,b,x,plain);
    if( err!=LIS_ERR_NOT_IMPLEMENTED )
    {
        fprintf(stderr,"lis_solve(USER,-p ilu) returned %d, expected %d\n",(int)err,(int)LIS_ERR_NOT_IMPLEMENTED);
        goto test_fail;
    }
    lis_solver_destroy(plain); plain = NULL;

    err = lis_solver_create(&solver); if( err ) goto api_fail;
    err = lis_solver_set_option("-i gmres -p ilu -fill 0 -scale none -tol 1.0e-12 -maxiter 200 -maxiter_noimp 0 -print none",solver); if( err ) goto api_fail;
    err = lis_solver_set_matrix(C,solver); if( err ) goto api_fail;
    err = lis_precon_psd_create(solver,&precon); if( err ) goto api_fail;
    err = lis_precon_psd_update(solver,precon); if( err ) goto api_fail;

    shell.matvec_calls = 0;
    lis_vector_set_all(0.0,x);
    err = lis_solve_kernel(U,b,x,solver,precon); if( err ) goto api_fail;
    err = lis_solver_get_status(solver,&status); if( err ) goto api_fail;
    if( status!=LIS_SUCCESS ) { fprintf(stderr,"first kernel solve status=%d\n",(int)status); goto test_fail; }
    if( shell.matvec_calls==0 ) { fprintf(stderr,"USER matvec was not used by first kernel solve\n"); goto test_fail; }
    lis_vector_axpy(-1.0,exact,x); lis_vector_nrm2(x,&nrm2);
    if( nrm2>1.0e-8 ) { fprintf(stderr,"first solution error %.15e\n",(double)nrm2); goto test_fail; }

    calls_after_first = shell.matvec_calls;
    lis_vector_set_all(0.0,x);
    err = lis_solve_kernel(U,b,x,solver,precon); if( err ) goto api_fail;
    err = lis_solver_get_status(solver,&status); if( err ) goto api_fail;
    if( status!=LIS_SUCCESS ) { fprintf(stderr,"second kernel solve status=%d\n",(int)status); goto test_fail; }
    if( shell.matvec_calls<=calls_after_first ) { fprintf(stderr,"USER matvec was not used by second kernel solve\n"); goto test_fail; }
    lis_vector_axpy(-1.0,exact,x); lis_vector_nrm2(x,&nrm2);
    if( nrm2>1.0e-8 ) { fprintf(stderr,"second solution error %.15e\n",(double)nrm2); goto test_fail; }

    printf("LIS_PREBUILT_USER_STAGE5D PASSED matvec_calls=%lu\n", shell.matvec_calls);
    lis_precon_destroy(precon);
    lis_solver_destroy(solver);
    lis_vector_destroy(exact); lis_vector_destroy(x); lis_vector_destroy(b);
    lis_matrix_destroy(U); lis_matrix_destroy(C);
    lis_finalize();
    return 0;

api_fail:
    fprintf(stderr,"API failure: %d\n",(int)err);
test_fail:
    if(precon) lis_precon_destroy(precon);
    if(solver) lis_solver_destroy(solver);
    if(plain) lis_solver_destroy(plain);
    if(exact) lis_vector_destroy(exact);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(U) lis_matrix_destroy(U);
    if(C) lis_matrix_destroy(C);
    lis_finalize();
    return 1;
}
