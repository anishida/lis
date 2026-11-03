/* Integration regression: native Krylov selected through LIS eigensolver API. */
#include "lis.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 12

typedef struct {
    LIS_SCALAR diagonal[N];
} user_diag;

static LIS_INT apply_diagonal(void *opaque,
                               const LIS_SCALAR *x, LIS_SCALAR *y)
{
    user_diag *d = (user_diag *)opaque;
    LIS_INT i;
    for (i = 0; i < N; ++i) y[i] = d->diagonal[i] * x[i];
    return LIS_SUCCESS;
}

static LIS_INT make_matrix(user_diag *data, LIS_MATRIX *A)
{
    LIS_INT err;
    *A = NULL;
    err = lis_matrix_create(LIS_COMM_WORLD, A);
    if (err) return err;
    err = lis_matrix_set_size(*A, 0, N);
    if (err) return err;
    return lis_matrix_set_user(*A, data, apply_diagonal, apply_diagonal);
}

static int solve_check(LIS_MATRIX A, LIS_MATRIX B,
                       const char *options, double first, double second)
{
    LIS_VECTOR x = NULL;
    LIS_ESOLVER es = NULL;
    LIS_SCALAR eig0, eig1;
    LIS_INT err, status;
    int fail = 1;
    const char *phase = "allocate vector";

    err = lis_vector_duplicate(A, &x);
    if (err) goto cleanup;
    phase = "initialize vector";
    err = lis_vector_set_all(1.0, x);
    if (err) goto cleanup;
    phase = "create eigensolver";
    err = lis_esolver_create(&es);
    if (err) goto cleanup;
    phase = "set options";
    err = lis_esolver_set_option((char *)options, es);
    if (err) goto cleanup;
    phase = "solve eigenproblem";
    err = B ? lis_gesolve(A, B, x, &eig0, es) :
              lis_esolve(A, x, &eig0, es);
    if (err) goto cleanup;
    phase = "get status";
    err = lis_esolver_get_status(es, &status);
    if (err || status != LIS_SUCCESS) goto cleanup;
    phase = "get second eigenvalue";
    err = lis_esolver_get_specific_evalue(es, 1, &eig1);
    if (err) goto cleanup;
    phase = "verify eigenvalues";
    if (fabs((double)eig0 - first) > 1.0e-5 ||
        fabs((double)eig1 - second) > 1.0e-5) {
        fprintf(stderr,"LIS Krylov unexpected Ritz values %g %g\n",
                (double)eig0, (double)eig1);
        goto cleanup;
    }
    printf("PASS LIS API: eigenvalues %.10g %.10g\n",
           (double)eig0, (double)eig1);
    fail = 0;
cleanup:
    if (fail) fprintf(stderr,"FAIL LIS API at %s: error=%lld status=%lld options=%s\n",phase,(long long)err,(long long)(es?es->retcode:-1),options);
    if (es) lis_esolver_destroy(es);
    if (x) lis_vector_destroy(x);
    return fail;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A = NULL,B = NULL;
    user_diag a,b;
    int rc = EXIT_FAILURE;
    LIS_INT i,err;
    err = lis_initialize(&argc,&argv);
    if (err) return EXIT_FAILURE;
    for (i=0;i<N;i++) {
        a.diagonal[i] = (double)(i + 2) * (double)(i + 1);
        b.diagonal[i] = (double)(i + 1);
    }
    err = make_matrix(&a,&A);
    if(err)goto cleanup;
    err = make_matrix(&b,&B);
    if(err)goto cleanup;

    if(solve_check(A,NULL,
         "-e krylov -ss 2 -nev 2 -etol 1e-7 -emaxiter 700 -kwhich LA -eprint none",
         156.0,132.0)) goto cleanup;
    if(solve_check(A,B,
         "-e krylov -ss 2 -nev 2 -etol 1e-7 -emaxiter 700 -shift -0.5 -kwhich LM -eprint none",
         2.0,3.0)) goto cleanup;

    rc=EXIT_SUCCESS;
cleanup:
    if(A)lis_matrix_destroy(A);
    if(B)lis_matrix_destroy(B);
    lis_finalize();
    return rc;
}
