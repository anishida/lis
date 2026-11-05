/* Extended Krylov API contract and allocation-limit regression.
 * Build under a 32-bit target to test pre-multiplication bounds directly.
 * No large buffer is allocated by this test.
 */
#include "lis_krylov.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int apply_diagonal(void *ctx, int n, const double *x, double *y)
{
    int i;
    (void)ctx;
    for (i = 0; i < n; ++i) y[i] = (double)(i + 1) * x[i];
    return 0;
}

int main(void)
{
    lis_krylov_problem p = {0};
    lis_krylov_options options;
    lis_krylov_stats stats;
    lis_krylov_ext_options ext;
    lis_krylov_ext_stats extra;
    lis_krylov_workspace *w = NULL;
    double values[2], vectors[16], residuals[2];
    int rc;

    lis_krylov_workspace_destroy(NULL);
    lis_krylov_options_default(&options);
    lis_krylov_ext_options_default(&ext);
    lis_krylov_ext_stats_default(&extra);
    p.n = 8;
    p.apply_a = apply_diagonal;
    options.nev = 2;
    options.ncv = 5;
    options.tol = 1e-8;
    options.maxiter = 100;

    /* The output stats structure has an independent version/size
     * check. Invalid memory layouts must be rejected before use. */
    extra.version = LIS_KRYLOV_EXT_VERSION + 1u;
    rc = lis_krylov_symmetric_ex(&p, &options, &ext, NULL,
                                 values, vectors, residuals, &stats, &extra);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT) {
        fprintf(stderr, "FAIL invalid statistics version: %d\n", rc);
        return 1;
    }
    lis_krylov_ext_stats_default(&extra);
    extra.struct_size--;
    rc = lis_krylov_symmetric_ex(&p, &options, &ext, NULL,
                                 values, vectors, residuals, &stats, &extra);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT) {
        fprintf(stderr, "FAIL invalid statistics size: %d\n", rc);
        return 1;
    }
    lis_krylov_ext_stats_default(&extra);
    ext.version++;
    w = (lis_krylov_workspace *)(uintptr_t)1;
    rc = lis_krylov_workspace_create(&p, &options, &ext, &w);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT || w != NULL) {
        fprintf(stderr, "FAIL invalid workspace version: %d\n", rc);
        return 1;
    }
    lis_krylov_ext_options_default(&ext);
    ext.struct_size++;
    rc = lis_krylov_symmetric_ex(&p, &options, &ext, NULL,
                                 values, vectors, residuals, &stats, &extra);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT) {
        fprintf(stderr, "FAIL invalid options size: %d\n", rc);
        return 1;
    }

    /* On ILP32, n * ncv = 2^34 and n * (ncv/2) = 2^33.
     * Both products wrap to zero in size_t. Previously the solver
     * computed these products BEFORE its overflow check. The corrected
     * allocator must reject this shape without trying to allocate. */
    if (sizeof(size_t) == 4) {
        p.n = 131072;
        options.nev = 1;
        options.ncv = 131072;
        lis_krylov_ext_options_default(&ext);
        w = (lis_krylov_workspace *)(uintptr_t)1;
        rc = lis_krylov_workspace_create(&p, &options, &ext, &w);
        if (rc != LIS_KRYLOV_BAD_ARGUMENT || w != NULL) {
            fprintf(stderr, "FAIL ILP32 overflow detection: %d\n", rc);
            return 1;
        }
        puts("PASS ILP32 workspace multiplication overflow rejected");
    } else {
        puts("SKIP ILP32 overflow case (64-bit target)");
    }
    puts("PASS versioned API validation and NULL workspace cleanup");
    return 0;
}
