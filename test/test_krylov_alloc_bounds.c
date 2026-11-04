/* Regression for checked Krylov workspace products on 32-bit size_t.
 * Run with gcc -m32 where the required 32-bit runtime is available.
 */
#include "lis_krylov.h"
#include <stdint.h>
#include <stdio.h>

static int apply_identity(void *context, int n,
                          const double *x, double *y)
{
    int i;
    (void)context;
    for (i = 0; i < n; ++i) y[i] = x[i];
    return 0;
}

int main(void)
{
    lis_krylov_problem problem = {0};
    lis_krylov_options options;
    lis_krylov_stats stats;
    double eigenvalue[1], eigenvector[1], residual[1];
    int rc;

    lis_krylov_options_default(&options);
    problem.n = 16;
    problem.transform = LIS_KRYLOV_STANDARD;
    problem.apply_a = apply_identity;
    options.nev = 1;
    options.ncv = 1; /* invalid even on platforms with 64-bit size_t */
    rc = lis_krylov_symmetric(&problem, &options, eigenvalue,
                              eigenvector, residual, &stats);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT) {
        fprintf(stderr, "FAIL ordinary capacity validation: %d\n", rc);
        return 1;
    }

    /* On ILP32, 131072 * 131072 = 2^34 and
     * 131072 * (131072/2) = 2^33. Both products wrap to zero in
     * 32-bit size_t. Rejection must precede all attempted allocations.
     * Do not run this huge shape with 64-bit size_t. */
    if (sizeof(size_t) == 4) {
        problem.n = 131072;
        options.nev = 1;
        options.ncv = 131072;
        rc = lis_krylov_symmetric(&problem, &options, eigenvalue,
                                  eigenvector, residual, &stats);
        if (rc != LIS_KRYLOV_BAD_ARGUMENT) {
            fprintf(stderr, "FAIL 32-bit multiplication overflow: %d\n",
                    rc);
            return 1;
        }
        puts("PASS ILP32 workspace allocation overflow rejected");
    } else {
        puts("SKIP ILP32 overflow case on 64-bit size_t");
    }
    puts("PASS Krylov allocation bounds");
    return 0;
}
