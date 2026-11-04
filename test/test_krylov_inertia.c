/* Inertia-based interval counting with independent analytical reference. */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>

typedef struct {
    double eigenvalues[8];
    double metric[8];
    int calls;
    int fail;
    int invalid;
} inertia_context;

static int diagonal_inertia(void *v, int n, double shift,
                            int *negative, int *zero)
{
    inertia_context *c = (inertia_context *)v;
    int i;
    if (n != 8) return 1;
    ++c->calls;
    if (c->fail) return 1;
    *negative = 0;
    *zero = 0;
    for (i = 0; i < n; ++i) {
        double diagonal = c->metric[i] * (c->eigenvalues[i] - shift);
        if (diagonal < 0.0) ++*negative;
        else if (diagonal == 0.0) ++*zero;
    }
    if (c->invalid) *negative += n;
    return 0;
}

static int check_interval(inertia_context *c, double lower,
                          double upper, int expected)
{
    int result = -1;
    int rc = lis_krylov_count_interval(8, c, diagonal_inertia,
                                       lower, upper, &result);
    if (rc != LIS_KRYLOV_SUCCESS || result != expected) {
        fprintf(stderr,
                "FAIL interval (%g,%g): status=%d count=%d expected=%d\n",
                lower, upper, rc, result, expected);
        return 1;
    }
    printf("PASS interval (%g,%g): %d eigenvalues\n",
           lower, upper, result);
    return 0;
}

int main(void)
{
    inertia_context c = {
        {-4.0, -2.0, 0.0, 1.0, 1.0, 3.0, 6.0, 9.0},
        {1.0, 2.0, 0.5, 3.0, 1.0, 4.0, 1.0, 0.25},
        0, 0, 0
    };
    int count, rc;

    if (check_interval(&c, -1.5, 2.5, 3) ||
        check_interval(&c, -3.0, 4.0, 5) ||
        check_interval(&c, -20.0, 20.0, 8) ||
        check_interval(&c, 10.0, 11.0, 0))
        return 1;

    rc = lis_krylov_count_interval(8, &c, diagonal_inertia,
                                   0.0, 2.0, &count);
    if (rc != LIS_KRYLOV_SINGULAR_SHIFT) return 1;
    rc = lis_krylov_count_interval(8, &c, diagonal_inertia,
                                   2.0, 2.0, &count);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT) return 1;
    rc = lis_krylov_count_interval(8, &c, diagonal_inertia,
                                   NAN, 2.0, &count);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT) return 1;

    c.fail = 1;
    rc = lis_krylov_count_interval(8, &c, diagonal_inertia,
                                   -1.0, 2.0, &count);
    if (rc != LIS_KRYLOV_CALLBACK_ERROR) return 1;
    c.fail = 0;
    c.invalid = 1;
    rc = lis_krylov_count_interval(8, &c, diagonal_inertia,
                                   -1.0, 2.0, &count);
    if (rc != LIS_KRYLOV_BREAKDOWN) return 1;

    puts("PASS singular boundaries, bad arguments and callback errors");
    return 0;
}
