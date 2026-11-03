/* Analytical regression for clustered and repeated generalized eigenvalues.
 * The positive-semidefinite B metric has rank 48 and a 48-dimensional
 * kernel. For each block, the finite eigenvalue is exactly l[i].
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 96
#define RANK 48
#define NEV 4

typedef struct {
    double eigenvalues[RANK];
    double sigma;
    int solves;
} test_context;

static int apply_a(void *v, int n, const double *x, double *y)
{
    test_context *c = (test_context *)v;
    int i;
    if (n != N) return 1;
    for (i = 0; i < RANK; ++i) {
        y[i] = (c->eigenvalues[i] + 0.09) * x[i] + 0.3 * x[RANK+i];
        y[RANK+i] = 0.3 * x[i] + x[RANK+i];
    }
    return 0;
}

static int apply_b(void *v, int n, const double *x, double *y)
{
    int i;
    (void)v;
    if (n != N) return 1;
    for (i = 0; i < RANK; ++i) {
        y[i] = x[i];
        y[RANK+i] = 0.0;
    }
    return 0;
}

static int solve_shifted(void *v, int n, const double *x, double *y)
{
    test_context *c = (test_context *)v;
    int i;
    if (n != N) return 1;
    ++c->solves;
    for (i = 0; i < RANK; ++i) {
        double a = c->eigenvalues[i] + 0.09 - c->sigma;
        double det = a - 0.09;
        if (det == 0.0) return 1;
        y[i] = (x[i] - 0.3 * x[RANK+i]) / det;
        y[RANK+i] = (a * x[RANK+i] - 0.3 * x[i]) / det;
    }
    return 0;
}

static int run_case(const char *name, int repeated, int ncv,
                    unsigned seed, int hard_locking)
{
    test_context c = {{0}, 0.9, 0};
    lis_krylov_problem problem = {0};
    lis_krylov_options options;
    lis_krylov_stats stats = {0};
    double values[NEV], vectors[N*NEV], residuals[NEV], start[N] = {0};
    int i, j, k, status;

    for (i = 0; i < RANK; ++i)
        c.eigenvalues[i] = i < NEV ?
                          (1.0 + (repeated ? 0.0 : 1e-5 * i)) :
                          (20.0 + 0.1 * i);
    for (i = RANK; i < N; ++i) start[i] = 1.0 + 0.03 * i;

    problem.n = N;
    problem.ctx = &c;
    problem.transform = LIS_KRYLOV_SHIFT_INVERT_B;
    problem.sigma = c.sigma;
    problem.apply_a = apply_a;
    problem.apply_b = apply_b;
    problem.solve = solve_shifted;
    lis_krylov_options_default(&options);
    options.nev = NEV;
    options.ncv = ncv;
    options.maxiter = 600;
    options.tol = 1e-9;
    options.initial_vector = start;
    options.random_seed = seed;
    options.allow_semidefinite_metric = 1;
    options.hard_locking = hard_locking;

    status = lis_krylov_symmetric(&problem, &options, values, vectors,
                                  residuals, &stats);
    if (status != LIS_KRYLOV_SUCCESS || stats.nconv != NEV ||
        stats.nreturned != NEV) {
        fprintf(stderr, "FAIL %s status=%d converged=%d returned=%d "
                "cycles=%d restarts=%d solves=%d maxres=%g\n",
                name, status, stats.nconv, stats.nreturned,
                stats.iterations, stats.restarts, stats.solve_count,
                stats.max_residual);
        return 1;
    }
    if (stats.restarts < 1 || stats.solve_count != c.solves) {
        fprintf(stderr, "FAIL %s restart/solve accounting\n", name);
        return 1;
    }
    for (i = 0; i < NEV; ++i) {
        double finite_norm = 0.0, target_energy = 0.0;
        double expected = 1.0 + (repeated ? 0.0 : 1e-5 * i);
        if (fabs(values[i] - expected) > 1e-7 ||
            residuals[i] > options.tol) {
            fprintf(stderr, "FAIL %s eigenpair %d value %.15g residual %g\n",
                    name, i, values[i], residuals[i]);
            return 1;
        }
        for (k = 0; k < RANK; ++k) {
            double z = vectors[(size_t)i*N+k];
            finite_norm += z*z;
            if (k < NEV) target_energy += z*z;
        }
        if (fabs(finite_norm - 1.0) > 1e-6 ||
            target_energy < 0.999999) {
            fprintf(stderr, "FAIL %s invariant subspace / metric norm\n", name);
            return 1;
        }
        for (j = 0; j < i; ++j) {
            double overlap = 0.0;
            for (k = 0; k < RANK; ++k)
                overlap += vectors[(size_t)i*N+k] *
                           vectors[(size_t)j*N+k];
            if (fabs(overlap) > 1e-6) {
                fprintf(stderr, "FAIL %s B-orthogonality %d/%d: %g\n",
                        name, i, j, overlap);
                return 1;
            }
        }
    }
    printf("PASS %s ncv=%d cycles=%d restarts=%d solves=%d maxres=%g\n",
           name, ncv, stats.iterations, stats.restarts,
           stats.solve_count, stats.max_residual);
    return 0;
}

int main(void)
{
    return run_case("clustered ncv=8", 0, 8, 2, 0) ||
           run_case("clustered ncv=12", 0, 12, 5, 0) ||
           run_case("multiple ncv=8", 1, 8, 1, 0) ||
           run_case("multiple ncv=12", 1, 12, 6, 0) ||
           run_case("hardlock clustered ncv=8", 0, 8, 2, 1) ||
           run_case("hardlock clustered ncv=12", 0, 12, 5, 1) ||
           run_case("hardlock multiple ncv=8", 1, 8, 1, 1) ||
           run_case("hardlock multiple ncv=12", 1, 12, 6, 1);
}
