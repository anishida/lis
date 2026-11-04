/* Opt-in block-enriched restart regression.
 * Generic operator-only symmetric generalized problem, with
 * positive-semidefinite metric and a known clustered finite spectrum.
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
    int fail_solve;
} fixture;

static int apply_a(void *opaque, int n, const double *x, double *y)
{
    fixture *c = (fixture *)opaque;
    int i;
    if (n != N) return -1;
    for (i = 0; i < RANK; ++i) {
        y[i] = (c->eigenvalues[i] + 0.09)*x[i] + 0.3*x[RANK+i];
        y[RANK+i] = 0.3*x[i] + x[RANK+i];
    }
    return 0;
}

static int apply_b(void *opaque, int n, const double *x, double *y)
{
    int i;
    (void)opaque;
    if (n != N) return -1;
    for (i = 0; i < RANK; ++i) {
        y[i] = x[i];
        y[RANK+i] = 0.0;
    }
    return 0;
}

static int solve_shifted(void *opaque, int n, const double *x, double *y)
{
    fixture *c = (fixture *)opaque;
    int i;
    if (n != N || c->fail_solve) return -1;
    ++c->solves;
    for (i = 0; i < RANK; ++i) {
        double d = c->eigenvalues[i] + 0.09 - c->sigma;
        double det = d - 0.09;
        if (det == 0.0) return -1;
        y[i] = (x[i] - 0.3*x[RANK+i]) / det;
        y[RANK+i] = (d*x[RANK+i] - 0.3*x[i]) / det;
    }
    return 0;
}

static void setup(fixture *c, lis_krylov_problem *p,
                  lis_krylov_options *o, double start[N],
                  int repeated, int hardlocking)
{
    int i;
    memset(c, 0, sizeof(*c));
    c->sigma = 0.9;
    for (i = 0; i < RANK; ++i)
        c->eigenvalues[i] = i < NEV ?
                            1.0 + (repeated ? 0.0 : 1.0e-5*i) :
                            20.0 + 0.1*i;
    memset(p, 0, sizeof(*p));
    p->n = N;
    p->transform = LIS_KRYLOV_SHIFT_INVERT_B;
    p->sigma = c->sigma;
    p->ctx = c;
    p->apply_a = apply_a;
    p->apply_b = apply_b;
    p->solve = solve_shifted;
    memset(start, 0, N*sizeof(*start));
    for (i = RANK; i < N; ++i)
        start[i] = 1.0 + 0.03*i;
    lis_krylov_options_default(o);
    o->nev = NEV;
    o->ncv = 12;
    o->maxiter = 600;
    o->tol = 1e-9;
    /* Analytically safe bound for the symmetric block operator's ||A||. */
    o->residual_scale = 30.0;
    o->initial_vector = start;
    o->random_seed = 6;
    o->allow_semidefinite_metric = 1;
    o->hard_locking = hardlocking;
}

/* Same checks for all modes: values, true generalized residual,
 * B-orthogonality, normalized vectors and inverse-call accounting. */
static int check(const char *name, const fixture *c,
                 const lis_krylov_options *o,
                 const double values[NEV], const double vectors[N*NEV],
                 const double residuals[NEV], const lis_krylov_stats *s,
                 int repeated)
{
    int i, j, k;
    if (s->nconv != NEV || s->nreturned != NEV ||
        s->solve_count != c->solves) {
        fprintf(stderr, "FAIL %s conv=%d returned=%d solves=%d/%d\n",
                name, s->nconv, s->nreturned, s->solve_count, c->solves);
        return 1;
    }
    for (j = 0; j < NEV; ++j) {
        double bnorm = 0.0, norm_r = 0.0, norm_ax = 0.0;
        double norm_bx = 0.0, norm_x = 0.0, denom, relative;
        double expected = 1.0 + (repeated ? 0.0 : 1.0e-5*j);
        double ax[N], bx[N];
        if (fabs(values[j] - expected) > 1.0e-7 ||
            residuals[j] > o->tol) {
            fprintf(stderr, "FAIL %s eigenvalue %d value=%.15g "
                    "residual=%g\n", name, j, values[j], residuals[j]);
            return 1;
        }
        apply_a((void *)c, N, vectors+(size_t)j*N, ax);
        apply_b((void *)c, N, vectors+(size_t)j*N, bx);
        for (k = 0; k < N; ++k) {
            double x = vectors[(size_t)j*N+k];
            double r = ax[k] - values[j]*bx[k];
            norm_r += r*r;
            norm_ax += ax[k]*ax[k];
            norm_bx += bx[k]*bx[k];
            norm_x += x*x;
        }
        for (k = 0; k < RANK; ++k) {
            double x = vectors[(size_t)j*N+k];
            bnorm += x*x;
        }
        /* Compute the published backward error independently from the
         * returned eigenvector. By fixing residual_scale above, the
         * denominator is known and need not imitate the solver's
         * history-dependent operator-norm estimator. There is no
         * arbitrary loosening of the requested residual tolerance. */
        denom = fmax(sqrt(norm_ax) + fabs(values[j])*sqrt(norm_bx),
                     o->residual_scale * sqrt(norm_x));
        relative = denom > 0.0 ? sqrt(norm_r)/denom : sqrt(norm_r);
        if (fabs(bnorm - 1.0) > 1e-6 ||
            relative > o->tol * (1.0 + 1e-5) ||
            fabs(relative - residuals[j]) > 1e-11) {
            fprintf(stderr,
                    "FAIL %s true residual/metric norm %d "
                    "metric=%.9g relative=%.9g reported=%.9g\n",
                    name, j, bnorm, relative, residuals[j]);
            return 1;
        }
        for (i = 0; i < j; ++i) {
            double overlap = 0.0;
            for (k = 0; k < RANK; ++k)
                overlap += vectors[(size_t)i*N+k] *
                           vectors[(size_t)j*N+k];
            if (fabs(overlap) > 1e-6) {
                fprintf(stderr, "FAIL %s orthogonality %d/%d\n", name,i,j);
                return 1;
            }
        }
    }
    return 0;
}

static int test_case(int repeated, int hardlocking, int block)
{
    fixture c;
    lis_krylov_problem p;
    lis_krylov_options o;
    lis_krylov_stats s, scalar_stats;
    double initial[N], values[NEV], vectors[N*NEV], residuals[NEV];
    double scalar_values[NEV], scalar_vectors[N*NEV], scalar_res[NEV];
    int rc, i;
    setup(&c, &p, &o, initial, repeated, hardlocking);
    if (block == 1) {
        rc = lis_krylov_symmetric(&p, &o, scalar_values,
                                  scalar_vectors, scalar_res, &scalar_stats);
        if (rc != LIS_KRYLOV_SUCCESS) {
            fprintf(stderr, "FAIL scalar baseline rc=%d\n", rc);
            return 1;
        }
        c.solves = 0;
    }
    rc = lis_krylov_symmetric_block_restart(
        &p, &o, block, values, vectors, residuals, &s);
    if (rc != LIS_KRYLOV_SUCCESS ||
        check("block restart", &c, &o, values, vectors,
              residuals, &s, repeated)) {
        fprintf(stderr, "FAIL block=%d repeated=%d hardlock=%d rc=%d "
                "cycles=%d restarts=%d\n", block, repeated, hardlocking,
                rc, s.iterations, s.restarts);
        return 1;
    }
    if (block == 1) {
        if (s.iterations != scalar_stats.iterations ||
            s.restarts != scalar_stats.restarts ||
            s.solve_count != scalar_stats.solve_count ||
            s.metric_count != scalar_stats.metric_count ||
            s.op_count != scalar_stats.op_count)
            return 1;
        for (i = 0; i < NEV; ++i)
            if (values[i] != scalar_values[i] ||
                residuals[i] != scalar_res[i])
                return 1;
    }
    printf("PASS block=%d repeated=%d hardlock=%d "
           "cycles=%d restarts=%d solves=%d worst=%.3g\n",
           block, repeated, hardlocking, s.iterations,
           s.restarts, s.solve_count, s.max_residual);
    return 0;
}

static int invalid_cases(void)
{
    fixture c;
    lis_krylov_problem p;
    lis_krylov_options o;
    lis_krylov_stats s;
    double start[N], values[NEV], vectors[N*NEV], residuals[NEV];
    setup(&c, &p, &o, start, 0, 0);
    if (lis_krylov_symmetric_block_restart(&p, &o, 0, values, vectors,
          residuals, &s) != LIS_KRYLOV_BAD_ARGUMENT ||
        lis_krylov_symmetric_block_restart(&p, &o, -1, values, vectors,
          residuals, &s) != LIS_KRYLOV_BAD_ARGUMENT ||
        lis_krylov_symmetric_block_restart(&p, &o, o.ncv+1, values, vectors,
          residuals, &s) != LIS_KRYLOV_BAD_ARGUMENT) return 1;
    c.fail_solve = 1;
    if (lis_krylov_symmetric_block_restart(&p, &o, 3, values, vectors,
          residuals, &s) != LIS_KRYLOV_CALLBACK_ERROR ||
        s.nreturned != 0 || s.nconv != 0) return 1;
    puts("PASS invalid batch size and callback failure");
    return 0;
}

int main(void)
{
    return test_case(0, 0, 1) ||
           test_case(1, 1, 1) ||
           test_case(0, 0, 3) ||
           test_case(1, 0, 3) ||
           test_case(0, 1, 3) ||
           test_case(1, 1, 3) ||
           invalid_cases();
}
