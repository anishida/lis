/* Native thick-restarted symmetric Rayleigh-Ritz subspace eigensolver.
 * Copyright (C) 2026, contributors to LIS.
 * Redistribution and use are governed by the LIS project BSD license.
 */
#include "lis_krylov.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* A direction is considered linearly dependent only when double-pass
 * orthogonalization removes essentially all of its input metric norm. */
#define KRYLOV_ORTH_REJECT 1.0e-12
#define KRYLOV_DGKS_ETA 0.717
/* Avoid OpenMP fork/join overhead on small vectors. Fixed chunk boundaries
 * make multi-threaded reductions independent of thread scheduling. */
#define KRYLOV_OMP_MIN_N 32768
#define KRYLOV_OMP_CHUNKS 32
#define KRYLOV_DIRECTION_REJECT (128.0 * DBL_EPSILON)
#define KRYLOV_MAX_SWEEPS 80
#define KRYLOV_JACOBI_EPS (16.0 * DBL_EPSILON)

struct krylov_workspace {
    const lis_krylov_problem *problem;
    const lis_krylov_options *options;
    lis_krylov_stats *stats;
    int n;
    int ncv;
    double norm_a;
    double *basis;
    double *metric_basis;
    double *images;
    double *projected;
    double *dense_work;
    double *ritz;
    double *eval;
    double *scratch;
    double *x;
    double *metric_x;
    double *rhs;
    double *check;
    int *order;
    int *locked;
    int nlocked;
    double *locked_basis;
    double *locked_metric;
    double *locked_values;
    double *locked_residuals;
};

struct lis_krylov_workspace {
    struct krylov_workspace buffers;
    int n;
    int ncv;
    int nev;
    int transform;
    int hard_locking;
};

static double compensated_dot(int n, const double *a, const double *b)
{
#ifdef _OPENMP
    if (n >= KRYLOV_OMP_MIN_N && omp_get_max_threads() > 1) {
        double partial[KRYLOV_OMP_CHUNKS];
        double total = 0.0, correction = 0.0;
        int chunk;
        #pragma omp parallel for schedule(static)
        for (chunk = 0; chunk < KRYLOV_OMP_CHUNKS; ++chunk) {
            const int begin = (int)(((long long)n * chunk) /
                                    KRYLOV_OMP_CHUNKS);
            const int end = (int)(((long long)n * (chunk + 1)) /
                                  KRYLOV_OMP_CHUNKS);
            double sum = 0.0, correction_local = 0.0;
            int k;
            for (k = begin; k < end; ++k) {
                double item = a[k] * b[k] - correction_local;
                double next = sum + item;
                correction_local = (next - sum) - item;
                sum = next;
            }
            partial[chunk] = sum;
        }
        for (chunk = 0; chunk < KRYLOV_OMP_CHUNKS; ++chunk) {
            double item = partial[chunk] - correction;
            double next = total + item;
            correction = (next - total) - item;
            total = next;
        }
        return total;
    }
#endif
    double sum = 0.0;
    double compensation = 0.0;
    int i;
    for (i = 0; i < n; ++i) {
        double value = a[i] * b[i] - compensation;
        double next = sum + value;
        compensation = (next - sum) - value;
        sum = next;
    }
    return sum;
}

static double stable_norm(int n, const double *x)
{
#ifdef _OPENMP
    if (n >= KRYLOV_OMP_MIN_N && omp_get_max_threads() > 1) {
        double scales[KRYLOV_OMP_CHUNKS];
        double squares[KRYLOV_OMP_CHUNKS];
        double global_scale = 0.0, global_squares = 1.0;
        int chunk;
        #pragma omp parallel for schedule(static)
        for (chunk = 0; chunk < KRYLOV_OMP_CHUNKS; ++chunk) {
            const int begin = (int)(((long long)n * chunk) /
                                    KRYLOV_OMP_CHUNKS);
            const int end = (int)(((long long)n * (chunk + 1)) /
                                  KRYLOV_OMP_CHUNKS);
            double scale = 0.0, local_squares = 1.0;
            int k;
            for (k = begin; k < end; ++k) {
                const double value = fabs(x[k]);
                if (!isfinite(value)) {
                    scale = INFINITY;
                    break;
                }
                if (value == 0.0) continue;
                if (scale < value) {
                    const double ratio = scale / value;
                    local_squares = 1.0 + local_squares * ratio * ratio;
                    scale = value;
                } else {
                    const double ratio = value / scale;
                    local_squares += ratio * ratio;
                }
            }
            scales[chunk] = scale;
            squares[chunk] = local_squares;
        }
        for (chunk = 0; chunk < KRYLOV_OMP_CHUNKS; ++chunk) {
            const double scale = scales[chunk];
            if (!isfinite(scale)) return INFINITY;
            if (scale == 0.0) continue;
            if (global_scale < scale) {
                const double ratio = global_scale / scale;
                global_squares = squares[chunk] +
                                 global_squares * ratio * ratio;
                global_scale = scale;
            } else {
                const double ratio = scale / global_scale;
                global_squares += squares[chunk] * ratio * ratio;
            }
        }
        return global_scale == 0.0 ? 0.0 :
               global_scale * sqrt(global_squares);
    }
#endif
    double scale = 0.0;
    double square_sum = 1.0;
    int i;
    for (i = 0; i < n; ++i) {
        double a = fabs(x[i]);
        if (!isfinite(a)) return INFINITY;
        if (a == 0.0) continue;
        if (scale < a) {
            double ratio = scale / a;
            square_sum = 1.0 + square_sum * ratio * ratio;
            scale = a;
        } else {
            double ratio = a / scale;
            square_sum += ratio * ratio;
        }
    }
    return scale == 0.0 ? 0.0 : scale * sqrt(square_sum);
}

static int all_finite(int n, const double *x)
{
    int i;
    for (i = 0; i < n; ++i)
        if (!isfinite(x[i])) return 0;
    return 1;
}

/* Vector kernels keep the original sequential order of projections.
 * Parallel updates finish before the next metric dot product. */
static void vector_axpy(int n, double alpha, const double *y, double *x)
{
    int i;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(n >= KRYLOV_OMP_MIN_N)
#endif
    for (i = 0; i < n; ++i) x[i] += alpha * y[i];
}

static void vector_scale(int n, double scale, double *x)
{
    int i;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(n >= KRYLOV_OMP_MIN_N)
#endif
    for (i = 0; i < n; ++i) x[i] *= scale;
}

static void vector_scale_pair(int n, double scale,
                              double *x, double *metric_x)
{
    int i;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(n >= KRYLOV_OMP_MIN_N)
#endif
    for (i = 0; i < n; ++i) {
        x[i] *= scale;
        metric_x[i] *= scale;
    }
}

static int apply_metric(struct krylov_workspace *s,
                        const double *x, double *y)
{
    const lis_krylov_problem *p = s->problem;
    lis_krylov_matvec fn;
    if (p->transform == LIS_KRYLOV_STANDARD) {
        memcpy(y, x, (size_t)s->n * sizeof(double));
    } else {
        fn = p->transform == LIS_KRYLOV_SHIFT_INVERT_A ?
             p->apply_a : p->apply_b;
        ++s->stats->metric_count;
        if (fn(p->ctx, s->n, x, y)) return -1;
    }
    return all_finite(s->n, y) ? 0 : -2;
}

/* Validate an approximate shifted inverse only when explicitly requested.
 * This costs two extra matrix-vector products per solve. */
static int verify_inverse(struct krylov_workspace *s,
                          const double *rhs, const double *solution)
{
    const lis_krylov_problem *p = s->problem;
    double numerator, denominator;
    int i;

    if (s->options->inner_tol == 0.0) return 0;
    if (p->apply_a(p->ctx, s->n, solution, s->check)) return -1;
    if (p->apply_b(p->ctx, s->n, solution, s->rhs)) return -1;
    for (i = 0; i < s->n; ++i)
        s->check[i] -= p->sigma * s->rhs[i] + rhs[i];
    numerator = stable_norm(s->n, s->check);
    denominator = stable_norm(s->n, rhs);
    if (!isfinite(numerator) || !isfinite(denominator)) return -2;
    if (numerator > s->options->inner_tol * denominator)
        return -3;
    return 0;
}

/* The metric image is already computed during orthogonalization. Reusing
 * it avoids an unnecessary matrix-vector product per Krylov extension. */
static int apply_operator(struct krylov_workspace *s,
                          const double *x, const double *metric_x,
                          double *y)
{
    const lis_krylov_problem *p = s->problem;
    int err;
    ++s->stats->op_count;
    if (p->transform == LIS_KRYLOV_STANDARD) {
        err = p->apply_a(p->ctx, s->n, x, y);
    } else {
        ++s->stats->solve_count;
        err = p->solve(p->ctx, s->n, metric_x, y);
        if (!err && all_finite(s->n, y))
            err = verify_inverse(s, metric_x, y);
    }
    if (err) return err == -3 ? -3 : -1;
    return all_finite(s->n, y) ? 0 : -2;
}

static int orthogonalize(struct krylov_workspace *s,
                         int count, double *x, double *metric_x)
{
    double initial, square, norm;
    int pass, j, err;

    norm = stable_norm(s->n, x);
    if (!isfinite(norm)) return -2;
    if (norm == 0.0) return 1;
    vector_scale(s->n, 1.0 / norm, x);

    err = apply_metric(s, x, metric_x);
    if (err) return err;
    initial = compensated_dot(s->n, x, metric_x);
    if (!isfinite(initial)) return -2;
    if (initial <= 0.0) {
        if (s->options->allow_semidefinite_metric && initial == 0.0)
            return 1;
        return -2;
    }

    /* DGKS selective reorthogonalization in the chosen metric.  The
     * norm-loss test compares squared metric norms to avoid square roots.
     * Recompute B*x after each pass; reusing a heavily cancelled metric
     * image can give an inaccurate norm on semidefinite B. */
    {
        double previous_square = initial;
        for (pass = 0; pass < 3; ++pass) {
            for (j = 0; j < s->nlocked; ++j) {
                const double *v = s->locked_basis + (size_t)j * s->n;
                const double *bv = s->locked_metric + (size_t)j * s->n;
                double coefficient = compensated_dot(s->n, bv, x);
                if (!isfinite(coefficient)) return -2;
                vector_axpy(s->n, -coefficient, v, x);
            }
            for (j = 0; j < count; ++j) {
                const double *v = s->basis + (size_t)j * s->n;
                const double *bv = s->metric_basis + (size_t)j * s->n;
                double coefficient = compensated_dot(s->n, bv, x);
                if (!isfinite(coefficient)) return -2;
                vector_axpy(s->n, -coefficient, v, x);
            }
            err = apply_metric(s, x, metric_x);
            if (err) return err;
            square = compensated_dot(s->n, x, metric_x);
            if (!isfinite(square)) return -2;
            /* Tiny negative roundoff is dependence only for a PSD B. */
            if (square <= initial * KRYLOV_ORTH_REJECT *
                          KRYLOV_ORTH_REJECT) {
                if (square >= 0.0 ||
                    (s->options->allow_semidefinite_metric &&
                     square >= -64.0 * DBL_EPSILON * initial))
                    return 1;
                return -2;
            }
            if (square >= KRYLOV_DGKS_ETA * KRYLOV_DGKS_ETA *
                          previous_square)
                break;
            previous_square = square;
        }
    }
    {
        double inverse_norm = 1.0 / sqrt(square);
        vector_scale_pair(s->n, inverse_norm, x, metric_x);
    }
    for (j = 0; j < s->nlocked; ++j) {
        const double *bv = s->locked_metric + (size_t)j * s->n;
        if (fabs(compensated_dot(s->n, bv, x)) > 1e-8)
            return 1;
    }
    for (j = 0; j < count; ++j) {
        const double *bv = s->metric_basis + (size_t)j * s->n;
        if (fabs(compensated_dot(s->n, bv, x)) > 1.0e-8)
            return 1;
    }
    return 0;
}

static int append_vector(struct krylov_workspace *s, int count, double *x);

/* For a positive-semidefinite metric B, the generalized problem has
 * finite modes only in the nonzero part of the transformed spectrum.
 * Start and replenishment directions in range((A-sigma B)^-1 B) so that
 * components of ker(B) are not introduced as arbitrary Ritz directions.
 * This still requires a nonsingular shifted operator and a nonnegative
 * B-inner product on the generated Krylov subspace. */
static int append_range_seed(struct krylov_workspace *s,
                             int count, double *x)
{
    int err;
    if (s->options->allow_semidefinite_metric) {
        err = apply_metric(s, x, s->metric_x);
        if (err) return err;
        if (stable_norm(s->n, s->metric_x) == 0.0) return 1;
        /* scratch has at least n elements, and does not overlap the
         * scratch vectors used by optional inner-solve certification. */
        err = apply_operator(s, x, s->metric_x, s->scratch);
        if (err) return err;
        memcpy(x, s->scratch, (size_t)s->n * sizeof(double));
    }
    return append_vector(s, count, x);
}

static int append_vector(struct krylov_workspace *s,
                         int count, double *x)
{
    int err = orthogonalize(s, count, x, s->metric_x);
    if (err) return err;
    memcpy(s->basis + (size_t)count * s->n, x,
           (size_t)s->n * sizeof(double));
    if (s->metric_basis != s->basis)
        memcpy(s->metric_basis + (size_t)count * s->n, s->metric_x,
               (size_t)s->n * sizeof(double));
    return apply_operator(s, x, s->metric_x,
                          s->images + (size_t)count * s->n);
}

static void seed_vector(int n, double *x, unsigned seed)
{
    unsigned state = seed ^ 0x6d2b79f5u;
    int i;
    for (i = 0; i < n; ++i) {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        x[i] = (double)(state & 0xffffu) / 32768.0 - 1.0;
    }
}

/* Dense cyclic Jacobi eigendecomposition of symmetric projected matrix.
 * The transformation is stable for modest restarted subspace sizes. */
static int diagonalize(int m, double *a, double *q,
                        double *eigenvalues, int stride)
{
    int sweep, p, r, i;
    for (i = 0; i < m; ++i)
        for (p = 0; p < m; ++p)
            q[(size_t)i * stride + p] = (i == p);

    for (sweep = 0; sweep < KRYLOV_MAX_SWEEPS; ++sweep) {
        double largest = 0.0;
        double diagonal_norm = 0.0;
        double pivot_threshold;
        for (p = 0; p < m; ++p) {
            diagonal_norm += fabs(a[(size_t)p * stride + p]);
            for (r = p + 1; r < m; ++r) {
                double off = fabs(a[(size_t)p * stride + r]);
                if (off > largest) largest = off;
            }
        }
        if (largest <= KRYLOV_JACOBI_EPS *
                       (diagonal_norm + DBL_MIN)) break;

        /* Thresholded cyclic Jacobi: prioritize significant pivots in
         * the first few sweeps. Skipping tiny rotations avoids needless
         * updates of m basis and Ritz-vector entries per pair. Later
         * sweeps use the original machine-precision cutoff so clusters
         * and near-multiple eigenvalues still receive full refinement.
         * A zero threshold from sweep 3 onward prevents stagnation. */
        pivot_threshold = sweep < 3 ? 0.05 * largest : 0.0;
        for (p = 0; p < m; ++p) {
            for (r = p + 1; r < m; ++r) {
                double apr = a[(size_t)p * stride + r];
                double app = a[(size_t)p * stride + p];
                double arr = a[(size_t)r * stride + r];
                double tau, angle, cosine, sine;
                if (fabs(apr) <= pivot_threshold ||
                    fabs(apr) <= DBL_EPSILON *
                                 (fabs(app) + fabs(arr) + DBL_MIN))
                    continue;
                tau = (arr - app) / (2.0 * apr);
                angle = copysign(1.0, tau) /
                        (fabs(tau) + hypot(1.0, tau));
                cosine = 1.0 / sqrt(1.0 + angle * angle);
                sine = angle * cosine;
                a[(size_t)p * stride + p] = app - angle * apr;
                a[(size_t)r * stride + r] = arr + angle * apr;
                a[(size_t)p * stride + r] = 0.0;
                a[(size_t)r * stride + p] = 0.0;
                for (i = 0; i < m; ++i) {
                    double vip, vir;
                    if (i != p && i != r) {
                        vip = a[(size_t)i * stride + p];
                        vir = a[(size_t)i * stride + r];
                        a[(size_t)i * stride + p] =
                            a[(size_t)p * stride + i] =
                            cosine * vip - sine * vir;
                        a[(size_t)i * stride + r] =
                            a[(size_t)r * stride + i] =
                            sine * vip + cosine * vir;
                    }
                    vip = q[(size_t)i * stride + p];
                    vir = q[(size_t)i * stride + r];
                    q[(size_t)i * stride + p] =
                        cosine * vip - sine * vir;
                    q[(size_t)i * stride + r] =
                        sine * vip + cosine * vir;
                }
            }
        }
    }
    if (sweep == KRYLOV_MAX_SWEEPS) return -1;
    for (i = 0; i < m; ++i)
        eigenvalues[i] = a[(size_t)i * stride + i];
    return 0;
}

static int earlier(double a, double b, const char *which)
{
    if (which[0] == 'L' && which[1] == 'M') return fabs(a) > fabs(b);
    if (which[0] == 'S' && which[1] == 'M') return fabs(a) < fabs(b);
    if (which[0] == 'L' && which[1] == 'A') return a > b;
    return a < b;
}

static void sort_ritz(int m, const double *values, int *order,
                      const char *which)
{
    int i, j;
    for (i = 0; i < m; ++i) order[i] = i;
    for (i = 1; i < m; ++i) {
        int key = order[i];
        j = i - 1;
        while (j >= 0 && earlier(values[key], values[order[j]], which)) {
            order[j + 1] = order[j];
            --j;
        }
        order[j + 1] = key;
    }
}

static void rotate_column(struct krylov_workspace *s,
                          const double *basis, int m,
                          int idx, double *output)
{
    int i, j;
#ifdef _OPENMP
    if (s->n >= KRYLOV_OMP_MIN_N && omp_get_max_threads() > 1) {
        #pragma omp parallel for schedule(static)
        for (i = 0; i < s->n; ++i) {
            double value = 0.0;
            int col;
            for (col = 0; col < m; ++col) {
                const double coeff = s->ritz[(size_t)col * s->ncv + idx];
                value += coeff * basis[(size_t)col * s->n + i];
            }
            output[i] = value;
        }
        return;
    }
#endif
    for (i = 0; i < s->n; ++i) output[i] = 0.0;
    for (j = 0; j < m; ++j) {
        double c = s->ritz[(size_t)j * s->ncv + idx];
        const double *col = basis + (size_t)j * s->n;
        for (i = 0; i < s->n; ++i) output[i] += c * col[i];
    }
}

static int eigenpair_residual(struct krylov_workspace *s,
                              int m, int index, double mu,
                              double *lambda, double *relative)
{
    const lis_krylov_problem *p = s->problem;
    double norm_a, norm_b, norm_r, norm_x, denom;
 
    rotate_column(s, s->basis, m, index, s->x);
    if (p->transform == LIS_KRYLOV_STANDARD) {
        *lambda = mu;
    } else if (p->transform == LIS_KRYLOV_SHIFT_INVERT_B) {
        if (fabs(mu) < DBL_MIN) return 1;
        *lambda = p->sigma + 1.0 / mu;
    } else {
        /* A Ritz pole too close to 1 cannot be mapped reliably:
         * the corresponding generalized value may be infinite. */
        if (fabs(mu - 1.0) <= 64.0 * DBL_EPSILON *
                                fmax(1.0, fabs(mu))) return 1;
        *lambda = p->sigma * mu / (mu - 1.0);
    }
    if (!isfinite(*lambda)) return 1;

    if (p->apply_a(p->ctx, s->n, s->x, s->metric_x)) return -1;
    norm_a = stable_norm(s->n, s->metric_x);
    if (p->transform == LIS_KRYLOV_STANDARD) {
        memcpy(s->rhs, s->x, (size_t)s->n * sizeof(double));
    } else if (p->apply_b(p->ctx, s->n, s->x, s->rhs)) {
        return -1;
    }
    norm_b = stable_norm(s->n, s->rhs);
    norm_x = stable_norm(s->n, s->x);
    if (!isfinite(norm_a) || !isfinite(norm_b) ||
        !isfinite(norm_x)) return -2;

    if (s->options->residual_scale == 0.0 && norm_x > 0.0 &&
        norm_a / norm_x > s->norm_a)
        s->norm_a = norm_a / norm_x;

    vector_axpy(s->n, -*lambda, s->rhs, s->metric_x);
    norm_r = stable_norm(s->n, s->metric_x);
    denom = norm_a + fabs(*lambda) * norm_b;
    if (s->norm_a * norm_x > denom)
        denom = s->norm_a * norm_x;
    if (!isfinite(denom) || !isfinite(norm_r)) return -2;
    *relative = denom > 0.0 ? norm_r / denom : norm_r;
    return isfinite(*relative) ? 0 : -2;
}

/* Update the projected matrix only for new columns. Projection entries
 * for retained Ritz vectors are restored on restart from the Ritz values. */
static void update_projected(struct krylov_workspace *s,
                             int previous, int m)
{
    int i, j;
    for (j = previous; j < m; ++j) {
        for (i = 0; i <= j; ++i) {
            double left = compensated_dot(s->n,
                               s->metric_basis + (size_t)i * s->n,
                               s->images + (size_t)j * s->n);
            double right = compensated_dot(s->n,
                               s->metric_basis + (size_t)j * s->n,
                               s->images + (size_t)i * s->n);
            double a = 0.5 * (left + right);
            s->projected[(size_t)i * s->ncv + j] = a;
            s->projected[(size_t)j * s->ncv + i] = a;
        }
    }
}

static void restart_basis(struct krylov_workspace *s,
                          int m, int retain)
{
    double *sets[3];
    int count = 0, k, j, i;
    sets[count++] = s->basis;
    if (s->metric_basis != s->basis)
        sets[count++] = s->metric_basis;
    sets[count++] = s->images;
    for (k = 0; k < count; ++k) {
        for (j = 0; j < retain; ++j)
            rotate_column(s, sets[k], m, s->order[j],
                          s->scratch + (size_t)j * s->n);
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) \
            if((size_t)retain * (size_t)s->n >= KRYLOV_OMP_MIN_N)
#endif
        for (j = 0; j < retain; ++j)
            memcpy(sets[k] + (size_t)j * s->n,
                   s->scratch + (size_t)j * s->n,
                   (size_t)s->n * sizeof(double));
    }
    for (i = 0; i < retain; ++i) {
        for (j = 0; j < retain; ++j)
            s->projected[(size_t)i * s->ncv + j] =
                (i == j ? s->eval[s->order[i]] : 0.0);
    }
}

void lis_krylov_options_default(lis_krylov_options *o)
{
    if (!o) return;
    memset(o, 0, sizeof(*o));
    o->tol = 1.0e-10;
    o->maxiter = 300;
    o->which = "LM";
}

static void release_workspace(struct krylov_workspace *s)
{
    free(s->basis);
    if (s->metric_basis != s->basis) free(s->metric_basis);
    free(s->images);
    free(s->projected);
    free(s->dense_work);
    free(s->ritz);
    free(s->eval);
    free(s->scratch);
    free(s->x);
    free(s->metric_x);
    free(s->rhs);
    free(s->check);
    free(s->order);
    free(s->locked);
    free(s->locked_basis);
    if (s->locked_metric != s->locked_basis) free(s->locked_metric);
    free(s->locked_values);
    free(s->locked_residuals);
}


/* All allocation is centralized so a one-shot solve and a reusable
 * workspace use precisely the same storage and alignment assumptions. */
static int allocate_workspace_buffers(struct krylov_workspace *s,
                                      int n, int ncv, int nev,
                                      int transform, int hard_locking)
{
    size_t basis_size, dense_size, scratch_size;
    if (!s || n < 3 || nev < 1 || nev >= n || ncv <= nev ||
        ncv > n || n > INT_MAX / 8 || ncv > INT_MAX / 8)
        return LIS_KRYLOV_BAD_ARGUMENT;
    basis_size = (size_t)n * (size_t)ncv;
    dense_size = (size_t)ncv * (size_t)ncv;
    scratch_size = (size_t)n *
                   (size_t)(nev > ncv / 2 ? nev : ncv / 2);
    if (basis_size > SIZE_MAX / sizeof(double) ||
        dense_size > SIZE_MAX / sizeof(double) ||
        scratch_size > SIZE_MAX / sizeof(double))
        return LIS_KRYLOV_BAD_ARGUMENT;
    s->n = n;
    s->ncv = ncv;
#define ALLOC(PTR, COUNT, TYPE)                                     \
    do {                                                            \
        (PTR) = (TYPE *)calloc((COUNT), sizeof(TYPE));              \
        if (!(PTR)) goto failed;                                    \
    } while (0)
    ALLOC(s->basis, basis_size, double);
    if (transform == LIS_KRYLOV_STANDARD)
        s->metric_basis = s->basis;
    else
        ALLOC(s->metric_basis, basis_size, double);
    ALLOC(s->images, basis_size, double);
    ALLOC(s->projected, dense_size, double);
    ALLOC(s->dense_work, dense_size, double);
    ALLOC(s->ritz, dense_size, double);
    ALLOC(s->eval, ncv, double);
    ALLOC(s->scratch, scratch_size, double);
    ALLOC(s->x, n, double);
    ALLOC(s->metric_x, n, double);
    ALLOC(s->rhs, n, double);
    ALLOC(s->check, n, double);
    ALLOC(s->order, ncv, int);
    ALLOC(s->locked, ncv, int);
    if (hard_locking) {
        ALLOC(s->locked_basis, (size_t)nev * n, double);
        if (transform == LIS_KRYLOV_STANDARD)
            s->locked_metric = s->locked_basis;
        else
            ALLOC(s->locked_metric, (size_t)nev * n, double);
        ALLOC(s->locked_values, nev, double);
        ALLOC(s->locked_residuals, nev, double);
    }
#undef ALLOC
    return LIS_KRYLOV_SUCCESS;
failed:
#undef ALLOC
    release_workspace(s);
    memset(s, 0, sizeof(*s));
    return LIS_KRYLOV_NO_MEMORY;
}

int lis_krylov_workspace_create(const lis_krylov_problem *p,
                                const lis_krylov_options *o,
                                lis_krylov_workspace **result)
{
    lis_krylov_workspace *w;
    int ncv, rc;
    if (result) *result = NULL;
    if (!p || !o || !result || p->n < 3 ||
        o->nev < 1 || o->nev >= p->n || o->ncv < 0 ||
        o->maxiter < 0 || o->hard_locking < 0 || o->hard_locking > 1 ||
        p->transform < LIS_KRYLOV_STANDARD ||
        p->transform > LIS_KRYLOV_SHIFT_INVERT_A)
        return LIS_KRYLOV_BAD_ARGUMENT;
    ncv = o->ncv ? o->ncv :
          (o->nev < (p->n - 8) / 3 ? 3 * o->nev + 8 : p->n);
    if (ncv <= o->nev || ncv > p->n || p->n > INT_MAX / 8 ||
        ncv > INT_MAX / 8)
        return LIS_KRYLOV_BAD_ARGUMENT;
    w = (lis_krylov_workspace *)calloc(1, sizeof(*w));
    if (!w) return LIS_KRYLOV_NO_MEMORY;
    rc = allocate_workspace_buffers(&w->buffers, p->n, ncv, o->nev,
                                    p->transform, o->hard_locking);
    if (rc) { free(w); return rc; }
    w->n = p->n;
    w->ncv = ncv;
    w->nev = o->nev;
    w->transform = p->transform;
    w->hard_locking = o->hard_locking;
    *result = w;
    return LIS_KRYLOV_SUCCESS;
}

void lis_krylov_workspace_destroy(lis_krylov_workspace *workspace)
{
    if (!workspace) return;
    release_workspace(&workspace->buffers);
    free(workspace);
}

static int lis_krylov_symmetric_core(lis_krylov_workspace *reuse,
                         const lis_krylov_problem *p,
                         const lis_krylov_options *o,
                         double *values, double *vectors, double *residuals,
                         lis_krylov_stats *stats)
{
    struct krylov_workspace s = {0};
    size_t dense_size;
    int m = 0, previous = 0, iter, ncv, nev, maxiter, rc;
    int i, j;
    double tol;
    unsigned seed;
    const char *which;

    if (stats) memset(stats, 0, sizeof(*stats));
    if (!p || !o || !values || !vectors || !residuals || !stats ||
        p->n < 3 || o->nev < 1 || o->nev >= p->n || !p->apply_a ||
        p->transform < LIS_KRYLOV_STANDARD ||
        p->transform > LIS_KRYLOV_SHIFT_INVERT_A ||
        (p->transform != LIS_KRYLOV_STANDARD &&
         (!p->apply_b || !p->solve)) ||
        (p->transform == LIS_KRYLOV_SHIFT_INVERT_A && p->sigma == 0.0) ||
        !isfinite(p->sigma) || !isfinite(o->tol) || o->tol < 0.0 ||
        (o->hard_locking != 0 && o->hard_locking != 1) ||
        (o->allow_semidefinite_metric != 0 &&
         o->allow_semidefinite_metric != 1) ||
        (o->allow_semidefinite_metric &&
         p->transform != LIS_KRYLOV_SHIFT_INVERT_B) ||
        !isfinite(o->inner_tol) || o->inner_tol < 0.0 ||
        !isfinite(o->residual_scale) || o->residual_scale < 0.0 ||
        o->maxiter < 0)
        return LIS_KRYLOV_BAD_ARGUMENT;
    which = o->which ? o->which : "LM";
    if (strcmp(which, "LM") && strcmp(which, "SM") &&
        strcmp(which, "LA") && strcmp(which, "SA"))
        return LIS_KRYLOV_BAD_ARGUMENT;

    nev = o->nev;
    maxiter = o->maxiter ? o->maxiter : 300;
    tol = o->tol ? o->tol : 1.0e-10;
    ncv = o->ncv ? o->ncv :
          (nev < (p->n - 8) / 3 ? 3 * nev + 8 : p->n);
    if (ncv <= nev || ncv > p->n || p->n > INT_MAX / 8 ||
        ncv > INT_MAX / 8) return LIS_KRYLOV_BAD_ARGUMENT;
    dense_size = (size_t)ncv * (size_t)ncv;
    if (reuse) {
        if (reuse->n != p->n || reuse->ncv != ncv ||
            reuse->nev != nev || reuse->transform != p->transform ||
            reuse->hard_locking != o->hard_locking)
            return LIS_KRYLOV_BAD_ARGUMENT;
        s = reuse->buffers;
    } else {
        rc = allocate_workspace_buffers(&s, p->n, ncv, nev,
                                        p->transform, o->hard_locking);
        if (rc != LIS_KRYLOV_SUCCESS) return rc;
    }
    s.problem = p;
    s.options = o;
    s.stats = stats;
    s.norm_a = o->residual_scale;
    s.nlocked = 0;
    seed = o->random_seed ? o->random_seed : 1u;
    rc = LIS_KRYLOV_NO_CONVERGENCE;

    if (o->initial_vector) {
        memcpy(s.x, o->initial_vector, (size_t)s.n * sizeof(double));
    } else {
        seed_vector(s.n, s.x, seed);
    }
    if (o->residual_scale == 0.0) {
        double norm_x = stable_norm(s.n, s.x);
        double norm_a;
        if (!isfinite(norm_x) || norm_x == 0.0) {
            rc = LIS_KRYLOV_BAD_ARGUMENT;
            goto done;
        }
        if (p->apply_a(p->ctx, s.n, s.x, s.rhs)) {
            rc = LIS_KRYLOV_CALLBACK_ERROR;
            goto done;
        }
        norm_a = stable_norm(s.n, s.rhs);
        if (!isfinite(norm_a)) {
            rc = LIS_KRYLOV_BREAKDOWN;
            goto done;
        }
        s.norm_a = norm_a / norm_x;
    }

    /* In semidefinite mode an explicit start may lie in ker(B).  Try
     * reproducible range-projected replacement directions instead of
     * reporting spurious breakdown or dividing by a zero metric norm. */
    {
        int attempt;
        for (attempt = 0; attempt <
             (o->allow_semidefinite_metric ? 16 : 1); ++attempt) {
            if (attempt > 0)
                seed_vector(s.n, s.x, seed + (unsigned)attempt);
            i = append_range_seed(&s, m, s.x);
            if (i == 0) break;
            if (i < 0) {
                rc = i == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                     i == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                     LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
        }
        if (i != 0) {
            rc = o->allow_semidefinite_metric ?
                 LIS_KRYLOV_NO_CONVERGENCE : LIS_KRYLOV_BREAKDOWN;
            goto done;
        }
    }
    m = 1;

    for (iter = 0; iter < maxiter; ++iter) {
        int old_m = m, extensions = 0;
        stats->iterations = iter + 1;
        update_projected(&s, previous, m);
        memcpy(s.dense_work, s.projected, dense_size * sizeof(double));
        if (diagonalize(m, s.dense_work, s.ritz, s.eval, ncv)) {
            rc = LIS_KRYLOV_BREAKDOWN;
            goto done;
        }
        sort_ritz(m, s.eval, s.order, which);
        if (o->hard_locking && s.nlocked < nev) {
            int captured = 0;
            /* Freeze at most one verified pair per cycle. New directions
             * are deflated against all frozen metric vectors. */
            for (j = 0; j < m && j < nev - s.nlocked; ++j) {
                int idx = s.order[j];
                double lambda = NAN, rel = INFINITY;
                int err = eigenpair_residual(&s, m, idx, s.eval[idx],
                                             &lambda, &rel);
                if (err < 0) {
                    rc = err == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                         LIS_KRYLOV_BREAKDOWN;
                    goto done;
                }
                if (err != 0 || rel > tol) continue;
                err = apply_metric(&s, s.x, s.metric_x);
                if (err) {
                    rc = err == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                         LIS_KRYLOV_BREAKDOWN;
                    goto done;
                }
                memcpy(s.locked_basis + (size_t)s.nlocked * s.n, s.x,
                       (size_t)s.n * sizeof(double));
                if (s.locked_metric != s.locked_basis)
                    memcpy(s.locked_metric + (size_t)s.nlocked * s.n,
                           s.metric_x, (size_t)s.n * sizeof(double));
                s.locked_values[s.nlocked] = lambda;
                s.locked_residuals[s.nlocked] = rel;
                ++s.nlocked;
                captured = 1;
                break;
            }
            if (captured) {
                int trial, ok = 0;
                stats->nreturned = stats->nconv = s.nlocked;
                stats->max_residual = 0.0;
                for (j = 0; j < s.nlocked; ++j) {
                    values[j] = s.locked_values[j];
                    residuals[j] = s.locked_residuals[j];
                    if (residuals[j] > stats->max_residual)
                        stats->max_residual = residuals[j];
                    memcpy(vectors + (size_t)j * s.n,
                           s.locked_basis + (size_t)j * s.n,
                           (size_t)s.n * sizeof(double));
                }
                if (s.nlocked == nev) {
                    rc = LIS_KRYLOV_SUCCESS;
                    goto done;
                }
                /* Purge the active Ritz subspace and restart with a new
                 * direction orthogonal to the permanently locked set. */
                m = previous = 0;
                for (trial = 0; trial < 32; ++trial) {
                    int attempt;
                    seed_vector(s.n, s.x,
                                seed + (unsigned)iter * 577u +
                                (unsigned)trial + 1001u);
                    attempt = append_range_seed(&s, 0, s.x);
                    if (attempt < 0) {
                        rc = attempt == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                             attempt == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                             LIS_KRYLOV_BREAKDOWN;
                        goto done;
                    }
                    if (attempt == 0) { m = 1; ok = 1; break; }
                }
                if (!ok) {
                    rc = LIS_KRYLOV_NO_CONVERGENCE;
                    goto done;
                }
                ++stats->restarts;
                continue;
            }
        }
        stats->nreturned = s.nlocked;
        stats->nconv = s.nlocked;
        stats->max_residual = 0.0;
        for (j = 0; j < s.nlocked; ++j) {
            values[j] = s.locked_values[j];
            residuals[j] = s.locked_residuals[j];
            memcpy(vectors + (size_t)j * s.n,
                   s.locked_basis + (size_t)j * s.n,
                   (size_t)s.n * sizeof(double));
            if (residuals[j] > stats->max_residual)
                stats->max_residual = residuals[j];
        }
        memset(s.locked, 0, (size_t)m * sizeof(int));
        for (j = 0; j < m && stats->nreturned < nev; ++j) {
            int idx = s.order[j];
            double lambda = NAN, rel = INFINITY;
            int err = eigenpair_residual(&s, m, idx, s.eval[idx],
                                         &lambda, &rel);
            if (err < 0) {
                rc = err == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                     LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
            /* Do not expose a nonrepresentable original eigenvalue
             * as a valid partial result. */
            if (err == 1) continue;
            s.locked[j] = (rel <= tol);
            if (s.locked[j]) ++stats->nconv;
            if (rel > stats->max_residual) stats->max_residual = rel;
            values[stats->nreturned] = lambda;
            residuals[stats->nreturned] = rel;
            memcpy(vectors + (size_t)stats->nreturned * s.n, s.x,
                   (size_t)s.n * sizeof(double));
            ++stats->nreturned;
        }
        if (stats->nconv == nev) {
            rc = LIS_KRYLOV_SUCCESS;
            goto done;
        }
        if (m >= ncv) {
            int keep = m / 2;
            if (keep < nev) keep = nev;
            if (keep >= m) keep = m - 1;
            restart_basis(&s, m, keep);
            ++stats->restarts;
            m = keep;
            previous = keep;
            /* A degeneracy cannot be resolved by residual expansion alone
             * when the retained subspace already approximates an invariant
             * eigenspace. Inject a reproducible independent search direction
             * after restarting to discover other vectors of a cluster. */
            for (i = 0; i < 4 && m < ncv; ++i) {
                int result;
                seed_vector(s.n, s.x,
                            seed + (unsigned)iter * 997u + (unsigned)i + 17u);
                result = append_range_seed(&s, m, s.x);
                if (result < 0) {
                    rc = result == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                         result == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                         LIS_KRYLOV_BREAKDOWN;
                    goto done;
                }
                if (result == 0) { ++m; break; }
            }
            continue;
        }

        for (j = 0; j < old_m && j < nev - s.nlocked && m < ncv; ++j) {
            int idx = s.order[j];
            double mu = s.eval[idx];
            double norm_w, norm_v, norm_r, threshold;
            if (s.locked[j]) continue;
            rotate_column(&s, s.images, old_m, idx, s.x);
            rotate_column(&s, s.basis, old_m, idx, s.metric_x);
            norm_w = stable_norm(s.n, s.x);
            norm_v = stable_norm(s.n, s.metric_x);
            for (i = 0; i < s.n; ++i)
                s.x[i] -= mu * s.metric_x[i];
            norm_r = stable_norm(s.n, s.x);
            threshold = KRYLOV_DIRECTION_REJECT *
                        (norm_w + fabs(mu) * norm_v);
            if (norm_r <= threshold) continue;
            i = append_vector(&s, m, s.x);
            if (i < 0) {
                rc = i == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                     i == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                     LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
            if (i == 0) { ++m; ++extensions; }
        }
        previous = old_m;
        if (!extensions && m < ncv) {
            int trial;
            for (trial = 0; trial < 6 && !extensions; ++trial) {
                seed_vector(s.n, s.x,
                            seed + (unsigned)iter * 7u + (unsigned)trial + 1u);
                i = append_range_seed(&s, m, s.x);
                if (i < 0) {
                    rc = i == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                         i == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                         LIS_KRYLOV_BREAKDOWN;
                    goto done;
                }
                if (i == 0) { ++m; ++extensions; }
            }
            if (!extensions) {
                /* A semidefinite metric can have rank smaller than nev.
                 * Preserve finite Ritz approximations and report that
                 * the requested number of pairs was not obtained. */
                rc = o->allow_semidefinite_metric ?
                     LIS_KRYLOV_NO_CONVERGENCE : LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
        }
    }
    rc = LIS_KRYLOV_NO_CONVERGENCE;

done:
    /* Locked pairs are accumulated across restarts. Restore the spectral
     * ordering of the returned original eigenpairs before publishing them. */
    if (o && o->hard_locking && s.scratch &&
        (rc == LIS_KRYLOV_SUCCESS || rc == LIS_KRYLOV_NO_CONVERGENCE)) {
        for (i = 1; i < stats->nreturned; ++i) {
            for (j = i; j > 0; --j) {
                double left = values[j], right = values[j-1];
                double mu_left, mu_right, tmp;
                if (p->transform == LIS_KRYLOV_STANDARD) {
                    mu_left = left;
                    mu_right = right;
                } else if (p->transform == LIS_KRYLOV_SHIFT_INVERT_B) {
                    mu_left = 1.0 / (left - p->sigma);
                    mu_right = 1.0 / (right - p->sigma);
                } else {
                    mu_left = left / (left - p->sigma);
                    mu_right = right / (right - p->sigma);
                }
                if (!earlier(mu_left, mu_right, which)) break;
                tmp = values[j-1];
                values[j-1] = values[j];
                values[j] = tmp;
                tmp = residuals[j-1];
                residuals[j-1] = residuals[j];
                residuals[j] = tmp;
                memcpy(s.scratch, vectors + (size_t)(j-1) * s.n,
                       (size_t)s.n * sizeof(double));
                memcpy(vectors + (size_t)(j-1) * s.n,
                       vectors + (size_t)j * s.n,
                       (size_t)s.n * sizeof(double));
                memcpy(vectors + (size_t)j * s.n, s.scratch,
                       (size_t)s.n * sizeof(double));
            }
        }
    }
    if (rc != LIS_KRYLOV_SUCCESS && rc != LIS_KRYLOV_NO_CONVERGENCE) {
        stats->nreturned = 0;
        stats->nconv = 0;
    }
    if (!reuse) release_workspace(&s);
    return rc;
}

int lis_krylov_symmetric(const lis_krylov_problem *p,
                         const lis_krylov_options *o,
                         double *values, double *vectors, double *residuals,
                         lis_krylov_stats *stats)
{
    return lis_krylov_symmetric_core(NULL, p, o, values, vectors,
                                     residuals, stats);
}

int lis_krylov_symmetric_reuse(lis_krylov_workspace *workspace,
                               const lis_krylov_problem *p,
                               const lis_krylov_options *o,
                               double *values, double *vectors,
                               double *residuals, lis_krylov_stats *stats)
{
    if (!workspace) return LIS_KRYLOV_BAD_ARGUMENT;
    return lis_krylov_symmetric_core(workspace, p, o, values, vectors,
                                     residuals, stats);
}

/* For SPD B, Sylvester's inertia law converts symmetric indefinite
 * factorizations of A-shift*B into an eigenvalue count. The callback owns
 * all factorization and pivot-sign decisions; the eigensolver has no direct
 * factorization dependency. The interval is open and its endpoints must
 * not coincide with eigenvalues (i.e. no zero inertia pivots). */
int lis_krylov_count_interval(int n, void *ctx,
                              lis_krylov_inertia_callback inertia,
                              double lower, double upper, int *count)
{
    int neg_lo = 0, neg_hi = 0, zero_lo = 0, zero_hi = 0;
    if (count) *count = 0;
    if (n <= 0 || !inertia || !count ||
        !isfinite(lower) || !isfinite(upper) || !(lower < upper))
        return LIS_KRYLOV_BAD_ARGUMENT;
    if (inertia(ctx, n, lower, &neg_lo, &zero_lo) ||
        inertia(ctx, n, upper, &neg_hi, &zero_hi))
        return LIS_KRYLOV_CALLBACK_ERROR;
    if (neg_lo < 0 || neg_lo > n || neg_hi < 0 || neg_hi > n ||
        zero_lo < 0 || zero_lo > n - neg_lo ||
        zero_hi < 0 || zero_hi > n - neg_hi ||
        neg_hi < neg_lo)
        return LIS_KRYLOV_BREAKDOWN;
    if (zero_lo || zero_hi) return LIS_KRYLOV_SINGULAR_SHIFT;
    *count = neg_hi - neg_lo;
    return LIS_KRYLOV_SUCCESS;
}

/* An independently counted symmetric window with a midpoint shift.
 * For SPD B, transformed eigenvalues inside the symmetric open interval
 * outrank those outside by magnitude. All requested pairs must converge
 * and lie strictly inside the window before success can be claimed. */
int lis_krylov_solve_interval(const lis_krylov_problem *p,
                              const lis_krylov_options *options,
                              lis_krylov_inertia_callback inertia,
                              double lower, double upper,
                              double *values, double *vectors,
                              double *residuals, lis_krylov_stats *stats,
                              int *expected_count)
{
    lis_krylov_options o;
    double center;
    int count = 0, rc, i;
    if (stats) memset(stats, 0, sizeof(*stats));
    if (expected_count) *expected_count = 0;
    if (!p || !options || !inertia || !stats || !expected_count ||
        !values || !vectors || !residuals ||
        p->transform != LIS_KRYLOV_SHIFT_INVERT_B ||
        options->allow_semidefinite_metric ||
        !isfinite(lower) || !isfinite(upper) || !(lower < upper))
        return LIS_KRYLOV_BAD_ARGUMENT;
    center = lower + 0.5 * (upper - lower);
    if (!isfinite(center) || p->sigma != center)
        return LIS_KRYLOV_BAD_ARGUMENT;
    rc = lis_krylov_count_interval(p->n, p->ctx, inertia,
                                   lower, upper, &count);
    if (rc != LIS_KRYLOV_SUCCESS) return rc;
    *expected_count = count;
    if (!count) return LIS_KRYLOV_SUCCESS;
    if (options->nev < count || count >= p->n ||
        (options->ncv && options->ncv <= count))
        return LIS_KRYLOV_BAD_ARGUMENT;
    o = *options;
    o.nev = count;
    o.which = "LM";
    rc = lis_krylov_symmetric(p, &o, values, vectors, residuals, stats);
    if (rc != LIS_KRYLOV_SUCCESS) return rc;
    if (stats->nreturned != count || stats->nconv != count)
        return LIS_KRYLOV_NO_CONVERGENCE;
    for (i = 0; i < count; ++i)
        if (!(values[i] > lower && values[i] < upper))
            return LIS_KRYLOV_NO_CONVERGENCE;
    return LIS_KRYLOV_SUCCESS;
}
