/* Repeated-solve microbenchmark for optional workspace allocations.
 * Diagonal symmetric shift-invert problem; operator callbacks are exact.
 * Measures full solve calls, not factorization/assembly and not a
 * representative production sparse workload. Outputs trial spread.
 *
 * cc -O3 -std=c99 -Iinclude src/esolver/lis_krylov.c \
 *    test/bench_krylov_workspace.c -lm -o bench_workspace
 * ./bench_workspace [repetitions]
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 48
#define NEV 3
#define TRIALS 5

typedef struct { double sigma; int solves; } fixture;

static int apply_a(void *ptr, int n, const double *x, double *y)
{
    int i;
    (void)ptr;
    if (n != N) return -1;
    for (i = 0; i < n; ++i) y[i] = (1.0 + 0.25*i)*x[i];
    return 0;
}
static int apply_b(void *ptr, int n, const double *x, double *y)
{
    int i;
    (void)ptr;
    if (n != N) return -1;
    for (i = 0; i < n; ++i) y[i] = x[i];
    return 0;
}
static int solve(void *ptr, int n, const double *x, double *y)
{
    fixture *f = (fixture *)ptr;
    int i;
    if (n != N) return -1;
    for (i = 0; i < n; ++i)
        y[i] = x[i]/(1.0 + 0.25*i - f->sigma);
    ++f->solves;
    return 0;
}

static int compare_double(const void *left, const void *right)
{
    double a = *(const double *)left;
    double b = *(const double *)right;
    return (a > b) - (a < b);
}

static int batch(int use_workspace, int count,
                 lis_krylov_problem *problem,
                 const lis_krylov_options *options,
                 lis_krylov_workspace *w, fixture *f, double *seconds)
{
    double values[NEV], vectors[N*NEV], residuals[NEV];
    lis_krylov_stats stats;
    clock_t begin = clock();
    int i, j, rc;
    for (i = 0; i < count; ++i) {
        f->sigma = 0.05 + 0.00001*(i%500);
        problem->sigma = f->sigma;
        if (use_workspace) {
            rc = lis_krylov_symmetric_ex(problem, options, NULL, w,
                     values, vectors, residuals, &stats, NULL);
        } else {
            rc = lis_krylov_symmetric(problem, options,
                     values, vectors, residuals, &stats);
        }
        if (rc != LIS_KRYLOV_SUCCESS || stats.nconv != NEV) {
            fprintf(stderr,"FAIL %s solve %d rc=%d conv=%d\n",
                    use_workspace ? "reuse":"one-shot",i,rc,stats.nconv);
            return 1;
        }
        for (j = 0; j < NEV; ++j)
            if (fabs(values[j]-(1.0+0.25*j))>1e-7 ||
                residuals[j]>options->tol) {
                fprintf(stderr,"FAIL output j=%d value=%.9g res=%.3g\n",
                        j,values[j],residuals[j]);
                return 1;
            }
    }
    *seconds = (double)(clock()-begin)/CLOCKS_PER_SEC;
    return 0;
}
int main(int argc, char **argv)
{
    fixture f = {0};
    lis_krylov_problem p = {0};
    lis_krylov_options o;
    lis_krylov_workspace *w = NULL;
    double samples[2][TRIALS];
    int repeats = argc > 1 ? atoi(argv[1]): 200;
    int pass, t;
    if (repeats < 1 || repeats > 100000) return 2;
    p.n=N; p.transform=LIS_KRYLOV_SHIFT_INVERT_B;
    p.ctx=&f; p.apply_a=apply_a; p.apply_b=apply_b;
    p.solve=solve;
    lis_krylov_options_default(&o);
    o.nev=NEV; o.ncv=16; o.maxiter=400; o.tol=1e-8;
    if (lis_krylov_workspace_create(&p,&o,NULL,&w) ||
        !w) return 1;
    /* Warm-up both modes before recording trials. */
    { double ignored;
      if(batch(0,3,&p,&o,w,&f,&ignored) ||
         batch(1,3,&p,&o,w,&f,&ignored)) return 1;
    }
    for (t=0;t<TRIALS;++t) {
        /* Alternate order so thermal/cache drift is not one-sided. */
        for (pass=0;pass<2;++pass) {
            int mode=(pass+t)%2;
            if (batch(mode,repeats,&p,&o,w,&f,&samples[mode][t]))
                return 1;
        }
    }
    puts("mode,repeats,trials,median_cpu_s,min_cpu_s,max_cpu_s,median_us_per_solve");
    for(pass=0;pass<2;++pass) {
        qsort(samples[pass],TRIALS,sizeof(double),compare_double);
        printf("%s,%d,%d,%.6f,%.6f,%.6f,%.3f\n",
               pass?"reused":"one_shot",repeats,TRIALS,
               samples[pass][TRIALS/2],samples[pass][0],
               samples[pass][TRIALS-1],
               samples[pass][TRIALS/2]*1e6/repeats);
    }
    lis_krylov_workspace_destroy(w);
    return 0;
}
