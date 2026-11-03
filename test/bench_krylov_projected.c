/* Standalone microbenchmark of the internal dense Jacobi projection.
 * This deliberately includes the implementation to exercise the exact
 * static diagonalize() implementation without exposing a new public API.
 *
 * cc -O3 -std=c99 -Iinclude test/bench_krylov_projected.c -lm -o bench
 * ./bench
 *
 * Compare revisions by defining KRYLOV_SOURCE to an alternate C file.
 */
#ifndef KRYLOV_SOURCE
#define KRYLOV_SOURCE "../src/esolver/lis_krylov.c"
#endif
#include KRYLOV_SOURCE
#include <stdio.h>
#include <time.h>

static void make_projected(int n, double *a)
{
    int i, j;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j) {
            double v = i == j ? 0.2*(double)(i + 1) :
                  0.45*sin(0.19*(double)(i+1)*(double)(j+1)) /
                  sqrt((double)(i+j+2));
            a[(size_t)i*n+j] = v;
        }
}

static int verify(int n, const double *orig,
                  const double *eigenvecs, const double *eigenvalues)
{
    double residual_square = 0.0, matrix_square = 0.0, orth_error = 0.0;
    int i, j, k;
    for (i = 0; i < n*n; ++i)
        matrix_square += orig[i]*orig[i];
    for (k = 0; k < n; ++k) {
        for (i = 0; i < n; ++i) {
            double left = 0.0;
            for (j = 0; j < n; ++j)
                left += orig[(size_t)i*n+j] *
                        eigenvecs[(size_t)j*n+k];
            left -= eigenvalues[k]*eigenvecs[(size_t)i*n+k];
            residual_square += left*left;
        }
        for (j = 0; j < n; ++j) {
            double overlap = 0.0;
            for (i = 0; i < n; ++i)
                overlap += eigenvecs[(size_t)i*n+j] *
                           eigenvecs[(size_t)i*n+k];
            orth_error = fmax(orth_error, fabs(overlap - (k == j)));
        }
    }
    if (sqrt(residual_square / fmax(matrix_square, 1.0)) > 1e-9 ||
        orth_error > 1e-10) {
        fprintf(stderr, "FAIL n=%d relative_residual=%g orth_error=%g\n",
                n, sqrt(residual_square / fmax(matrix_square, 1.0)),
                orth_error);
        return -1;
    }
    return 0;
}

static int bench(int n, int repeats)
{
    size_t elements = (size_t)n*n;
    double *original = (double *)malloc(elements*sizeof(double));
    double *a = (double *)malloc(elements*sizeof(double));
    double *q = (double *)malloc(elements*sizeof(double));
    double *eigenvalues = (double *)malloc((size_t)n*sizeof(double));
    clock_t begin, finish;
    double seconds;
    int i, rc = 0;

    if (!original || !a || !q || !eigenvalues) {
        free(original); free(a); free(q); free(eigenvalues);
        return 1;
    }
    make_projected(n, original);
    begin = clock();
    for (i = 0; i < repeats; ++i) {
        memcpy(a, original, elements*sizeof(double));
        rc = diagonalize(n, a, q, eigenvalues, n);
        if (rc) break;
    }
    finish = clock();
    if (!rc) rc = verify(n, original, q, eigenvalues);
    seconds = (double)(finish-begin) / CLOCKS_PER_SEC;
    printf("%d,%d,%.8f,%d\n", n, repeats, seconds, rc);
    free(original); free(a); free(q); free(eigenvalues);
    return !!rc;
}

int main(void)
{
    puts("n,repeats,cpu_seconds,failed");
    return bench(24, 1600) || bench(48, 480) ||
           bench(72, 120) || bench(96, 40);
}
