/* Optional serial double-precision static condensation core for LIS.
 * No dependency on LIS, SuperLU, BLAS, LAPACK or other numerical libraries.
 * All sparse rows must contain strictly ascending, unique column indices.
 */
#ifndef LIS_CONDENSE_CORE_H
#define LIS_CONDENSE_CORE_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t n, nnz;
    const size_t *rowptr; /* length n+1 */
    const size_t *colind; /* length nnz */
    const double *values; /* length nnz */
} lis_sc_csr;

typedef struct {
    size_t nblocks;
    const size_t *block_ptr; /* length nblocks+1, block_ptr[0]=0 */
    const size_t *block_dof; /* length block_ptr[nblocks] */
} lis_sc_blocks;

typedef enum {
    LIS_SC_OK = 0,
    LIS_SC_BAD_ARGUMENT,
    LIS_SC_BAD_CSR,
    LIS_SC_BAD_MAP,
    LIS_SC_COUPLED_BLOCKS,
    LIS_SC_SINGULAR_BLOCK,
    LIS_SC_NONFINITE,
    LIS_SC_SIZE,
    LIS_SC_OUT_OF_MEMORY,
    LIS_SC_PATTERN_MISMATCH,
    LIS_SC_UNSUPPORTED
} lis_sc_status;

typedef enum {
    LIS_SC_FULL = 0,
    LIS_SC_NUMERIC = 1, /* reuse Schur sparsity + block mapping; refactor numerical values */
    LIS_SC_REUSE_ALL = 2
} lis_sc_reuse_mode;

typedef struct lis_sc_plan lis_sc_plan;

#ifdef LIS_SC_TEST_ALLOC_FAIL
/* Test-only allocation failure injection. Not compiled in production LIS. */
void lis_sc_test_fail_after(long allocation_index);
#endif

/* Transactional build: returns a fresh, owned candidate without changing base.
 * FULL and NUMERIC return non-NULL *candidate on success.
 * REUSE_ALL returns NULL on success if the full operator and metadata match;
 * this means the base is reusable and must remain alive.
 * On error *candidate is always NULL and base remains valid.
 */
lis_sc_status lis_sc_stage(lis_sc_plan **candidate, const lis_sc_plan *base,
                           const lis_sc_csr *a, const lis_sc_blocks *blocks,
                           lis_sc_reuse_mode mode, double relative_pivot_tolerance);

/* Atomic replacement: on any error, *slot is unchanged, including its matrix. */
lis_sc_status lis_sc_build(lis_sc_plan **slot, const lis_sc_csr *a,
                            const lis_sc_blocks *blocks,
                            lis_sc_reuse_mode mode, double relative_pivot_tolerance);
void lis_sc_free(lis_sc_plan *plan);
const lis_sc_csr *lis_sc_matrix(const lis_sc_plan *plan); /* borrowed */
size_t lis_sc_physical_dofs(const lis_sc_plan *plan);
size_t lis_sc_internal_dofs(const lis_sc_plan *plan);
const size_t *lis_sc_physical_to_full(const lis_sc_plan *plan); /* borrowed */

/* Vector lengths: b_full[n], g[n_phys], u[n_phys], x_full[n]. */
lis_sc_status lis_sc_rhs(const lis_sc_plan *plan, const double *b_full,
                          double *g);
lis_sc_status lis_sc_reconstruct(const lis_sc_plan *plan, const double *u,
                                  const double *b_full, double *x_full);
lis_sc_status lis_sc_true_residual(const lis_sc_plan *plan, const double *b_full,
                                    const double *x_full, double *relres);

/* Modes are contiguous vectors, mode_count vectors of length n (input) or n_phys (output). */
lis_sc_status lis_sc_project_modes(const lis_sc_plan *plan, size_t mode_count,
                                    const double *full_modes, double *physical_modes);
lis_sc_status lis_sc_project_labels(const lis_sc_plan *plan,
                                     const long *full_labels, long *physical_labels);
const char *lis_sc_status_string(lis_sc_status status);

#ifdef __cplusplus
}
#endif
#endif
