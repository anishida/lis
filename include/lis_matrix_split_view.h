/* Non-owning split-triangular external matrix view. */
#ifndef LIS_MATRIX_SPLIT_VIEW_H
#define LIS_MATRIX_SPLIT_VIEW_H

#include <stdint.h>
#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct LIS_MATRIX_SPLIT_TRIANGULAR_VIEW_STRUCT
    *LIS_MATRIX_SPLIT_TRIANGULAR_VIEW;

typedef struct LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY_STRUCT
    *LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY;

typedef struct LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW_STRUCT
    *LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW;

/*
 * Create a serial LIS_MATRIX_USER view over split triangular storage.
 *
 * ptr/index describe the strictly lower triangle by columns. For column j,
 * every index entry must identify a row i>j and row indices within that
 * column must be strictly increasing. origin selects 0- or 1-based ptr/index
 * storage and must be either 0 or 1.
 *
 * diag and lower_values are required. upper_values==NULL mirrors
 * lower_values (symmetric storage); otherwise upper_values[k] supplies the
 * opposite entry A(j,i) corresponding to lower_values[k] == A(i,j).
 *
 * The application retains ownership of ptr/index/value arrays and must keep
 * them alive while the view is used. The returned view owns only compact
 * integer row metadata. Destroy the LIS_MATRIX before destroying the view.
 *
 * Row access is read-only and thread-safe as long as the application does not
 * modify or replace the external numerical arrays concurrently with matrix
 * operations. Use lis_matrix_split_triangular_view_set_values() only between
 * such operations.
 */
LIS_INT lis_matrix_create_split_triangular_view(
    LIS_Comm comm,
    LIS_INT n,
    LIS_INT nnz,
    const LIS_INT ptr[],
    const LIS_INT index[],
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[],
    LIS_INT origin,
    LIS_MATRIX *A,
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW *view);

/*
 * Fixed-width external-index variant for embedding applications whose sparse
 * structure uses signed 64-bit indices independently of the LIS build. The
 * structure is range-checked and narrowed only into temporary LIS_INT arrays;
 * those temporary arrays are released after the ordinary split-view row
 * metadata has been built. Numerical arrays remain non-owning and zero-copy.
 */
LIS_INT lis_matrix_create_split_triangular_view_i64(
    LIS_Comm comm,
    int64_t n,
    int64_t nnz,
    const int64_t ptr[],
    const int64_t index[],
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[],
    int64_t origin,
    LIS_MATRIX *A,
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW *view);

/*
 * Build immutable expanded-row metadata once and share it between multiple
 * split-triangular matrices that have the same sparse structure.
 *
 * The topology stores only LIS-owned integer row metadata. The input ptr/index
 * arrays are validated during construction and are not retained afterwards.
 * Topologies are reference-counted: every shared view retains one reference,
 * so the caller may destroy its own topology handle immediately after all
 * desired views have been created.
 */
LIS_INT lis_matrix_create_split_triangular_topology(
    LIS_INT n,
    LIS_INT nnz,
    const LIS_INT ptr[],
    const LIS_INT index[],
    LIS_INT origin,
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY *topology);

LIS_INT lis_matrix_create_split_triangular_topology_i64(
    int64_t n,
    int64_t nnz,
    const int64_t ptr[],
    const int64_t index[],
    int64_t origin,
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY *topology);

/*
 * Create a serial LIS_MATRIX_USER backed by an existing immutable topology.
 * Numerical arrays remain non-owning and may differ between views that share
 * the same topology. upper_values==NULL selects symmetric mirrored storage.
 */
LIS_INT lis_matrix_create_split_triangular_shared_view(
    LIS_Comm comm,
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology,
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[],
    LIS_MATRIX *A,
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW *view);

/* Replace only the numerical pointers of a shared-topology view. */
LIS_INT lis_matrix_split_triangular_shared_view_set_values(
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view,
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[]);

void lis_matrix_split_triangular_shared_view_destroy(
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view);

/*
 * Borrow external real64 numerical storage without allocating a LIS_SCALAR
 * copy. The same immutable topology (including an i64-imported topology)
 * can back multiple real64 and native LIS_SCALAR matrices.
 *
 * Values are converted only on row/diagonal/matvec reads. The imaginary
 * component is zero for complex LIS. In wider-precision configurations the
 * source remains IEEE double precision; importing does not create precision.
 *
 * Arrays belong to the caller, must outlive all use of the view, and must
 * not be modified while an LIS operation is active. For updates in place,
 * call set_values_real64 even when the addresses do not change, so generation
 * and CSR caches are invalidated. Symmetric/asymmetric mode is immutable,
 * and real64 and native LIS_SCALAR setters cannot be mixed.
 * Destroy LIS_MATRIX before destroying the returned view.
 */
LIS_INT lis_matrix_create_split_triangular_shared_view_real64(
    LIS_Comm comm,
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology,
    const double diag[],
    const double lower_values[],
    const double upper_values[],
    LIS_MATRIX *A,
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW *view);

LIS_INT lis_matrix_split_triangular_shared_view_set_values_real64(
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view,
    const double diag[],
    const double lower_values[],
    const double upper_values[]);

void lis_matrix_split_triangular_topology_destroy(
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology);

/* Replace only the non-owning numerical array pointers and publish a value
 * generation change to the associated LIS matrix. Symmetric/nonsymmetric mode
 * is fixed at creation time and cannot be changed by this call. */
LIS_INT lis_matrix_split_triangular_view_set_values(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view,
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[]);

/*
 * Atomically replace the strictly-lower structure and all external numerical
 * array pointers while keeping the matrix size n unchanged. nnz and the
 * symmetric/nonsymmetric storage mode may change. New row metadata is built
 * and validated before the active view is modified. On success the associated
 * matrix structure generation is advanced; on validation/allocation failure
 * the previous view remains usable and unchanged.
 */
LIS_INT lis_matrix_split_triangular_view_rebind(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view,
    LIS_INT nnz,
    const LIS_INT ptr[],
    const LIS_INT index[],
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[],
    LIS_INT origin);

/* Checked signed-64-bit structure variant of the transactional rebind. */
LIS_INT lis_matrix_split_triangular_view_rebind_i64(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view,
    int64_t nnz,
    const int64_t ptr[],
    const int64_t index[],
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[],
    int64_t origin);

void lis_matrix_split_triangular_view_destroy(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view);

#ifdef __cplusplus
}
#endif

#endif
