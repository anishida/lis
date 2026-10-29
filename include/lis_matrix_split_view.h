/* Non-owning split-triangular external matrix view. */
#ifndef LIS_MATRIX_SPLIT_VIEW_H
#define LIS_MATRIX_SPLIT_VIEW_H

#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct LIS_MATRIX_SPLIT_TRIANGULAR_VIEW_STRUCT
    *LIS_MATRIX_SPLIT_TRIANGULAR_VIEW;

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

void lis_matrix_split_triangular_view_destroy(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view);

#ifdef __cplusplus
}
#endif

#endif
