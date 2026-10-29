/* Generic cached CSR materialization for row-backed matrices. */
#ifndef LIS_MATRIX_MATERIALIZE_H
#define LIS_MATRIX_MATERIALIZE_H

#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

#define LIS_MATRIX_MATERIALIZE_REUSE   0
#define LIS_MATRIX_MATERIALIZE_VALUES  1
#define LIS_MATRIX_MATERIALIZE_REBUILD 2

typedef struct LIS_MATRIX_CSR_CACHE_STRUCT *LIS_MATRIX_CSR_CACHE;

/*
 * Return a cache-owned CSR materialization of source.
 *
 * The first call builds the CSR. Subsequent calls reuse it when the source
 * identity and effective structure/value generations are unchanged. A pure
 * value change refreshes the existing CSR allocation; a structural change
 * rebuilds it transactionally. Nested LIS_MATRIX_OPERATOR operands are
 * included in the generation stamp, so operand updates are observed even
 * when the top-level operator object itself is unchanged.
 *
 * The returned CSR pointer is borrowed from cache and remains valid until the
 * next successful rebuild or lis_matrix_csr_cache_destroy(). The initial
 * implementation is serial and requires row-accessible source matrices.
 */
LIS_INT lis_matrix_materialize_csr_cached(
    LIS_MATRIX source,
    LIS_MATRIX_CSR_CACHE *cache,
    LIS_MATRIX *csr,
    LIS_INT *action);

void lis_matrix_csr_cache_destroy(
    LIS_MATRIX_CSR_CACHE cache);

#ifdef __cplusplus
}
#endif

#endif
