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

#define LIS_MATRIX_POLICY_DIRECT       0
#define LIS_MATRIX_POLICY_CACHED_CSR   1
#define LIS_MATRIX_POLICY_ACTION_DIRECT (-1)

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
 * implementation requires a single-process row distribution. On OpenMP
 * builds, sufficiently large USER/OPERATOR matrices are built using
 * independent per-thread row sources only when every USER leaf has explicitly
 * opted in to concurrent row callback access. Refresh stays serial after a
 * full index-by-index pattern check. Other callbacks and smaller inputs keep
 * the serial path. No matrix values or callback state may be mutated while
 * materialization is in progress.
 */
LIS_INT lis_matrix_materialize_csr_cached(
    LIS_MATRIX source,
    LIS_MATRIX_CSR_CACHE *cache,
    LIS_MATRIX *csr,
    LIS_INT *action);


/*
 * Select the matrix representation for a configured solver.
 *
 * This policy is applied to LIS_MATRIX_USER and LIS_MATRIX_OPERATOR sources.
 * Native CSR and other explicit formats are returned unchanged so existing
 * LIS format/conversion behavior is preserved.
 *
 * USER/OPERATOR matrices are passed through directly when the configured
 * preconditioner can consume their available callbacks. Otherwise, a cached
 * CSR view is selected. The cache tracks source generations and reports
 * REUSE, VALUES, or REBUILD through action.
 *
 * selected is borrowed. For DIRECT it aliases source; for CACHED_CSR it is
 * owned by cache and remains valid until a successful rebuild or cache
 * destruction. action is LIS_MATRIX_POLICY_ACTION_DIRECT for the direct path.
 *
 * The function does not change solver configuration and does not call
 * lis_solve(); the caller must use the returned selected matrix.
 *
 * This is an advanced, opt-in helper declared in this separate header,
 * not an implicit part of lis_solve(). Neither the source nor the
 * returned selected matrix should be destroyed by the helper's caller.
 * For CACHED_CSR, the cache owns the matrix and may replace its pointer
 * on a subsequent structural update; no borrowed pointer may be retained
 * beyond a successful rebuild or explicit cache destruction.
 * On failed materialization, selected remains NULL and errors propagate.
 * The cached materialization fallback is currently serial-only.
 */
LIS_INT lis_matrix_prepare_solver_matrix(
    LIS_MATRIX source,
    LIS_SOLVER solver,
    LIS_MATRIX_CSR_CACHE *cache,
    LIS_MATRIX *selected,
    LIS_INT *policy,
    LIS_INT *action);

void lis_matrix_csr_cache_destroy(
    LIS_MATRIX_CSR_CACHE cache);

#ifdef __cplusplus
}
#endif

#endif
