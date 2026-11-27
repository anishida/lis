#ifndef PACKED_LDLT_H
#define PACKED_LDLT_H
#include <stddef.h>
#include <stdint.h>
/* Real double, zero-based indices, factor ordering, unit lower triangular L.
   For each panel [first,first+width), rows are strictly increasing external
   row indices. Values are column-major: strict lower in-panel entries, then
   the external rows, for each column. No padding. */
typedef struct {
    uint32_t first, width, nrows;
    const uint32_t *rows;
    const double *values;
    size_t nvalues;
} p55_panel_input;
typedef struct p55_factor p55_factor;
enum { P55_OK=0, P55_INVALID=1, P55_NOMEM=2, P55_SINGULAR=3, P55_RANGE=4 };
/* kind: 1 scalar; 2 start of pair; 0 continuation. off[k] is D[k,k+1]
   for kind[k]==2. Pairs must be within one panel; L[k+1,k] must be zero.
   Inputs are copied. On failure *out is NULL. */
int p55_pack(p55_factor **out, uint32_t n, size_t np, const p55_panel_input *p,
             const uint8_t *kind, const double *diag, const double *off);
/* In-place solve L D L^T x = rhs. No allocation. On failure x may be changed.
   Caller must apply any external ordering/permutation before/after this call. */
int p55_solve(const p55_factor *f, double *x);
void p55_free(p55_factor *f);
size_t p55_block_count(const p55_factor *f);
size_t p55_index_bytes(const p55_factor *f);
#endif
