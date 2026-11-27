#ifndef PIVOT_LDLT_H
#define PIVOT_LDLT_H
#include "lis_adapter.h"
/* Symbolic panel: contiguous columns and sorted external row list. */
typedef struct { uint32_t first,width,nrows; const uint32_t *rows; } p58_symbolic;
typedef struct { size_t scalar, pairs, swaps; } p58_stats;
enum { P58_NO_PIVOT=9, P58_PATTERN=10 };
/* Dense reference numeric factorization, row-major real symmetric A.
   max_n is mandatory and limits quadratic workspace; caller supplies A[n*n].
   tau in (0,1], scalar threshold and reciprocal multiplier bound for pairs.
   Swaps are confined to one panel and preserve earlier panel row membership.
   Reject any nonzero L outside symbolic pattern. Transactional on failure;
   stats, if non-NULL, updated only on success. */
int p58_factorize(p56_context **slot,uint32_t n,const double *a,size_t np,
 const p58_symbolic *sym,double tau,uint32_t max_n,p58_stats *stats);
#endif
