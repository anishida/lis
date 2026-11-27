#ifndef MULTIFRONTAL_LDLT_H
#define MULTIFRONTAL_LDLT_H
#include "delayed_ldlt.h"
/* Lower CSC, strictly increasing rows >= column. Missing diagonal is zero.
   parent[k] > k or UINT32_MAX (root); each structural offdiagonal row must
   be an ancestor of its column. One original column per front. */
typedef struct {
 uint32_t max_front; size_t max_work_bytes,max_factor_slots,max_delays;
 long double max_growth; /* 0 disables growth limit */
} p61_limits;
typedef struct {
 size_t fronts,peak_front,peak_work_bytes,factor_slots,pairs,scalars,swaps,delays;
 long double growth,max_multiplier;
} p61_stats;
/* Workspace budget counts requested bytes of allocations by this module,
   including temporary factor packing, excluding caller input, allocator
   metadata, existing context and final p56_rebuild allocations. */
int p61_factorize(p56_context **slot,uint32_t n,const size_t *col,
 const uint32_t *row,const double *val,const uint32_t *parent,
 const p61_limits *lim,p61_stats *stats);
#endif
