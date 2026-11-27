#ifndef SYMBOLIC_LDLT_H
#define SYMBOLIC_LDLT_H
#include "multifrontal_ldlt.h"
/* Lower CSC structure (including explicit zeros), natural ordering.
   Counts include diagonal even when omitted. Supernodes are strict,
   contiguous, and diagnostic only; numerical stage still uses one column/front. */
typedef struct {
 size_t max_work_bytes,max_slots; /* slots excludes unit diagonal */
} p62_limits;
typedef struct {
 uint32_t n; uint32_t *parent,*count,*super;
 size_t nsuper,roots,slots,max_front,peak_work_bytes;
} p62_symbolic;
/* *out must initially be NULL or owned by p62_analyze. Unchanged on error.
   nnz is the actual length of row; col has n+1 entries. */
int p62_analyze(p62_symbolic **out,uint32_t n,size_t nnz,const size_t *col,
 const uint32_t *row,const p62_limits *lim);
void p62_free(p62_symbolic *s);
/* Automatic symbolic + existing numerical stage; separate budgets, no
   global peak RAM guarantee. Both slot and stats unchanged on failure. */
int p62_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,
 const uint32_t *row,const double *val,const p62_limits *slim,
 const p61_limits *nlim,p61_stats *stats);
#endif
