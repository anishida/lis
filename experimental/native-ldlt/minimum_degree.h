#ifndef MINIMUM_DEGREE_H
#define MINIMUM_DEGREE_H
#include "symbolic_ldlt.h"
/* Exact explicit-graph minimum degree, deterministic ties by original index.
   Budgets apply separately to ordering and CSC permutation, not total RSS. */
typedef struct {size_t max_work_bytes,max_slots,max_pair_checks;} p63_limits;
typedef struct {uint32_t n,*perm,*inv;size_t slots,max_front,pair_checks,peak_work_bytes;} p63_order;
typedef struct {p63_order summary;size_t permutation_bytes;p61_stats numeric;} p63_stats;
/* summary.perm and summary.inv in stats are always NULL. */
void p63_free(p63_order *o);
int p63_analyze(p63_order **out,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const p63_limits *lim);
/* Maps original -> ordered CSC -> numerically pivoted factors -> original.
   Both context and stats unchanged on failure. */
int p63_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p63_limits *ol,const p62_limits *sl,const p61_limits *nl,p63_stats *stats);
/* Quotient graph: max_pair_checks bounds adjacency/element visits instead.
   summary.pair_checks reports these visits, not explicit clique pair checks. */
int p63q_analyze(p63_order **out,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const p63_limits *lim);
int p63q_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p63_limits *ol,const p62_limits *sl,const p61_limits *nl,uint32_t width,p63_stats *stats);
#endif
