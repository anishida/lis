#ifndef QMD_PLAN_H
#define QMD_PLAN_H
#include "minimum_degree.h"
/* Application-owned immutable symbolic plan. Pattern comparison is exact.
   Each construction phase has a separate budget: plan bytes (including
   temporary permutation records), QMD analysis, symbolic analysis.
   Numerical refresh redoes dynamic pivoting; it does not freeze L/D.
   No global cache. Serialize factor replacement against apply/destroy. */
typedef struct p71_plan p71_plan;
void p71_plan_free(p71_plan *plan);
int p71_plan_create(p71_plan **out,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const p63_limits *ol,const p62_limits *sl,uint32_t width,size_t max_plan_bytes);
/* The input dimensions and complete CSC pattern must match the plan.
   New factor and stats are committed only on success. */
int p71_refactor(const p71_plan *plan,p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p61_limits *nl,p63_stats *stats);
#endif
