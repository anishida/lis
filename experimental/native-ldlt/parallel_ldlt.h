#ifndef PARALLEL_LDLT_H
#define PARALLEL_LDLT_H
#include "supernodal_ldlt.h"
typedef struct {p61_stats numeric;unsigned threads_used;size_t levels,max_ready;} p67_stats;
/* Level-scheduled independent fronts, 1..64 requested threads. Build with
   -fopenmp. Caller serializes context replacement. Scratch budget shared
   across workers; excludes runtime and final adapter allocation. */
int p67_numeric(p56_context **slot,uint32_t n,const size_t *col,const uint32_t *row,const double *val,const uint32_t *parent,const uint32_t *bounds,size_t ng,const p61_limits *lim,unsigned threads,p67_stats *stats);
int p67_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p62_limits *sl,const p61_limits *nl,uint32_t max_width,unsigned threads,p67_stats *stats);
int p67_ordered_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *cp,const uint32_t *ri,const double *v,const p64_limits *ol,const p62_limits *sl,const p61_limits *nl,p64_method method,uint32_t max_width,unsigned threads,p64_stats *stats,p67_stats *parallel);
#endif
