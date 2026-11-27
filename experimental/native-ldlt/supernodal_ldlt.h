#ifndef SUPERNODAL_LDLT_H
#define SUPERNODAL_LDLT_H
#include "nested_dissection.h"
/* Internal numeric entry. bounds partitions contiguous parent chains;
   NULL retains one original column per front. n+1 CSC offsets required. */
int p65_numeric(p56_context **slot,uint32_t n,const size_t *col,const uint32_t *row,const double *val,const uint32_t *parent,const uint32_t *bounds,size_t ng,const p61_limits *lim,p61_stats *stats);
/* Strict symbolic supernodes, split at max_width. Width 1 is legacy path.
   Symbolic limits include temporary boundaries; numeric limits independent. */
int p65_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p62_limits *sl,const p61_limits *nl,uint32_t max_width,p61_stats *stats);
int p65_ordered_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p64_limits *ol,const p62_limits *sl,const p61_limits *nl,p64_method method,uint32_t max_width,p64_stats *stats);
#endif
