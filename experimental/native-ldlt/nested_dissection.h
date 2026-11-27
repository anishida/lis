#ifndef NESTED_DISSECTION_H
#define NESTED_DISSECTION_H
#include "minimum_degree.h"
typedef enum {P64_NATURAL=0,P64_MD=1,P64_ND=2,P64_AUTO=3} p64_method;
typedef struct {size_t max_work_bytes,max_slots,max_pair_checks; uint32_t leaf_size,max_depth,max_front;} p64_limits;
typedef struct {size_t slots,front,work_bytes;long double work;int status;} p64_estimate;
typedef struct {p64_method selected;p64_estimate candidates[3];p61_stats numeric;} p64_stats;
/* Factor order -> original. All limits are per phase, not total process RAM.
   In AUTO, candidates failing a resource limit are skipped. No numeric retry.
   Output context and stats unchanged on error. */
int p64_order(uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,
 const p64_limits *lim,const p62_limits *sl,p64_method method,p63_order **out,p64_stats *stats);
int p64_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,
 const p64_limits *lim,const p62_limits *sl,const p61_limits *nl,p64_method method,p64_stats *stats);
#endif
