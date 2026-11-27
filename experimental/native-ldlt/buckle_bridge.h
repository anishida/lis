#ifndef BUCKLE_BRIDGE_H
#define BUCKLE_BRIDGE_H
#include "parallel_ldlt.h"
/* CSR real binary64; arrays borrowed for the duration of a call. Sorted,
   strictly increasing columns per row. Full storage requires exact symmetry,
   including structure. Missing entries are zero. base is 0 or 1. */
typedef enum { P69_FULL=0, P69_LOWER=1, P69_UPPER=2 } p69_storage;
typedef struct {uint32_t n; size_t nnz; const size_t *ptr; const uint32_t *index; const double *value; unsigned base; p69_storage storage;} p69_csr;
typedef struct {size_t max_conversion_bytes; p64_limits ordering; p62_limits symbolic; p61_limits numeric; p64_method method; uint32_t width; unsigned threads;} p69_options;
typedef struct {size_t input_k,input_g,shifted_nnz,conversion_bytes; p64_stats ordering; p67_stats parallel;} p69_stats;
/* K - sigma*G. G may be NULL only for sigma == 0. Both matrices can use
   different bases/storage. No symmetry averaging, no discarded nonzeros.
   On failure slot and stats are unchanged. Caller serializes context changes.
   Conversion budget counts requested CSC arrays, not other phases. */
/* Validate CSR without allocating or changing input. */
int p69_validate_csr(const p69_csr *a);
void p69_defaults(p69_options *o);
int p69_factorize(p56_context **slot,const p69_csr *k,const p69_csr *g,double sigma,const p69_options *o,p69_stats *stats);
#endif
