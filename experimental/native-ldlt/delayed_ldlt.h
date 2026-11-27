#ifndef DELAYED_LDLT_H
#define DELAYED_LDLT_H
#include "bk_ldlt.h"
/* Dense reference, sequential panels (not an elimination tree). All limits
   mandatory except max_growth=0 disables growth limit. Zero delays/pending/
   extra slots means none allowed. max_front caps eligible+pending at entry.
   max_extra_slots counts allocated L doubles versus original symbolic slots,
   not set difference of original nonzero positions or total process memory.
   Output and stats unchanged on failure. */
typedef struct {uint32_t max_n,max_front;size_t max_delays,max_pending,max_extra_slots;long double max_growth;} p60_limits;
typedef struct {size_t direct,secondary,alternate,pairs,swaps,delays,peak_pending,peak_front,baseline_slots,factor_slots,panels;long double growth,max_multiplier;} p60_stats;
enum {P60_LIMIT=12};
int p60_factorize(p56_context **slot,uint32_t n,const double *a,size_t np,const p58_symbolic *sym,const p60_limits *limits,p60_stats *stats);
#endif
