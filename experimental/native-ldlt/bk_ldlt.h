#ifndef BK_LDLT_H
#define BK_LDLT_H
#include "pivot_ldlt.h"
/* Classical lower Bunch–Kaufman selection on the full active Schur matrix.
   Swaps restricted to current symbolic panel; no fallback or delayed pivot.
   max_growth=0 disables rejection; otherwise >=1. Statistics on success only. */
typedef struct { size_t direct, secondary, alternate, pairs, swaps; long double growth, max_multiplier; } p59_stats;
enum { P59_GROWTH=11 };
int p59_factorize(p56_context **slot,uint32_t n,const double *a,size_t np,
 const p58_symbolic *sym,uint32_t max_n,long double max_growth,p59_stats *stats);
#endif
