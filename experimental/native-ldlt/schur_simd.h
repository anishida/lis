#ifndef SCHUR_SIMD_H
#define SCHUR_SIMD_H
#include <stddef.h>
typedef enum {P66_AUTO=0,P66_REFERENCE=1,P66_AVX2=2,P66_AVX512=3} p66_mode;
/* Per-thread selection; returns 0 for invalid/unavailable mode. */
int p66_available(p66_mode mode);
int p66_set_mode(p66_mode mode);
p66_mode p66_effective_mode(void);
const char *p66_mode_name(p66_mode mode);
/* Internal triangular update: pivot width must be 1 or 2.
   Scratch may be partly modified on failure. */
typedef int (*p66_kernel)(double *,size_t,size_t,size_t,long double *);
p66_kernel p66_select_kernel(void);
/* Active symmetric front stored in its lower triangle only. */
p66_kernel p66_select_kernel_lower(void);
#endif
