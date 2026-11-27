#ifndef LIS_LDLT_USER_BRIDGE_H
#define LIS_LDLT_USER_BRIDGE_H
#include "lis.h"
/* Register once after lis_initialize, before creating/using solvers.
 * Experimental application-linked bridge; not an installed LIS API. */
LIS_INT lis_ldlt_experimental_register(void);
/* Only ldltexpr owns a reusable plan. Explicit refresh for an application-
 * owned preconditioner before lis_solve_kernel; not automatic sequence
 * refresh. Failure retains the old factor, which remains stale. Caller
 * serializes matrix changes, refresh and apply. */
LIS_INT lis_ldlt_experimental_refresh(LIS_SOLVER solver,LIS_PRECON precon);
/* Per-calling-thread opt-in. Requires OpenMP-enabled prototype for n>1. */
LIS_INT lis_ldlt_experimental_set_threads(unsigned n);
#endif
