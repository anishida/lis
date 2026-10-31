/* Generic exact-nullspace solution lift construction. */
#ifndef LIS_NULLSPACE_LIFT_H
#define LIS_NULLSPACE_LIFT_H

#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct LIS_NULLSPACE_LIFT_DIAGNOSTICS_STRUCT
{
    LIS_INT nullspace_dim;
    LIS_INT active_count;
    LIS_REAL rhs_norm;
    LIS_REAL max_relative_component;
    LIS_REAL lift_norm;
}
LIS_NULLSPACE_LIFT_DIAGNOSTICS;

/*
 * Build a bounded solution lift from the solver's exact right-nullspace.
 *
 * The exact basis stored by LIS is orthonormal. For every basis vector q,
 * the coefficient q^H b is measured relative to ||b||_2. Components whose
 * relative magnitude is greater than tolerance contribute
 *
 *     (q^H b / alpha) q
 *
 * to lift. The matrix and RHS are not modified. The caller owns lift and may
 * add it to a compatible solution after the projected solve.
 */
LIS_INT lis_solver_build_nullspace_lift(
    LIS_SOLVER solver,
    LIS_VECTOR rhs,
    LIS_REAL alpha,
    LIS_REAL tolerance,
    LIS_VECTOR lift,
    LIS_NULLSPACE_LIFT_DIAGNOSTICS *diagnostics);

#ifdef __cplusplus
}
#endif

#endif
