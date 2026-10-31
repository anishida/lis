/* Generic effective solver/preconditioner capability queries. */

#ifndef __LIS_SOLVER_CAPABILITIES_H__
#define __LIS_SOLVER_CAPABILITIES_H__

#include "lis.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Effective configuration capabilities. */
#define LIS_SOLVER_CAP_FLEXIBLE_PRECON   0x01
#define LIS_SOLVER_CAP_RECYCLING         0x02
#define LIS_SOLVER_CAP_NEAR_NULLSPACE    0x04
#define LIS_SOLVER_CAP_NESTED_PRECON     0x08

/*
 * Return the effective leaf preconditioner selected by the solver.
 *
 * For an ordinary preconditioner, *precon_type is LIS_OPTIONS_PRECON and
 * *nested is LIS_FALSE.  For the supported one-level Hybrid composition,
 * *precon_type is LIS_OPTIONS_PPRECON and *nested is LIS_TRUE.
 *
 * Recursive Hybrid composition is intentionally not interpreted here because
 * the current Hybrid implementation rejects it during construction.
 */
extern LIS_INT lis_solver_get_effective_preconditioner(
    LIS_SOLVER solver,
    LIS_INT *precon_type,
    LIS_INT *nested);

/*
 * Return a bit mask describing capabilities of the currently configured
 * solver/preconditioner composition.  The query is configuration-only: it
 * does not require lis_solve_setup() and does not allocate solver state.
 */
extern LIS_INT lis_solver_get_effective_capabilities(
    LIS_SOLVER solver,
    LIS_INT *capabilities);

#ifdef __cplusplus
}
#endif

#endif
