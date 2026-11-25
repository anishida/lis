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
 * Recursive Hybrid and Hybrid with MFAMG_ILUTP are rejected, matching
 * the current Hybrid construction restrictions.
 */
extern LIS_INT lis_solver_get_effective_preconditioner(
    LIS_SOLVER solver,
    LIS_INT *precon_type,
    LIS_INT *nested);

/*
 * Return a bit mask describing capabilities of the currently configured
 * solver/preconditioner composition.  The query is configuration-only: it
 * does not require lis_solve_setup() and does not allocate solver state.
 * NEAR_NULLSPACE indicates that the configured MFAMG leaf consumes a basis;
 * it does not certify solver/precision compatibility or matrix eligibility.
 */
extern LIS_INT lis_solver_get_effective_capabilities(
    LIS_SOLVER solver,
    LIS_INT *capabilities);

/*
 * Generic preflight for exact and near-nullspace settings.  The returned
 * bitmask is zero when the current solver profile is compatible.
 * Query-only: neither nullspace vectors nor numerical values are changed.
 * An application may report the conflicting settings before a large setup.
 */
#define LIS_NULLSPACE_CONFIG_EXACT_SOLVER     0x01
#define LIS_NULLSPACE_CONFIG_EXACT_SCALE      0x02
#define LIS_NULLSPACE_CONFIG_EXACT_PRECISION  0x04
/* Unsupported outer solver / effective preconditioner combination.
 * The current solve kernel checks the top-level preconditioner type;
 * nested Hybrid+MFAMG is not treated as direct MFAMG here.
 * MFAMG_ILUTP retains the flexible outer solver requirement. */
#define LIS_NULLSPACE_CONFIG_NEAR_COMBINATION 0x08
#define LIS_NULLSPACE_CONFIG_NEAR_PRECISION   0x10

extern LIS_INT lis_solver_check_nullspace_configuration(
    LIS_SOLVER solver,
    LIS_INT *issues);

/* Read-only queries for application-owned configuration copied by LIS.
 * A zero dimension/count denotes an absent basis/partition. */
extern LIS_INT lis_solver_get_near_nullspace_dim(
    LIS_SOLVER solver,LIS_INT *dim);
extern LIS_INT lis_solver_get_mfamg_contact_configuration(
    LIS_SOLVER solver,LIS_INT *enabled,LIS_INT *nlocal);
/* Compare label values and length, not pointer identity. For an empty
 * partition pass nlocal=0 and NULL. Does not reset persistent state.
 * Invalid arguments leave the output unchanged. */
extern LIS_INT lis_solver_compare_mfamg_contact_partition(
    LIS_SOLVER solver,LIS_INT nlocal,const LIS_INT *partition,LIS_INT *equal);

#ifdef __cplusplus
}
#endif

#endif
