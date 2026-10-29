/* Transactional solver option-profile and runtime helpers. */

#ifndef __LIS_SOLVER_OPTIONS_H__
#define __LIS_SOLVER_OPTIONS_H__

#include "lis.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * Replace the complete text-configurable solver option profile.
 *
 * Unlike lis_solver_set_option(), which updates only options present in the
 * supplied string, this routine starts from LIS defaults, applies the string,
 * and then installs the resulting options/parameters on the target solver.
 *
 * Application-owned solver configuration that is not part of the option text
 * (callbacks, exact/left/near nullspaces, sequence policy, thread override and
 * half-storage policy) is preserved.  Solver work, sequence preconditioner and
 * recycle state are invalidated before the new profile becomes active.
 */
extern LIS_INT lis_solver_set_option_replace(
    const char *text,
    LIS_SOLVER solver);

/*
 * Replace the complete text-configurable solver profile only when
 * its effective options or parameters differ from the active profile.
 * An equivalent replacement preserves reusable solver state.
 */
extern LIS_INT lis_solver_set_option_replace_if_changed(
    const char *text,
    LIS_SOLVER solver,
    LIS_INT *changed);

/*
 * Typed solve-time controls.  These values are read for each solve and do not
 * change solver work dimensions or invalidate persistent preconditioner state.
 * Tolerance must be finite and non-negative; zero is accepted to retain
 * the underlying LIS solver-option semantics. Maximum iterations may be
 * zero (zero permitted iterations) but cannot be negative.
 */
extern LIS_INT lis_solver_set_tolerance(
    LIS_SOLVER solver,
    LIS_REAL tolerance);

extern LIS_INT lis_solver_get_tolerance(
    LIS_SOLVER solver,
    LIS_REAL *tolerance);

extern LIS_INT lis_solver_set_maxiter(
    LIS_SOLVER solver,
    LIS_INT maxiter);

extern LIS_INT lis_solver_get_maxiter(
    LIS_SOLVER solver,
    LIS_INT *maxiter);

/*
 * Compact snapshot of diagnostics already exposed by individual getters.
 * Time fields intentionally use double to match lis_solver_get_timeex().
 * Before any solve, values reflect the initialized solver state. After a
 * failed solve, status is the most recent return status and other fields
 * are whatever the underlying solver getters currently report; they do
 * not assert successful convergence.
 */
typedef struct
{
    LIS_INT status;
    LIS_INT iterations;
    LIS_INT iterations_double;
    LIS_INT iterations_quad;
    LIS_REAL residual_norm;
    double time;
    double iteration_time;
    double preconditioner_time;
    double preconditioner_create_time;
    double preconditioner_iteration_time;
} LIS_SOLVE_REPORT;

extern LIS_INT lis_solver_get_report(
    LIS_SOLVER solver,
    LIS_SOLVE_REPORT *report);

#ifdef __cplusplus
}
#endif

#endif
