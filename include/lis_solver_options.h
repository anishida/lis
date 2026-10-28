/* Transactional solver option-profile helpers. */

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

#ifdef __cplusplus
}
#endif

#endif
