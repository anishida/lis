/* Generic selection of converged eigensolver pairs. */
#ifndef LIS_ESOLVER_SELECTION_H
#define LIS_ESOLVER_SELECTION_H

#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Sort accepted pairs by the real part of the eigenvalue. */
#define LIS_ESOLVER_SELECT_ASCENDING_REAL 0x01
/* Keep only pairs whose eigenvalue has a strictly positive real part. */
#define LIS_ESOLVER_SELECT_POSITIVE_REAL  0x02

typedef struct LIS_ESOLVER_SELECTED_PAIR_STRUCT
{
    LIS_INT mode;
    LIS_SCALAR value;
    LIS_REAL residual;
}
LIS_ESOLVER_SELECTED_PAIR;

/*
 * Select converged eigenpairs from the first scan_count stored results.
 *
 * A result is accepted when its eigenvalue and residual are finite and
 * residual <= residual_tolerance. Optional flags can sort the accepted
 * pairs by ascending real eigenvalue and/or reject non-positive real
 * eigenvalues. The function writes at most requested_count entries to
 * pairs and reports the total number of accepted results through
 * available_count when that pointer is non-NULL.
 *
 * scan_count must not exceed the number of results available from the
 * eigensolver. Access errors are propagated to the caller.
 */
LIS_INT lis_esolver_select_converged_pairs(
    LIS_ESOLVER esolver,
    LIS_INT scan_count,
    LIS_INT requested_count,
    LIS_REAL residual_tolerance,
    LIS_INT flags,
    LIS_ESOLVER_SELECTED_PAIR pairs[],
    LIS_INT *available_count);

#ifdef __cplusplus
}
#endif

#endif
