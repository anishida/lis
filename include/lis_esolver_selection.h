/* Generic selection of converged eigensolver pairs. */
#ifndef LIS_ESOLVER_SELECTION_H
#define LIS_ESOLVER_SELECTION_H

#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

#define LIS_ESOLVER_SELECT_ASCENDING_REAL 0x01
#define LIS_ESOLVER_SELECT_POSITIVE_REAL  0x02

typedef struct LIS_ESOLVER_SELECTED_PAIR_STRUCT
{
    LIS_INT mode;
    LIS_SCALAR value;
    LIS_REAL residual;
}
LIS_ESOLVER_SELECTED_PAIR;

typedef LIS_INT (*LIS_ESOLVER_SELECTED_PAIR_VISITOR)(
    void *context,
    LIS_INT output_index,
    LIS_INT source_mode,
    LIS_SCALAR value,
    LIS_REAL residual,
    LIS_VECTOR vector);

LIS_INT lis_esolver_select_converged_pairs(
    LIS_ESOLVER esolver,
    LIS_INT scan_count,
    LIS_INT requested_count,
    LIS_REAL residual_tolerance,
    LIS_INT flags,
    LIS_ESOLVER_SELECTED_PAIR pairs[],
    LIS_INT *available_count);

LIS_INT lis_esolver_copy_selected_eigenpairs(
    LIS_ESOLVER esolver,
    LIS_INT scan_count,
    LIS_INT requested_count,
    LIS_REAL residual_tolerance,
    LIS_INT flags,
    LIS_SCALAR values[],
    LIS_VECTOR vectors[],
    LIS_INT *available_count);

/*
 * Select converged pairs and visit them one at a time using one caller-owned
 * work vector. Before each callback, work_vector contains the selected
 * eigenvector for source_mode. This avoids allocating O(n*requested_count)
 * output vectors in embedding applications that stream results into their own
 * storage.
 */
LIS_INT lis_esolver_visit_converged_pairs(
    LIS_ESOLVER esolver,
    LIS_INT scan_count,
    LIS_INT requested_count,
    LIS_REAL residual_tolerance,
    LIS_INT flags,
    LIS_VECTOR work_vector,
    LIS_ESOLVER_SELECTED_PAIR_VISITOR visitor,
    void *context,
    LIS_INT *available_count);

#ifdef __cplusplus
}
#endif

#endif
