/* Generic graph queries derived from matrix entries. */

#ifndef __LIS_MATRIX_GRAPH_H__
#define __LIS_MATRIX_GRAPH_H__

#include "lis.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * Compute weakly connected components of the matrix graph.
 *
 * Every off-diagonal entry whose algebraically combined magnitude is greater
 * than drop_tolerance creates an undirected edge.  Duplicate column entries in
 * one row are combined before the threshold test, matching matrix semantics.
 *
 * labels must contain space for A->n entries.  Returned labels are consecutive
 * integers in [0,*component_count).  CSR, row-backed LIS_MATRIX_USER and
 * explicit formats convertible to CSR are supported.  The initial API is
 * intentionally serial; distributed matrices return LIS_ERR_NOT_IMPLEMENTED.
 */
extern LIS_INT lis_matrix_connected_components(
    LIS_MATRIX A,
    LIS_REAL drop_tolerance,
    LIS_INT *labels,
    LIS_INT *component_count);

#ifdef __cplusplus
}
#endif

#endif
