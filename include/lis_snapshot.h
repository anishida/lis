/* Optional, solver-independent Matrix Market snapshot exporter.
 * No FEA model, mesh, DOF numbering or contact assumptions are made.
 * This advanced API is explicitly invoked by the embedding application.
 */
#ifndef LIS_SNAPSHOT_H
#define LIS_SNAPSHOT_H

#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Export the exact current assembled linear operator and optional vectors.
 *
 * Creates <prefix>_A.mtx, optionally _b.mtx, _x0.mtx, _x.mtx and _modes.mtx.
 * Modes are an n-by-basis_count Matrix Market array (column-major).
 * The distinct modes suffix avoids b/B collisions on Windows.
 *
 * Serial real-valued matrices only in this first version. Explicit CSR,
 * other convertible explicit formats, and row-backed USER/OPERATOR sources
 * are supported; unsupported/matrix-free sources without row access fail.
 * This does NOT solve or modify the operator, run FEA, or create a DOF map.
 *
 * All arguments are borrowed. A failed call returns nonzero; files from
 * earlier successful writes may remain, so discard the entire prefix set
 * on error. Use a unique prefix for each Newton/nonlinear snapshot.
 */
LIS_INT lis_snapshot_export_mm(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x0,
    LIS_VECTOR x,
    LIS_INT basis_count,
    LIS_VECTOR basis[],
    const char *prefix);

/*
 * Optional physical/entity DOF labels from the embedding application.
 *
 * The equation index is zero-based and enumerated in LIS active equation
 * order. Pointers returned by the callback are borrowed until it returns;
 * node and dof must be nonempty strings; body may be NULL/empty.  Labels
 * can identify physical nodes, edges, multiplier entities, components,
 * thermomechanical fields, or arbitrary application DOF types.
 *
 * Output: <prefix>_dofs.csv with columns eq,node,dof,body.
 * CSV quoting/escaping is handled by LIS. No specific FEA package, mesh
 * numbering, element type or DOFs per node is assumed.
 *
 * This is opt-in and independent of lis_snapshot_export_mm(); it never
 * modifies A. Serial equations only in this first implementation.
 */
typedef LIS_INT (*LIS_SNAPSHOT_DOF_MAP_CALLBACK)(
    void *user_data,
    LIS_INT active_equation,
    const char **node,
    const char **dof,
    const char **body);

LIS_INT lis_snapshot_export_dof_map_csv(
    LIS_MATRIX A,
    LIS_SNAPSHOT_DOF_MAP_CALLBACK callback,
    void *user_data,
    const char *prefix);

#ifdef __cplusplus
}
#endif

#endif
