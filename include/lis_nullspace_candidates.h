/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Generic exact/near-nullspace candidate classification.
 */
#ifndef LIS_NULLSPACE_CANDIDATES_H
#define LIS_NULLSPACE_CANDIDATES_H

#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Bind accepted near-null directions to the solver. */
#define LIS_NULLSPACE_CANDIDATE_BIND_NEAR          0x01
/*
 * Also test right-exact directions against A^H and bind those that are
 * exact on the left. This flag does not accept a separate left-candidate
 * set; it only derives a left basis from the same candidates supplied for
 * right-side classification.
 */
#define LIS_NULLSPACE_CANDIDATE_BIND_LEFT          0x02
/* Include exact directions in the near-nullspace basis as well. */
#define LIS_NULLSPACE_CANDIDATE_NEAR_INCLUDE_EXACT 0x04

typedef struct LIS_NULLSPACE_CANDIDATE_DIAGNOSTICS_STRUCT
{
    LIS_INT input_count;
    LIS_INT independent_count;
    LIS_INT exact_count;
    LIS_INT near_count;
    LIS_INT left_exact_count;

    LIS_REAL operator_scale;
    LIS_REAL max_exact_defect;
    LIS_REAL max_near_defect;
}
LIS_NULLSPACE_CANDIDATE_DIAGNOSTICS;

typedef struct LIS_COMPONENT_NULLSPACE_DIAGNOSTICS_STRUCT
{
    LIS_INT component_count;
    LIS_INT localized_candidate_count;
    LIS_NULLSPACE_CANDIDATE_DIAGNOSTICS classifier;
}
LIS_COMPONENT_NULLSPACE_DIAGNOSTICS;

/*
 * Classify application-provided candidate vectors into exact and near
 * nullspace directions.
 *
 * The candidates are orthogonalized within each group, transformed by A,
 * and rotated by a small one-sided Jacobi SVD-like process so that their
 * scaled defects can be tested independently. groups may be NULL, in which
 * case all candidates belong to one group. Otherwise groups[i] identifies
 * the independent candidate group for candidates[i].
 *
 * operator_scale is an application-supplied positive scale for A. The
 * classification defects are ||A q||_2 / operator_scale. Therefore
 * exact_tolerance and near_tolerance are dimensionless relative tolerances,
 * not absolute matrix norms. near_tolerance must be at least
 * exact_tolerance.
 *
 * Exact directions are bound with lis_solver_set_nullspace(). When
 * LIS_NULLSPACE_CANDIDATE_BIND_NEAR is set, accepted near directions are
 * bound with lis_solver_set_near_nullspace(). When
 * LIS_NULLSPACE_CANDIDATE_BIND_LEFT is set, each right-exact direction is
 * additionally tested with ||A^H q||_2 / operator_scale and the accepted
 * subset is bound as the left nullspace.
 */
LIS_INT lis_solver_classify_nullspace_candidates(
    LIS_MATRIX A,
    LIS_SOLVER solver,
    LIS_INT nvec,
    LIS_VECTOR candidates[],
    const LIS_INT groups[],
    LIS_REAL operator_scale,
    LIS_REAL exact_tolerance,
    LIS_REAL near_tolerance,
    LIS_INT flags,
    LIS_NULLSPACE_CANDIDATE_DIAGNOSTICS *diagnostics);

/*
 * Split global candidates over weakly connected matrix components and pass
 * the resulting nonzero component-local candidates through the generic
 * classifier above.  This is useful when applications know a small global
 * candidate basis but disconnected matrix components must receive independent
 * nullspace directions.
 *
 * The initial implementation is intentionally serial, matching
 * lis_matrix_connected_components().
 */
LIS_INT lis_solver_classify_component_nullspace_candidates(
    LIS_MATRIX A,
    LIS_SOLVER solver,
    LIS_INT nvec,
    LIS_VECTOR candidates[],
    LIS_REAL graph_drop_tolerance,
    LIS_REAL operator_scale,
    LIS_REAL exact_tolerance,
    LIS_REAL near_tolerance,
    LIS_INT flags,
    LIS_COMPONENT_NULLSPACE_DIAGNOSTICS *diagnostics);

/* Return the dimension of the exact right-nullspace currently bound. */
LIS_INT lis_solver_get_nullspace_dim(
    LIS_SOLVER solver,
    LIS_INT *nvec);

/* Copy one exact right-nullspace vector into an application vector. */
LIS_INT lis_solver_copy_nullspace_vector(
    LIS_SOLVER solver,
    LIS_INT index,
    LIS_VECTOR vector);

#ifdef __cplusplus
}
#endif

#endif
