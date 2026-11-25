/* Optional bridge from public LIS Matrix/Vector APIs to portable double CSR core.
 * Serial real-double only. No automatic activation; no changes to MFAMG defaults.
 */
#ifndef LIS_CONDENSE_LIS_H
#define LIS_CONDENSE_LIS_H
#include "lis.h"
#include "lis_condense_core.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct lis_sc_lis_session lis_sc_lis_session;
/* Creates a fully owned condensed matrix. Outputs allocated session on success only. */
lis_sc_status lis_sc_lis_create(LIS_MATRIX A,LIS_INT nblocks,
                                 const LIS_INT *block_ptr,const LIS_INT *block_dof,
                                 double pivot_tolerance,lis_sc_lis_session **out);
/* Atomic update of an existing LIS session, with caller-supplied full matrix.
 * FULL: allows any new CSR sparsity, dimensions and block mapping.
 * NUMERIC: requires identical full CSR pattern, block map and pivot tolerance.
 * REUSE_ALL: requires bitwise identical full CSR values as well.
 * The new LIS Schur matrix is assembled *before* committing the new plan.
 * A failing update leaves lis_sc_lis_matrix(session) and every core object intact.
 * The caller MUST destroy/rebuild any solver or preconditioner built on the old
 * Schur matrix before calling this API; old borrowed matrix pointers expire
 * after a successful FULL or NUMERIC update.
 */
lis_sc_status lis_sc_lis_update(lis_sc_lis_session *session, LIS_MATRIX A,
                                 LIS_INT nblocks, const LIS_INT *block_ptr,
                                 const LIS_INT *block_dof, double pivot_tolerance,
                                 lis_sc_reuse_mode mode);
void lis_sc_lis_destroy(lis_sc_lis_session *session);
/* Matrix remains owned by session and MUST NOT be destroyed by the caller. */
LIS_MATRIX lis_sc_lis_matrix(const lis_sc_lis_session *session);
lis_sc_status lis_sc_lis_rhs(const lis_sc_lis_session *session,
                              LIS_VECTOR b_full,LIS_VECTOR g_reduced);
lis_sc_status lis_sc_lis_reconstruct(const lis_sc_lis_session *session,
                                      LIS_VECTOR u_reduced,LIS_VECTOR b_full,
                                      LIS_VECTOR x_full);
lis_sc_status lis_sc_lis_residual(const lis_sc_lis_session *session,
                                   LIS_VECTOR b_full,LIS_VECTOR x_full,double *rr);
/* Caller supplies element counts for both buffers; invalid lengths leave the
 * destination untouched. These are counts of LIS_INT, not byte lengths. */
lis_sc_status lis_sc_lis_labels(const lis_sc_lis_session *session,
                                 const LIS_INT *full_labels,size_t full_label_count,
                                 LIS_INT *physical_labels,size_t physical_label_count);
/* Project caller-owned full-space LIS vectors to preallocated reduced vectors.
 * The function does not allocate or take ownership of either vector array.
 * No assumption is made about the number or physical meaning of the modes.
 */
lis_sc_status lis_sc_lis_modes(const lis_sc_lis_session *session,
                                LIS_INT nvec,
                                LIS_VECTOR const full_modes[],
                                LIS_VECTOR const reduced_modes[]);
/* Managed serial real-double solve on the owned Schur operator.
 * The session owns condensation, solver and work vectors. It uses sequence
 * solves for repeated RHS. Configure the borrowed solver, but never destroy
 * it or attach an externally owned preconditioner. The borrowed low-level
 * condensation is for read-only mode/label projection; never update or
 * destroy its matrix independently.
 * FULL updates require explicit clearing of dimension-dependent basis and
 * partition metadata (REBIND_REQUIRED); NUMERIC retains geometric settings
 * and the reduced warm start, but requires clearing exact/left nullspaces.
 * Reclassify mathematical nullspaces against the changed operator.
 * Updates stage new plan, CSR and workspace before retiring old solver state.
 * FULL and NUMERIC retire the preconditioner; REUSE_ALL preserves it.
 * The plan owns a snapshot of the full operator. A caller modification
 * becomes visible only after a successful update; certification is against
 * that snapshot, not against a subsequently modified caller matrix.
 * REUSE_ALL never retires the solver. On solve failure x_full is unchanged.
 * Only unscaled, default-precision, non-transposed CSR solves are supported.
 */
typedef struct lis_sc_solver_session lis_sc_solver_session;
typedef struct {
    lis_sc_status status;
    LIS_INT solver_called;
    LIS_INT api_return;
    LIS_INT solver_status;
    LIS_INT iterations;
    LIS_INT full_residual_valid;
    LIS_INT solution_committed;
    double full_relative_residual;
} lis_sc_solver_report;
lis_sc_status lis_sc_solver_session_create(
    LIS_MATRIX A,LIS_INT nblocks,const LIS_INT *block_ptr,
    const LIS_INT *block_dof,double pivot_tolerance,lis_sc_solver_session **out);
lis_sc_status lis_sc_solver_session_update(
    lis_sc_solver_session *session,LIS_MATRIX A,LIS_INT nblocks,
    const LIS_INT *block_ptr,const LIS_INT *block_dof,double pivot_tolerance,
    lis_sc_reuse_mode mode);
LIS_SOLVER lis_sc_solver_session_get_solver(lis_sc_solver_session *session);
const lis_sc_lis_session *lis_sc_solver_session_get_condensation(
    const lis_sc_solver_session *session);
/* The previous successful reduced iterate is available to the configured
 * solver as an initial guess when -initx_zeros false is selected. */
lis_sc_status lis_sc_solver_session_solve(
    lis_sc_solver_session *session,LIS_VECTOR b_full,LIS_VECTOR x_full,
    double full_residual_tolerance);
/* Report exists after a solve attempt, and is invalidated on matrix update. */
lis_sc_status lis_sc_solver_session_get_report(
    const lis_sc_solver_session *session,lis_sc_solver_report *report);
void lis_sc_solver_session_destroy(lis_sc_solver_session *session);

#ifdef __cplusplus
}
#endif
#endif
