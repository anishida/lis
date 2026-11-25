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
#ifdef __cplusplus
}
#endif
#endif
