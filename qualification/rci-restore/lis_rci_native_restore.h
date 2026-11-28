/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifndef LIS_RCI_NATIVE_RESTORE_H
#define LIS_RCI_NATIVE_RESTORE_H
#include "lis_rci_restore.h"
#include "lis_rci_bridge.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct lis_rci_native_restore lis_rci_native_restore;
/* Borrow assembled A, optional B, and optional prepared LIS_SOLVER with an
 * already-created preconditioner for A. The owner must keep them alive and
 * immutable during the adapter lifetime. Native preconditioners may require
 * other prepared solver state; the caller is responsible for preparation.
 * Generation changes invalidate this adapter permanently (recreate it).
 * Serial/OpenMP real-double, no MPI. Use matching LIS headers and flags.
 */
lis_rci_native_restore *lis_rci_native_restore_create(LIS_MATRIX A,LIS_MATRIX B,LIS_SOLVER prepared);
void lis_rci_native_restore_destroy(lis_rci_native_restore *);
int lis_rci_native_restore_generations(void *,uint64_t *,uint64_t *);
int lis_rci_native_restore_apply(void *,uint32_t,size_t,const double *,double *);
int lis_rci_native_restore_solve(lis_rci_native_restore *,lis_rci_context *,const double *,double *,lis_rci_info *);
#ifdef __cplusplus
}
#endif
#endif
