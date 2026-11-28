/* AI-assisted continuation reconstruction from the actual native LIS sources.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifndef LIS_RCI_NATIVE_LINEAR_RESTORE_H
#define LIS_RCI_NATIVE_LINEAR_RESTORE_H
#include "lis_rci_restore.h"
#include "lis.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Qualification-only native continuation. The solver is already prepared
 * by the normal lis_solve_kernel setup, including preconditioner, vectors,
 * convergence policy and history. It is exclusively borrowed until terminal.
 * No full solve setup/teardown is performed here. All native vector algebra
 * and original stopping/monitor policies are preserved. Callback buffers
 * are native LIS_VECTOR handles, not owned by this session.
 * Request codes A/PRECON use native LIS_INT completion codes. Native
 * status is available separately; cancellation can leave a partial iterate.
 * The classical driver preserves native behavior; strict external mode
 * additionally requires every output entry to be filled and finite.
 * No MPI or double-double solver mode. Native scalar precision is determined by LIS.
 */
typedef struct lis_rci_native_linear_restore lis_rci_native_linear_restore;
typedef struct {
 uint64_t token;
 uint32_t operation;
 LIS_MATRIX A;
 LIS_SOLVER solver;
 LIS_VECTOR in,out;
} lis_rci_native_linear_request;
lis_rci_native_linear_restore *lis_rci_native_linear_restore_create(void);
void lis_rci_native_linear_restore_destroy(lis_rci_native_linear_restore *);
int lis_rci_native_linear_restore_start(lis_rci_native_linear_restore *,LIS_SOLVER,int strict,lis_rci_native_linear_request *);
/* Explicit method for an already prepared native solver. The method must
 * match the prepared work layout, including native FGCRO-DR -> FGMRES
 * fallback. It does not alter solver->options. */
int lis_rci_native_linear_restore_start_method(lis_rci_native_linear_restore *,LIS_SOLVER,LIS_INT method,int strict,lis_rci_native_linear_request *);
int lis_rci_native_linear_restore_advance(lis_rci_native_linear_restore *,uint64_t,LIS_INT,lis_rci_native_linear_request *);
int lis_rci_native_linear_restore_cancel(lis_rci_native_linear_restore *,uint64_t);
int lis_rci_native_linear_restore_status(const lis_rci_native_linear_restore *,LIS_INT *,uint64_t *,uint64_t *);
LIS_INT lis_rci_native_linear_restore_classic(LIS_SOLVER);
/* Qualification dispatch hook, per thread, called after native setup. */
typedef LIS_INT (*lis_rci_native_linear_test_driver)(LIS_SOLVER,LIS_INT method);
void lis_rci_native_linear_test_select(int reference,lis_rci_native_linear_test_driver driver);
#ifdef __cplusplus
}
#endif
#endif
