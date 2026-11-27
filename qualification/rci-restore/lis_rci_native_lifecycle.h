/* Experimental full-lifecycle linear RCI, not an installed API. */
#ifndef LIS_RCI_NATIVE_LIFECYCLE_H
#define LIS_RCI_NATIVE_LIFECYCLE_H
#include "lis_rci_native_session.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct lis_rci_native_lifecycle lis_rci_native_lifecycle;
lis_rci_native_lifecycle *lis_rci_native_lifecycle_create(void);
void lis_rci_native_lifecycle_destroy(lis_rci_native_lifecycle *);
/* Borrow assembled unscaled CSR A, b, x and a normal configured solver.
 * NONE/Jacobi, no scaling/conversion/nullspace/ADDS/thread override.
 * Synchronous native setup may allocate; all iterative A/P and final
 * true residual A requests are external. x is published only on DONE.
 * Same calling thread and OpenMP setting for the entire active solve.
 * A/b/x/solver must remain alive, exclusive and unchanged while borrowed.
 * The matrix generation is checked on advance. Worker must stop before
 * cancellation/destroy. No implicit use of the old array RCI API. */
int lis_rci_native_lifecycle_start(lis_rci_native_lifecycle *,LIS_MATRIX,LIS_VECTOR,LIS_VECTOR,LIS_SOLVER,lis_rci_native_session_request *);
int lis_rci_native_lifecycle_advance(lis_rci_native_lifecycle *,uint64_t,LIS_INT,lis_rci_native_session_request *);
int lis_rci_native_lifecycle_cancel(lis_rci_native_lifecycle *,uint64_t);
int lis_rci_native_lifecycle_status(const lis_rci_native_lifecycle *,LIS_INT *);
#ifdef __cplusplus
}
#endif
#endif
