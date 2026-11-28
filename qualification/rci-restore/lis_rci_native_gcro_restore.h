/* AI-assisted continuation reconstruction from the actual native LIS sources.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifndef LIS_RCI_NATIVE_GCRO_RESTORE_H
#define LIS_RCI_NATIVE_GCRO_RESTORE_H
#include "lis_rci_native_linear_restore.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Borrow a prepared native GCRO-DR solver and all of its vectors exclusively.
 * Native real-double, non-MPI, no double-double mode. Unlike the separate
 * transactional generic recycle component, native failure semantics clear
 * the possibly modified recycle basis. Cancellation unwinds the original
 * failure/cleanup path without executing additional callbacks. */
typedef struct lis_rci_native_gcro_restore lis_rci_native_gcro_restore;
lis_rci_native_gcro_restore *lis_rci_native_gcro_restore_create(void);
void lis_rci_native_gcro_restore_destroy(lis_rci_native_gcro_restore *);
int lis_rci_native_gcro_restore_start(lis_rci_native_gcro_restore *,LIS_SOLVER,int,lis_rci_native_linear_request *);
int lis_rci_native_gcro_restore_advance(lis_rci_native_gcro_restore *,uint64_t,LIS_INT,lis_rci_native_linear_request *);
int lis_rci_native_gcro_restore_cancel(lis_rci_native_gcro_restore *,uint64_t);
int lis_rci_native_gcro_restore_status(const lis_rci_native_gcro_restore *,LIS_INT *,uint64_t *,uint64_t *);
typedef LIS_INT (*lis_rci_native_gcro_test_driver)(LIS_SOLVER);
void lis_rci_native_gcro_test_select(int,lis_rci_native_gcro_test_driver);
#ifdef __cplusplus
}
#endif
#endif
