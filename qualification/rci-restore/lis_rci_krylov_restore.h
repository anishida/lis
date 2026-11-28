/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifndef LIS_RCI_KRYLOV_RESTORE_H
#define LIS_RCI_KRYLOV_RESTORE_H
#include "lis_rci_restore.h"
#include "lis_krylov.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Qualification-only C11 continuation of the actual native symmetric core.
 * No worker, stack switching or replay. Problem/options/extension, optional
 * workspace, and result arrays are borrowed and must remain valid, exclusive
 * and unchanged until completion. Non-null native callback pointers are
 * required by the original argument validation, but never invoked by RCI.
 * SHIFTED_SOLVE uses the native convention (A-sigma B)^-1, NOT the FEAST
 * projector convention (sigma B-A)^-1. Caller prepares its inverse.
 * DONE guarantees native SUCCESS; MAXITER represents NO_CONVERGENCE or
 * BUDGET_EXHAUSTED (inspect native status). Native partial results retained
 * only on NO_CONVERGENCE. Cancel/error invalidates result counts.
 */
typedef struct lis_rci_krylov_restore lis_rci_krylov_restore;
typedef struct {
 uint64_t token;
 uint32_t operation;
 int64_t n;
 const double *in;
 double *out;
 double sigma;
} lis_rci_krylov_restore_request;
lis_rci_krylov_restore *lis_rci_krylov_restore_create(void);
void lis_rci_krylov_restore_destroy(lis_rci_krylov_restore *);
int lis_rci_krylov_restore_start(lis_rci_krylov_restore *,const lis_krylov_problem64 *,const lis_krylov_options *,
 const lis_krylov_ext_options *,lis_krylov_workspace *,const lis_krylov_budget *,
 double *,double *,double *,lis_krylov_stats *,lis_krylov_ext_stats *,lis_rci_krylov_restore_request *);
int lis_rci_krylov_restore_advance(lis_rci_krylov_restore *,uint64_t,int,lis_rci_krylov_restore_request *);
int lis_rci_krylov_restore_cancel(lis_rci_krylov_restore *,uint64_t);
int lis_rci_krylov_restore_status(const lis_rci_krylov_restore *,int *,lis_krylov_budget_stats *);
/* Blocking driver over the very same continuation core, using the native
 * callback ABI. Returns native LIS_KRYLOV_* status. The existing classical
 * entry points are not redirected by this qualification experiment. */
int lis_rci_krylov_restore_classic64(const lis_krylov_problem64 *,const lis_krylov_options *,
 const lis_krylov_ext_options *,lis_krylov_workspace *,const lis_krylov_budget *,
 double *,double *,double *,lis_krylov_stats *,lis_krylov_ext_stats *,lis_krylov_budget_stats *);
#ifdef LIS_RCI_KRYLOV_NATIVE_LIBRARY
/* Test-only direct original reference, available in the patched native build. */
int lis_rci_krylov_restore_reference64(const lis_krylov_problem64 *,const lis_krylov_options *,
 const lis_krylov_ext_options *,lis_krylov_workspace *,double *,double *,double *,
 lis_krylov_stats *,lis_krylov_ext_stats *);
#endif
#ifdef __cplusplus
}
#endif
#endif
