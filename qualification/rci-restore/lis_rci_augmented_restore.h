/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifndef LIS_RCI_AUGMENTED_RESTORE_H
#define LIS_RCI_AUGMENTED_RESTORE_H
#include "lis_rci_recycle_restore.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Experimental fixed/flexible augmented GMRES sequence. Reuses the existing
 * RCI GMRES/FGMRES recurrence one cycle at a time. All large arrays allocated
 * at creation; native harmonic selection may allocate small dense workspace.
 * Selector sees column-major D^T D and D^T F (D=A Z, F=U or physical Z);
 * coefficient columns represent selected real directions or conjugate pairs.
 * Failure is nonfatal and preserves committed recycle state. */
typedef int (*lis_rci_harmonic_select)(void *,size_t,size_t,const double *,const double *,double *,size_t *);
int lis_rci_restore_harmonic_native(void *,size_t,size_t,const double *,const double *,double *,size_t *);
typedef struct {
 size_t n,restart,recycle,max_cycles;
 double rtol,atol,rank_tol;
 int flexible;
 lis_rci_harmonic_select select;
 void *user;
} lis_rci_augmented_restore_options;
typedef struct lis_rci_augmented_restore lis_rci_augmented_restore;
typedef struct {
 size_t cycles,iterations,matvecs,precons,builds,updates,rejects,recycle_uses,recycle_rank;
 long double residual;
 int status,callback_code;
} lis_rci_augmented_restore_info;
lis_rci_augmented_restore *lis_rci_augmented_restore_create(const lis_rci_augmented_restore_options *);
void lis_rci_augmented_restore_destroy(lis_rci_augmented_restore *);
int lis_rci_augmented_restore_start(lis_rci_augmented_restore *,const double *b,const double *x0,uint64_t a,uint64_t p,lis_rci_recycle_request *);
int lis_rci_augmented_restore_advance(lis_rci_augmented_restore *,uint64_t token,uint64_t a,uint64_t p,int callback_code,lis_rci_recycle_request *);
int lis_rci_augmented_restore_cancel(lis_rci_augmented_restore *,uint64_t token);
int lis_rci_augmented_restore_get_info(const lis_rci_augmented_restore *,lis_rci_augmented_restore_info *);
const double *lis_rci_augmented_restore_solution(const lis_rci_augmented_restore *);
int lis_rci_augmented_restore_basis(const lis_rci_augmented_restore *,lis_rci_recycle_basis *);
#ifdef __cplusplus
}
#endif
#endif
