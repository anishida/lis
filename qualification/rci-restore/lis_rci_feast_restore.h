/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifndef LIS_RCI_FEAST_RESTORE_H
#define LIS_RCI_FEAST_RESTORE_H
#include "lis_rci_restore.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Real symmetric A and SPD B (or identity), real-double storage, conjugate
 * quadrature. The outer machine uses split-complex handles for A/B and the
 * contour projector. The small projected solve is supplied by the caller
 * (e.g. lis_eigen_projected_solve). It overwrites row-major H and returns
 * eigenvectors as columns of the row-major coefficient matrix.
 * Only DONE publishes results. Stable count/trace and physical residual
 * tests do not certify the absence of missing eigenvalues in the interval.
 * No classical FEAST options or ABI compatibility are implied. */
typedef int (*lis_rci_feast_dense)(void *,size_t,double *,double *,double *);
typedef struct {
 size_t n,columns,nodes,max_cycles;
 int generalized;
 double lo,hi,residual_tol,trace_tol,rank_tol,imag_tol;
 const lis_rci_complex_coeff *shift,*weight;
 lis_rci_feast_dense projected_solve;
 void *user;
} lis_rci_feast_restore_options;
typedef struct lis_rci_feast_restore lis_rci_feast_restore;
typedef struct {
 size_t count,cycles,rank;
 const double *values,*vectors,*residuals; /* vectors column-major n*count */
} lis_rci_feast_restore_result;
lis_rci_feast_restore *lis_rci_feast_restore_create(const lis_rci_feast_restore_options *);
void lis_rci_feast_restore_destroy(lis_rci_feast_restore *);
/* initial is n*columns real values in column-major order, copied at start. */
int lis_rci_feast_restore_start(lis_rci_feast_restore *,const double *initial,uint64_t a,uint64_t b,lis_rci_opaque_request *);
int lis_rci_feast_restore_advance(lis_rci_feast_restore *,uint64_t token,uint64_t a,uint64_t b,int callback_code,lis_rci_opaque_request *);
int lis_rci_feast_restore_cancel(lis_rci_feast_restore *,uint64_t token);
int lis_rci_feast_restore_error(const lis_rci_feast_restore *,int *);
int lis_rci_feast_restore_get_result(const lis_rci_feast_restore *,lis_rci_feast_restore_result *);
#ifdef __cplusplus
}
#endif
#endif
