/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifndef LIS_RCI_RECYCLE_RESTORE_H
#define LIS_RCI_RECYCLE_RESTORE_H
#include "lis_rci_restore.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Recycle-space component for later augmented Krylov integration.
 * Fixed mode: U, ZU=M^-1 U, C=A ZU. Flexible: W, physical Z, C=A Z.
 * All three bases receive the same orthogonalization transforms.
 * Storage is real double, column-major, owned by the session.
 * New bases are staged; errors/cancel preserve the last committed basis.
 * Does not select harmonic Ritz coefficients or implement a Krylov method. */
typedef struct lis_rci_recycle_restore lis_rci_recycle_restore;
typedef struct {
 size_t n,count;
 uint64_t operator_generation,precon_generation;
 const double *domain,*physical,*image;
} lis_rci_recycle_basis;
typedef struct {
 uint64_t token,operator_generation,precon_generation;
 uint32_t operation; /* LIS_RCI_OP_PRECON or LIS_RCI_OP_A */
 size_t n;
 const double *in;
 double *out;
} lis_rci_recycle_request;
lis_rci_recycle_restore *lis_rci_recycle_restore_create(size_t n,size_t capacity,int flexible,double rank_tol);
void lis_rci_recycle_restore_destroy(lis_rci_recycle_restore *);
/* Fixed mode ignores physical (must be NULL); flexible requires physical.
 * Inputs have n*count values and are copied. */
int lis_rci_recycle_restore_begin(lis_rci_recycle_restore *,size_t count,const double *domain,const double *physical,
 uint64_t a,uint64_t p,lis_rci_recycle_request *);
/* Refresh C against a changed A; fixed mode reapplies M^-1.
 * A changed preconditioner generation in flexible mode is conservatively
 * rejected, without discarding the previous basis. Use begin for a new basis.
 */
int lis_rci_recycle_restore_refresh(lis_rci_recycle_restore *,uint64_t a,uint64_t p,lis_rci_recycle_request *);
int lis_rci_recycle_restore_advance(lis_rci_recycle_restore *,uint64_t token,uint64_t a,uint64_t p,int callback_code,lis_rci_recycle_request *);
int lis_rci_recycle_restore_cancel(lis_rci_recycle_restore *,uint64_t token);
int lis_rci_recycle_restore_error(const lis_rci_recycle_restore *,int *);
int lis_rci_recycle_restore_get(const lis_rci_recycle_restore *,lis_rci_recycle_basis *);
/* Apply x += Z*(C^T r), r -= C*(C^T r). The caller supplies true r=b-Ax.
 * Apply the same projection to a candidate direction w,z,az, with az=A z:
 * w -= W*(C^T az), z -= Z*(C^T az), az -= C*(C^T az).
 * Disjoint arrays required; result is staged and committed only if finite. */
int lis_rci_recycle_restore_project(lis_rci_recycle_restore *,uint64_t a,uint64_t p,double *x,double *r);
int lis_rci_recycle_restore_direction(lis_rci_recycle_restore *,uint64_t a,uint64_t p,double *w,double *z,double *az);
/* Commit a known Arnoldi relation without extra MATVEC/PRECON. The caller
 * guarantees image=A*physical and, in fixed mode, physical=M^-1*domain.
 * All arrays are copied before any transformation. Rank failure preserves old basis. */
int lis_rci_recycle_restore_commit_images(lis_rci_recycle_restore *,size_t count,const double *domain,const double *physical,const double *image,uint64_t a,uint64_t p);
#ifdef __cplusplus
}
#endif
#endif
