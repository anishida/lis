/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#ifndef LIS_RCI_LIS_BRIDGE_H
#define LIS_RCI_LIS_BRIDGE_H
#include "lis_rci.h"
#include <lis.h>
#include <lis_abi.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Check client headers against the linked library before ANY typed LIS call.
 * This does not translate layouts; pass compile-time sizes from the client. */
int lis_rci_bridge_check_abi(unsigned int_bytes,unsigned scalar_bytes,unsigned real_bytes,unsigned scalar_kind);
#define LIS_RCI_CHECK_CLIENT_ABI() lis_rci_bridge_check_abi((unsigned)sizeof(LIS_INT),(unsigned)sizeof(LIS_SCALAR),(unsigned)sizeof(LIS_REAL),LIS_ABI_SCALAR_REAL)
/* Non-MPI real-double LIS (serial or OpenMP). A is borrowed, assembled and
 * must remain valid for the bridge lifetime and unchanged during solve.
 * Do not use the same bridge concurrently. LIS calls inherit LIS thread-safety
 * restrictions. BADARG includes invalid input or allocation failure at create. */
typedef struct lis_rci_bridge lis_rci_bridge;
/* jacobi=0 identity, jacobi=1 diagonal inverse; CG requires positive diagonal. */
int lis_rci_bridge_create(LIS_MATRIX A,const lis_rci_options *opt,int jacobi,lis_rci_bridge **out);
void lis_rci_bridge_destroy(lis_rci_bridge *p);
/* b/x are arrays of opt.n doubles. x initial guess and result (on successful solve only). */
int lis_rci_bridge_solve(lis_rci_bridge *p,const double *b,double *x,lis_rci_info *info);
#ifdef __cplusplus
}
#endif
#endif
