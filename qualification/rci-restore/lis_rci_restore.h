/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifndef LIS_RCI_RESTORE_PROTOCOL_H
#define LIS_RCI_RESTORE_PROTOCOL_H
#include "lis_rci.h"
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Separate experimental ABI; never extend an existing public struct in place.
 * A handle is an application-defined identity, not necessarily a host pointer.
 * The protocol never dereferences handles. Input handles are read-only.
 * The caller owns all vector objects until the session finishes. */
#define LIS_RCI_OPAQUE_VERSION 1u
typedef void *lis_rci_handle;
typedef struct { long double real, imag; } lis_rci_complex_coeff;
enum lis_rci_opaque_operation {
 LIS_RCI_OP_A=101, LIS_RCI_OP_B, LIS_RCI_OP_PRECON,
 LIS_RCI_OP_SHIFTED_SOLVE, LIS_RCI_OP_ZERO, LIS_RCI_OP_AXPBY
};
typedef struct {
 uint32_t version, size, operation, flags;
 uint64_t token, generation_a, generation_b;
 size_t columns;
 const lis_rci_handle *in, *in2, *out;
 lis_rci_complex_coeff shift, alpha, beta;
} lis_rci_opaque_request;
/* Return DONE on valid request, BADARG otherwise. bytes is the available
 * buffer size; validate before reading fields past the version/size prefix.
 * Reject unsupported flags, null handles and all output/input aliases.
 * Distinct handles must also refer to disjoint storage (caller obligation).
 * SHIFTED_SOLVE means (shift*B-A) X=R, with B=I in standard problems;
 * the backend is configured for the same problem as the projector.
 * AXPBY means out=alpha*in+beta*in2. Complex representation is backend-owned. */
int lis_rci_opaque_validate(const lis_rci_opaque_request *request,size_t bytes);
typedef int (*lis_rci_opaque_execute)(void *user,const lis_rci_opaque_request *request);
typedef int (*lis_rci_generation_read)(void *user,uint64_t *a,uint64_t *b);
typedef struct {
 void *user;
 lis_rci_opaque_execute block, single;
 lis_rci_generation_read generations;
} lis_rci_opaque_backend;
/* Prefer a block callback, otherwise dispatch each column with the same token.
 * No automatic retry on callback errors; earlier columns may be modified.
 * generation checks occur before and after every callback.
 * callback_code receives the first nonzero backend code. */
int lis_rci_opaque_dispatch(const lis_rci_opaque_backend *backend,
 const lis_rci_opaque_request *request,size_t bytes,int *callback_code);

/* Non-MPI contour projector. Caller supplies disjoint handles for five
 * column arrays: Q, R, X, Y0, Y1. All workspace and quadrature data are
 * allocated/copied at creation; start/advance do not allocate.
 * Q is read-only, R/X/Y0/Y1 are scratch. No result published before DONE.
 * For standard problems R=Q (B is identity); for generalized R=BQ.
 * Y=sum_j weight[j]*(shift[j]*B-A)^-1 R.
 * Neither the projector nor the dispatcher implements the linear solve.
 */
typedef struct lis_rci_projector_restore lis_rci_projector_restore;
typedef struct {
 size_t columns,nodes;
 int generalized;
 const lis_rci_complex_coeff *shift,*weight;
 const lis_rci_handle *q,*r,*x,*y0,*y1;
} lis_rci_projector_restore_options;
lis_rci_projector_restore *lis_rci_projector_restore_create(const lis_rci_projector_restore_options *);
void lis_rci_projector_restore_destroy(lis_rci_projector_restore *);
int lis_rci_projector_restore_start(lis_rci_projector_restore *,uint64_t a,uint64_t b,lis_rci_opaque_request *);
int lis_rci_projector_restore_advance(lis_rci_projector_restore *,uint64_t token,uint64_t a,uint64_t b,int callback_code,lis_rci_opaque_request *);
int lis_rci_projector_restore_cancel(lis_rci_projector_restore *,uint64_t token);
int lis_rci_projector_restore_error(const lis_rci_projector_restore *,int *callback_code);
const lis_rci_handle *lis_rci_projector_restore_result(const lis_rci_projector_restore *);

/* Synchronous nested runner of the existing real-double linear core.
 * Uses borrowed double-array handles only for this adapter, not a new solver.
 * Parent request must be a single-column shifted solve. Caller owns/validates
 * its actual parent context and must advance it only after successful return.
 * The backend interprets OP_A as the parent shifted operator.
 * Child failure leaves output unchanged; child_code reports backend error.
 */
int lis_rci_nested_restore_run(lis_rci_context *child,
 const lis_rci_opaque_request *parent,size_t bytes,size_t n,
 const lis_rci_opaque_backend *backend,int *child_code);

/* Split-complex adapter: each opaque handle identifies two disjoint real
 * arrays of length n. Only this backend dereferences those vector handles.
 * A and B are real operators; B is identity for standard problems.
 * Shifted solves use a persistent doubled-real core (BiCGSTAB/GMRES/FGMRES),
 * identity preconditioning, and real MATVEC callbacks. No new solver kernel.
 * Apply must completely fill out, must not modify in, and return 0 on success.
 * This backend is not reentrant. A user-defined callback code is retained in
 * split_error; backend operation failures are otherwise returned as RCI status.
 */
typedef struct { size_t n; double *real,*imag; } lis_rci_split_vector;
typedef int (*lis_rci_real_apply)(void *,uint32_t,size_t,const double *,double *);
typedef struct lis_rci_split_restore lis_rci_split_restore;
lis_rci_split_restore *lis_rci_split_restore_create(size_t n,int generalized,
 const lis_rci_options *child_options,void *user,lis_rci_real_apply apply,
 lis_rci_generation_read generations);
void lis_rci_split_restore_destroy(lis_rci_split_restore *);
int lis_rci_split_restore_backend(lis_rci_split_restore *,lis_rci_opaque_backend *);
int lis_rci_split_restore_error(const lis_rci_split_restore *,int *callback_code);
#ifdef __cplusplus
}
#endif
#endif
