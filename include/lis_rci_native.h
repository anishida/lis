/* Experimental native RCI API; proposed under LIS COPYING, subject to contributor confirmation.
 * Native LIS scalar, serial/OpenMP, no MPI or double-double solver mode. Not a stable ABI. */
#ifndef LIS_RCI_NATIVE_H
#define LIS_RCI_NATIVE_H
#include "lis.h"
#include "lis_rci_status.h"
#include "lis_krylov.h"
#include "lis_abi.h"
#ifdef __cplusplus
extern "C" {
#endif
#define LIS_RCI_NATIVE_API_VERSION 1u
enum { LIS_NR_OP_A=101,LIS_NR_OP_B=102,LIS_NR_OP_PRECON=103,LIS_NR_OP_SHIFTED_SOLVE=104 };
enum { LIS_NR_LINEAR=1,LIS_NR_GCRO=2,LIS_NR_FGCRO=3,LIS_NR_FEAST=4,LIS_NR_KRYLOV=5 };
enum { LIS_NR_VECTOR=1,LIS_NR_SPLIT=2,LIS_NR_RAW64=3 };
enum { LIS_NR_STALE=-9 };
/* Same operation codes as the opaque Generic RCI protocol. Vector arrays
 * are borrowed native handles; raw64 arrays are borrowed double buffers.
 * Shift sign explicitly identifies sigma B-A (+1) versus A-sigma B (-1).
 * Scalar/block dispatch is preserved even for a one-column block.
 * The solver/eigensolver and operands are exclusively borrowed until terminal.
 * Context kind is fixed, tokens survive restarts, wrong tokens are retryable.
 * Output descriptors are unchanged on nonconsuming errors. */
typedef struct {
 uint64_t token;
 uint32_t operation,representation;
 LIS_INT count;
 int block,shift_sign;
 LIS_MATRIX A,B;
 LIS_SOLVER solver;
 LIS_VECTOR *in,*out;
 LIS_VECTOR in_real,in_imag,out_real,out_imag;
 LIS_REAL zr,zi;
 int64_t n64;
 const double *raw_in;
 double *raw_out;
 double sigma;
} lis_rci_native_session_request;

typedef struct lis_rci_native_context lis_rci_native_context;
/* Pass client-compiled ABI descriptor and request size, before using typed objects.
 * Valid kind: LIS_NR_LINEAR (all six linear solvers), LIS_NR_FEAST, LIS_NR_KRYLOV.
 * NULL means unsupported kind, ABI mismatch or allocation failure. */
lis_rci_native_context *lis_rci_native_create(int kind,uint32_t version,size_t request_size,const LIS_ABI_INFO *client,size_t client_size);
void lis_rci_native_destroy(lis_rci_native_context *);
/* Linear: assembled unscaled CSR, NONE/Jacobi, no conversion, nullspace,
 * ADDS or thread override. FEAST: Hermitian CSR A and optional SPD CSR B,
 * no conversion, automatic retry or explicit transform. Native setup is
 * synchronous; iterative A/B/P/shift and final linear residual requests
 * are external. x and FEAST value are committed only on DONE.
 * Krylov: fixed-double callback ABI; caller-owned result storage and
 * native partial-result semantics. All borrowed operands/options must
 * remain valid, exclusive and unchanged until terminal.
 * Same calling thread and OpenMP settings; no concurrent LIS calls.
 * Worker must stop accessing request buffers before cancel/destroy. */
int lis_rci_native_start_linear(lis_rci_native_context *,LIS_MATRIX,LIS_VECTOR,LIS_VECTOR,LIS_SOLVER,lis_rci_native_session_request *);
int lis_rci_native_start_feast(lis_rci_native_context *,LIS_MATRIX,LIS_MATRIX,LIS_VECTOR,LIS_SCALAR *,LIS_ESOLVER,lis_rci_native_session_request *);
int lis_rci_native_start_krylov(lis_rci_native_context *,const lis_krylov_problem64 *,const lis_krylov_options *,const lis_krylov_ext_options *,lis_krylov_workspace *,const lis_krylov_budget *,double *,double *,double *,lis_krylov_stats *,lis_krylov_ext_stats *,lis_rci_native_session_request *);
int lis_rci_native_advance(lis_rci_native_context *,uint64_t,LIS_INT,lis_rci_native_session_request *);
int lis_rci_native_cancel(lis_rci_native_context *,uint64_t);
int lis_rci_native_pending(const lis_rci_native_context *,lis_rci_native_session_request *);
int lis_rci_native_status(const lis_rci_native_context *,LIS_INT *native);
/* Optional synchronous LIS executor for vector and split requests only.
 * Completes output but does not advance. Returns native callback code.
 * For RAW64 the application must execute its own callback. */
LIS_INT lis_rci_native_execute_default(lis_rci_native_context *,uint64_t);
#ifdef __cplusplus
}
#endif
#endif
