/* AI-assisted native continuation candidate, LIS COPYING terms proposed
 * subject to contributor confirmation. Qualification only. */
#ifndef LIS_RCI_NATIVE_FEAST_H
#define LIS_RCI_NATIVE_FEAST_H
#include "lis.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct lis_rci_native_feast lis_rci_native_feast;
enum { NF_DONE=0,NF_A=1,NF_B=2,NF_SHIFT=4,NF_ERROR=-1,NF_PROTOCOL=-2,
 NF_CANCELLED=-3,NF_NONFINITE=-4,NF_STATE=-5,NF_BADARG=-6,NF_STALE=-7 };
typedef struct {
 uint64_t token;
 int operation;
 LIS_MATRIX A,B;
 LIS_INT count;
 int block; /* Preserve the native scalar vs block dispatch even for count=1. */
 LIS_VECTOR *in,*out;
 LIS_REAL zr,zi;
 LIS_VECTOR in_real,in_imag,out_real,out_imag;
 LIS_SOLVER solver;
} lis_rci_native_feast_request;
/* Borrow the fully prepared native eigensolver until terminal.
 * A/B requests use count vectors (including scalar count=1).
 * Shift requests solve (zr+i*zi)B-A in split-complex representation.
 * NULL B means identity. In strict mode all requested output entries
 * must be finite and written. A/B identity and generation counters must
 * remain unchanged until terminal; destruction is forbidden while borrowed.
 * Cancellation exposes no certified eigenpairs.
 * This is not an installed API. Real-double non-MPI only. */
lis_rci_native_feast *lis_rci_native_feast_create(void);
void lis_rci_native_feast_destroy(lis_rci_native_feast *);
int lis_rci_native_feast_start(lis_rci_native_feast *,LIS_ESOLVER,int,lis_rci_native_feast_request *);
int lis_rci_native_feast_advance(lis_rci_native_feast *,uint64_t,LIS_INT,lis_rci_native_feast_request *);
int lis_rci_native_feast_cancel(lis_rci_native_feast *,uint64_t);
LIS_INT lis_rci_native_feast_native(const lis_rci_native_feast *);
typedef LIS_INT (*lis_rci_native_feast_test_driver)(LIS_ESOLVER);
void lis_rci_native_feast_test_select(int,lis_rci_native_feast_test_driver);
#ifdef __cplusplus
}
#endif
#endif
