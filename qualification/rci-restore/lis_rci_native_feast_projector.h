/* AI-assisted native continuation candidate; LIS COPYING terms proposed,
 * subject to contributor confirmation. Qualification-only, real double. */
#ifndef LIS_RCI_NATIVE_FEAST_PROJECTOR_H
#define LIS_RCI_NATIVE_FEAST_PROJECTOR_H
#include "lis.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct lis_rci_native_feast_projector lis_rci_native_feast_projector;
enum { FP_DONE=0, FP_B=1, FP_SHIFT=2, FP_ERROR=-1, FP_PROTOCOL=-2,
       FP_CANCELLED=-3, FP_NONFINITE=-4, FP_STATE=-5, FP_BADARG=-6 };
typedef struct {
 uint64_t token;
 int operation;
 LIS_MATRIX A,B;
 LIS_REAL zr,zi;
 LIS_VECTOR in_real,in_imag,out_real,out_imag;
 LIS_SOLVER solver;
} lis_rci_native_feast_projector_request;
/* A/B, solver and input/output vectors borrowed exclusively until terminal.
 * The output buffers of requests are write-only; external strict mode
 * poisons them to detect incomplete replies. Shift callbacks must solve
 * (zr+i*zi)B-A, with B=NULL meaning identity. On error the projector output
 * may be partial, as in the original LIS function. */
lis_rci_native_feast_projector *lis_rci_native_feast_projector_create(void);
void lis_rci_native_feast_projector_destroy(lis_rci_native_feast_projector *);
int lis_rci_native_feast_projector_start(lis_rci_native_feast_projector *,
 LIS_MATRIX,LIS_MATRIX,LIS_REAL,LIS_REAL,LIS_REAL,LIS_INT,LIS_VECTOR,LIS_VECTOR,
 LIS_VECTOR,LIS_VECTOR,LIS_SOLVER,int,lis_rci_native_feast_projector_request *);
int lis_rci_native_feast_projector_advance(lis_rci_native_feast_projector *,uint64_t,LIS_INT,lis_rci_native_feast_projector_request *);
int lis_rci_native_feast_projector_cancel(lis_rci_native_feast_projector *,uint64_t);
LIS_INT lis_rci_native_feast_projector_native(const lis_rci_native_feast_projector *);
#ifdef __cplusplus
}
#endif
#endif
