/* Experimental full FEAST lifecycle. Not an installed ABI. */
#ifndef LIS_RCI_NATIVE_FEAST_LIFECYCLE_H
#define LIS_RCI_NATIVE_FEAST_LIFECYCLE_H
#include "lis_rci_native_session.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct lis_rci_native_feast_lifecycle lis_rci_native_feast_lifecycle;
lis_rci_native_feast_lifecycle *lis_rci_native_feast_lifecycle_create(void);
void lis_rci_native_feast_lifecycle_destroy(lis_rci_native_feast_lifecycle *);
/* Borrow A, optional SPD B, output x, eigenvalue and configured FEAST solver
 * exclusively until terminal. Only native scalar Hermitian CSR, no conversion, retry
 * or explicit spectral transform. A/B/shift requests are all external.
 * The result arrays inside esolver are scratch while active; never inspect
 * them until terminal. x and *value are committed only on DONE; failure
 * leaves nev zero and invalidates current results. Same thread/settings. */
int lis_rci_native_feast_lifecycle_start(lis_rci_native_feast_lifecycle *,LIS_MATRIX,LIS_MATRIX,LIS_VECTOR,LIS_SCALAR *,LIS_ESOLVER,lis_rci_native_session_request *);
int lis_rci_native_feast_lifecycle_advance(lis_rci_native_feast_lifecycle *,uint64_t,LIS_INT,lis_rci_native_session_request *);
int lis_rci_native_feast_lifecycle_cancel(lis_rci_native_feast_lifecycle *,uint64_t);
int lis_rci_native_feast_lifecycle_status(const lis_rci_native_feast_lifecycle *,LIS_INT *);
#ifdef __cplusplus
}
#endif
#endif
