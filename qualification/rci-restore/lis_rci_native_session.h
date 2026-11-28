/* Candidate typed native Generic RCI session; not an installed ABI. */
#ifndef LIS_RCI_NATIVE_SESSION_H
#define LIS_RCI_NATIVE_SESSION_H
#include "lis_rci_restore.h"
#include "lis.h"
#include "lis_krylov.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct lis_rci_native_session lis_rci_native_session;
#include "lis_rci_native.h"
lis_rci_native_session *lis_rci_native_session_create(int kind);
void lis_rci_native_session_destroy(lis_rci_native_session *);
int lis_rci_native_session_start_linear(lis_rci_native_session *,LIS_SOLVER,LIS_INT method,lis_rci_native_session_request *);
int lis_rci_native_session_start_feast(lis_rci_native_session *,LIS_ESOLVER,lis_rci_native_session_request *);
int lis_rci_native_session_start_krylov(lis_rci_native_session *,const lis_krylov_problem64 *,const lis_krylov_options *,const lis_krylov_ext_options *,lis_krylov_workspace *,const lis_krylov_budget *,double *,double *,double *,lis_krylov_stats *,lis_krylov_ext_stats *,lis_rci_native_session_request *);
int lis_rci_native_session_advance(lis_rci_native_session *,uint64_t,LIS_INT,lis_rci_native_session_request *);
int lis_rci_native_session_cancel(lis_rci_native_session *,uint64_t);
int lis_rci_native_session_pending(const lis_rci_native_session *,lis_rci_native_session_request *);
int lis_rci_native_session_status(const lis_rci_native_session *,LIS_INT *native,uint64_t *a,uint64_t *p,lis_krylov_budget_stats *);
#ifdef __cplusplus
}
#endif
#endif
