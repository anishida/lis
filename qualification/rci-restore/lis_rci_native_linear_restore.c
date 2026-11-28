/* AI-assisted continuation reconstruction from the actual native LIS sources.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_native_linear_restore.h"
#include "lislib.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#if defined(USE_MPI) || defined(USE_QUAD_PRECISION)
#error "native linear continuation qualification is non-MPI without double-double solver mode"
#endif
#include "linear_continuation.inc"

static int scalar_finite(LIS_SCALAR v)
{
#ifdef _COMPLEX
 return isfinite(creal(v))&&isfinite(cimag(v));
#else
 return isfinite(v);
#endif
}

static int map_native(LIS_INT rc)
{
 if(rc==LIS_SUCCESS)return LIS_RCI_DONE;
 if(rc==LIS_MAXITER)return LIS_RCI_MAXITER;
 if(rc==LIS_BREAKDOWN)return LIS_RCI_BREAKDOWN;
 if(rc==LIS_ERR_ILL_ARG)return LIS_RCI_BADARG;
 if(rc==LIS_ERR_OUT_OF_MEMORY)return LIS_RCI_BADARG;
 return LIS_RCI_CALLBACK;
}
static int linear_return(lis_rci_native_linear_restore *K,int parent,LIS_INT rc)
{
 K->active=parent;K->retval=rc;
 if(parent)return L_STEP_AGAIN;
 K->native_code=rc;K->status=map_native(rc);
 memset(&K->request,0,sizeof K->request);return L_STEP_STOP;
}
static int abort_linear(lis_rci_native_linear_restore *K,int status)
{
 if(K->gm.h){lis_free(K->gm.h);K->gm.h=NULL;}
 if(K->fgm.h){lis_free(K->fgm.h);K->fgm.h=NULL;}
 K->active=0;K->status=status;K->native_code=LIS_BREAKDOWN;
 if(K->solver){K->solver->retcode=LIS_BREAKDOWN;K->solver->resid=NAN;}
 memset(&K->request,0,sizeof K->request);return L_STEP_STOP;
}
static int linear_suspend(lis_rci_native_linear_restore *K,uint32_t op,LIS_MATRIX A,LIS_VECTOR in,LIS_VECTOR out)
{
 if(K->token==UINT64_MAX)return abort_linear(K,LIS_RCI_PROTOCOL);
 if(K->strict){
  if(!in||!out||in==out||in->n!=out->n)return abort_linear(K,LIS_RCI_BADARG);
  for(LIS_INT i=0;i<in->n;i++)if(!scalar_finite(in->value[i]))return abort_linear(K,LIS_RCI_NONFINITE);
  for(LIS_INT i=0;i<out->n;i++)out->value[i]=
#ifdef _COMPLEX
   (LIS_SCALAR)(NAN+I*NAN);
#else
   (LIS_SCALAR)NAN;
#endif
 }
 if(op==LIS_RCI_OP_A)K->matvecs++;else K->precons++;
 K->request=(lis_rci_native_linear_request){++K->token,op,A,K->solver,in,out};
 K->status=(int)op;return L_STEP_YIELD;
}
static int drive(lis_rci_native_linear_restore *K,lis_rci_native_linear_request *q)
{
 int step=L_STEP_AGAIN;while(step==L_STEP_AGAIN)step=linear_step(K);
 *q=K->request;return K->status;
}
lis_rci_native_linear_restore *lis_rci_native_linear_restore_create(void)
{
 lis_rci_native_linear_restore *K=calloc(1,sizeof(*K));if(K)K->status=LIS_RCI_STATE;return K;
}
void lis_rci_native_linear_restore_destroy(lis_rci_native_linear_restore *K)
{if(K){if(K->active)(void)abort_linear(K,LIS_RCI_CANCELLED);free(K);}}
int lis_rci_native_linear_restore_start_method(lis_rci_native_linear_restore *K,LIS_SOLVER solver,LIS_INT method,int strict,lis_rci_native_linear_request *q)
{
 if(!K||!q||!solver||!solver->A||!solver->b||!solver->x||!solver->work||
    (strict!=0&&strict!=1))return LIS_RCI_BADARG;
 if(K->active)return LIS_RCI_STATE;
 if(method!=LIS_SOLVER_CG&&method!=LIS_SOLVER_BICGSTAB&&method!=LIS_SOLVER_GMRES&&method!=LIS_SOLVER_FGMRES)return LIS_RCI_BADARG;
 uint64_t token=K->token;memset(K,0,sizeof(*K));K->token=token;K->solver=solver;K->strict=strict;
 K->method=method;
 if(method==LIS_SOLVER_CG){K->active=2;K->cg.solver=solver;}
 else if(method==LIS_SOLVER_BICGSTAB){K->active=3;K->bi.solver=solver;}
 else if(method==LIS_SOLVER_GMRES){K->active=4;K->gm.solver=solver;}
 else{K->active=5;K->fgm.solver=solver;}
 return drive(K,q);
}
int lis_rci_native_linear_restore_start(lis_rci_native_linear_restore *K,LIS_SOLVER solver,int strict,lis_rci_native_linear_request *q)
{
 if(!solver)return LIS_RCI_BADARG;
 return lis_rci_native_linear_restore_start_method(K,solver,solver->options[LIS_OPTIONS_SOLVER],strict,q);
}
int lis_rci_native_linear_restore_advance(lis_rci_native_linear_restore *K,uint64_t token,LIS_INT code,lis_rci_native_linear_request *q)
{
 if(!K||!q)return LIS_RCI_BADARG;
 if(!K->active||!K->request.operation)return LIS_RCI_STATE;
 if(token!=K->token)return LIS_RCI_PROTOCOL;
 if(K->strict&&!code){
  LIS_VECTOR out=K->request.out;
  for(LIS_INT i=0;i<out->n;i++)if(!scalar_finite(out->value[i])){
   (void)abort_linear(K,LIS_RCI_NONFINITE);*q=K->request;return K->status;
  }
 }
 K->retval=code;memset(&K->request,0,sizeof K->request);return drive(K,q);
}
int lis_rci_native_linear_restore_cancel(lis_rci_native_linear_restore *K,uint64_t token)
{
 if(!K)return LIS_RCI_BADARG;
 if(!K->active)return LIS_RCI_STATE;
 if(token!=K->token)return LIS_RCI_PROTOCOL;
 (void)abort_linear(K,LIS_RCI_CANCELLED);return LIS_RCI_CANCELLED;
}
int lis_rci_native_linear_restore_status(const lis_rci_native_linear_restore *K,LIS_INT *native,uint64_t *a,uint64_t *p)
{
 if(!K||!native||!a||!p)return LIS_RCI_BADARG;
 *native=K->native_code;*a=K->matvecs;*p=K->precons;return K->status;
}
static LIS_INT classic_method(LIS_SOLVER solver,LIS_INT method)
{
 lis_rci_native_linear_restore *K=lis_rci_native_linear_restore_create();if(!K)return LIS_ERR_OUT_OF_MEMORY;
 lis_rci_native_linear_request q;
 int st=lis_rci_native_linear_restore_start_method(K,solver,method,0,&q);
 while(st>0){
  LIS_INT e=q.operation==LIS_RCI_OP_A?lis_matvec(q.A,q.in,q.out):lis_psolve(q.solver,q.in,q.out);
  st=lis_rci_native_linear_restore_advance(K,q.token,e,&q);
 }
 LIS_INT native=K->native_code;
 if(st==LIS_RCI_BADARG&&K->status==LIS_RCI_STATE)native=LIS_ERR_ILL_ARG;
 lis_rci_native_linear_restore_destroy(K);return native;
}
LIS_INT lis_rci_native_linear_restore_classic(LIS_SOLVER solver)
{
 return solver?classic_method(solver,solver->options[LIS_OPTIONS_SOLVER]):LIS_ERR_ILL_ARG;
}
/* The hook is limited to this disposable qualification build. */
static _Thread_local int reference_mode;
static _Thread_local lis_rci_native_linear_test_driver external_driver;
void lis_rci_native_linear_test_select(int reference,lis_rci_native_linear_test_driver driver)
{reference_mode=reference;external_driver=driver;}
extern LIS_INT lis_cg_native_reference(LIS_SOLVER);
extern LIS_INT lis_bicgstab_native_reference(LIS_SOLVER);
extern LIS_INT lis_gmres_native_reference(LIS_SOLVER);
extern LIS_INT lis_fgmres_native_reference(LIS_SOLVER);
static LIS_INT dispatch_linear(LIS_SOLVER solver,LIS_INT method,LIS_INT (*reference)(LIS_SOLVER))
{
 if(reference_mode)return reference(solver);
 if(external_driver)return external_driver(solver,method);
 return classic_method(solver,method);
}
LIS_INT lis_cg(LIS_SOLVER s){return dispatch_linear(s,LIS_SOLVER_CG,lis_cg_native_reference);}
LIS_INT lis_bicgstab(LIS_SOLVER s){return dispatch_linear(s,LIS_SOLVER_BICGSTAB,lis_bicgstab_native_reference);}
LIS_INT lis_gmres(LIS_SOLVER s){return dispatch_linear(s,LIS_SOLVER_GMRES,lis_gmres_native_reference);}
LIS_INT lis_fgmres(LIS_SOLVER s){return dispatch_linear(s,LIS_SOLVER_FGMRES,lis_fgmres_native_reference);}
