/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_krylov_restore.h"
/* Compile original native source in the same TU, sharing all non-suspending
 * numerical kernels and workspace helpers. Not linked into liblis yet. */
#ifndef LIS_RCI_KRYLOV_NATIVE_LIBRARY
#include "../../src/esolver/lis_krylov.c"
#endif
#include "krylov_continuation.inc"

static int map_status(int rc)
{
 switch(rc){
 case LIS_KRYLOV_SUCCESS:return LIS_RCI_DONE;
 case LIS_KRYLOV_BAD_ARGUMENT:return LIS_RCI_BADARG;
 case LIS_KRYLOV_CALLBACK_ERROR:return LIS_RCI_CALLBACK;
 case LIS_KRYLOV_NO_CONVERGENCE:case LIS_KRYLOV_BUDGET_EXHAUSTED:return LIS_RCI_MAXITER;
 default:return LIS_RCI_BREAKDOWN;
 }
}
static int k_return(lis_rci_krylov_restore *K,int parent,int rc)
{
 K->active=parent;K->retval=rc;
 if(parent)return K_STEP_AGAIN;
 K->native_code=rc;K->status=map_status(rc);
 memset(&K->request,0,sizeof K->request);
 return K_STEP_STOP;
}
static int k_abort(lis_rci_krylov_restore *K,int status,int native)
{
 struct frame_lis_krylov_symmetric_core *F=&K->f_lis_krylov_symmetric_core;
 if(K->active){
  if(!F->reuse)release_workspace(&F->s);
  if(F->stats)F->stats->nreturned=F->stats->nconv=0;
 }
 K->active=0;K->status=status;K->native_code=native;
 memset(&K->request,0,sizeof K->request);
 return K_STEP_STOP;
}
static int k_suspend(lis_rci_krylov_restore *K,uint32_t op,int64_t n,const double *in,double *out)
{
 uint64_t *counter,limit;
 if(op==LIS_RCI_OP_A){counter=&K->counts.apply_a;limit=K->budget.max_apply_a;}
 else if(op==LIS_RCI_OP_B){counter=&K->counts.apply_b;limit=K->budget.max_apply_b;}
 else{counter=&K->counts.solve;limit=K->budget.max_solve;}
 if((limit&&*counter>=limit)||(K->budget.max_total&&K->counts.total>=K->budget.max_total)||
    K->counts.total>=INT_MAX/2){
  K->counts.exhausted=1;return k_abort(K,LIS_RCI_MAXITER,LIS_KRYLOV_BUDGET_EXHAUSTED);
 }
 if(K->token==UINT64_MAX)return k_abort(K,LIS_RCI_PROTOCOL,LIS_KRYLOV_CALLBACK_ERROR);
 ++*counter;++K->counts.total;
 /* Preserve the native numerical path and require complete output. All
  * transformed callback sites use nonoverlapping input/output arrays. */
 for(int64_t i=0;i<n;i++)out[i]=NAN;
 K->request=(lis_rci_krylov_restore_request){++K->token,op,n,in,out,K->f_lis_krylov_symmetric_core.p->sigma};
 K->status=(int)op;return K_STEP_YIELD;
}
static int run(lis_rci_krylov_restore *K,lis_rci_krylov_restore_request *q)
{
 int step=K_STEP_AGAIN;
 while(step==K_STEP_AGAIN)step=k_step(K);
 *q=K->request;return K->status;
}
lis_rci_krylov_restore *lis_rci_krylov_restore_create(void)
{
 lis_rci_krylov_restore *K=calloc(1,sizeof(*K));if(K)K->status=LIS_RCI_STATE;return K;
}
void lis_rci_krylov_restore_destroy(lis_rci_krylov_restore *K)
{if(K){if(K->active)(void)k_abort(K,LIS_RCI_CANCELLED,LIS_KRYLOV_CALLBACK_ERROR);free(K);}}
int lis_rci_krylov_restore_start(lis_rci_krylov_restore *K,const lis_krylov_problem64 *p,const lis_krylov_options *o,
 const lis_krylov_ext_options *ext,lis_krylov_workspace *reuse,const lis_krylov_budget *budget,
 double *values,double *vectors,double *residuals,lis_krylov_stats *stats,lis_krylov_ext_stats *ext_stats,lis_rci_krylov_restore_request *q)
{
 if(!K||!q)return LIS_RCI_BADARG;
 if(K->active)return LIS_RCI_STATE;
 if(budget&&(budget->struct_size<sizeof(*budget)||budget->version!=LIS_KRYLOV_BUDGET_VERSION))return LIS_RCI_BADARG;
 uint64_t token=K->token;memset(K,0,sizeof(*K));K->token=token;
 if(budget)K->budget=*budget;
 struct frame_lis_krylov_symmetric_core *F=&K->f_lis_krylov_symmetric_core;
 F->reuse=reuse;F->p=p;F->o=o;F->ext=ext;F->values=values;F->vectors=vectors;
 F->residuals=residuals;F->stats=stats;F->ext_stats=ext_stats;K->active=9;
 return run(K,q);
}
int lis_rci_krylov_restore_advance(lis_rci_krylov_restore *K,uint64_t token,int code,lis_rci_krylov_restore_request *q)
{
 if(!K||!q)return LIS_RCI_BADARG;
 if(!K->active||!K->request.operation)return LIS_RCI_STATE;
 if(token!=K->token)return LIS_RCI_PROTOCOL;
 K->retval=code;memset(&K->request,0,sizeof K->request);
 return run(K,q);
}
int lis_rci_krylov_restore_cancel(lis_rci_krylov_restore *K,uint64_t token)
{
 if(!K)return LIS_RCI_BADARG;
 if(!K->active)return LIS_RCI_STATE;
 if(token!=K->token)return LIS_RCI_PROTOCOL;
 (void)k_abort(K,LIS_RCI_CANCELLED,LIS_KRYLOV_CALLBACK_ERROR);
 return LIS_RCI_CANCELLED;
}
int lis_rci_krylov_restore_status(const lis_rci_krylov_restore *K,int *native,lis_krylov_budget_stats *counts)
{
 if(!K||!native||!counts)return LIS_RCI_BADARG;
 *native=K->native_code;*counts=K->counts;return K->status;
}

int lis_rci_krylov_restore_classic64(const lis_krylov_problem64 *p,const lis_krylov_options *o,
 const lis_krylov_ext_options *ext,lis_krylov_workspace *reuse,const lis_krylov_budget *budget,
 double *values,double *vectors,double *residuals,lis_krylov_stats *stats,lis_krylov_ext_stats *extra,
 lis_krylov_budget_stats *counts)
{
 if(counts)memset(counts,0,sizeof(*counts));
 if(budget&&(budget->struct_size<sizeof(*budget)||budget->version!=LIS_KRYLOV_BUDGET_VERSION))
  return LIS_KRYLOV_BAD_ARGUMENT;
 lis_rci_krylov_restore *K=lis_rci_krylov_restore_create();
 if(!K)return LIS_KRYLOV_NO_MEMORY;
 lis_rci_krylov_restore_request q;
 int st=lis_rci_krylov_restore_start(K,p,o,ext,reuse,budget,values,vectors,residuals,stats,extra,&q);
 while(st>0){
  lis_krylov_matvec64 fn=q.operation==LIS_RCI_OP_A?p->apply_a:
   q.operation==LIS_RCI_OP_B?p->apply_b:p->solve;
  int e=fn(p->ctx,q.n,q.in,q.out);
  st=lis_rci_krylov_restore_advance(K,q.token,e,&q);
 }
 int native=K->native_code;
 if(counts)*counts=K->counts;
 lis_rci_krylov_restore_destroy(K);
 return native;
}
