#include <limits.h>
/* Public facade: numerical recurrences remain in the qualified native cores. */
#include "lis_rci_native.h"
#include "lis_rci_native_lifecycle.h"
#include "lis_rci_native_feast_lifecycle.h"
struct lis_rci_native_context {
 int kind,status;
 union {lis_rci_native_lifecycle *linear;lis_rci_native_feast_lifecycle *feast;lis_rci_native_session *krylov;} core;
 lis_rci_native_session_request pending;
};
lis_rci_native_context *lis_rci_native_create(int kind,uint32_t version,size_t bytes,const LIS_ABI_INFO *client,size_t client_size)
{
 LIS_ABI_INFO linked;
 if(version!=LIS_RCI_NATIVE_API_VERSION||bytes!=sizeof(lis_rci_native_session_request)||
    lis_abi_get_info(&linked,sizeof linked)||
    !lis_abi_is_layout_compatible(&linked,sizeof linked,client,client_size,LIS_ABI_LAYOUT_FLAGS))return NULL;
 lis_rci_native_context *s=calloc(1,sizeof(*s));if(!s)return NULL;
 s->kind=kind;s->status=LIS_RCI_STATE;
 switch(kind){
 case LIS_NR_LINEAR:if((s->core.linear=lis_rci_native_lifecycle_create()))return s;break;
 case LIS_NR_FEAST:if((s->core.feast=lis_rci_native_feast_lifecycle_create()))return s;break;
 case LIS_NR_KRYLOV:if((s->core.krylov=lis_rci_native_session_create(kind)))return s;break;
 }
 free(s);return NULL;
}
void lis_rci_native_destroy(lis_rci_native_context *s)
{
 if(!s)return;
 switch(s->kind){
 case LIS_NR_LINEAR:lis_rci_native_lifecycle_destroy(s->core.linear);break;
 case LIS_NR_FEAST:lis_rci_native_feast_lifecycle_destroy(s->core.feast);break;
 case LIS_NR_KRYLOV:lis_rci_native_session_destroy(s->core.krylov);break;
 }
 free(s);
}
static int public_publish(lis_rci_native_context *s,int st,lis_rci_native_session_request *q)
{
 s->status=st;
 if(st<=0)memset(&s->pending,0,sizeof s->pending);
 if(q)*q=s->pending;
 return st;
}
int lis_rci_native_start_linear(lis_rci_native_context *s,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_SOLVER solver,lis_rci_native_session_request *q)
{
 if(!s||!q||s->kind!=LIS_NR_LINEAR)return LIS_RCI_BADARG;
 if(s->status>0)return LIS_RCI_STATE;
 return public_publish(s,lis_rci_native_lifecycle_start(s->core.linear,A,b,x,solver,&s->pending),q);
}
int lis_rci_native_start_feast(lis_rci_native_context *s,LIS_MATRIX A,LIS_MATRIX B,LIS_VECTOR x,LIS_SCALAR *v,LIS_ESOLVER e,lis_rci_native_session_request *q)
{
 if(!s||!q||s->kind!=LIS_NR_FEAST)return LIS_RCI_BADARG;
 if(s->status>0)return LIS_RCI_STATE;
 return public_publish(s,lis_rci_native_feast_lifecycle_start(s->core.feast,A,B,x,v,e,&s->pending),q);
}
int lis_rci_native_start_krylov(lis_rci_native_context *s,const lis_krylov_problem64 *p,const lis_krylov_options *o,const lis_krylov_ext_options *ext,lis_krylov_workspace *w,const lis_krylov_budget *b,double *v,double *x,double *r,lis_krylov_stats *stats,lis_krylov_ext_stats *es,lis_rci_native_session_request *q)
{
 if(!s||!q||s->kind!=LIS_NR_KRYLOV)return LIS_RCI_BADARG;
 if(s->status>0)return LIS_RCI_STATE;
 return public_publish(s,lis_rci_native_session_start_krylov(s->core.krylov,p,o,ext,w,b,v,x,r,stats,es,&s->pending),q);
}
int lis_rci_native_advance(lis_rci_native_context *s,uint64_t token,LIS_INT code,lis_rci_native_session_request *q)
{
 if(!s||!q)return LIS_RCI_BADARG;
 if(s->status<=0)return LIS_RCI_STATE;
 if(token!=s->pending.token)return LIS_RCI_PROTOCOL;
 int st;
 switch(s->kind){
 case LIS_NR_LINEAR:st=lis_rci_native_lifecycle_advance(s->core.linear,token,code,&s->pending);break;
 case LIS_NR_FEAST:st=lis_rci_native_feast_lifecycle_advance(s->core.feast,token,code,&s->pending);break;
 case LIS_NR_KRYLOV:
  if(code<INT_MIN||code>INT_MAX)return LIS_RCI_BADARG;
  if(code){lis_rci_native_session_cancel(s->core.krylov,token);st=LIS_RCI_CALLBACK;break;}
  for(int64_t i=0;i<s->pending.n64;i++)if(!isfinite(s->pending.raw_out[i])){
   lis_rci_native_session_cancel(s->core.krylov,token);
   return public_publish(s,LIS_RCI_NONFINITE,q);
  }
  st=lis_rci_native_session_advance(s->core.krylov,token,code,&s->pending);break;
 default:return LIS_RCI_BADARG;
 }
 return public_publish(s,st,q);
}
int lis_rci_native_cancel(lis_rci_native_context *s,uint64_t token)
{
 if(!s)return LIS_RCI_BADARG;
 if(s->status<=0)return LIS_RCI_STATE;
 if(token!=s->pending.token)return LIS_RCI_PROTOCOL;
 int st;
 switch(s->kind){
 case LIS_NR_LINEAR:st=lis_rci_native_lifecycle_cancel(s->core.linear,token);break;
 case LIS_NR_FEAST:st=lis_rci_native_feast_lifecycle_cancel(s->core.feast,token);break;
 case LIS_NR_KRYLOV:st=lis_rci_native_session_cancel(s->core.krylov,token);break;
 default:return LIS_RCI_BADARG;
 }
 return public_publish(s,st,NULL);
}
int lis_rci_native_pending(const lis_rci_native_context *s,lis_rci_native_session_request *q)
{if(!s||!q)return LIS_RCI_BADARG;*q=s->pending;return s->status>0?s->status:LIS_RCI_STATE;}
int lis_rci_native_status(const lis_rci_native_context *s,LIS_INT *native)
{
 if(!s||!native)return LIS_RCI_BADARG;*native=0;
 uint64_t a,p;
 switch(s->kind){
 case LIS_NR_LINEAR:lis_rci_native_lifecycle_status(s->core.linear,native);break;
 case LIS_NR_FEAST:lis_rci_native_feast_lifecycle_status(s->core.feast,native);break;
 case LIS_NR_KRYLOV:lis_rci_native_session_status(s->core.krylov,native,&a,&p,NULL);break;
 }
 return s->status;
}
LIS_INT lis_rci_native_execute_default(lis_rci_native_context *s,uint64_t token)
{
 if(!s||s->status<=0||token!=s->pending.token)return LIS_ERR_ILL_ARG;
 lis_rci_native_session_request *q=&s->pending;
 if(q->representation==LIS_NR_SPLIT){
  if(q->operation!=LIS_NR_OP_SHIFTED_SOLVE||q->shift_sign!=1)return LIS_ERR_ILL_ARG;
  lis_vector_set_all(0,q->out_real);lis_vector_set_all(0,q->out_imag);
  return lis_efeast_shifted_solve_split(q->A,q->B,q->zr,q->zi,q->in_real,q->in_imag,q->out_real,q->out_imag,q->solver);
 }
 if(q->representation!=LIS_NR_VECTOR||!q->in||!q->out||q->count<=0)return LIS_ERR_ILL_ARG;
 if(q->operation==LIS_NR_OP_PRECON)return lis_psolve(q->solver,q->in[0],q->out[0]);
 if(q->operation!=LIS_NR_OP_A&&q->operation!=LIS_NR_OP_B)return LIS_ERR_ILL_ARG;
 return q->block?lis_matvec_block(q->A,q->count,q->in,q->out):lis_matvec(q->A,q->in[0],q->out[0]);
}
