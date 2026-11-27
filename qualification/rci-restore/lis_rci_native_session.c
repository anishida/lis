/* Qualification-only dispatcher; no numerical recurrence here. */
#include "lis_rci_native_session.h"
#include "lis_rci_native_linear_restore.h"
#include "lis_rci_native_gcro_restore.h"
#include "lis_rci_native_fgcro_restore.h"
#include "lis_rci_native_feast.h"
#include "lis_rci_krylov_restore.h"
#include <limits.h>
struct lis_rci_native_session {
 int kind,status;
 union {lis_rci_native_linear_restore *linear;lis_rci_native_gcro_restore *gcro;
  lis_rci_native_fgcro_restore *fgcro;lis_rci_native_feast *feast;lis_rci_krylov_restore *krylov;} core;
 lis_rci_native_linear_request l;
 lis_rci_native_feast_request f;
 lis_rci_krylov_restore_request k;
 lis_rci_native_session_request q;
};
lis_rci_native_session *lis_rci_native_session_create(int kind)
{
 lis_rci_native_session *s=calloc(1,sizeof(*s));if(!s)return NULL;
 s->kind=kind;s->status=LIS_RCI_STATE;
 switch(kind){
 case LIS_NR_LINEAR:if(s->core.linear=lis_rci_native_linear_restore_create())return s;break;
 case LIS_NR_GCRO:if(s->core.gcro=lis_rci_native_gcro_restore_create())return s;break;
 case LIS_NR_FGCRO:if(s->core.fgcro=lis_rci_native_fgcro_restore_create())return s;break;
 case LIS_NR_FEAST:if(s->core.feast=lis_rci_native_feast_create())return s;break;
 case LIS_NR_KRYLOV:if(s->core.krylov=lis_rci_krylov_restore_create())return s;break;
 }
 free(s);return NULL;
}
void lis_rci_native_session_destroy(lis_rci_native_session *s)
{
 if(!s)return;
 switch(s->kind){
 case LIS_NR_LINEAR:lis_rci_native_linear_restore_destroy(s->core.linear);break;
 case LIS_NR_GCRO:lis_rci_native_gcro_restore_destroy(s->core.gcro);break;
 case LIS_NR_FGCRO:lis_rci_native_fgcro_restore_destroy(s->core.fgcro);break;
 case LIS_NR_FEAST:lis_rci_native_feast_destroy(s->core.feast);break;
 case LIS_NR_KRYLOV:lis_rci_krylov_restore_destroy(s->core.krylov);break;
 }
 free(s);
}
static int nr_feast_status(lis_rci_native_session *s,int st)
{
 switch(st){
 case NF_DONE:return LIS_RCI_DONE;
 case NF_A:return LIS_RCI_OP_A;
 case NF_B:return LIS_RCI_OP_B;
 case NF_SHIFT:return LIS_RCI_OP_SHIFTED_SOLVE;
 case NF_PROTOCOL:return LIS_RCI_PROTOCOL;
 case NF_CANCELLED:return LIS_RCI_CANCELLED;
 case NF_NONFINITE:return LIS_RCI_NONFINITE;
 case NF_STATE:return LIS_RCI_STATE;
 case NF_BADARG:return LIS_RCI_BADARG;
 case NF_STALE:return LIS_NR_STALE;
 default:{LIS_INT n=lis_rci_native_feast_native(s->core.feast);
 return n==LIS_MAXITER?LIS_RCI_MAXITER:n==LIS_BREAKDOWN?LIS_RCI_BREAKDOWN:LIS_RCI_CALLBACK;}
 }
}
static int nr_publish(lis_rci_native_session *s,int st,lis_rci_native_session_request *out)
{
 if(s->kind==LIS_NR_FEAST)st=nr_feast_status(s,st);
 memset(&s->q,0,sizeof(s->q));
 s->status=st;
 if(st>0){
  s->q.operation=(uint32_t)st;
  if(s->kind==LIS_NR_KRYLOV){
   s->q.token=s->k.token;s->q.representation=LIS_NR_RAW64;s->q.count=1;
   s->q.n64=s->k.n;s->q.raw_in=s->k.in;s->q.raw_out=s->k.out;s->q.sigma=s->k.sigma;
   if(st==LIS_RCI_OP_SHIFTED_SOLVE)s->q.shift_sign=-1;
  }else if(s->kind==LIS_NR_FEAST){
   s->q.token=s->f.token;s->q.A=s->f.A;s->q.B=s->f.B;s->q.solver=s->f.solver;
   if(st==LIS_RCI_OP_SHIFTED_SOLVE){
    s->q.representation=LIS_NR_SPLIT;s->q.count=1;s->q.shift_sign=1;
    s->q.zr=s->f.zr;s->q.zi=s->f.zi;s->q.in_real=s->f.in_real;s->q.in_imag=s->f.in_imag;
    s->q.out_real=s->f.out_real;s->q.out_imag=s->f.out_imag;
   }else{
    s->q.representation=LIS_NR_VECTOR;s->q.count=s->f.count;s->q.block=s->f.block;
    s->q.in=s->f.in;s->q.out=s->f.out;
   }
  }else{
   s->q.token=s->l.token;s->q.representation=LIS_NR_VECTOR;s->q.count=1;
   s->q.A=s->l.A;s->q.solver=s->l.solver;s->q.in=&s->l.in;s->q.out=&s->l.out;
  }
 }
 if(out)*out=s->q;return st;
}
int lis_rci_native_session_start_linear(lis_rci_native_session *s,LIS_SOLVER solver,LIS_INT method,lis_rci_native_session_request *q)
{
 if(!s||!solver||!q)return LIS_RCI_BADARG;
 if(s->status>0)return LIS_RCI_STATE;
 int st;
 switch(s->kind){
 case LIS_NR_LINEAR:st=lis_rci_native_linear_restore_start_method(s->core.linear,solver,method,1,&s->l);break;
 case LIS_NR_GCRO:if(method!=LIS_SOLVER_GCRODR)return LIS_RCI_BADARG;
  st=lis_rci_native_gcro_restore_start(s->core.gcro,solver,1,&s->l);break;
 case LIS_NR_FGCRO:if(method!=LIS_SOLVER_FGCRODR)return LIS_RCI_BADARG;
  st=lis_rci_native_fgcro_restore_start(s->core.fgcro,solver,1,&s->l);break;
 default:return LIS_RCI_BADARG;
 }
 return nr_publish(s,st,q);
}
int lis_rci_native_session_start_feast(lis_rci_native_session *s,LIS_ESOLVER e,lis_rci_native_session_request *q)
{
 if(!s||!e||!q||s->kind!=LIS_NR_FEAST)return LIS_RCI_BADARG;
 if(s->status>0)return LIS_RCI_STATE;
 return nr_publish(s,lis_rci_native_feast_start(s->core.feast,e,1,&s->f),q);
}
int lis_rci_native_session_start_krylov(lis_rci_native_session *s,const lis_krylov_problem64 *p,const lis_krylov_options *o,const lis_krylov_ext_options *ext,lis_krylov_workspace *w,const lis_krylov_budget *b,double *v,double *x,double *r,lis_krylov_stats *stats,lis_krylov_ext_stats *es,lis_rci_native_session_request *q)
{
 if(!s||!q||s->kind!=LIS_NR_KRYLOV)return LIS_RCI_BADARG;
 if(s->status>0)return LIS_RCI_STATE;
 return nr_publish(s,lis_rci_krylov_restore_start(s->core.krylov,p,o,ext,w,b,v,x,r,stats,es,&s->k),q);
}
int lis_rci_native_session_advance(lis_rci_native_session *s,uint64_t token,LIS_INT code,lis_rci_native_session_request *q)
{
 if(!s||!q)return LIS_RCI_BADARG;
 if(s->status<=0)return LIS_RCI_STATE;
 if(token!=s->q.token)return LIS_RCI_PROTOCOL;
 int st;
 switch(s->kind){
 case LIS_NR_LINEAR:st=lis_rci_native_linear_restore_advance(s->core.linear,token,code,&s->l);break;
 case LIS_NR_GCRO:st=lis_rci_native_gcro_restore_advance(s->core.gcro,token,code,&s->l);break;
 case LIS_NR_FGCRO:st=lis_rci_native_fgcro_restore_advance(s->core.fgcro,token,code,&s->l);break;
 case LIS_NR_FEAST:st=lis_rci_native_feast_advance(s->core.feast,token,code,&s->f);break;
 case LIS_NR_KRYLOV:
  if(code<INT_MIN||code>INT_MAX)return LIS_RCI_BADARG;
  st=lis_rci_krylov_restore_advance(s->core.krylov,token,(int)code,&s->k);break;
 default:return LIS_RCI_BADARG;
 }
 return nr_publish(s,st,q);
}
int lis_rci_native_session_cancel(lis_rci_native_session *s,uint64_t token)
{
 if(!s)return LIS_RCI_BADARG;
 if(s->status<=0)return LIS_RCI_STATE;
 if(token!=s->q.token)return LIS_RCI_PROTOCOL;
 int st;
 switch(s->kind){
 case LIS_NR_LINEAR:st=lis_rci_native_linear_restore_cancel(s->core.linear,token);break;
 case LIS_NR_GCRO:st=lis_rci_native_gcro_restore_cancel(s->core.gcro,token);break;
 case LIS_NR_FGCRO:st=lis_rci_native_fgcro_restore_cancel(s->core.fgcro,token);break;
 case LIS_NR_FEAST:st=lis_rci_native_feast_cancel(s->core.feast,token);break;
 case LIS_NR_KRYLOV:st=lis_rci_krylov_restore_cancel(s->core.krylov,token);break;
 default:return LIS_RCI_BADARG;
 }
 return nr_publish(s,st,NULL);
}
int lis_rci_native_session_pending(const lis_rci_native_session *s,lis_rci_native_session_request *q)
{
 if(!s||!q)return LIS_RCI_BADARG;*q=s->q;
 return s->status>0?s->status:LIS_RCI_STATE;
}
int lis_rci_native_session_status(const lis_rci_native_session *s,LIS_INT *native,uint64_t *a,uint64_t *p,lis_krylov_budget_stats *budget)
{
 if(!s||!native||!a||!p)return LIS_RCI_BADARG;
 *a=*p=0;*native=0;if(budget)memset(budget,0,sizeof(*budget));
 switch(s->kind){
 case LIS_NR_LINEAR:lis_rci_native_linear_restore_status(s->core.linear,native,a,p);break;
 case LIS_NR_GCRO:lis_rci_native_gcro_restore_status(s->core.gcro,native,a,p);break;
 case LIS_NR_FGCRO:lis_rci_native_fgcro_restore_status(s->core.fgcro,native,a,p);break;
 case LIS_NR_FEAST:*native=lis_rci_native_feast_native(s->core.feast);break;
 case LIS_NR_KRYLOV:{int n=0;lis_krylov_budget_stats local;lis_rci_krylov_restore_status(s->core.krylov,&n,budget?budget:&local);*native=n;break;}
 }
 return s->status;
}
