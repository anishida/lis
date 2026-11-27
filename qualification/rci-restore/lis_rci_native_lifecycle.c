/* Included in the native solver TU in a disposable qualification build. */
#include "lis_rci_native_lifecycle.h"
#include "native_lifecycle_frame.inc"
#define LIFE_YIELD ((LIS_INT)-10000)
#define LIFE_AGAIN ((LIS_INT)-10001)
struct lis_rci_native_lifecycle {
 struct life_frame f;
 lis_rci_native_session *child;
 lis_rci_native_session_request inner,q;
 LIS_VECTOR output,stage,one_in,one_out;
 LIS_PRECON owned_precon;
 LIS_UNSIGNED_INT ag,vg;
 LIS_INT completion,native,threads;
 int phase,status,forced,active;
 uint64_t token;
};
static LIS_INT life_step(lis_rci_native_lifecycle *);
static int life_finite(LIS_SCALAR x)
{
#ifdef _COMPLEX
 return isfinite(creal(x))&&isfinite(cimag(x));
#else
 return isfinite(x);
#endif
}
static LIS_SCALAR life_poison(void)
{
#ifdef _COMPLEX
 return (LIS_SCALAR)(NAN+I*NAN);
#else
 return (LIS_SCALAR)NAN;
#endif
}
static int life_map(LIS_INT rc)
{
 return rc==LIS_SUCCESS?LIS_RCI_DONE:rc==LIS_MAXITER?LIS_RCI_MAXITER:
 rc==LIS_BREAKDOWN?LIS_RCI_BREAKDOWN:rc==LIS_ERR_ILL_ARG?LIS_RCI_BADARG:LIS_RCI_CALLBACK;
}
static LIS_INT life_core_result(lis_rci_native_lifecycle *K,int st)
{
 if(st>0){K->q=K->inner;return LIFE_YIELD;}
 uint64_t a,p;LIS_INT native=0;
 lis_rci_native_session_status(K->child,&native,&a,&p,NULL);
 K->completion=native;
 if(st!=LIS_RCI_DONE&&st!=LIS_RCI_MAXITER&&st!=LIS_RCI_BREAKDOWN){
  K->forced=st;if(!native)K->completion=LIS_BREAKDOWN;
 }
 return LIFE_AGAIN;
}
static LIS_INT life_core(lis_rci_native_lifecycle *K)
{
 K->phase=1;
 int st=lis_rci_native_session_start_linear(K->child,K->f.solver,K->f.nsolver,&K->inner);
 return life_core_result(K,st);
}
static LIS_INT life_residual(lis_rci_native_lifecycle *K)
{
 K->phase=2;
 K->one_in=K->f.xx;K->one_out=K->f.t;
 for(LIS_INT i=0;i<K->one_out->n;i++)K->one_out->value[i]=life_poison();
 K->q=(lis_rci_native_session_request){.operation=LIS_RCI_OP_A,.representation=LIS_NR_VECTOR,.count=1,.A=K->f.A,.solver=K->f.solver,.in=&K->one_in,.out=&K->one_out};
 return LIFE_YIELD;
}
#include "native_lifecycle_step.inc"
static int life_finish(lis_rci_native_lifecycle *K,LIS_INT rc)
{
 LIS_SOLVER s=K->f.solver;
 if(rc!=LIS_SUCCESS){s->retcode=rc;lis_solver_work_destroy(s);}
 K->native=s->retcode;
 if(K->forced){
  s->retcode=K->native=K->completion?K->completion:LIS_BREAKDOWN;s->resid=NAN;
  if(K->f.nsolver==LIS_SOLVER_GCRODR)lis_gcrodr_clear_recycle(s);
  if(K->f.nsolver==LIS_SOLVER_FGCRODR)lis_fgcrodr_clear_recycle(s);
 }
 K->status=K->forced?K->forced:life_map(K->native);
 if(K->status==LIS_RCI_DONE)lis_vector_copy(K->stage,K->output);
 if(s->precon==K->owned_precon)s->precon=NULL;
 lis_precon_destroy(K->owned_precon);K->owned_precon=NULL;
 lis_vector_destroy(K->stage);K->stage=NULL;
 /* Do not leave dangling aliases to the kernel's or wrapper's scratch. */
 s->x=NULL;s->xx=NULL;
 lis_rci_native_session_destroy(K->child);K->child=NULL;
 K->active=0;memset(&K->q,0,sizeof K->q);
 return K->status;
}
static int life_drive(lis_rci_native_lifecycle *K,LIS_INT rc,lis_rci_native_session_request *q)
{
 while(rc==LIFE_AGAIN)rc=life_step(K);
 if(rc!=LIFE_YIELD){int st=life_finish(K,rc);if(q)*q=K->q;return st;}
 /* Token space is practically inexhaustible, but never wrap. */
 if(K->token==UINT64_MAX){
  if(K->phase==1)lis_rci_native_session_cancel(K->child,K->inner.token);
  K->forced=LIS_RCI_PROTOCOL;K->completion=LIS_BREAKDOWN;
  return life_drive(K,life_step(K),q);
 }
 K->q.token=++K->token;K->status=(int)K->q.operation;
 if(q)*q=K->q;return K->status;
}
lis_rci_native_lifecycle *lis_rci_native_lifecycle_create(void)
{lis_rci_native_lifecycle *K=calloc(1,sizeof(*K));if(K)K->status=LIS_RCI_STATE;return K;}
void lis_rci_native_lifecycle_destroy(lis_rci_native_lifecycle *K)
{if(K){if(K->active)lis_rci_native_lifecycle_cancel(K,K->q.token);free(K);}}
int lis_rci_native_lifecycle_start(lis_rci_native_lifecycle *K,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_SOLVER s,lis_rci_native_session_request *q)
{
 if(!K||!A||!b||!x||!s||!q)return LIS_RCI_BADARG;
 if(K->active)return LIS_RCI_STATE;
 LIS_INT method=s->options[LIS_OPTIONS_SOLVER],pre=s->options[LIS_OPTIONS_PRECON];
 int kind=method==LIS_SOLVER_GCRODR?LIS_NR_GCRO:method==LIS_SOLVER_FGCRODR?LIS_NR_FGCRO:LIS_NR_LINEAR;
 if(kind==LIS_NR_LINEAR&&method!=LIS_SOLVER_CG&&method!=LIS_SOLVER_BICGSTAB&&method!=LIS_SOLVER_GMRES&&method!=LIS_SOLVER_FGMRES)return LIS_RCI_BADARG;
 if(A->matrix_type!=LIS_MATRIX_CSR||A->n<=0||!A->ptr||!A->value||!b->value||!x->value||b==x||A->n!=b->n||A->n!=x->n||
 A->is_scaled||b->is_scaled||s->setup||s->num_threads||s->nullspace_dim||s->near_nullspace_dim||s->left_nullspace_dim||
 s->options[LIS_OPTIONS_SCALE]||s->options[LIS_OPTIONS_STORAGE]||s->options[LIS_OPTIONS_ADDS]||
 s->options[LIS_OPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE||(pre!=LIS_PRECON_TYPE_NONE&&pre!=LIS_PRECON_TYPE_JACOBI))return LIS_RCI_BADARG;
 LIS_INT saved,changed;
 if(lis_solver_thread_scope_begin(s,&saved,&changed))return LIS_RCI_BADARG;
 lis_solver_thread_scope_end(saved,changed);
 uint64_t token=K->token;memset(K,0,sizeof(*K));K->token=token;K->threads=saved;K->status=LIS_RCI_STATE;
 K->child=lis_rci_native_session_create(kind);if(!K->child)return LIS_RCI_BADARG;
 LIS_INT rc=lis_vector_duplicate(A,&K->stage);
 if(rc){lis_rci_native_session_destroy(K->child);K->child=NULL;return LIS_RCI_BADARG;}
 lis_vector_copy(x,K->stage);K->output=x;
 s->A=A;
 rc=lis_precon_create(s,&K->owned_precon);
 if(rc){lis_vector_destroy(K->stage);lis_rci_native_session_destroy(K->child);K->stage=NULL;K->child=NULL;s->retcode=rc;return life_map(rc);}
 K->f.A=A;K->f.b=b;K->f.x=K->stage;K->f.solver=s;K->f.precon=K->owned_precon;
 K->ag=A->structure_generation;K->vg=A->value_generation;K->active=1;
 return life_drive(K,life_step(K),q);
}
static int life_interrupt(lis_rci_native_lifecycle *K,int reason,LIS_INT code,lis_rci_native_session_request *q)
{
 K->forced=reason;K->completion=code;
 if(K->phase==1)lis_rci_native_session_cancel(K->child,K->inner.token);
 return life_drive(K,life_step(K),q);
}
int lis_rci_native_lifecycle_advance(lis_rci_native_lifecycle *K,uint64_t token,LIS_INT code,lis_rci_native_session_request *q)
{
 if(!K||!q)return LIS_RCI_BADARG;if(!K->active)return LIS_RCI_STATE;if(token!=K->token)return LIS_RCI_PROTOCOL;
 if(K->f.A->structure_generation!=K->ag||K->f.A->value_generation!=K->vg)return life_interrupt(K,LIS_NR_STALE,LIS_BREAKDOWN,q);
#ifdef _OPENMP
 if(omp_get_max_threads()!=K->threads)return life_interrupt(K,LIS_RCI_STATE,LIS_BREAKDOWN,q);
#endif
 if(code)return life_interrupt(K,LIS_RCI_CALLBACK,code,q);
 if(K->phase==1){
  int st=lis_rci_native_session_advance(K->child,K->inner.token,code,&K->inner);
  return life_drive(K,life_core_result(K,st),q);
 }
 if(code){K->forced=LIS_RCI_CALLBACK;K->completion=code;}
 else{
  K->completion=0;
  for(LIS_INT i=0;i<K->one_out->n;i++)if(!life_finite(K->one_out->value[i])){K->forced=LIS_RCI_NONFINITE;K->completion=LIS_BREAKDOWN;break;}
 }
 return life_drive(K,life_step(K),q);
}
int lis_rci_native_lifecycle_cancel(lis_rci_native_lifecycle *K,uint64_t token)
{if(!K)return LIS_RCI_BADARG;if(!K->active)return LIS_RCI_STATE;if(token!=K->token)return LIS_RCI_PROTOCOL;return life_interrupt(K,LIS_RCI_CANCELLED,LIS_BREAKDOWN,NULL);}
int lis_rci_native_lifecycle_status(const lis_rci_native_lifecycle *K,LIS_INT *native)
{if(!K||!native)return LIS_RCI_BADARG;*native=K->native;return K->status;}
