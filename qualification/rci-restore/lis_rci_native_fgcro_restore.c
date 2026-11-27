/* AI-assisted continuation reconstruction from the actual native LIS sources.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_native_fgcro_restore.h"
#if defined(USE_MPI)||defined(USE_QUAD_PRECISION)
#error "native FGCRO continuation candidate requires non-MPI without quad solver mode"
#endif
#include "fgcro_continuation.inc"
static LIS_SCALAR native_poison_value(void)
{
#ifdef _COMPLEX
 return (LIS_SCALAR)(NAN + I*NAN);
#else
 return (LIS_SCALAR)NAN;
#endif
}
static int native_scalar_finite(LIS_SCALAR x)
{
#ifdef _COMPLEX
 return isfinite(creal(x))&&isfinite(cimag(x));
#else
 return isfinite(x);
#endif
}
/* forced_status uses the otherwise unused generated method field. */
static int linear_return(lis_rci_native_fgcro_restore *K,int parent,LIS_INT rc)
{
 K->active=parent;K->retval=rc;
 if(parent)return L_STEP_AGAIN;
 if(K->method){ lis_fgcrodr_clear_recycle(K->solver); rc=LIS_BREAKDOWN; K->solver->retcode=rc; }
 K->native_code=rc;
 K->status=K->method?K->method:rc==LIS_SUCCESS?LIS_RCI_DONE:rc==LIS_MAXITER?LIS_RCI_MAXITER:rc==LIS_BREAKDOWN?LIS_RCI_BREAKDOWN:LIS_RCI_CALLBACK;
 memset(&K->request,0,sizeof K->request);
 return L_STEP_STOP;
}
static int linear_suspend(lis_rci_native_fgcro_restore *K,uint32_t op,LIS_MATRIX A,LIS_VECTOR in,LIS_VECTOR out)
{
 if(K->token==UINT64_MAX&&!K->method)K->method=LIS_RCI_PROTOCOL;
 if(K->strict&&!K->method){
  if(!in||!out||in==out||in->n!=out->n)K->method=LIS_RCI_BADARG;
  else for(LIS_INT i=0;i<in->n;i++)if(!native_scalar_finite(in->value[i])){K->method=LIS_RCI_NONFINITE;break;}
 }
 /* Resume the original failing call and unwind all original cleanup labels. */
 if(K->method){K->retval=LIS_BREAKDOWN;return L_STEP_AGAIN;}
 if(K->strict)for(LIS_INT i=0;i<out->n;i++)out->value[i]=native_poison_value();
 if(op==LIS_RCI_OP_A)K->matvecs++;else K->precons++;
 K->request=(lis_rci_native_linear_request){++K->token,op,A,K->solver,in,out};
 K->status=(int)op;return L_STEP_YIELD;
}
static int drive(lis_rci_native_fgcro_restore *K,lis_rci_native_linear_request *q)
{
 int step=L_STEP_AGAIN;
 while(step==L_STEP_AGAIN)step=linear_step(K);
 if(q)*q=K->request;
 return K->status;
}
lis_rci_native_fgcro_restore *lis_rci_native_fgcro_restore_create(void)
{
 lis_rci_native_fgcro_restore *K=calloc(1,sizeof(*K));
 if(K)K->status=LIS_RCI_STATE;
 return K;
}
void lis_rci_native_fgcro_restore_destroy(lis_rci_native_fgcro_restore *K)
{
 if(K){
  if(K->active)(void)lis_rci_native_fgcro_restore_cancel(K,K->token);
  free(K);
 }
}
int lis_rci_native_fgcro_restore_start(lis_rci_native_fgcro_restore *K,LIS_SOLVER s,int strict,lis_rci_native_linear_request *q)
{
 if(!K||!q||!s||!s->A||!s->b||!s->x||!s->work||s->options[LIS_OPTIONS_SOLVER]!=LIS_SOLVER_FGCRODR||(s->fgcrodr_state&&s->fgcrodr_state->target_dim==0)||(strict!=0&&strict!=1))return LIS_RCI_BADARG;
 if(K->active)return LIS_RCI_STATE;
 uint64_t token=K->token;memset(K,0,sizeof(*K));K->token=token;
 K->solver=s;K->strict=strict;K->fgcro.solver=s;K->active=10;
 return drive(K,q);
}
int lis_rci_native_fgcro_restore_advance(lis_rci_native_fgcro_restore *K,uint64_t token,LIS_INT code,lis_rci_native_linear_request *q)
{
 if(!K||!q)return LIS_RCI_BADARG;
 if(!K->active||!K->request.operation)return LIS_RCI_STATE;
 if(token!=K->token)return LIS_RCI_PROTOCOL;
 if(K->strict&&!code){
  for(LIS_INT i=0;i<K->request.out->n;i++)if(!native_scalar_finite(K->request.out->value[i])){
   K->method=LIS_RCI_NONFINITE;code=LIS_BREAKDOWN;break;
  }
 }
 K->retval=code;memset(&K->request,0,sizeof K->request);return drive(K,q);
}
int lis_rci_native_fgcro_restore_cancel(lis_rci_native_fgcro_restore *K,uint64_t token)
{
 if(!K)return LIS_RCI_BADARG;
 if(!K->active||!K->request.operation)return LIS_RCI_STATE;
 if(token!=K->token)return LIS_RCI_PROTOCOL;
 K->method=LIS_RCI_CANCELLED;K->retval=LIS_BREAKDOWN;
 memset(&K->request,0,sizeof K->request);
 return drive(K,NULL);
}
int lis_rci_native_fgcro_restore_status(const lis_rci_native_fgcro_restore *K,LIS_INT *native,uint64_t *a,uint64_t *p)
{
 if(!K||!native||!a||!p)return LIS_RCI_BADARG;
 *native=K->native_code;*a=K->matvecs;*p=K->precons;return K->status;
}
static _Thread_local int fgcro_reference;
static _Thread_local lis_rci_native_fgcro_test_driver fgcro_driver;
void lis_rci_native_fgcro_test_select(int reference,lis_rci_native_fgcro_test_driver driver)
{fgcro_reference=reference;fgcro_driver=driver;}
LIS_INT lis_fgcrodr(LIS_SOLVER s)
{
 if(fgcro_reference || (s && s->fgcrodr_state && s->fgcrodr_state->target_dim==0))return lis_fgcrodr_native_reference(s);
 if(fgcro_driver)return fgcro_driver(s);
 lis_rci_native_fgcro_restore *K=lis_rci_native_fgcro_restore_create();
 if(!K)return LIS_ERR_OUT_OF_MEMORY;
 lis_rci_native_linear_request q;
 int st=lis_rci_native_fgcro_restore_start(K,s,0,&q);
 while(st>0){
  LIS_INT e=q.operation==LIS_RCI_OP_A?lis_matvec(q.A,q.in,q.out):lis_psolve(q.solver,q.in,q.out);
  st=lis_rci_native_fgcro_restore_advance(K,q.token,e,&q);
 }
 LIS_INT native=K->native_code;
 if(st==LIS_RCI_BADARG&&K->status==LIS_RCI_STATE)native=LIS_ERR_ILL_ARG;
 lis_rci_native_fgcro_restore_destroy(K);return native;
}

/* Keep private state private: differential harness receives only a copy. */
LIS_INT lis_rci_native_fgcro_test_basis(LIS_SOLVER s,LIS_SCALAR *out,size_t capacity)
{
 if(!s||!s->fgcrodr_state||!out)return LIS_ERR_ILL_ARG;
 LIS_FGCRODR_STATE r=s->fgcrodr_state;
 if(!r->valid||r->dim<=0||!s->A||s->A->n<=0)return LIS_ERR_ILL_ARG;
 size_t n=(size_t)s->A->n,k=(size_t)r->dim;
 if(k>capacity/n/3)return LIS_ERR_ILL_ARG;
 for(size_t i=0;i<k;i++){
  memcpy(out+i*n,r->w[i]->value,n*sizeof(*out));
  memcpy(out+(k+i)*n,r->z[i]->value,n*sizeof(*out));
  memcpy(out+(2*k+i)*n,r->c[i]->value,n*sizeof(*out));
 }
 return LIS_SUCCESS;
}
