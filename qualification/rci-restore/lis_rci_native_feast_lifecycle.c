/* Included in native esolver TU only in qualification checkout. */
#include "lis_rci_native_feast_lifecycle.h"
#include "native_feast_lifecycle_frame.inc"
#define EF_YIELD ((LIS_INT)-10000)
#define EF_AGAIN ((LIS_INT)-10001)
struct lis_rci_native_feast_lifecycle {
 struct ef_frame f;
 lis_rci_native_session *child;
 lis_rci_native_session_request inner,q;
 LIS_VECTOR output,stage;
 LIS_SCALAR *value,stage_value;
 LIS_INT completion,native;
 int phase,status,forced,active;
 uint64_t token;
};
static LIS_INT ef_step(lis_rci_native_feast_lifecycle *);
static int ef_map(LIS_INT n)
{return n==LIS_SUCCESS?LIS_RCI_DONE:n==LIS_MAXITER?LIS_RCI_MAXITER:n==LIS_BREAKDOWN?LIS_RCI_BREAKDOWN:n==LIS_ERR_ILL_ARG?LIS_RCI_BADARG:LIS_RCI_CALLBACK;}
static LIS_INT ef_core_result(lis_rci_native_feast_lifecycle *K,int st)
{
 if(st>0){K->q=K->inner;return EF_YIELD;}
 uint64_t a,p;LIS_INT native=0;
 lis_rci_native_session_status(K->child,&native,&a,&p,NULL);
 K->completion=native;
 if(st!=LIS_RCI_DONE){K->forced=st;if(!native)K->completion=LIS_BREAKDOWN;}
 return EF_AGAIN;
}
static LIS_INT ef_core(lis_rci_native_feast_lifecycle *K)
{
 K->phase=1;
 return ef_core_result(K,lis_rci_native_session_start_feast(K->child,K->f.esolver,&K->inner));
}
#include "native_feast_lifecycle_step.inc"
static int ef_finish(lis_rci_native_feast_lifecycle *K,LIS_INT rc)
{
 LIS_ESOLVER e=K->f.esolver;
 K->native=rc?rc:e->retcode;
 K->status=K->forced?K->forced:ef_map(K->native);
 if(K->status==LIS_RCI_DONE){lis_vector_copy(K->stage,K->output);*K->value=K->stage_value;}
 else {e->nev=0;e->retcode=K->native;}
 lis_vector_destroy(K->stage);K->stage=NULL;
 lis_rci_native_session_destroy(K->child);K->child=NULL;
 K->active=0;memset(&K->q,0,sizeof K->q);return K->status;
}
static int ef_drive(lis_rci_native_feast_lifecycle *K,LIS_INT rc,lis_rci_native_session_request *q)
{
 while(rc==EF_AGAIN)rc=ef_step(K);
 if(rc!=EF_YIELD){int st=ef_finish(K,rc);if(q)*q=K->q;return st;}
 if(K->token==UINT64_MAX){
  lis_rci_native_session_cancel(K->child,K->inner.token);K->forced=LIS_RCI_PROTOCOL;K->completion=LIS_BREAKDOWN;
  return ef_drive(K,ef_step(K),q);
 }
 K->q.token=++K->token;K->status=(int)K->q.operation;
 if(q)*q=K->q;return K->status;
}
lis_rci_native_feast_lifecycle *lis_rci_native_feast_lifecycle_create(void)
{lis_rci_native_feast_lifecycle *K=calloc(1,sizeof(*K));if(K)K->status=LIS_RCI_STATE;return K;}
void lis_rci_native_feast_lifecycle_destroy(lis_rci_native_feast_lifecycle *K)
{if(K){if(K->active)lis_rci_native_feast_lifecycle_cancel(K,K->token);free(K);}}
int lis_rci_native_feast_lifecycle_start(lis_rci_native_feast_lifecycle *K,LIS_MATRIX A,LIS_MATRIX B,LIS_VECTOR x,LIS_SCALAR *value,LIS_ESOLVER e,lis_rci_native_session_request *q)
{
 if(!K||!A||!x||!value||!e||!q)return LIS_RCI_BADARG;
 if(K->active)return LIS_RCI_STATE;
 if(A->matrix_type!=LIS_MATRIX_CSR||A->n<=0||!A->ptr||!A->value||!x->value||x->n!=A->n||
 (B&&(B->matrix_type!=LIS_MATRIX_CSR||B->n!=A->n||!B->ptr||!B->value))||
 e->options[LIS_EOPTIONS_ESOLVER]!=LIS_ESOLVER_FEAST||e->options[LIS_EOPTIONS_STORAGE]||
 e->options[LIS_EOPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE||e->transform!=LIS_EIGEN_TRANSFORM_AUTO||
 e->subspace_retry_max||e->options[LIS_EOPTIONS_SUBSPACE]<=0||e->options[LIS_EOPTIONS_SUBSPACE]>A->n||
 e->options[LIS_EOPTIONS_MODE]<0||e->options[LIS_EOPTIONS_MODE]>=e->options[LIS_EOPTIONS_SUBSPACE]||
 e->options[LIS_EOPTIONS_MAXITER]<0)return LIS_RCI_BADARG;
 if(e->requested_nev_active&&(e->requested_nev<0||e->requested_nev>e->options[LIS_EOPTIONS_SUBSPACE]))return LIS_RCI_BADARG;
 uint64_t token=K->token;memset(K,0,sizeof(*K));K->token=token;K->status=LIS_RCI_STATE;
 K->child=lis_rci_native_session_create(LIS_NR_FEAST);if(!K->child)return LIS_RCI_BADARG;
 if(lis_vector_duplicate(A,&K->stage)){lis_rci_native_session_destroy(K->child);K->child=NULL;return LIS_RCI_BADARG;}
 lis_vector_copy(x,K->stage);K->output=x;K->value=value;
 K->f.A=A;K->f.B=B;K->f.x=K->stage;K->f.evalue0=&K->stage_value;K->f.esolver=e;
 if(e->requested_nev_active)e->nev=e->requested_nev;
 e->subspace_retry_attempts=0;e->subspace_retry_final_size=e->options[LIS_EOPTIONS_SUBSPACE];
 K->active=1;return ef_drive(K,ef_step(K),q);
}
int lis_rci_native_feast_lifecycle_advance(lis_rci_native_feast_lifecycle *K,uint64_t token,LIS_INT code,lis_rci_native_session_request *q)
{
 if(!K||!q)return LIS_RCI_BADARG;if(!K->active)return LIS_RCI_STATE;if(token!=K->token)return LIS_RCI_PROTOCOL;
 if(code){
  lis_rci_native_session_cancel(K->child,K->inner.token);
  K->forced=LIS_RCI_CALLBACK;K->completion=code;
  return ef_drive(K,ef_step(K),q);
 }
 return ef_drive(K,ef_core_result(K,lis_rci_native_session_advance(K->child,K->inner.token,code,&K->inner)),q);
}
int lis_rci_native_feast_lifecycle_cancel(lis_rci_native_feast_lifecycle *K,uint64_t token)
{
 if(!K)return LIS_RCI_BADARG;if(!K->active)return LIS_RCI_STATE;if(token!=K->token)return LIS_RCI_PROTOCOL;
 return ef_drive(K,ef_core_result(K,lis_rci_native_session_cancel(K->child,K->inner.token)),NULL);
}
int lis_rci_native_feast_lifecycle_status(const lis_rci_native_feast_lifecycle *K,LIS_INT *native)
{if(!K||!native)return LIS_RCI_BADARG;*native=K->native;return K->status;}
