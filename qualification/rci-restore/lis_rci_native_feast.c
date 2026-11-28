/* Qualification-only native full FEAST continuation. */
#include "lis_rci_native_feast.h"
#if defined(USE_MPI)||defined(USE_QUAD_PRECISION)
#error "native FEAST continuation candidate requires non-MPI without quad solver mode"
#endif
#include "feast_native_continuation.inc"
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
static int nf_return(lis_rci_native_feast *K,int parent,LIS_INT rc)
{
 K->active=parent;K->completion=rc;
 if(parent)return 3;
 if(K->forced){rc=LIS_BREAKDOWN;K->esolver->nev=0;K->esolver->retcode=rc;}
 K->native=rc;K->status=K->forced?K->forced:rc==LIS_SUCCESS?NF_DONE:NF_ERROR;
 memset(&K->request,0,sizeof K->request);
 return K->status;
}
static int nf_finite(LIS_VECTOR x)
{
 if(!x||!x->value)return 0;
 for(LIS_INT i=0;i<x->n;i++)if(!native_scalar_finite(x->value[i]))return 0;
 return 1;
}
static void nf_poison(LIS_VECTOR x){for(LIS_INT i=0;i<x->n;i++)x->value[i]=native_poison_value();}
static int nf_force(lis_rci_native_feast *K,int status)
{K->forced=status;K->completion=LIS_BREAKDOWN;return 3;}
static int nf_generation_ok(const lis_rci_native_feast *K)
{
 if(K->esolver->A!=K->owner_a||K->esolver->B!=K->owner_b)return 0;
 if(K->owner_a&&(K->owner_a->structure_generation!=K->as||K->owner_a->value_generation!=K->av))return 0;
 if(K->owner_b&&(K->owner_b->structure_generation!=K->bs||K->owner_b->value_generation!=K->bv))return 0;
 return 1;
}
static int nf_block(lis_rci_native_feast *K,int op,LIS_MATRIX A,LIS_INT n,LIS_VECTOR *in,LIS_VECTOR *out)
{
 if(K->forced)return nf_force(K,K->forced);
 if(K->strict&&!nf_generation_ok(K))return nf_force(K,NF_STALE);
 if(K->token==UINT64_MAX)return nf_force(K,NF_PROTOCOL);
 if(K->strict){
  if(n<=0||!in||!out)return nf_force(K,NF_BADARG);
  for(LIS_INT i=0;i<n;i++)if(!nf_finite(in[i])||!out[i]||in[i]->n!=out[i]->n)return nf_force(K,NF_NONFINITE);
  for(LIS_INT i=0;i<n;i++)for(LIS_INT j=0;j<n;j++)if(in[i]==out[j]||(i!=j&&out[i]==out[j]))return nf_force(K,NF_BADARG);
  for(LIS_INT i=0;i<n;i++)nf_poison(out[i]);
 }
 K->request=(lis_rci_native_feast_request){.token=++K->token,.operation=op,.A=A,.count=n,.block=1,.in=in,.out=out};
 return K->status=op;
}
static int nf_one(lis_rci_native_feast *K,int op,LIS_MATRIX A,LIS_VECTOR in,LIS_VECTOR out)
{
 K->one_in=in;K->one_out=out;
 int st=nf_block(K,op,A,1,&K->one_in,&K->one_out);
 if(st==op)K->request.block=0;
 return st;
}
static int nf_shift(lis_rci_native_feast *K,LIS_MATRIX A,LIS_MATRIX B,LIS_REAL zr,LIS_REAL zi,
 LIS_VECTOR ir,LIS_VECTOR ii,LIS_VECTOR or,LIS_VECTOR oi,LIS_SOLVER s)
{
 if(K->forced)return nf_force(K,K->forced);
 if(K->strict&&!nf_generation_ok(K))return nf_force(K,NF_STALE);
 if(K->token==UINT64_MAX)return nf_force(K,NF_PROTOCOL);
 if(K->strict){
  if(!nf_finite(ir)||!nf_finite(ii)||!or||!oi||or==oi||or==ir||or==ii||oi==ir||oi==ii)return nf_force(K,NF_NONFINITE);
  nf_poison(or);nf_poison(oi);
 }
 K->request=(lis_rci_native_feast_request){.token=++K->token,.operation=NF_SHIFT,.A=A,.B=B,.zr=zr,.zi=zi,
 .in_real=ir,.in_imag=ii,.out_real=or,.out_imag=oi,.solver=s};
 return K->status=NF_SHIFT;
}
static int nf_drive(lis_rci_native_feast *K,lis_rci_native_feast_request *q)
{
 int st;do{st=nf_step(K);}while(st==3);
 if(q)*q=K->request;
 return st;
}
lis_rci_native_feast *lis_rci_native_feast_create(void)
{
 lis_rci_native_feast *K=calloc(1,sizeof(*K));if(K)K->status=NF_STATE;return K;
}
void lis_rci_native_feast_destroy(lis_rci_native_feast *K)
{if(K){if(K->active)(void)lis_rci_native_feast_cancel(K,K->token);free(K);}}
int lis_rci_native_feast_start(lis_rci_native_feast *K,LIS_ESOLVER e,int strict,lis_rci_native_feast_request *q)
{
 if(!K||!e||!q||(strict!=0&&strict!=1))return NF_BADARG;
 if(K->active)return NF_STATE;
 uint64_t token=K->token;memset(K,0,sizeof(*K));K->token=token;
 K->esolver=e;K->feast.esolver=e;K->strict=strict;K->active=9;
 K->owner_a=e->A;K->owner_b=e->B;
 if(e->A){K->as=e->A->structure_generation;K->av=e->A->value_generation;}
 if(e->B){K->bs=e->B->structure_generation;K->bv=e->B->value_generation;}
 return nf_drive(K,q);
}
int lis_rci_native_feast_advance(lis_rci_native_feast *K,uint64_t token,LIS_INT rc,lis_rci_native_feast_request *q)
{
 if(!K||!q)return NF_BADARG;
 if(!K->active||!K->request.operation)return NF_STATE;
 if(token!=K->token)return NF_PROTOCOL;
 if(K->strict&&!nf_generation_ok(K)){K->forced=NF_STALE;rc=LIS_BREAKDOWN;}
 if(K->strict&&!rc){
  int ok=1;
  if(K->request.operation==NF_SHIFT)ok=nf_finite(K->request.out_real)&&nf_finite(K->request.out_imag);
  else for(LIS_INT i=0;i<K->request.count;i++)if(!nf_finite(K->request.out[i])){ok=0;break;}
  if(!ok){K->forced=NF_NONFINITE;rc=LIS_BREAKDOWN;}
 }
 K->completion=rc;memset(&K->request,0,sizeof K->request);return nf_drive(K,q);
}
int lis_rci_native_feast_cancel(lis_rci_native_feast *K,uint64_t token)
{
 if(!K)return NF_BADARG;
 if(!K->active||!K->request.operation)return NF_STATE;
 if(token!=K->token)return NF_PROTOCOL;
 K->forced=NF_CANCELLED;K->completion=LIS_BREAKDOWN;memset(&K->request,0,sizeof K->request);return nf_drive(K,NULL);
}
LIS_INT lis_rci_native_feast_native(const lis_rci_native_feast *K){return K?K->native:LIS_ERR_ILL_ARG;}
static _Thread_local int nf_reference;
static _Thread_local lis_rci_native_feast_test_driver nf_driver;
void lis_rci_native_feast_test_select(int ref,lis_rci_native_feast_test_driver driver){nf_reference=ref;nf_driver=driver;}
LIS_INT lis_efeast(LIS_ESOLVER e)
{
 if(nf_reference)return lis_efeast_native_reference(e);
 if(nf_driver)return nf_driver(e);
 lis_rci_native_feast *K=lis_rci_native_feast_create();if(!K)return LIS_OUT_OF_MEMORY;
 lis_rci_native_feast_request q;
 int st=lis_rci_native_feast_start(K,e,0,&q);
 while(st>0){
  LIS_INT rc=q.operation==NF_SHIFT?lis_efeast_shifted_solve_split(q.A,q.B,q.zr,q.zi,q.in_real,q.in_imag,q.out_real,q.out_imag,q.solver):
    (!q.block?lis_matvec(q.A,q.in[0],q.out[0]):lis_matvec_block(q.A,q.count,q.in,q.out));
  st=lis_rci_native_feast_advance(K,q.token,rc,&q);
 }
 LIS_INT rc=st==NF_BADARG?LIS_ERR_ILL_ARG:lis_rci_native_feast_native(K);
 lis_rci_native_feast_destroy(K);return rc;
}
