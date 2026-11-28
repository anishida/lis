/* AI-assisted native continuation candidate, qualification only. */
#include "lis_rci_native_feast_projector.h"
#if defined(USE_MPI)||defined(USE_QUAD_PRECISION)||defined(_COMPLEX)||defined(_LONG__DOUBLE)
#error "native FEAST projector continuation requires real double non-MPI"
#endif
#include "feast_projector_continuation.inc"

static int fp_finish(lis_rci_native_feast_projector *F,LIS_INT rc)
{
 F->active=0;F->native=rc;
 F->status=F->forced?F->forced:rc==LIS_SUCCESS?FP_DONE:FP_ERROR;
 memset(&F->request,0,sizeof F->request);
 return F->status;
}
static int fp_finite(LIS_VECTOR x)
{
 if(!x||!x->value)return 0;
 for(LIS_INT i=0;i<x->n;i++)if(!isfinite(x->value[i]))return 0;
 return 1;
}
static int fp_suspend(lis_rci_native_feast_projector *F,int op,LIS_MATRIX A,LIS_MATRIX B,
 LIS_REAL zr,LIS_REAL zi,LIS_VECTOR ir,LIS_VECTOR ii,LIS_VECTOR or,LIS_VECTOR oi,LIS_SOLVER solver)
{
 if(F->token==UINT64_MAX){F->forced=FP_PROTOCOL;F->completion=LIS_BREAKDOWN;return 3;}
 if(F->strict){
  if(!fp_finite(ir)||(ii&&!fp_finite(ii))){F->forced=FP_NONFINITE;F->completion=LIS_BREAKDOWN;return 3;}
  for(LIS_INT i=0;i<or->n;i++)or->value[i]=NAN;
  if(oi)for(LIS_INT i=0;i<oi->n;i++)oi->value[i]=NAN;
 }
 F->request=(lis_rci_native_feast_projector_request){++F->token,op,A,B,zr,zi,ir,ii,or,oi,solver};
 return F->status=op;
}
static int fp_drive(lis_rci_native_feast_projector *F,lis_rci_native_feast_projector_request *q)
{
 int st;
 do{st=fp_step(F);}while(st==3);
 if(q)*q=F->request;
 return st;
}
lis_rci_native_feast_projector *lis_rci_native_feast_projector_create(void)
{
 lis_rci_native_feast_projector *F=calloc(1,sizeof(*F));
 if(F)F->status=FP_STATE;
 return F;
}
void lis_rci_native_feast_projector_destroy(lis_rci_native_feast_projector *F)
{
 if(F){if(F->active)(void)lis_rci_native_feast_projector_cancel(F,F->token);free(F);}
}
int lis_rci_native_feast_projector_start(lis_rci_native_feast_projector *F,
 LIS_MATRIX A,LIS_MATRIX B,LIS_REAL emin,LIS_REAL emax,LIS_REAL ratio,LIS_INT npoints,
 LIS_VECTOR xr,LIS_VECTOR xi,LIS_VECTOR qr,LIS_VECTOR qi,LIS_SOLVER solver,int strict,
 lis_rci_native_feast_projector_request *q)
{
 if(!F||!q||(strict!=0&&strict!=1))return FP_BADARG;
 if(F->active)return FP_STATE;
 uint64_t token=F->token;memset(F,0,sizeof(*F));F->token=token;
 F->active=1;F->strict=strict;
 F->A=A;F->B=B;F->emin=emin;F->emax=emax;F->ratio=ratio;F->npoints=npoints;
 F->x_real=xr;F->x_imag=xi;F->q_real=qr;F->q_imag=qi;F->solver=solver;
 return fp_drive(F,q);
}
int lis_rci_native_feast_projector_advance(lis_rci_native_feast_projector *F,uint64_t token,LIS_INT rc,lis_rci_native_feast_projector_request *q)
{
 if(!F||!q)return FP_BADARG;
 if(!F->active||!F->request.operation)return FP_STATE;
 if(token!=F->token)return FP_PROTOCOL;
 if(F->strict&&!rc&&(!fp_finite(F->request.out_real)||(F->request.out_imag&&!fp_finite(F->request.out_imag)))){
  F->forced=FP_NONFINITE;rc=LIS_BREAKDOWN;
 }
 F->completion=rc;memset(&F->request,0,sizeof F->request);
 return fp_drive(F,q);
}
int lis_rci_native_feast_projector_cancel(lis_rci_native_feast_projector *F,uint64_t token)
{
 if(!F)return FP_BADARG;
 if(!F->active||!F->request.operation)return FP_STATE;
 if(token!=F->token)return FP_PROTOCOL;
 F->forced=FP_CANCELLED;F->completion=LIS_BREAKDOWN;memset(&F->request,0,sizeof F->request);
 return fp_drive(F,NULL);
}
LIS_INT lis_rci_native_feast_projector_native(const lis_rci_native_feast_projector *F)
{return F?F->native:LIS_ERR_ILL_ARG;}
LIS_INT lis_efeast_projector_apply_split(LIS_MATRIX A,LIS_MATRIX B,LIS_REAL emin,LIS_REAL emax,
 LIS_REAL ratio,LIS_INT npoints,LIS_VECTOR xr,LIS_VECTOR xi,LIS_VECTOR qr,LIS_VECTOR qi,LIS_SOLVER solver)
{
 lis_rci_native_feast_projector *F=lis_rci_native_feast_projector_create();
 if(!F)return LIS_OUT_OF_MEMORY;
 lis_rci_native_feast_projector_request q;
 int st=lis_rci_native_feast_projector_start(F,A,B,emin,emax,ratio,npoints,xr,xi,qr,qi,solver,0,&q);
 while(st>0){
  LIS_INT e=q.operation==FP_B?lis_matvec(q.A,q.in_real,q.out_real):
   lis_efeast_shifted_solve_split(q.A,q.B,q.zr,q.zi,q.in_real,q.in_imag,q.out_real,q.out_imag,q.solver);
  st=lis_rci_native_feast_projector_advance(F,q.token,e,&q);
 }
 LIS_INT rc=lis_rci_native_feast_projector_native(F);
 lis_rci_native_feast_projector_destroy(F);
 return rc;
}
