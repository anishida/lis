/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#include "lis_rci.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
/* These optimizations invalidate the nonfinite-output protocol. */
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0)
#error "experimental RCI requires IEEE nonfinite checks (disable fast/finite-only math)"
#endif

enum phase {
 IDLE, B_NORM, INITIAL, VERIFY, RES_NORM,
 CG_P, CG_RZ, CG_A, CG_PV, CG_NORM,
 BI_P, BI_RHO_FIRST, BI_A, BI_SV, BI_NORM_S, BI_S, BI_T,
 BI_TT, BI_TS, BI_NORM_R, BI_RHO_NEXT,
 GM_P, GM_A, GM_BEFORE, GM_DOT, GM_AFTER,
 RES_VECTOR, CG_P_READY, CG_X_READY, CG_R_READY,
 BI_SHADOW_READY, BI_P_READY, BI_S_READY, BI_SHORT_READY,
 BI_X1_READY, BI_X2_READY, BI_R_READY, BI_TMP_READY, BI_NEW_P_READY,
 GM_V_READY, GM_ORTH_READY, GM_NEXT_V_READY, GM_UPDATE_READY, END
};
struct lis_rci_context {
 lis_rci_options o; lis_rci_info info; enum phase phase; int verify_reason;
 double *mem,*b,*x,*r,*z,*p,*v,*s,*t,*shadow,*ax,*V,*Z,*H,*cs,*sn,*g,*y,*out,*scratch,*lin_dest;
 size_t j, orth_i, update_i, update_k; long double tol,rho,alpha,omega,work,before,lin_a,lin_b,residual_wide;
 int fresh, orth_pass, external, update_reason, happy;
 uint64_t token;
 lis_rci_request_ex pending;
 lis_rci_stats stats;
};
static int mul(size_t a,size_t b,size_t *z){if(b && a>SIZE_MAX/b)return 0;*z=a*b;return 1;}
static int add(size_t a,size_t b,size_t *z){if(a>SIZE_MAX-b)return 0;*z=a+b;return 1;}
/* Avoid serialized OpenMP regions below the existing parallel threshold.
 * Keep the numerical reduction order within each selected backend. */
#define RCI_OMP_MIN 4096u
static int finite_vec(const double *a,size_t n)
{
#ifdef _OPENMP
 if(n>=RCI_OMP_MIN){
  int ok=1;
#pragma omp parallel for reduction(&:ok)
  for(size_t i=0;i<n;i++)ok &= isfinite(a[i])!=0;
  return ok;
 }
#endif
 for(size_t i=0;i<n;i++)if(!isfinite(a[i]))return 0;
 return 1;
}
static long double dot(const double *a,const double *b,size_t n)
{
 long double s=0;
#ifdef _OPENMP
 if(n>=RCI_OMP_MIN){
#pragma omp parallel for reduction(+:s)
  for(size_t i=0;i<n;i++)s+=(long double)a[i]*b[i];
  return s;
 }
#endif
 for(size_t i=0;i<n;i++)s+=(long double)a[i]*b[i];
 return s;
}
static long double norm(const double *a,size_t n)
{
 /* Scale before squaring, including when long double has narrow range. */
 double m=0;long double s=0;
#ifdef _OPENMP
 if(n>=RCI_OMP_MIN){
#pragma omp parallel for reduction(max:m)
  for(size_t i=0;i<n;i++)if(fabs(a[i])>m)m=fabs(a[i]);
  if(m==0)return 0;
#pragma omp parallel for reduction(+:s)
  for(size_t i=0;i<n;i++){long double v=(long double)a[i]/m;s+=v*v;}
  return (long double)m*sqrtl(s);
 }
#endif
 for(size_t i=0;i<n;i++)if(fabs(a[i])>m)m=fabs(a[i]);
 if(m==0)return 0;
 for(size_t i=0;i<n;i++){long double v=(long double)a[i]/m;s+=v*v;}
 return (long double)m*sqrtl(s);
}
/* Do not cast a finite long-double intermediate outside the double range.
 * Existing request validation rejects the resulting nonfinite value. */
static double narrow(long double value)
{
 if(value>DBL_MAX)return INFINITY;
 if(value< -DBL_MAX)return -INFINITY;
 return (double)value;
}
static void copy(double *a,const double *b,size_t n){memcpy(a,b,n*sizeof(double));}
static void axpy(double *a,const double *b,long double k,size_t n)
{
#ifdef _OPENMP
 if(n>=RCI_OMP_MIN){
#pragma omp parallel for
  for(size_t i=0;i<n;i++)a[i]=narrow((long double)a[i]+k*b[i]);
  return;
 }
#endif
 for(size_t i=0;i<n;i++)a[i]=narrow((long double)a[i]+k*b[i]);
}
/* Request tokens are context-local and never reset on context reuse. */
static int end(lis_rci_context *c,int status,lis_rci_request_ex *q)
{
 c->phase=END; c->info.status=status; c->out=NULL; c->lin_dest=NULL;
 c->pending=(lis_rci_request_ex){0};
 *q=(lis_rci_request_ex){0}; return status;
}
static int request(lis_rci_context *c,enum phase p,int status,
 const double *in,const double *in2,double *out,lis_rci_request_ex *q)
{
 if(!finite_vec(in,c->o.n)||(in2&&in2!=in&&!finite_vec(in2,c->o.n)))
  return end(c,LIS_RCI_NONFINITE,q);
 if(c->token==UINT64_MAX)return end(c,LIS_RCI_PROTOCOL,q);
 if(out)for(size_t i=0;i<c->o.n;i++)out[i]=NAN;
 c->phase=p;c->out=out;c->info.status=status;
 *q=(lis_rci_request_ex){++c->token,status,c->o.n,in,in2,out};
 c->pending=*q;c->stats.issued[status-1]++;c->stats.last_token=c->token;
 if(status==LIS_RCI_MATVEC)c->info.matvecs++;
 if(status==LIS_RCI_PRECON)c->info.precons++;
 return status;
}
static int req(lis_rci_context *c,enum phase p,int status,
 const double *in,double *out,lis_rci_request_ex *q)
{return request(c,p,status,in,NULL,out,q);}
static int scalar_req(lis_rci_context *c,enum phase p,const double *in,
 const double *in2,lis_rci_request_ex *q)
{return request(c,p,in2?LIS_RCI_DOT:LIS_RCI_NORM,in,in2,NULL,q);}
/* Full mode stages vector output separately so a callback never overwrites
 * an input operand, even when the algorithm updates a vector in place. */
static int linear_req(lis_rci_context *c,enum phase p,int op,double *dest,
 const double *a,const double *b,long double alpha,long double beta,
 lis_rci_request_ex *q)
{
 if(!isfinite(alpha)||!isfinite(beta)||(op==LIS_RCI_DIVIDE&&alpha<=0))
  return end(c,LIS_RCI_NONFINITE,q);
 c->lin_dest=dest;c->lin_a=alpha;c->lin_b=beta;
 return request(c,p,op,a,b,c->scratch,q);
}
static int combine(lis_rci_context *c,enum phase p,double *dest,
 const double *a,const double *b,long double alpha,long double beta,lis_rci_request_ex *q)
{return linear_req(c,p,LIS_RCI_AXPBY,dest,a,b,alpha,beta,q);}
static int divide_req(lis_rci_context *c,enum phase p,double *dest,const double *a,
 long double divisor,lis_rci_request_ex *q)
{return linear_req(c,p,LIS_RCI_DIVIDE,dest,a,NULL,divisor,0,q);}
static int verify(lis_rci_context *c,int reason,lis_rci_request_ex *q)
{c->verify_reason=reason;return req(c,VERIFY,LIS_RCI_MATVEC,c->x,c->ax,q);}
lis_rci_context *lis_rci_create(const lis_rci_options*o){
 if(!o||!o->n||!o->maxiter||o->method<LIS_RCI_CG||o->method>LIS_RCI_FGMRES||!isfinite(o->rtol)||!isfinite(o->atol)||o->rtol<0||o->atol<0||(!o->rtol&&!o->atol))return NULL;
 size_t n=o->n,m=0,total,tmp,mp=0;if(!mul(n,11,&total))return NULL;
 if(o->method>=LIS_RCI_GMRES){if(!o->restart)return NULL;m=o->restart;if(m>o->maxiter)m=o->maxiter;if(m>n)m=n;if(!add(m,1,&mp)||!mul(mp,n,&tmp)||!add(total,tmp,&total)||!mul(m,n,&tmp)||!add(total,tmp,&total)||!mul(mp,m,&tmp)||!add(total,tmp,&total)||!mul(m,4,&tmp)||!add(total,tmp,&total)||!add(total,1,&total))return NULL;}
 if(!mul(total,sizeof(double),&tmp))return NULL;
 lis_rci_context*c=calloc(1,sizeof(*c));if(!c)return NULL;c->mem=calloc(total,sizeof(double));if(!c->mem){free(c);return NULL;}c->o=*o;c->o.restart=m;c->info.status=LIS_RCI_STATE;c->info.residual=NAN;c->residual_wide=NAN;c->stats.error_rank=UINT32_MAX;
 double *d=c->mem;
 c->b=d;d+=n;c->x=d;d+=n;c->r=d;d+=n;c->z=d;d+=n;c->p=d;d+=n;c->v=d;d+=n;c->s=d;d+=n;c->t=d;d+=n;c->shadow=d;d+=n;c->ax=d;d+=n;c->scratch=d;d+=n;
 if(m){c->V=d;d+=mp*n;c->Z=d;d+=m*n;c->H=d;d+=mp*m;c->cs=d;d+=m;c->sn=d;d+=m;c->g=d;d+=mp;c->y=d;}
 return c;
}
void lis_rci_destroy(lis_rci_context*c){if(c){free(c->mem);free(c);}}
const double*lis_rci_solution(const lis_rci_context*c){return c?c->x:NULL;}
int lis_rci_get_info(const lis_rci_context*c,lis_rci_info*i){if(!c||!i)return LIS_RCI_BADARG;*i=c->info;return 0;}
int lis_rci_get_residual(const lis_rci_context *c,long double *residual)
{
 if(!c||!residual)return LIS_RCI_BADARG;
 *residual=c->residual_wide;
 return LIS_RCI_DONE;
}
static int cycle(lis_rci_context *c,long double rn,lis_rci_request_ex *q)
{
 size_t n=c->o.n,m=c->o.restart;
 if(rn<=c->tol)return end(c,LIS_RCI_DONE,q);
 if(c->info.iterations>=c->o.maxiter)return end(c,LIS_RCI_MAXITER,q);
 c->fresh=1;c->rho=c->alpha=c->omega=1;
 if(c->o.method==LIS_RCI_CG)return req(c,CG_P,LIS_RCI_PRECON,c->r,c->z,q);
 if(c->o.method==LIS_RCI_BICGSTAB){
  if(c->external==2)return combine(c,BI_SHADOW_READY,c->shadow,c->r,NULL,1,0,q);
  copy(c->shadow,c->r,n);copy(c->p,c->r,n);
  return req(c,BI_P,LIS_RCI_PRECON,c->p,c->z,q);
 }
 if(rn>DBL_MAX)return end(c,LIS_RCI_NONFINITE,q);
 c->j=0;memset(c->H,0,(m+1)*m*sizeof(double));
 memset(c->g,0,(m+1)*sizeof(double));c->g[0]=(double)rn;
 if(c->external==2)return divide_req(c,GM_V_READY,c->V,c->r,rn,q);
 for(size_t i=0;i<n;i++)c->V[i]=narrow((long double)c->r[i]/rn);
 return req(c,GM_P,LIS_RCI_PRECON,c->V,c->Z,q);
}
static int start_impl(lis_rci_context *c,const double *b,const double *x0,
 int external,lis_rci_request_ex *q)
{
 if(!c||!b||!q)return LIS_RCI_BADARG;
 if(c->phase!=IDLE&&c->phase!=END)return LIS_RCI_STATE;
 size_t n=c->o.n;
 c->info=(lis_rci_info){.residual=NAN};c->residual_wide=NAN;c->external=external;c->verify_reason=0;
 c->stats=(lis_rci_stats){.last_token=c->token,.error_rank=UINT32_MAX};
 if(!finite_vec(b,n)||(x0&&!finite_vec(x0,n)))return end(c,LIS_RCI_NONFINITE,q);
 memmove(c->b,b,n*sizeof(double));
 if(x0)memmove(c->x,x0,n*sizeof(double));else memset(c->x,0,n*sizeof(double));
 return scalar_req(c,B_NORM,c->b,NULL,q);
}
static int gm_update(lis_rci_context*c,size_t k,int reason,lis_rci_request_ex*q){size_t n=c->o.n,m=c->o.restart;
 for(size_t ii=k;ii>0;ii--){size_t i=ii-1;long double a=c->g[i];for(size_t j=i+1;j<k;j++)a-=(long double)c->H[j*(m+1)+i]*c->y[j];double d=c->H[i*(m+1)+i];if(d==0||!isfinite(d))return verify(c,LIS_RCI_BREAKDOWN,q);c->y[i]=narrow(a/d);if(!isfinite(c->y[i]))return verify(c,LIS_RCI_BREAKDOWN,q);}
 if(c->external==2){
  c->update_i=0;c->update_k=k;c->update_reason=reason;
  return combine(c,GM_UPDATE_READY,c->x,c->x,c->Z,1,c->y[0],q);
 }
 for(size_t j=0;j<k;j++)axpy(c->x,c->Z+j*n,c->y[j],n);return verify(c,reason,q);}
static int gm_finish_column(lis_rci_context *c,lis_rci_request_ex *q)
{
 size_t j=c->j,m=c->o.restart,n=c->o.n;double *h=c->H+j*(m+1);
 for(size_t i=0;i<j;i++){double u=h[i],w=h[i+1];h[i]=c->cs[i]*u+c->sn[i]*w;h[i+1]=-c->sn[i]*u+c->cs[i]*w;}
 double d=hypot(h[j],h[j+1]);
 if(!isfinite(d)||d==0){if(j)return gm_update(c,j,LIS_RCI_BREAKDOWN,q);return verify(c,LIS_RCI_BREAKDOWN,q);}
 c->cs[j]=h[j]/d;c->sn[j]=h[j+1]/d;h[j]=d;h[j+1]=0;
 c->g[j+1]=-c->sn[j]*c->g[j];c->g[j]*=c->cs[j];c->info.iterations++;
 if(c->happy||fabs(c->g[j+1])<=c->tol||j+1==m||c->info.iterations>=c->o.maxiter)
  return gm_update(c,j+1,c->happy?LIS_RCI_BREAKDOWN:LIS_RCI_MAXITER,q);
 c->j++;return req(c,GM_P,LIS_RCI_PRECON,c->V+c->j*n,c->Z+c->j*n,q);
}
static int advance_impl(lis_rci_context *c,uint64_t token,int code,
 long double scalar,lis_rci_request_ex *q)
{
 if(!c||!q)return LIS_RCI_BADARG;
 if(c->phase==IDLE||c->phase==END)return LIS_RCI_STATE;
 /* A stale or duplicate completion must leave the outstanding request intact. */
 if(token!=c->token)return LIS_RCI_PROTOCOL;
 c->stats.completed[c->info.status-1]++;
 if(code){c->info.callback_code=code;return end(c,LIS_RCI_CALLBACK,q);}
 if(c->info.status==LIS_RCI_DOT||c->info.status==LIS_RCI_NORM){
  if(!isfinite(scalar)||(c->info.status==LIS_RCI_NORM&&scalar<0))
   return end(c,LIS_RCI_NONFINITE,q);
 }else if(!finite_vec(c->out,c->o.n))return end(c,LIS_RCI_NONFINITE,q);
 if(c->info.status==LIS_RCI_AXPBY||c->info.status==LIS_RCI_DIVIDE){
  copy(c->lin_dest,c->scratch,c->o.n);c->lin_dest=NULL;
 }
 size_t n=c->o.n,m=c->o.restart;long double a,b,rn;
 switch(c->phase){
 case B_NORM:
  c->tol=fmaxl(c->o.atol,(long double)c->o.rtol*scalar);
  if(!isfinite(c->tol))return end(c,LIS_RCI_NONFINITE,q);
  return req(c,INITIAL,LIS_RCI_MATVEC,c->x,c->ax,q);
 case INITIAL:case VERIFY:
  if(c->external==2)return combine(c,RES_VECTOR,c->r,c->b,c->ax,1,-1,q);
  for(size_t i=0;i<n;i++)c->r[i]=c->b[i]-c->ax[i];
  /* fall through */
 case RES_VECTOR:return scalar_req(c,RES_NORM,c->r,NULL,q);
 case RES_NORM:
  c->residual_wide=scalar;
  /* Preserve the actual convergence decision even if the legacy double
   * diagnostic cannot represent this norm. Avoid an out-of-range cast. */
  c->info.residual=scalar>DBL_MAX?INFINITY:(double)scalar;
  if(scalar<=c->tol)return end(c,LIS_RCI_DONE,q);
  if(c->verify_reason==LIS_RCI_BREAKDOWN)return end(c,LIS_RCI_BREAKDOWN,q);
  return cycle(c,scalar,q);
 case CG_P:return scalar_req(c,CG_RZ,c->r,c->z,q);
 case CG_RZ:
  a=scalar;if(a<=0)return verify(c,LIS_RCI_BREAKDOWN,q);
  if(c->external==2){
   b=c->fresh?0:a/c->rho;c->fresh=0;c->rho=a;
   return combine(c,CG_P_READY,c->p,c->z,b==0?NULL:c->p,1,b,q);
  }
  if(c->fresh){copy(c->p,c->z,n);c->fresh=0;}
  else{b=a/c->rho;for(size_t i=0;i<n;i++)c->p[i]=narrow((long double)c->z[i]+b*c->p[i]);}
  c->rho=a;
  /* fall through */
 case CG_P_READY:return req(c,CG_A,LIS_RCI_MATVEC,c->p,c->v,q);
 case CG_A:return scalar_req(c,CG_PV,c->p,c->v,q);
 case CG_PV:
  if(scalar<=0)return verify(c,LIS_RCI_BREAKDOWN,q);
  b=c->rho/scalar;if(!isfinite(b))return verify(c,LIS_RCI_BREAKDOWN,q);
  if(c->external==2){c->work=b;return combine(c,CG_X_READY,c->x,c->x,c->p,1,b,q);}
  axpy(c->x,c->p,b,n);axpy(c->r,c->v,-b,n);
  goto cg_residual;
 case CG_X_READY:return combine(c,CG_R_READY,c->r,c->r,c->v,1,-c->work,q);
 case CG_R_READY:
 cg_residual:
  c->info.iterations++;
  if(!finite_vec(c->x,n))return end(c,LIS_RCI_NONFINITE,q);
  return scalar_req(c,CG_NORM,c->r,NULL,q);
 case CG_NORM:
  if(scalar<=c->tol||c->info.iterations>=c->o.maxiter)return verify(c,LIS_RCI_MAXITER,q);
  return req(c,CG_P,LIS_RCI_PRECON,c->r,c->z,q);
 case BI_SHADOW_READY:return combine(c,BI_P_READY,c->p,c->r,NULL,1,0,q);
 case BI_P_READY:return req(c,BI_P,LIS_RCI_PRECON,c->p,c->z,q);
 case BI_P:
  if(c->fresh)return scalar_req(c,BI_RHO_FIRST,c->shadow,c->r,q);
  return req(c,BI_A,LIS_RCI_MATVEC,c->z,c->v,q);
 case BI_RHO_FIRST:
  c->rho=scalar;c->fresh=0;if(!scalar)return verify(c,LIS_RCI_BREAKDOWN,q);
  return req(c,BI_A,LIS_RCI_MATVEC,c->z,c->v,q);
 case BI_A:return scalar_req(c,BI_SV,c->shadow,c->v,q);
 case BI_SV:
  if(!scalar)return verify(c,LIS_RCI_BREAKDOWN,q);
  c->alpha=c->rho/scalar;if(!isfinite(c->alpha))return verify(c,LIS_RCI_BREAKDOWN,q);
  if(c->external==2)return combine(c,BI_S_READY,c->s,c->r,c->v,1,-c->alpha,q);
  copy(c->s,c->r,n);axpy(c->s,c->v,-c->alpha,n);
  /* fall through */
 case BI_S_READY:
  return scalar_req(c,BI_NORM_S,c->s,NULL,q);
 case BI_NORM_S:
  if(scalar<=c->tol){
   if(c->external==2)return combine(c,BI_SHORT_READY,c->x,c->x,c->z,1,c->alpha,q);
   axpy(c->x,c->z,c->alpha,n);c->info.iterations++;return verify(c,LIS_RCI_MAXITER,q);
  }
  return req(c,BI_S,LIS_RCI_PRECON,c->s,c->ax,q);
 case BI_SHORT_READY:c->info.iterations++;return verify(c,LIS_RCI_MAXITER,q);
 case BI_S:return req(c,BI_T,LIS_RCI_MATVEC,c->ax,c->t,q);
 case BI_T:return scalar_req(c,BI_TT,c->t,c->t,q);
 case BI_TT:
  if(scalar<=0)return verify(c,LIS_RCI_BREAKDOWN,q);
  c->work=scalar;return scalar_req(c,BI_TS,c->t,c->s,q);
 case BI_TS:
  c->omega=scalar/c->work;if(!c->omega||!isfinite(c->omega))return verify(c,LIS_RCI_BREAKDOWN,q);
  if(c->external==2)return combine(c,BI_X1_READY,c->x,c->x,c->z,1,c->alpha,q);
  axpy(c->x,c->z,c->alpha,n);axpy(c->x,c->ax,c->omega,n);
  copy(c->r,c->s,n);axpy(c->r,c->t,-c->omega,n);
  goto bi_residual;
 case BI_X1_READY:return combine(c,BI_X2_READY,c->x,c->x,c->ax,1,c->omega,q);
 case BI_X2_READY:return combine(c,BI_R_READY,c->r,c->s,c->t,1,-c->omega,q);
 case BI_R_READY:
 bi_residual:
  c->info.iterations++;
  if(!finite_vec(c->x,n))return end(c,LIS_RCI_NONFINITE,q);
  return scalar_req(c,BI_NORM_R,c->r,NULL,q);
 case BI_NORM_R:
  if(scalar<=c->tol||c->info.iterations>=c->o.maxiter)return verify(c,LIS_RCI_MAXITER,q);
  return scalar_req(c,BI_RHO_NEXT,c->shadow,c->r,q);
 case BI_RHO_NEXT:
  if(!scalar||!c->omega)return verify(c,LIS_RCI_BREAKDOWN,q);
  b=(scalar/c->rho)*(c->alpha/c->omega);
  if(!isfinite(b))return verify(c,LIS_RCI_BREAKDOWN,q);
  if(c->external==2){
   c->rho=scalar;c->work=b;
   return combine(c,BI_TMP_READY,c->s,c->p,c->v,1,-c->omega,q);
  }
  for(size_t i=0;i<n;i++)c->p[i]=narrow((long double)c->r[i]+b*((long double)c->p[i]-c->omega*c->v[i]));
  c->rho=scalar;return req(c,BI_P,LIS_RCI_PRECON,c->p,c->z,q);
 case BI_TMP_READY:return combine(c,BI_NEW_P_READY,c->p,c->r,c->s,1,c->work,q);
 case BI_NEW_P_READY:return req(c,BI_P,LIS_RCI_PRECON,c->p,c->z,q);
 case GM_V_READY:return req(c,GM_P,LIS_RCI_PRECON,c->V,c->Z,q);
 case GM_P:return req(c,GM_A,LIS_RCI_MATVEC,c->Z+c->j*n,c->v,q);
 case GM_A:return scalar_req(c,GM_BEFORE,c->v,NULL,q);
 case GM_BEFORE:
  c->before=scalar;c->orth_i=0;c->orth_pass=0;
  return scalar_req(c,GM_DOT,c->V,c->v,q);
 case GM_DOT:{
  double *h=c->H+c->j*(m+1);size_t i=c->orth_i;
  double d=narrow(scalar);
  if(!isfinite(d))return end(c,LIS_RCI_NONFINITE,q);
  h[i]+=d;
  if(!isfinite(h[i]))return end(c,LIS_RCI_NONFINITE,q);
  if(c->external==2)return combine(c,GM_ORTH_READY,c->v,c->v,c->V+i*n,1,-scalar,q);
  axpy(c->v,c->V+i*n,-scalar,n);
 }
 /* fall through */
 case GM_ORTH_READY:{
  if(++c->orth_i>c->j){c->orth_i=0;c->orth_pass++;}
  if(c->orth_pass<2)return scalar_req(c,GM_DOT,c->V+c->orth_i*n,c->v,q);
  return scalar_req(c,GM_AFTER,c->v,NULL,q);
 }
 case GM_AFTER:{
  size_t j=c->j;double *h=c->H+j*(m+1);rn=scalar;
  if(rn>DBL_MAX)return end(c,LIS_RCI_NONFINITE,q);
  h[j+1]=(double)rn;c->happy=rn<=DBL_EPSILON*c->before;
  if(!c->happy){
   if(c->external==2)return divide_req(c,GM_NEXT_V_READY,c->V+(j+1)*n,c->v,rn,q);
   for(size_t i=0;i<n;i++)c->V[(j+1)*n+i]=narrow((long double)c->v[i]/rn);
  }
  return gm_finish_column(c,q);
 }
 case GM_NEXT_V_READY:return gm_finish_column(c,q);
 case GM_UPDATE_READY:
  if(++c->update_i<c->update_k)
   return combine(c,GM_UPDATE_READY,c->x,c->x,c->Z+c->update_i*n,1,c->y[c->update_i],q);
  return verify(c,c->update_reason,q);
 default:return end(c,LIS_RCI_STATE,q);
 }
}
int lis_rci_start_ex(lis_rci_context *c,const double *b,const double *x0,lis_rci_request_ex *q)
{return start_impl(c,b,x0,1,q);}
int lis_rci_advance_ex(lis_rci_context *c,uint64_t token,int code,long double scalar,lis_rci_request_ex *q)
{
 if(!c||!q)return LIS_RCI_BADARG;
 if(c->external!=1)return LIS_RCI_STATE;
 return advance_impl(c,token,code,scalar,q);
}
/* Full-mode output uses a new structure; the old request ABI is unchanged. */
static void full_request(const lis_rci_context *c,int s,const lis_rci_request_ex *e,lis_rci_request_full *q)
{
 *q=(lis_rci_request_full){0};q->request=*e;
 if(s==LIS_RCI_AXPBY||s==LIS_RCI_DIVIDE){q->alpha=c->lin_a;q->beta=c->lin_b;}
}
int lis_rci_start_full(lis_rci_context *c,const double *b,const double *x0,lis_rci_request_full *q)
{
 if(!q)return LIS_RCI_BADARG;
 lis_rci_request_ex e={0};int s=start_impl(c,b,x0,2,&e);
 if(s!=LIS_RCI_BADARG&&s!=LIS_RCI_STATE)full_request(c,s,&e,q);
 return s;
}
int lis_rci_advance_full(lis_rci_context *c,uint64_t token,int code,long double scalar,lis_rci_request_full *q)
{
 if(!c||!q)return LIS_RCI_BADARG;
 if(c->external!=2)return LIS_RCI_STATE;
 lis_rci_request_ex e={0};int s=advance_impl(c,token,code,scalar,&e);
 /* A wrong token must not alter the caller's outstanding descriptor. */
 if(s!=LIS_RCI_BADARG&&s!=LIS_RCI_STATE&&(s!=LIS_RCI_PROTOCOL||c->phase==END))
  full_request(c,s,&e,q);
 return s;
}
/* Pending inspection/cancellation do not consume a numerical completion. */
int lis_rci_get_pending(const lis_rci_context *c,lis_rci_request_full *q)
{
 if(!c||!q)return LIS_RCI_BADARG;
 if(c->phase==IDLE||c->phase==END){*q=(lis_rci_request_full){0};return LIS_RCI_STATE;}
 full_request(c,c->info.status,&c->pending,q);
 return c->info.status;
}
int lis_rci_cancel(lis_rci_context *c,uint64_t token)
{
 if(!c)return LIS_RCI_BADARG;
 if(c->phase==IDLE||c->phase==END)return LIS_RCI_STATE;
 if(token!=c->token)return LIS_RCI_PROTOCOL;
 c->stats.cancelled=1;
 lis_rci_request_ex q;
 return end(c,LIS_RCI_CANCELLED,&q);
}
int lis_rci_get_stats(const lis_rci_context *c,lis_rci_stats *stats)
{
 if(!c||!stats)return LIS_RCI_BADARG;
 *stats=c->stats;return LIS_RCI_DONE;
}
/* Checked completion with the optional transport-neutral consensus record.
 * The tag is the pending operation. Protocol disagreement takes precedence
 * over callback errors; the selected callback error is still recorded. */
int lis_rci_advance_consensus(lis_rci_context *c,uint64_t token,
 const uint64_t words[5],long double scalar,lis_rci_request_full *q)
{
 if(!c||!words||!q)return LIS_RCI_BADARG;
 if(c->external==0||c->phase==IDLE||c->phase==END)return LIS_RCI_STATE;
 if(token!=c->token)return LIS_RCI_PROTOCOL;
 lis_rci_consensus_result result;
 int s=lis_rci_consensus_unpack(words,&result);
 if(s)return s;
 c->stats.error_rank=result.error_rank;
 if(result.token_mismatch||result.tag_mismatch||
    words[0]!=token||words[2]!=(uint64_t)c->info.status){
  c->stats.completed[c->info.status-1]++;
  c->info.callback_code=result.error;
  lis_rci_request_ex e;
  s=end(c,LIS_RCI_PROTOCOL,&e);full_request(c,s,&e,q);return s;
 }
 lis_rci_request_ex e={0};
 s=advance_impl(c,token,(int)result.error,scalar,&e);
 full_request(c,s,&e,q);return s;
}
/* Compatibility entry points complete scalar reductions with the local backend. */
static int local_continue(lis_rci_context *c,int status,lis_rci_request_ex *e,lis_rci_request *q)
{
 while(status==LIS_RCI_DOT||status==LIS_RCI_NORM){
  long double s=status==LIS_RCI_DOT?dot(e->in,e->in2,e->n):norm(e->in,e->n);
  status=advance_impl(c,e->token,0,s,e);
 }
 *q=(lis_rci_request){e->in,e->out,e->n};return status;
}
int lis_rci_start(lis_rci_context *c,const double *b,const double *x0,lis_rci_request *q)
{
 if(!q)return LIS_RCI_BADARG;
 lis_rci_request_ex e={0};int s=start_impl(c,b,x0,0,&e);
 if(s==LIS_RCI_BADARG||s==LIS_RCI_STATE)return s;
 return local_continue(c,s,&e,q);
}
int lis_rci_advance(lis_rci_context *c,int code,lis_rci_request *q)
{
 if(!c||!q)return LIS_RCI_BADARG;
 if(c->external||c->phase==IDLE||c->phase==END)return LIS_RCI_STATE;
 lis_rci_request_ex e={0};int s=advance_impl(c,c->token,code,0,&e);
 return local_continue(c,s,&e,q);
}
/* Transport-independent five-word maximum-reduction protocol. */
int lis_rci_consensus_pack(uint64_t token,uint64_t tag,uint32_t rank,int32_t error,uint64_t words[5])
{
 if(!words)return LIS_RCI_BADARG;
 words[0]=token;words[1]=~token;words[2]=tag;words[3]=~tag;
 words[4]=error?((uint64_t)(UINT32_MAX-rank)<<32)|(uint32_t)error:0;
 return LIS_RCI_DONE;
}
int lis_rci_consensus_merge(uint64_t dst[5],const uint64_t src[5])
{
 if(!dst||!src)return LIS_RCI_BADARG;
 for(size_t i=0;i<5;i++)if(src[i]>dst[i])dst[i]=src[i];
 return LIS_RCI_DONE;
}
int lis_rci_consensus_unpack(const uint64_t words[5],lis_rci_consensus_result *result)
{
 if(!words||!result||~words[1]>words[0]||~words[3]>words[2])return LIS_RCI_BADARG;
 uint32_t code=(uint32_t)words[4];
 if(words[4]&&!code)return LIS_RCI_BADARG;
 *result=(lis_rci_consensus_result){
  words[0]!=~words[1],words[2]!=~words[3],
  words[4]?UINT32_MAX-(uint32_t)(words[4]>>32):UINT32_MAX,
  (int32_t)(code<=INT32_MAX?(int64_t)code:(int64_t)code-INT64_C(4294967296))
 };
 return LIS_RCI_DONE;
}

/* Internal qualification-only size accessor for the nested adapter. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((visibility("hidden")))
#endif
size_t lis_rci_restore_dimension(const lis_rci_context *c){return c?c->o.n:0;}
