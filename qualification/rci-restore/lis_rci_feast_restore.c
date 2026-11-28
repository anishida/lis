/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_feast_restore.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
enum phase { F_IDLE,F_PROJECT,F_BQ,F_AQ,F_AX,F_BX,F_DONE,F_FAILED };
enum groups { Q,R,X,Y0,Y1,BQ,AQ,RITZ,AX,BX,NGROUPS };
struct lis_rci_feast_restore {
 lis_rci_feast_restore_options o;
 lis_rci_projector_restore *p;
 lis_rci_split_vector *v;
 lis_rci_handle *h;
 double *mem,*dense,*coeff,*eval,*res,*out,*out_eval,*out_res;
 size_t rank,cycles,count,previous_count;
 long double previous_trace;
 uint64_t a,b,token,inner_token;
 int code,status,inner_active;
 enum phase phase;
};
static int mul(size_t a,size_t b,size_t *out)
{if(b&&a>SIZE_MAX/b)return 0;*out=a*b;return 1;}
static int store(double *out,long double x)
{if(!isfinite(x)||fabsl(x)>DBL_MAX)return 0;*out=(double)x;return 1;}
static lis_rci_split_vector *vec(lis_rci_feast_restore *s,size_t g,size_t j)
{return s->v+g*s->o.columns+j;}
static lis_rci_handle *handles(lis_rci_feast_restore *s,size_t g)
{return s->h+g*s->o.columns;}
static int end(lis_rci_feast_restore *s,int status,lis_rci_opaque_request *q)
{
 if(s->inner_active){(void)lis_rci_projector_restore_cancel(s->p,s->inner_token);s->inner_active=0;}
 s->phase=status==0?F_DONE:F_FAILED;s->status=status;
 if(q)memset(q,0,sizeof(*q));
 return status;
}
static int publish(lis_rci_feast_restore *s,enum phase phase,lis_rci_opaque_request *v,lis_rci_opaque_request *q)
{
 if(s->token==UINT64_MAX)return end(s,LIS_RCI_PROTOCOL,q);
 v->token=++s->token;s->phase=phase;s->status=(int)v->operation;*q=*v;return s->status;
}
static int delegate(lis_rci_feast_restore *s,int status,lis_rci_opaque_request *v,lis_rci_opaque_request *q)
{
 if(status<=0)return end(s,status?status:LIS_RCI_STATE,q);
 s->inner_token=v->token;s->inner_active=1;
 return publish(s,F_PROJECT,v,q);
}
static int issue(lis_rci_feast_restore *s,enum phase phase,uint32_t op,size_t from,size_t to,size_t columns,lis_rci_opaque_request *q)
{
 lis_rci_opaque_request v={0};v.version=LIS_RCI_OPAQUE_VERSION;v.size=sizeof v;v.operation=op;
 v.generation_a=s->a;v.generation_b=s->b;v.columns=columns;v.in=handles(s,from);v.out=handles(s,to);
 return publish(s,phase,&v,q);
}
static int begin_projection(lis_rci_feast_restore *s,lis_rci_opaque_request *q)
{
 lis_rci_opaque_request v;int status=lis_rci_projector_restore_start(s->p,s->a,s->b,&v);
 return delegate(s,status,&v,q);
}
lis_rci_feast_restore *lis_rci_feast_restore_create(const lis_rci_feast_restore_options *o)
{
 size_t mn=0,mm=0,vn=0,plane=0;
 if(!o||!o->n||!o->columns||o->columns>o->n||!o->nodes||!o->max_cycles||
    !o->projected_solve||!o->shift||!o->weight||(o->generalized!=0&&o->generalized!=1)||
    !isfinite(o->lo)||!isfinite(o->hi)||o->lo>=o->hi||
    !isfinite(o->residual_tol)||o->residual_tol<=0||
    !isfinite(o->trace_tol)||o->trace_tol<=0||
    !isfinite(o->rank_tol)||o->rank_tol<=0||o->rank_tol>=1||
    !isfinite(o->imag_tol)||o->imag_tol<=0||
    !mul(o->n,o->columns,&mn)||!mul(o->columns,o->columns,&mm)||
    !mul(NGROUPS,o->columns,&vn)||!mul(2*NGROUPS,mn,&plane)||
    plane>SIZE_MAX/sizeof(double)||mn>SIZE_MAX/sizeof(double)||
    mm>SIZE_MAX/sizeof(double)||vn>SIZE_MAX/sizeof(lis_rci_split_vector)||
    vn>SIZE_MAX/sizeof(lis_rci_handle)||o->columns>SIZE_MAX/sizeof(double))
  return NULL;
 lis_rci_feast_restore *s=calloc(1,sizeof(*s));if(!s)return NULL;
 s->o=*o;s->status=LIS_RCI_STATE;
 s->v=calloc(vn,sizeof(*s->v));s->h=calloc(vn,sizeof(*s->h));s->mem=calloc(plane,sizeof(double));
 s->dense=malloc(mm*sizeof(double));s->coeff=malloc(mm*sizeof(double));
 s->eval=malloc(o->columns*sizeof(double));s->res=malloc(o->columns*sizeof(double));
 s->out=malloc(mn*sizeof(double));s->out_eval=malloc(o->columns*sizeof(double));
 s->out_res=malloc(o->columns*sizeof(double));
 if(!s->v||!s->h||!s->mem||!s->dense||!s->coeff||!s->eval||!s->res||
    !s->out||!s->out_eval||!s->out_res){lis_rci_feast_restore_destroy(s);return NULL;}
 for(size_t i=0;i<vn;i++){s->v[i]=(lis_rci_split_vector){o->n,s->mem+2*i*o->n,s->mem+(2*i+1)*o->n};s->h[i]=s->v+i;}
 lis_rci_projector_restore_options p={o->columns,o->nodes,o->generalized,o->shift,o->weight,
  handles(s,Q),handles(s,R),handles(s,X),handles(s,Y0),handles(s,Y1)};
 s->p=lis_rci_projector_restore_create(&p);
 if(!s->p){lis_rci_feast_restore_destroy(s);return NULL;}
 /* The projector has copied the quadrature; do not retain borrowed pointers. */
 s->o.shift=s->o.weight=NULL;
 return s;
}
void lis_rci_feast_restore_destroy(lis_rci_feast_restore *s)
{
 if(!s)return;
 lis_rci_projector_restore_destroy(s->p);free(s->v);free(s->h);free(s->mem);
 free(s->dense);free(s->coeff);free(s->eval);free(s->res);free(s->out);free(s->out_eval);free(s->out_res);free(s);
}
int lis_rci_feast_restore_start(lis_rci_feast_restore *s,const double *initial,uint64_t a,uint64_t b,lis_rci_opaque_request *q)
{
 if(!s||!initial||!q)return LIS_RCI_BADARG;
 if(s->phase!=F_IDLE&&s->phase!=F_DONE&&s->phase!=F_FAILED)return LIS_RCI_STATE;
 s->a=a;s->b=b;s->code=0;s->cycles=s->count=s->rank=0;s->previous_count=SIZE_MAX;s->previous_trace=0;
 size_t n=s->o.n;
 for(size_t j=0;j<s->o.columns;j++){
  lis_rci_split_vector *v=vec(s,Q,j);
  for(size_t i=0;i<n;i++)if(!isfinite(initial[j*n+i]))return end(s,LIS_RCI_NONFINITE,q);
  memcpy(v->real,initial+j*n,n*sizeof(double));memset(v->imag,0,n*sizeof(double));
 }
 return begin_projection(s,q);
}
static int real_vector(const lis_rci_split_vector *v,double tol)
{
 double re=0,im=0;
 for(size_t i=0;i<v->n;i++){
  if(!isfinite(v->real[i])||!isfinite(v->imag[i]))return 0;
  re=fmax(re,fabs(v->real[i]));im=fmax(im,fabs(v->imag[i]));
 }
 return (long double)im<=(long double)tol*(1+(long double)re);
}
static long double dot(const double *a,const double *b,size_t n)
{long double r=0;for(size_t i=0;i<n;i++)r+=(long double)a[i]*b[i];return r;}
static int b_orthonormalize(lis_rci_feast_restore *s)
{
 size_t n=s->o.n,m=s->o.columns,rank=0;
 for(size_t j=0;j<m;j++){
  lis_rci_split_vector *v=vec(s,Q,j),*bv=vec(s,BQ,j);
  if(!real_vector(bv,s->o.imag_tol))return LIS_RCI_NONFINITE;
  long double original=dot(v->real,bv->real,n);
  if(!isfinite(original)||original<0)return LIS_RCI_BREAKDOWN;
  if(original==0)continue;
  for(int pass=0;pass<2;pass++)for(size_t k=0;k<rank;k++){
   lis_rci_split_vector *u=vec(s,Q,k),*bu=vec(s,BQ,k);
   long double t=dot(u->real,bv->real,n);
   if(!isfinite(t))return LIS_RCI_NONFINITE;
   for(size_t i=0;i<n;i++){
    if(!store(v->real+i,(long double)v->real[i]-t*u->real[i])||
       !store(bv->real+i,(long double)bv->real[i]-t*bu->real[i]))return LIS_RCI_NONFINITE;
   }
  }
  long double rr=dot(v->real,bv->real,n);
  if(!isfinite(rr))return LIS_RCI_NONFINITE;
  long double threshold=(long double)s->o.rank_tol*s->o.rank_tol*original;
  if(rr<=threshold)continue;
  long double d=sqrtl(rr);
  lis_rci_split_vector *u=vec(s,Q,rank),*bu=vec(s,BQ,rank);
  for(size_t i=0;i<n;i++){
   if(!store(u->real+i,(long double)v->real[i]/d)||!store(bu->real+i,(long double)bv->real[i]/d))return LIS_RCI_NONFINITE;
   u->imag[i]=bu->imag[i]=0;
  }
  rank++;
 }
 s->rank=rank;
 return rank?0:LIS_RCI_BREAKDOWN;
}
static int rayleigh_ritz(lis_rci_feast_restore *s)
{
 size_t n=s->o.n,r=s->rank;
 for(size_t j=0;j<r;j++)if(!real_vector(vec(s,AQ,j),s->o.imag_tol))return LIS_RCI_NONFINITE;
 for(size_t i=0;i<r;i++)for(size_t j=0;j<r;j++){
  long double t=dot(vec(s,Q,i)->real,vec(s,AQ,j)->real,n);
  if(!store(s->dense+i*r+j,t))return LIS_RCI_NONFINITE;
 }
 for(size_t i=0;i<r;i++)for(size_t j=i+1;j<r;j++){
  long double a=s->dense[i*r+j],b=s->dense[j*r+i];
  if(fabsl(a-b)>1e-9L*(1+fabsl(a)+fabsl(b)))return LIS_RCI_BREAKDOWN;
  s->dense[i*r+j]=s->dense[j*r+i]=(double)(.5L*a+.5L*b);
 }
 for(size_t i=0;i<r;i++)s->eval[i]=NAN;
 for(size_t i=0;i<r*r;i++)s->coeff[i]=NAN;
 int code=s->o.projected_solve(s->o.user,r,s->dense,s->eval,s->coeff);
 if(code){s->code=code;return LIS_RCI_CALLBACK;}
 for(size_t i=0;i<r;i++)if(!isfinite(s->eval[i]))return LIS_RCI_NONFINITE;
 for(size_t i=0;i<r*r;i++)if(!isfinite(s->coeff[i]))return LIS_RCI_NONFINITE;
 for(size_t j=0;j<r;j++){
  lis_rci_split_vector *v=vec(s,RITZ,j);
  long double norm=0;
  for(size_t i=0;i<n;i++){
   long double a=0;for(size_t k=0;k<r;k++)a+=(long double)vec(s,Q,k)->real[i]*s->coeff[k*r+j];
   if(!store(v->real+i,a))return LIS_RCI_NONFINITE;
   norm=hypotl(norm,a);v->imag[i]=0;
  }
  if(norm==0||!isfinite(norm))return LIS_RCI_BREAKDOWN;
 }
 return 0;
}
static int verify_ritz(lis_rci_feast_restore *s,lis_rci_opaque_request *q)
{
 size_t n=s->o.n,r=s->rank,count=0;long double trace=0;
 int converged=1;
 for(size_t j=0;j<r;j++){
  lis_rci_split_vector *a=vec(s,AX,j),*b=vec(s,BX,j);
  if(!real_vector(a,s->o.imag_tol)||!real_vector(b,s->o.imag_tol))return end(s,LIS_RCI_NONFINITE,q);
  long double nr=0,na=0,nb=0,lambda=s->eval[j];
  for(size_t i=0;i<n;i++){
   nr=hypotl(nr,(long double)a->real[i]-lambda*b->real[i]);
   na=hypotl(na,a->real[i]);nb=hypotl(nb,b->real[i]);
  }
  long double den=na+fabsl(lambda)*nb;
  if(!isfinite(nr)||!isfinite(den)||!store(s->res+j,den>0?nr/den:nr))
   return end(s,LIS_RCI_NONFINITE,q);
  if(lambda>=s->o.lo&&lambda<=s->o.hi){
   if(s->res[j]>s->o.residual_tol)converged=0;
   trace+=lambda;count++;
  }
 }
 s->cycles++;s->count=count;
 if(count&&converged&&count==s->previous_count&&
    isfinite(trace)&&fabsl(trace-s->previous_trace)<=(long double)s->o.trace_tol*(1+fabsl(trace))){
  size_t k=0;
  for(size_t j=0;j<r;j++)if(s->eval[j]>=s->o.lo&&s->eval[j]<=s->o.hi){
   s->out_eval[k]=s->eval[j];s->out_res[k]=s->res[j];
   memcpy(s->out+k*n,vec(s,RITZ,j)->real,n*sizeof(double));k++;
  }
  for(size_t i=0;i<count;i++)for(size_t j=i+1;j<count;j++)if(s->out_eval[j]<s->out_eval[i]){
   double t=s->out_eval[i];s->out_eval[i]=s->out_eval[j];s->out_eval[j]=t;
   t=s->out_res[i];s->out_res[i]=s->out_res[j];s->out_res[j]=t;
   for(size_t k2=0;k2<n;k2++){t=s->out[i*n+k2];s->out[i*n+k2]=s->out[j*n+k2];s->out[j*n+k2]=t;}
  }
  return end(s,0,q);
 }
 s->previous_count=count;s->previous_trace=trace;
 if(s->cycles>=s->o.max_cycles)return end(s,LIS_RCI_MAXITER,q);
 for(size_t j=0;j<s->o.columns;j++){
  lis_rci_split_vector *v=vec(s,Q,j);
  if(j<r)memcpy(v->real,vec(s,RITZ,j)->real,n*sizeof(double));else memset(v->real,0,n*sizeof(double));
  memset(v->imag,0,n*sizeof(double));
 }
 return begin_projection(s,q);
}
int lis_rci_feast_restore_advance(lis_rci_feast_restore *s,uint64_t token,uint64_t a,uint64_t b,int code,lis_rci_opaque_request *q)
{
 if(!s||!q)return LIS_RCI_BADARG;
 if(s->phase==F_IDLE||s->phase==F_DONE||s->phase==F_FAILED)return LIS_RCI_STATE;
 if(token!=s->token)return LIS_RCI_PROTOCOL;
 if(a!=s->a||b!=s->b)return end(s,LIS_RCI_PROTOCOL,q);
 if(code){s->code=code;return end(s,LIS_RCI_CALLBACK,q);}
 int e=0;
 switch(s->phase){
 case F_PROJECT:{
  lis_rci_opaque_request v;
  e=lis_rci_projector_restore_advance(s->p,s->inner_token,a,b,0,&v);
  if(e>0)return delegate(s,e,&v,q);
  s->inner_active=0;if(e)return end(s,e,q);
  const lis_rci_handle *res=lis_rci_projector_restore_result(s->p);
  if(!res)return end(s,LIS_RCI_STATE,q);
  for(size_t j=0;j<s->o.columns;j++){
   lis_rci_split_vector *src=res[j],*dst=vec(s,Q,j);
   if(!real_vector(src,s->o.imag_tol))return end(s,LIS_RCI_NONFINITE,q);
   memcpy(dst->real,src->real,s->o.n*sizeof(double));memset(dst->imag,0,s->o.n*sizeof(double));
  }
  return issue(s,F_BQ,LIS_RCI_OP_B,Q,BQ,s->o.columns,q);
 }
 case F_BQ:
  e=b_orthonormalize(s);if(e)return end(s,e,q);
  return issue(s,F_AQ,LIS_RCI_OP_A,Q,AQ,s->rank,q);
 case F_AQ:
  e=rayleigh_ritz(s);if(e)return end(s,e,q);
  return issue(s,F_AX,LIS_RCI_OP_A,RITZ,AX,s->rank,q);
 case F_AX:return issue(s,F_BX,LIS_RCI_OP_B,RITZ,BX,s->rank,q);
 case F_BX:return verify_ritz(s,q);
 default:return end(s,LIS_RCI_STATE,q);
 }
}
int lis_rci_feast_restore_cancel(lis_rci_feast_restore *s,uint64_t token)
{
 if(!s)return LIS_RCI_BADARG;
 if(s->phase==F_IDLE||s->phase==F_DONE||s->phase==F_FAILED)return LIS_RCI_STATE;
 if(token!=s->token)return LIS_RCI_PROTOCOL;
 return end(s,LIS_RCI_CANCELLED,NULL);
}
int lis_rci_feast_restore_error(const lis_rci_feast_restore *s,int *code)
{if(!s||!code)return LIS_RCI_BADARG;*code=s->code;return s->status;}
int lis_rci_feast_restore_get_result(const lis_rci_feast_restore *s,lis_rci_feast_restore_result *r)
{
 if(!s||!r)return LIS_RCI_BADARG;
 memset(r,0,sizeof(*r));
 if(s->phase!=F_DONE)return LIS_RCI_STATE;
 *r=(lis_rci_feast_restore_result){s->count,s->cycles,s->rank,s->out_eval,s->out,s->out_res};return 0;
}
