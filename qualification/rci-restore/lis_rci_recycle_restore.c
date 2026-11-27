/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_recycle_restore.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
struct lis_rci_recycle_restore {
 size_t n,k,stride,committed,requested,j,rank;
 int flexible,active,stage,code,status;
 double tol,*mem,*scratch;
 uint64_t token,a,p,next_a,next_p;
};
static int finite_vec(const double *v,size_t n)
{if(!v)return 0;for(size_t i=0;i<n;i++)if(!isfinite(v[i]))return 0;return 1;}
static int store(double *out,long double v)
{if(!isfinite(v)||fabsl(v)>DBL_MAX)return 0;*out=(double)v;return 1;}
static long double dot(const double *a,const double *b,size_t n)
{long double v=0;for(size_t i=0;i<n;i++)v+=(long double)a[i]*b[i];return v;}
static long double norm(const double *a,size_t n)
{long double v=0;for(size_t i=0;i<n;i++)v=hypotl(v,a[i]);return v;}
static double *base(lis_rci_recycle_restore *s,int bank,size_t which)
{return s->mem+((size_t)bank*3+which)*s->stride;}
static double *staged(lis_rci_recycle_restore *s,size_t which,size_t j)
{return base(s,1-s->active,which)+j*s->n;}
static int fail(lis_rci_recycle_restore *s,int status,lis_rci_recycle_request *q)
{s->stage=0;s->status=status;if(q)memset(q,0,sizeof(*q));return status;}
static int issue(lis_rci_recycle_restore *s,uint32_t op,lis_rci_recycle_request *q)
{
 if(s->token==UINT64_MAX)return fail(s,LIS_RCI_PROTOCOL,q);
 double *in=staged(s,op==LIS_RCI_OP_PRECON?0:1,s->j);
 double *out=staged(s,op==LIS_RCI_OP_PRECON?1:2,s->j);
 for(size_t i=0;i<s->n;i++)out[i]=NAN;
 *q=(lis_rci_recycle_request){++s->token,s->next_a,s->next_p,op,s->n,in,out};
 s->stage=(int)op;s->status=(int)op;return s->status;
}
lis_rci_recycle_restore *lis_rci_recycle_restore_create(size_t n,size_t k,int flexible,double tol)
{
 if(!n||!k||k>n||n>SIZE_MAX/k||n*k>SIZE_MAX/6/sizeof(double)||
    n>SIZE_MAX/3/sizeof(double)||(flexible!=0&&flexible!=1)||!isfinite(tol)||tol<=0||tol>=1)return NULL;
 lis_rci_recycle_restore *s=calloc(1,sizeof(*s));if(!s)return NULL;
 s->n=n;s->k=k;s->stride=n*k;s->flexible=flexible;s->tol=tol;s->status=LIS_RCI_STATE;
 s->mem=malloc(6*s->stride*sizeof(double));s->scratch=malloc(3*n*sizeof(double));
 if(!s->mem||!s->scratch){lis_rci_recycle_restore_destroy(s);return NULL;}
 return s;
}
void lis_rci_recycle_restore_destroy(lis_rci_recycle_restore *s)
{if(s){free(s->mem);free(s->scratch);free(s);}}
int lis_rci_recycle_restore_begin(lis_rci_recycle_restore *s,size_t k,const double *w,const double *z,
 uint64_t a,uint64_t p,lis_rci_recycle_request *q)
{
 if(!s||!q||!k||k>s->k||!w||(s->flexible?!z:z!=NULL))return LIS_RCI_BADARG;
 if(s->stage)return LIS_RCI_STATE;
 if(!finite_vec(w,k*s->n)||(s->flexible&&!finite_vec(z,k*s->n)))return LIS_RCI_NONFINITE;
 s->requested=k;s->j=s->rank=0;s->next_a=a;s->next_p=p;s->code=0;
 memmove(staged(s,0,0),w,k*s->n*sizeof(double));
 if(s->flexible)memmove(staged(s,1,0),z,k*s->n*sizeof(double));
 return issue(s,s->flexible?LIS_RCI_OP_A:LIS_RCI_OP_PRECON,q);
}
int lis_rci_recycle_restore_refresh(lis_rci_recycle_restore *s,uint64_t a,uint64_t p,lis_rci_recycle_request *q)
{
 if(!s||!q)return LIS_RCI_BADARG;
 if(s->stage||!s->committed)return LIS_RCI_STATE;
 if(s->flexible&&s->p!=p)return LIS_RCI_PROTOCOL;
 return lis_rci_recycle_restore_begin(s,s->committed,base(s,s->active,0),
  s->flexible?base(s,s->active,1):NULL,a,p,q);
}
static int process_column(lis_rci_recycle_restore *s)
{
 size_t n=s->n,j=s->j;
 double *v[3]={staged(s,0,j),staged(s,1,j),staged(s,2,j)};
 long double orig=norm(v[2],n);if(!isfinite(orig))return LIS_RCI_NONFINITE;
 if(orig==0)return 0;
 for(int pass=0;pass<2;pass++)for(size_t k=0;k<s->rank;k++){
  long double a=dot(staged(s,2,k),v[2],n);if(!isfinite(a))return LIS_RCI_NONFINITE;
  for(size_t b=0;b<3;b++)for(size_t i=0;i<n;i++)
   if(!store(v[b]+i,(long double)v[b][i]-a*staged(s,b,k)[i]))return LIS_RCI_NONFINITE;
 }
 long double d=norm(v[2],n);if(!isfinite(d))return LIS_RCI_NONFINITE;
 if(d<=(long double)s->tol*orig)return 0;
 for(size_t b=0;b<3;b++){
  double *dst=staged(s,b,s->rank);
  for(size_t i=0;i<n;i++)if(!store(dst+i,(long double)v[b][i]/d))return LIS_RCI_NONFINITE;
 }
 s->rank++;return 0;
}
int lis_rci_recycle_restore_advance(lis_rci_recycle_restore *s,uint64_t token,uint64_t a,uint64_t p,int code,lis_rci_recycle_request *q)
{
 if(!s||!q)return LIS_RCI_BADARG;
 if(!s->stage)return LIS_RCI_STATE;
 if(token!=s->token)return LIS_RCI_PROTOCOL;
 if(a!=s->next_a||p!=s->next_p)return fail(s,LIS_RCI_PROTOCOL,q);
 if(code){s->code=code;return fail(s,LIS_RCI_CALLBACK,q);}
 if(!finite_vec(staged(s,s->stage==LIS_RCI_OP_PRECON?1:2,s->j),s->n))return fail(s,LIS_RCI_NONFINITE,q);
 if(s->stage==LIS_RCI_OP_PRECON)return issue(s,LIS_RCI_OP_A,q);
 int e=process_column(s);if(e)return fail(s,e,q);
 if(++s->j<s->requested)return issue(s,s->flexible?LIS_RCI_OP_A:LIS_RCI_OP_PRECON,q);
 if(!s->rank)return fail(s,LIS_RCI_BREAKDOWN,q);
 s->active=1-s->active;s->committed=s->rank;s->a=a;s->p=p;
 return fail(s,0,q);
}
int lis_rci_recycle_restore_cancel(lis_rci_recycle_restore *s,uint64_t token)
{
 if(!s)return LIS_RCI_BADARG;
 if(!s->stage)return LIS_RCI_STATE;
 if(token!=s->token)return LIS_RCI_PROTOCOL;
 return fail(s,LIS_RCI_CANCELLED,NULL);
}
int lis_rci_recycle_restore_error(const lis_rci_recycle_restore *s,int *code)
{if(!s||!code)return LIS_RCI_BADARG;*code=s->code;return s->status;}
int lis_rci_recycle_restore_get(const lis_rci_recycle_restore *s,lis_rci_recycle_basis *b)
{
 if(!s||!b)return LIS_RCI_BADARG;
 memset(b,0,sizeof(*b));if(!s->committed)return LIS_RCI_STATE;
 const double *w=s->mem+(size_t)s->active*3*s->stride;
 *b=(lis_rci_recycle_basis){s->n,s->committed,s->a,s->p,w,w+s->stride,w+2*s->stride};
 return 0;
}
static int ready(lis_rci_recycle_restore *s,uint64_t a,uint64_t p)
{
 if(!s||!s->committed||s->stage)return LIS_RCI_STATE;
 return s->a==a&&s->p==p?0:LIS_RCI_PROTOCOL;
}
int lis_rci_recycle_restore_project(lis_rci_recycle_restore *s,uint64_t a,uint64_t p,double *x,double *r)
{
 int e=ready(s,a,p);if(e)return e;
 size_t n=s->n;if(x==r||!finite_vec(x,n)||!finite_vec(r,n))return LIS_RCI_BADARG;
 double *xx=s->scratch,*rr=xx+n;
 memcpy(xx,x,n*sizeof(double));memcpy(rr,r,n*sizeof(double));
 for(size_t j=0;j<s->committed;j++){
  const double *z=base(s,s->active,1)+j*n,*c=base(s,s->active,2)+j*n;
  long double t=dot(c,rr,n);if(!isfinite(t))return LIS_RCI_NONFINITE;
  for(size_t i=0;i<n;i++)if(!store(xx+i,(long double)xx[i]+t*z[i])||
     !store(rr+i,(long double)rr[i]-t*c[i]))return LIS_RCI_NONFINITE;
 }
 memcpy(x,xx,n*sizeof(double));memcpy(r,rr,n*sizeof(double));return 0;
}
int lis_rci_recycle_restore_direction(lis_rci_recycle_restore *s,uint64_t a,uint64_t p,double *w,double *z,double *az)
{
 int e=ready(s,a,p);if(e)return e;
 size_t n=s->n;if(w==z||w==az||z==az||!finite_vec(w,n)||!finite_vec(z,n)||!finite_vec(az,n))return LIS_RCI_BADARG;
 double *arr[3]={w,z,az},*dst[3]={s->scratch,s->scratch+n,s->scratch+2*n};
 for(size_t k=0;k<3;k++)memcpy(dst[k],arr[k],n*sizeof(double));
 for(size_t j=0;j<s->committed;j++){
  long double t=dot(base(s,s->active,2)+j*n,dst[2],n);
  if(!isfinite(t))return LIS_RCI_NONFINITE;
  for(size_t k=0;k<3;k++)for(size_t i=0;i<n;i++)
   if(!store(dst[k]+i,(long double)dst[k][i]-t*base(s,s->active,k)[j*n+i]))return LIS_RCI_NONFINITE;
 }
 for(size_t k=0;k<3;k++)memcpy(arr[k],dst[k],n*sizeof(double));
 return 0;
}

int lis_rci_recycle_restore_commit_images(lis_rci_recycle_restore *s,size_t k,const double *w,const double *z,const double *az,uint64_t a,uint64_t p)
{
 if(!s||!k||k>s->k||!w||!z||!az)return LIS_RCI_BADARG;
 if(s->stage)return LIS_RCI_STATE;
 size_t size=k*s->n;
 if(!finite_vec(w,size)||!finite_vec(z,size)||!finite_vec(az,size))return LIS_RCI_NONFINITE;
 memmove(staged(s,0,0),w,size*sizeof(double));
 memmove(staged(s,1,0),z,size*sizeof(double));
 memmove(staged(s,2,0),az,size*sizeof(double));
 s->rank=0;s->code=0;
 for(s->j=0;s->j<k;s->j++){
  int e=process_column(s);if(e)return fail(s,e,NULL);
 }
 if(!s->rank)return fail(s,LIS_RCI_BREAKDOWN,NULL);
 s->active=1-s->active;s->committed=s->rank;s->a=a;s->p=p;
 return fail(s,0,NULL);
}
