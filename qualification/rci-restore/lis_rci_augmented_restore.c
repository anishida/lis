/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_augmented_restore.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
enum aug_phase { A_IDLE,A_PROJECT,A_CHILD_A,A_CHILD_P,A_DIRECTION_A,A_DONE };
struct lis_rci_augmented_restore {
 lis_rci_augmented_restore_options o;
 lis_rci_augmented_restore_info info;
 lis_rci_context *child;
 lis_rci_recycle_restore *recycle;
 lis_rci_request child_q;
 lis_rci_recycle_basis old;
 enum aug_phase phase;
 size_t m,p,used,oldk;
 uint64_t token,a,pc;
 int child_active;
 double *mem,*b,*x,*r,*tmp,*w,*z,*d,*neww,*newz,*newd,*gram,*overlap,*coef;
};
static int mul(size_t a,size_t b,size_t *z)
{if(b&&a>SIZE_MAX/b)return 0;*z=a*b;return 1;}
static int finite_vec(const double *v,size_t n)
{if(!v)return 0;for(size_t i=0;i<n;i++)if(!isfinite(v[i]))return 0;return 1;}
static int put(double *out,long double x)
{if(!isfinite(x)||fabsl(x)>DBL_MAX)return 0;*out=(double)x;return 1;}
static long double dot(const double *x,const double *y,size_t n)
{long double d=0;for(size_t i=0;i<n;i++)d+=(long double)x[i]*y[i];return d;}
static void copy(double *d,const double *s,size_t n){memmove(d,s,n*sizeof(double));}
static int finish(lis_rci_augmented_restore *s,int status,lis_rci_recycle_request *q)
{
 if(s->child_active){
  lis_rci_request_full pending;
  if(lis_rci_get_pending(s->child,&pending)>0)(void)lis_rci_cancel(s->child,pending.request.token);
  s->child_active=0;
 }
 s->phase=A_DONE;s->info.status=status;
 lis_rci_recycle_basis b;
 s->info.recycle_rank=lis_rci_recycle_restore_get(s->recycle,&b)?0:b.count;
 if(q)memset(q,0,sizeof(*q));
 return status;
}
static int issue(lis_rci_augmented_restore *s,enum aug_phase phase,uint32_t op,const double *in,double *out,lis_rci_recycle_request *q)
{
 if(s->token==UINT64_MAX)return finish(s,LIS_RCI_PROTOCOL,q);
 if(!finite_vec(in,s->o.n))return finish(s,LIS_RCI_NONFINITE,q);
 for(size_t i=0;i<s->o.n;i++)out[i]=NAN;
 s->phase=phase;s->info.status=(int)op;
 *q=(lis_rci_recycle_request){++s->token,s->a,s->pc,op,s->o.n,in,out};
 if(op==LIS_RCI_OP_A)s->info.matvecs++;else s->info.precons++;
 return (int)op;
}
lis_rci_augmented_restore *lis_rci_augmented_restore_create(const lis_rci_augmented_restore_options *o)
{
 size_t n,m,k,p,large,dense,total;
 if(!o||!o->n||!o->restart||!o->recycle||o->recycle>=o->restart||o->restart>o->n||
    o->restart>512||!o->max_cycles||o->max_cycles>SIZE_MAX/o->restart||
    (o->flexible!=0&&o->flexible!=1)||!o->select)return NULL;
 n=o->n;m=o->restart;k=o->recycle;p=m+k;
 if(!mul(n,4+3*p+3*k,&large)||!mul(p,p,&dense)||dense>SIZE_MAX/3||
    large>SIZE_MAX-3*dense)return NULL;
 total=large+3*dense;if(total>SIZE_MAX/sizeof(double))return NULL;
 lis_rci_augmented_restore *s=calloc(1,sizeof(*s));if(!s)return NULL;
 s->o=*o;s->m=m;s->p=p;s->info.status=LIS_RCI_STATE;s->info.residual=NAN;
 lis_rci_options child={n,m,m,o->rtol,o->atol,o->flexible?LIS_RCI_FGMRES:LIS_RCI_GMRES};
 s->child=lis_rci_create(&child);s->recycle=lis_rci_recycle_restore_create(n,k,o->flexible,o->rank_tol);
 s->mem=calloc(total,sizeof(double));
 if(!s->child||!s->recycle||!s->mem){lis_rci_augmented_restore_destroy(s);return NULL;}
 s->b=s->mem;s->x=s->b+n;s->r=s->x+n;s->tmp=s->r+n;
 s->w=s->tmp+n;s->z=s->w+p*n;s->d=s->z+p*n;
 s->neww=s->d+p*n;s->newz=s->neww+k*n;s->newd=s->newz+k*n;
 s->gram=s->newd+k*n;s->overlap=s->gram+dense;s->coef=s->overlap+dense;
 return s;
}
void lis_rci_augmented_restore_destroy(lis_rci_augmented_restore *s)
{if(s){lis_rci_destroy(s->child);lis_rci_recycle_restore_destroy(s->recycle);free(s->mem);free(s);}}
/* Selection uses the augmented physical Arnoldi relation, without reapplying
 * either A or M. The same coefficient matrix transforms all three bases. */
static void update_recycle(lis_rci_augmented_restore *s)
{
 size_t n=s->o.n,p=s->oldk+s->used,k=s->o.recycle;
 if(!s->used||!p)return;
 if(k>p)k=p;
 const double *f=s->o.flexible?s->z:s->w;
 for(size_t j=0;j<p;j++)for(size_t i=0;i<p;i++){
  if(!put(s->gram+i+j*p,dot(s->d+i*n,s->d+j*n,n))||
     !put(s->overlap+i+j*p,dot(s->d+i*n,f+j*n,n))){s->info.rejects++;return;}
 }
 for(size_t i=0;i<p*k;i++)s->coef[i]=NAN;
 size_t selected=0;
 if(s->o.select(s->o.user,p,k,s->gram,s->overlap,s->coef,&selected)||
    !selected||selected>k||!finite_vec(s->coef,p*selected)){s->info.rejects++;return;}
 double *src[3]={s->w,s->z,s->d},*dst[3]={s->neww,s->newz,s->newd};
 for(size_t g=0;g<3;g++)for(size_t j=0;j<selected;j++)for(size_t i=0;i<n;i++){
  long double a=0;
  for(size_t l=0;l<p;l++)a+=(long double)src[g][l*n+i]*s->coef[l+j*p];
  if(!put(dst[g]+j*n+i,a)){s->info.rejects++;return;}
 }
 if(lis_rci_recycle_restore_commit_images(s->recycle,selected,s->neww,s->newz,s->newd,s->a,s->pc)){
  s->info.rejects++;return;
 }
 if(s->oldk)s->info.updates++;else s->info.builds++;
}
static int child_result(lis_rci_augmented_restore *,int,lis_rci_recycle_request *);
static int begin_cycle(lis_rci_augmented_restore *s,lis_rci_recycle_request *q)
{
 s->oldk=s->used=0;
 if(!lis_rci_recycle_restore_get(s->recycle,&s->old)&&
    s->old.operator_generation==s->a&&s->old.precon_generation==s->pc){
  s->oldk=s->old.count;size_t len=s->oldk*s->o.n;
  copy(s->w,s->old.domain,len);copy(s->z,s->old.physical,len);copy(s->d,s->old.image,len);
  s->info.recycle_uses++;
  return issue(s,A_PROJECT,LIS_RCI_OP_A,s->x,s->tmp,q);
 }
 int status=lis_rci_start(s->child,s->b,s->x,&s->child_q);
 s->child_active=status>0;return child_result(s,status,q);
}
static int child_result(lis_rci_augmented_restore *s,int status,lis_rci_recycle_request *q)
{
 if(status==LIS_RCI_PRECON){
  if(s->used>=s->m)return finish(s,LIS_RCI_STATE,q);
  size_t off=(s->oldk+s->used)*s->o.n;
  copy(s->w+off,s->child_q.in,s->o.n);
  return issue(s,A_CHILD_P,LIS_RCI_OP_PRECON,s->child_q.in,s->z+off,q);
 }
 if(status==LIS_RCI_MATVEC)
  return issue(s,A_CHILD_A,LIS_RCI_OP_A,s->child_q.in,s->child_q.out,q);
 s->child_active=0;
 lis_rci_info ci;lis_rci_get_info(s->child,&ci);
 s->info.iterations+=ci.iterations;s->info.cycles++;
 lis_rci_get_residual(s->child,&s->info.residual);
 if(status!=LIS_RCI_DONE&&status!=LIS_RCI_MAXITER)return finish(s,status,q);
 copy(s->x,lis_rci_solution(s->child),s->o.n);
 update_recycle(s);
 if(status==LIS_RCI_DONE||s->info.cycles>=s->o.max_cycles)return finish(s,status,q);
 return begin_cycle(s,q);
}
int lis_rci_augmented_restore_start(lis_rci_augmented_restore *s,const double *b,const double *x0,uint64_t a,uint64_t p,lis_rci_recycle_request *q)
{
 if(!s||!b||!q)return LIS_RCI_BADARG;
 if(s->phase!=A_IDLE&&s->phase!=A_DONE)return LIS_RCI_STATE;
 s->info=(lis_rci_augmented_restore_info){.residual=NAN};s->a=a;s->pc=p;
 if(!finite_vec(b,s->o.n)||(x0&&!finite_vec(x0,s->o.n)))return finish(s,LIS_RCI_NONFINITE,q);
 copy(s->b,b,s->o.n);if(x0)copy(s->x,x0,s->o.n);else memset(s->x,0,s->o.n*sizeof(double));
 return begin_cycle(s,q);
}
int lis_rci_augmented_restore_advance(lis_rci_augmented_restore *s,uint64_t token,uint64_t a,uint64_t p,int code,lis_rci_recycle_request *q)
{
 if(!s||!q)return LIS_RCI_BADARG;
 if(s->phase==A_IDLE||s->phase==A_DONE)return LIS_RCI_STATE;
 if(token!=s->token)return LIS_RCI_PROTOCOL;
 if(a!=s->a||p!=s->pc)return finish(s,LIS_RCI_PROTOCOL,q);
 if(code){s->info.callback_code=code;return finish(s,LIS_RCI_CALLBACK,q);}
 size_t n=s->o.n,off=(s->oldk+s->used)*n;
 int status,e;
 switch(s->phase){
 case A_PROJECT:
  if(!finite_vec(s->tmp,n))return finish(s,LIS_RCI_NONFINITE,q);
  for(size_t i=0;i<n;i++)if(!put(s->r+i,(long double)s->b[i]-s->tmp[i]))return finish(s,LIS_RCI_NONFINITE,q);
  e=lis_rci_recycle_restore_project(s->recycle,a,p,s->x,s->r);if(e)return finish(s,e,q);
  status=lis_rci_start(s->child,s->b,s->x,&s->child_q);s->child_active=status>0;
  return child_result(s,status,q);
 case A_CHILD_A:
  status=lis_rci_advance(s->child,0,&s->child_q);return child_result(s,status,q);
 case A_CHILD_P:
  if(!finite_vec(s->z+off,n))return finish(s,LIS_RCI_NONFINITE,q);
  return issue(s,A_DIRECTION_A,LIS_RCI_OP_A,s->z+off,s->d+off,q);
 case A_DIRECTION_A:
  if(!finite_vec(s->d+off,n))return finish(s,LIS_RCI_NONFINITE,q);
  if(s->oldk){
   e=lis_rci_recycle_restore_direction(s->recycle,a,p,s->w+off,s->z+off,s->d+off);
   if(e)return finish(s,e,q);
  }
  copy(s->child_q.out,s->z+off,n);
  status=lis_rci_advance(s->child,0,&s->child_q);
  if(status!=LIS_RCI_MATVEC)return finish(s,status<0?status:LIS_RCI_STATE,q);
  /* GM_P -> GM_A: complete with A*z already computed and projected.
   * The child has not seen or requested any second operator evaluation. */
  if(memcmp(s->child_q.in,s->z+off,n*sizeof(double)))return finish(s,LIS_RCI_STATE,q);
  copy(s->child_q.out,s->d+off,n);s->used++;
  status=lis_rci_advance(s->child,0,&s->child_q);return child_result(s,status,q);
 default:return finish(s,LIS_RCI_STATE,q);
 }
}
int lis_rci_augmented_restore_cancel(lis_rci_augmented_restore *s,uint64_t token)
{
 if(!s)return LIS_RCI_BADARG;
 if(s->phase==A_IDLE||s->phase==A_DONE)return LIS_RCI_STATE;
 if(token!=s->token)return LIS_RCI_PROTOCOL;
 return finish(s,LIS_RCI_CANCELLED,NULL);
}
int lis_rci_augmented_restore_get_info(const lis_rci_augmented_restore *s,lis_rci_augmented_restore_info *i)
{if(!s||!i)return LIS_RCI_BADARG;*i=s->info;return 0;}
const double *lis_rci_augmented_restore_solution(const lis_rci_augmented_restore *s)
{return s&&s->phase==A_DONE&&s->info.status==0?s->x:NULL;}
int lis_rci_augmented_restore_basis(const lis_rci_augmented_restore *s,lis_rci_recycle_basis *b)
{return s?lis_rci_recycle_restore_get(s->recycle,b):LIS_RCI_BADARG;}
