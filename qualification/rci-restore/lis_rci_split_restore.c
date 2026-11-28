/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_restore.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct lis_rci_split_restore {
 size_t n;
 int generalized,status,code;
 void *user;
 lis_rci_real_apply apply;
 lis_rci_generation_read generations;
 lis_rci_context *child;
 double *mem,*rhs,*solution,*ar,*ai,*br,*bi;
 uint64_t a,b;
 lis_rci_complex_coeff shift;
};
static int vector_ok(const lis_rci_split_vector *v,size_t n)
{
 return v&&v->n==n&&v->real&&v->imag&&v->real!=v->imag;
}
static int alias(const lis_rci_split_vector *a,const lis_rci_split_vector *b)
{
 return a->real==b->real||a->real==b->imag||a->imag==b->real||a->imag==b->imag;
}
static int split_gens(void *user,uint64_t *a,uint64_t *b)
{
 lis_rci_split_restore *s=user;
 return s->generations(s->user,a,b);
}
static int check_gen(lis_rci_split_restore *s)
{
 uint64_t a=0,b=0;
 int e=split_gens(s,&a,&b);
 if(e){s->code=e;return LIS_RCI_CALLBACK;}
 return a==s->a&&b==s->b?0:LIS_RCI_PROTOCOL;
}
static int real_apply(lis_rci_split_restore *s,uint32_t op,const double *in,double *out)
{
 int e=check_gen(s);if(e)return e;
 if(op==LIS_RCI_OP_B&&!s->generalized)memcpy(out,in,s->n*sizeof(double));
 else{
  for(size_t i=0;i<s->n;i++)out[i]=NAN;
  e=s->apply(s->user,op,s->n,in,out);
  if(e){s->code=e;return LIS_RCI_CALLBACK;}
 }
 for(size_t i=0;i<s->n;i++)if(!isfinite(out[i]))return LIS_RCI_NONFINITE;
 return check_gen(s);
}
static int store(double *out,long double v)
{
 if(!isfinite(v)||fabsl(v)>DBL_MAX)return LIS_RCI_NONFINITE;
 *out=(double)v;return 0;
}
static int shifted_action(void *user,const lis_rci_opaque_request *q)
{
 lis_rci_split_restore *s=user;
 const double *in=q->in[0];double *out=q->out[0];size_t n=s->n;
 if(q->operation==LIS_RCI_OP_PRECON){
  memcpy(out,in,2*n*sizeof(double));return 0;
 }
 if(q->operation!=LIS_RCI_OP_A)return LIS_RCI_BADARG;
 int e=real_apply(s,LIS_RCI_OP_A,in,s->ar);if(e)return e;
 e=real_apply(s,LIS_RCI_OP_A,in+n,s->ai);if(e)return e;
 e=real_apply(s,LIS_RCI_OP_B,in,s->br);if(e)return e;
 e=real_apply(s,LIS_RCI_OP_B,in+n,s->bi);if(e)return e;
 for(size_t i=0;i<n;i++){
  long double r=s->shift.real*s->br[i]-s->shift.imag*s->bi[i]-s->ar[i];
  long double im=s->shift.imag*s->br[i]+s->shift.real*s->bi[i]-s->ai[i];
  if(store(out+i,r)||store(out+n+i,im))return LIS_RCI_NONFINITE;
 }
 return 0;
}
static int solve_shift(lis_rci_split_restore *s,const lis_rci_opaque_request *q,
 const lis_rci_split_vector *in,lis_rci_split_vector *out)
{
 size_t n=s->n;s->shift=q->shift;
 memcpy(s->rhs,in->real,n*sizeof(double));
 memcpy(s->rhs+n,in->imag,n*sizeof(double));
 lis_rci_handle rhs=s->rhs,sol=s->solution;
 lis_rci_opaque_request parent=*q;
 /* The child's real vector is the concatenation of real and imaginary planes.
  * The backend owns the actual complex shift in s->shift. */
 parent.shift=(lis_rci_complex_coeff){0,0};
 parent.in=&rhs;parent.out=&sol;
 lis_rci_opaque_backend b={s,NULL,shifted_action,split_gens};
 int code=0;
 int e=lis_rci_nested_restore_run(s->child,&parent,sizeof parent,2*n,&b,&code);
 if(e==LIS_RCI_CALLBACK){
  if(s->code)return LIS_RCI_CALLBACK;
  /* shifted_action returns RCI statuses for internal validation failures. */
  if(code==LIS_RCI_NONFINITE||code==LIS_RCI_PROTOCOL||code==LIS_RCI_BADARG)return code;
  s->code=code;
 }
 if(e)return e;
 e=check_gen(s);if(e)return e;
 memcpy(out->real,s->solution,n*sizeof(double));
 memcpy(out->imag,s->solution+n,n*sizeof(double));
 return 0;
}
static int split_single(void *user,const lis_rci_opaque_request *q)
{
 lis_rci_split_restore *s=user;
 s->code=0;s->status=LIS_RCI_BADARG;
 if(lis_rci_opaque_validate(q,sizeof(*q))||q->columns!=1)return LIS_RCI_BADARG;
 lis_rci_split_vector *out=q->out[0];
 const lis_rci_split_vector *a=q->in?q->in[0]:NULL,*b=q->in2?q->in2[0]:NULL;
 if(!vector_ok(out,s->n)||(a&&(!vector_ok(a,s->n)||alias(out,a)))||
    (b&&(!vector_ok(b,s->n)||alias(out,b))))return LIS_RCI_BADARG;
 s->a=q->generation_a;s->b=q->generation_b;
 int e=check_gen(s);if(e){s->status=e;return s->code?s->code:e;}
 size_t n=s->n;
 switch(q->operation){
 case LIS_RCI_OP_SHIFTED_SOLVE:e=solve_shift(s,q,a,out);break;
 case LIS_RCI_OP_ZERO:
  memset(out->real,0,n*sizeof(double));memset(out->imag,0,n*sizeof(double));break;
 case LIS_RCI_OP_A:case LIS_RCI_OP_B:
  e=real_apply(s,q->operation,a->real,s->ar);
  if(!e)e=real_apply(s,q->operation,a->imag,s->ai);
  if(!e){memcpy(out->real,s->ar,n*sizeof(double));memcpy(out->imag,s->ai,n*sizeof(double));}
  break;
 case LIS_RCI_OP_PRECON:
  memcpy(out->real,a->real,n*sizeof(double));memcpy(out->imag,a->imag,n*sizeof(double));break;
 case LIS_RCI_OP_AXPBY:
  for(size_t i=0;i<n;i++){
   long double r=q->alpha.real*a->real[i]-q->alpha.imag*a->imag[i]
      +q->beta.real*b->real[i]-q->beta.imag*b->imag[i];
   long double im=q->alpha.imag*a->real[i]+q->alpha.real*a->imag[i]
      +q->beta.imag*b->real[i]+q->beta.real*b->imag[i];
   if(store(s->ar+i,r)||store(s->ai+i,im)){e=LIS_RCI_NONFINITE;break;}
  }
  if(!e){memcpy(out->real,s->ar,n*sizeof(double));memcpy(out->imag,s->ai,n*sizeof(double));}
  break;
 default:e=LIS_RCI_BADARG;
 }
 if(!e)e=check_gen(s);
 s->status=e;return s->code?s->code:e;
}
lis_rci_split_restore *lis_rci_split_restore_create(size_t n,int generalized,
 const lis_rci_options *o,void *user,lis_rci_real_apply apply,lis_rci_generation_read generations)
{
 if(!n||n>SIZE_MAX/8/sizeof(double)||(generalized!=0&&generalized!=1)||!o||
    o->n!=2*n||o->method<LIS_RCI_BICGSTAB||o->method>LIS_RCI_FGMRES||!apply||!generations)
  return NULL;
 lis_rci_split_restore *s=calloc(1,sizeof(*s));if(!s)return NULL;
 s->n=n;s->generalized=generalized;s->user=user;s->apply=apply;s->generations=generations;
 s->status=LIS_RCI_STATE;s->child=lis_rci_create(o);s->mem=malloc(8*n*sizeof(double));
 if(!s->child||!s->mem){lis_rci_split_restore_destroy(s);return NULL;}
 s->rhs=s->mem;s->solution=s->rhs+2*n;s->ar=s->solution+2*n;
 s->ai=s->ar+n;s->br=s->ai+n;s->bi=s->br+n;
 return s;
}
void lis_rci_split_restore_destroy(lis_rci_split_restore *s)
{if(s){lis_rci_destroy(s->child);free(s->mem);free(s);}}
int lis_rci_split_restore_backend(lis_rci_split_restore *s,lis_rci_opaque_backend *b)
{
 if(!s||!b)return LIS_RCI_BADARG;
 *b=(lis_rci_opaque_backend){s,NULL,split_single,split_gens};return 0;
}
int lis_rci_split_restore_error(const lis_rci_split_restore *s,int *code)
{
 if(!s||!code)return LIS_RCI_BADARG;
 *code=s->code;return s->status;
}
