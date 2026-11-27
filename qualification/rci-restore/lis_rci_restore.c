/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_restore.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int coeff_ok(lis_rci_complex_coeff c)
{ return isfinite(c.real)&&isfinite(c.imag); }
static int handles_ok(const lis_rci_handle *a,size_t n)
{
 if(!a)return 0;
 for(size_t i=0;i<n;i++)if(!a[i])return 0;
 return 1;
}
static int outputs_ok(const lis_rci_handle *out,const lis_rci_handle *a,
 const lis_rci_handle *b,size_t n)
{
 for(size_t i=0;i<n;i++){
  for(size_t j=0;j<n;j++){
   if((a&&out[i]==a[j])||(b&&out[i]==b[j])||(i!=j&&out[i]==out[j]))return 0;
  }
 }
 return 1;
}
int lis_rci_opaque_validate(const lis_rci_opaque_request *q,size_t bytes)
{
 uint32_t prefix[2];
 if(!q||bytes<sizeof prefix)return LIS_RCI_BADARG;
 memcpy(prefix,q,sizeof prefix);
 if(prefix[0]!=LIS_RCI_OPAQUE_VERSION||prefix[1]!=sizeof(*q)||bytes<sizeof(*q))
  return LIS_RCI_BADARG;
 if(!q->token||q->flags||!q->columns||q->columns>SIZE_MAX/sizeof(lis_rci_handle)
    ||!coeff_ok(q->shift)||!coeff_ok(q->alpha)||!coeff_ok(q->beta)
    ||!handles_ok(q->out,q->columns))return LIS_RCI_BADARG;
 switch(q->operation){
 case LIS_RCI_OP_ZERO:
  if(q->in||q->in2)return LIS_RCI_BADARG;
  break;
 case LIS_RCI_OP_AXPBY:
  if(!handles_ok(q->in,q->columns)||!handles_ok(q->in2,q->columns))return LIS_RCI_BADARG;
  break;
 case LIS_RCI_OP_A:case LIS_RCI_OP_B:case LIS_RCI_OP_PRECON:case LIS_RCI_OP_SHIFTED_SOLVE:
  if(!handles_ok(q->in,q->columns)||q->in2)return LIS_RCI_BADARG;
  break;
 default:return LIS_RCI_BADARG;
 }
 return outputs_ok(q->out,q->in,q->in2,q->columns)?LIS_RCI_DONE:LIS_RCI_BADARG;
}
static int generation_check(const lis_rci_opaque_backend *b,
 const lis_rci_opaque_request *q,int *code)
{
 uint64_t a=0,c=0;
 int e=b->generations(b->user,&a,&c);
 if(e){*code=e;return LIS_RCI_CALLBACK;}
 return a==q->generation_a&&c==q->generation_b?LIS_RCI_DONE:LIS_RCI_PROTOCOL;
}
int lis_rci_opaque_dispatch(const lis_rci_opaque_backend *b,
 const lis_rci_opaque_request *q,size_t bytes,int *code)
{
 if(!code)return LIS_RCI_BADARG;
 *code=0;
 if(!b||!b->generations||(!b->block&&!b->single)||lis_rci_opaque_validate(q,bytes))
  return LIS_RCI_BADARG;
 size_t n=b->block?1:q->columns;
 for(size_t i=0;i<n;i++){
  int s=generation_check(b,q,code);
  if(s)return s;
  lis_rci_opaque_request one=*q;
  if(!b->block){
   one.columns=1;one.out=q->out+i;
   if(q->in)one.in=q->in+i;
   if(q->in2)one.in2=q->in2+i;
  }
  int e=(b->block?b->block:b->single)(b->user,&one);
  if(e){*code=e;return LIS_RCI_CALLBACK;}
  s=generation_check(b,q,code);
  if(s)return s;
 }
 return LIS_RCI_DONE;
}

enum stage { P_IDLE,P_ZERO,P_B,P_SOLVE,P_ACCUM,P_DONE,P_FAILED };
struct lis_rci_projector_restore {
 size_t m,p,j,active;
 int generalized,status,code;
 enum stage stage;
 uint64_t token,a,b;
 lis_rci_complex_coeff *nodes;
 lis_rci_handle *handles;
};
static lis_rci_handle *block(lis_rci_projector_restore *p,size_t i)
{ return p->handles+i*p->m; }
static int finish(lis_rci_projector_restore *p,int s,lis_rci_opaque_request *q)
{
 p->stage=s==LIS_RCI_DONE?P_DONE:P_FAILED;p->status=s;
 if(q)memset(q,0,sizeof(*q));
 return s;
}
static int issue(lis_rci_projector_restore *p,enum stage stage,lis_rci_opaque_request *q)
{
 if(p->token==UINT64_MAX)return finish(p,LIS_RCI_PROTOCOL,q);
 lis_rci_opaque_request v={0};
 v.version=LIS_RCI_OPAQUE_VERSION;v.size=sizeof(v);v.token=++p->token;
 v.generation_a=p->a;v.generation_b=p->b;v.columns=p->m;
 switch(stage){
 case P_ZERO:v.operation=LIS_RCI_OP_ZERO;v.out=block(p,3);break;
 case P_B:v.operation=LIS_RCI_OP_B;v.in=block(p,0);v.out=block(p,1);break;
 case P_SOLVE:v.operation=LIS_RCI_OP_SHIFTED_SOLVE;
  v.shift=p->nodes[p->j];v.in=block(p,p->generalized?1:0);v.out=block(p,2);break;
 case P_ACCUM:v.operation=LIS_RCI_OP_AXPBY;
  v.in=block(p,3+p->active);v.in2=block(p,2);v.out=block(p,4-p->active);
  v.alpha.real=1;v.beta=p->nodes[p->p+p->j];break;
 default:return finish(p,LIS_RCI_STATE,q);
 }
 p->stage=stage;p->status=(int)v.operation;*q=v;return p->status;
}
lis_rci_projector_restore *lis_rci_projector_restore_create(const lis_rci_projector_restore_options *o)
{
 if(!o||!o->columns||!o->nodes||(o->generalized!=0&&o->generalized!=1)
    ||o->columns>SIZE_MAX/5/sizeof(lis_rci_handle)
    ||o->nodes>SIZE_MAX/2/sizeof(lis_rci_complex_coeff)||!o->shift||!o->weight)
  return NULL;
 const lis_rci_handle *arr[]={o->q,o->r,o->x,o->y0,o->y1};
 for(size_t k=0;k<5;k++)if(!handles_ok(arr[k],o->columns))return NULL;
 for(size_t k=0;k<o->nodes;k++)if(!coeff_ok(o->shift[k])||!coeff_ok(o->weight[k]))return NULL;
 lis_rci_projector_restore *p=calloc(1,sizeof(*p));if(!p)return NULL;
 p->m=o->columns;p->p=o->nodes;p->generalized=o->generalized;p->status=LIS_RCI_STATE;
 p->nodes=malloc(2*p->p*sizeof(*p->nodes));p->handles=malloc(5*p->m*sizeof(*p->handles));
 if(!p->nodes||!p->handles){lis_rci_projector_restore_destroy(p);return NULL;}
 memcpy(p->nodes,o->shift,p->p*sizeof(*p->nodes));
 memcpy(p->nodes+p->p,o->weight,p->p*sizeof(*p->nodes));
 for(size_t k=0;k<5;k++)memcpy(block(p,k),arr[k],p->m*sizeof(*p->handles));
 if(!outputs_ok(p->handles,NULL,NULL,5*p->m)){lis_rci_projector_restore_destroy(p);return NULL;}
 return p;
}
void lis_rci_projector_restore_destroy(lis_rci_projector_restore *p)
{ if(p){free(p->handles);free(p->nodes);free(p);} }
int lis_rci_projector_restore_start(lis_rci_projector_restore *p,uint64_t a,uint64_t b,lis_rci_opaque_request *q)
{
 if(!p||!q)return LIS_RCI_BADARG;
 if(p->stage!=P_IDLE&&p->stage!=P_DONE&&p->stage!=P_FAILED)return LIS_RCI_STATE;
 p->a=a;p->b=b;p->j=0;p->active=0;p->code=0;
 return issue(p,P_ZERO,q);
}
int lis_rci_projector_restore_advance(lis_rci_projector_restore *p,uint64_t token,uint64_t a,uint64_t b,int code,lis_rci_opaque_request *q)
{
 if(!p||!q)return LIS_RCI_BADARG;
 if(p->stage==P_IDLE||p->stage==P_DONE||p->stage==P_FAILED)return LIS_RCI_STATE;
 if(token!=p->token)return LIS_RCI_PROTOCOL;
 if(a!=p->a||b!=p->b)return finish(p,LIS_RCI_PROTOCOL,q);
 if(code){p->code=code;return finish(p,LIS_RCI_CALLBACK,q);}
 switch(p->stage){
 case P_ZERO:return issue(p,p->generalized?P_B:P_SOLVE,q);
 case P_B:return issue(p,P_SOLVE,q);
 case P_SOLVE:return issue(p,P_ACCUM,q);
 case P_ACCUM:
  p->active=1-p->active;
  if(++p->j==p->p)return finish(p,LIS_RCI_DONE,q);
  return issue(p,P_SOLVE,q);
 default:return finish(p,LIS_RCI_STATE,q);
 }
}
int lis_rci_projector_restore_cancel(lis_rci_projector_restore *p,uint64_t token)
{
 if(!p)return LIS_RCI_BADARG;
 if(p->stage==P_IDLE||p->stage==P_DONE||p->stage==P_FAILED)return LIS_RCI_STATE;
 if(token!=p->token)return LIS_RCI_PROTOCOL;
 return finish(p,LIS_RCI_CANCELLED,NULL);
}
int lis_rci_projector_restore_error(const lis_rci_projector_restore *p,int *code)
{
 if(!p||!code)return LIS_RCI_BADARG;
 *code=p->code;return p->status;
}
const lis_rci_handle *lis_rci_projector_restore_result(const lis_rci_projector_restore *p)
{ return p&&p->stage==P_DONE?p->handles+(3+p->active)*p->m:NULL; }

int lis_rci_nested_restore_run(lis_rci_context *child,const lis_rci_opaque_request *parent,
 size_t bytes,size_t n,const lis_rci_opaque_backend *backend,int *code)
{
 if(!code)return LIS_RCI_BADARG;
 *code=0;
 /* Current real-double solver cannot represent complex shifted systems. */
 if(!child||!n||n>SIZE_MAX/sizeof(double)||!backend||!backend->generations
    ||(!backend->block&&!backend->single)||lis_rci_opaque_validate(parent,bytes)
    ||parent->operation!=LIS_RCI_OP_SHIFTED_SOLVE||parent->columns!=1||parent->shift.imag!=0)
  return LIS_RCI_BADARG;
 int s=generation_check(backend,parent,code);
 if(s)return s;
 /* The child dimension must be checked BEFORE reading the parent's data.
  * The existing core has no public dimension accessor; use the separately
  * supplied child dimension accessor in this reconstruction. */
 extern size_t lis_rci_restore_dimension(const lis_rci_context *);
 if(lis_rci_restore_dimension(child)!=n)return LIS_RCI_BADARG;
 lis_rci_request q;
 s=lis_rci_start(child,parent->in[0],NULL,&q);
 while(s==LIS_RCI_MATVEC||s==LIS_RCI_PRECON){
  lis_rci_request_full pending;
  if(lis_rci_get_pending(child,&pending)!=s)return LIS_RCI_STATE;
  lis_rci_handle in=(void *)q.in,out=q.out;
  lis_rci_opaque_request r=*parent;
  r.operation=s==LIS_RCI_MATVEC?LIS_RCI_OP_A:LIS_RCI_OP_PRECON;
  r.token=pending.request.token;r.in=&in;r.in2=NULL;r.out=&out;
  int e=0;int d=lis_rci_opaque_dispatch(backend,&r,sizeof r,&e);
  if(d){
   if(d==LIS_RCI_CALLBACK){*code=e;(void)lis_rci_advance(child,e,&q);}
   else (void)lis_rci_cancel(child,pending.request.token);
   return d;
  }
  s=lis_rci_advance(child,0,&q);
 }
 if(s!=LIS_RCI_DONE)return s;
 s=generation_check(backend,parent,code);if(s)return s;
 memcpy(parent->out[0],lis_rci_solution(child),n*sizeof(double));
 return LIS_RCI_DONE;
}
