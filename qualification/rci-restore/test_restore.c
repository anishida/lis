/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_restore.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#define C(x) do{if(!(x)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#x);return 1;}}while(0)
typedef struct { long double complex v[3]; } vec;
typedef struct { uint64_t a,b; unsigned calls,fail,change; int generr, generalized; } env;
static int gens(void *u,uint64_t *a,uint64_t *b)
{env *e=u;*a=e->a;*b=e->b;return e->generr;}
static long double complex z(lis_rci_complex_coeff c){return c.real+c.imag*I;}
static int perform(void *u,const lis_rci_opaque_request *q)
{
 env *e=u;e->calls++;
 if(e->fail==e->calls)return -73;
 if(e->change==e->calls)e->a++;
 for(size_t j=0;j<q->columns;j++){
  vec *o=q->out[j];const vec *a=q->in?q->in[j]:NULL,*b=q->in2?q->in2[j]:NULL;
  for(size_t k=0;k<3;k++){
   long double d=(long double)k+2,B=e->generalized?1+k:1;
   switch(q->operation){
   case LIS_RCI_OP_ZERO:o->v[k]=0;break;
   case LIS_RCI_OP_B:o->v[k]=B*a->v[k];break;
   case LIS_RCI_OP_SHIFTED_SOLVE:o->v[k]=a->v[k]/(z(q->shift)*B-d);break;
   case LIS_RCI_OP_AXPBY:o->v[k]=z(q->alpha)*a->v[k]+z(q->beta)*b->v[k];break;
   default:return -99;
   }
  }
 }
 return 0;
}
static int projector(int general,int use_block)
{
 vec v[10]={0};lis_rci_handle h[10];
 for(size_t i=0;i<10;i++)h[i]=&v[i];
 for(size_t j=0;j<2;j++)for(size_t k=0;k<3;k++)v[j].v[k]=(1+j+k)+(2*j-k*.5L)*I;
 lis_rci_complex_coeff shift[]={{1,2},{4,-3},{6,.5L}},weight[]={{.5L,.2L},{-.3L,.7L},{.2L,-.1L}};
 lis_rci_projector_restore_options o={2,3,general,shift,weight,h,h+2,h+4,h+6,h+8};
 lis_rci_projector_restore *p=lis_rci_projector_restore_create(&o);C(p);
 env e={.a=UINT64_C(9007199254740995),.b=8,.generalized=general};
 lis_rci_opaque_backend b={&e,use_block?perform:NULL,use_block?NULL:perform,gens};
 lis_rci_opaque_request q;int code=0;uint64_t last=0;unsigned requests=0;
 for(unsigned repeat=0;repeat<2;repeat++){
  int s=lis_rci_projector_restore_start(p,e.a,e.b,&q);
  C(s==LIS_RCI_OP_ZERO&&lis_rci_projector_restore_result(p)==NULL);
  while(s>0){
   C(q.token>last);last=q.token;requests++;
   C(lis_rci_projector_restore_advance(p,q.token-1,e.a,e.b,0,&q)==LIS_RCI_PROTOCOL);
   C(q.token==last);C(lis_rci_opaque_dispatch(&b,&q,sizeof q,&code)==0&&code==0);
   s=lis_rci_projector_restore_advance(p,q.token,e.a,e.b,0,&q);
  }
  C(s==0);const lis_rci_handle *result=lis_rci_projector_restore_result(p);C(result);
  for(size_t j=0;j<2;j++)for(size_t k=0;k<3;k++){
   long double B=general?1+k:1;
   long double complex want=0;
   for(size_t i=0;i<3;i++)want+=z(weight[i])*B*v[j].v[k]/(z(shift[i])*B-(k+2));
   C(cabsl(((vec*)result[j])->v[k]-want)<1e-15L);
  }
 }
 C(requests==2u*(7u+(unsigned)general));C(e.calls==requests*(use_block?1:2));
 /* Inject callback failure, stale generation and cancellation at every phase. */
 for(unsigned stop=0;stop<7+(unsigned)general;stop++)for(int fault=0;fault<3;fault++){
  int s=lis_rci_projector_restore_start(p,e.a,e.b,&q);C(s>0);
  for(unsigned i=0;i<stop;i++){
   C(lis_rci_opaque_dispatch(&b,&q,sizeof q,&code)==0);
   s=lis_rci_projector_restore_advance(p,q.token,e.a,e.b,0,&q);C(s>0);
  }
  if(fault==0){e.fail=e.calls+1;C(lis_rci_opaque_dispatch(&b,&q,sizeof q,&code)==LIS_RCI_CALLBACK&&code==-73);e.fail=0;
   C(lis_rci_projector_restore_advance(p,q.token,e.a,e.b,code,&q)==LIS_RCI_CALLBACK);
   C(lis_rci_projector_restore_error(p,&code)==LIS_RCI_CALLBACK&&code==-73);
  }else if(fault==1){e.change=e.calls+1;C(lis_rci_opaque_dispatch(&b,&q,sizeof q,&code)==LIS_RCI_PROTOCOL);e.change=0;
   C(lis_rci_projector_restore_advance(p,q.token,e.a,e.b,0,&q)==LIS_RCI_PROTOCOL);
  }else{C(lis_rci_projector_restore_cancel(p,q.token-1)==LIS_RCI_PROTOCOL);C(lis_rci_projector_restore_cancel(p,q.token)==LIS_RCI_CANCELLED);}
  C(!lis_rci_projector_restore_result(p));
 }
 lis_rci_projector_restore_destroy(p);
 h[9]=h[0];C(!lis_rci_projector_restore_create(&o));
 return 0;
}
static int real_apply(void *u,const lis_rci_opaque_request *q)
{
 env *e=u;e->calls++;if(e->fail==e->calls)return -73;
 double *out=q->out[0];const double *in=q->in[0];
 for(size_t i=0;i<3;i++)out[i]=q->operation==LIS_RCI_OP_A?(double)(q->shift.real-(i+1))*in[i]:in[i];
 return 0;
}
static int nested(void)
{
 env e={.a=7,.b=11};lis_rci_opaque_backend b={&e,NULL,real_apply,gens};
 double rhs[]={2,3,5},out[]={-1,-1,-1};lis_rci_handle ih=rhs,oh=out;
 lis_rci_opaque_request q={.version=1,.size=sizeof(q),.operation=LIS_RCI_OP_SHIFTED_SOLVE,.token=UINT64_MAX-7,.generation_a=7,.generation_b=11,.columns=1,.in=&ih,.out=&oh,.shift={6,0}};
 int code=0;
 C(!lis_rci_opaque_validate(&q,sizeof q));C(lis_rci_opaque_validate(&q,7)==LIS_RCI_BADARG);
 q.version=2;C(lis_rci_opaque_validate(&q,sizeof q)==LIS_RCI_BADARG);q.version=1;
 q.size--;C(lis_rci_opaque_validate(&q,sizeof q)==LIS_RCI_BADARG);q.size++;
 q.flags=1;C(lis_rci_opaque_validate(&q,sizeof q)==LIS_RCI_BADARG);q.flags=0;
 q.token=0;C(lis_rci_opaque_validate(&q,sizeof q)==LIS_RCI_BADARG);q.token=UINT64_MAX-7;
 q.out=&ih;C(lis_rci_opaque_validate(&q,sizeof q)==LIS_RCI_BADARG);q.out=&oh;
 for(int method=LIS_RCI_CG;method<=LIS_RCI_FGMRES;method++){
  lis_rci_options o={3,20,3,1e-12,0,(lis_rci_method)method};
  lis_rci_context *p=lis_rci_create(&o);C(p);
  C(lis_rci_nested_restore_run(p,&q,sizeof q,2,&b,&code)==LIS_RCI_BADARG);
  C(!lis_rci_nested_restore_run(p,&q,sizeof q,3,&b,&code));
  for(size_t i=0;i<3;i++)C(fabs(out[i]-rhs[i]/(5-i))<1e-10);
  unsigned requests=e.calls;e.fail=e.calls+1;out[0]=out[1]=out[2]=-1;
  C(lis_rci_nested_restore_run(p,&q,sizeof q,3,&b,&code)==LIS_RCI_CALLBACK&&code==-73);
  C(e.calls==requests+1&&out[0]==-1&&out[1]==-1&&out[2]==-1);e.fail=0;
  C(!lis_rci_nested_restore_run(p,&q,sizeof q,3,&b,&code));C(q.token==UINT64_MAX-7);
  e.a++;C(lis_rci_nested_restore_run(p,&q,sizeof q,3,&b,&code)==LIS_RCI_PROTOCOL);e.a--;
  lis_rci_destroy(p);
 }
 return 0;
}
int main(void)
{
 for(int i=0;i<2;i++)for(int j=0;j<2;j++)C(!projector(i,j));
 C(!nested());
 puts("RESTORE PASS: 8 block/column standard/generalized complex projections, 90 phase faults, 4 nested methods, descriptor/generation guards");
 return 0;
}
