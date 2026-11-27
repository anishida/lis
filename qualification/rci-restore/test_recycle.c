/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_recycle_restore.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define N 5
#define K 3
#define C(e) do{if(!(e)){fprintf(stderr,"RECYCLE FAIL %d: %s\n",__LINE__,#e);return 1;}}while(0)
static void op(uint32_t code,const double *in,double *out,int changed)
{
 for(int i=0;i<N;i++){
  if(code==LIS_RCI_OP_PRECON)out[i]=in[i]/(2+i+changed);
  else out[i]=(4+i+changed)*in[i]+(i?-.7*in[i-1]:0)+(i+1<N?.2*in[i+1]:0);
 }
}
static int drive(lis_rci_recycle_restore *s,int status,lis_rci_recycle_request *q,int changed,uint64_t *last)
{
 int count=0;
 while(status>0){
  if(q->token<=*last||++count>2*K)return 1;
  *last=q->token;
  if(lis_rci_recycle_restore_advance(s,q->token-1,q->operator_generation,q->precon_generation,0,q)!=LIS_RCI_PROTOCOL)return 1;
  op(q->operation,q->in,q->out,changed);
  status=lis_rci_recycle_restore_advance(s,q->token,q->operator_generation,q->precon_generation,0,q);
 }
 return status;
}
static int check_basis(lis_rci_recycle_restore *s,int flex,int changed)
{
 lis_rci_recycle_basis b;C(!lis_rci_recycle_restore_get(s,&b));C(b.count==2);
 double out[N];
 for(size_t j=0;j<b.count;j++){
  op(LIS_RCI_OP_A,b.physical+j*N,out,changed);
  for(int i=0;i<N;i++)C(fabs(out[i]-b.image[j*N+i])<1e-12);
  if(!flex){
   op(LIS_RCI_OP_PRECON,b.domain+j*N,out,changed);
   for(int i=0;i<N;i++)C(fabs(out[i]-b.physical[j*N+i])<1e-12);
  }
  for(size_t k=0;k<b.count;k++){
   double d=0;for(int i=0;i<N;i++)d+=b.image[j*N+i]*b.image[k*N+i];
   C(fabs(d-(j==k))<1e-12);
  }
 }
 return 0;
}
static int run(int flex)
{
 double w[K*N]={1,.3,0,.4,0, .2,1,1,0,.3, 1,.3,0,.4,0};
 double z[K*N];for(int j=0;j<K;j++)for(int i=0;i<N;i++)z[j*N+i]=w[j*N+i]/(3+i+(j==1?.5:0));
 lis_rci_recycle_restore *s=lis_rci_recycle_restore_create(N,K,flex,1e-12);C(s);
 lis_rci_recycle_request q;uint64_t last=0;int status;
 status=lis_rci_recycle_restore_begin(s,K,w,flex?z:NULL,7,9,&q);
 C(status==(flex?LIS_RCI_OP_A:LIS_RCI_OP_PRECON));
 C(!drive(s,status,&q,0,&last));C(!check_basis(s,flex,0));
 lis_rci_recycle_basis b;C(!lis_rci_recycle_restore_get(s,&b));
 double saved[6*N];memcpy(saved,b.domain,2*N*sizeof(double));memcpy(saved+2*N,b.physical,2*N*sizeof(double));memcpy(saved+4*N,b.image,2*N*sizeof(double));
 for(int rep=0;rep<2;rep++){
  double exact[N],x[N],r[N],rhs[N],ax[N];
  for(int i=0;i<N;i++){exact[i]=sin(i+rep);x[i]=.1;}
  op(LIS_RCI_OP_A,exact,rhs,0);op(LIS_RCI_OP_A,x,ax,0);
  for(int i=0;i<N;i++)r[i]=rhs[i]-ax[i];
  C(!lis_rci_recycle_restore_project(s,7,9,x,r));
  op(LIS_RCI_OP_A,x,ax,0);
  for(int i=0;i<N;i++)C(fabs(r[i]-(rhs[i]-ax[i]))<1e-12);
  for(size_t j=0;j<b.count;j++){double d=0;for(int i=0;i<N;i++)d+=r[i]*b.image[j*N+i];C(fabs(d)<1e-12);}
  double v[N],p[N],av[N];for(int i=0;i<N;i++)v[i]=cos(i+.3);
  op(LIS_RCI_OP_PRECON,v,p,0);op(LIS_RCI_OP_A,p,av,0);
  C(!lis_rci_recycle_restore_direction(s,7,9,v,p,av));op(LIS_RCI_OP_A,p,ax,0);
  for(int i=0;i<N;i++)C(fabs(av[i]-ax[i])<1e-12);
  if(!flex){op(LIS_RCI_OP_PRECON,v,ax,0);for(int i=0;i<N;i++)C(fabs(ax[i]-p[i])<1e-12);}
 }
 /* Each injected error preserves all three committed arrays byte-for-byte. */
 for(int refresh=0;refresh<2;refresh++)for(int stop=0;stop<(flex?1:2)*(refresh?2:K);stop++)for(int fault=0;fault<4;fault++){
  status=refresh?lis_rci_recycle_restore_refresh(s,8,9,&q):lis_rci_recycle_restore_begin(s,K,w,flex?z:NULL,8,9,&q);C(status>0);
  for(int i=0;i<stop;i++){op(q.operation,q.in,q.out,1);status=lis_rci_recycle_restore_advance(s,q.token,8,9,0,&q);C(status>0);}
  if(fault==0){C(lis_rci_recycle_restore_advance(s,q.token,8,9,71,&q)==LIS_RCI_CALLBACK);}
  else if(fault==1){C(lis_rci_recycle_restore_advance(s,q.token,9,9,0,&q)==LIS_RCI_PROTOCOL);}
  else if(fault==2){C(lis_rci_recycle_restore_cancel(s,q.token)==LIS_RCI_CANCELLED);}
  else{q.out[0]=0;C(lis_rci_recycle_restore_advance(s,q.token,8,9,0,&q)==LIS_RCI_NONFINITE);}
  C(!lis_rci_recycle_restore_get(s,&b)&&b.operator_generation==7&&b.precon_generation==9);
  C(!memcmp(saved,b.domain,2*N*sizeof(double))&&!memcmp(saved+2*N,b.physical,2*N*sizeof(double))&&!memcmp(saved+4*N,b.image,2*N*sizeof(double)));
 }
 if(flex)C(lis_rci_recycle_restore_refresh(s,8,10,&q)==LIS_RCI_PROTOCOL);
 status=lis_rci_recycle_restore_refresh(s,8,flex?9:10,&q);C(status>0);C(!drive(s,status,&q,1,&last));C(!check_basis(s,flex,1));
 C(!lis_rci_recycle_restore_get(s,&b)&&b.operator_generation==8);
 lis_rci_recycle_restore_destroy(s);return 0;
}
int main(void)
{
 C(!run(0));C(!run(1));
 puts("RECYCLE PASS: fixed U/ZU/C and flexible W/Z/C, rank drop, two RHS, direction projection, refresh, 60 transactional phase faults");
 return 0;
}
