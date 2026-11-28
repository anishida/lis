#include "lis_rci.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"RCI lifecycle line %d: %s\n",__LINE__,#c);return 1;}}while(0)
static int exercise(int method,int mode)
{
 lis_rci_options o={3,100,3,1e-12,1e-14,(lis_rci_method)method};
 lis_rci_context *c=lis_rci_create(&o);CHECK(c);
 lis_rci_request_full f;lis_rci_stats stats;lis_rci_info info;
 CHECK(lis_rci_get_pending(c,&f)==LIS_RCI_STATE&&f.request.token==0);
 CHECK(lis_rci_cancel(c,0)==LIS_RCI_STATE);
 double b[]={2,4,6};uint64_t prior=0;
 for(int rep=0;rep<2;rep++){
  int s;
  if(mode==0){lis_rci_request q;s=lis_rci_start(c,b,NULL,&q);}
  else if(mode==1){lis_rci_request_ex q;s=lis_rci_start_ex(c,b,NULL,&q);}
  else s=lis_rci_start_full(c,b,NULL,&f);
  CHECK(s==(mode?LIS_RCI_NORM:LIS_RCI_MATVEC));
  CHECK(lis_rci_get_pending(c,&f)==s&&f.request.token>prior);
  uint64_t token=f.request.token;CHECK(lis_rci_cancel(c,token-1)==LIS_RCI_PROTOCOL);
  CHECK(lis_rci_get_pending(c,&f)==s&&f.request.token==token);
  if(rep==0){
   CHECK(lis_rci_cancel(c,token)==LIS_RCI_CANCELLED);
   CHECK(lis_rci_get_info(c,&info)==0&&info.status==LIS_RCI_CANCELLED);
   CHECK(lis_rci_get_stats(c,&stats)==0&&stats.cancelled==1&&stats.last_token==token);
   CHECK(lis_rci_get_pending(c,&f)==LIS_RCI_STATE&&f.request.token==0);
   prior=token;continue;
  }
  int count=0;
  while(s>0){
   CHECK(++count<1000&&lis_rci_get_pending(c,&f)==s);
   lis_rci_request_ex *q=&f.request;uint64_t t=q->token;long double value=0;
   CHECK(t>prior);prior=t;
   if(q->out)for(size_t i=0;i<3;i++){
    if(s==LIS_RCI_MATVEC)q->out[i]=2*q->in[i];
    else if(s==LIS_RCI_PRECON)q->out[i]=q->in[i];
    else if(s==LIS_RCI_AXPBY)q->out[i]=(double)(f.alpha*q->in[i]+(q->in2?f.beta*q->in2[i]:0));
    else {CHECK(s==LIS_RCI_DIVIDE);q->out[i]=(double)(q->in[i]/f.alpha);}
   }else{
    for(size_t i=0;i<3;i++)value+=(long double)q->in[i]*(q->in2?q->in2[i]:q->in[i]);
    if(s==LIS_RCI_NORM)value=sqrtl(value);
   }
   if(!mode){lis_rci_request legacy;s=lis_rci_advance(c,0,&legacy);}
   else {
    uint64_t words[5],other[5];
    CHECK(!lis_rci_consensus_pack(t,(uint64_t)s,0,0,words));
    CHECK(!lis_rci_consensus_pack(t,(uint64_t)s,9,0,other));
    CHECK(!lis_rci_consensus_merge(words,other));
    CHECK(lis_rci_advance_consensus(c,t-1,words,value,&f)==LIS_RCI_PROTOCOL);
    CHECK(f.request.token==t);
    s=lis_rci_advance_consensus(c,t,words,value,&f);
   }
  }
  CHECK(s==LIS_RCI_DONE);
  const double *x=lis_rci_solution(c);CHECK(fabs(x[0]-1)<1e-12&&fabs(x[1]-2)<1e-12&&fabs(x[2]-3)<1e-12);
  CHECK(lis_rci_get_stats(c,&stats)==0&&!stats.cancelled);
  CHECK(lis_rci_get_info(c,&info)==0);
  CHECK(stats.issued[0]==info.matvecs&&stats.issued[1]==info.precons);
  for(size_t i=0;i<6;i++)CHECK(stats.issued[i]==stats.completed[i]);
  CHECK(stats.issued[3]>0&&stats.issued[0]>0);
 }
 if(mode){
  for(int mismatch=0;mismatch<3;mismatch++){
   CHECK(lis_rci_start_full(c,b,NULL,&f)==LIS_RCI_NORM);
   uint64_t t=f.request.token,w[5],other[5];
   CHECK(!lis_rci_consensus_pack(t,LIS_RCI_NORM,5,-18,w));
   CHECK(!lis_rci_consensus_pack(t+(mismatch==1),LIS_RCI_NORM+(mismatch==2),2,73,other));
   CHECK(!lis_rci_consensus_merge(w,other));
   CHECK(lis_rci_advance_consensus(c,t,w,0,&f)==(mismatch?LIS_RCI_PROTOCOL:LIS_RCI_CALLBACK));
   CHECK(lis_rci_get_info(c,&info)==0&&info.callback_code==73);
   CHECK(lis_rci_get_stats(c,&stats)==0&&stats.error_rank==2&&stats.completed[3]==1);
   CHECK(f.request.token==0);
  }
 }
 lis_rci_destroy(c);return 0;
}
/* Verify terminal states with all three drivers, independently of convergence
 * on the normal test matrices. */
static int edge(int method,int mode,int scenario)
{
 lis_rci_options o={3,scenario==2?1:20,3,1e-12,1e-14,(lis_rci_method)method};
 lis_rci_context *c=lis_rci_create(&o);CHECK(c);
 double b[3]={2,4,8},x0[3]={0,0,0};
 if(scenario==0)b[0]=b[1]=b[2]=0;
 if(scenario==1)x0[0]=x0[1]=x0[2]=1;
 lis_rci_request q;lis_rci_request_ex ex;lis_rci_request_full f;
 int s=mode==0?lis_rci_start(c,b,x0,&q):mode==1?lis_rci_start_ex(c,b,x0,&ex):lis_rci_start_full(c,b,x0,&f);
 int requests=0;
 while(s>0){
  CHECK(++requests<2000&&lis_rci_get_pending(c,&f)==s);
  lis_rci_request_ex *p=&f.request;
  long double scalar=0;
  if(s==LIS_RCI_MATVEC){
   for(size_t i=0;i<3;i++)p->out[i]=scenario==3?0:(double)(2u<<i)*p->in[i];
  }else if(s==LIS_RCI_PRECON){
   for(size_t i=0;i<3;i++)p->out[i]=p->in[i];
  }else if(s==LIS_RCI_DOT||s==LIS_RCI_NORM){
   for(size_t i=0;i<3;i++)scalar+=(long double)p->in[i]*(p->in2?p->in2[i]:p->in[i]);
   if(s==LIS_RCI_NORM)scalar=sqrtl(scalar);
  }else{
   for(size_t i=0;i<3;i++)
    p->out[i]=(double)(s==LIS_RCI_DIVIDE?(long double)p->in[i]/f.alpha:
       f.alpha*p->in[i]+(p->in2?f.beta*p->in2[i]:0));
  }
  uint64_t token=p->token;
  if(mode==0)s=lis_rci_advance(c,0,&q);
  else if(mode==1)s=lis_rci_advance_ex(c,token,0,scalar,&ex);
  else s=lis_rci_advance_full(c,token,0,scalar,&f);
 }
 int expected=scenario<2?LIS_RCI_DONE:scenario==2?LIS_RCI_MAXITER:LIS_RCI_BREAKDOWN;
 lis_rci_info info;lis_rci_stats stats;
 CHECK(s==expected&&!lis_rci_get_info(c,&info)&&info.status==s);
 CHECK(!lis_rci_get_stats(c,&stats));
 CHECK(lis_rci_get_pending(c,&f)==LIS_RCI_STATE&&!f.request.token);
 for(size_t i=0;i<6;i++)CHECK(stats.issued[i]==stats.completed[i]);
 if(scenario<2)CHECK(info.iterations==0&&info.residual==0&&info.matvecs==1&&info.precons==0);
 if(scenario==2)CHECK(info.iterations==1&&info.residual>0);
 if(scenario==3)CHECK(info.residual>0);
 lis_rci_destroy(c);return 0;
}
int main(void)
{
 for(int m=1;m<=4;m++)for(int mode=0;mode<3;mode++)CHECK(!exercise(m,mode));
 for(int m=1;m<=4;m++)for(int mode=0;mode<3;mode++)for(int e=0;e<4;e++)CHECK(!edge(m,mode,e));
 puts("RCI terminal states 48 method/mode/scenario cases PASS");
 puts("RCI lifecycle 12 method/mode combinations and checked consensus PASS");return 0;
}
