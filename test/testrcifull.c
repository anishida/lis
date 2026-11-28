#include "lis_rci.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"full RCI line %d: %s\n",__LINE__,#c);return 1;}}while(0)
enum {N=7};
static int run(int method,int fault,int fault_type,int *requests)
{
 lis_rci_options o={N,200,4,1e-11,1e-14,(lis_rci_method)method};
 lis_rci_context *c=lis_rci_create(&o);CHECK(c);
 double b[N],x0[N];unsigned seen=0;uint64_t last=0;
 for(int rep=1;rep<=2;rep++){
  for(size_t i=0;i<N;i++){
   b[i]=(4+i*.2)*(i+1)*rep+(i? (method==LIS_RCI_CG?1.0:-1.0)*(double)i*rep:0)+(i+1<N?(i+2)*rep:0);
   x0[i]=rep==2?.25:0;
  }
  lis_rci_request_full f;int s=lis_rci_start_full(c,b,x0,&f),count=0;
  while(s>0){
   lis_rci_request_ex *q=&f.request;
   CHECK(++count<10000&&q->operation==s&&q->n==N&&q->in&&q->token>last);
   seen|=1u<<s;
   if(rep==1&&count==fault){
    last=q->token;
    if(fault_type==3){
     CHECK(lis_rci_cancel(c,last-1)==LIS_RCI_PROTOCOL);
     CHECK(lis_rci_cancel(c,last)==LIS_RCI_CANCELLED);
     lis_rci_stats st;CHECK(!lis_rci_get_stats(c,&st)&&st.cancelled);
     CHECK(lis_rci_get_pending(c,&f)==LIS_RCI_STATE&&f.request.token==0);
    }
    else if(fault_type==1){CHECK(lis_rci_advance_full(c,q->token,91,0,&f)==LIS_RCI_CALLBACK);}
    else if(q->out){q->out[0]=0;CHECK(lis_rci_advance_full(c,q->token,0,0,&f)==LIS_RCI_NONFINITE);}
    else {CHECK(lis_rci_advance_full(c,q->token,0,NAN,&f)==LIS_RCI_NONFINITE);}
    lis_rci_info info;CHECK(!lis_rci_get_info(c,&info)&&info.status==(fault_type==3?LIS_RCI_CANCELLED:fault_type==1?LIS_RCI_CALLBACK:LIS_RCI_NONFINITE));
    if(fault_type==1)CHECK(info.callback_code==91);
    /* Reuse after the injected error: the next start must issue a newer token. */
    CHECK(lis_rci_start_full(c,b,x0,&f)==LIS_RCI_NORM&&f.request.token>last);
    CHECK(lis_rci_advance_full(c,f.request.token,42,0,&f)==LIS_RCI_CALLBACK);
    lis_rci_destroy(c);return 0;
   }
   lis_rci_request_full saved=f;
   CHECK(lis_rci_advance_full(c,last,0,0,&f)==LIS_RCI_PROTOCOL);
   CHECK(f.request.token==saved.request.token&&f.request.in==saved.request.in&&f.request.out==saved.request.out&&f.alpha==saved.alpha&&f.beta==saved.beta);
   lis_rci_request_ex ex;lis_rci_request legacy;
   CHECK(lis_rci_advance_ex(c,q->token,0,0,&ex)==LIS_RCI_STATE);
   CHECK(lis_rci_advance(c,0,&legacy)==LIS_RCI_STATE);
   long double v=0;
   if(s==LIS_RCI_DOT||s==LIS_RCI_NORM){
    CHECK(!q->out);
    for(size_t i=0;i<N;i++)v+=(long double)q->in[i]*(q->in2?q->in2[i]:q->in[i]);
    if(s==LIS_RCI_NORM)v=sqrtl(v);
   }else{
    CHECK(q->out&&q->out!=q->in&&q->out!=q->in2);
    for(size_t i=0;i<N;i++){
     if(s==LIS_RCI_MATVEC)q->out[i]=(4+i*.2)*q->in[i]+(i?(method==LIS_RCI_CG?1:-1)*q->in[i-1]:0)+(i+1<N?q->in[i+1]:0);
     else if(s==LIS_RCI_PRECON)q->out[i]=q->in[i]/((4+i*.2)*(method==LIS_RCI_FGMRES?1+.1*((i+count)%3):1));
     else if(s==LIS_RCI_AXPBY)q->out[i]=(double)(f.alpha*q->in[i]+(q->in2?f.beta*q->in2[i]:0));
     else {CHECK(s==LIS_RCI_DIVIDE&&f.alpha>0&&!q->in2);q->out[i]=(double)((long double)q->in[i]/f.alpha);}
    }
   }
   last=q->token;s=lis_rci_advance_full(c,last,0,v,&f);
  }
  CHECK(!fault&&s==LIS_RCI_DONE&&f.request.token==0&&f.request.out==NULL);
  unsigned expected=method>=LIS_RCI_GMRES?126u:62u;
  CHECK(seen==expected);
  const double *x=lis_rci_solution(c);
  long double rn=0,bn=0;
  for(size_t i=0;i<N;i++){
   CHECK(fabs(x[i]-(i+1)*rep)<1e-8);
   long double r=b[i]-(4+i*.2)*x[i]-(i?(method==LIS_RCI_CG?1:-1)*x[i-1]:0)-(i+1<N?x[i+1]:0);
   rn+=r*r;bn+=(long double)b[i]*b[i];
  }
  CHECK(sqrtl(rn/bn)<1.05e-11L);
  if(rep==1)*requests=count;
 }
 lis_rci_destroy(c);return 0;
}
int main(void)
{
 int total=0;
 for(int m=1;m<=4;m++){
  int n=0;CHECK(!run(m,0,0,&n));CHECK(n>0);
  for(int i=1;i<=n;i++){CHECK(!run(m,i,1,&n));CHECK(!run(m,i,2,&n));CHECK(!run(m,i,3,&n));total+=3;}
 }
 printf("RCI full algebra: 8 solves and %d fault/cancellation injections PASS\n",total);
 return 0;
}
