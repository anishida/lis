#include "lis_rci.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"extended RCI line %d: %s\n",__LINE__,#c);return 1;}}while(0)
static long double scalar(const lis_rci_request_ex *q)
{
 long double s=0;
 for(size_t i=0;i<q->n;i++)s+=(long double)q->in[i]*(q->in2?q->in2[i]:q->in[i]);
 return q->operation==LIS_RCI_NORM?sqrtl(s):s;
}
static int test_solver(int method)
{
 lis_rci_options o={5,100,5,1e-12,1e-14,(lis_rci_method)method};
 lis_rci_context *c=lis_rci_create(&o);CHECK(c);
 double a[5][5]={{4,1,0,0,0},{1,5,1,0,0},{0,1,6,1,0},{0,0,1,7,1},{0,0,0,1,8}};
 if(method!=LIS_RCI_CG)a[1][0]=-1;
 uint64_t last=0;double b[5],x0[5]={0};lis_rci_request_ex q;
 for(int rep=1;rep<=2;rep++){
  for(size_t i=0;i<5;i++){b[i]=0;for(size_t j=0;j<5;j++)b[i]+=a[i][j]*(j+1)*rep;}
  int st=lis_rci_start_ex(c,b,x0,&q),calls=0,seen=0;
  while(st>0){
   CHECK(++calls<5000&&q.operation==st&&q.n==5&&q.in&&q.token>last);
   seen|=1<<st;
   lis_rci_request_ex save=q;
   CHECK(lis_rci_advance_ex(c,last,0,0,&q)==LIS_RCI_PROTOCOL);
   CHECK(q.token==save.token&&q.in==save.in&&q.in2==save.in2&&q.out==save.out);
   lis_rci_request legacy;
   CHECK(lis_rci_advance(c,0,&legacy)==LIS_RCI_STATE);
   long double value=0;
   if(st==LIS_RCI_MATVEC||st==LIS_RCI_PRECON){
    CHECK(q.out&&!q.in2);
    for(size_t i=0;i<5;i++){
     if(st==LIS_RCI_MATVEC){q.out[i]=0;for(size_t j=0;j<5;j++)q.out[i]+=a[i][j]*q.in[j];}
     else q.out[i]=q.in[i]/(a[i][i]*(method==LIS_RCI_FGMRES?1+.1*((i+calls)%3):1));
    }
   }else{
    CHECK(!q.out&&((st==LIS_RCI_DOT&&q.in2)||(st==LIS_RCI_NORM&&!q.in2)));
    value=scalar(&q);
   }
   last=q.token;st=lis_rci_advance_ex(c,last,0,value,&q);
  }
  CHECK(st==LIS_RCI_DONE&&seen==30&&q.token==0&&q.in==NULL);
  const double *x=lis_rci_solution(c);
  long double r2=0,b2=0;
  for(size_t i=0;i<5;i++){
   CHECK(fabs(x[i]-(i+1)*rep)<1e-9);
   long double rr=b[i];for(size_t j=0;j<5;j++)rr-=a[i][j]*x[j];
   r2+=rr*rr;b2+=(long double)b[i]*b[i];
  }
  CHECK(sqrtl(r2)<=1.05e-12L*sqrtl(b2));
 }
 CHECK(lis_rci_start_ex(c,b,NULL,&q)==LIS_RCI_NORM&&q.token>last);
 CHECK(lis_rci_advance_ex(c,q.token,42,0,&q)==LIS_RCI_CALLBACK);
 lis_rci_info info;CHECK(!lis_rci_get_info(c,&info)&&info.callback_code==42);
 CHECK(lis_rci_start_ex(c,b,NULL,&q)==LIS_RCI_NORM);
 CHECK(lis_rci_advance_ex(c,q.token,0,-1,&q)==LIS_RCI_NONFINITE);
 CHECK(lis_rci_start_ex(c,b,NULL,&q)==LIS_RCI_NORM);
 CHECK(lis_rci_advance_ex(c,q.token,0,NAN,&q)==LIS_RCI_NONFINITE);
 /* Legacy mode still works after an extended-mode failure. */
 lis_rci_request l;CHECK(lis_rci_start(c,b,NULL,&l)==LIS_RCI_MATVEC);
 CHECK(lis_rci_advance_ex(c,1,0,0,&q)==LIS_RCI_STATE);
 CHECK(lis_rci_advance(c,77,&l)==LIS_RCI_CALLBACK);
 lis_rci_destroy(c);return 0;
}
static uint64_t rng=UINT64_C(123456789);
static uint64_t next(void){rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;return rng;}
static int test_consensus(void)
{
 for(unsigned n=1;n<=64;n++)for(unsigned t=0;t<128;t++){
  uint64_t w[64][5],f[5],rev[5],tag=next(),token=next();
  uint32_t rank=UINT32_MAX;int32_t error=0;int tm=0,nm=0;
  for(unsigned i=0;i<n;i++){
   uint64_t a=token^(i&&t%3==0&&i==n-1?UINT64_C(1):0);
   uint64_t b=tag^(i&&t%5==0&&i==n-1?UINT64_C(1):0);
   int32_t e=(t%7==0&&i==n-1)?INT32_MIN:((t%11==0&&i>0)?(int32_t)i:0);
   uint32_t k=i==n-1?UINT32_MAX:i;
   if(e&&(!error||k<rank)){rank=k;error=e;}
   tm|=a!=token;nm|=b!=tag;
   CHECK(!lis_rci_consensus_pack(a,b,k,e,w[i]));
  }
  memcpy(f,w[0],sizeof f);memcpy(rev,w[n-1],sizeof rev);
  for(unsigned i=1;i<n;i++){CHECK(!lis_rci_consensus_merge(f,w[i]));CHECK(!lis_rci_consensus_merge(rev,w[n-1-i]));}
  CHECK(!memcmp(f,rev,sizeof f));lis_rci_consensus_result result;
  CHECK(!lis_rci_consensus_unpack(f,&result));
  CHECK(result.token_mismatch==tm&&result.tag_mismatch==nm&&result.error==error&&result.error_rank==rank);
 }
 uint64_t a[5],b[5];lis_rci_consensus_result r;
 CHECK(!lis_rci_consensus_pack(UINT64_C(9007199254740992),0,9,5,a));
 CHECK(!lis_rci_consensus_pack(UINT64_C(9007199254740993),UINT64_MAX,2,-7,b));
 CHECK(!lis_rci_consensus_merge(a,b)&&!lis_rci_consensus_unpack(a,&r));
 CHECK(r.token_mismatch&&r.tag_mismatch&&r.error==-7&&r.error_rank==2);
 puts("RCI consensus 8192 scenarios + exact uint64 edge PASS");return 0;
}
int main(void)
{
 for(int m=LIS_RCI_CG;m<=LIS_RCI_FGMRES;m++)CHECK(!test_solver(m));
 CHECK(!test_consensus());
 puts("RCI external reductions, tokens, failure recovery PASS");return 0;
}
