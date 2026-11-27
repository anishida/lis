/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
/* Adapted from the archived independent model tests, not original RC3 code.
 * Exercise the actual public LIS consensus functions against an independent oracle. */
#include "lis_rci.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
typedef struct {uint64_t token,bits;int32_t error;uint32_t rank;} input;
static uint64_t s=UINT64_C(0xd1b54a32d192ed03);
static uint64_t rng(void){s^=s>>12;s^=s<<25;s^=s>>27;return s*UINT64_C(2685821657736338717);}
static const uint64_t edges[]={0,1,UINT64_MAX,UINT64_MAX-1,UINT64_C(1)<<63,(UINT64_C(1)<<63)-1,(UINT64_C(1)<<53),(UINT64_C(1)<<53)+1,UINT64_C(0x7ff0000000000000),UINT64_C(0xfff0000000000000),UINT64_C(0x7ff8000000000001),UINT64_C(0x7ff8000000000002),UINT64_C(0x0000000000000001),UINT64_C(0x8000000000000001)};
static const int32_t errs[]={0,1,-1,INT32_MIN,INT32_MAX,42,-42};
static unsigned long checks;
#define C(e) do{if(!(e)){fprintf(stderr,"AUDIT FAIL scenario=%lu line=%d: %s\n",checks,__LINE__,#e);exit(1);}}while(0)
static void compare(const uint64_t w[5],const input *a,size_t n)
{
 lis_rci_consensus_result got;int tm=0,bm=0;int32_t err=0;uint32_t rank=UINT32_MAX;
 for(size_t i=0;i<n;i++){tm|=a[i].token!=a[0].token;bm|=a[i].bits!=a[0].bits;if(a[i].error&&(!err||a[i].rank<rank)){err=a[i].error;rank=a[i].rank;}}
 C(!lis_rci_consensus_unpack(w,&got));C(got.token_mismatch==tm&&got.tag_mismatch==bm&&got.error==err&&got.error_rank==rank);
}
static void test(input *a,size_t n)
{
 uint64_t f[5]={0},r[5]={0},t[64][5];
 for(size_t i=0;i<n;i++){C(!lis_rci_consensus_pack(a[i].token,a[i].bits,a[i].rank,a[i].error,t[i]));C(!lis_rci_consensus_merge(f,t[i]));}
 for(size_t i=n;i>0;i--)C(!lis_rci_consensus_merge(r,t[i-1]));
 for(size_t stride=1;stride<n;stride*=2)for(size_t i=0;i+stride<n;i+=2*stride)C(!lis_rci_consensus_merge(t[i],t[i+stride]));
 compare(f,a,n);compare(r,a,n);compare(t[0],a,n);checks++;
}
int main(void)
{
 input a[64];size_t ne=sizeof edges/sizeof *edges,ns=sizeof errs/sizeof *errs;
 for(size_t i=0;i<ne;i++)for(size_t j=0;j<ne;j++)for(size_t e=0;e<ns;e++)for(size_t f=0;f<ns;f++){a[0]=(input){edges[i],edges[j],errs[e],0};a[1]=(input){edges[j],edges[i],errs[f],1};test(a,2);}
 for(size_t i=0;i<ns;i++){a[0]=(input){0,0,errs[i],UINT32_MAX};test(a,1);}
 for(unsigned k=0;k<100000;k++){
  size_t n=1+rng()%64;uint64_t t=rng(),b=rng();
  for(size_t i=0;i<n;i++)a[i]=(input){k%4==0?t:rng(),k%3==0?b:rng(),k%5==0?0:errs[rng()%ns],(uint32_t)i};
  if(k%7==0){size_t i=rng()%n;a[i].token=t^1;a[i].bits=b^1;}
  for(size_t i=n;i>1;i--){size_t j=rng()%i;input tmp=a[i-1];a[i-1]=a[j];a[j]=tmp;}
  test(a,n);
 }
 uint64_t zero[5]={0},w[5];lis_rci_consensus_result result;
 C(lis_rci_consensus_unpack(zero,&result)==LIS_RCI_BADARG);
 C(lis_rci_consensus_pack(1,2,3,0,NULL)==LIS_RCI_BADARG);
 C(lis_rci_consensus_merge(NULL,zero)==LIS_RCI_BADARG);
 C(!lis_rci_consensus_pack(1,2,3,0,w));w[4]=UINT64_C(1)<<32;
 C(lis_rci_consensus_unpack(w,&result)==LIS_RCI_BADARG);
 printf("AUDIT ACTUAL CONSENSUS %lu scenarios, forward/reverse/tree + malformed records PASS\n",checks);
 return 0;
}
