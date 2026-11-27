#include "qmd_plan.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"FAIL %d %s\n",__LINE__,#x);exit(1);}}while(0)
static unsigned tests;
static void check(p56_context *f,const size_t *c,const uint32_t *r,const double *a,unsigned n){
 double b[8]={0},x[8],scratch[8];
 for(unsigned j=0;j<n;j++)for(size_t z=c[j];z<c[j+1];z++){unsigned i=r[z];b[i]+=a[z]*(j+1);if(i!=j)b[j]+=a[z]*(i+1);}
 CHECK(!p56_apply(f,b,x,scratch,n));for(unsigned j=0;j<n;j++)CHECK(fabs(x[j]-(j+1))<1e-11);tests++;
}
int main(void){
 size_t c[]={0,4,7,9,10};uint32_t r[]={0,1,2,3,1,2,3,2,3,3};
 double a[]={0,4,0,0,0,0,0,2,0,-3};double b[]={5,1,0.125,0,-2,0.25,0,4,0.5,-3};double zero[10]={0};
 p63_limits ol={1<<20,10000,100000};p62_limits sl={1<<20,10000};p61_limits nl={100,1<<20,10000,10000,0};p71_plan *p=NULL;p56_context *f=NULL,*fresh=NULL;
 for(unsigned width=1;width<=4;width++){
  CHECK(!p71_plan_create(&p,4,10,c,r,&ol,&sl,width,1<<20));
  CHECK(!p71_refactor(p,&f,4,10,c,r,a,&nl,NULL));check(f,c,r,a,4);
  CHECK(!p63q_factorize(&fresh,4,10,c,r,a,&ol,&sl,&nl,width,NULL));check(fresh,c,r,a,4);
  CHECK(!p71_refactor(p,&f,4,10,c,r,b,&nl,NULL));check(f,c,r,b,4);
  CHECK(!p63q_factorize(&fresh,4,10,c,r,b,&ol,&sl,&nl,width,NULL));check(fresh,c,r,b,4);
  p56_context *old=f;p71_plan *op=p;p63_stats st={0};st.numeric.delays=123;
  CHECK(p71_refactor(p,&f,4,10,c,r,zero,&nl,&st)!=0&&f==old&&st.numeric.delays==123);check(f,c,r,b,4);
  uint32_t bad[10];memcpy(bad,r,sizeof r);bad[1]=2;CHECK(p71_refactor(p,&f,4,10,c,bad,b,&nl,NULL)==P55_INVALID&&f==old);check(f,c,r,b,4);
  p61_limits tiny=nl;tiny.max_work_bytes=1;CHECK(p71_refactor(p,&f,4,10,c,r,b,&tiny,NULL)==P60_LIMIT&&f==old);check(f,c,r,b,4);
  CHECK(p71_plan_create(&p,4,10,c,r,&ol,&sl,width,1)==P60_LIMIT&&p==op);tests++;
 }
 p71_plan_free(p);p56_destroy(f);p56_destroy(fresh);printf("QMD plan %u/%u PASS\n",tests,tests);return 0;
}
