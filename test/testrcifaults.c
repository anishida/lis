/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#include "lis_rci.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
static unsigned checks,failed;
#define CHECK(x) do{checks++;if(!(x)){fprintf(stderr,"FAIL line %d %s\n",__LINE__,#x);failed++;}}while(0)
static void fill(int st,lis_rci_request q){for(size_t i=0;i<q.n;i++)q.out[i]=st==LIS_RCI_PRECON?q.in[i]/4.: (4.+(i%7)*.1)*q.in[i]-(i?q.in[i-1]:0)-(i+1<q.n?q.in[i+1]*.9:0);}
/* CG must use symmetric A. */
static void fill_method(lis_rci_method m,int st,lis_rci_request q){if(m!=LIS_RCI_CG||st!=LIS_RCI_MATVEC){fill(st,q);return;}for(size_t i=0;i<q.n;i++)q.out[i]=(4.+(i%7)*.1)*q.in[i]-(i?q.in[i-1]:0)-(i+1<q.n?q.in[i+1]:0);}
static int run(lis_rci_context*c,lis_rci_method m,double*b,int at,int type,int*count){lis_rci_request q;int s=lis_rci_start(c,b,NULL,&q),k=0;while(s>0&&k<2000){k++;if(k==at&&type==1){s=lis_rci_advance(c,-123,&q);break;}if(k==at&&type==2){for(size_t i=0;i+1<q.n;i++)q.out[i]=0;s=lis_rci_advance(c,0,&q);break;}fill_method(m,s,q);s=lis_rci_advance(c,0,&q);}*count=k;return s;}
int main(void){double b[29];for(size_t i=0;i<29;i++)b[i]=sin((double)i)+1;for(int m=1;m<=4;m++){lis_rci_options o={29,200,5,1e-10,1e-13,(lis_rci_method)m};lis_rci_context*c=lis_rci_create(&o);CHECK(c!=NULL);if(!c)continue;int count,k;CHECK(run(c,o.method,b,0,0,&count)==LIS_RCI_DONE);CHECK(count>3);for(int i=1;i<=count;i++){CHECK(run(c,o.method,b,i,1,&k)==LIS_RCI_CALLBACK);CHECK(k==i);lis_rci_info inf;lis_rci_get_info(c,&inf);CHECK(inf.callback_code==-123);CHECK(run(c,o.method,b,i,2,&k)==LIS_RCI_NONFINITE);CHECK(k==i);CHECK(run(c,o.method,b,0,0,&k)==LIS_RCI_DONE);CHECK(k==count);}printf("method=%d injection_positions=%d\n",m,count);lis_rci_destroy(c);}printf("FAULTS checks=%u failures=%u\n",checks,failed);return failed?1:0;}
