/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
/* Regression of the observed wide-residual/finite-check audit findings. */
#include "lis_rci.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#define C(e) do{if(!(e)){fprintf(stderr,"RCI RANGE line %d: %s\n",__LINE__,#e);return 1;}}while(0)
/* A legitimate finite wide DOT can exceed the dense double Hessenberg
 * storage range. This input reaches that branch without lying about DOT. */
static int arnoldi_range(void)
{
#if LDBL_MAX_EXP > DBL_MAX_EXP
 lis_rci_options o={2,5,2,1e-12,0,LIS_RCI_GMRES};
 lis_rci_context *c=lis_rci_create(&o);C(c);
 double b[]={1,1};lis_rci_request_ex q;
 int s=lis_rci_start_ex(c,b,NULL,&q);
 unsigned steps=0;
 while(s>0&&steps++<20){
  long double scalar=0;
  if(s==LIS_RCI_NORM)scalar=hypotl(q.in[0],q.in[1]);
  else if(s==LIS_RCI_DOT)scalar=(long double)q.in[0]*q.in2[0]+(long double)q.in[1]*q.in2[1];
  else if(s==LIS_RCI_PRECON){for(size_t i=0;i<2;i++)q.out[i]=q.in[i];}
  else if(s==LIS_RCI_MATVEC){
   /* A = DBL_MAX * [[1,1/4],[1/4,1]], finite for these inputs. */
   for(size_t i=0;i<2;i++){
    long double v=(long double)DBL_MAX*((long double)q.in[i]+.25L*q.in[1-i]);
    C(fabsl(v)<=DBL_MAX);q.out[i]=(double)v;
   }
  }
  else {lis_rci_destroy(c);return 1;}
  s=lis_rci_advance_ex(c,q.token,0,scalar,&q);
 }
 C(s==LIS_RCI_NONFINITE);
 lis_rci_destroy(c);
#endif
 return 0;
}
int main(void)
{
 C(!arnoldi_range());
 lis_rci_options o={2,10,2,1,0,LIS_RCI_CG};
 lis_rci_context *c=lis_rci_create(&o);C(c);
 double b[]={DBL_MAX,DBL_MAX};lis_rci_request q;lis_rci_info info;long double wide;
 C(!lis_rci_get_residual(c,&wide)&&isnan(wide));
 int s=lis_rci_start(c,b,NULL,&q);
#if LDBL_MAX_EXP > DBL_MAX_EXP
 C(s==LIS_RCI_MATVEC);
 for(size_t i=0;i<2;i++)q.out[i]=q.in[i];
 C(lis_rci_advance(c,0,&q)==LIS_RCI_DONE);
 C(!lis_rci_get_info(c,&info)&&isinf(info.residual)&&info.residual>0);
 C(!lis_rci_get_residual(c,&wide)&&isfinite(wide)&&wide>DBL_MAX);
 C(fabsl(wide/DBL_MAX-sqrtl(2.L))<1e-14L);
#else
 C(s==LIS_RCI_NONFINITE);
#endif
 double bad[]={NAN,1};
 C(lis_rci_start(c,bad,NULL,&q)==LIS_RCI_NONFINITE);
 C(!lis_rci_get_residual(c,&wide)&&isnan(wide));
 double normal[]={1,2};
 C(lis_rci_start(c,normal,NULL,&q)==LIS_RCI_MATVEC);
 q.out[0]=0;
 C(lis_rci_advance(c,0,&q)==LIS_RCI_NONFINITE);
 C(!lis_rci_get_residual(c,&wide)&&isnan(wide));
 C(lis_rci_get_residual(NULL,&wide)==LIS_RCI_BADARG);
 C(lis_rci_get_residual(c,NULL)==LIS_RCI_BADARG);
 C(lis_rci_start(c,normal,NULL,&q)==LIS_RCI_MATVEC);
 for(size_t i=0;i<2;i++)q.out[i]=0;
 C(lis_rci_advance(c,0,&q)==LIS_RCI_DONE);
 C(!lis_rci_get_residual(c,&wide)&&fabsl(wide-sqrtl(5.L))<1e-14L);
 lis_rci_destroy(c);
 puts("RCI WIDE RESIDUAL and nonfinite guards PASS");return 0;
}
