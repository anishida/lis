/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#include <lis_rci.h>
#include <cmath>
#include <cstdio>
int main(){
 lis_rci_options o{3,10,3,1e-12,1e-15,LIS_RCI_FGMRES};lis_rci_context*c=lis_rci_create(&o);if(!c)return 1;
 double b[]={2.,4.,6.};lis_rci_request q{};int st=lis_rci_start(c,b,nullptr,&q);
 while(st>0){for(size_t i=0;i<q.n;i++)q.out[i]=(st==LIS_RCI_MATVEC?2.:1.)*q.in[i];st=lis_rci_advance(c,0,&q);}
 const double*x=lis_rci_solution(c);int ok=st==LIS_RCI_DONE&&std::abs(x[0]-1)<1e-12&&std::abs(x[1]-2)<1e-12&&std::abs(x[2]-3)<1e-12;lis_rci_destroy(c);
 std::puts(ok?"INSTALLED_CXX_CONSUMER PASS":"INSTALLED_CXX_CONSUMER FAIL");return ok?0:1;
}
