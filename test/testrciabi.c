/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#include <lis_rci_bridge.h>
#include <stdio.h>
int main(void){
 if(LIS_RCI_CHECK_CLIENT_ABI()!=LIS_RCI_DONE)return 1;
 if(lis_rci_bridge_check_abi(sizeof(LIS_INT)==4?8:4,8,8,LIS_ABI_SCALAR_REAL)!=LIS_RCI_BADARG)return 2;
 if(lis_rci_bridge_check_abi(sizeof(LIS_INT),16,8,LIS_ABI_SCALAR_REAL)!=LIS_RCI_BADARG)return 3;
 puts("RCI_CLIENT_ABI 3/3 PASS");return 0;
}
