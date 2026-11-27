/* Only the fixed-width ABI query is used with baseline headers.
 * Never pass an old-layout LIS object into the new library. */
#include "lis.h"
#include "lis_abi.h"
#include <stdio.h>
#include <stddef.h>
int main(void)
{
    LIS_ABI_INFO local,linked;
    int same;
    lis_abi_get_client_info(&local);
    if(lis_abi_get_info(&linked,sizeof linked)) return 1;
    same=lis_abi_is_layout_compatible(&linked,sizeof linked,&local,sizeof local,LIS_ABI_LAYOUT_FLAGS);
    printf("client_solver=%u linked_solver=%u client_precon=%u linked_precon=%u compatible=%d\n",
        local.solver_struct_bytes,linked.solver_struct_bytes,
        local.precon_struct_bytes,linked.precon_struct_bytes,same);
    printf("offsets=%zu,%zu,%zu,%zu,%zu\n",
        offsetof(struct LIS_SOLVER_STRUCT,A),offsetof(struct LIS_SOLVER_STRUCT,precon),
        offsetof(struct LIS_SOLVER_STRUCT,options),offsetof(struct LIS_PRECON_STRUCT,L),
        offsetof(struct LIS_PRECON_STRUCT,user_data));
#ifdef BASELINE_CLIENT
    if(same || local.solver_struct_bytes>=linked.solver_struct_bytes ||
       local.precon_struct_bytes>=linked.precon_struct_bytes) return 1;
    puts("BASELINE_CLIENT_REJECTED=PASS");
#else
    if(!same) return 1;
    puts("MATCHING_CLIENT_ACCEPTED=PASS");
#endif
    return 0;
}
