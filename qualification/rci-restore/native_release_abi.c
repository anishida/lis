/* Client compiled against baseline headers; run with both shared libraries. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <stdio.h>
#include <string.h>
int main(void)
{
 LIS_ABI_INFO linked,local;
 memset(&linked,0,sizeof linked);lis_abi_get_client_info(&local);
 if(lis_abi_get_info(&linked,sizeof linked))return 1;
 if(!lis_abi_is_layout_compatible(&linked,sizeof linked,&local,sizeof local,LIS_ABI_LAYOUT_FLAGS))return 2;
 const unsigned char *bytes=(const unsigned char *)&linked;
 for(size_t i=0;i<sizeof linked;i++)printf("%02x",(unsigned)bytes[i]);
 puts("");
 return 0;
}
