/* Fixed-layout runtime build ABI query for embedding applications. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif
#include <string.h>
#include "lislib.h"
#include "lis_abi.h"

int lis_abi_get_info(LIS_ABI_INFO *info, size_t info_size)
{
    LIS_ABI_INFO result;
    size_t copy_size;

    /* The v1 prefix must never change, even if later versions grow. */
    typedef char lis_abi_v1_prefix_size_check[
        (offsetof(LIS_ABI_INFO,pointer_bytes)==LIS_ABI_INFO_V1_SIZE)?1:-1];
    (void)sizeof(lis_abi_v1_prefix_size_check);
    if(info==NULL || info_size<LIS_ABI_INFO_V1_SIZE)
        return LIS_ABI_ERR_ARGUMENT;

    /* Prefix assertions also guard future ABI descriptor extensions. */
    typedef char lis_abi_v2_prefix_size_check[
        (offsetof(LIS_ABI_INFO,comm_bytes)==LIS_ABI_INFO_V2_SIZE)?1:-1];
    typedef char lis_abi_v3_prefix_size_check[
        (offsetof(LIS_ABI_INFO,esolver_struct_bytes)+sizeof(uint32_t)
         ==LIS_ABI_INFO_V3_SIZE)?1:-1];
    (void)sizeof(lis_abi_v2_prefix_size_check);
    (void)sizeof(lis_abi_v3_prefix_size_check);
    memset(&result,0,sizeof(result));
    lis_abi_get_client_info(&result);
    copy_size=info_size<sizeof(result) ? info_size : sizeof(result);
    memcpy(info,&result,copy_size);
    return LIS_ABI_SUCCESS;
}

int lis_abi_is_compatible(const LIS_ABI_INFO *info,
    uint32_t lis_int_bytes, uint32_t lis_scalar_bytes,
    uint32_t lis_real_bytes, uint32_t scalar_kind)
{
    if(info==NULL || info->struct_size<LIS_ABI_INFO_V1_SIZE ||
       info->version==0)
        return 0;
    return info->lis_int_bytes==lis_int_bytes &&
           info->lis_scalar_bytes==lis_scalar_bytes &&
           info->lis_real_bytes==lis_real_bytes &&
           info->scalar_kind==scalar_kind;
}

int lis_abi_is_layout_compatible(const LIS_ABI_INFO *linked,size_t linked_size,
    const LIS_ABI_INFO *client,size_t client_size,uint32_t flag_mask)
{
    if(!linked || !client ||
       linked_size<LIS_ABI_INFO_V3_SIZE || client_size<LIS_ABI_INFO_V3_SIZE ||
       linked->version<3u || client->version<3u ||
       linked->struct_size<LIS_ABI_INFO_V3_SIZE ||
       client->struct_size<LIS_ABI_INFO_V3_SIZE ||
       linked->reserved || client->reserved)
        return 0;
    return lis_abi_is_compatible(linked,client->lis_int_bytes,
             client->lis_scalar_bytes,client->lis_real_bytes,client->scalar_kind) &&
           (linked->flags & flag_mask)==(client->flags & flag_mask) &&
           linked->pointer_bytes==client->pointer_bytes &&
           linked->endianness==client->endianness &&
           (linked->endianness==LIS_ABI_ENDIAN_LITTLE ||
            linked->endianness==LIS_ABI_ENDIAN_BIG) &&
           linked->comm_bytes==client->comm_bytes &&
           linked->vector_struct_bytes==client->vector_struct_bytes &&
           linked->matrix_struct_bytes==client->matrix_struct_bytes &&
           linked->precon_struct_bytes==client->precon_struct_bytes &&
           linked->solver_struct_bytes==client->solver_struct_bytes &&
           linked->esolver_struct_bytes==client->esolver_struct_bytes;
}
