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
    const uint16_t byte_order_probe=1;
    size_t copy_size;

    /* The v1 prefix must never change, even if later versions grow. */
    typedef char lis_abi_v1_prefix_size_check[
        (offsetof(LIS_ABI_INFO,pointer_bytes)==LIS_ABI_INFO_V1_SIZE)?1:-1];
    (void)sizeof(lis_abi_v1_prefix_size_check);
    if(info==NULL || info_size<LIS_ABI_INFO_V1_SIZE)
        return LIS_ABI_ERR_ARGUMENT;

    memset(&result,0,sizeof(result));
    result.version=LIS_ABI_INFO_VERSION;
    result.struct_size=(uint32_t)sizeof(result);
    result.lis_int_bytes=(uint32_t)sizeof(LIS_INT);
    result.lis_scalar_bytes=(uint32_t)sizeof(LIS_SCALAR);
    result.lis_real_bytes=(uint32_t)sizeof(LIS_REAL);
#ifdef _COMPLEX
    result.scalar_kind=LIS_ABI_SCALAR_COMPLEX;
#else
    result.scalar_kind=LIS_ABI_SCALAR_REAL;
#endif
#ifdef USE_MPI
    result.flags|=LIS_ABI_FLAG_MPI;
#endif
#ifdef _OPENMP
    result.flags|=LIS_ABI_FLAG_OPENMP;
#endif
#ifdef USE_QUAD_PRECISION
    result.flags|=LIS_ABI_FLAG_QUAD_SOLVER;
#endif
#ifdef _LONG__DOUBLE
    result.flags|=LIS_ABI_FLAG_LONG_DOUBLE;
#endif
#if defined(_COMPLEX) && defined(_LONG__DOUBLE) && defined(HAVE_QUADMATH_H)
    result.flags|=LIS_ABI_FLAG_FLOAT128;
#endif
    result.pointer_bytes=(uint32_t)sizeof(void *);
    result.endianness=*((const unsigned char *)&byte_order_probe) ?
        LIS_ABI_ENDIAN_LITTLE : LIS_ABI_ENDIAN_BIG;
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
