#include "lis_abi.h"
#include "lis.h"
#include <stdio.h>
#include <string.h>

int main(void)
{
    LIS_ABI_INFO info,old;
    union { LIS_ABI_INFO aligned; unsigned char bytes[sizeof(LIS_ABI_INFO)+8]; } extended;
    LIS_ABI_INFO short_result;
    int fail=0;
    memset(&info,0x5a,sizeof(info));
    old=info;
    if(lis_abi_get_info(NULL,sizeof(info))==0 ||
       lis_abi_get_info(&info,LIS_ABI_INFO_V1_SIZE-1)==LIS_ABI_SUCCESS ||
       memcmp(&old,&info,sizeof(info))!=0)
        fail=1;
    if(lis_abi_get_info(&info,sizeof(info))!=0 ||
       info.version!=LIS_ABI_INFO_VERSION ||
       info.struct_size!=sizeof(info) ||
       info.lis_int_bytes!=sizeof(LIS_INT) ||
       info.lis_scalar_bytes!=sizeof(LIS_SCALAR) ||
       info.lis_real_bytes!=sizeof(LIS_REAL) ||
       info.reserved!=0 ||
       info.pointer_bytes!=sizeof(void *) ||
       (info.endianness!=LIS_ABI_ENDIAN_LITTLE &&
        info.endianness!=LIS_ABI_ENDIAN_BIG))
        fail=1;
#ifdef _COMPLEX
    if(info.scalar_kind!=LIS_ABI_SCALAR_COMPLEX) fail=1;
#else
    if(info.scalar_kind!=LIS_ABI_SCALAR_REAL) fail=1;
#endif
#ifdef USE_MPI
    if(!(info.flags & LIS_ABI_FLAG_MPI)) fail=1;
#else
    if(info.flags & LIS_ABI_FLAG_MPI) fail=1;
#endif
#ifdef _OPENMP
    if(!(info.flags & LIS_ABI_FLAG_OPENMP)) fail=1;
#else
    if(info.flags & LIS_ABI_FLAG_OPENMP) fail=1;
#endif
#ifdef USE_QUAD_PRECISION
    if(!(info.flags & LIS_ABI_FLAG_QUAD_SOLVER)) fail=1;
#else
    if(info.flags & LIS_ABI_FLAG_QUAD_SOLVER) fail=1;
#endif
#ifdef _LONG__DOUBLE
    if(!(info.flags & LIS_ABI_FLAG_LONG_DOUBLE)) fail=1;
#else
    if(info.flags & LIS_ABI_FLAG_LONG_DOUBLE) fail=1;
#endif
#if defined(_COMPLEX) && defined(_LONG__DOUBLE) && defined(HAVE_QUADMATH_H)
    if(!(info.flags & LIS_ABI_FLAG_FLOAT128)) fail=1;
#else
    if(info.flags & LIS_ABI_FLAG_FLOAT128) fail=1;
#endif
    /* A future-library-style v2 result must fit an old v1 buffer. */
    memset(&short_result,0x5a,sizeof(short_result));
    if(lis_abi_get_info(&short_result,LIS_ABI_INFO_V1_SIZE)!=LIS_ABI_SUCCESS ||
       short_result.version!=LIS_ABI_INFO_VERSION ||
       short_result.struct_size!=sizeof(LIS_ABI_INFO) ||
       short_result.lis_int_bytes!=sizeof(LIS_INT) ||
       short_result.pointer_bytes!=UINT32_C(0x5a5a5a5a) ||
       short_result.endianness!=UINT32_C(0x5a5a5a5a))
        fail=1;
    /* A longer client buffer must have its unknown tail untouched. */
    memset(extended.bytes,0x5a,sizeof(extended.bytes));
    if(lis_abi_get_info(&extended.aligned,sizeof(extended.bytes))!=LIS_ABI_SUCCESS ||
       extended.bytes[sizeof(LIS_ABI_INFO)]!=0x5a ||
       extended.bytes[sizeof(extended.bytes)-1]!=0x5a)
        fail=1;
    if(!lis_abi_is_compatible(&info,
        (uint32_t)sizeof(LIS_INT),(uint32_t)sizeof(LIS_SCALAR),
        (uint32_t)sizeof(LIS_REAL),info.scalar_kind) ||
       lis_abi_is_compatible(&info,
        (uint32_t)sizeof(LIS_INT)+1,(uint32_t)sizeof(LIS_SCALAR),
        (uint32_t)sizeof(LIS_REAL),info.scalar_kind) ||
       lis_abi_is_compatible(NULL,1,1,1,LIS_ABI_SCALAR_REAL))
        fail=1;
    printf("LIS_ABI_INFO=%s\n",fail?"FAIL":"PASS");
    return fail;
}
