/* ABI-stable runtime description of the linked LIS library.
 *
 * This header intentionally does NOT include lis.h.  A client may call
 * lis_abi_get_info() before making any calls that use LIS_INT/LIS_SCALAR,
 * even if its generated lis_config.h differs from the linked library.
 */
#ifndef LIS_ABI_H
#define LIS_ABI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * The first 32 bytes form the immutable v1 prefix. A future library may
 * append fields without preventing clients compiled with older headers from
 * querying that prefix.
 */
#define LIS_ABI_INFO_VERSION 3u
#define LIS_ABI_INFO_V1_SIZE UINT32_C(32)
#define LIS_ABI_INFO_V2_SIZE UINT32_C(40)
#define LIS_ABI_INFO_V3_SIZE UINT32_C(64)
#define LIS_ABI_SUCCESS 0
#define LIS_ABI_ERR_ARGUMENT 1
#define LIS_ABI_ENDIAN_UNKNOWN 0u
#define LIS_ABI_ENDIAN_LITTLE 1u
#define LIS_ABI_ENDIAN_BIG 2u
#define LIS_ABI_SCALAR_REAL 0u
#define LIS_ABI_SCALAR_COMPLEX 1u

#define LIS_ABI_FLAG_MPI          UINT32_C(0x00000001)
#define LIS_ABI_FLAG_OPENMP       UINT32_C(0x00000002)
#define LIS_ABI_FLAG_QUAD_SOLVER  UINT32_C(0x00000004)
#define LIS_ABI_FLAG_LONG_DOUBLE  UINT32_C(0x00000008)
#define LIS_ABI_FLAG_FLOAT128     UINT32_C(0x00000010)

typedef struct LIS_ABI_INFO_STRUCT
{
    uint32_t version;
    uint32_t struct_size;
    uint32_t lis_int_bytes;
    uint32_t lis_scalar_bytes;
    uint32_t lis_real_bytes;
    uint32_t scalar_kind;
    uint32_t flags;
    uint32_t reserved; /* must be zero */
    /* Added in v2; unavailable to v1 callers passing 32 bytes. */
    uint32_t pointer_bytes;
    uint32_t endianness;
    /* Added in v3. Public struct sizes, not an offset/layout fingerprint. */
    uint32_t comm_bytes;
    uint32_t vector_struct_bytes;
    uint32_t matrix_struct_bytes;
    uint32_t precon_struct_bytes;
    uint32_t solver_struct_bytes;
    uint32_t esolver_struct_bytes;
} LIS_ABI_INFO;

/*
 * Returns LIS_ABI_SUCCESS or LIS_ABI_ERR_ARGUMENT; this independent API
 * deliberately does not depend on LIS_INT or LIS_ERR_* definitions.
 *
 * At least LIS_ABI_INFO_V1_SIZE bytes (32) are required. On success, writes
 * min(info_size, linked_library_struct_size) bytes. The leading 32 bytes are
 * stable, with version and struct_size describing the library's full result.
 * A v1 client can therefore query a v2 library with its 32-byte buffer.
 * Extra caller bytes after the copied prefix are NOT modified.
 * On invalid arguments, the caller buffer is not modified.
 *
 * This reports the linked library's ABI, it does not translate layouts.
 * The library returns local build capabilities, not MPI communicator state.
 * Flags: QUAD_SOLVER means USE_QUAD_PRECISION; FLOAT128 indicates the
 * complex+long-double+quadmath build combination. Flags are informational.
 * Neither this query nor the compatibility helper validates every platform
 * calling convention, Fortran ABI, or external MPI configuration.
 */
int lis_abi_get_info(LIS_ABI_INFO *info, size_t info_size);

/* Compare the essential scalar ABI; nonzero means compatible.
 * The client supplies its OWN expected widths and scalar kind.
 * This helper does not inspect feature flags, pointer ABI or endianness.
 */
int lis_abi_is_compatible(const LIS_ABI_INFO *info,
    uint32_t lis_int_bytes, uint32_t lis_scalar_bytes,
    uint32_t lis_real_bytes, uint32_t scalar_kind);

/* Compare v3 client and linked-library descriptions. Supply the actual
 * buffer lengths, not merely the struct_size reported by the library.
 * Returns zero if either buffer lacks v3 data. Exact known layout fields
 * and caller-selected flags must match. This is a conservative screening
 * check, not a proof of identical offsets, MPI implementation or calling
 * convention. The OpenMP flag may be omitted from flag_mask.
 *
 * External client example (matching installed include and lib paths):
 *   #define LIS_USE_INSTALLED_CONFIG
 *   #include <lis.h>
 *   LIS_ABI_INFO local, linked;
 *   lis_abi_get_client_info(&local);
 *   if(lis_abi_get_info(&linked,sizeof(linked)) ||
 *      !lis_abi_is_layout_compatible(&linked,sizeof(linked),
 *             &local,sizeof(local),LIS_ABI_LAYOUT_FLAGS))
 *       ... reject before using LIS objects ...
 */
#define LIS_ABI_LAYOUT_FLAGS (LIS_ABI_FLAG_MPI | LIS_ABI_FLAG_QUAD_SOLVER | \
                              LIS_ABI_FLAG_LONG_DOUBLE | LIS_ABI_FLAG_FLOAT128)
int lis_abi_is_layout_compatible(const LIS_ABI_INFO *linked,size_t linked_size,
    const LIS_ABI_INFO *client,size_t client_size,uint32_t flag_mask);

#ifdef __cplusplus
}
#endif
#endif
