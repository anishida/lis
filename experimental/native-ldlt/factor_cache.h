#ifndef FACTOR_CACHE_H
#define FACTOR_CACHE_H
#include "lis_adapter.h"
/* Caller supplies a 32-byte key binding the matrix, ordering, factor options
   and algorithm/build compatibility. CRC detects accidental corruption,
   not malicious modification. POSIX regular files, IEEE binary64. */
enum { P57_IO=5, P57_CORRUPT=6, P57_MISS=7, P57_LIMIT=8 };
int p57_save(const p56_context *ctx, const char *path, const uint8_t key[32]);
/* max_bytes and max_n are mandatory nonzero resource limits. On any failure
   *slot is unchanged. Caller serializes load against applies/rebuild/destroy. */
int p57_load(p56_context **slot, const char *path, const uint8_t key[32],
             uint64_t max_bytes, uint32_t max_n);
#endif
