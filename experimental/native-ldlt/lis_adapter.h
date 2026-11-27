#ifndef LIS_ADAPTER_H
#define LIS_ADAPTER_H
#include "packed_ldlt.h"
/* Standalone real-double adapter. Not part of the LIS public ABI. */
typedef struct p56_context p56_context;
/* factor_to_original[k] is the original index of factor coordinate k.
   NULL means identity. Copies permutation and factors; no caller-owned
   pointers retained. On failure *slot is unchanged. Initialize *slot=NULL.
   Caller must serialize rebuild/destroy against applies. */
int p56_rebuild(p56_context **slot, uint32_t n, size_t np,
                const p55_panel_input *panels, const uint8_t *kind,
                const double *diag, const double *off,
                const uint32_t *factor_to_original);
void p56_destroy(p56_context *ctx);
uint32_t p56_size(const p56_context *ctx);
size_t p56_block_count(const p56_context *ctx);
/* Scratch must hold n doubles and be disjoint from input/output.
   Input/output may overlap, including in-place. Output unchanged on failure.
   No allocation, shared mutable state or global state during apply.
   Concurrent applies require separate scratch/output and immutable context. */
int p56_apply(const p56_context *ctx, const double *rhs, double *out,
              double *scratch, size_t scratch_count);
#endif
