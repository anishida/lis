/* Experimental application-side reference driver; not an installed LIS API. */
#ifndef LIS_ANCHOR_PILOT_H
#define LIS_ANCHOR_PILOT_H
#include "lis.h"

typedef struct lis_anchor_pilot lis_anchor_pilot;
typedef struct {
    LIS_INT raw_status, pilot_iterations, continuation_iterations, fresh_iterations;
    int accepted, continued, fallback, factorizations;
    double relative_residual, ratio, setup_seconds;
} lis_anchor_pilot_report;

/* All retention modes use one FGMRES solve, with a projected checkpoint
 * at iteration 20 and the configured restart length of 60.
 * Serial real-double, unscaled assembled CSR, one active thread.
 * P is borrowed and must stay alive and unchanged until destroy. A/b/x
 * share its DOF ordering. The driver owns the solver and factors, starts
 * from zero for each RHS, and never modifies x unless physically accepted.
 * Native numeric statuses are recorded separately from certified acceptance.
 * API/argument failures are not retried. Recreate after a hard solve error. */
LIS_INT lis_anchor_pilot_create(LIS_MATRIX P, lis_anchor_pilot **out);
/* Memory-constrained variant: retire the anchor before a fresh fallback,
 * and rebuild it lazily on the next solve. Never holds two factors.
 * This trades extra anchor setups for lower peak memory. */
LIS_INT lis_anchor_pilot_create_low_memory(LIS_MATRIX P, lis_anchor_pilot **out);
/* Optional rolling mode: retire the old factor before fallback and retain
 * the fresh factor only after physical acceptance and output commit.
 * Keep both matrices alive during solve. After EVERY return, query the
 * borrowed current anchor; only that matrix must remain alive and immutable
 * until the next solve or destroy. On failure the anchor identity is unchanged
 * (its factor may be absent). The context never owns any matrix. */
LIS_INT lis_anchor_pilot_create_rolling(LIS_MATRIX P, lis_anchor_pilot **out);
LIS_MATRIX lis_anchor_pilot_current_anchor(const lis_anchor_pilot *ctx);
LIS_INT lis_anchor_pilot_solve(lis_anchor_pilot *ctx, LIS_MATRIX A,
    LIS_VECTOR b, LIS_VECTOR x, lis_anchor_pilot_report *report);
/* Matched fresh-factorization reference, same options and certification. */
LIS_INT lis_anchor_fresh_solve(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,
    lis_anchor_pilot_report *report);
void lis_anchor_pilot_destroy(lis_anchor_pilot *ctx);
#endif
