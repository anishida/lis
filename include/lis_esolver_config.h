/* Generic typed eigensolver configuration helpers.
 *
 * Setters validate the complete active configuration before publication.
 * Their order matters for dependent options: to leave block Lanczos, first
 * set block_size=1, then change eigensolver type; to enter block Lanczos,
 * set the Lanczos type and adequate subspace before increasing block_size.
 * Reduce mode/block_size/nev constraints before shrinking the subspace.
 *
 * Type/subspace/block/transform settings may invalidate work vectors, while
 * scalar setters preserve work. Previously computed eigenpairs and result
 * arrays should be treated as results of the prior solve, not as results
 * corresponding to a newly configured solver; run a new solve after changes.
 */
#ifndef LIS_ESOLVER_CONFIG_H
#define LIS_ESOLVER_CONFIG_H

#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

LIS_INT lis_esolver_set_type(
    LIS_ESOLVER esolver,
    LIS_INT type);
LIS_INT lis_esolver_get_type(
    LIS_ESOLVER esolver,
    LIS_INT *type);

LIS_INT lis_esolver_set_tolerance(
    LIS_ESOLVER esolver,
    LIS_REAL tolerance);
LIS_INT lis_esolver_get_tolerance(
    LIS_ESOLVER esolver,
    LIS_REAL *tolerance);

LIS_INT lis_esolver_set_maxiter(
    LIS_ESOLVER esolver,
    LIS_INT maxiter);
LIS_INT lis_esolver_get_maxiter(
    LIS_ESOLVER esolver,
    LIS_INT *maxiter);

LIS_INT lis_esolver_set_subspace_size(
    LIS_ESOLVER esolver,
    LIS_INT subspace_size);
LIS_INT lis_esolver_get_subspace_size(
    LIS_ESOLVER esolver,
    LIS_INT *subspace_size);

/*
 * Requested eigenpair count (-nev): 0 retains the historical default of
 * using the configured subspace size; a positive count must fit the current
 * subspace. Set subspace_size first when increasing both values.
 *
 * The requested count is stored separately from esolver->nev, which some
 * eigensolvers overwrite with their result count. The getter remains stable
 * across solves, including retries and executions that return fewer pairs.
 * A request configured using -nev is also preserved across solves. Historical
 * text parsing permits -nev before -ss and validates their final relationship
 * when solving; this typed setter validates against the active subspace now.
 * Legacy callers that directly modify esolver->nev remain supported if they
 * have not activated typed/text requested-nev configuration.
 */
LIS_INT lis_esolver_set_requested_nev(
    LIS_ESOLVER esolver,
    LIS_INT requested_nev);
LIS_INT lis_esolver_get_requested_nev(
    LIS_ESOLVER esolver,
    LIS_INT *requested_nev);

LIS_INT lis_esolver_set_mode(
    LIS_ESOLVER esolver,
    LIS_INT mode);
LIS_INT lis_esolver_get_mode(
    LIS_ESOLVER esolver,
    LIS_INT *mode);

LIS_INT lis_esolver_set_block_size(
    LIS_ESOLVER esolver,
    LIS_INT block_size);
LIS_INT lis_esolver_get_block_size(
    LIS_ESOLVER esolver,
    LIS_INT *block_size);

LIS_INT lis_esolver_set_transform(
    LIS_ESOLVER esolver,
    LIS_INT transform);
LIS_INT lis_esolver_get_transform(
    LIS_ESOLVER esolver,
    LIS_INT *transform);

/*
 * Configure the real and imaginary parts used by -shift and -shift_im.
 * A nonzero imaginary part is rejected by real-scalar builds.
 */
LIS_INT lis_esolver_set_shift(
    LIS_ESOLVER esolver,
    LIS_REAL real_shift,
    LIS_REAL imaginary_shift);
LIS_INT lis_esolver_get_shift(
    LIS_ESOLVER esolver,
    LIS_REAL *real_shift,
    LIS_REAL *imaginary_shift);

#ifdef __cplusplus
}
#endif

#endif
