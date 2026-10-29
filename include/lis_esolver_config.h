/* Generic typed eigensolver configuration helpers. */
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
