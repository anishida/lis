/* Mutable lifecycle helpers for non-owning matrix operators. */
#ifndef LIS_MATRIX_OPERATOR_H
#define LIS_MATRIX_OPERATOR_H

#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Replace the scalar coefficients of a linear-combination operator
 * alpha*A + beta*B and publish a value-generation change.
 */
LIS_INT lis_matrix_operator_set_coefficients(
    LIS_MATRIX C,
    LIS_SCALAR alpha,
    LIS_SCALAR beta);

/*
 * Transactionally replace the operands of an existing linear-combination
 * operator while preserving its scalar coefficients. The replacement
 * operands must use the same external row distribution as C. A new operator
 * work vector is prepared before the active operands are changed. On success
 * the matrix structure generation is advanced.
 */
LIS_INT lis_matrix_operator_rebind(
    LIS_MATRIX C,
    LIS_MATRIX A,
    LIS_MATRIX B);

#ifdef __cplusplus
}
#endif

#endif
