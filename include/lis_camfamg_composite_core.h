/* Callback-based, non-mutating multiplicative preconditioner composition.
 * Actual LIS preconditioner handles are bound in an adapter, not here.
 */
#ifndef LIS_CAMFAMG_COMPOSITE_CORE_H
#define LIS_CAMFAMG_COMPOSITE_CORE_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef int (*lis_cam_apply_fn)(void *user, const double *rhs, double *out);

typedef struct {
    size_t n;
    lis_cam_apply_fn matvec; void *matrix_ctx;
    lis_cam_apply_fn ca_apply; void *ca_ctx;
    lis_cam_apply_fn ilutp_apply; void *ilutp_ctx;
    /* If true, use safe min-residual damping of CA step; otherwise alpha=1.
     * Optional damping generally makes preconditioning nonlinear -> FGMRES.
     */
    int safe_damping;
} lis_cam_composite;

typedef enum {LIS_CAM_OK=0,LIS_CAM_INVALID=1,LIS_CAM_NONFINITE=2,
              LIS_CAM_ALLOCATION=3,LIS_CAM_CALLBACK_ERROR=4} lis_cam_status;

lis_cam_status lis_cam_apply_multiplicative(const lis_cam_composite *ops,
                                             const double *r, double *z);
#ifdef __cplusplus
}
#endif
#endif
