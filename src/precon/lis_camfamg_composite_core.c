#include "lis_camfamg_composite_core.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

lis_cam_status lis_cam_apply_multiplicative(const lis_cam_composite *ops,
                                             const double *r,double *z){
    size_t i,n;double *buf=NULL,*z1,*az1,*r1,*z2,*out;
    double alpha=1.0;
    if(!ops||!r||!z||!ops->matvec||!ops->ca_apply||!ops->ilutp_apply||!ops->n)
        return LIS_CAM_INVALID;
    n=ops->n;
    if(n>SIZE_MAX/(5*sizeof(double)))return LIS_CAM_ALLOCATION;
    buf=(double*)calloc(5*n,sizeof(double));
    if(!buf)return LIS_CAM_ALLOCATION;
    z1=buf;az1=buf+n;r1=buf+2*n;z2=buf+3*n;out=buf+4*n;
    for(i=0;i<n;i++)if(!isfinite(r[i])){free(buf);return LIS_CAM_NONFINITE;}
    if(ops->ca_apply(ops->ca_ctx,r,z1)!=0){free(buf);return LIS_CAM_CALLBACK_ERROR;}
    for(i=0;i<n;i++)if(!isfinite(z1[i])){free(buf);return LIS_CAM_NONFINITE;}
    if(ops->matvec(ops->matrix_ctx,z1,az1)!=0){free(buf);return LIS_CAM_CALLBACK_ERROR;}
    for(i=0;i<n;i++)if(!isfinite(az1[i])){free(buf);return LIS_CAM_NONFINITE;}
    if(ops->safe_damping){
        /* alpha minimizes ||r-alpha A z1||_2.  Avoid scale-dependent 1e-30 cutoffs. */
        double scale=0,num=0,den=0;
        for(i=0;i<n;i++){double t=fabs(az1[i]);if(t>scale)scale=t;}
        if(scale>0){
            for(i=0;i<n;i++){double v=az1[i]/scale;num+=r[i]*v;den+=v*v;}
            if(den>0)alpha=(num/den)/scale;
            else alpha=0;
        } else alpha=0;
        /* No clipping to a matrix-/snapshot-specific bound. */
        if(!isfinite(alpha)){free(buf);return LIS_CAM_NONFINITE;}
    }
    for(i=0;i<n;i++){
        r1[i]=r[i]-alpha*az1[i];
        if(!isfinite(r1[i])){free(buf);return LIS_CAM_NONFINITE;}
    }
    if(ops->ilutp_apply(ops->ilutp_ctx,r1,z2)!=0){free(buf);return LIS_CAM_CALLBACK_ERROR;}
    for(i=0;i<n;i++){
        out[i]=alpha*z1[i]+z2[i];
        if(!isfinite(out[i])){free(buf);return LIS_CAM_NONFINITE;}
    }
    for(i=0;i<n;i++)z[i]=out[i];
    free(buf);
    return LIS_CAM_OK;
}
