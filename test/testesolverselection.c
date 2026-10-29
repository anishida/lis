/* Focused regression for typed eigensolver configuration and pair selection. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lis.h"
#include "lis_esolver_config.h"
#include "lis_esolver_selection.h"

static int fail(const char *msg, int err)
{
    fprintf(stderr,"FAIL: %s (err=%d)\n",msg,err);
    return 1;
}

int main(int argc, char *argv[])
{
    LIS_ESOLVER esolver = NULL;
    LIS_ESOLVER_SELECTED_PAIR pairs[4];
    LIS_REAL tolerance,real_shift,imaginary_shift;
    LIS_INT type,maxiter,subspace_size,mode,block_size,transform;
    LIS_INT available = 0;
    LIS_INT err;

    err = lis_initialize(&argc,&argv);
    if(err) return fail("lis_initialize",err);

    err = lis_esolver_create(&esolver);
    if(err) return fail("lis_esolver_create",err);

    /* Configure the common embedding controls without parsing option text. */
    err = lis_esolver_set_type(esolver,LIS_ESOLVER_LI);
    if(err) return fail("set type",err);
    err = lis_esolver_set_subspace_size(esolver,4);
    if(err) return fail("set subspace",err);
    err = lis_esolver_set_mode(esolver,1);
    if(err) return fail("set mode",err);
    err = lis_esolver_set_block_size(esolver,2);
    if(err) return fail("set block size",err);
    err = lis_esolver_set_tolerance(esolver,(LIS_REAL)1.0e-9);
    if(err) return fail("set tolerance",err);
    err = lis_esolver_set_maxiter(esolver,123);
    if(err) return fail("set maxiter",err);
    err = lis_esolver_set_transform(
        esolver,LIS_EIGEN_TRANSFORM_SHIFT_INVERT);
    if(err) return fail("set transform",err);
    err = lis_esolver_set_shift(esolver,(LIS_REAL)2.5,(LIS_REAL)0.0);
    if(err) return fail("set shift",err);

    err = lis_esolver_get_type(esolver,&type);
    if(err) return fail("get type",err);
    err = lis_esolver_get_subspace_size(esolver,&subspace_size);
    if(err) return fail("get subspace",err);
    err = lis_esolver_get_mode(esolver,&mode);
    if(err) return fail("get mode",err);
    err = lis_esolver_get_block_size(esolver,&block_size);
    if(err) return fail("get block size",err);
    err = lis_esolver_get_tolerance(esolver,&tolerance);
    if(err) return fail("get tolerance",err);
    err = lis_esolver_get_maxiter(esolver,&maxiter);
    if(err) return fail("get maxiter",err);
    err = lis_esolver_get_transform(esolver,&transform);
    if(err) return fail("get transform",err);
    err = lis_esolver_get_shift(esolver,&real_shift,&imaginary_shift);
    if(err) return fail("get shift",err);

    if(type!=LIS_ESOLVER_LI ||
       subspace_size!=4 || mode!=1 || block_size!=2 ||
       tolerance!=(LIS_REAL)1.0e-9 || maxiter!=123 ||
       transform!=LIS_EIGEN_TRANSFORM_SHIFT_INVERT ||
       real_shift!=(LIS_REAL)2.5 || imaginary_shift!=(LIS_REAL)0.0)
        return fail("typed configuration round trip",LIS_BREAKDOWN);

    /* Invalid typed updates must not publish partial configuration. */
    err = lis_esolver_set_type(esolver,LIS_ESOLVER_CG);
    if(err==LIS_SUCCESS) return fail("invalid block/type accepted",LIS_BREAKDOWN);
    err = lis_esolver_set_tolerance(esolver,(LIS_REAL)-1.0);
    if(err==LIS_SUCCESS) return fail("invalid tolerance accepted",LIS_BREAKDOWN);
    err = lis_esolver_set_maxiter(esolver,0);
    if(err==LIS_SUCCESS) return fail("invalid maxiter accepted",LIS_BREAKDOWN);
    err = lis_esolver_set_subspace_size(esolver,1);
    if(err==LIS_SUCCESS) return fail("invalid subspace accepted",LIS_BREAKDOWN);
    err = lis_esolver_set_mode(esolver,4);
    if(err==LIS_SUCCESS) return fail("invalid mode accepted",LIS_BREAKDOWN);
    err = lis_esolver_set_block_size(esolver,5);
    if(err==LIS_SUCCESS) return fail("invalid block size accepted",LIS_BREAKDOWN);
    err = lis_esolver_set_transform(esolver,99);
    if(err==LIS_SUCCESS) return fail("invalid transform accepted",LIS_BREAKDOWN);
    err = lis_esolver_set_shift(esolver,(LIS_REAL)NAN,(LIS_REAL)0.0);
    if(err==LIS_SUCCESS) return fail("invalid shift accepted",LIS_BREAKDOWN);

#ifndef _COMPLEX
    err = lis_esolver_set_shift(esolver,(LIS_REAL)2.5,(LIS_REAL)0.25);
    if(err==LIS_SUCCESS)
        return fail("imaginary shift accepted in real build",LIS_BREAKDOWN);
#else
    err = lis_esolver_set_shift(esolver,(LIS_REAL)2.5,(LIS_REAL)0.25);
    if(err) return fail("complex shift",err);
    err = lis_esolver_get_shift(esolver,&real_shift,&imaginary_shift);
    if(err || real_shift!=(LIS_REAL)2.5 || imaginary_shift!=(LIS_REAL)0.25)
        return fail("complex shift round trip",err ? err : LIS_BREAKDOWN);
    err = lis_esolver_set_shift(esolver,(LIS_REAL)2.5,(LIS_REAL)0.0);
    if(err) return fail("reset complex shift",err);
#endif

    lis_esolver_get_type(esolver,&type);
    lis_esolver_get_subspace_size(esolver,&subspace_size);
    lis_esolver_get_mode(esolver,&mode);
    lis_esolver_get_block_size(esolver,&block_size);
    lis_esolver_get_tolerance(esolver,&tolerance);
    lis_esolver_get_maxiter(esolver,&maxiter);
    lis_esolver_get_transform(esolver,&transform);
    lis_esolver_get_shift(esolver,&real_shift,&imaginary_shift);

    if(type!=LIS_ESOLVER_LI ||
       subspace_size!=4 || mode!=1 || block_size!=2 ||
       tolerance!=(LIS_REAL)1.0e-9 || maxiter!=123 ||
       transform!=LIS_EIGEN_TRANSFORM_SHIFT_INVERT ||
       real_shift!=(LIS_REAL)2.5 || imaginary_shift!=(LIS_REAL)0.0)
        return fail("failed update changed active configuration",LIS_BREAKDOWN);

    /* Dependent setters require an order that preserves valid intermediate states. */
    err = lis_esolver_set_block_size(esolver,1);
    if(err) return fail("decrease block before changing type",err);
    err = lis_esolver_set_type(esolver,LIS_ESOLVER_CG);
    if(err) return fail("change type after reducing block",err);
    err = lis_esolver_set_type(esolver,LIS_ESOLVER_LI);
    if(err) return fail("restore Lanczos type",err);
    err = lis_esolver_set_block_size(esolver,2);
    if(err) return fail("restore block size",err);

    /* Non-finite tolerances and shifts must not affect published values. */
    err = lis_esolver_set_tolerance(esolver,(LIS_REAL)NAN);
    if(err==LIS_SUCCESS)
        return fail("non-finite tolerance accepted",LIS_BREAKDOWN);
    err = lis_esolver_set_tolerance(esolver,(LIS_REAL)INFINITY);
    if(err==LIS_SUCCESS)
        return fail("infinite tolerance accepted",LIS_BREAKDOWN);
    err = lis_esolver_set_shift(esolver,(LIS_REAL)INFINITY,(LIS_REAL)0.0);
    if(err==LIS_SUCCESS)
        return fail("infinite shift accepted",LIS_BREAKDOWN);

    lis_esolver_get_tolerance(esolver,&tolerance);
    lis_esolver_get_shift(esolver,&real_shift,&imaginary_shift);
    if(tolerance!=(LIS_REAL)1.0e-9 ||
       real_shift!=(LIS_REAL)2.5 || imaginary_shift!=(LIS_REAL)0.0)
        return fail("non-finite update changed state",LIS_BREAKDOWN);

    /* Cover every boundary of the typed setters and reject partial state. */
    if(lis_esolver_set_type(esolver,0)==LIS_SUCCESS ||
       lis_esolver_set_type(esolver,LIS_ESOLVER_LEN+1)==LIS_SUCCESS ||
       lis_esolver_set_tolerance(esolver,(LIS_REAL)0.0)==LIS_SUCCESS ||
       lis_esolver_set_tolerance(esolver,(LIS_REAL)-INFINITY)==LIS_SUCCESS ||
       lis_esolver_set_maxiter(esolver,-1)==LIS_SUCCESS ||
       lis_esolver_set_subspace_size(esolver,0)==LIS_SUCCESS ||
       lis_esolver_set_subspace_size(esolver,2)==LIS_SUCCESS ||
       lis_esolver_set_mode(esolver,-1)==LIS_SUCCESS ||
       lis_esolver_set_block_size(esolver,0)==LIS_SUCCESS ||
       lis_esolver_set_block_size(esolver,5)==LIS_SUCCESS ||
       lis_esolver_set_transform(esolver,-1)==LIS_SUCCESS ||
       lis_esolver_set_shift(esolver,(LIS_REAL)-INFINITY,(LIS_REAL)0.0)==LIS_SUCCESS ||
       lis_esolver_set_shift(esolver,(LIS_REAL)2.5,(LIS_REAL)INFINITY)==LIS_SUCCESS ||
       lis_esolver_set_shift(esolver,(LIS_REAL)2.5,(LIS_REAL)NAN)==LIS_SUCCESS)
        return fail("additional invalid configuration accepted",LIS_BREAKDOWN);

    /* Null outputs and objects must be rejected, not dereferenced. */
    if(lis_esolver_get_type(esolver,NULL)==LIS_SUCCESS ||
       lis_esolver_get_tolerance(esolver,NULL)==LIS_SUCCESS ||
       lis_esolver_get_maxiter(esolver,NULL)==LIS_SUCCESS ||
       lis_esolver_get_subspace_size(esolver,NULL)==LIS_SUCCESS ||
       lis_esolver_get_mode(esolver,NULL)==LIS_SUCCESS ||
       lis_esolver_get_block_size(esolver,NULL)==LIS_SUCCESS ||
       lis_esolver_get_transform(esolver,NULL)==LIS_SUCCESS ||
       lis_esolver_get_shift(esolver,NULL,&imaginary_shift)==LIS_SUCCESS ||
       lis_esolver_get_shift(esolver,&real_shift,NULL)==LIS_SUCCESS ||
       lis_esolver_get_type(NULL,&type)==LIS_SUCCESS)
        return fail("invalid getter arguments accepted",LIS_BREAKDOWN);

    /* Verify rollback also for every additional rejected update. */
    if(lis_esolver_get_type(esolver,&type) ||
       lis_esolver_get_subspace_size(esolver,&subspace_size) ||
       lis_esolver_get_mode(esolver,&mode) ||
       lis_esolver_get_block_size(esolver,&block_size) ||
       lis_esolver_get_tolerance(esolver,&tolerance) ||
       lis_esolver_get_maxiter(esolver,&maxiter) ||
       lis_esolver_get_transform(esolver,&transform) ||
       lis_esolver_get_shift(esolver,&real_shift,&imaginary_shift) ||
       type!=LIS_ESOLVER_LI || subspace_size!=4 || mode!=1 ||
       block_size!=2 || tolerance!=(LIS_REAL)1.0e-9 ||
       maxiter!=123 || transform!=LIS_EIGEN_TRANSFORM_SHIFT_INVERT ||
       real_shift!=(LIS_REAL)2.5 || imaginary_shift!=(LIS_REAL)0.0)
        return fail("extended invalid updates modified configuration",LIS_BREAKDOWN);

    /* Existing converged-pair selection remains unchanged. */
    esolver->nev = 4;

    esolver->evalue = (LIS_SCALAR *)calloc(6,sizeof(LIS_SCALAR));
    esolver->resid = (LIS_REAL *)calloc(6,sizeof(LIS_REAL));
    if(esolver->evalue==NULL || esolver->resid==NULL)
        return fail("allocation",LIS_OUT_OF_MEMORY);

    esolver->evalue[0] = (LIS_SCALAR)3.0;
    esolver->evalue[1] = (LIS_SCALAR)(-1.0);
    esolver->evalue[2] = (LIS_SCALAR)2.0;
    esolver->evalue[3] = (LIS_SCALAR)1.0;

    esolver->resid[0] = (LIS_REAL)1.0e-8;
    esolver->resid[1] = (LIS_REAL)1.0e-8;
    esolver->resid[2] = (LIS_REAL)1.0e-3;
    esolver->resid[3] = (LIS_REAL)1.0e-9;

    /*
     * Scalar controls after results exist must not modify stored eigenpair
     * arrays. These results remain from the previous configuration/solve,
     * not newly computed pairs for the updated tolerances or transforms.
     */
    {
        LIS_SCALAR *old_values=esolver->evalue;
        LIS_REAL *old_residuals=esolver->resid;
        LIS_INT old_nev=esolver->nev;
        err=lis_esolver_set_tolerance(esolver,(LIS_REAL)2.0e-9);
        if(err)return fail("post-result tolerance update",err);
        err=lis_esolver_set_transform(esolver,LIS_EIGEN_TRANSFORM_AUTO);
        if(err)return fail("post-result transform update",err);
        err=lis_esolver_set_shift(esolver,(LIS_REAL)1.5,(LIS_REAL)0.0);
        if(err)return fail("post-result shift update",err);
        if(esolver->evalue!=old_values || esolver->resid!=old_residuals ||
           esolver->nev!=old_nev ||
           esolver->evalue[0]!=(LIS_SCALAR)3.0 ||
           esolver->resid[3]!=(LIS_REAL)1.0e-9)
            return fail("post-result configuration mutated saved results",LIS_BREAKDOWN);
    }

    err = lis_esolver_select_converged_pairs(
        esolver,4,2,(LIS_REAL)1.0e-6,
        LIS_ESOLVER_SELECT_ASCENDING_REAL,
        pairs,&available);
    if(err) return fail("algebraic selection",err);

    if(available!=3 ||
       (double)pairs[0].value!=-1.0 ||
       (double)pairs[1].value!=1.0 ||
       pairs[0].mode!=1 || pairs[1].mode!=3)
        return fail("algebraic ordering",LIS_BREAKDOWN);

    err = lis_esolver_select_converged_pairs(
        esolver,4,2,(LIS_REAL)1.0e-6,
        LIS_ESOLVER_SELECT_ASCENDING_REAL |
        LIS_ESOLVER_SELECT_POSITIVE_REAL,
        pairs,&available);
    if(err) return fail("positive selection",err);

    if(available!=2 ||
       (double)pairs[0].value!=1.0 ||
       (double)pairs[1].value!=3.0 ||
       pairs[0].mode!=3 || pairs[1].mode!=0)
        return fail("positive ordering",LIS_BREAKDOWN);

    available = -1;
    err = lis_esolver_select_converged_pairs(
        esolver,6,1,(LIS_REAL)1.0e-6,
        LIS_ESOLVER_SELECT_ASCENDING_REAL,
        pairs,&available);
    if(err==LIS_SUCCESS)
        return fail("out-of-range scan error propagation",LIS_BREAKDOWN);

    printf("PASS: typed eigensolver configuration and pair selection\n");

    lis_esolver_destroy(esolver);
    lis_finalize();
    return 0;
}
