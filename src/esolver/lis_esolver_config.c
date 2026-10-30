/* Generic typed eigensolver configuration helpers. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>

#include "lislib.h"
#include "lis_esolver_config.h"

#undef __FUNC__
#define __FUNC__ "lis_esolver_config_check"
static LIS_INT
lis_esolver_config_check(LIS_ESOLVER esolver)
{
    if(esolver==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver is undefined\n");
        return LIS_ERR_ILL_ARG;
    }

    return LIS_SUCCESS;
}

static LIS_INT
lis_esolver_config_reset_work(LIS_ESOLVER esolver)
{
    if(esolver->work==NULL)
        return LIS_SUCCESS;

    return lis_esolver_work_destroy(esolver);
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_type"
LIS_INT
lis_esolver_set_type(
    LIS_ESOLVER esolver,
    LIS_INT type)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;

    if(type<1 || type>LIS_ESOLVER_LEN)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid eigensolver type\n");
        return LIS_ERR_ILL_ARG;
    }

    if(esolver->block_size>1 &&
       type!=LIS_ESOLVER_LI &&
       type!=LIS_ESOLVER_GLI)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "block eigensolver mode requires Lanczos\n");
        return LIS_ERR_ILL_ARG;
    }

    if(esolver->options[LIS_EOPTIONS_ESOLVER]==type)
    {
        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
    }

    err = lis_esolver_config_reset_work(esolver);
    if(err) return err;

    esolver->options[LIS_EOPTIONS_ESOLVER] = type;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_type"
LIS_INT
lis_esolver_get_type(
    LIS_ESOLVER esolver,
    LIS_INT *type)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;
    if(type==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver type output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *type = esolver->options[LIS_EOPTIONS_ESOLVER];

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_tolerance"
LIS_INT
lis_esolver_set_tolerance(
    LIS_ESOLVER esolver,
    LIS_REAL tolerance)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;

    if(!isfinite((double)tolerance) || tolerance<=0.0)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid eigensolver tolerance\n");
        return LIS_ERR_ILL_ARG;
    }

    esolver->params[LIS_EPARAMS_RESID-LIS_EOPTIONS_LEN] =
        (LIS_SCALAR)tolerance;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_tolerance"
LIS_INT
lis_esolver_get_tolerance(
    LIS_ESOLVER esolver,
    LIS_REAL *tolerance)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;
    if(tolerance==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver tolerance output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *tolerance = (LIS_REAL)
        esolver->params[LIS_EPARAMS_RESID-LIS_EOPTIONS_LEN];

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_maxiter"
LIS_INT
lis_esolver_set_maxiter(
    LIS_ESOLVER esolver,
    LIS_INT maxiter)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;

    if(maxiter<1)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "eigensolver maximum iteration count must be positive\n");
        return LIS_ERR_ILL_ARG;
    }

    esolver->options[LIS_EOPTIONS_MAXITER] = maxiter;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_maxiter"
LIS_INT
lis_esolver_get_maxiter(
    LIS_ESOLVER esolver,
    LIS_INT *maxiter)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;
    if(maxiter==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver maxiter output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *maxiter = esolver->options[LIS_EOPTIONS_MAXITER];

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_subspace_size"
LIS_INT
lis_esolver_set_subspace_size(
    LIS_ESOLVER esolver,
    LIS_INT subspace_size)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;

    if(subspace_size<1 ||
       esolver->options[LIS_EOPTIONS_MODE]>=subspace_size ||
       esolver->block_size>subspace_size ||
       (esolver->nev>0 && esolver->nev>subspace_size))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "invalid eigensolver subspace size for active configuration\n");
        return LIS_ERR_ILL_ARG;
    }

    if(esolver->options[LIS_EOPTIONS_SUBSPACE]==subspace_size)
    {
        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
    }

    err = lis_esolver_config_reset_work(esolver);
    if(err) return err;

    esolver->options[LIS_EOPTIONS_SUBSPACE] = subspace_size;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_subspace_size"
LIS_INT
lis_esolver_get_subspace_size(
    LIS_ESOLVER esolver,
    LIS_INT *subspace_size)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;
    if(subspace_size==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver subspace output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *subspace_size = esolver->options[LIS_EOPTIONS_SUBSPACE];

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

/*
 * -nev shares the existing esolver->nev slot with the result count of some
 * eigensolver algorithms. This setter configures the requested count before
 * the next solve, without changing any public LIS_EOPTIONS_* identifiers.
 */
#undef __FUNC__
#define __FUNC__ "lis_esolver_set_requested_nev"
LIS_INT
lis_esolver_set_requested_nev(
    LIS_ESOLVER esolver,
    LIS_INT requested_nev)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;

    if(requested_nev<0 ||
       requested_nev>esolver->options[LIS_EOPTIONS_SUBSPACE])
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "requested eigenpair count must be zero or fit the subspace\n");
        return LIS_ERR_ILL_ARG;
    }

    if(esolver->nev==requested_nev)
    {
        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
    }

    err = lis_esolver_config_reset_work(esolver);
    if(err) return err;

    esolver->nev = requested_nev;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_requested_nev"
LIS_INT
lis_esolver_get_requested_nev(
    LIS_ESOLVER esolver,
    LIS_INT *requested_nev)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;
    if(requested_nev==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "requested eigenpair count output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *requested_nev = esolver->nev;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_mode"
LIS_INT
lis_esolver_set_mode(
    LIS_ESOLVER esolver,
    LIS_INT mode)
{
    LIS_INT err;
    LIS_INT subspace_size;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;

    subspace_size = esolver->options[LIS_EOPTIONS_SUBSPACE];
    if(mode<0 || mode>=subspace_size)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "eigensolver mode must be within the active subspace\n");
        return LIS_ERR_ILL_ARG;
    }

    esolver->options[LIS_EOPTIONS_MODE] = mode;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_mode"
LIS_INT
lis_esolver_get_mode(
    LIS_ESOLVER esolver,
    LIS_INT *mode)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;
    if(mode==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver mode output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *mode = esolver->options[LIS_EOPTIONS_MODE];

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_block_size"
LIS_INT
lis_esolver_set_block_size(
    LIS_ESOLVER esolver,
    LIS_INT block_size)
{
    LIS_INT err;
    LIS_INT type;
    LIS_INT subspace_size;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;

    type = esolver->options[LIS_EOPTIONS_ESOLVER];
    subspace_size = esolver->options[LIS_EOPTIONS_SUBSPACE];

    if(block_size<1 || block_size>subspace_size)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "invalid eigensolver block size\n");
        return LIS_ERR_ILL_ARG;
    }

    if(block_size>1 &&
       type!=LIS_ESOLVER_LI &&
       type!=LIS_ESOLVER_GLI)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "block eigensolver mode requires Lanczos\n");
        return LIS_ERR_ILL_ARG;
    }

    if(esolver->block_size==block_size)
    {
        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
    }

    err = lis_esolver_config_reset_work(esolver);
    if(err) return err;

    esolver->block_size = block_size;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_block_size"
LIS_INT
lis_esolver_get_block_size(
    LIS_ESOLVER esolver,
    LIS_INT *block_size)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;
    if(block_size==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver block-size output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *block_size = esolver->block_size;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_transform"
LIS_INT
lis_esolver_set_transform(
    LIS_ESOLVER esolver,
    LIS_INT transform)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;

    if(transform!=LIS_EIGEN_TRANSFORM_AUTO &&
       transform!=LIS_EIGEN_TRANSFORM_SHIFT_INVERT &&
       transform!=LIS_EIGEN_TRANSFORM_BUCKLING)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "invalid eigensolver spectral transform\n");
        return LIS_ERR_ILL_ARG;
    }

    esolver->transform = transform;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_transform"
LIS_INT
lis_esolver_get_transform(
    LIS_ESOLVER esolver,
    LIS_INT *transform)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;
    if(transform==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver transform output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *transform = esolver->transform;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_shift"
LIS_INT
lis_esolver_set_shift(
    LIS_ESOLVER esolver,
    LIS_REAL real_shift,
    LIS_REAL imaginary_shift)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;

    if(!isfinite((double)real_shift) ||
       !isfinite((double)imaginary_shift))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid eigensolver shift\n");
        return LIS_ERR_ILL_ARG;
    }

#ifndef _COMPLEX
    if(imaginary_shift!=0.0)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "imaginary eigensolver shift requires complex scalars\n");
        return LIS_ERR_ILL_ARG;
    }
#endif

    esolver->params[LIS_EPARAMS_SHIFT-LIS_EOPTIONS_LEN] =
        (LIS_SCALAR)real_shift;
    esolver->params[LIS_EPARAMS_SHIFT_IM-LIS_EOPTIONS_LEN] =
        (LIS_SCALAR)imaginary_shift;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_shift"
LIS_INT
lis_esolver_get_shift(
    LIS_ESOLVER esolver,
    LIS_REAL *real_shift,
    LIS_REAL *imaginary_shift)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_esolver_config_check(esolver);
    if(err) return err;
    if(real_shift==NULL || imaginary_shift==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"eigensolver shift output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *real_shift = (LIS_REAL)
        esolver->params[LIS_EPARAMS_SHIFT-LIS_EOPTIONS_LEN];
    *imaginary_shift = (LIS_REAL)
        esolver->params[LIS_EPARAMS_SHIFT_IM-LIS_EOPTIONS_LEN];

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}
