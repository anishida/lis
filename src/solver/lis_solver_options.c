/* Generic solver option-profile helpers. */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <string.h>
#include <math.h>
#include "lislib.h"
#include "lis_solver_options.h"

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_replace"
LIS_INT
lis_solver_set_option_replace(
    const char *text,
    LIS_SOLVER solver)
{
    LIS_SOLVER fresh;
    LIS_INT err;
    LIS_INT i;
    char *copy;
    size_t len;

    LIS_DEBUG_FUNC_IN;

    if(solver==NULL || text==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"solver or option text is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    fresh = NULL;
    copy = NULL;

    err = lis_solver_create(&fresh);
    if(err) return err;

    len = strlen(text);
    copy = (char *)lis_malloc(len+1,"lis_solver_set_option_replace::copy");
    if(copy==NULL)
    {
        lis_solver_destroy(fresh);
        LIS_SETERR_MEM(len+1);
        return LIS_ERR_OUT_OF_MEMORY;
    }
    memcpy(copy,text,len+1);

    /* Parse first so invalid input leaves the active solver unchanged. */
    err = lis_solver_set_option(copy,fresh);
    lis_free(copy);
    if(err)
    {
        lis_solver_destroy(fresh);
        return err;
    }

    /* Allocate FGCRO-DR's private configuration state before any profile
     * reset or mutation, so allocation failure leaves the profile intact. */
    err = lis_fgcrodr_profile_prepare(solver,fresh);
    if(err)
    {
        lis_solver_destroy(fresh);
        return err;
    }

    /*
     * A complete option-profile change invalidates solve-dependent sequence
     * state.  lis_solver_reset_sequence() intentionally preserves application
     * callbacks, nullspace configuration, sequence policy and thread override.
     */
    err = lis_solver_reset_sequence(solver);
    if(err)
    {
        lis_solver_destroy(fresh);
        return err;
    }

    /*
     * Copy only the public solver profile arrays.  Do not bulk-copy solver
     * internals: application-owned callbacks, nullspace objects and persistent
     * sequence policy belong to the existing solver object.
     */
    for(i=0;i<LIS_OPTIONS_LEN;i++)
        solver->options[i] = fresh->options[i];

    for(i=0;i<LIS_PARAMS_LEN;i++)
        solver->params[i] = fresh->params[i];

    /* This parsed option is stored outside the public options[] array. */
    solver->mfamg_block_size = fresh->mfamg_block_size;
    solver->mfamg_fine_agg_factor = fresh->mfamg_fine_agg_factor;
    solver->mfamg_vcycle_smooth_steps = fresh->mfamg_vcycle_smooth_steps;
    solver->mfamg_deep_sibling_merge = fresh->mfamg_deep_sibling_merge;
    solver->ilut_rcm = fresh->ilut_rcm;
    solver->ilut_row_scale = fresh->ilut_row_scale;
    solver->ilutp_enabled = fresh->ilutp_enabled;
    solver->ilutp_threshold = fresh->ilutp_threshold;
    solver->ilutp_full = fresh->ilutp_full;
    solver->ilutp_full_limit = fresh->ilutp_full_limit;
    solver->ilut_equil = fresh->ilut_equil;

    /* GCRO-DR configuration is stored outside options[] / params[]. */
    solver->recycle_target_dim = fresh->recycle_target_dim;
    solver->recycle_keep_threshold = fresh->recycle_keep_threshold;
    solver->recycle_reject_threshold = fresh->recycle_reject_threshold;

    /* Preflight already completed: copy cannot fail or allocate here. */
    lis_fgcrodr_profile_copy(solver,fresh);

    solver->precision = solver->options[LIS_OPTIONS_PRECISION];
    solver->setup = LIS_FALSE;
    solver->retcode = LIS_SUCCESS;
    solver->iter = 0;
    solver->iter2 = 0;
    solver->resid = (LIS_REAL)0.0;

    lis_solver_destroy(fresh);

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_replace_if_changed"
LIS_INT
lis_solver_set_option_replace_if_changed(
    const char *text,
    LIS_SOLVER solver,
    LIS_INT *changed)
{
    LIS_SOLVER fresh;
    LIS_INT err;
    LIS_INT i;
    LIS_INT different;
    char *copy;
    size_t len;

    LIS_DEBUG_FUNC_IN;

    if(changed)
        *changed = LIS_FALSE;

    if(solver==NULL || text==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"solver or option text is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    fresh = NULL;
    copy = NULL;
    different = LIS_FALSE;

    err = lis_solver_create(&fresh);
    if(err)
        return err;

    len = strlen(text);
    copy = (char *)lis_malloc(
        len+1,
        "lis_solver_set_option_replace_if_changed::copy");
    if(copy==NULL)
    {
        lis_solver_destroy(fresh);
        LIS_SETERR_MEM(len+1);
        return LIS_ERR_OUT_OF_MEMORY;
    }

    memcpy(copy,text,len+1);

    err = lis_solver_set_option(copy,fresh);
    lis_free(copy);
    if(err)
    {
        lis_solver_destroy(fresh);
        return err;
    }

    for(i=0;i<LIS_OPTIONS_LEN;i++)
    {
        if(solver->options[i]!=fresh->options[i])
        {
            different = LIS_TRUE;
            break;
        }
    }

    if(!different)
    {
        for(i=0;i<LIS_PARAMS_LEN;i++)
        {
            if(solver->params[i]!=fresh->params[i])
            {
                different = LIS_TRUE;
                break;
            }
        }
    }

    if(!different && solver->mfamg_block_size!=fresh->mfamg_block_size)
        different = LIS_TRUE;
    if(!different && solver->mfamg_fine_agg_factor!=fresh->mfamg_fine_agg_factor)
        different = LIS_TRUE;
    if(!different && solver->mfamg_vcycle_smooth_steps!=fresh->mfamg_vcycle_smooth_steps)
        different = LIS_TRUE;
    if(!different && solver->mfamg_deep_sibling_merge!=fresh->mfamg_deep_sibling_merge)
        different = LIS_TRUE;
    if(!different && (solver->ilut_rcm!=fresh->ilut_rcm ||
                      solver->ilut_row_scale!=fresh->ilut_row_scale ||
                      solver->ilutp_enabled!=fresh->ilutp_enabled ||
                      solver->ilutp_threshold!=fresh->ilutp_threshold ||
                      solver->ilutp_full!=fresh->ilutp_full ||
                      solver->ilutp_full_limit!=fresh->ilutp_full_limit ||
                      solver->ilut_equil!=fresh->ilut_equil))
        different = LIS_TRUE;

    /* Include solver-specific GCRO-DR profile fields in comparisons. */
    if(!different &&
       (solver->recycle_target_dim!=fresh->recycle_target_dim ||
        solver->recycle_keep_threshold!=fresh->recycle_keep_threshold ||
        solver->recycle_reject_threshold!=fresh->recycle_reject_threshold))
        different = LIS_TRUE;

    if(!different && !lis_fgcrodr_profile_equal(solver,fresh))
        different = LIS_TRUE;

    lis_solver_destroy(fresh);

    if(!different)
    {
        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
    }

    err = lis_solver_set_option_replace(text,solver);
    if(err)
        return err;

    if(changed)
        *changed = LIS_TRUE;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_set_tolerance"
LIS_INT
lis_solver_set_tolerance(
    LIS_SOLVER solver,
    LIS_REAL tolerance)
{
    LIS_DEBUG_FUNC_IN;

    if(solver==NULL || !isfinite((double)tolerance) || tolerance<0.0)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "solver is NULL or tolerance is invalid\n");
        return LIS_ERR_ILL_ARG;
    }

    solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN] =
        (LIS_SCALAR)tolerance;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_get_tolerance"
LIS_INT
lis_solver_get_tolerance(
    LIS_SOLVER solver,
    LIS_REAL *tolerance)
{
    LIS_DEBUG_FUNC_IN;

    if(solver==NULL || tolerance==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "solver or tolerance output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *tolerance = (LIS_REAL)
        solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_set_maxiter"
LIS_INT
lis_solver_set_maxiter(
    LIS_SOLVER solver,
    LIS_INT maxiter)
{
    LIS_DEBUG_FUNC_IN;

    if(solver==NULL || maxiter<0)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "solver is NULL or maximum iteration count is negative\n");
        return LIS_ERR_ILL_ARG;
    }

    solver->options[LIS_OPTIONS_MAXITER] = maxiter;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_get_maxiter"
LIS_INT
lis_solver_get_maxiter(
    LIS_SOLVER solver,
    LIS_INT *maxiter)
{
    LIS_DEBUG_FUNC_IN;

    if(solver==NULL || maxiter==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "solver or maximum-iteration output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *maxiter = solver->options[LIS_OPTIONS_MAXITER];

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_get_report"
LIS_INT
lis_solver_get_report(
    LIS_SOLVER solver,
    LIS_SOLVE_REPORT *report)
{
    LIS_SOLVE_REPORT tmp;
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    if(solver==NULL || report==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "solver or report output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    memset(&tmp,0,sizeof(tmp));

    err = lis_solver_get_status(solver,&tmp.status);
    if(err) return err;

    err = lis_solver_get_iterex(
        solver,
        &tmp.iterations,
        &tmp.iterations_double,
        &tmp.iterations_quad);
    if(err) return err;

    err = lis_solver_get_residualnorm(solver,&tmp.residual_norm);
    if(err) return err;

    err = lis_solver_get_timeex(
        solver,
        &tmp.time,
        &tmp.iteration_time,
        &tmp.preconditioner_time,
        &tmp.preconditioner_create_time,
        &tmp.preconditioner_iteration_time);
    if(err) return err;

    *report = tmp;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}
