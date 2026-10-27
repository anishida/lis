/* Generic solver option-profile helpers. */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <string.h>
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
