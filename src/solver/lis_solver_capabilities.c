/* Generic solver/preconditioner capability queries. */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include "lislib.h"
#include "lis_solver_capabilities.h"

#undef __FUNC__
#define __FUNC__ "lis_solver_get_effective_preconditioner"
LIS_INT
lis_solver_get_effective_preconditioner(
    LIS_SOLVER solver,
    LIS_INT *precon_type,
    LIS_INT *nested)
{
    LIS_INT outer;
    LIS_INT inner;

    if(solver==NULL || precon_type==NULL || nested==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid effective preconditioner query\n");
        return LIS_ERR_ILL_ARG;
    }

    outer = solver->options[LIS_OPTIONS_PRECON];
    if(outer!=LIS_PRECON_TYPE_HYBRID)
    {
        *precon_type = outer;
        *nested = LIS_FALSE;
        return LIS_SUCCESS;
    }

    inner = solver->options[LIS_OPTIONS_PPRECON];
    if(inner==LIS_PRECON_TYPE_HYBRID)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "recursive Hybrid preconditioner composition is not supported\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }

    *precon_type = inner;
    *nested = LIS_TRUE;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_effective_capabilities"
LIS_INT
lis_solver_get_effective_capabilities(
    LIS_SOLVER solver,
    LIS_INT *capabilities)
{
    LIS_INT nsolver;
    LIS_INT precon_type;
    LIS_INT nested;
    LIS_INT caps;
    LIS_INT err;

    if(solver==NULL || capabilities==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid solver capability query\n");
        return LIS_ERR_ILL_ARG;
    }

    nsolver = solver->options[LIS_OPTIONS_SOLVER];
    err = lis_solver_get_effective_preconditioner(
        solver,&precon_type,&nested);
    if(err) return err;

    caps = 0;

    if(nsolver==LIS_SOLVER_FGMRES ||
       nsolver==LIS_SOLVER_FGCRODR)
    {
        caps |= LIS_SOLVER_CAP_FLEXIBLE_PRECON;
    }

    if(nsolver==LIS_SOLVER_GCRODR ||
       nsolver==LIS_SOLVER_FGCRODR)
    {
        caps |= LIS_SOLVER_CAP_RECYCLING;
    }

    if(precon_type==LIS_PRECON_TYPE_MFAMG)
    {
        caps |= LIS_SOLVER_CAP_NEAR_NULLSPACE;
    }

    if(nested)
    {
        caps |= LIS_SOLVER_CAP_NESTED_PRECON;
    }

    *capabilities = caps;
    return LIS_SUCCESS;
}

/*
 * Configuration-only nullspace preflight, mirroring the compatibility
 * conditions consumed by lis_solve_kernel().  Do not silently change the
 * requested algorithm or discard any application nullspace vectors.
 */
#undef __FUNC__
#define __FUNC__ "lis_solver_check_nullspace_configuration"
LIS_INT
lis_solver_check_nullspace_configuration(
    LIS_SOLVER solver,
    LIS_INT *issues)
{
    LIS_INT nsolver,precon_type,scale,precision,flags=0;

    if(solver==NULL || issues==NULL)
        return LIS_ERR_ILL_ARG;

    nsolver=solver->options[LIS_OPTIONS_SOLVER];
    precon_type=solver->options[LIS_OPTIONS_PRECON];
    scale=solver->options[LIS_OPTIONS_SCALE];
    precision=solver->options[LIS_OPTIONS_PRECISION];

    if(solver->nullspace_dim>0)
    {
        if(nsolver!=LIS_SOLVER_GMRES &&
           nsolver!=LIS_SOLVER_FGMRES &&
           nsolver!=LIS_SOLVER_GCRODR &&
           nsolver!=LIS_SOLVER_FGCRODR)
            flags |= LIS_NULLSPACE_CONFIG_EXACT_SOLVER;
        if(scale!=LIS_SCALE_NONE)
            flags |= LIS_NULLSPACE_CONFIG_EXACT_SCALE;
        if(precision!=LIS_PRECISION_DEFAULT)
            flags |= LIS_NULLSPACE_CONFIG_EXACT_PRECISION;
    }

    if(solver->near_nullspace_dim>0)
    {
        if(nsolver!=LIS_SOLVER_FGMRES &&
           nsolver!=LIS_SOLVER_FGCRODR &&
           precon_type!=LIS_PRECON_TYPE_MFAMG)
            flags |= LIS_NULLSPACE_CONFIG_NEAR_COMBINATION;
        if(precision!=LIS_PRECISION_DOUBLE)
            flags |= LIS_NULLSPACE_CONFIG_NEAR_PRECISION;
    }

    *issues=flags;
    return LIS_SUCCESS;
}
