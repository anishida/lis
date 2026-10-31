/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Generic exact-nullspace solution lift construction.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#include <string.h>

#include "lislib.h"
#include "lis_nullspace_lift.h"


static LIS_INT
lis_solver_nullspace_lift_vector_compatible(
    LIS_VECTOR a,
    LIS_VECTOR b)
{
    if(a==NULL || b==NULL)
        return LIS_FALSE;

    return
        a->label==LIS_LABEL_VECTOR &&
        b->label==LIS_LABEL_VECTOR &&
        a->precision==b->precision &&
        a->n==b->n &&
        a->gn==b->gn &&
        a->np==b->np &&
        a->pad==b->pad &&
        a->origin==b->origin &&
        a->my_rank==b->my_rank &&
        a->nprocs==b->nprocs &&
        a->is==b->is &&
        a->ie==b->ie;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_build_nullspace_lift"
LIS_INT
lis_solver_build_nullspace_lift(
    LIS_SOLVER solver,
    LIS_VECTOR rhs,
    LIS_REAL alpha,
    LIS_REAL tolerance,
    LIS_VECTOR lift,
    LIS_NULLSPACE_LIFT_DIAGNOSTICS *diagnostics)
{
    LIS_INT i,err;
    LIS_INT active_count = 0;
    LIS_REAL rhs_norm = 0.0;
    LIS_REAL lift_norm = 0.0;
    LIS_REAL max_relative_component = 0.0;

    LIS_DEBUG_FUNC_IN;

    if(diagnostics)
        memset(diagnostics,0,sizeof(*diagnostics));

    if(solver==NULL || rhs==NULL || lift==NULL || rhs==lift ||
       !isfinite((double)alpha) || alpha<=0.0 ||
       !isfinite((double)tolerance) || tolerance<0.0)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "invalid nullspace lift arguments\n");
        return LIS_ERR_ILL_ARG;
    }

    if(!lis_solver_nullspace_lift_vector_compatible(rhs,lift))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "nullspace lift and RHS layouts are incompatible\n");
        return LIS_ERR_ILL_ARG;
    }

    err = lis_vector_set_all((LIS_SCALAR)0.0,lift);
    if(err)
        return err;

    if(solver->nullspace==NULL || solver->nullspace_dim<=0)
    {
        if(diagnostics)
            diagnostics->nullspace_dim = 0;
        return LIS_SUCCESS;
    }

    if(!lis_solver_nullspace_lift_vector_compatible(
            solver->nullspace[0],rhs))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "RHS layout is incompatible with exact nullspace\n");
        return LIS_ERR_ILL_ARG;
    }

    err = lis_vector_nrm2(rhs,&rhs_norm);
    if(err)
        return err;

    if(!isfinite((double)rhs_norm) || rhs_norm<0.0)
    {
        LIS_SETERR(LIS_BREAKDOWN,
            "nullspace lift RHS norm is non-finite\n");
        return LIS_BREAKDOWN;
    }

    for(i=0;i<solver->nullspace_dim;i++)
    {
        LIS_SCALAR coefficient;
        LIS_REAL magnitude;
        LIS_REAL relative;

        if(!lis_solver_nullspace_lift_vector_compatible(
                solver->nullspace[i],rhs))
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "exact nullspace contains incompatible vector layouts\n");
            return LIS_ERR_ILL_ARG;
        }

        err = lis_vector_dot(
            solver->nullspace[i],rhs,&coefficient);
        if(err)
            return err;

        magnitude = (LIS_REAL)fabs(coefficient);
        if(!isfinite((double)magnitude))
        {
            LIS_SETERR(LIS_BREAKDOWN,
                "nullspace lift coefficient is non-finite\n");
            return LIS_BREAKDOWN;
        }

        relative = rhs_norm>0.0 ? magnitude/rhs_norm : 0.0;
        if(relative>max_relative_component)
            max_relative_component = relative;

        if(rhs_norm>0.0 && relative>tolerance)
        {
            err = lis_vector_axpy(
                coefficient/(LIS_SCALAR)alpha,
                solver->nullspace[i],lift);
            if(err)
                return err;
            active_count++;
        }
    }

    err = lis_vector_nrm2(lift,&lift_norm);
    if(err)
        return err;

    if(!isfinite((double)lift_norm))
    {
        LIS_SETERR(LIS_BREAKDOWN,
            "nullspace lift norm is non-finite\n");
        return LIS_BREAKDOWN;
    }

    if(diagnostics)
    {
        diagnostics->nullspace_dim = solver->nullspace_dim;
        diagnostics->active_count = active_count;
        diagnostics->rhs_norm = rhs_norm;
        diagnostics->max_relative_component = max_relative_component;
        diagnostics->lift_norm = lift_norm;
    }

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}
