/* Generic nullspace basis query helpers. */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include "lislib.h"
#include "lis_nullspace_candidates.h"

#undef __FUNC__
#define __FUNC__ "lis_solver_get_nullspace_dim"
LIS_INT
lis_solver_get_nullspace_dim(
    LIS_SOLVER solver,
    LIS_INT *nvec)
{
    if(solver==NULL || nvec==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid nullspace dimension query\n");
        return LIS_ERR_ILL_ARG;
    }

    *nvec = solver->nullspace_dim;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_copy_nullspace_vector"
LIS_INT
lis_solver_copy_nullspace_vector(
    LIS_SOLVER solver,
    LIS_INT index,
    LIS_VECTOR vector)
{
    if(solver==NULL || vector==NULL ||
       index<0 || index>=solver->nullspace_dim ||
       solver->nullspace==NULL ||
       solver->nullspace[index]==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid nullspace vector query\n");
        return LIS_ERR_ILL_ARG;
    }

    return lis_vector_copy(solver->nullspace[index],vector);
}
