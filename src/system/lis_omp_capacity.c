/* Elastic OpenMP scratch-capacity management. */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <stddef.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lislib.h"
#include "lis_omp_capacity.h"

#undef __FUNC__
#define __FUNC__ "lis_omp_get_thread_capacity"
LIS_INT
lis_omp_get_thread_capacity(
    LIS_INT *capacity)
{
    if(capacity==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"thread-capacity output is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

#ifdef _OPENMP
#pragma omp critical(lis_omp_capacity_resize)
    {
        *capacity = lis_omp_thread_capacity;
    }
#else
    *capacity = 1;
#endif
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_omp_ensure_thread_capacity"
LIS_INT
lis_omp_ensure_thread_capacity(
    LIS_INT num_threads)
{
#ifdef _OPENMP
    LIS_INT err;
#endif

    if(num_threads<1)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"OpenMP thread count must be positive\n");
        return LIS_ERR_ILL_ARG;
    }

#ifndef _OPENMP
    if(num_threads==1) return LIS_SUCCESS;
    LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,"OpenMP support is not enabled\n");
    return LIS_ERR_NOT_IMPLEMENTED;
#else
    if(omp_in_parallel())
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "OpenMP scratch capacity cannot be resized inside a parallel region\n");
        return LIS_ERR_ILL_ARG;
    }

    err = LIS_SUCCESS;

    /*
     * Serialize the capacity check and any reallocation as one operation.
     * This protects the global pointer/capacity pair against concurrent
     * capacity-management calls.  The public API still requires callers not
     * to overlap this operation with an active LIS solve, because solves use
     * lis_vec_tmp without taking this management lock.
     */
#pragma omp critical(lis_omp_capacity_resize)
    {
        if(num_threads>lis_omp_thread_capacity)
        {
            LIS_SCALAR *tmp;
            size_t count;
            size_t bytes;

            if(lis_vec_tmp==NULL)
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "LIS must be initialized before growing OpenMP scratch capacity\n");
                err = LIS_ERR_ILL_ARG;
            }
            else if((size_t)num_threads >
                    ((size_t)-1)/(LIS_VEC_TMP_PADD*sizeof(LIS_SCALAR)))
            {
                LIS_SETERR(LIS_ERR_OUT_OF_MEMORY,
                    "OpenMP scratch capacity size overflow\n");
                err = LIS_ERR_OUT_OF_MEMORY;
            }
            else
            {
                count = (size_t)num_threads*(size_t)LIS_VEC_TMP_PADD;
                bytes = count*sizeof(LIS_SCALAR);

                tmp = (LIS_SCALAR *)lis_realloc(lis_vec_tmp,bytes);
                if(tmp==NULL)
                {
                    LIS_SETERR_MEM(bytes);
                    err = LIS_ERR_OUT_OF_MEMORY;
                }
                else
                {
                    lis_vec_tmp = tmp;
                    lis_omp_thread_capacity = num_threads;
                }
            }
        }
    }

    return err;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_num_threads_elastic"
LIS_INT
lis_solver_set_num_threads_elastic(
    LIS_SOLVER solver,
    LIS_INT num_threads)
{
    LIS_INT err;

    if(solver==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"solver is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    err = lis_omp_ensure_thread_capacity(num_threads);
    if(err) return err;

    return lis_solver_set_num_threads(solver,num_threads);
}
