/* Elastic OpenMP scratch-capacity management. */

#ifndef __LIS_OMP_CAPACITY_H__
#define __LIS_OMP_CAPACITY_H__

#include "lis.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * Ensure that LIS global OpenMP reduction scratch can serve at least
 * num_threads threads.  Capacity only grows.
 *
 * The call must be made outside an active OpenMP parallel region and must not
 * overlap an active LIS operation that may use the global reduction scratch.
 * Concurrent capacity-management calls are serialized internally.
 *
 * Calling omp_set_num_threads() directly does not resize LIS scratch.  An
 * embedding application that raises the OpenMP thread count after
 * lis_initialize() must call this function first, or use
 * lis_solver_set_num_threads_elastic().
 */
extern LIS_INT lis_omp_ensure_thread_capacity(LIS_INT num_threads);

/* Return the currently allocated OpenMP scratch capacity. */
extern LIS_INT lis_omp_get_thread_capacity(LIS_INT *capacity);

/*
 * Convenience wrapper for embedded applications: grow global LIS scratch when
 * necessary and then apply the existing per-solver thread override.
 * The same no-active-solve/no-active-parallel-region rule applies.
 */
extern LIS_INT lis_solver_set_num_threads_elastic(
    LIS_SOLVER solver,
    LIS_INT num_threads);

#ifdef __cplusplus
}
#endif

#endif
