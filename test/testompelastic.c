#include "lis.h"
#include "lis_omp_capacity.h"
#include <stdio.h>
#include <math.h>
#ifdef _COMPLEX
#include <complex.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

static double scalar_error(LIS_SCALAR value)
{
#ifdef _COMPLEX
    return (double)cabs(value);
#else
    return fabs((double)value);
#endif
}

int main(int argc, char **argv)
{
    LIS_SOLVER solver = NULL;
    LIS_VECTOR x = NULL, y = NULL;
    LIS_INT err, capacity, initial_capacity, target, threads;
    LIS_SCALAR dot;
    int failed = 0;
#ifdef _OPENMP
    LIS_INT parallel_err;
#endif

#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(1);
#endif

    err = lis_initialize(&argc,&argv);
    if(err) return 1;

    if(lis_omp_get_thread_capacity(NULL)!=LIS_ERR_ILL_ARG)
    {
        fprintf(stderr,"NULL capacity output was not rejected\n");
        failed = 1;
    }
    if(lis_omp_ensure_thread_capacity(0)!=LIS_ERR_ILL_ARG ||
       lis_omp_ensure_thread_capacity(-1)!=LIS_ERR_ILL_ARG)
    {
        fprintf(stderr,"invalid thread counts were not rejected\n");
        failed = 1;
    }
    if(lis_solver_set_num_threads_elastic(NULL,1)!=LIS_ERR_ILL_ARG)
    {
        fprintf(stderr,"NULL solver was not rejected\n");
        failed = 1;
    }

    err = lis_omp_get_thread_capacity(&capacity);
    if(err) failed = 1;
    initial_capacity = capacity;

#ifdef _OPENMP
    if(initial_capacity<1)
    {
        fprintf(stderr,"invalid initial capacity=%d\n",(int)initial_capacity);
        failed = 1;
    }

    target = initial_capacity < 4 ? 4 : initial_capacity + 1;

    if(!failed) err = lis_omp_ensure_thread_capacity(target);
    if(err) failed = 1;

    if(!failed) err = lis_omp_get_thread_capacity(&capacity);
    if(err || capacity!=target)
    {
        fprintf(stderr,"capacity did not grow exactly: err=%d initial=%d target=%d capacity=%d\n",
                (int)err,(int)initial_capacity,(int)target,(int)capacity);
        failed = 1;
    }

    /* Capacity is monotonic: smaller requests must not shrink the buffer. */
    if(!failed) err = lis_omp_ensure_thread_capacity(target-1);
    if(err) failed = 1;
    if(!failed) err = lis_omp_get_thread_capacity(&capacity);
    if(err || capacity!=target)
    {
        fprintf(stderr,"capacity unexpectedly shrank: err=%d capacity=%d target=%d\n",
                (int)err,(int)capacity,(int)target);
        failed = 1;
    }

    /* Resizing from an active OpenMP region is explicitly unsupported. */
    parallel_err = LIS_SUCCESS;
#pragma omp parallel num_threads(2) shared(parallel_err)
    {
#pragma omp single
        parallel_err = lis_omp_ensure_thread_capacity(target+1);
    }
    if(parallel_err!=LIS_ERR_ILL_ARG)
    {
        fprintf(stderr,"parallel-region resize was not rejected: err=%d\n",
                (int)parallel_err);
        failed = 1;
    }
    if(!failed) err = lis_omp_get_thread_capacity(&capacity);
    if(err || capacity!=target)
    {
        fprintf(stderr,"rejected resize changed capacity: err=%d capacity=%d target=%d\n",
                (int)err,(int)capacity,(int)target);
        failed = 1;
    }

    if(!failed) err = lis_solver_create(&solver);
    if(err) failed = 1;
    if(!failed) err = lis_solver_set_num_threads_elastic(solver,target);
    if(err) failed = 1;
    if(!failed) err = lis_solver_get_num_threads(solver,&threads);
    if(err || threads!=target)
    {
        fprintf(stderr,"solver elastic thread override mismatch: err=%d threads=%d target=%d\n",
                (int)err,(int)threads,(int)target);
        failed = 1;
    }

    if(!failed) err = lis_vector_create(LIS_COMM_WORLD,&x);
    if(err) failed = 1;
    if(!failed) err = lis_vector_set_size(x,0,100);
    if(err) failed = 1;
    if(!failed) err = lis_vector_duplicate(x,&y);
    if(err) failed = 1;
    if(!failed) err = lis_vector_set_all((LIS_SCALAR)1.0,x);
    if(err) failed = 1;
    if(!failed) err = lis_vector_set_all((LIS_SCALAR)1.0,y);
    if(err) failed = 1;

    if(!failed)
    {
        omp_set_num_threads((int)target);
        err = lis_vector_dot(x,y,&dot);
        if(err || scalar_error(dot-(LIS_SCALAR)100.0)>1.0e-12)
        {
            fprintf(stderr,"elastic-thread reduction failed: err=%d\n",(int)err);
            failed = 1;
        }
    }
#else
    if(capacity!=1)
    {
        fprintf(stderr,"non-OpenMP capacity must be 1, got %d\n",(int)capacity);
        failed = 1;
    }
    if(lis_omp_ensure_thread_capacity(1)!=LIS_SUCCESS ||
       lis_omp_ensure_thread_capacity(2)!=LIS_ERR_NOT_IMPLEMENTED)
    {
        fprintf(stderr,"non-OpenMP capacity behavior mismatch\n");
        failed = 1;
    }

    if(!failed) err = lis_solver_create(&solver);
    if(err) failed = 1;
    if(!failed && lis_solver_set_num_threads_elastic(solver,1)!=LIS_SUCCESS)
    {
        fprintf(stderr,"non-OpenMP elastic single-thread setting failed\n");
        failed = 1;
    }
    if(!failed && lis_solver_set_num_threads_elastic(solver,2)!=LIS_ERR_NOT_IMPLEMENTED)
    {
        fprintf(stderr,"non-OpenMP elastic multi-thread setting was not rejected\n");
        failed = 1;
    }
#endif

    if(y) lis_vector_destroy(y);
    if(x) lis_vector_destroy(x);
    if(solver) lis_solver_destroy(solver);
    lis_finalize();

    if(failed)
    {
        fprintf(stderr,"OMP_ELASTIC_CAPACITY FAILED\n");
        return 1;
    }

    printf("OMP_ELASTIC_CAPACITY PASSED\n");
    return 0;
}
