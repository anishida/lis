/* Verify that CSR transpose SpMV survives a failed OpenMP scratch
 * allocation and returns a correct answer instead of dereferencing NULL.
 * Compiles only with -DLIS_TEST_OOM in both LIS and this test program. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <math.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lislib.h"
#ifndef LIS_TEST_OOM
#error "Build this regression with -DLIS_TEST_OOM"
#endif


/* Real and complex off-diagonal entries intentionally differ: the
 * regression checks the conjugate transpose, not merely a diagonal A*x. */
static LIS_SCALAR test_diag(void)
{
#ifdef _COMPLEX
    return (LIS_SCALAR)(4.0+0.25*I);
#else
    return (LIS_SCALAR)4.0;
#endif
}

static LIS_SCALAR test_upper(void)
{
#ifdef _COMPLEX
    return (LIS_SCALAR)(-1.2+0.4*I);
#else
    return (LIS_SCALAR)-1.2;
#endif
}

static LIS_SCALAR test_lower(void)
{
#ifdef _COMPLEX
    return (LIS_SCALAR)(0.3-0.15*I);
#else
    return (LIS_SCALAR)0.3;
#endif
}

static int check_transpose_result(LIS_VECTOR x,LIS_VECTOR y,LIS_INT n)
{
    LIS_INT i;
    for(i=0;i<n;i++)
    {
        LIS_SCALAR expected=conj(test_diag())*x->value[i];
        if(i>0) expected+=conj(test_upper())*x->value[i-1];
        if(i+1<n) expected+=conj(test_lower())*x->value[i+1];
#ifdef _COMPLEX
        if(cabs(y->value[i]-expected)>1e-11)
#else
        if(fabs((double)(y->value[i]-expected))>1e-11)
#endif
        {
            fprintf(stderr,"CSR transpose incorrect at %d\n",(int)i);
            return 1;
        }
    }
    return 0;
}

static int check_case(int split,int threads)
{
    const LIS_INT n=64;
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_INT i,err=0;
    size_t calls=0,normal_calls=0,live_before=0,expected_calls=0;
    int bad=1;
#ifdef _OPENMP
    omp_set_num_threads(threads);
#else
    (void)threads;
#endif
#ifdef _OPENMP
    expected_calls=(size_t)(split || threads!=1 ? 1 : 0);
#else
    /* Without OpenMP both layouts execute allocation-free serial SpMV. */
    expected_calls=0;
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A); if(err)goto done;
    err=lis_matrix_set_size(A,0,n); if(err)goto done;
    for(i=0;i<n;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,test_diag(),A);
        if(err)goto done;
        if(i+1<n)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     test_upper(),A);
            if(err)goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,i+1,i,
                                     test_lower(),A);
            if(err)goto done;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR); if(err)goto done;
    err=lis_matrix_assemble(A); if(err)goto done;
    if(split){err=lis_matrix_split(A); if(err)goto done;}
    err=lis_vector_duplicate(A,&x); if(err)goto done;
    err=lis_vector_duplicate(A,&y); if(err)goto done;
    for(i=0;i<n;i++)
    {
#ifdef _COMPLEX
        x->value[i]=(LIS_SCALAR)(1.0+0.01*i+(0.2+0.005*i)*I);
#else
        x->value[i]=(LIS_SCALAR)(1.0+0.01*i);
#endif
    }

    /* Track the allocation count in the successful path as well:
     * one worker must bypass scratch, while two workers use exactly
     * one scratch allocation and free it before returning. */
    live_before=lis_test_oom_live();
    lis_test_oom_begin(0);
    err=lis_matvech(A,x,y);
    normal_calls=lis_test_oom_calls();
    lis_test_oom_end();
    if(err || normal_calls!=expected_calls ||
       lis_test_oom_live()!=live_before)goto done;
    if(check_transpose_result(x,y,n))goto done;
    lis_test_oom_begin(1);
    err=lis_matvech(A,x,y);
    calls=lis_test_oom_calls();
    lis_test_oom_end();
    if(err || calls!=expected_calls ||
       lis_test_oom_live()!=live_before)goto done;
    if(check_transpose_result(x,y,n))goto done;
    printf("CSR_TRANSPOSE_CASE_OK split=%d threads=%d normal_alloc=%llu oom_alloc=%llu\n",
           split,threads,(unsigned long long)normal_calls,
           (unsigned long long)calls);
    bad=0;
done:
    lis_test_oom_end();
    if(bad) fprintf(stderr,"CSR_TRANSPOSE_OOM FAIL split=%d threads=%d err=%d normal_alloc=%llu oom_alloc=%llu\n",
                    split,threads,(int)err,
                    (unsigned long long)normal_calls,
                    (unsigned long long)calls);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return bad;
}

int main(int argc,char **argv)
{
    int failed;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS)return 2;
#ifdef _OPENMP
    failed=check_case(0,2)|check_case(1,2)|
           check_case(0,1)|check_case(1,1);
#else
    failed=check_case(0,1)|check_case(1,1);
#endif
    lis_finalize();
#ifdef _OPENMP
    if(!failed)puts("CSR_TRANSPOSE_OPENMP_EXECUTED");
    if(!failed)puts("CSR_TRANSPOSE_OOM PASS");
#else
    /* Real serial numerical regression; do not claim an OpenMP/OOM check. */
    if(!failed)puts("CSR_TRANSPOSE_SERIAL_ONLY_PASS");
#endif
    return failed;
}
