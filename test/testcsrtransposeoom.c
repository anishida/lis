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

static int check_case(int split,int threads)
{
    const LIS_INT n=64;
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_INT i,err=0;
    size_t calls=0;
    int bad=1;
#ifdef _OPENMP
    omp_set_num_threads(threads);
#else
    (void)threads;
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A); if(err)goto done;
    err=lis_matrix_set_size(A,0,n); if(err)goto done;
    for(i=0;i<n;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)4.0,A);
        if(err)goto done;
        if(i+1<n)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)-1.2,A);
            if(err)goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,i+1,i,
                                     (LIS_SCALAR)0.3,A);
            if(err)goto done;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR); if(err)goto done;
    err=lis_matrix_assemble(A); if(err)goto done;
    if(split){err=lis_matrix_split(A); if(err)goto done;}
    err=lis_vector_duplicate(A,&x); if(err)goto done;
    err=lis_vector_duplicate(A,&y); if(err)goto done;
    for(i=0;i<n;i++) x->value[i]=(LIS_SCALAR)(1.0+0.01*i);

    /* Normal and allocation-failed calls must return the same y. */
    err=lis_matvech(A,x,y); if(err)goto done;
    for(i=0;i<n;i++)
    {
        double expect=4.0*(double)x->value[i];
        if(i>0) expect-=1.2*(double)x->value[i-1];
        if(i+1<n) expect+=0.3*(double)x->value[i+1];
        if(fabs((double)y->value[i]-expect)>1e-11)goto done;
    }
    lis_test_oom_begin(1);
    err=lis_matvech(A,x,y);
    calls=lis_test_oom_calls();
    lis_test_oom_end();
    if(err || (threads==1 && !split ? calls!=0 : calls==0))goto done;
    for(i=0;i<n;i++)
    {
        double expect=4.0*(double)x->value[i];
        if(i>0) expect-=1.2*(double)x->value[i-1];
        if(i+1<n) expect+=0.3*(double)x->value[i+1];
        if(fabs((double)y->value[i]-expect)>1e-11)
        {
            fprintf(stderr,"CSR OOM fallback wrong entry %d\\n",(int)i);
            goto done;
        }
    }
    bad=0;
done:
    lis_test_oom_end();
    if(bad) fprintf(stderr,"CSR_TRANSPOSE_OOM FAIL split=%d err=%d calls=%llu\\n",
                    split,(int)err,(unsigned long long)calls);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return bad;
}

int main(int argc,char **argv)
{
    int failed;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS)return 2;
    failed=check_case(0,2)|check_case(1,2)|check_case(0,1);
    lis_finalize();
    if(!failed)puts("CSR_TRANSPOSE_OOM PASS");
    return failed;
}
