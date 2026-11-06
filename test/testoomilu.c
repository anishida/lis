/* Native OOM tests exercising LIS ILU constructors and cleanup. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <stdint.h>
#include <stdio.h>

#ifndef LIS_TEST_OOM
#error "testoomilu requires --enable-oom-tests"
#endif

#define CHECK(expr) do { \
    if( !(expr) ) { \
        fprintf(stderr,"testoomilu:%d: %s failed\n",__LINE__,#expr); \
        return 1; \
    } \
} while(0)

static int test_create(size_t base)
{
    size_t i;
    LIS_MATRIX_ILU A;
    LIS_INT rc;
    for(i=1;i<=3;i++) {
        A=(LIS_MATRIX_ILU)(uintptr_t)1;
        lis_test_oom_begin(i);
        rc=lis_matrix_ilu_create(7,2,&A);
        lis_test_oom_end();
        CHECK(rc==LIS_OUT_OF_MEMORY && A==NULL);
        CHECK(lis_test_oom_calls()==i);
        CHECK(lis_test_oom_live()==base);
    }
    lis_test_oom_begin(0);
    rc=lis_matrix_ilu_create(7,2,&A);
    lis_test_oom_end();
    CHECK(rc==LIS_SUCCESS && A && A->nnz && A->index);
    CHECK(lis_test_oom_calls()==3);
    CHECK(A->n==7 && A->bs==2);
    CHECK(lis_matrix_ilu_destroy(A)==LIS_SUCCESS);
    CHECK(lis_test_oom_live()==base);
    puts("ILU constructor: three failure paths and success passed");
    return 0;
}

static int test_vr(size_t base)
{
    size_t i;
    LIS_MATRIX_ILU A;
    LIS_INT rc;
    for(i=1;i<=2;i++) {
        CHECK(lis_matrix_ilu_create(7,2,&A)==LIS_SUCCESS);
        lis_test_oom_begin(i);
        rc=lis_matrix_ilu_setVR(A);
        lis_test_oom_end();
        CHECK(rc==LIS_OUT_OF_MEMORY);
        CHECK(A->bsz==NULL && A->values==NULL);
        CHECK(lis_test_oom_calls()==i);
        CHECK(lis_test_oom_live()==base+3);
        CHECK(lis_matrix_ilu_setVR(A)==LIS_SUCCESS);
        CHECK(A->bsz && A->values);
        CHECK(lis_matrix_ilu_destroy(A)==LIS_SUCCESS);
        CHECK(lis_test_oom_live()==base);
    }
    puts("ILU variable blocks: two failure paths and retry passed");
    return 0;
}

int main(int argc,char **argv)
{
    size_t base;
    if( lis_initialize(&argc,&argv)!=LIS_SUCCESS ) return 1;
    base=lis_test_oom_live();
    if( test_create(base)!=0 ) return 1;
    if( test_vr(base)!=0 ) return 1;
    puts("NATIVE_ILU_OOM_TESTS_PASSED");
    return lis_finalize()!=LIS_SUCCESS;
}
