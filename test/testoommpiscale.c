/* Native MPI diagonal-growth OOM regression for the production helper. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <stdio.h>
#include <string.h>

#if !defined(LIS_TEST_OOM) || !defined(USE_MPI)
#error "testoommpiscale requires --enable-oom-tests --enable-mpi"
#endif

#define CHECK(expr) do { \
    if( !(expr) ) { \
        fprintf(stderr,"testoommpiscale:%d: %s failed\n",__LINE__,#expr); \
        return 1; \
    } \
} while(0)

int main(int argc,char **argv)
{
    struct LIS_MATRIX_STRUCT matrix;
    struct LIS_VECTOR_STRUCT diag;
    LIS_SCALAR *old;
    size_t base;
    if( lis_initialize(&argc,&argv)!=LIS_SUCCESS ) return 1;
    memset(&matrix,0,sizeof(matrix));
    memset(&diag,0,sizeof(diag));
    base=lis_test_oom_live();
    matrix.np=5;
    diag.np=2;
    diag.value=(LIS_SCALAR *)lis_malloc(2*sizeof(LIS_SCALAR),"testoommpiscale::value");
    CHECK(diag.value!=NULL);
    diag.value[0]=(LIS_SCALAR)2;
    diag.value[1]=(LIS_SCALAR)7;
    old=diag.value;

    lis_test_oom_begin(1);
    CHECK(lis_matrix_scale_expand_diagonal(&matrix,&diag)==LIS_OUT_OF_MEMORY);
    lis_test_oom_end();
    CHECK(lis_test_oom_calls()==1);
    CHECK(diag.value==old && diag.value[0]==(LIS_SCALAR)2);
    CHECK(diag.value[1]==(LIS_SCALAR)7);
    CHECK(lis_test_oom_live()==base+1);

    lis_test_oom_begin(0);
    CHECK(lis_matrix_scale_expand_diagonal(&matrix,&diag)==LIS_SUCCESS);
    lis_test_oom_end();
    CHECK(lis_test_oom_calls()==1);
    CHECK(diag.value!=NULL && diag.value[0]==(LIS_SCALAR)2);
    CHECK(diag.value[1]==(LIS_SCALAR)7);
    diag.value[4]=(LIS_SCALAR)9;
    CHECK(diag.value[4]==(LIS_SCALAR)9);

    matrix.np=2;
    lis_test_oom_begin(1);
    CHECK(lis_matrix_scale_expand_diagonal(&matrix,&diag)==LIS_SUCCESS);
    lis_test_oom_end();
    CHECK(lis_test_oom_calls()==0);
    lis_free(diag.value);
    CHECK(lis_test_oom_live()==base);
    puts("NATIVE_MPI_SCALE_OOM_TESTS_PASSED");
    return lis_finalize()!=LIS_SUCCESS;
}
