/* Native OOM regression for RCO allocation and partially successful realloc. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_matrix.h"
#include <stdio.h>
#include <string.h>

#ifndef LIS_TEST_OOM
#error "testoomrco requires --enable-oom-tests"
#endif

#define REQUIRE(expr) do { \
    if(!(expr)) { \
        fprintf(stderr,"testoomrco:%d: %s failed\n",__LINE__,#expr); \
        failed=1; goto done; \
    } \
} while(0)

static void release_rows(LIS_INT n,LIS_INT *row,
                         LIS_INT **index,LIS_SCALAR **value)
{
    LIS_INT i;
    for(i=0;i<n;i++)
    {
        if(index && index[i]) lis_free(index[i]);
        if(value && value[i]) lis_free(value[i]);
    }
    if(row) lis_free(row);
    if(index) lis_free(index);
    if(value) lis_free(value);
}

static int test_malloc_oom(void)
{
    LIS_INT counts[2]={1,1};
    LIS_INT *row=NULL,**index=NULL;
    LIS_SCALAR **value=NULL;
    LIS_INT err;
    size_t baseline=lis_test_oom_live(),fail_at,calls;
    int failed=0;

    /* Three outer arrays and two arrays for each of the two rows. */
    for(fail_at=1;fail_at<=7;fail_at++)
    {
        lis_test_oom_begin(fail_at);
        err=lis_matrix_malloc_rco(2,counts,&row,&index,&value);
        lis_test_oom_end();
        calls=lis_test_oom_calls();
        REQUIRE(err==LIS_OUT_OF_MEMORY && calls==fail_at);
        REQUIRE(row==NULL && index==NULL && value==NULL);
        REQUIRE(lis_test_oom_live()==baseline);
    }

    lis_test_oom_begin(0);
    err=lis_matrix_malloc_rco(2,counts,&row,&index,&value);
    lis_test_oom_end();
    REQUIRE(err==LIS_SUCCESS && lis_test_oom_calls()==7);
    REQUIRE(row && index && value && index[0] && value[0] &&
            index[1] && value[1]);

done:
    release_rows(2,row,index,value);
    if(lis_test_oom_live()!=baseline) failed=1;
    return failed;
}

static int test_realloc_oom(void)
{
    LIS_INT counts[1]={1};
    LIS_INT *row=NULL,**index=NULL,*old_index;
    LIS_SCALAR **value=NULL,*old_value;
    unsigned char bytes[sizeof(LIS_SCALAR)];
    LIS_INT err;
    size_t baseline,expected,fail_at;
    int failed=0;

    err=lis_matrix_malloc_rco(1,counts,&row,&index,&value);
    REQUIRE(err==LIS_SUCCESS && index && value &&
            index[0] && value[0]);
    index[0][0]=17;
    memset(&value[0][0],0x5a,sizeof(LIS_SCALAR));
    memcpy(bytes,&value[0][0],sizeof(bytes));
    baseline=lis_test_oom_live();

    for(fail_at=1;fail_at<=2;fail_at++)
    {
        old_index=index[0];
        old_value=value[0];
        lis_test_oom_begin(fail_at);
        err=lis_matrix_realloc_rco(0,4,&index,&value);
        lis_test_oom_end();
        expected=lis_test_oom_calls();

        REQUIRE(err==LIS_OUT_OF_MEMORY && expected==fail_at);
        /* On failure #2 the index realloc committed its new pointer.
         * The original values buffer must remain live and unchanged.
         */
        if(fail_at==1) REQUIRE(index[0]==old_index);
        REQUIRE(value[0]==old_value);
        REQUIRE(index[0] && value[0] && index[0][0]==17);
        REQUIRE(memcmp(&value[0][0],bytes,sizeof(bytes))==0);
        REQUIRE(lis_test_oom_live()==baseline);
    }

    lis_test_oom_begin(0);
    err=lis_matrix_realloc_rco(0,4,&index,&value);
    lis_test_oom_end();
    REQUIRE(err==LIS_SUCCESS && lis_test_oom_calls()==2);
    REQUIRE(index[0][0]==17);
    REQUIRE(memcmp(&value[0][0],bytes,sizeof(bytes))==0);

done:
    release_rows(1,row,index,value);
    if(err==LIS_SUCCESS && !failed && lis_test_oom_live()!=baseline-5)
    {
        /* One row = three outer arrays + two row arrays. */
        fprintf(stderr,"testoomrco: allocation ownership imbalance\n");
        failed=1;
    }
    return failed;
}

int main(int argc,char **argv)
{
    int failed=0;
    size_t baseline;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
    baseline=lis_test_oom_live();
    failed |= test_malloc_oom();
    failed |= test_realloc_oom();
    if(lis_test_oom_live()!=baseline) failed=1;
    if(lis_finalize()!=LIS_SUCCESS) failed=1;
    if(failed) fprintf(stderr,"RCO_NATIVE_OOM=FAILED\n");
    else puts("RCO_NATIVE_OOM=PASSED");
    return failed;
}
