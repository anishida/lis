/* Native OOM regression for the real BSR/BSC/VBR allocators. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <stdint.h>
#include <stdio.h>

#ifndef LIS_TEST_OOM
#error "testoomalloc requires --enable-oom-tests"
#endif

#define CHECK(expr) do { \
    if( !(expr) ) { \
        fprintf(stderr,"testoomalloc:%d: %s failed\n",__LINE__,#expr); \
        return 1; \
    } \
} while(0)

extern LIS_INT lis_matrix_malloc_bsr(LIS_INT,LIS_INT,LIS_INT,LIS_INT,
                                     LIS_INT **,LIS_INT **,LIS_SCALAR **);
extern LIS_INT lis_matrix_malloc_bsc(LIS_INT,LIS_INT,LIS_INT,LIS_INT,
                                     LIS_INT **,LIS_INT **,LIS_SCALAR **);
extern LIS_INT lis_matrix_malloc_vbr(LIS_INT,LIS_INT,LIS_INT,LIS_INT,LIS_INT,
                                     LIS_INT **,LIS_INT **,LIS_INT **,
                                     LIS_INT **,LIS_INT **,LIS_SCALAR **);

typedef LIS_INT (*block_fn)(LIS_INT,LIS_INT,LIS_INT,LIS_INT,
                            LIS_INT **,LIS_INT **,LIS_SCALAR **);

static int check_block(const char *name, block_fn alloc, size_t base)
{
    size_t i;
    LIS_INT *ptr,*index,rc;
    LIS_SCALAR *value;

    for(i=1;i<=3;i++) {
        ptr=(LIS_INT *)(uintptr_t)1;
        index=(LIS_INT *)(uintptr_t)1;
        value=(LIS_SCALAR *)(uintptr_t)1;
        lis_test_oom_begin(i);
        rc=alloc(12,2,3,5,&ptr,&index,&value);
        lis_test_oom_end();
        CHECK(rc==LIS_OUT_OF_MEMORY);
        CHECK(ptr==NULL && index==NULL && value==NULL);
        CHECK(lis_test_oom_calls()==i);
        CHECK(lis_test_oom_live()==base);
    }
    lis_test_oom_begin(0); /* Count all successful allocations. */
    rc=alloc(12,2,3,5,&ptr,&index,&value);
    lis_test_oom_end();
    CHECK(rc==LIS_SUCCESS && ptr && index && value);
    CHECK(lis_test_oom_calls()==3);
    CHECK(lis_test_oom_live()==base+3);
    ptr[0]=1; index[0]=2; value[0]=(LIS_SCALAR)3.5;
    CHECK(ptr[0]==1 && index[0]==2);
    lis_free2(3,ptr,index,value);
    CHECK(lis_test_oom_live()==base);
    printf("%s: three failure paths and success passed\n",name);
    return 0;
}

static int check_vbr(size_t base)
{
    size_t i;
    LIS_INT *row,*col,*ptr,*bptr,*index,rc;
    LIS_SCALAR *value;

    for(i=1;i<=6;i++) {
        row=(LIS_INT *)(uintptr_t)1;
        col=(LIS_INT *)(uintptr_t)1;
        ptr=(LIS_INT *)(uintptr_t)1;
        bptr=(LIS_INT *)(uintptr_t)1;
        index=(LIS_INT *)(uintptr_t)1;
        value=(LIS_SCALAR *)(uintptr_t)1;
        lis_test_oom_begin(i);
        rc=lis_matrix_malloc_vbr(12,48,6,4,5,
                                &row,&col,&ptr,&bptr,&index,&value);
        lis_test_oom_end();
        CHECK(rc==LIS_OUT_OF_MEMORY);
        CHECK(row==NULL && col==NULL && ptr==NULL);
        CHECK(bptr==NULL && index==NULL && value==NULL);
        CHECK(lis_test_oom_calls()==i);
        CHECK(lis_test_oom_live()==base);
    }
    lis_test_oom_begin(0);
    rc=lis_matrix_malloc_vbr(12,48,6,4,5,
                            &row,&col,&ptr,&bptr,&index,&value);
    lis_test_oom_end();
    CHECK(rc==LIS_SUCCESS);
    CHECK(row && col && ptr && bptr && index && value);
    CHECK(lis_test_oom_calls()==6);
    CHECK(lis_test_oom_live()==base+6);
    row[0]=1; col[0]=2; ptr[0]=3; bptr[0]=4; index[0]=5;
    value[0]=(LIS_SCALAR)6;
    CHECK(row[0]==1 && col[0]==2 && ptr[0]==3 && bptr[0]==4 && index[0]==5);
    lis_free2(6,row,col,ptr,bptr,index,value);
    CHECK(lis_test_oom_live()==base);
    puts("VBR: six failure paths and success passed");
    return 0;
}

int main(int argc, char **argv)
{
    size_t base;
    if( lis_initialize(&argc,&argv)!=LIS_SUCCESS ) return 1;
    base=lis_test_oom_live();
    if( check_block("BSR",lis_matrix_malloc_bsr,base)!=0 ) return 1;
    if( check_block("BSC",lis_matrix_malloc_bsc,base)!=0 ) return 1;
    if( check_vbr(base)!=0 ) return 1;
    puts("NATIVE_BLOCK_ALLOC_OOM_TESTS_PASSED");
    return lis_finalize()!=LIS_SUCCESS;
}
