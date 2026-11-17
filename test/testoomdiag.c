/* Native OOM tests for real diagonal block allocation and duplication. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifndef LIS_TEST_OOM
#error "testoomdiag requires --enable-oom-tests"
#endif

#define CHECK(expr) do { \
    if( !(expr) ) { \
        fprintf(stderr,"testoomdiag:%d: %s failed\n",__LINE__,#expr); \
        return 1; \
    } \
} while(0)

static int test_blocksize(LIS_MATRIX_DIAG D)
{
    size_t i,base=lis_test_oom_live();
    LIS_INT sizes[3]={2,3,1},rc;
    LIS_SCALAR *old_value=D->value;
    LIS_INT old_nr=D->nr,old_bn=D->bn;
    for(i=1;i<=5;i++) {
        lis_test_oom_begin(i);
        rc=lis_matrix_diag_set_blocksize(D,3,sizes);
        lis_test_oom_end();
        CHECK(rc==LIS_OUT_OF_MEMORY && lis_test_oom_calls()==i);
        CHECK(D->value==old_value && D->bns==NULL && D->v_value==NULL);
        CHECK(D->nr==old_nr && D->bn==old_bn);
        CHECK(lis_test_oom_live()==base);
    }
    lis_test_oom_begin(0);
    rc=lis_matrix_diag_set_blocksize(D,3,sizes);
    lis_test_oom_end();
    CHECK(rc==LIS_SUCCESS && lis_test_oom_calls()==5);
    CHECK(D->bns && D->v_value && D->nr==3);
    CHECK(D->value==NULL);
    CHECK(lis_test_oom_live()==base+4);
    for(i=0;i<3;i++) CHECK(D->bns[i]==sizes[i] && D->v_value[i]);
    puts("Diagonal blocksize: five failure paths, success passed");
    return 0;
}

static int test_duplicate(LIS_MATRIX_DIAG D)
{
    size_t i,base=lis_test_oom_live();
    LIS_MATRIX_DIAG out;
    LIS_INT rc;
    /* 1 struct + (in MPI builds) 1 ranges + 2 tables + 3 raw blocks. */
#ifdef USE_MPI
    const size_t failures=7;
#else
    const size_t failures=6;
#endif
    for(i=1;i<=failures;i++) {
        out=(LIS_MATRIX_DIAG)(uintptr_t)1;
        lis_test_oom_begin(i);
        rc=lis_matrix_diag_duplicate(D,&out);
        lis_test_oom_end();
        CHECK(rc==LIS_OUT_OF_MEMORY && out==NULL);
        CHECK(lis_test_oom_calls()==i);
        CHECK(lis_test_oom_live()==base);
    }
    lis_test_oom_begin(0);
    rc=lis_matrix_diag_duplicate(D,&out);
    lis_test_oom_end();
    CHECK(rc==LIS_SUCCESS && out && out->bns && out->v_value);
    CHECK(lis_test_oom_calls()==failures);
    CHECK(lis_matrix_diag_destroy(out)==LIS_SUCCESS);
    CHECK(lis_test_oom_live()==base);
    puts("Diagonal duplication: injected failures and success passed");
    return 0;
}

static int test_duplicate_matrix(void)
{
    LIS_INT rows[4]={0,2,5,6},ranges[2]={0,6};
    LIS_MATRIX A;
    LIS_MATRIX_DIAG out;
    LIS_INT rc;
    size_t i,base;
#ifdef USE_MPI
    const size_t vbr_allocs=7,block_allocs=3;
#else
    const size_t vbr_allocs=6,block_allocs=2;
#endif
    A=(LIS_MATRIX)lis_malloc(sizeof(*A),"native OOM test matrix");
    CHECK(A!=NULL);
    memset(A,0,sizeof(*A));
    A->status=LIS_MATRIX_VBR;
    A->matrix_type=LIS_MATRIX_VBR;
    A->n=6;
    A->np=6;
    A->nr=3;
    A->nprocs=1;
    A->row=rows;
    A->ranges=ranges;
    base=lis_test_oom_live();
    for(i=1;i<=vbr_allocs;i++) {
        out=(LIS_MATRIX_DIAG)(uintptr_t)1;
        lis_test_oom_begin(i);
        rc=lis_matrix_diag_duplicateM(A,&out);
        lis_test_oom_end();
        CHECK(rc==LIS_OUT_OF_MEMORY && out==NULL);
        CHECK(lis_test_oom_calls()==i && lis_test_oom_live()==base);
    }
    CHECK(lis_matrix_diag_duplicateM(A,&out)==LIS_SUCCESS);
    CHECK(out && out->bns && out->v_value);
    CHECK(lis_matrix_diag_destroy(out)==LIS_SUCCESS);
    CHECK(lis_test_oom_live()==base);
    A->bnr=2;
    A->bnc=2;
    A->matrix_type=LIS_MATRIX_BSR;
    for(i=0;i<2;i++) {
        size_t j;
        if(i==1) A->matrix_type=LIS_MATRIX_BSC;
        A->status=A->matrix_type;
        for(j=1;j<=block_allocs;j++) {
            out=(LIS_MATRIX_DIAG)(uintptr_t)1;
            lis_test_oom_begin(j);
            rc=lis_matrix_diag_duplicateM(A,&out);
            lis_test_oom_end();
            CHECK(rc==LIS_OUT_OF_MEMORY && out==NULL);
            CHECK(lis_test_oom_calls()==j && lis_test_oom_live()==base);
        }
        CHECK(lis_matrix_diag_duplicateM(A,&out)==LIS_SUCCESS);
        CHECK(lis_matrix_diag_destroy(out)==LIS_SUCCESS);
        CHECK(lis_test_oom_live()==base);
    }
    lis_free(A); /* rows and ranges belong to this stack fixture. */
    puts("Matrix diagonal duplication: VBR/BSR/BSC failures passed");
    return 0;
}

static int test_fixed_transition(LIS_MATRIX_DIAG D)
{
    size_t base=lis_test_oom_live();
    LIS_INT old_nr=D->nr,old_bn=D->bn,rc;
    LIS_INT *old_bns=D->bns;
    LIS_SCALAR **old_blocks=D->v_value;
    lis_test_oom_begin(1);
    rc=lis_matrix_diag_set_blocksize(D,2,NULL);
    lis_test_oom_end();
    CHECK(rc==LIS_OUT_OF_MEMORY);
    CHECK(D->nr==old_nr && D->bn==old_bn);
    CHECK(D->bns==old_bns && D->v_value==old_blocks);
    CHECK(lis_test_oom_live()==base);
    CHECK(lis_matrix_diag_set_blocksize(D,2,NULL)==LIS_SUCCESS);
    CHECK(D->bns==NULL && D->v_value==NULL && D->value!=NULL);
    CHECK(lis_test_oom_live()<base);
    puts("Diagonal transition: failed and successful variable-to-fixed passed");
    return 0;
}

int main(int argc,char **argv)
{
    LIS_MATRIX_DIAG D;
    size_t base;
    if( lis_initialize(&argc,&argv)!=LIS_SUCCESS ) return 1;
    base=lis_test_oom_live();
    CHECK(lis_matrix_diag_create(6,0,LIS_COMM_WORLD,&D)==LIS_SUCCESS);
    CHECK(test_blocksize(D)==0);
    CHECK(test_duplicate(D)==0);
    CHECK(test_duplicate_matrix()==0);
    CHECK(test_fixed_transition(D)==0);
    CHECK(lis_matrix_diag_destroy(D)==LIS_SUCCESS);
    CHECK(lis_test_oom_live()==base);
    puts("NATIVE_DIAGONAL_OOM_TESTS_PASSED");
    return lis_finalize()!=LIS_SUCCESS;
}
