#!/usr/bin/env python3
"""Allocation failure tests for the real lis_matrix_diag_set_blocksize().

Run from the LIS checkout:
    python3 test/testdiagsetblocksizeoom.py
    LIS_DIAG_BLOCK_ASAN=1 python3 test/testdiagsetblocksizeoom.py

Extracts the production function from src/matrix/lis_matrix_diag.c and
compiles it with a tracked fault-injecting allocator. Requires only
Python 3 and a C compiler; no configured LIS build or dependencies.
"""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def extract(source):
    start = source.index("LIS_INT lis_matrix_diag_set_blocksize(")
    end = source.index("\n#undef __FUNC__", start)
    return source[start:end]


PRELUDE = r'''
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef MOCK_LONG_LONG
typedef long long LIS_INT;
#else
typedef int LIS_INT;
#endif
typedef double LIS_SCALAR;
typedef struct LIS_MATRIX_DIAG_STRUCT {
    LIS_INT n,bn,nr;
    LIS_SCALAR *value,**v_value;
    LIS_INT *bns;
} *LIS_MATRIX_DIAG;

#define LIS_SUCCESS 0
#define LIS_OUT_OF_MEMORY 3
#define LIS_ERR_ILL_ARG 1
#define LIS_MATRIX_CHECK_NULL 1
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR(code,msg) ((void)(code),(void)(msg))
#define LIS_SETERR_MEM(n) ((void)(n))
#define _max(a,b) ((a)>(b)?(a):(b))

#define MAX_TRACKED 32
static void *tracked[MAX_TRACKED];
static int allocation_calls, fail_at, live_count;
static void *tracked_malloc(size_t bytes)
{
    int i;
    void *p;
    ++allocation_calls;
    if (allocation_calls==fail_at) return NULL;
    p=malloc(bytes ? bytes : 1);
    assert(p);
    for(i=0;i<MAX_TRACKED;i++)
        if(!tracked[i]) {
            tracked[i]=p;
            live_count++;
            return p;
        }
    abort();
}
static void tracked_free(void *p)
{
    int i;
    if(!p) return;
    for(i=0;i<MAX_TRACKED;i++)
        if(tracked[i]==p) {
            tracked[i]=NULL;
            live_count--;
            free(p);
            return;
        }
    fprintf(stderr,"invalid or duplicate free %p\n",p);
    abort();
}
static void *lis_malloc(size_t bytes,char *tag)
{
    (void)tag;
    return tracked_malloc(bytes);
}
static void lis_free(void *p) { tracked_free(p); }
static LIS_INT lis_matrix_diag_check(LIS_MATRIX_DIAG D,LIS_INT flag)
{
    (void)flag;
    return D ? LIS_SUCCESS : LIS_ERR_ILL_ARG;
}

/* Intercept only raw block malloc/free in the extracted LIS function. */
#define malloc(bytes) tracked_malloc(bytes)
#define free(ptr) tracked_free(ptr)
'''

TESTS = r'''
static void begin_case(LIS_MATRIX_DIAG D)
{
    memset(D,0,sizeof(*D));
    assert(live_count==0);
    D->n=6;
    D->bn=1;
    D->nr=6;
    fail_at=0;
    allocation_calls=0;
    D->value=(LIS_SCALAR *)lis_malloc(6*sizeof(LIS_SCALAR),"existing");
    assert(D->value);
    D->value[0]=8.25;
    assert(live_count==1);
    allocation_calls=0;
}

static void assert_intact(LIS_MATRIX_DIAG D,LIS_SCALAR *old)
{
    assert(D->value==old && D->value[0]==8.25);
    assert(!D->bns && !D->v_value);
    assert(D->bn==1 && D->nr==6);
    assert(live_count==1);
}

static void dispose(LIS_MATRIX_DIAG D)
{
    int i;
    if(D->v_value) {
        for(i=0;i<D->nr;i++) {
            if(D->v_value[i]) tracked_free(D->v_value[i]);
        }
        lis_free(D->v_value);
    }
    lis_free(D->bns);
    lis_free(D->value);
    memset(D,0,sizeof(*D));
    assert(live_count==0);
}

static void test_variable_block(void)
{
    struct LIS_MATRIX_DIAG_STRUCT D;
    LIS_INT sizes[3]={2,3,1};
    LIS_SCALAR *old;
    int failure,rc,i;
    for(failure=1;failure<=5;failure++) {
        begin_case(&D);
        old=D.value;
        fail_at=failure;
        rc=lis_matrix_diag_set_blocksize(&D,3,sizes);
        assert(rc==LIS_OUT_OF_MEMORY);
        assert(allocation_calls==failure);
        assert_intact(&D,old);
        dispose(&D);
    }

    begin_case(&D);
    rc=lis_matrix_diag_set_blocksize(&D,3,sizes);
    assert(rc==LIS_SUCCESS);
    assert(D.value==NULL && D.bns && D.v_value);
    assert(D.bn==3 && D.nr==3 && live_count==5);
    for(i=0;i<3;i++) {
        assert(D.bns[i]==sizes[i]);
        assert(D.v_value[i]);
        D.v_value[i][0]=(LIS_SCALAR)(i+1);
        assert(D.v_value[i][0]==(LIS_SCALAR)(i+1));
    }
    /* Reconfiguring an already variable-block diagonal remains a no-op. */
    assert(lis_matrix_diag_set_blocksize(&D,3,sizes)==LIS_SUCCESS);
    assert(live_count==5);
    {
        LIS_INT *old_bns=D.bns;
        LIS_SCALAR **old_table=D.v_value;
        LIS_SCALAR *old_block=D.v_value[0];
        /* Failing to switch representations leaves every old block intact. */
        fail_at=allocation_calls+1;
        assert(lis_matrix_diag_set_blocksize(&D,2,NULL)==LIS_OUT_OF_MEMORY);
        assert(D.bns==old_bns && D.v_value==old_table);
        assert(D.v_value[0]==old_block && D.v_value[0][0]==1.0);
        assert(D.value==NULL && D.bn==3 && D.nr==3 && live_count==5);
        fail_at=0;
        assert(lis_matrix_diag_set_blocksize(&D,2,NULL)==LIS_SUCCESS);
        assert(D.bns==NULL && D.v_value==NULL);
        assert(D.value && D.bn==2 && D.nr==3 && live_count==1);
    }
    dispose(&D);
    puts("variable-block: 5 failures, success and fixed-block switch passed");
}

static void test_fixed_block(void)
{
    struct LIS_MATRIX_DIAG_STRUCT D;
    LIS_SCALAR *old;
    begin_case(&D);
    old=D.value;
    fail_at=1;
    assert(lis_matrix_diag_set_blocksize(&D,2,NULL)==LIS_OUT_OF_MEMORY);
    assert_intact(&D,old);
    fail_at=0;
    assert(lis_matrix_diag_set_blocksize(&D,2,NULL)==LIS_SUCCESS);
    assert(D.value!=NULL && D.bn==2 && D.nr==3);
    assert(live_count==1);
    dispose(&D);
    puts("fixed-block: failure and success passed");
}

static void test_invalid_sizes(void)
{
    struct LIS_MATRIX_DIAG_STRUCT D;
    LIS_SCALAR *old;
    LIS_INT bad[3]={2,-1,5};
    begin_case(&D);
    old=D.value;
    assert(lis_matrix_diag_set_blocksize(&D,0,NULL)==LIS_ERR_ILL_ARG);
    assert(lis_matrix_diag_set_blocksize(&D,0,bad)==LIS_ERR_ILL_ARG);
    assert(lis_matrix_diag_set_blocksize(&D,3,bad)==LIS_ERR_ILL_ARG);
    assert(allocation_calls==0);
    assert_intact(&D,old);
    dispose(&D);
    puts("invalid sizes: preserved state");
}

int main(void)
{
    test_variable_block();
    test_fixed_block();
    test_invalid_sizes();
    puts("DIAG_SET_BLOCKSIZE_OOM_TESTS_PASSED");
    return 0;
}
'''


def main():
    src = (Path(sys.argv[1]) if len(sys.argv)>1 else
           Path(__file__).resolve().parent.parent / "src/matrix/lis_matrix_diag.c")
    func = extract(src.read_text(encoding="utf-8"))
    with TemporaryDirectory(prefix="lis-diag-set-block-") as tmp:
        c_file=Path(tmp)/"test_diag_set_block.c"
        c_file.write_text(PRELUDE+func+TESTS,encoding="utf-8")
        for mode,extra in (("index32", []),("index64", ["-DMOCK_LONG_LONG"])):
            exe=Path(tmp)/("test_diag_set_block_"+mode)
            flags=[os.environ.get("CC","cc"),"-std=c99","-O0","-g",
                   "-Wall","-Wextra","-Werror"]+extra
            if os.environ.get("LIS_DIAG_BLOCK_ASAN")=="1":
                flags+=["-fsanitize=address,undefined","-fno-omit-frame-pointer"]
            subprocess.run(flags+[str(c_file),"-o",str(exe)],check=True)
            subprocess.run([str(exe)],check=True)


if __name__=="__main__":
    main()
