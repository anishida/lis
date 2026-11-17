#!/usr/bin/env python3
"""Fault-injection regression for LIS BSR, BSC and VBR allocators.

Run: python3 test/testblockallocoom.py
     LIS_BLOCK_ALLOC_ASAN=1 python3 test/testblockallocoom.py

Compiles the actual functions extracted from the LIS source tree.
No configured build, YAML, or third-party Python packages are needed.
"""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

FUNCTIONS = [
    ("src/matrix/lis_matrix_bsr.c", "lis_matrix_malloc_bsr"),
    ("src/matrix/lis_matrix_bsc.c", "lis_matrix_malloc_bsc"),
    ("src/matrix/lis_matrix_vbr.c", "lis_matrix_malloc_vbr"),
]


def extract(source, name):
    start = source.index("LIS_INT " + name + "(")
    end = source.index("\n#undef __FUNC__", start)
    return source[start:end]


PRELUDE = r'''
#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef int LIS_INT;
typedef double LIS_SCALAR;
#define LIS_SUCCESS 0
#define LIS_OUT_OF_MEMORY 3
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR_MEM(n) ((void)(n))
#define MAX_TRACKED 32

static void *tracked[MAX_TRACKED];
static int fail_at, calls, live;

static void *lis_malloc(size_t size, char *tag)
{
    int i;
    void *p;
    (void)tag;
    if (++calls==fail_at) return NULL;
    p=malloc(size ? size : 1);
    assert(p);
    for(i=0;i<MAX_TRACKED;i++) {
        if(!tracked[i]) {
            tracked[i]=p;
            live++;
            return p;
        }
    }
    abort();
}

static void lis_free(void *p)
{
    int i;
    if(!p) return;
    for(i=0;i<MAX_TRACKED;i++) {
        if(tracked[i]==p) {
            tracked[i]=NULL;
            live--;
            free(p);
            return;
        }
    }
    fprintf(stderr,"invalid or double free %p\n",p);
    abort();
}

#define lis_free2(n, ...) do { \
    void *items[] = {__VA_ARGS__}; \
    int idx; \
    assert((int)(sizeof(items)/sizeof(items[0]))==(n)); \
    for(idx=0;idx<(n);idx++) lis_free(items[idx]); \
} while(0)

static void reset(void)
{
    int i;
    assert(live==0);
    for(i=0;i<MAX_TRACKED;i++) assert(!tracked[i]);
    calls=0;
    fail_at=0;
}
'''

TESTS = r'''
typedef LIS_INT (*block_alloc_fn)(LIS_INT,LIS_INT,LIS_INT,LIS_INT,
                                  LIS_INT **,LIS_INT **,LIS_SCALAR **);

static void test_block(const char *name, block_alloc_fn fn)
{
    LIS_INT *ptr,*index;
    LIS_SCALAR *value;
    int fail,rc;
    for(fail=1;fail<=3;fail++) {
        reset();
        fail_at=fail;
        ptr=(LIS_INT *)(uintptr_t)1;
        index=(LIS_INT *)(uintptr_t)1;
        value=(LIS_SCALAR *)(uintptr_t)1;
        rc=fn(12,2,3,5,&ptr,&index,&value);
        assert(rc==LIS_OUT_OF_MEMORY);
        assert(calls==fail);
        assert(!ptr && !index && !value && live==0);
    }
    reset();
    rc=fn(12,2,3,5,&ptr,&index,&value);
    assert(rc==LIS_SUCCESS && live==3 && ptr && index && value);
    ptr[0]=1; index[0]=2; value[0]=3.5;
    assert(ptr[0]==1 && index[0]==2 && value[0]==3.5);
    lis_free2(3,ptr,index,value);
    assert(live==0);
    printf("%s: three failure paths and success passed\n",name);
}

static void test_vbr(void)
{
    LIS_INT *row,*col,*ptr,*bptr,*bindex;
    LIS_SCALAR *value;
    int fail,rc;
    for(fail=1;fail<=6;fail++) {
        reset();
        fail_at=fail;
        row=(LIS_INT *)(uintptr_t)1;
        col=(LIS_INT *)(uintptr_t)1;
        ptr=(LIS_INT *)(uintptr_t)1;
        bptr=(LIS_INT *)(uintptr_t)1;
        bindex=(LIS_INT *)(uintptr_t)1;
        value=(LIS_SCALAR *)(uintptr_t)1;
        rc=lis_matrix_malloc_vbr(12,48,6,4,5,
                                 &row,&col,&ptr,&bptr,&bindex,&value);
        assert(rc==LIS_OUT_OF_MEMORY);
        assert(calls==fail);
        assert(!row && !col && !ptr && !bptr && !bindex && !value);
        assert(live==0);
    }
    reset();
    rc=lis_matrix_malloc_vbr(12,48,6,4,5,
                             &row,&col,&ptr,&bptr,&bindex,&value);
    assert(rc==LIS_SUCCESS && live==6);
    assert(row && col && ptr && bptr && bindex && value);
    row[0]=1; col[0]=2; ptr[0]=3; bptr[0]=4;
    bindex[0]=5; value[0]=6.5;
    assert(row[0]==1 && col[0]==2 && ptr[0]==3 && bptr[0]==4 &&
           bindex[0]==5 && value[0]==6.5);
    lis_free2(6,row,col,ptr,bptr,bindex,value);
    assert(live==0);
    puts("VBR: six failure paths and success passed");
}

int main(void)
{
    test_block("BSR",lis_matrix_malloc_bsr);
    test_block("BSC",lis_matrix_malloc_bsc);
    test_vbr();
    puts("BLOCK_MATRIX_ALLOC_OOM_TESTS_PASSED");
    return 0;
}
'''


def main():
    root = (Path(sys.argv[1]) if len(sys.argv)>1
            else Path(__file__).resolve().parent.parent)
    functions = "\n".join(
        extract((root/path).read_text(encoding="utf-8"), name)
        for path,name in FUNCTIONS
    )
    with TemporaryDirectory(prefix="lis-block-oom-") as tmp:
        source=Path(tmp)/"block_oom.c"
        binary=Path(tmp)/"block_oom"
        source.write_text(PRELUDE+functions+TESTS, encoding="utf-8")
        cc=[os.environ.get("CC","cc"),"-std=c99","-O0","-g",
            "-Wall","-Wextra","-Werror","-Wno-unused-parameter"]
        if os.environ.get("LIS_BLOCK_ALLOC_ASAN")=="1":
            cc+=["-fsanitize=address,undefined","-fno-omit-frame-pointer"]
        subprocess.run(cc+[str(source),"-o",str(binary)],check=True)
        subprocess.run([str(binary)],check=True)


if __name__=="__main__":
    main()
