#!/usr/bin/env python3
"""Fault-injection tests for the real LIS ILU creation and VR setup routines.

Run: python3 test/testilucreateoom.py
     LIS_ILU_ASAN=1 python3 test/testilucreateoom.py
     LIS_ILU_OPENMP=1 python3 test/testilucreateoom.py

The two functions are extracted directly from src/matrix/lis_matrix_ilu.c.
No configured LIS build or third-party packages are required.
"""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def extract(source, name):
    start = source.index("LIS_INT " + name + "(")
    end = source.index("\n#undef __FUNC__", start)
    return source[start:end]


PRELUDE = r'''
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int LIS_INT;
typedef double LIS_SCALAR;
typedef struct LIS_MATRIX_ILU_STRUCT {
    LIS_INT n, bs;
    LIS_INT *nnz, **index, *nnz_ma, *bsz;
    LIS_SCALAR **value, ***values;
} *LIS_MATRIX_ILU;

#define LIS_SUCCESS 0
#define LIS_OUT_OF_MEMORY 3
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR_MEM(size) ((void)(size))

#define MAX_TRACKED 32
static void *tracked[MAX_TRACKED];
static int fail_call, allocation_calls, live_count;

static void *lis_malloc(size_t size, char *tag)
{
    int i;
    void *p;
    (void)tag;
    if (++allocation_calls==fail_call) return NULL;
    p=malloc(size ? size : 1);
    assert(p);
    for(i=0;i<MAX_TRACKED;i++) {
        if(!tracked[i]) {
            tracked[i]=p;
            live_count++;
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
            live_count--;
            free(p);
            return;
        }
    }
    fprintf(stderr,"Invalid or double free %p\n",p);
    abort();
}

static void reset(void)
{
    int i;
    assert(live_count==0);
    for(i=0;i<MAX_TRACKED;i++) assert(!tracked[i]);
    fail_call=0;
    allocation_calls=0;
}
'''

TESTS = r'''
static void release_ilu(LIS_MATRIX_ILU A)
{
    lis_free(A->bsz);
    lis_free(A->values);
    lis_free(A->nnz);
    lis_free(A->index);
    lis_free(A);
}

static void test_create(void)
{
    int failure,rc,i;
    LIS_MATRIX_ILU A;

    for(failure=1;failure<=3;failure++) {
        reset();
        A=(LIS_MATRIX_ILU)(size_t)1;
        fail_call=failure;
        rc=lis_matrix_ilu_create(7,2,&A);
        assert(rc==LIS_OUT_OF_MEMORY);
        assert(allocation_calls==failure);
        assert(A==NULL && live_count==0);
    }

    reset();
    rc=lis_matrix_ilu_create(7,2,&A);
    assert(rc==LIS_SUCCESS && A && live_count==3);
    assert(A->n==7 && A->bs==2);
    assert(A->nnz && A->index);
    assert(!A->nnz_ma && !A->value && !A->values && !A->bsz);
    for(i=0;i<7;i++) assert(A->nnz[i]==0 && A->index[i]==NULL);
    release_ilu(A);
    assert(live_count==0);
    puts("ILU create: 3 failure paths and success passed");
}

static void test_vr(void)
{
    int failure,rc;
    LIS_MATRIX_ILU A;

    for(failure=1;failure<=2;failure++) {
        reset();
        rc=lis_matrix_ilu_create(7,2,&A);
        assert(rc==LIS_SUCCESS && live_count==3);
        fail_call=allocation_calls+failure;
        rc=lis_matrix_ilu_setVR(A);
        assert(rc==LIS_OUT_OF_MEMORY);
        assert(A->bsz==NULL && A->values==NULL);
        assert(live_count==3);

        fail_call=0;
        rc=lis_matrix_ilu_setVR(A);
        assert(rc==LIS_SUCCESS);
        assert(A->bsz && A->values && live_count==5);
        release_ilu(A);
        assert(live_count==0);
    }
    puts("ILU VR setup: 2 failure paths and retry passed");
}

int main(void)
{
    test_create();
    test_vr();
    puts("ILU_CREATE_OOM_TESTS_PASSED");
    return 0;
}
'''


def main():
    path = (Path(sys.argv[1]) if len(sys.argv)>1 else
            Path(__file__).resolve().parent.parent/"src/matrix/lis_matrix_ilu.c")
    text = path.read_text(encoding="utf-8")
    bodies = "\n".join(extract(text, name) for name in (
        "lis_matrix_ilu_create", "lis_matrix_ilu_setVR"))
    with TemporaryDirectory(prefix="lis-ilu-oom-") as tmp:
        src=Path(tmp)/"ilu_oom.c"
        exe=Path(tmp)/"ilu_oom"
        src.write_text(PRELUDE+bodies+TESTS,encoding="utf-8")
        cc=[os.environ.get("CC","cc"),"-std=c99","-O0","-g",
            "-Wall","-Wextra","-Werror"]
        if os.environ.get("LIS_ILU_ASAN")=="1":
            cc+=["-fsanitize=address,undefined","-fno-omit-frame-pointer"]
        if os.environ.get("LIS_ILU_OPENMP")=="1":
            cc+=["-fopenmp"]
        subprocess.run(cc+[str(src),"-o",str(exe)],check=True)
        subprocess.run([str(exe)],check=True)


if __name__=="__main__":
    main()
