#!/usr/bin/env python3
"""Fault-injection tests for the real LIS CSR/COO allocation routines.

Extracts the implementations directly from the two C source files and
compiles them with a deterministic tracked allocator. No LIS build, YAML
workflow, or third-party Python dependency is required.

    python3 test/testcsrcooalloc.py
    LIS_ALLOC_ASAN=1 python3 test/testcsrcooalloc.py
"""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def extract_function(source: str, name: str) -> str:
    signature = "LIS_INT " + name + "("
    start = source.index(signature)
    end = source.index("\n#undef __FUNC__", start)
    return source[start:end]


PRELUDE = r'''
#include <assert.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef int LIS_INT;
typedef double LIS_SCALAR;
#define LIS_SUCCESS 0
#define LIS_OUT_OF_MEMORY 3
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR_MEM(size) ((void)(size))

#define MAX_TRACKED 16
static void *tracked[MAX_TRACKED];
static int fail_call;
static int calls;
static int live_count;

static void *lis_malloc(size_t size, char *tag)
{
    int i;
    void *p;
    (void)tag;
    ++calls;
    if (calls==fail_call) return NULL;
    p=malloc(size ? size : 1);
    assert(p!=NULL);
    for (i=0; i<MAX_TRACKED; ++i) {
        if (!tracked[i]) {
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
    if (!p) return;
    for (i=0; i<MAX_TRACKED; ++i) {
        if (tracked[i]==p) {
            tracked[i]=NULL;
            --live_count;
            free(p);
            return;
        }
    }
    fprintf(stderr,"Invalid free %p\n",p);
    abort();
}

static void lis_free2(LIS_INT n, ...)
{
    va_list ap;
    int i;
    va_start(ap,n);
    for (i=0; i<n; ++i) lis_free(va_arg(ap,void *));
    va_end(ap);
}

static void reset(void)
{
    int i;
    assert(live_count==0);
    for (i=0; i<MAX_TRACKED; ++i) assert(tracked[i]==NULL);
    fail_call=0;
    calls=0;
}
'''

TESTS = r'''
static void test_csr(void)
{
    int failure,rc;
    LIS_INT *ptr,*index;
    LIS_SCALAR *value;

    for (failure=1; failure<=3; ++failure) {
        reset();
        fail_call=failure;
        ptr=(LIS_INT *)(uintptr_t)1;
        index=(LIS_INT *)(uintptr_t)1;
        value=(LIS_SCALAR *)(uintptr_t)1;
        rc=lis_matrix_malloc_csr(4,8,&ptr,&index,&value);
        assert(rc==LIS_OUT_OF_MEMORY);
        assert(calls==failure);
        assert(ptr==NULL && index==NULL && value==NULL);
        assert(live_count==0);
    }

    reset();
    rc=lis_matrix_malloc_csr(4,8,&ptr,&index,&value);
    assert(rc==LIS_SUCCESS && live_count==3);
    assert(ptr && index && value);
    ptr[0]=0;
    index[0]=2;
    value[0]=4.25;
    assert(ptr[0]==0 && index[0]==2 && value[0]==4.25);
    lis_free2(3,ptr,index,value);
    assert(live_count==0);
    puts("csr: three failure points and success passed");
}

static void test_coo(void)
{
    int failure,rc;
    LIS_INT *row,*col;
    LIS_SCALAR *value;

    for (failure=1; failure<=3; ++failure) {
        reset();
        fail_call=failure;
        row=(LIS_INT *)(uintptr_t)1;
        col=(LIS_INT *)(uintptr_t)1;
        value=(LIS_SCALAR *)(uintptr_t)1;
        rc=lis_matrix_malloc_coo(8,&row,&col,&value);
        assert(rc==LIS_OUT_OF_MEMORY);
        assert(calls==failure);
        assert(row==NULL && col==NULL && value==NULL);
        assert(live_count==0);
    }

    reset();
    rc=lis_matrix_malloc_coo(8,&row,&col,&value);
    assert(rc==LIS_SUCCESS && live_count==3);
    assert(row && col && value);
    row[0]=3;
    col[0]=1;
    value[0]=7.5;
    assert(row[0]==3 && col[0]==1 && value[0]==7.5);
    lis_free2(3,row,col,value);
    assert(live_count==0);
    puts("coo: three failure points and success passed");
}

int main(void)
{
    test_csr();
    test_coo();
    puts("CSR_COO_ALLOC_TESTS_PASSED");
    return 0;
}
'''


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    if len(sys.argv) == 2:
        root = Path(sys.argv[1])
    bodies = []
    for path, name in [
        ("src/matrix/lis_matrix_csr.c", "lis_matrix_malloc_csr"),
        ("src/matrix/lis_matrix_coo.c", "lis_matrix_malloc_coo"),
    ]:
        bodies.append(extract_function(
            (root/path).read_text(encoding="utf-8"), name
        ))
    with TemporaryDirectory(prefix="lis-csr-coo-alloc-") as directory:
        src = Path(directory)/"allocation_tests.c"
        exe = Path(directory)/"allocation_tests"
        src.write_text(PRELUDE+"\n".join(bodies)+TESTS, encoding="utf-8")
        cmd = [
            os.environ.get("CC", "cc"),
            "-std=c99", "-O0", "-g", "-Wall", "-Wextra", "-Werror",
        ]
        if os.environ.get("LIS_ALLOC_ASAN") == "1":
            cmd += [
                "-fsanitize=address,undefined", "-fno-omit-frame-pointer",
            ]
        subprocess.run(cmd+[str(src), "-o", str(exe)], check=True)
        subprocess.run([str(exe)], check=True)


if __name__ == "__main__":
    main()
