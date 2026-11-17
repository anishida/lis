#!/usr/bin/env python3
"""Standalone fault-injection checks for the *actual* LIS RCO allocation functions.

Run from an LIS checkout:
    python3 test/testrcoalloc.py
    LIS_RCO_ASAN=1 python3 test/testrcoalloc.py  # if sanitizer runtime is available

This script extracts two function bodies from src/matrix/lis_matrix_rco.c and
compiles them with an instrumented allocator. It does not require a configured
LIS build, and it never changes files in the checkout.
"""
import os
import pathlib
import subprocess
import sys
import tempfile


def extract(src: str, name: str) -> str:
    signature = "LIS_INT " + name + "("
    start = src.find(signature)
    if start < 0:
        raise RuntimeError(f"missing function {name}")
    stop = src.find("#undef __FUNC__", start)
    if stop < 0:
        raise RuntimeError(f"end marker missing for {name}")
    return src[start:stop].strip() + "\n"


PRELUDE = r'''
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef int LIS_INT;
typedef double LIS_SCALAR;
#define LIS_SUCCESS 0
#define LIS_OUT_OF_MEMORY 2
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR_MEM(bytes) ((void)(bytes))

#define MAX_PTRS 128
static void *tracked[MAX_PTRS];
static int fail_call = 0;
static int calls = 0;

static int slot_for(void *p)
{
    int i;
    for (i=0; i<MAX_PTRS; ++i) if (tracked[i]==p) return i;
    return -1;
}

static int live_count(void)
{
    int i, result=0;
    for (i=0; i<MAX_PTRS; ++i) result += tracked[i]!=NULL;
    return result;
}

static void track(void *p)
{
    int i;
    assert(p!=NULL);
    assert(slot_for(p)<0);
    for (i=0; i<MAX_PTRS; ++i)
        if (!tracked[i]) { tracked[i]=p; return; }
    abort();
}

static void *lis_malloc(size_t size, char *tag)
{
    void *p;
    (void)tag;
    ++calls;
    if (calls==fail_call) return NULL;
    p=malloc(size ? size : 1);
    assert(p!=NULL);
    track(p);
    return p;
}

static void *lis_realloc(void *p, size_t size)
{
    void *q;
    int slot=slot_for(p);
    assert(slot>=0);
    ++calls;
    if (calls==fail_call) return NULL;
    q=realloc(p,size ? size : 1);
    assert(q!=NULL);
    tracked[slot]=q;
    return q;
}

static void lis_free(void *p)
{
    int slot;
    if (!p) return;
    slot=slot_for(p);
    if (slot<0) { fprintf(stderr, "bad/double free: %p\n",p); abort(); }
    tracked[slot]=NULL;
    free(p);
}

static void lis_free2(LIS_INT count, ...)
{
    va_list ap;
    int i;
    va_start(ap,count);
    for (i=0; i<count; ++i) lis_free(va_arg(ap,void *));
    va_end(ap);
}
'''

POSTLUDE = r'''
static void reset_test(void)
{
    assert(live_count()==0);
    fail_call=0;
    calls=0;
}

static void release_rows(int n, LIS_INT *row, LIS_INT **index, LIS_SCALAR **value)
{
    int i;
    for (i=0; i<n; ++i) { lis_free(index[i]); lis_free(value[i]); }
    lis_free2(3,row,index,value);
}

static void test_malloc_rco(void)
{
    int i, result;
    LIS_INT nnz[] = {2,0,3};
    LIS_INT *row;
    LIS_INT **index;
    LIS_SCALAR **value;

    for (i=1; i<=7; ++i) {
        reset_test();
        fail_call=i;
        row=(LIS_INT *)(uintptr_t)1;
        index=(LIS_INT **)(uintptr_t)1;
        value=(LIS_SCALAR **)(uintptr_t)1;
        result=lis_matrix_malloc_rco(3,nnz,&row,&index,&value);
        assert(result==LIS_OUT_OF_MEMORY);
        assert(row==NULL && index==NULL && value==NULL);
        assert(live_count()==0);
    }

    reset_test();
    result=lis_matrix_malloc_rco(3,nnz,&row,&index,&value);
    assert(result==LIS_SUCCESS && live_count()==7);
    assert(row[0]==0 && row[1]==0 && row[2]==0);
    assert(index[0] && value[0] && index[1]==NULL && value[1]==NULL);
    assert(index[2] && value[2]);
    release_rows(3,row,index,value);
    assert(live_count()==0);

    reset_test();
    result=lis_matrix_malloc_rco(3,NULL,&row,&index,&value);
    assert(result==LIS_SUCCESS && live_count()==3);
    for (i=0; i<3; ++i) assert(index[i]==NULL && value[i]==NULL);
    release_rows(3,row,index,value);
    assert(live_count()==0);
    puts("malloc_rco: 7 injected allocation failures + 2 success paths passed");
}

static void test_realloc_rco(void)
{
    int phase, result;
    LIS_INT **index;
    LIS_SCALAR **value;
    for (phase=0; phase<3; ++phase) {
        reset_test();
        index=(LIS_INT **)lis_malloc(sizeof(*index),"test-index");
        value=(LIS_SCALAR **)lis_malloc(sizeof(*value),"test-value");
        index[0]=(LIS_INT *)lis_malloc(sizeof(LIS_INT),"index-row");
        value[0]=(LIS_SCALAR *)lis_malloc(sizeof(LIS_SCALAR),"value-row");
        index[0][0]=17;
        value[0][0]=3.25;
        if (phase<2) fail_call=calls+phase+1;
        result=lis_matrix_realloc_rco(0,4,&index,&value);
        assert(result==(phase<2 ? LIS_OUT_OF_MEMORY : LIS_SUCCESS));
        assert(index[0]!=NULL && value[0]!=NULL);
        assert(index[0][0]==17 && value[0][0]==3.25);
        assert(live_count()==4);
        lis_free(index[0]); lis_free(value[0]);
        lis_free(index); lis_free(value);
        assert(live_count()==0);
    }
    puts("realloc_rco: first/second realloc failures + success passed");
}

int main(void)
{
    test_malloc_rco();
    test_realloc_rco();
    puts("RCO_ALLOC_TESTS_PASSED");
    return 0;
}
'''


def main() -> None:
    source = pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else (
        pathlib.Path(__file__).resolve().parent.parent / "src/matrix/lis_matrix_rco.c"
    )
    original = source.read_text(encoding="utf-8")
    bodies = "\n".join(extract(original, name) for name in (
        "lis_matrix_malloc_rco", "lis_matrix_realloc_rco"
    ))
    with tempfile.TemporaryDirectory(prefix="lis-rco-check-") as tmp:
        cfile = pathlib.Path(tmp) / "rco_faults.c"
        binary = pathlib.Path(tmp) / "rco_faults"
        cfile.write_text(PRELUDE + bodies + POSTLUDE, encoding="utf-8")
        cc = [os.environ.get("CC", "cc"), "-std=c99", "-O0", "-g", "-Wall", "-Wextra", "-Werror"]
        if os.environ.get("LIS_RCO_ASAN") == "1":
            cc += ["-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
        cc += [str(cfile), "-o", str(binary)]
        print("Compiling RCO functions from", source, flush=True)
        subprocess.run(cc,check=True)
        subprocess.run([str(binary)],check=True)


if __name__ == "__main__":
    main()