#!/usr/bin/env python3
"""Deterministic fault injection for the actual MSR/ELL/DIA/JAD allocators.

Run from the LIS source tree:
    python3 test/testothermatrixalloc.py
    LIS_OTHER_ALLOC_ASAN=1 python3 test/testothermatrixalloc.py
    LIS_OTHER_ALLOC_OPENMP=1 python3 test/testothermatrixalloc.py

The test extracts the real LIS C function bodies. It needs only a C
compiler and Python 3, not a configured LIS build or any external library.
"""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

FUNCTIONS = [
    ("src/matrix/lis_matrix_msr.c", "lis_matrix_malloc_msr"),
    ("src/matrix/lis_matrix_ell.c", "lis_matrix_malloc_ell"),
    ("src/matrix/lis_matrix_dia.c", "lis_matrix_malloc_dia"),
    ("src/matrix/lis_matrix_jad.c", "lis_matrix_malloc_jad"),
]


def extract_function(text, name):
    start = text.index("LIS_INT " + name + "(")
    stop = text.index("\n#undef __FUNC__", start)
    return text[start:stop]


PRELUDE = r'''
#include <assert.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef int LIS_INT;
typedef double LIS_SCALAR;
#define LIS_SUCCESS 0
#define LIS_OUT_OF_MEMORY 3
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR_MEM(n) ((void)(n))

#define MAX_TRACKED 32
static void *tracked[MAX_TRACKED];
static int allocation_calls, failure_call, live_count;

static void *lis_malloc(size_t bytes, char *tag)
{
    int i;
    void *ptr;
    (void)tag;
    ++allocation_calls;
    if (allocation_calls == failure_call) return NULL;
    ptr = malloc(bytes ? bytes : 1);
    assert(ptr != NULL);
    for (i=0; i<MAX_TRACKED; i++) {
        if (!tracked[i]) {
            tracked[i] = ptr;
            ++live_count;
            return ptr;
        }
    }
    abort();
}

static void lis_free(void *ptr)
{
    int i;
    if (!ptr) return;
    for (i=0; i<MAX_TRACKED; i++) {
        if (tracked[i] == ptr) {
            tracked[i] = NULL;
            --live_count;
            free(ptr);
            return;
        }
    }
    fprintf(stderr, "Invalid or double free %p\n", ptr);
    abort();
}

static void lis_free2(LIS_INT n, ...)
{
    int i;
    va_list args;
    va_start(args, n);
    for (i=0; i<n; i++) lis_free(va_arg(args, void *));
    va_end(args);
}

static void reset(void)
{
    int i;
    assert(live_count == 0);
    for (i=0; i<MAX_TRACKED; i++) assert(tracked[i] == NULL);
    allocation_calls = 0;
    failure_call = 0;
}
'''

TESTS = r'''
static void check_pair(const char *name, LIS_INT (*allocate)(
    LIS_INT, LIS_INT, LIS_INT **, LIS_SCALAR **),
    LIS_INT a, LIS_INT b)
{
    int fail, result;
    LIS_INT *index;
    LIS_SCALAR *value;
    for (fail=1; fail<=2; fail++) {
        reset();
        failure_call=fail;
        index=(LIS_INT *)(uintptr_t)1;
        value=(LIS_SCALAR *)(uintptr_t)1;
        result=allocate(a,b,&index,&value);
        assert(result==LIS_OUT_OF_MEMORY);
        assert(allocation_calls==fail);
        assert(index==NULL && value==NULL);
        assert(live_count==0);
    }
    reset();
    result=allocate(a,b,&index,&value);
    assert(result==LIS_SUCCESS);
    assert(index!=NULL && value!=NULL);
    assert(live_count==2);
    index[0]=7; value[0]=2.75;
    assert(index[0]==7 && value[0]==2.75);
    lis_free2(2,index,value);
    assert(live_count==0);
    printf("%s: 2 OOM paths and success passed\n",name);
}

static void test_msr(void)
{
    int fail, result;
    LIS_INT *index;
    LIS_SCALAR *value;
    for (fail=1; fail<=2; fail++) {
        reset();
        failure_call=fail;
        index=(LIS_INT *)(uintptr_t)1;
        value=(LIS_SCALAR *)(uintptr_t)1;
        result=lis_matrix_malloc_msr(6,9,6,&index,&value);
        assert(result==LIS_OUT_OF_MEMORY);
        assert(allocation_calls==fail);
        assert(index==NULL && value==NULL);
        assert(live_count==0);
    }
    reset();
    result=lis_matrix_malloc_msr(6,9,6,&index,&value);
    assert(result==LIS_SUCCESS);
    assert(index && value && live_count==2);
    index[0]=3; value[0]=1.25;
    assert(index[0]==3 && value[0]==1.25);
    lis_free2(2,index,value);
    assert(live_count==0);
    puts("MSR: 2 OOM paths and success passed");
}

static void test_jad(void)
{
    int fail, result;
    LIS_INT *perm,*ptr,*index;
    LIS_SCALAR *value;
    for (fail=1; fail<=4; fail++) {
        reset();
        failure_call=fail;
        perm=(LIS_INT *)(uintptr_t)1;
        ptr=(LIS_INT *)(uintptr_t)1;
        index=(LIS_INT *)(uintptr_t)1;
        value=(LIS_SCALAR *)(uintptr_t)1;
        result=lis_matrix_malloc_jad(6,12,3,&perm,&ptr,&index,&value);
        assert(result==LIS_OUT_OF_MEMORY);
        assert(allocation_calls==fail);
        assert(perm==NULL && ptr==NULL && index==NULL && value==NULL);
        assert(live_count==0);
    }
    reset();
    result=lis_matrix_malloc_jad(6,12,3,&perm,&ptr,&index,&value);
    assert(result==LIS_SUCCESS);
    assert(perm && ptr && index && value && live_count==4);
    perm[0]=2; ptr[0]=3; index[0]=4; value[0]=5.5;
    assert(perm[0]==2 && ptr[0]==3 && index[0]==4 && value[0]==5.5);
    lis_free2(4,perm,ptr,index,value);
    assert(live_count==0);
    puts("JAD: 4 OOM paths and success passed");
}

int main(void)
{
    test_msr();
    check_pair("ELL", lis_matrix_malloc_ell, 6,3);
    check_pair("DIA", lis_matrix_malloc_dia, 6,3);
    test_jad();
    puts("OTHER_MATRIX_ALLOC_TESTS_PASSED");
    return 0;
}
'''


def main():
    root = (Path(sys.argv[1]) if len(sys.argv)>1
            else Path(__file__).resolve().parent.parent)
    funcs = []
    for path,name in FUNCTIONS:
        funcs.append(extract_function((root/path).read_text(encoding="utf-8"),name))
    with TemporaryDirectory(prefix="lis-other-matrix-alloc-") as folder:
        source=Path(folder)/"other_matrix_alloc.c"
        binary=Path(folder)/"other_matrix_alloc"
        source.write_text(PRELUDE+"\n".join(funcs)+TESTS,encoding="utf-8")
        flags=[os.environ.get("CC","cc"),"-std=c99","-O0","-g",
               "-Wall","-Wextra","-Werror","-Wno-unused-parameter"]
        if os.environ.get("LIS_OTHER_ALLOC_ASAN") == "1":
            flags+=["-fsanitize=address,undefined","-fno-omit-frame-pointer"]
        if os.environ.get("LIS_OTHER_ALLOC_OPENMP") == "1":
            flags+=["-fopenmp"]
        subprocess.run(flags+[str(source),"-o",str(binary)],check=True)
        subprocess.run([str(binary)],check=True)


if __name__=="__main__":
    main()
