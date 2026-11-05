#!/usr/bin/env python3
"""Fault-injection regression for the actual lis_vector_reuse implementation.

Run from a LIS checkout:
    python3 test/testvectorreuseoom.py
    LIS_VECTOR_REUSE_ASAN=1 python3 test/testvectorreuseoom.py

No LIS build or external Python packages are needed.
"""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def extract(source):
    start = source.index("LIS_INT lis_vector_reuse(")
    end = source.index("\n#undef __FUNC__", start)
    return source[start:end]


PRELUDE = r'''
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

typedef int LIS_INT;
typedef double LIS_SCALAR;
typedef struct vector {
    LIS_INT np, precision, status, is_copy;
    LIS_SCALAR *value, *value_lo, *work;
} *LIS_VECTOR;

#define LIS_TRUE 1
#define LIS_FALSE 0
#define LIS_SUCCESS 0
#define LIS_VECTOR_NULL 0
#define LIS_VECTOR_ASSEMBLED 1
#define LIS_VECTOR_CHECK_NULL 1
#define LIS_OUT_OF_MEMORY 3
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR_MEM(size) ((void)(size))

#define MAX_TRACKED 16
static void *tracked[MAX_TRACKED];
static int calls, fail_call, live_count;

static LIS_INT lis_vector_check(LIS_VECTOR v,LIS_INT what)
{
    (void)what;
    return v ? LIS_SUCCESS : 1;
}

static void *lis_malloc(size_t size,char *tag)
{
    int i;
    void *p;
    (void)tag;
    if(++calls==fail_call) return NULL;
    p=malloc(size ? size : 1);
    assert(p);
    for(i=0;i<MAX_TRACKED;i++)
        if(!tracked[i]) {
            tracked[i]=p;
            ++live_count;
            return p;
        }
    abort();
}

static void lis_free(void *p)
{
    int i;
    if(!p) return;
    for(i=0;i<MAX_TRACKED;i++)
        if(tracked[i]==p) {
            tracked[i]=NULL;
            --live_count;
            free(p);
            return;
        }
    fprintf(stderr,"invalid/double free: %p\n",p);
    abort();
}
'''

TESTS = r'''
static void reset(void)
{
    int i;
    assert(live_count==0);
    for(i=0;i<MAX_TRACKED;i++) assert(!tracked[i]);
    calls=0;
    fail_call=0;
}

static void check_null(LIS_VECTOR v,LIS_INT precision)
{
    assert(v->np==7 && v->precision==precision);
    assert(v->status==LIS_VECTOR_NULL && v->is_copy==LIS_FALSE);
    assert(v->value==NULL && v->value_lo==NULL && v->work==NULL);
}

static void test_reuse(LIS_INT precision)
{
    struct vector storage;
    LIS_VECTOR v;
    LIS_SCALAR *old_value, *old_work;
    int failure, total=precision ? 2 : 1;

    for(failure=1;failure<=total;failure++) {
        reset();
        storage.np=7;
        storage.precision=precision;
        storage.status=LIS_VECTOR_NULL;
        storage.is_copy=LIS_FALSE;
        storage.value=NULL;
        storage.value_lo=NULL;
        storage.work=NULL;
        v=&storage;

        fail_call=failure;
        assert(lis_vector_reuse(&v)==LIS_OUT_OF_MEMORY);
        assert(v==&storage);
        assert(calls==failure);
        assert(live_count==0);
        check_null(v,precision);

        fail_call=0;
        assert(lis_vector_reuse(&v)==LIS_SUCCESS);
        assert(v==&storage);
        assert(v->status==LIS_VECTOR_ASSEMBLED && v->is_copy==LIS_TRUE);
        assert(v->value);
        assert(v->value_lo==(precision ? v->value+8 : NULL));
        assert(!!v->work==!!precision);
        assert(live_count==total);

        old_value=v->value;
        old_work=v->work;
        calls=0;
        assert(lis_vector_reuse(&v)==LIS_SUCCESS);
        assert(calls==0);
        assert(v->value==old_value && v->work==old_work);

        lis_free(v->work);
        lis_free(v->value);
        v->value=NULL;
        v->value_lo=NULL;
        v->work=NULL;
        assert(live_count==0);
    }
    printf("%s precision: all OOM/retry cases passed\n",
           precision ? "extended" : "default");
}

int main(void)
{
    test_reuse(0);
    test_reuse(1);
    puts("VECTOR_REUSE_OOM_PASSED");
    return 0;
}
'''


def main():
    source = (Path(sys.argv[1]) if len(sys.argv)>1 else
              Path(__file__).resolve().parent.parent/"src/vector/lis_vector.c")
    func = extract(source.read_text(encoding="utf-8"))
    with TemporaryDirectory(prefix="lis-vector-reuse-") as tmp:
        cfile=Path(tmp)/"vector_reuse_oom.c"
        binary=Path(tmp)/"vector_reuse_oom"
        cfile.write_text(PRELUDE+func+TESTS,encoding="utf-8")
        cc=[os.environ.get("CC","cc"),"-std=c99","-O0","-g",
            "-Wall","-Wextra","-Werror"]
        if os.environ.get("LIS_VECTOR_REUSE_ASAN")=="1":
            cc+=["-fsanitize=address,undefined","-fno-omit-frame-pointer"]
        subprocess.run(cc+[str(cfile),"-o",str(binary)],check=True)
        subprocess.run([str(binary)],check=True)


if __name__=="__main__":
    main()
