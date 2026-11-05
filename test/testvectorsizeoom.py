#!/usr/bin/env python3
"""Fault-injection regression for the actual lis_vector_set_size() C function.

Extracts the function from src/vector/lis_vector.c and compiles it with a
tracked allocator. Tests do not need a configured LIS build.

  python3 test/testvectorsizeoom.py
  LIS_VECTOR_ASAN=1 python3 test/testvectorsizeoom.py
"""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def extract_function(source):
    start = source.index("LIS_INT lis_vector_set_size(")
    end = source.index("\n#undef __FUNC__", start)
    return source[start:end]


PRELUDE = r'''
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int LIS_INT;
typedef int LIS_Comm;
typedef double LIS_SCALAR;
typedef struct lis_vector {
    LIS_INT precision, comm;
    LIS_INT *ranges;
    LIS_SCALAR *value, *value_lo, *work;
    LIS_INT is_copy, value_external, status;
    LIS_INT n, gn, np, my_rank, nprocs, is, ie;
} *LIS_VECTOR;

#define LIS_TRUE 1
#define LIS_FALSE 0
#define LIS_VECTOR_NULL 0
#define LIS_VECTOR_ASSEMBLED 1
#define LIS_SUCCESS 0
#define LIS_ERR_ILL_ARG 1
#define LIS_OUT_OF_MEMORY 3
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR2(code,fmt,a,b) ((void)(code),(void)(fmt),(void)(a),(void)(b))
#define LIS_SETERR_MEM(size) ((void)(size))

#define TRACKED_MAX 32
static void *tracked[TRACKED_MAX];
static int fail_call, calls, live_count;

static void *lis_malloc(size_t size, char *tag)
{
    void *p;
    int i;
    (void)tag;
    if (++calls==fail_call) return NULL;
    p=malloc(size ? size : 1);
    assert(p!=NULL);
    for (i=0;i<TRACKED_MAX;i++) {
        if (tracked[i]==NULL) {
            tracked[i]=p;
            ++live_count;
            return p;
        }
    }
    abort();
}

static void lis_free(void *p)
{
    int i;
    if (!p) return;
    for (i=0;i<TRACKED_MAX;i++) {
        if (tracked[i]==p) {
            tracked[i]=NULL;
            --live_count;
            free(p);
            return;
        }
    }
    fprintf(stderr,"invalid or duplicate free: %p\n",p);
    abort();
}

/* Exercise the ranges ownership path used by MPI builds. */
static LIS_INT lis_ranges_create(LIS_Comm comm, LIS_INT *local_n,
    LIS_INT *global_n, LIS_INT **ranges, LIS_INT *is, LIS_INT *ie,
    LIS_INT *nprocs, LIS_INT *my_rank)
{
    LIS_INT *r;
    (void)comm;
    *ranges=NULL;
    r=(LIS_INT *)lis_malloc(2*sizeof(LIS_INT),"mock::ranges");
    if (!r) return LIS_OUT_OF_MEMORY;
    if (*local_n==0) *local_n=*global_n;
    if (*global_n==0) *global_n=*local_n;
    *is=0;
    *ie=*local_n;
    *nprocs=1;
    *my_rank=0;
    r[0]=0;
    r[1]=*local_n;
    *ranges=r;
    return LIS_SUCCESS;
}
'''

TESTS = r'''
static void reset(void)
{
    int i;
    assert(live_count==0);
    for (i=0;i<TRACKED_MAX;i++) assert(tracked[i]==NULL);
    fail_call=0;
    calls=0;
}

static void check_null_state(LIS_VECTOR v, LIS_INT precision)
{
    assert(v->precision==precision);
    assert(v->status==LIS_VECTOR_NULL);
    assert(v->ranges==NULL && v->value==NULL);
    assert(v->value_lo==NULL && v->work==NULL);
    assert(v->n==0 && v->gn==0);
}

static void cleanup(LIS_VECTOR v)
{
    lis_free(v->work);
    lis_free(v->value);
    lis_free(v->ranges);
    memset(v,0,sizeof(*v));
    assert(live_count==0);
}

static void test_precision(LIS_INT precision)
{
    int failure, highest_failure = precision ? 3 : 2;
    struct lis_vector v;
    LIS_INT err,i;

    for (failure=1;failure<=highest_failure;failure++) {
        reset();
        memset(&v,0,sizeof(v));
        v.precision=precision;
        fail_call=failure;
        err=lis_vector_set_size(&v,0,5);
        assert(err==LIS_OUT_OF_MEMORY);
        assert(calls==failure);
        assert(live_count==0);
        check_null_state(&v,precision);

        /* A failed sizing attempt must not poison a later retry. */
        fail_call=0;
        err=lis_vector_set_size(&v,0,5);
        assert(err==LIS_SUCCESS);
        assert(v.status==LIS_VECTOR_ASSEMBLED);
        assert(v.n==5 && v.gn==5 && v.np==5);
        assert(v.ranges && v.ranges[0]==0 && v.ranges[1]==5);
        assert(v.value && v.is_copy==LIS_TRUE);
        assert(v.value_external==LIS_FALSE);
        if(precision) {
            assert(v.value_lo==v.value+6);
            assert(v.work);
        } else {
            assert(v.value_lo==NULL && v.work==NULL);
        }
        for(i=0;i<5;i++) {
            assert(v.value[i]==0.0);
            if(precision) assert(v.value_lo[i]==0.0);
        }
        assert(live_count==(precision ? 3 : 2));
        cleanup(&v);
    }

    reset();
    memset(&v,0,sizeof(v));
    v.precision=precision;
    err=lis_vector_set_size(&v,-1,3);
    assert(err==LIS_ERR_ILL_ARG && calls==0);
    err=lis_vector_set_size(&v,4,3);
    assert(err==LIS_ERR_ILL_ARG && calls==0);
    check_null_state(&v,precision);
    puts(precision ? "quad sizing failures/retries passed" :
                     "default sizing failures/retries passed");
}

int main(void)
{
    test_precision(0);
    test_precision(1);
    puts("VECTOR_SET_SIZE_OOM_TESTS_PASSED");
    return 0;
}
'''


def main():
    src = (Path(sys.argv[1]) if len(sys.argv)==2 else
           Path(__file__).resolve().parent.parent/"src/vector/lis_vector.c")
    body = extract_function(src.read_text(encoding="utf-8"))
    with TemporaryDirectory(prefix="lis-vector-size-oom-") as tmp:
        cfile = Path(tmp)/"vector_size_oom.c"
        executable = Path(tmp)/"vector_size_oom"
        cfile.write_text(PRELUDE+body+TESTS,encoding="utf-8")
        cmd = [os.environ.get("CC","cc"),"-std=c99","-O0","-g",
               "-Wall","-Wextra","-Werror"]
        if os.environ.get("LIS_VECTOR_ASAN") == "1":
            cmd += ["-fsanitize=address,undefined","-fno-omit-frame-pointer"]
        subprocess.run(cmd+[str(cfile),"-o",str(executable)],check=True)
        subprocess.run([str(executable)],check=True)


if __name__ == "__main__":
    main()
