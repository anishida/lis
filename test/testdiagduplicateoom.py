#!/usr/bin/env python3
"""Fault-injection checks for the real LIS variable-block diagonal duplicators.

Run from a LIS checkout:
    python3 test/testdiagduplicateoom.py
    LIS_DIAG_ASAN=1 python3 test/testdiagduplicateoom.py

The test compiles four function bodies directly from lis_matrix_diag.c.
It does not require a configured LIS build or third-party Python packages.
"""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

FUNCTIONS = (
    "lis_matrix_diag_init",
    "lis_matrix_diag_destroy",
    "lis_matrix_diag_duplicate",
    "lis_matrix_diag_duplicateM",
)


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

#ifdef MOCK_LONG_LONG
typedef long long LIS_INT;
#else
typedef int LIS_INT;
#endif
typedef int LIS_Comm;
typedef double LIS_SCALAR;
typedef struct LIS_MATRIX_DIAG_STRUCT {
    LIS_SCALAR *value, *work, **v_value;
    LIS_INT *bns, *ptr, *ranges;
    LIS_INT bn, nr, n, gn, np, my_rank, nprocs, is, ie, origin;
    LIS_Comm comm;
} *LIS_MATRIX_DIAG;

typedef struct LIS_MATRIX_STRUCT {
    LIS_INT *ranges, *row;
    LIS_INT n, gn, np, nr, bnr, bnc, nprocs, my_rank, is, ie, origin;
    LIS_INT matrix_type;
    LIS_Comm comm;
} *LIS_MATRIX;

#define LIS_SUCCESS 0
#define LIS_OUT_OF_MEMORY 3
#define LIS_MATRIX_CHECK_ALL 1
#define LIS_MATRIX_BSR 7
#define LIS_MATRIX_BSC 8
#define LIS_MATRIX_VBR 9
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR_MEM(n) ((void)(n))
#define _max(a,b) ((a)>(b)?(a):(b))

#define MAX_TRACKED 48
static void *tracked[MAX_TRACKED];
static int alloc_calls, fail_call, live_count;

static void *tracked_malloc(size_t bytes)
{
    int i;
    void *p;
    alloc_calls++;
    if(alloc_calls==fail_call) return NULL;
    p = malloc(bytes ? bytes : 1);
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
    fprintf(stderr,"invalid/double free %p\n",p);
    abort();
}

static void *lis_malloc(size_t bytes, char *tag)
{
    (void)tag;
    return tracked_malloc(bytes);
}
static void lis_free(void *p) { tracked_free(p); }
static void lis_free2(LIS_INT n, void *a, void *b)
{
    assert(n==2);
    lis_free(a);
    lis_free(b);
}
static LIS_INT lis_matrix_diag_check(LIS_MATRIX_DIAG D, LIS_INT mode)
{
    (void)mode;
    return D ? LIS_SUCCESS : 1;
}
static LIS_INT lis_matrix_check(LIS_MATRIX A, LIS_INT mode)
{
    (void)mode;
    return A ? LIS_SUCCESS : 1;
}

/* Route the raw block malloc/free through the same fault injector. */
#define malloc(bytes) tracked_malloc(bytes)
#define free(ptr) tracked_free(ptr)
'''

TESTS = r'''
static void reset(void)
{
    int i;
    assert(live_count==0);
    for(i=0;i<MAX_TRACKED;i++) assert(!tracked[i]);
    alloc_calls=0;
    fail_call=0;
}

static void check_duplicate(LIS_MATRIX_DIAG Din, int vbr, int allocs)
{
    LIS_MATRIX_DIAG out;
    int fail,i;
    for(fail=1;fail<=allocs;fail++) {
        reset();
        fail_call=fail;
        out=(LIS_MATRIX_DIAG)(size_t)1;
        assert(lis_matrix_diag_duplicate(Din,&out)==LIS_OUT_OF_MEMORY);
        assert(alloc_calls==fail);
        assert(out==NULL);
        assert(live_count==0);
    }
    reset();
    out=NULL;
    assert(lis_matrix_diag_duplicate(Din,&out)==LIS_SUCCESS);
    assert(out && out->nr==Din->nr);
    if(vbr) {
        assert(out->bns && out->v_value);
        for(i=0;i<out->nr;i++) {
            assert(out->bns[i]==Din->bns[i]);
            assert(out->v_value[i]);
        }
    } else assert(out->value);
    lis_matrix_diag_destroy(out);
    assert(live_count==0);
}

static void check_duplicate_matrix(LIS_MATRIX Ain, int vbr, int allocs)
{
    LIS_MATRIX_DIAG out;
    int fail,i;
    for(fail=1;fail<=allocs;fail++) {
        reset();
        fail_call=fail;
        out=(LIS_MATRIX_DIAG)(size_t)1;
        assert(lis_matrix_diag_duplicateM(Ain,&out)==LIS_OUT_OF_MEMORY);
        assert(alloc_calls==fail);
        assert(out==NULL);
        assert(live_count==0);
    }
    reset();
    out=NULL;
    assert(lis_matrix_diag_duplicateM(Ain,&out)==LIS_SUCCESS);
    assert(out);
    if(vbr) {
        assert(out->nr==Ain->nr && out->bns && out->v_value);
        for(i=0;i<out->nr;i++) {
            assert(out->bns[i]==Ain->row[i+1]-Ain->row[i]);
            assert(out->v_value[i]);
        }
    } else assert(out->value);
    lis_matrix_diag_destroy(out);
    assert(live_count==0);
}

int main(void)
{
    struct LIS_MATRIX_DIAG_STRUCT diag;
    struct LIS_MATRIX_STRUCT matrix;
    LIS_INT sizes[3]={2,3,1};
    LIS_INT rows[4]={0,2,5,6};
    LIS_INT ranges[2]={0,6};
    memset(&diag,0,sizeof(diag));
    memset(&matrix,0,sizeof(matrix));
    diag.nr=3;
    diag.bn=3;
    diag.n=6;
    diag.nprocs=1;
    diag.bns=sizes;
    diag.ranges=ranges;
    matrix.nr=3;
    matrix.n=6;
    matrix.np=6;
    matrix.nprocs=1;
    matrix.row=rows;
    matrix.ranges=ranges;
    matrix.matrix_type=LIS_MATRIX_VBR;

    check_duplicate(&diag,1,
#ifdef USE_MPI
        7
#else
        6
#endif
    );
    check_duplicate_matrix(&matrix,1,
#ifdef USE_MPI
        7
#else
        6
#endif
    );
    diag.bns=NULL;
    check_duplicate(&diag,0,
#ifdef USE_MPI
        3
#else
        2
#endif
    );
    matrix.matrix_type=LIS_MATRIX_BSR;
    matrix.bnr=2;
    matrix.bnc=2;
    check_duplicate_matrix(&matrix,0,
#ifdef USE_MPI
        3
#else
        2
#endif
    );
    matrix.matrix_type=LIS_MATRIX_BSC;
    check_duplicate_matrix(&matrix,0,
#ifdef USE_MPI
        3
#else
        2
#endif
    );

    puts("DIAG_DUPLICATE_OOM_TESTS_PASSED");
    return 0;
}
'''


def main():
    source = (Path(sys.argv[1]) if len(sys.argv)>1 else
              Path(__file__).resolve().parent.parent /
              "src/matrix/lis_matrix_diag.c")
    text = source.read_text(encoding="utf-8")
    bodies = "\n".join(extract(text, name) for name in FUNCTIONS)
    with TemporaryDirectory(prefix="lis-diag-oom-") as tmp:
        src=Path(tmp)/"diag_oom.c"
        exe=Path(tmp)/"diag_oom"
        src.write_text(PRELUDE+bodies+TESTS,encoding="utf-8")
        for mode,defines in (("serial32",[]),("serial64",["-DMOCK_LONG_LONG"]),
                             ("mpi32",["-DUSE_MPI"]),
                             ("mpi64",["-DUSE_MPI","-DMOCK_LONG_LONG"])):
            executable=Path(tmp)/(mode+"_diag_oom")
            cc=[os.environ.get("CC","cc"),"-std=c99","-O0","-g",
                "-Wall","-Wextra","-Werror","-Wno-unused-parameter",
                "-Wno-unused-but-set-variable"]+defines
            if os.environ.get("LIS_DIAG_ASAN")=="1":
                cc+=["-fsanitize=address,undefined","-fno-omit-frame-pointer"]
            subprocess.run(cc+[str(src),"-o",str(executable)],check=True)
            subprocess.run([str(executable)],check=True)


if __name__=="__main__":
    main()
