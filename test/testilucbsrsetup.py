#!/usr/bin/env python3
"""Fault-inject the real ILUC BSR setup and cleanup C statements.

Extracts the exact initialization/allocations/error handler from
lis_precon_create_iluc_bsr(), plus the real BSR allocator. The numerical
factorization body is intentionally not tested here.

    python3 test/testilucbsrsetup.py
    LIS_ILUC_BSR_ASAN=1 python3 test/testilucbsrsetup.py

This is a focused unit test; full LIS integration testing is separate.
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


def extract_setup(source):
    body = extract(source, "lis_precon_create_iluc_bsr")
    cutoff = "\n\tfor(i=0;i<nr;i++) iw[i] = 0;"
    if body.count(cutoff) != 1:
        raise RuntimeError("ILUC BSR setup boundary changed")
    prefix = body[:body.index(cutoff)]
    label = "\nsetup_failed:\n"
    if body.count(label) != 1:
        raise RuntimeError("ILUC BSR failure handler changed")
    tail = body[body.index(label):]
    return prefix + r'''
    /* The numerical factorization is outside this setup regression. */
    precon->L=L;
    precon->U=U;
    precon->WD=D;
    lis_free2(12,tmp,w,z,iw,wc,wl,iz,zc,zl,ptr,index,value);
    return LIS_SUCCESS;
''' + tail


PRELUDE = r'''
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int LIS_INT;
typedef int LIS_Comm;
typedef double LIS_SCALAR;
typedef double LIS_REAL;
typedef struct matrix {
    LIS_INT n,nr,bnr,bnnz;
    struct matrix *L;
} *LIS_MATRIX;
typedef struct ilu {
    LIS_INT bnnz;
    LIS_SCALAR *value;
} *LIS_MATRIX_ILU;
typedef struct diag {
    LIS_INT n;
} *LIS_MATRIX_DIAG;
typedef struct solver {
    LIS_MATRIX A;
    LIS_REAL params[3];
} *LIS_SOLVER;
typedef struct precon {
    LIS_MATRIX_ILU L,U;
    LIS_MATRIX_DIAG WD;
} *LIS_PRECON;

#define LIS_COMM_WORLD 0
#define LIS_OPTIONS_LEN 0
#define LIS_PARAMS_DROP 0
#define LIS_PARAMS_RATE 1
#define LIS_PARAMS_GAMMA 2
#define LIS_SUCCESS 0
#define LIS_OUT_OF_MEMORY 3
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR_MEM(n) ((void)(n))
#define MAX_TRACKED 48

static void *tracked[MAX_TRACKED];
static int fail_at,allocations,live;

static void *lis_malloc(size_t size,const char *tag)
{
    void *p;
    int i;
    (void)tag;
    if(++allocations==fail_at) return NULL;
    p=malloc(size ? size : 1);
    assert(p);
    for(i=0;i<MAX_TRACKED;i++)
        if(!tracked[i]) {
            tracked[i]=p;
            ++live;
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
            --live;
            free(p);
            return;
        }
    fprintf(stderr,"invalid or duplicate free %p\n",p);
    abort();
}

static void lis_free2_impl(LIS_INT n, void *const items[])
{
    LIS_INT i;
    for(i=0;i<n;i++) lis_free(items[i]);
}
#define lis_free2(n,...) do { \
    void *const items[] = {__VA_ARGS__}; \
    assert((sizeof(items)/sizeof(items[0]))==(size_t)(n)); \
    lis_free2_impl(n,items); \
} while(0)

static LIS_INT lis_matrix_ilu_create(LIS_INT n,LIS_INT bs,
                                    LIS_MATRIX_ILU *out)
{
    (void)n;
    (void)bs;
    *out=(LIS_MATRIX_ILU)lis_malloc(sizeof(**out),"mock ILU");
    if(!*out) return LIS_OUT_OF_MEMORY;
    (*out)->value=NULL;
    (*out)->bnnz=2;
    return LIS_SUCCESS;
}
static LIS_INT lis_matrix_ilu_setCR(LIS_MATRIX_ILU A)
{
    LIS_SCALAR *v=(LIS_SCALAR *)lis_malloc(2*sizeof(LIS_SCALAR),"mock CR");
    if(!v) return LIS_OUT_OF_MEMORY;
    A->value=v;
    return LIS_SUCCESS;
}
static LIS_INT lis_matrix_ilu_destroy(LIS_MATRIX_ILU A)
{
    lis_free(A->value);
    lis_free(A);
    return LIS_SUCCESS;
}
static LIS_INT lis_matrix_diag_duplicateM(LIS_MATRIX A,
                                          LIS_MATRIX_DIAG *out)
{
    *out=(LIS_MATRIX_DIAG)lis_malloc(sizeof(**out),"mock diagonal");
    if(!*out) return LIS_OUT_OF_MEMORY;
    (*out)->n=A->n;
    return LIS_SUCCESS;
}
static LIS_INT lis_matrix_diag_destroy(LIS_MATRIX_DIAG D)
{
    lis_free(D);
    return LIS_SUCCESS;
}
'''

TESTS = r'''
static void check_clean(void)
{
    int i;
    assert(live==0);
    for(i=0;i<MAX_TRACKED;i++) assert(tracked[i]==NULL);
}

int main(void)
{
    struct matrix matrix,lower;
    struct solver solver;
    struct precon precon;
    int failure,success_allocations;
    LIS_INT result;

    memset(&matrix,0,sizeof(matrix));
    memset(&lower,0,sizeof(lower));
    matrix.n=4;
    matrix.nr=2;
    matrix.bnr=2;
    matrix.bnnz=2;
    matrix.L=&lower;
    lower.bnnz=2;
    memset(&solver,0,sizeof(solver));
    solver.A=&matrix;
    solver.params[LIS_PARAMS_DROP]=0.0;
    solver.params[LIS_PARAMS_RATE]=1.0;
    solver.params[LIS_PARAMS_GAMMA]=1.0;

    for(failure=1;failure<=17;failure++) {
        check_clean();
        allocations=0;
        fail_at=failure;
        memset(&precon,0,sizeof(precon));
        result=lis_precon_create_iluc_bsr(&solver,&precon);
        assert(result==LIS_OUT_OF_MEMORY);
        assert(allocations==failure);
        assert(precon.L==NULL && precon.U==NULL && precon.WD==NULL);
        check_clean();
    }

    allocations=0;
    fail_at=0;
    memset(&precon,0,sizeof(precon));
    result=lis_precon_create_iluc_bsr(&solver,&precon);
    assert(result==LIS_SUCCESS);
    success_allocations=allocations;
    assert(success_allocations==17);
    assert(precon.L && precon.U && precon.WD);
    assert(live==5);
    lis_matrix_ilu_destroy(precon.L);
    lis_matrix_ilu_destroy(precon.U);
    lis_matrix_diag_destroy(precon.WD);
    check_clean();

    puts("ILUC_BSR_SETUP_OOM_TESTS_PASSED");
    return 0;
}
'''


def main():
    root = (Path(sys.argv[1]) if len(sys.argv)>1 else
            Path(__file__).resolve().parent.parent)
    iluc = (root/"src/precon/lis_precon_iluc.c").read_text(encoding="utf-8")
    bsr = (root/"src/matrix/lis_matrix_bsr.c").read_text(encoding="utf-8")
    source = (PRELUDE + "\n" + extract(bsr,"lis_matrix_malloc_bsr") +
              "\n" + extract_setup(iluc) + "\n" + TESTS)
    with TemporaryDirectory(prefix="lis-iluc-bsr-oom-") as tmp:
        cfile=Path(tmp)/"test_iluc_bsr_oom.c"
        exe=Path(tmp)/"test_iluc_bsr_oom"
        cfile.write_text(source,encoding="utf-8")
        cc=[os.environ.get("CC","cc"),"-std=c99","-O0","-g",
            "-Wall","-Wextra","-Werror",
            "-Wno-unused-variable","-Wno-unused-but-set-variable",
            "-Wno-unused-parameter"]
        if os.environ.get("LIS_ILUC_BSR_ASAN")=="1":
            cc+=["-fsanitize=address,undefined","-fno-omit-frame-pointer"]
        subprocess.run(cc+[str(cfile),"-o",str(exe)],check=True)
        subprocess.run([str(exe)],check=True)


if __name__=="__main__":
    main()
