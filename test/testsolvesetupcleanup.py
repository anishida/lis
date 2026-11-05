#!/usr/bin/env python3
"""Exercise the actual lis_solve_setup() source with deterministic failures.

This standalone test extracts the function from src/solver/lis_solver.c,
compiles it with minimal mocked LIS operations and checks that every
allocated temporary vector is released. No LIS build is required.

Run: python3 test/testsolvesetupcleanup.py
"""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def extract_setup(text):
    start = text.index("LIS_INT lis_solve_setup(")
    end = text.index("\n#undef __FUNC__", start)
    return text[start:end]


PRELUDE = r'''
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

typedef int LIS_INT;
typedef void *LIS_MATRIX;
typedef struct vector {
    int magic;
} *LIS_VECTOR;
typedef struct solver {
    LIS_INT setup,retcode;
} *LIS_SOLVER;

#define LIS_TRUE 1
#define LIS_FALSE 0
#define LIS_SUCCESS 0
#define LIS_ERR_ILL_ARG 1
#define LIS_OUT_OF_MEMORY 3
#define LIS_DEBUG_FUNC_IN ((void)0)
#define LIS_DEBUG_FUNC_OUT ((void)0)
#define LIS_SETERR(code,message) ((void)(code), (void)(message))

static int fail_duplicate=0;
static int duplicate_calls=0;
static int live_vectors=0;
static int destroy_calls=0;
static int solve_calls=0;
static int work_cleanup_calls=0;
static LIS_INT solve_result=0;

static LIS_INT lis_vector_duplicate(LIS_MATRIX A,LIS_VECTOR *out)
{
    (void)A;
    *out=NULL;
    duplicate_calls++;
    if(duplicate_calls==fail_duplicate) return LIS_OUT_OF_MEMORY;
    *out=(LIS_VECTOR)malloc(sizeof(**out));
    assert(*out!=NULL);
    (*out)->magic=12345;
    live_vectors++;
    return LIS_SUCCESS;
}

static LIS_INT lis_vector_destroy(LIS_VECTOR v)
{
    if(v) {
        assert(v->magic==12345);
        v->magic=0;
        free(v);
        live_vectors--;
        destroy_calls++;
    }
    return LIS_SUCCESS;
}

static LIS_INT lis_solver_work_destroy(LIS_SOLVER solver)
{
    (void)solver;
    work_cleanup_calls++;
    return LIS_SUCCESS;
}

static LIS_INT lis_solve(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_SOLVER solver)
{
    (void)A;
    assert(b!=NULL && x!=NULL);
    assert(solver->setup==LIS_TRUE);
    solve_calls++;
    return solve_result;
}
'''

TESTS = r'''
static void reset(void)
{
    assert(live_vectors==0);
    fail_duplicate=0;
    duplicate_calls=0;
    destroy_calls=0;
    solve_calls=0;
    work_cleanup_calls=0;
    solve_result=LIS_SUCCESS;
}

int main(void)
{
    struct solver s;
    LIS_MATRIX A=(LIS_MATRIX)&s;
    LIS_INT rc;

    reset();
    s.setup=LIS_FALSE;
    s.retcode=42;
    rc=lis_solve_setup(NULL,&s);
    assert(rc==LIS_ERR_ILL_ARG);
    rc=lis_solve_setup(A,NULL);
    assert(rc==LIS_ERR_ILL_ARG);
    assert(duplicate_calls==0);

    reset();
    s.setup=LIS_FALSE;
    s.retcode=42;
    fail_duplicate=1;
    rc=lis_solve_setup(A,&s);
    assert(rc==LIS_OUT_OF_MEMORY);
    assert(s.retcode==LIS_OUT_OF_MEMORY);
    assert(s.setup==LIS_FALSE);
    assert(duplicate_calls==1 && destroy_calls==0);
    assert(live_vectors==0 && solve_calls==0);

    reset();
    s.setup=LIS_FALSE;
    s.retcode=42;
    fail_duplicate=2;
    rc=lis_solve_setup(A,&s);
    assert(rc==LIS_OUT_OF_MEMORY);
    assert(s.retcode==LIS_OUT_OF_MEMORY);
    assert(s.setup==LIS_FALSE);
    assert(duplicate_calls==2 && destroy_calls==1);
    assert(live_vectors==0 && solve_calls==0);

    reset();
    s.setup=LIS_FALSE;
    s.retcode=42;
    solve_result=LIS_ERR_ILL_ARG;
    rc=lis_solve_setup(A,&s);
    assert(rc==LIS_ERR_ILL_ARG);
    assert(s.retcode==LIS_ERR_ILL_ARG);
    assert(s.setup==LIS_FALSE);
    assert(duplicate_calls==2 && destroy_calls==2);
    assert(work_cleanup_calls==1 && live_vectors==0 && solve_calls==1);

    /* A failed setup must not poison the next setup attempt. */
    reset();
    s.setup=LIS_FALSE;
    s.retcode=42;
    rc=lis_solve_setup(A,&s);
    assert(rc==LIS_SUCCESS);
    assert(s.setup==LIS_TRUE);
    assert(duplicate_calls==2 && destroy_calls==2);
    assert(work_cleanup_calls==0 && live_vectors==0 && solve_calls==1);

    reset();
    s.setup=LIS_TRUE;
    solve_result=LIS_ERR_ILL_ARG;
    rc=lis_solve_setup(A,&s);
    assert(rc==LIS_ERR_ILL_ARG);
    assert(s.setup==LIS_TRUE);
    assert(destroy_calls==2 && live_vectors==0);

    puts("SOLVE_SETUP_CLEANUP_PASSED");
    return 0;
}
'''


def main():
    path = (Path(sys.argv[1]) if len(sys.argv)>1
            else Path(__file__).resolve().parent.parent/"src/solver/lis_solver.c")
    body = extract_setup(path.read_text(encoding="utf-8"))
    with TemporaryDirectory(prefix="lis-setup-cleanup-") as tmp:
        source = Path(tmp)/"setup_cleanup.c"
        executable = Path(tmp)/"setup_cleanup"
        source.write_text(PRELUDE + body + TESTS, encoding="utf-8")
        args = [
            os.environ.get("CC","cc"), "-std=c99", "-O0", "-g",
            "-Wall", "-Wextra", "-Werror",
        ]
        if os.environ.get("LIS_SETUP_ASAN") == "1":
            args += ["-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
        subprocess.run(args + [str(source), "-o", str(executable)], check=True)
        subprocess.run([str(executable)],check=True)


if __name__ == "__main__":
    main()
