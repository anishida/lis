#!/usr/bin/env python3
"""Test the actual MPI allocation snippet in lis_matrix_scale.

Usage: python3 test/testmpiscalerealloc.py
       LIS_MPI_SCALE_ASAN=1 python3 test/testmpiscalerealloc.py

Extracts the source statement block and compiles it with a controllable
allocation failure; no MPI runtime, third-party packages or YAML required.
"""
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def extract_mpi_growth(text):
    start_fn = text.index("LIS_INT lis_matrix_scale(")
    start = text.index("#ifdef USE_MPI\n\t\tif( A->np>D->np )", start_fn)
    end = text.index("\n#endif", start) + len("\n#endif")
    return text[start:end]


PRELUDE = r'''
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

typedef double LIS_SCALAR;
typedef int LIS_INT;
typedef struct matrix { LIS_INT np; void *commtable; } *LIS_MATRIX;
typedef struct vector { LIS_INT np; LIS_SCALAR *value; } *LIS_VECTOR;
#define LIS_OUT_OF_MEMORY 3
#define LIS_SUCCESS 0
#define LIS_SETERR_MEM(n) ((void)(n),++error_reports)

static int realloc_calls, force_failure, send_calls, error_reports;
static void *lis_realloc(void *p,size_t size)
{
    realloc_calls++;
    if(force_failure) return NULL;
    return realloc(p,size);
}
static void lis_send_recv(void *table,LIS_SCALAR *value)
{
    (void)table;
    assert(value);
    send_calls++;
    value[0] += 1.0;
}
static LIS_INT run_mpi_growth(LIS_MATRIX A,LIS_VECTOR D)
{
    LIS_SCALAR *d = D->value;
'''

TESTS = r'''
    return LIS_SUCCESS;
}

int main(void)
{
    struct matrix a = {5,NULL};
    struct vector d = {2,NULL};
    LIS_SCALAR *old;
    d.value=(LIS_SCALAR *)malloc(2*sizeof(LIS_SCALAR));
    assert(d.value);
    d.value[0]=2.5;
    d.value[1]=7.25;
    old=d.value;

    force_failure=1;
    assert(run_mpi_growth(&a,&d)==LIS_OUT_OF_MEMORY);
    assert(realloc_calls==1 && send_calls==0 && error_reports==1);
    assert(d.value==old && d.value[0]==2.5 && d.value[1]==7.25);

    force_failure=0;
    assert(run_mpi_growth(&a,&d)==LIS_SUCCESS);
    assert(realloc_calls==2 && send_calls==1);
    assert(d.value && d.value[0]==3.5 && d.value[1]==7.25);
    d.value[4]=9.0;
    assert(d.value[4]==9.0);
    free(d.value);

    a.np=2;
    d.value=(LIS_SCALAR *)malloc(2*sizeof(LIS_SCALAR));
    assert(d.value);
    d.value[0]=1.0;
    d.value[1]=4.0;
    old=d.value;
    assert(run_mpi_growth(&a,&d)==LIS_SUCCESS);
    assert(realloc_calls==2 && send_calls==2);
    assert(d.value==old && d.value[0]==2.0 && d.value[1]==4.0);
    free(d.value);

    puts("MPI_SCALE_REALLOC_TESTS_PASSED");
    return 0;
}
'''


def main():
    source = (Path(sys.argv[1]) if len(sys.argv) > 1 else
              Path(__file__).resolve().parent.parent / "src/matrix/lis_matrix_ops.c")
    body = extract_mpi_growth(source.read_text(encoding="utf-8"))
    with TemporaryDirectory(prefix="lis-mpi-scale-") as tmp:
        cfile = Path(tmp) / "mpi_realloc.c"
        executable = Path(tmp) / "mpi_realloc"
        cfile.write_text(PRELUDE + body + TESTS, encoding="utf-8")
        cmd = [os.environ.get("CC", "cc"), "-std=c99", "-O0", "-g", "-DUSE_MPI",
               "-Wall", "-Wextra", "-Werror"]
        if os.environ.get("LIS_MPI_SCALE_ASAN") == "1":
            cmd += ["-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
        subprocess.run(cmd + [str(cfile), "-o", str(executable)], check=True)
        subprocess.run([str(executable)], check=True)


if __name__ == "__main__":
    main()
