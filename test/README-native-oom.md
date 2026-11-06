# Native out-of-memory regression tests

The opt-in native fault injector tests the **compiled LIS implementation**,
not function bodies extracted into temporary C files. Normal production builds
do not contain the injection hooks.

## Build and run (serial)

```sh
autoreconf -fiv
./configure --enable-oom-tests --enable-test \
    --disable-fortran --disable-f90 --disable-mpi --disable-omp
make -j2
make check
```

The opt-in `make check` suite includes:

| C test | Scope |
| --- | --- |
| `testoomalloc.c` | BSR/BSC/VBR allocators, failure at each allocation and success |
| `testoomilu.c` | ILU create, variable-block setup, retry and cleanup |
| `testoomdiag.c` | Variable/fixed block transitions and diagonal duplication including VBR/BSR/BSC |
| `testoomiluc.c` | Native BSR ILUC preconditioner and injected early OOM failures |
| `testoommpiscale.c` | MPI-only diagonal buffer growth helper and failed `lis_realloc` |

For a selected test after building:

```sh
make -C test testoomdiag
./test/testoomdiag
```

## MPI

```sh
./configure MPICC=mpicc --enable-oom-tests --enable-test --enable-mpi \
    --disable-fortran --disable-f90 --disable-omp
make -j2
make -C test testoommpiscale
./test/testoommpiscale
```

This MPI test requires an MPI-capable build, but runs locally with one process
to isolate the allocation-failure path. It does **not** prove correctness of
distributed communication with two or more ranks.

## Injector semantics

The test-only functions are declared in `include/lis_system.h` when
`LIS_TEST_OOM` is enabled.

- `lis_test_oom_begin(N)` resets the call count and causes logical
  allocation number `N` to return `NULL`. With `N == 0`, no allocation
  is failed and calls are counted.
- `lis_test_oom_end()` disables injection without erasing the last count.
- `lis_test_oom_calls()` reports the number of logical attempts while armed.
- `lis_test_oom_live()` counts tracked LIS blocks plus test-wrapped raw
  diagonal blocks; compare it against a baseline established **after**
  persistent test fixtures are initialized.

The injector covers `lis_malloc`, `lis_calloc`, and `lis_realloc`.
Three raw diagonal-block allocation sites in `src/matrix/lis_matrix_diag.c`
are wrapped using `LIS_DIAG_BLOCK_MALLOC/FREE`: outside this test build
these macros expand to the original libc `malloc/free` operations.
The allocator's *own* bookkeeping `malloc` calls are not counted or failed.

**Do not arm the injector concurrently across OpenMP threads.** Its state
and the LIS allocation tracker are process-global and not synchronized.

## Coverage and migration

The previously merged source-extracting Python tests (`testblockallocoom.py`,
`testilucbsrsetup.py`, `testilucreateoom.py`,
`testdiagduplicateoom.py`, `testmpiscalerealloc.py`,
`testdiagsetblocksizeoom.py`) are retained for parity checking.
The native tests cover the real library, but do not reproduce every
configuration from the Python harnesses (including mocked 64-bit MPI paths).
Keep the old scripts until these configurations are independently validated.

Sanitizer builds, multi-rank MPI, alternative index widths and Windows
integration should be run separately before deleting the legacy tests.
