# Retained-anchor FGMRES reference driver

Opt-in application-side experiment; no new installed API, solver fields, or default solver options. Supports unscaled assembled serial real-double CSR with one active thread and identical DOF ordering across a sequence. MPI, complex and extended precision are not supported by this driver.

The caller lends an immutable anchor matrix P; the driver owns the solver/factors. Each solve starts from zero. A single native FGMRES invocation has a 100-iteration budget and restart length 60. At completed iteration 20, continue in the current Arnoldi cycle only if the finite projected residual ratio r20/r10 <= 0.25 (r10 positive). Early convergence may be accepted. Independently certify ||b-Ax||/||b|| <= 1e-8 (absolute norm for zero RHS); otherwise attempt at most one fresh ILUTP factorization of A and certify again. API errors are not retried. Commit the output only after acceptance. Preserve raw native status separately. Recreate the context after a hard solve error.

Retention policies:
- retain: preserve P and its factor across a fresh fallback;
- low-memory: release the old factor before fallback, lazily rebuild immutable P on the next solve;
- rolling: release the old factor before fallback and promote the new matrix/factor only after successful physical acceptance and output copy. Query current_anchor after every call; keep both matrices alive during solve. The driver never owns matrices.

Low-memory and rolling bound simultaneous live anchor/fallback factors, not process RSS. All modes use the same checkpoint, fixed experimental profile and certification.

## Build and run

Configure a serial real-double static build with OpenMP enabled, then make. For example:

```sh
./configure --disable-mpi --enable-omp --disable-fortran --disable-f90 --disable-saamg --disable-test --disable-shared --enable-static CFLAGS="-O2 -DNDEBUG"
make -j2
cc -O2 -fopenmp -DHAVE_CONFIG_H -Iinclude tools/lis_anchor_pilot.c tools/bench_anchor_pilot.c src/.libs/liblis.a -lm -o bench-anchor-pilot
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE ./bench-anchor-pilot pilot-single-pass A1.mtx b1.mtx A2.mtx b2.mtx
```

Other modes: fresh, pilot-low-memory, pilot-rolling. Inputs are Matrix Market CSR-compatible matrices and one-column real array RHS. Optional LIS_ANCHOR_AUDIT_DIR points to an existing directory for exact round-trip solution output.

The dedicated native regression workflow covers checkpoint equivalence, continuation, fallback, ownership, generation invalidation, and five injected boundary failures. It is not an exhaustive native allocation-failure sweep. Historical experiments and controls remain in the research fork; performance depends on the matrix sequence and is not guaranteed. Original nonlinear contact histories and multi-thread/MPI support require separate validation.
