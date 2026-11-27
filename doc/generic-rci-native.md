# Experimental native Generic RCI

The existing array RCI API is unchanged. The opt-in `lis_rci_native.h`
facade exposes six linear methods, FEAST and Krylov under one typed context.

## Build

Run `./configure --enable-rci --disable-mpi --disable-quad --enable-shared`,
then `make -j2` and `make install`. This candidate requires Python 3,
Node.js and GCC/Clang hidden symbol visibility. Disabled builds and installed
clients do not require the generators. Generation uses pinned local sources
in build directories, without modifying the original source or accessing
the network.

## Contract

Native LIS real/complex scalar, non-MPI Serial/OpenMP; linear unscaled CSR with NONE/Jacobi,
FEAST standard/generalized CSR without retries or conversion. Krylov retains
its fixed-double callback ABI. Check the client/library ABI descriptor before
typed calls. Create a context with API version and request size, start the
appropriate method and execute requests while status is positive. Complete
each request with its token and native callback status. DONE commits linear
and FEAST solution outputs; Krylov retains native partial-result semantics.
Wrong tokens do not consume requests. Keep operands and solver alive and
exclusive; stop the backend before cancellation or destruction. All calls
must use the same thread and OpenMP settings. Shift signs are explicit:
FEAST sigma B-A, Krylov A-sigma B. See test_native_installed.c in
qualification/rci-restore for a complete client.

The native continuation is reconstructed from pinned LIS source, not recovered
RC3 code. This is a candidate, not a stable ABI. The native-only scalar interface passed the installed-client matrix in CI 46004501846 (six scalar configurations). MPI and double-double
solver mode are excluded. Maintainer review and contributor licensing
confirmation remain necessary.

A nonzero external completion code is a callback failure even when its
numeric value overlaps a native MAXITER or BREAKDOWN result. No additional
certification callback is launched after that failure. The public driver
rejects unwritten/nonfinite raw Krylov outputs. Forced linear termination
invalidates recycling state. Installed-client tests cover these distinctions,
wrong tokens, stale generations and recovery after interruption.

## Native scalar configuration

`--enable-rci=native` enables only the typed native API, preserving the existing
real-double array interface's configuration contract. Real double, extended long double, binary128 and their three complex counterparts passed installed-client smoke and fault tests.
Krylov remains fixed-real-double. The complex FGCRO generation applies the
qualified Hermitian projection-order correction to its private build copy.

## Release gate

The final qualification branch runs a same-commit installed-client matrix,
legacy enabled/disabled regression, exact baseline exports and ABI descriptor
comparison, native classic benchmark, sanitizer and source-tarball rebuild.
The new interface remains experimental. The reference PR is not modified.
No claim of final qualification is made until the complete gate is green.
