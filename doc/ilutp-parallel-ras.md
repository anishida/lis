# Parallel local ILUTP using independent OpenMP RAS domains (experimental)

## Scope and honest algorithm distinction

This feature runs **independent native ILUTP factorizations and solves in
parallel across RAS subdomains**. It does NOT parallelize the dependent pivot
sequence of one global ILUTP/LU factor. The original whole-matrix `-ilutp 1`
continues to reject `omp_get_max_threads() > 1`, as before.

The preconditioner is restricted additive Schwarz (RAS), using contiguous
nonoverlapping owner cores with optional sparsity-graph overlap. Local factors
may use dynamic column pivoting (`-ilutp 1`), optional RCM, row scaling and
bounded complete fill (`-ilutp_full 1`). Each worker temporarily sets its
**nested OpenMP thread budget to 1** so that the serial ILUTP pivoting and
triangular solve remain correct; it restores the caller's thread count when
its local work finishes. Up to min(`OMP_NUM_THREADS`, `-ilut_ras_domains`)
local factorizations and applies are handled concurrently.

This is opt-in and alters the preconditioner mathematically; it does not
promise the convergence of global complete-fill LU. Cross-domain couplings
are approximated by RAS, not included in a global factor. The Krylov solver
still operates on the original matrix A and original right-hand side b.

## Windows / CalculiX contact trial

From a Windows batch file, before launching CalculiX:

```bat
set "OMP_NUM_THREADS=8"
set "CCX_LIS_NL_CONTACT_OPTIONS=-i fgmres -restart 120 -p ilut -ilutp 1 -ilutp_threshold 0.5 -ilut_ras_domains 8 -ilut_ras_overlap 1 -ilut_schur 0 -ilut_drop 0.005 -ilut_rate 12 -ilu_pivot_tol 1e-9 -scale none -ilut_rcm 1 -ilut_row_scale 1 -omp_num_threads 8 -maxiter 1000"
```

An alternative bounded complete-fill local LU trial is:

```bat
set "OMP_NUM_THREADS=8"
set "CCX_LIS_NL_CONTACT_OPTIONS=-i fgmres -restart 120 -p ilut -ilutp 1 -ilutp_full 1 -ilutp_threshold 0.1 -ilutp_full_limit 60000000 -ilut_ras_domains 8 -ilut_ras_overlap 1 -ilut_schur 0 -scale none -ilut_rcm 1 -ilut_row_scale 1 -omp_num_threads 8 -maxiter 1000"
```

The full-fill limit applies **to each local factorization**, not as a total
RAM limit; each domain may allocate a large amount of memory. Prefer the
limited-fill trial first. The default `-ilut_schur 0` is deliberate: the
experimental dense approximate interface Schur correction is limited to
256 interface DOFs and is not suitable by default for large contact systems.

## Correctness and thread model

* Input: serial, real-double explicit CSR only. Not distributed MPI RAS.
* Local ILUTP: one OpenMP thread **per factor**, concurrently across domains.
  Calling the unrestricted whole-matrix ILUTP at 8 threads is still rejected.
* Construct and apply must not run simultaneously on the same preconditioner.
  The same RAS preconditioner instance cannot be applied concurrently.
* Each local RHS, local solve vector, pivot map and factor storage is
  exclusively owned by its subdomain; RAS writes only to disjoint owner cores.
* A failed domain setup finishes all workers before generic cleanup; apply
  collects error codes after the barrier instead of returning from a worker.
* ILUTP does not support FP16 local factor storage. Adjoint RAS is still
  explicitly unsupported. Overlap values 0..2, domains 2..16.
* The existing core partitioning is contiguous, not a graph partitioner.
  Performance and convergence may be sensitive to the matrix ordering.
* Local factorization scheduling is dynamic (one domain per work item)
  to balance unequal fill-in; apply uses static scheduling.
* The global allocation registry is protected against concurrent mutations.
  Large raw allocations and free() occur outside the registry lock, but
  tracked reallocations and metadata updates remain serialized.
* The shared error-printing/debug subsystem and other legacy code have not
  been proven fully race-free under arbitrary concurrent failing operations.
  Full end-to-end race analysis (e.g. Helgrind/TSan with runtime-compatible
  OpenMP instrumentation) and real contact-system benchmarks remain review
  tasks. The new native OOM+OpenMP CI job injects allocation failures while
  eight local factors are being constructed and verifies that all partial
  subdomain allocations are released after the worker join. Passing OpenMP regression tests is not a proof of race freedom.
* RCM and row scaling operate on **each extracted local matrix**. They do not
  implement whole-system two-sided Equil in this branch.

## Regression / validation

`test/testhalfilut.c` extends the manufactured local pivot test to
32x32 matrices containing sixteen off-diagonal 2x2 pivot blocks with zero
diagonal and known solution. Eight disjoint four-row RAS domains are tested
at 1, 2, 4, and 8 OpenMP threads. The test verifies the original-system
relative residual and exact-solution error below 1e-10, and ensures the
outer OpenMP thread budget is restored. It retains the negative test that
full fill without ILUTP is rejected.

A second 32x32 regression deliberately couples different owner cores using
nonsymmetric off-diagonal links, with graph overlap 1 and 2, local RCM and
row scaling, and both dropped ILUTP and bounded complete-fill ILUTP.
For each mode it checks the original-system residual and requires an exactly
identical preconditioner application result for 1/2/4/8 threads and repeated
calls, including reads from overlapping imported RHS rows.

Run, for example:

```sh
autoreconf -fiv
./configure --enable-test --enable-omp --disable-fortran --disable-f90 --disable-mpi --disable-saamg
make -j2
OMP_DYNAMIC=FALSE OMP_NUM_THREADS=8 ./test/testhalfilut --parallel-ras-ilutp
make -j2 check
```

An FGMRES solve of a **real contact system** must be benchmarked separately.
Compare the same original A/b on the same physical residual, factor setup
wall time, solve wall time, memory high-water mark and convergence. This
feature has no demonstrated end-to-end nonlinear CalculiX contact benefit
until such cases have been tested. In particular it does **not** fix
mechanical Newton divergence, rank deficiency or an invalid tangent matrix.

## Relationship to upstream

This change is stacked on
[3rav/lis#63](https://github.com/3rav/lis/pull/63) (experimental RAS/Schur,
based on native ILUTP [anishida/lis#229](https://github.com/anishida/lis/pull/229)).
It must be reviewed/merged **after** the dependencies. It leaves the native
ILUTP global one-thread guard intact and adds no new public solver flags.
