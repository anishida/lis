# Linear solver USER callback failures

This change concerns the standard-precision linear solvers and the
existing application callback ABI. It is not new MPI, QUAD or SWITCH
support.

## Status and output

A nonzero operator or preconditioner callback result must stop the
affected recurrence before the failed output is used. The callback code
is preserved in the solver's status; the existing public API convention
is unchanged: the caller must check both `lis_solve` and
`lis_solver_get_status`. The solve driver skips final residual certification after a non-successful
iteration. It must not call an application operator again after a failure.
A successful public API return does not imply a successful solve.

- On failure before a usable residual is available, `resid` is NaN.
- Otherwise the last available recurrence residual estimate is retained
  where available. It is diagnostic, not a certification of the returned
  iterate. If a newly updated iterate could not be checked, the
  affected restart paths report NaN.
- `iter` describes solver progress at exit according to the individual
  method's existing iteration convention. Initial-residual failure is
  iteration zero; it is not a uniform count of completed outer steps
  across all algorithms.
- A failed solve does not promise a converged physical solution. Most
  methods retain the last available or partial physical iterate. In
  BiCGSTAB(l), a failed preconditioned correction is not a physical
  iterate, so the saved initial iterate is restored and residual is NaN.
  Do not accept any returned vector as a solution after a nonzero status.

The initial residual helper already records callback errors and the
solve driver clears the status before each execution. For the default relative-to-initial-residual condition,
a finite nonzero initial norm is already normalized to one: the
additional diagnostic initialization reuses that fact rather than
performing another norm reduction. Alternate convergence conditions
may require an extra norm. No extra history entry or print operation
is added by the diagnostic initialization.

## BiCGSTAB(l) cleanup and recycling

BiCGSTAB(l) frees tau on callback failure and no-improvement exits.
Its numerical initial-iterate reconstruction correction is separate.
GCRO-DR and FGCRO-DR discard an affected recycle relation on callback
failure before returning. A subsequent solve can build a fresh relation.
Ordinary exhaustion of the iteration budget preserves the recycle relation. An early no-improvement exit from GCRO-DR takes the conservative failure path and discards it.

NaN diagnostics require IEEE-compatible floating-point semantics.
Compiler options assuming all values are finite (such as -ffinite-math-only,
-Ofast or fast floating-point models) can invalidate NaN tests and are
not supported by this diagnostic contract.

## Precision and MPI boundaries

CG and MINRES already have callback checks. The tests for these methods
also verify the shared certification guard. QUAD/SWITCH implementations
are not extended here. Direct
matrix-free USER/OPERATOR solves are restricted to default precision
by the existing capability check.

The shared internal callback-origin scope uses thread-local storage in serial builds where supported. The eigensolver driver activates it; it does not change the public linear API or add a public status field. Nested application calls are isolated and only the returned callback status crosses the application boundary. A row-parallel matvec reports a worker failure to its caller after the worksharing region; other rows already scheduled in that same operation may still execute. No subsequent solver operation may consume the failed result.

The new callback tests are serial. This patch does not collectively
agree on rank-local callback error codes. For distributed application
callbacks, all participating ranks must follow a compatible collective
schedule and agree on success/failure before returning; a rank-local
error alone can leave peers blocked in a later collective. The same
requirement applies to distributed preconditioner callbacks. Do not
interpret this patch as support for rank-asymmetric failures.

## Regression and CI

The ten native tests and their Automake targets are all in this
change; no dependency on the eigensolver or serial-contract PRs.
Fault injection covers selected initial and iteration applications
and selected restart/reconstruction paths, not an exhaustive proof
of every callback site or persistent recycle transition. Recycle tests
include a previously populated space, invalidation on failure and retry.

The workflow has two targeted PR profiles and four profiles on
matching develop pushes. The explicit paths intentionally exclude
test/Makefile.am and generated Makefile.in, following the CI budget
policy. Unrelated test registrations must not launch this
specialized suite. Test sources and the workflow itself trigger it.
The existing baseline handles build-system registration changes.
OpenMP is enabled but this suite uses one thread; dedicated complex
and multithreaded testing belongs to the separate serial-contract
package.
