# Experimental reverse communication interface (proposal)

Applications with an application-owned operator or preconditioner can drive
CG, BiCGSTAB, restarted GMRES or FGMRES without embedding their callbacks in
the existing LIS solver dispatch. The optional bridge runs the same interface
with an assembled LIS matrix and identity or Jacobi preconditioning. No existing
solver is replaced. This is a proposal for maintainer review, not a release or
a commitment to API stability.

Build with `./configure --enable-rci --disable-mpi --disable-quad --disable-complex --disable-longdouble`; add `--enable-omp` and/or `--enable-longlong` if desired.
Default builds do not install RCI headers or export RCI symbols. The experimental
module requires C11 and non-MPI real-double LIS. Unsupported enabled build
combinations fail during configure.

## Usage

Include `lis_rci.h` and create a context with the dimension, iteration limit,
restart length (GMRES/FGMRES), relative and absolute tolerances and method.
Call `lis_rci_start` with b and an optional initial guess. While it returns
MATVEC or PRECON, fill every entry of the request output with A times the input
or the application of the preconditioner, respectively. Call `lis_rci_advance`
with zero on success or an application error code on failure. Do not change the
input or reuse request pointers after advancing. DONE indicates that an
independent matvec verified the residual; negative status indicates failure.
The solution and diagnostic buffers remain owned by the context. Reuse a
finished context with another right-hand side of the same dimension.
See `test/testrci.c` for an executable example.

CG requires a symmetric positive definite operator and fixed SPD preconditioner.
BiCGSTAB and GMRES require fixed linear preconditioning; FGMRES allows changes.
GMRES and FGMRES share a right-preconditioned Arnoldi implementation storing
preconditioned basis vectors. Convergence is not guaranteed for arbitrary
operators and preconditioners. The bridge supports only identity and Jacobi.
No MPI, complex, long double, quad or other platforms are qualified.
A single context/bridge must not be used concurrently. Independent generic
contexts contain no shared mutable state; LIS bridge calls inherit LIS
thread-safety constraints. Allocation occurs at creation; callbacks own
their allocations. The diagnostic residual is absolute, not relative,
and can be stale on failure.

For the bridge, compile against matching installed LIS headers and flags,
including `-D_LONG__LONG` when appropriate; check `LIS_RCI_CHECK_CLIENT_ABI()`
before creating typed LIS objects. This check screens scalar/index widths and
is not a full layout proof. Keep the borrowed assembled matrix alive until
destroying the bridge; never mutate it during a solve.

## Origin and qualification

The sources are an independent AI-assisted implementation of standard iterative
methods developed in this contribution, not imported from the unavailable prior
RC3 or copied from an external RCI package. The proposed terms are those in the
repository's COPYING; the contributor must confirm the license before merge.
No third-party matrix data, opaque source archive or external numerical
dependency is included in this PR. Historical qualification (including
CHOLMOD-assisted structural cases and failed weak-preconditioner trials) is
retained in the contributor's fork, not part of upstream's test obligations.
The single CI workflow runs deterministic tests and compares existing
ABI exports against the current merge base. See the PR for evidence.

## Extended RCI and protocol

The original two-request API remains available. The extended interface
(`lis_rci_start_ex`, `lis_rci_advance_ex`) additionally externalizes every
dot product and Euclidean norm, including RHS norm and final residual norm.
All vectors are full-length real-double vectors. DOT supplies two read-only
inputs; NORM supplies one. Supply a finite long-double scalar on completion
(nonnegative for NORM); vector requests use their existing output arrays.
The caller is responsible for computing correct scalars, just as for MATVEC.
Each request has a nonzero uint64_t token, monotonic within one context across
restarts and solves. A wrong token returns PROTOCOL without consuming the
outstanding request. Tokens do not protect writes to expired pointers.
Do not mix legacy and extended entry points within a solve. The legacy API
executes the scalar requests internally using the existing local backend.
No vectors are allocated for individual requests.

An optional, transport-independent consensus helper packs five uint64_t words:
`[token, ~token, tag, ~tag, error_word]`. The tag is an opaque numeric
identifier, compared bitwise; it is not floating point equality. Componentwise
unsigned maximum combines records. An error word is zero for success;
otherwise its high 32 bits encode the inverse rank and its low 32 bits encode
the signed error. After merging a nonempty set with unique ranks, unpack reports
token/tag disagreement and the error of the lowest failing rank separately.
The caller chooses precedence. This is a pure local primitive, not an MPI
collective or distributed solver. It does not compare arbitrary long-double
scalar values, and it does not establish compatibility with the unavailable
original RC3 protocol/API.

`test/testrciext.c` drives all four solvers with external scalar operations,
tests stale tokens and reuse, and checks 8192 deterministic consensus scenarios
and a token difference above 2^53. The installed-library CI runs both legacy
and extended APIs.

## Full vector-algebra mode

Use `lis_rci_start_full` and `lis_rci_advance_full` to additionally receive
AXPBY (out = alpha*in + beta*in2, with an optional second operand) and
DIVIDE (out = in/alpha). These operations cover iterative vector updates,
residual formation, direction updates, Arnoldi orthogonalization and basis
normalization. Evaluate the formula in long double and convert to double.
The output is a separate context-owned staging vector; a completed finite
result is copied into the actual destination only after token/status validation.
One extra length-n vector is reserved at context creation, with no per-request
allocation. The small dense Hessenberg and Givens operations, input snapshots
and scratch commits remain internal. BiCGSTAB full mode uses two AXPBYs for
its direction recurrence, so intermediate rounding may differ from legacy.
The earlier two modes continue to perform these operations internally.
This is not an opaque-vector/device-resident or distributed backend interface.
`test/testrcifull.c` exercises all methods, context reuse, mode isolation,
stale tokens and both callback failures and incomplete/nonfinite outputs at
every request position of a deterministic solve.

## Request lifecycle and checked consensus

`lis_rci_get_pending` inspects the active request in any of the three modes.
`lis_rci_cancel` cancels the current token once the application has stopped
all workers accessing its request buffers; cancelled solves can reuse their
context. Cancellation does not roll back partial iterates. A stale token is
retryable and cannot consume or cancel the live request.
`lis_rci_get_stats` exposes per-solve issued/accepted-completion counts for
each operation, last token, cancellation state and consensus error rank.
Legacy internal DOT/NORM operations are included. Counts reset at start,
while tokens remain monotonic.

`lis_rci_advance_consensus` is an optional completion path for Extended/Full
mode. Pack the request token and operation as the tag, merge the records,
then submit the result and the computed scalar (if applicable). Agreement
on the pending token and operation is mandatory. A disagreement terminates
with PROTOCOL; protocol error takes precedence over a selected callback
error, whose code and rank remain available for diagnostics. A local wrong
token or malformed record does not consume the pending request. The caller
must provide truthful records, unique ranks, correct computation and any
transport; this checks neither arbitrary scalar agreement nor cross-context
identity. All vectors are still local full-length arrays. No MPI.

## Cross-mode qualification

The large native LIS test drives every method through Legacy, Extended and Full
against the same 8192-row assembled matrix and checks an independent residual,
forward error and agreement with the legacy solution. CG uses an SPD operator;
the other methods use a nonsymmetric operator. FGMRES varies its preconditioner
by preconditioner call count. This adds 24 cross-mode solves per test invocation.
The lifecycle test covers zero RHS, exact initial guess, iteration limit and
breakdown for all method/mode combinations. Full-mode fault injection also
cancels at every request position and verifies context reuse after cancellation.

## Local execution overhead

For small vectors (n < 4096) the built-in kernels in OpenMP builds use a direct
serial loop instead of entering an inactive OpenMP region. Larger vectors
retain the existing parallel reductions. A DOT request with identical operands
validates the input once rather than twice. NaN output poisoning, input/output
finite checks, token validation and public operation counters are preserved.
These changes do not alter the request protocol or expose opaque/device
vectors. Timing depends on machine, compiler, method and callback cost.

## Retained earlier qualification sources

The earlier reconstructed integrated candidate included three additional
standalone tests. Their source bodies are retained as `test/testrcisweep.c`
(previously `test_rci.c`), `test/testrcifaults.c` (`test_faults.c`) and
`test/testrcipoisson.c` (`test_poisson.c`). These restore the earlier size,
restart, preconditioner and initial-guess sweep, fault injection at every
legacy request, and scaled two-dimensional Poisson systems. They are built and
executed alongside the later three-mode tests. No external matrix files or
MPI dependencies are required. The archived structural matrix drivers and
industrial qualification remain on the fork's archive and qualification
branches. The archived implementation was itself a reconstruction, not the
unavailable original RC3.

## Numerical contract and recovered audit

RCI requires nonfinite-aware floating-point compilation. GCC/Clang-style
fast-math and finite-math-only options are rejected by configure and by the
RCI translation units; these optimizations invalidate NaN output poisoning.
Other compilers must preserve IEEE nonfinite checks. The legacy double
residual diagnostic explicitly reports +Inf if the verified absolute norm
exceeds DBL_MAX. The additive `lis_rci_get_residual` accessor returns the
stored long-double norm (NaN before verification), without changing the
convergence decision or any existing public structure.

The recovered independent consensus oracle now tests the actual exported
functions on 109611 scenarios (forward, reverse and tree reductions).
The native bridge test checks invalid dimensions, zero and negative Jacobi
diagonals, nonfinite input, unchanged output on failure and recovery. This
adapts earlier mock tests; mock-injected LIS matvec errors and allocation
failures remain separate qualification work. See `generic-rci-abi.md` for
the proposed compatibility policy. No MPI or opaque device vectors are added.

## CI and review gates

The experimental RCI workflow covers the four index/thread combinations,
with ASan/UBSan on each combination as well as release builds. It tests
installed-library C/C++ consumers, unsupported configure options, and all
legacy LIS symbols and headers. Exported RCI symbols must match the public
function declarations. On develop after merge the compatibility job compares
the default and opt-in libraries from the same revision; PRs also compare
the default upstream merge base. No external matrix download is required.
The upstream maintainer controls which checks are required for merge; all
RCI jobs and existing upstream checks should pass before review readiness.
The industrial benchmarks are informational and remain in the fork.

Wide intermediates are range-checked before narrowing to double; a finite DOT
outside the representable range of the dense GMRES Hessenberg storage is
reported as NONFINITE rather than silently overflowing a cast. This does not
extend the vector or dense storage representation beyond real double.


## Scope relative to the earlier Generic RCI design

The recovered historical design reports eight algorithms: CG, BiCGSTAB,
GMRES, FGMRES, GCRO-DR, FGCRO-DR, FEAST and Krylov eigensolver. This PR
implements only the first four, as an independent experimental subset.
Historical reports and artifact names do not establish that the original
source files have been recovered or that their results apply to this PR.

The broader non-MPI design is not superseded by this subset:

| Historical design component | Current implementation |
| --- | --- |
| Four basic linear methods | Independent real-double RCI |
| Versioned size-checked request descriptor, opaque vector handles | Not implemented |
| A/B operators, block RHS, complex shifts, split-complex | Not implemented |
| Nested contexts, operator generation checks, reusable workspace | Partial (local context reuse only) |
| Native LIS_SOLVER and lis_psolve adapter | Partial (matrix + identity/Jacobi only) |
| GCRO-DR U/ZU/C and FGCRO-DR W/Z/C, harmonic recycling | Not implemented |
| Persistent sequences with transactional recycle commit | Not implemented |
| FEAST contour projection and outer B-metric Rayleigh-Ritz iteration | Not implemented |
| Stackless Krylov with A/B/shifted solve continuations | Not implemented |
| Shared native classic/RCI solver kernels and policies | Not implemented |
| Complex, extended and binary128 precision | Not implemented |
| MPI and collective transport | Explicitly outside current scope |

The historical migration approach is to reuse the actual native LIS kernels,
not to add competing implementations of GCRO-DR, FGCRO-DR, FEAST or Krylov.
Preserve existing native workspace, recycling keys, convergence criteria,
residual history, error policies, ABI and operation counts. The planned
sequence is (1) recover exact historical sources and API layout where possible,
(2) versioned opaque backend and operator lifetime/dispatch contract,
(3) native shared linear kernels and recycling, (4) shifted/block projection,
FEAST and stackless Krylov, (5) full build and differential qualification.
Avoid advertising any unimplemented algorithm or numerical configuration.

Tracked build inputs for this subset were regenerated with Autoconf 2.72 and
Automake 1.18.1, retaining the existing Libtool macro/helper set. The clean
checkout serial/OpenMP build and upstream tracked Automake comparison passed
in the fork's targeted generator qualification (run 45966766386).
