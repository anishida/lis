# Reconstructed historical RCI foundation (qualification only)

Based on the supplied historical transcript (v0.42, v0.46 and v0.51).
Not recovered original source; no claim of its ABI or numerical equivalence.

- Separate size/version-checked opaque handle request with A/B, PRECON,
  shifted solve, zero and complex AXPBY; generation and token guards.
- Block dispatch, with deterministic per-column fallback and error handling.
- Non-MPI contour projection, with complex shifts and weights and caller-owned
  vector storage, using the same machine for standard and generalized cases.
  The caller implements the backend, including shifted solves. No result is
  published on error; scratch can be modified and is not rolled back.
- Synchronous nested real-double solver adapter uses existing four core methods;
  it does not implement complex shifted solves or a native LIS preconditioner.
- Tests: independent diagonal complex reference, block/column dispatch, reuse,
  error/generation/cancellation injection at each phase, nested four methods,
  version, size and aliasing checks.

No public headers installed, no configure/Autotools changes. The source change
adds an internal qualification dimension accessor, to avoid reading beyond an
incompatible child context. Keep this branch separate from the upstream PR until
the new API design and qualification have been reviewed.

Not yet shared classic solver kernels, GCRO-DR/FGCRO-DR,
stackless Krylov or extended precision LIS builds. Complex
coefficients/opaque vectors are not equivalent to full complex LIS support.

The split-complex adapter now implements (zB-A) with two real planes and a
persistent doubled-real linear context. It supports BiCGSTAB, GMRES, FGMRES
with identity preconditioning, not CG (complex shifts need not be SPD).
Standard/generalized non-diagonal 3x3 tests compare the actual nested solves
and two-RHS contour projection against independent complex Gaussian
elimination. Input and output planes must be disjoint; the caller owns them.
Callback failure preserves the current column's output; earlier columns
may already have changed. No full complex LIS ABI support is claimed.

A separate native real-double adapter borrows assembled A/B matrices and a
prepared LIS_SOLVER/preconditioner. It uses lis_matvec and lis_psolve,
checks both structure/value generations and invalidates on replacement of the
prepared preconditioner. The adapter does not prepare arbitrary native
preconditioners itself; the caller must meet their workspace and lifecycle
requirements. Qualification tests use native Jacobi and SSOR, with all four
linear methods in Serial/OpenMP and INT32/INT64. Native bridge integration
does not imply a shared classic numerical recurrence or native FEAST.

The experimental FEAST outer continuation uses the split-complex contour
projector, B-metric twice-MGS rank reduction, native projected eigensolve
callback, Ritz reconstruction and explicit independent A/B residual
verification. It requires real symmetric A and SPD B, conjugate quadrature,
and a suitable starting subspace. It returns only converged interval values
after stable count/trace; it does not certify completeness or an empty
interval, and is not integrated into the classical FEAST policy. Every outer
request has a context-local monotonic token independent of nested projector
tokens. The native test uses lis_eigen_projected_solve and actual assembled
LIS A/B matrices (known spectra 1,2,5,8); tests repeat, cancellation, callback
error, generation invalidation and independently check physical residuals.

The recycle-space continuation stores fixed GCRO U/ZU/C and flexible
FGCRO W/Z/C separately (Z authoritative), refreshes images with A and
fixed preconditioning, performs twice-MGS with identical transforms,
drops rank-deficient candidates, and commits via a two-bank transaction.
Cancellation, incomplete output, generation mismatch or callback failure
keeps the previous committed basis. Residual and direction projection
preserve A*Z=C and orthogonality to C. Tests cover two RHS, changed operator,
fixed/flexible policy, and all request-phase faults. This is a building
block, NOT complete GCRO-DR/FGCRO-DR: native harmonic Ritz selection,
automatic integration into the augmented Arnoldi/restart recurrence and
differential native qualification remain open.

The augmented linear sequence now drives the existing GMRES/FGMRES RCI
recurrence cycle-by-cycle. It projects residuals and preconditioned
directions against a compatible committed recycle basis, retains the
U/ZU/C or W/Z/C augmented Arnoldi relation, and uses the actual native
GCRO-DR small QR/extraction/inverse-iteration routines for harmonic
selection. It commits transformed physical/image bases without new A/M
calls. Selection failure is nonfatal; stale operator/preconditioner
generations are conservatively not reused. The native dense selector can
allocate and is qualification-only. This is an experimental augmented
sequence, not yet differential equivalence to classical GCRO-DR/FGCRO-DR
or shared native classic recurrence. Initial coarse projection incurs one
extra true MATVEC before child verification; no duplicated Arnoldi MATVEC.

The native symmetric Krylov continuation is generated from the exact audited
upstream source blob f275ffbb5cf3ab0a6757d5e59eaefd9dd9f9b74c. Nine
functions and 28 suspend sites are transformed into persistent C11 frames;
all synchronous numerical helpers remain the original functions. The
generator deliberately rejects any other source blob and CI checks the
committed generated output. No native callback invocation, worker thread,
stack switching or replay in the generated path. The native shift convention
is A-sigma B. Qualification compares against lis_krylov_symmetric_ex64,
including callback input/output byte traces, results, native/extended stats,
errors, cancel and budget/workspace reuse. This is not yet integrated into
the installed library or the native classical entry point.

The qualification-only synchronous driver now invokes the same stackless
session through the existing native callback ABI. Three-way comparison
covers the original classical implementation, external RCI, and the
blocking driver. An interleaved-session probe cancels one context while
the other continues; a 32768-vector test crosses the native OpenMP
threshold with four threads. No existing classical entry point is
redirected yet.

A disposable native integration job uses qualify_krylov_native.py to patch
the exact audited source in a fresh checkout, redirect the existing
lis_krylov_symmetric_ex64 entry point to the blocking continuation driver,
and retain the original core as a differential reference. It links an
external RCI consumer against the actual rebuilt liblis.so and runs all
17 existing upstream Krylov regression executables. The patch is never
committed to native source; these are qualification builds, not changes
to the default upstream PR or installed headers.

The source-specific native linear generator targets actual CG, BiCGSTAB,
GMRES and FGMRES plus the shared initial-residual helper. It exposes native
LIS_VECTOR A/PRECON requests at the real call sites, leaving reductions,
convergence, residual histories and work vectors under native LIS policy.
A disposable qualification patch redirects the four existing classical
symbols through a blocking continuation driver while retaining original
reference functions. The normal lis_solve setup/teardown is unchanged.
Strict external mode validates output; classic mode preserves native
behavior. Callback wait time is included in native preconditioner timing.
Cancellation may leave the native iterate partially updated. This is not
yet a replacement for the existing generic real-double RCI API, nor an
installed native interface. Precision and native-policy coverage require
qualification rather than an assumed equivalence.

Native dispatch must use the called function's method, not the solver option:
FGCRO-DR intentionally calls lis_fgmres directly for its zero-recycle
fallback. A full-regression failure exposed the distinction. The driver
now passes an explicit method without mutating solver options, including
the test-only external dispatcher. Full regression is rerun after this fix.

Native GCRO-DR qualification now generates continuations from the actual
GCRO-DR main recurrence, initial residual, recycle refresh, policy assessment
and harmonic update. The original helper functions remain in the native
translation unit; no callback is hidden inside those four transformed
GCRO-DR routines. Cancellation/incomplete output unwinds the original
error path, which may invalidate the native recycle basis (distinct from
the transactional generic experiment). This is a real-double-only
qualification patch, not an installed API. FGCRO-DR native continuation
and FEAST native continuation remain separate open tasks.

The native GCRO differential harness also injects cancellation, incomplete
output and a native callback error at each observed request of a cold and
warm solve, including refresh and harmonic-update calls. It checks that
native failure semantics invalidate any possibly partial recycle relation.

Native positive-k FGCRO-DR continuation candidate: the pinned generator
covers nine native FGCRO-DR functions plus the shared initial-residual
helper (19 expected continuation sites, including child calls). It preserves
W/Z/C ownership and the original harmonic build/update and replacement
routines. The native zero-recycle FGMRES fallback is left unchanged; the
external positive-k session rejects zero recycle dimension rather than
silently executing hidden callbacks. This is not a public or installed API.
The qualification workflow compares original, blocking, and external
paths, including private W/Z/C snapshots and independent relation checks,
and exercises cancellation and incomplete output on cold/warm solves.
CI success and native callback-error injection are not yet established.

Native FEAST qualification: the isolated split-projector run 45992283834 passed all three CI jobs (including ASan/UBSan). The full native candidate generates nine continuation frames from three pinned source blobs and exposes scalar A/B, block A/B, and split-complex shifted solves. The original outer driver, filtering, rank reduction, Rayleigh–Ritz and certification code are retained as a differential reference. The two private synchronous block helpers are exposed only in the disposable qualification checkout, not the installed public API. Strict mode poisons output vectors, validates completion and rejects stale tokens. Qualification compares standard/generalized real-double eigenpairs and diagnostics bitwise, independently computes physical residuals, and exercises failure/cancellation/incomplete output at selected boundaries including every request type. Full integration and qualification remain pending; MPI is outside scope.

The full FEAST strict session additionally snapshots the borrowed A/B identities and structure/value generations. A stale generation is rejected before accepting a reply; the original error path unwinds scratch and publishes no certified eigenpairs. Tests inject generation counter changes at the first and final requests. This does not make destruction or concurrent mutation of borrowed native objects safe. The native dispatcher reports numerical failure through esolver->retcode, even when lis_gesolve returns success for dispatch; qualification asserts that contract.

All eight separate native core candidates have passed qualification in real double (full FEAST run 45992950618). The combined qualification applies linear, GCRO, FGCRO, Krylov and FEAST transformations in separate Python processes to one clean checkout, builds a single liblis.so, links every native differential harness to that same library, and runs the full upstream test suite in serial and OpenMP. Sanitized combined differential tests are separate. A baseline-versus-combined dynamic symbol and SONAME check detects missing symbols but does not constitute a complete ABI/layout qualification; temporary helper exports must be removed or hidden before release. This remains an experimental branch, not an installed public RCI interface. Combined results are pending.

Installed API candidate: run prepare_native_installed_api.py on a pristine pinned checkout, then configure --enable-rci and build/install. The public header and ten facade symbols expose the full linear/FEAST lifecycle and native Krylov under one context. Existing classic and array RCI ABI stays separate. The disposable ELF export allowlist contains baseline symbols plus precisely the ten proposed public symbols. Disabled build retains original dispatch. Installed-client smoke uses only installed headers and shared library. This is not yet permanent source/distribution integration: preparation and ELF symbol policy are qualification-only, pending CI and maintainer review. Real-double supported subset only; scalar-extended experiments remain separate.
