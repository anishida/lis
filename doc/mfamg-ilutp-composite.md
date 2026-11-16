# Experimental MFAMG + native ILUTP residual composition

This opt-in *built-in* preconditioner combines an MFAMG V-cycle with a
native pivoted ILUTP correction. It is selected with
`-p mfamg_ilutp` in the current development branch. By default ILUTP uses
one global serial factor; an explicit `-ilut_ras_overlap 0..2` option
selects **domain-parallel, local native ILUTP through restricted additive
Schwarz (RAS)**.

## Algorithm

For input residual `r`, the composite applies

```text
z1 = MFAMG(r)
r1 = r - A*z1             (same assembled A used by the Krylov solver)
z2 = ILUTP(r1)
z  = z1 + z2
```

Both numerical preconditioners are built during setup and reused for each
application. Failures are returned to the caller; there is **no** silent
ILUTP-only or MFAMG-only fallback. This is multiplicative residual
correction, *not* ILUTP pivoting inside MFAMG blocks and *not* recursive
`Hybrid`. The outer FGMRES still verifies convergence on the original
system.

## Supported matrices and execution modes

- The operator is a real-double, non-split CSR matrix **within one LIS/MPI
  process** (`nprocs=1`); no distributed matrix or inter-rank exchanges are
  supported. This serial *matrix ownership* restriction does **not** prohibit
  multiple OpenMP threads operating on independent local RAS domains.
- `-scale none`; at least one caller-supplied MFAMG near-nullspace vector.
- **Global mode (default):** one native ILUTP factor with a serial pivot
  sequence. Set `OMP_NUM_THREADS=1` **before** creating the solver.
- **Local RAS mode (opt-in):** set `-ilut_ras_overlap 0..2` and
  `-ilut_ras_domains 2..16`. OpenMP workers create/apply independent native
  ILUTP factors on overlapping subdomains. The local pivot sequences remain
  serial per domain; parallelism is **across domains**. Domain count must not
  exceed the matrix dimension. The RAS correction is **not equivalent** to
  the global ILUTP factorization. Without OpenMP, the same RAS path is
  serial and provides no multithreaded speedup.
- No distributed MPI ranks, complex arithmetic, quad arithmetic, or
  Additive Schwarz wrapper.
- No transpose/adjoint `psolveh` implementation.
- No `-hybrid_p mfamg_ilutp` nesting. Select `-p mfamg_ilutp` directly.
- No optional RAS Schur interface correction (`-ilut_schur 1` is rejected).
  The RAS path uses restricted owner-core updates; no MPI communication.
  Merely setting `-ilut_ras_domains` without enabling overlap fails.
  The existing global-ILUTP path remains unchanged.
- The composite is **not reentrant**: one preconditioner instance contains shared
  apply-time work vectors and must not be called concurrently from different
  threads. Independent instances can be created separately.
- Numeric PSD/sequence preconditioner refresh is not implemented in this
  first version. Each apply checks the original matrix pointer and the
  structure/value generation recorded at setup, rejecting stale factors
  **only when the generation counters are updated through LIS APIs**.
  Applications that modify `A->ptr`, `A->index` or `A->value` directly,
  without calling `lis_matrix_mark_structure_changed()` or
  `lis_matrix_mark_values_changed()` as appropriate, bypass this safeguard.
  Never mutate the CSR storage during a preconditioner apply.
  Rebuild with ordinary `lis_solve` after each changed Newton/contact
  matrix; do not cache this preconditioner across matrix updates.

## Examples (LIS option strings)

```text
-i fgmres -p mfamg_ilutp -scale none
-mfamg_block_size 2 -mfamg_vcycle_smooth_steps 1
-ilutp_threshold 0.5 -ilut_equil 1 -ilut_rcm 1
-ilut_drop 0.005 -ilut_rate 8 -restart 80
-maxiter 1000 -tol 1e-8
```

For domain-parallel correction (requires an OpenMP-enabled build):
```sh
export OMP_NUM_THREADS=4
```
```text
-i fgmres -p mfamg_ilutp -scale none
-mfamg_block_size 1 -mfamg_vcycle_smooth_steps 1
-ilutp_threshold 0.5 -ilut_equil 1 -ilut_rcm 1
-ilut_ras_overlap 1 -ilut_ras_domains 4
-ilut_drop 0.005 -ilut_rate 8 -restart 120
-maxiter 1000 -tol 1e-8
```

The composite selects native ILUTP automatically, even without an explicit
`-ilutp 1` option. A strictly positive `-ilutp_threshold` (at most 1) is
mandatory; threshold zero would select legacy ILUT in limited-fill mode and
is therefore rejected. Explicit Schur options also fail setup. The ILUTP
threshold, optional complete-fill cap, RCM, row scale, Equil, and drop/rate
are inherited from the outer configuration. **RAS full-fill limits apply per
local domain**. `-ilut_equil 1` is forwarded into each local native ILUTP
factor. More RAS domains may improve parallel scaling but can reduce
preconditioner quality and increase memory usage.
MFAMG inherits its block size, fine-aggregate factor, smoothing passes,
deep-level sibling merge, near-nullspace basis, and optional FP16 storage.
ILUTP factors themselves remain in full precision. Existing `-p mfamg`
and `-p ilut` behavior is unchanged.

The host application must call `lis_solver_set_near_nullspace()` with
compatible vectors. A constant vector is useful for diffusion-like models;
contact/mechanical systems often require physically meaningful rigid-body
modes or other application-selected near-nullspace vectors. This code does
not invent those vectors.

## Validation

Run `autoreconf -fiv` and `./configure --enable-test`, then
`make -j2 && make -j2 check`. The new `testmfamgilutp` regression exercises
`lis_solve` for MFAMG and the composite on a deterministic real nonsymmetric
CSR system, checks exact-solution error, and repeatedly creates/applies/destroys
its children, including application with the input and output vectors aliased.
Negative tests check missing near-nullspace data, Additive Schwarz, non-none
scaling, Hybrid nesting, zero ILUTP threshold, invalid RAS settings, requested
Schur mode, multi-rank MPI, multiple OpenMP threads with global ILUTP, and
invalidation on an in-place matrix generation change. In OpenMP builds
`./test/testmfamgilutp --parallel-ras` additionally tests actual two-/four-
worker teams, overlaps 0/1/2, bounded complete-fill local ILUTP, repeated
FGMRES solves, original-system residuals and deterministic solutions.
For a fixed four-domain configuration with overlap 1, the test also compares
the entire final solution **bitwise across two and four OpenMP workers**
(`memcmp`), not merely across successive solves at the same thread count.
An additional RAS-specific ASan/UBSan/LSan job exercises one, two and four
workers; the LLVM Archer/ThreadSanitizer job checks the composite RAS path.
A separate CI job covers complex, quad-build and 64-bit LIS_INT combinations.

The source Automake file is `src/precon/Makefile.am` / `test/Makefile.am`;
tracked `Makefile.in` files are kept in sync to support builds without
`autoreconf`. CI tests both Autotools-regenerated builds and builds using
tracked `Makefile.in`. Regeneration with the project's maintained Automake
version and a minimized diff must be checked before release.

The synthetic tests **do not** establish race freedom for arbitrary
concurrent applications of one instance or safe reuse when matrix storage
is changed without notifying LIS. On production matrices, compare *the same*
A, b, near-nullspace, FGMRES restart, iteration limit and tolerance, and
report **true** residual, factor setup, total wall time and peak RSS.
The global variant builds an additional full-matrix ILUTP factorization,
which can be expensive. The RAS variant builds **one local ILUTP factor per
domain**, potentially reducing individual setup cost but increasing aggregate
memory as overlapping subdomains are replicated. Its `-ilutp_full_limit`
applies **per domain**, not globally. Increasing the number of domains may
also weaken the preconditioner; fewer iterations or more concurrent threads
do not necessarily mean faster end-to-end solves.

## External numerical validation for review (2026-11-16)

The following are **anonymized aggregates from frozen linear Newton/contact
systems** evaluated with native LIS. The matrix and RHS files, meshes,
near-nullspace vectors and displacement solutions remain private and are
**not** part of this repository. These are *single-run observations*, not
Cal-culiX end-to-end nonlinear speedups or controlled performance claims.
The full-system MFAMG experiments used the separately developed,
**experimental** `LIS_MFAMG_SPARSE_TERMINAL=1` facility, because the
unmodified MFAMG hierarchy otherwise stalled on one tested full system.
This is a numerical limitation, not a configuration a user should have to
assume for every input.

For an anonymous difficult 60,265-unknown real-double CSR system, with
six supplied near-nullspace candidates, FGMRES `restart=120`,
`maxiter=500`, `tol=1e-6`, identical native ILUTP parameters and
true residual verified against the original A and b:

| Preconditioner | Threads | Iterations | True relative residual | Setup + solve, seconds | Outcome |
| --- | ---: | ---: | ---: | ---: | --- |
| ILUTP alone | 1 | 154 | 9.93e-7 | 24.30 | converged |
| MFAMG + global ILUTP | 1 | 114 | 6.05e-7 | 25.98 | converged |
| MFAMG alone | 4 | limit 500 | 4.25e-6 | 18.96 | not converged |
| MFAMG + 4-domain RAS ILUTP | 4 | 227 | 8.65e-7 | 44.32 | converged |

The one- and four-thread runtimes are **not a controlled speedup
comparison**. For the two serial ILUTP configurations the composition
reduced Krylov iterations, **but did not reduce total time**. The tested
simple contiguous-row RAS partition was slower and costlier in memory.
With another difficult system, the same `restart=120` configuration
converged in 104 iterations / 14.03 s using ILUTP alone versus 91
iterations / 20.83 s using the global composite.

On a separate difficult system, **standalone native ILUTP, the global
composite and the local RAS composite all returned `LIS_BREAKDOWN`
during incomplete factorization**. A zero pivot emerged after sparse
elimination and pivoting, not from a zero original diagonal. Explicit
pivot regularization permitted setup but caused FGMRES breakdown; it
is not a numerical fix. No silent pivot regularization or fallback has
been introduced.

These results support **reviewing a deliberately limited experimental
preconditioner**, not describing it as production-ready or universally
robust. The separate, model-informed C3D8I static-condensation prototype
is **not part of this PR** and does not constrain the general LIS
preconditioner to a particular finite-element type.

## Ready-for-review boundaries

- Scope is opt-in, real-double, single-LIS-rank, unsplit CSR; multiple
  OpenMP threads are allowed only with explicitly selected RAS correction.
- Caller supplies physically appropriate near-nullspace vectors.
  Numerical convergence is **not guaranteed**; native ILUTP breakdown
  is propagated without a hidden fallback.
- A single composite instance is not reentrant. Rebuild for changed
  matrix generations. Direct edits to CSR arrays bypass generation checks.
- USERDEF numeric ID changes from 13 to 14. The `LIS_PRECON_STRUCT`
  ABI layout is unchanged, but external consumers must recompile.
- The separate C3D8I reduction remains outside this PR.
- The CI matrix covers Linux/Windows, serial/OpenMP, MPI rejection,
  complex/quad and 64-bit integer builds, ASan/UBSan/LSan, TSan/Archer
  and injected-OOM paths. A passing CI is **not** evidence that this
  preconditioner solves every nonlinear contact system successfully.
- Maintainer review/approval and any requested changes remain necessary
  before merging into `develop`.

## ABI / rollout

This development change appends one built-in preconditioner ID (13),
moving `LIS_PRECON_TYPE_USERDEF` from 13 to 14. The internal composite
state reuses the existing `LIS_PRECON_STRUCT::user_data` pointer only for
the new built-in type; the preconditioner struct size and layout are
**unchanged**. Registered user-defined preconditioners keep their prior
`user_data` and `user_destroy` semantics. To prevent applications from
overwriting or obtaining the composite's private owned state, public
`lis_precon_set_user_data()` and `lis_precon_get_user_data()` explicitly
return `LIS_ERR_NOT_IMPLEMENTED` for `MFAMG_ILUTP`; these APIs remain
unchanged for registered user-defined preconditioners. The public `LIS_PRECONNAME_MAX` is deliberately kept
at 10 to avoid changing the `LIS_PRECON_REGISTER_STRUCT` name-buffer size;
the human-readable built-in label is shortened to `MFAMG+P`. Search external
code for hardcoded preconditioner ID 13 before upgrading. Recompile code
against a consistent LIS header and binary; code with hardcoded
user-defined preconditioner IDs must be updated.
This PR is experimental and may be made **Ready for Review** when the
relevant CI checks pass and numerical limitations are documented; it
must not be treated as production-ready or merged without maintainer
review. The optional RAS composition has a wider validation surface.
Splitting it into a follow-up PR is reasonable if maintainers prefer
the global composite to be reviewed independently.
