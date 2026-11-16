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

These tests **do not** establish a contact-model convergence advantage,
iteration speedup, race freedom under arbitrary concurrent calls, or safe
reuse across changing matrix generations. Benchmark on saved *identical*
CalculiX/contact `A,b,near-nullspace` snapshots with the same outer
FGMRES settings. Compare true physical residual, iterations, factor setup
time, total wall time, and peak resident memory against standalone ILUTP
and standalone MFAMG. The global variant builds an additional full-matrix ILUTP factorization,
which can be expensive. The RAS variant builds **one local ILUTP factor per
domain**, potentially reducing individual setup cost but increasing aggregate
memory as overlapping subdomains are replicated. Its `-ilutp_full_limit`
applies **per domain**, not globally. Increasing the number of domains may
also weaken the preconditioner; fewer iterations or more concurrent threads
do not necessarily mean faster end-to-end solves.

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
This PR is experimental and should not be merged without CI and real
contact-matrix verification. The optional RAS composition is a separate
numerical mode with a wider validation surface. Splitting it into a
follow-up PR is reasonable if maintainers prefer the global composite to
be reviewed independently.
