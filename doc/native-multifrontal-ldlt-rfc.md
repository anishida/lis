# RFC: opt-in native multifrontal LDLT preconditioner

Status: draft proposal with the standalone RC1 numerical prototype imported under `experimental/native-ldlt/`. The prototype is not yet registered or built as a LIS preconditioner.

## Motivation

Symmetric indefinite systems occur in structural mechanics and shift-invert eigensolvers. A reusable pivoted LDLT factorization can provide a strong preconditioner to existing GMRES/FGMRES and generalized eigensolvers. This follows the opt-in native complete-fill LU precedent (#229), without changing the iterative focus of LIS.

## Proposed scope

- Serial, real double precision, assembled symmetric CSR input initially.
- Permuted factorization P A P^T = L D L^T, with 1x1/2x2 pivot blocks, threshold pivoting and delayed pivots.
- Sparse symbolic analysis, elimination tree, strict supernodes, and multifrontal numerical elimination.
- Fill-reducing ordering, bounded allocation and factor reuse; preserve caller matrix and reject unsupported input explicitly.
- Portable C implementation and optional runtime-dispatched SIMD after portability tests.
- Apply the factor via the existing preconditioner lifecycle and solve hooks; verify ownership, cleanup, and factor invalidation on matrix changes. Do not use an old factor after a failed numerical update.
- No required external linear algebra library, no default solver changes, no public structure layout changes.
- A possible option is `-i fgmres -p ldlt`; exact option registration is subject to maintainer review.

A symmetric indefinite inverse is not necessarily positive definite. Do not advertise unrestricted use with CG or MINRES; the initial supported iterative path is GMRES/FGMRES. Exact factorization may give one outer iteration on suitable problems, but this is not a general guarantee.

## Deferred features

Distributed MPI factorization, complex and extended precision, experimental OpenMP factorization (unresolved ThreadSanitizer findings in the prototype), persistent disk cache, and CalculiX-specific integration. Unsupported configurations must be rejected rather than silently misinterpreted. No changes to default behavior.

## Validation plan

- Compare against dense reference systems with known solutions: SPD, indefinite, 2x2 pivots, delayed pivots, singular inputs and permutation-sensitive systems.
- Verify physical residual, factor reuse and invalidation, allocation limits, failure cleanup and matrix ownership.
- Run existing LIS regression tests, ASan/UBSan, real and complex build compatibility, 32/64-bit index builds, OpenMP and MPI build compatibility (unsupported runtime modes explicitly rejected).
- Benchmark factor setup, repeated apply, total wall time and peak RSS against native ILUTP/complete-fill LU on matched symmetric matrices.
- Validate actual shift-invert / buckling matrices independently before claiming application-level performance or eigenpair certification.

## Current evidence

A separate C11 research prototype (stages 55–70) exists with standalone numerical and sanitizer regression logs. These are not upstream LIS CI results. The standalone source import is complete. Native adapter registration, build integration and upstream testing remain open. No speedup, memory reduction or production readiness is claimed for LIS or CalculiX.

## Question for maintainers

Would an opt-in, bounded, serial symmetric-indefinite factorization preconditioner fit LIS's scope? Feedback on the API surface, naming, and reuse lifecycle is welcome before integrating the imported numerical implementation into LIS.

## First LIS integration checkpoint (application-linked)

The experimental bridge registers `ldltexp` using the existing USERDEF
preconditioner registration API; it does not add an installed header,
change public layouts, or register a built-in preconditioner. It uses the
serial symbolic/numerical multifrontal path (natural ordering), checks
exact symmetry and sorted duplicate-free assembled CSR, and owns a factor
and apply workspace per preconditioner instance. Matrix ownership remains
with the caller. The caller must serialize operations and register once
before solving; no concurrent apply on the same instance is supported.
Generation changes reject a stale factor. No numerical reuse, MPI, complex,
quad, CG/MINRES, or OpenMP numerical factorization is supported at this
checkpoint. Allocation caps cover specific module allocations, not total RSS.
The dedicated workflow links the prototype into LIS and tests physical
residuals on indefinite 2x2 pivot systems with GMRES and FGMRES, plus rejection
of nonsymmetry, singularity and CG. These are pending CI qualification,
not a production-readiness or performance claim. Native built-in registration,
release distribution, larger independent tests and a maintainer API decision
remain future work.

## Second integration checkpoint: bounded quotient minimum degree

The experimental `ldltexp` adapter now uses a quotient graph to represent
elimination cliques implicitly, with exact lazy minimum-degree selection
and deterministic smallest-index ties. The earlier explicit-graph routine
is retained as a small-graph oracle. Strict symbolic supernodes are split
at width 32. Ordering counts adjacency/element visits rather than pair
checks; the two counters are not directly comparable. The bridge input
caps are 200,000 DOF and 10 million full CSR entries. Phase-local ordering,
symbolic and numeric workspace caps are 512, 256 and 1536 MiB respectively;
64 million factor slots, 2 billion ordering visits and maximum numeric
front 8192. These are experimental bounds, not an RSS guarantee. No LIS
public headers, installed ABI, or default solver options change.

Local release and ASan/UBSan tests compare exact permutations, fill and
front size against the explicit graph on random and paired indefinite
cases, including numerical solves and resource rejection. A private
152,592-DOF CalculiX-derived shifted linear system completed ordering
(21,108,420 symbolic off-diagonal L slots, maximum front 1485); a full
numeric solve was not established in the preliminary observation window.
No full-size LDLT success or application-level eigenpair claim is made.

## Third integration checkpoint: SIMD Schur update fast path

For ordinary finite, normal SIMD lanes, the rank-one/rank-two update now
uses vector range and cancellation checks, contiguous lower-triangle stores,
symmetric mirrored stores, and a vector maximum reduction. Mixed or unsafe
vectors retain the existing long-double fallback, as do triangular tails.
The runtime mode selection, scalar reference, pivot policy, workspace and
factor budgets are unchanged. The full-matrix local experiment uses AVX2.

On the same private 152,592-DOF shifted CalculiX linear system and the same
QMD/width-32 path, preliminary single-run setup time decreased from 68.815459
to 46.4829149 seconds; solve time was 0.132695913 seconds, peak RSS 744492 KiB,
22,172,255 numeric factor slots and 15,473 delayed pivots. The physical
relative residual (7.2414170947189105e-15), relative forward error
(4.8939654672876622e-12), and one FGMRES iteration were unchanged in the
recorded precision. These are local observations, not statistically
qualified performance or buckling-eigenpair results.

The local SIMD regression passed 1330/1330 in GCC release and Clang
ASan/UBSan with AVX2 and AVX-512 available. A GCC 14 ASan crash was reproduced
in both the previous and modified AVX-512 paths on the local VM; its cause
is not established. The dedicated CI adds runtime SIMD regression using
GCC release and Clang ASan/UBSan, alongside existing portable bridge and
QMD checks. Runtime SIMD availability depends on the runner hardware.
