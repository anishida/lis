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
