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

## Fourth integration checkpoint: finite maximum fast path

In numeric front assembly, pivot selection and multiplier statistics,
all maximum operands are finite nonnegative magnitudes after input and
range checks. Use an inline comparison instead of repeated general-purpose
long-double libm maximum calls. This does not change pivot acceptance,
delayed-pivot policy, arithmetic precision, range checks or resource limits.

A preliminary local single-run observation on the same 152,592-DOF
CalculiX-derived shifted linear system decreased setup from 46.4829149 s
to 40.4585161 s. Numeric factor slots (22,172,255), delayed pivots (15,473),
physical relative residual (7.2414170947189105e-15), forward error
(4.8939654672876622e-12) and one FGMRES iteration were unchanged in the
recorded precision. Peak RSS was 744636 KiB. The comparison is not a
statistically qualified timing result; the complete buckling eigensolver
was not tested. Local GCC release and Clang ASan/UBSan regressions each
passed multifrontal 543/543, supernodal 783/783 and SIMD 1330/1330.
The dedicated CI now runs the first two suites in addition to the existing
bridge, QMD and SIMD checks.

A separate local power-of-two symmetric equilibration experiment gave
41.0708539 s and 15,416 delays, with physical residual 7.44e-15. The
reduction in delayed pivots was negligible, so this scaling policy is
not included in the PR.

## Fifth integration checkpoint: tiled Schur updates

The AVX2 and AVX-512 paths pack immutable pivot couplings into fixed
64-column tiles and reuse contiguous SIMD loads in place of repeated
strided gathers. The tile width is divisible by both vector lengths,
preserving vector group boundaries and each lane's arithmetic. The
reference kernel, long-double fallback, pivot acceptance, growth checks
and factor budgets are unchanged. Scratch storage is two 64-double arrays
(1 KiB per active kernel). This is per-pivot cache tiling, not a delayed
multi-pivot BLAS-3 update.

On the same private 152,592-DOF shifted CalculiX system, sequential local
single-run A/B observations with the same native LIS build gave 40.9835579 s
(before) and 38.1748061 s (tiled) setup, a 6.9% reduction. Peak RSS was
744592 and 744524 KiB; both runs used 22,172,255 numeric factor slots,
15,473 delayed pivots, one FGMRES iteration, physical relative residual
7.2414170947189105e-15 and forward error 4.8939654672876622e-12 in the
recorded precision. No statistically qualified speedup or eigenpair claim
is made. The added tests exercise multiple tile boundaries, rank-one
and rank-two updates and all SIMD remainders at nonzero pivot offsets.
Local GCC release and Clang ASan/UBSan passed SIMD 1382/1382,
multifrontal 543/543 and supernodal 783/783.

## Sixth integration checkpoint: packed symmetric child updates

Child-to-parent contributions now retain only the packed lower triangle
rather than an entire square symmetric update. Checked triangular-size
arithmetic precedes the allocation, and the assembly order, active
front storage, pivot policy and accumulation precision are unchanged.
The pivot column maximum is read from the contiguous symmetric active
row instead of a strided column. No additional public API or ABI changes.

A 32-child, 24-variable shared-separator fixture completed with 108,848
bytes of numeric-module workspace instead of 179,504 bytes (39.4% lower).
The regression constrains workspace below the lower bound for the old
unpacked child buffers alone, 147,456 bytes. The preceding implementation
returns P60_LIMIT with that constraint; the packed implementation solves
and verifies the same system. Existing resource-failure atomicity tests
are retained.

On the same private 152,592-DOF CalculiX-derived shifted system, the
preliminary local setup observation was 37.165884 s (previous 38.1748061 s),
solve 0.135046005 s, peak RSS 716832 KiB (previous 744524 KiB), and numeric
workspace peak 481113856 bytes (previous 493144800). The recorded
physical residual (7.2414170947189105e-15), forward error
(4.8939654672876622e-12), factor slots (22,172,255), delays (15,473)
and one FGMRES iteration are unchanged. Timing observations are not
statistically qualified; they do not establish eigenpair accuracy.
Local GCC release and Clang ASan/UBSan passed multifrontal 543/543
(including the new shared-separator budget test), supernodal 783/783,
and SIMD 1382/1382. The number of workspace-failure cases depends on
the measured peak, so aggregate counts can change with storage layout.

## Seventh integration checkpoint: lower-triangle active front

The active numeric front now treats the lower triangle as authoritative.
Schur updates avoid redundant upper-triangle stores, and assembly and
symmetric interchanges use the corresponding lower-triangle addresses.
Eliminated pivot rows still hold their multipliers; the allocation shape
and workspace accounting are unchanged. The mirrored kernel remains
available as an internal reference. Regression tests poison the unused
active upper triangle with NaNs and check bitwise equality of the lower
triangle and peak magnitude against the mirrored kernel in each available
mode. The pivot policy and arithmetic range/cancellation checks are unchanged.

Local GCC release and Clang ASan/UBSan each passed SIMD 2732/2732,
multifrontal 543/543 and supernodal 783/783. A private full-size shifted
CalculiX-derived system (152,592 DOF) completed in 31.8586659 s setup and
0.149282932 s solve, versus the previous 36.334445 s single-run AVX2
observation. Peak RSS was 713396 KiB. Numeric factor slots (22,172,255),
delays (15,473), one FGMRES iteration, physical relative residual
(7.2414170947189105e-15) and forward error (4.8939654672876622e-12)
were unchanged in the recorded precision. Compilation/regression work
overlapped part of the new timing; this is preliminary, not a controlled
multi-run speedup measurement. Forcing AVX-512 in the previous implementation
gave 53.1606159 s setup and a certified physical residual of 7.43e-15;
the conservative AVX2 default is unchanged. No full eigenpair claim is made.

## Eighth integration checkpoint: compile-time SIMD specialization

The existing AVX2 and AVX-512 kernels now force compile-time specialization
for pivot width (1 or 2) and storage policy (mirrored or lower-only).
GCC -O2 previously emitted a common implementation despite the inline
specifier; the new wrappers remove runtime dispatch from the hot loops.
No valid pivot update arithmetic, range/cancellation fallback, growth
policy, or allocation bounds change. Signed structural-zero cases are
added to SIMD regression.

Local exact A/B comparison against the previous implementation passed
2088/2088 cases for both available instruction sets and storage policies.
GCC release and Clang ASan/UBSan passed SIMD 2988/2988, multifrontal
543/543 and supernodal 783/783 each. The portable SIMD-disabled build
passed 612/612. A separate five-repetition, alternating-order 1024-by-1024
AVX2 microbenchmark (32 repeated updates per sample, bitwise equal results)
recorded median times of 0.028237334 versus 0.026548096 s for rank one
and 0.049245926 versus 0.046496657 s for rank two (6.0% and 5.6% lower).

On the private 152,592-DOF shifted CalculiX-derived linear system,
sequential single-run setup was 32.3367729 s (control) and 31.425508 s
(specialized), peak RSS 713260 versus 713152 KiB. Both used 22,172,255
numeric factor slots, 15,473 delayed pivots, one FGMRES iteration,
physical relative residual 7.2414170947189105e-15 and forward error
4.8939654672876622e-12 in the recorded precision. The full-size timing
is a preliminary observation; it is not a multi-run speedup claim.

Two alternative local structural-zero SIMD fast paths gave 33.8219349 s
and 32.406774 s, respectively. Neither demonstrated a useful full-size
benefit; neither is included in this change. The new tests retain coverage
of signed zeros without modifying the existing numerical fallback.
