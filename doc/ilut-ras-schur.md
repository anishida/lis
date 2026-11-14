# Experimental serial RAS-ILUT and approximate interface Schur

LIS ILUT now has a **separate opt-in** serial restricted additive Schwarz (RAS) variant, inspired by the decomposition architecture of hypre. It does not depend on hypre or use hypre code.

## Configuration

```
-i fgmres -p ilut -scale none -ilut_ras_domains 4 \\
-ilut_ras_overlap 1 -ilut_schur 1 -ilu_pivot_tol 1e-12
```

- `-ilut_ras_overlap 0|1|2` turns on RAS (default: off). Overlap is the number of sparsity-graph edges added beyond each core partition.
- `-ilut_ras_domains 2..16` selects the count of deterministic, contiguous partitions (default: 2). Each subdomain uses the normal LIS ILUT factorization on its extracted local matrix, including the outer `-precon_half` policy.
- `-ilut_schur 0|1` adds an optional **diagonal-approximate interface Schur correction** (default: off). The interface is rows coupled across the original nonoverlapping partitions.

RAS applies the local factor solve on each overlapping subdomain but injects corrections into that subdomain's **nonoverlapping core** only. For the Schur variant, LIS computes

\[\widetilde S=A_{\Gamma\Gamma}-A_{\Gamma I}\operatorname{diag}(A_{II})^{-1}A_{I\Gamma}\]

and uses a small pivoted dense LU factorization to correct the interface residual. It then applies diagonal-approximate interior back-substitution. This is a **one-step approximate Schur correction**, not hypre's distributed GMRES-on-Schur algorithm and not a general-purpose exact Schur complement.

## Exact-interface correctness witnesses (8×8 and 9×9)

Reproducible tests, `ILUT_RAS_EXACT_INTERFACE` (sizes 8 and 9), demonstrate the numerical
purpose of the Schur step **without relying on a Krylov solver**.

Construct `A = 4I` of order 8 or 9 and add symmetric off-diagonal entries
`-1` at (2,3), (3,4), (4,5). Partition into owner domains
`[0,4)` and `[4,n)` (n=8 or 9). The interface is `Γ={3,4}`; all
other rows are interior, with `A_II=4I`. Thus
`diag(A_II)^(-1)` is the **exact** inverse, and the computed
approximate Schur operator equals the true block Schur complement:

```text
S = [[3.75, -1.00],
     [-1.00, 3.75]]
```

Manufacture `x_true=e3+e4` and `b=A*x_true`. Apply the native
LIS RAS-ILUT preconditioner once (no outer iterative solve):

| Option (two domains, overlap 0) | True relative residual |
|---|---:|
| `-ilut_schur 0` | 0.2319003617 |
| `-ilut_schur 1` | zero up to floating-point roundoff |

Both regressions assert an uncorrected residual between 0.20 and 0.26,
then `||Ax-b||/||b|| < 1e-10` and
`max_i |x_i - x_true_i| < 1e-10` for the Schur version. It runs
in the existing `testhalfilut` suite and checks the **actual native
preconditioner application**, rather than just a separate dense-matrix
reimplementation. The table is a mathematically evaluated reference;
the updated native CI must pass independently.

The 9×9 case additionally catches inconsistent interface ownership when
subdomain sizes differ (four rows versus five). The implementation now
uses the same quotient/remainder partition boundaries for core
extraction and owner mapping. These are intentionally exact-special-case
correctness witnesses.
No corresponding Schur exactness, convergence advantage, or FEA
performance benefit is claimed when `A_II` is not diagonal.

## Coupled-interior witness: inexact yet useful Schur correction

The same owner partitions and manufactured solution are also tested
with the symmetric tridiagonal matrix `A=tridiag(-1,4,-1)` for n=8
and n=9. Here `A_II` contains **nonzero interior off-diagonals**,
so the diagonal inverse is genuinely approximate.

Independent double-precision block calculation gives:

| Size | Schur off, relative residual | Schur on, relative residual |
|---|---:|---:|
| 8×8 | 0.2314968813 | 0.0210451710 |
| 9×9 | 0.2314959114 | 0.0210450829 |

The new native regression checks that a **single** Schur-corrected
preconditioner application reduces the true residual by at least 5×
and achieves a relative residual below 0.05. It does not assert exact
solution recovery for this coupled case. This is a limited example
of a useful interface correction, **not** a claim about Krylov iteration
counts, wall time, or real-world contact mechanics.

## Synthetic contact-KKT regression

`testhalfilut` now assembles an indefinite **48×48 stabilized saddle-point system** with a 36×36 spring-stiffness block, twelve displacement-difference constraints, and a negative 0.1 multiplier diagonal. It compares the legacy ILUT control with the opt-in variant, requiring convergence to a **true relative residual below 1e-7** and printing Krylov iterations and residuals. This is an algebraic contact-style regression, **not** a production nonlinear FEA benchmark and not evidence of a speedup.

## Integration and review scope

Upstream [#224](https://github.com/anishida/lis/pull/224) was merged
into `develop` on 2026-11-13. This PR adds the experimental
RAS/Schur path on top of that already-merged ILUT robustness work.

Related PRs: [#226](https://github.com/anishida/lis/pull/226)
(numerical tests), [#227](https://github.com/anishida/lis/pull/227)
(FEA benchmark metrics), [#229](https://github.com/anishida/lis/pull/229)
(ILUTP), [#230](https://github.com/anishida/lis/pull/230)
(two-sided equilibration). PR #229/#230 modify some shared ILUT
and solver-option code; integration must be verified before merging
these features together. This branch enables genuine local ILUTP in RAS when `-ilutp 1` is set, by forwarding the pivoting and full-fill policies to each subdomain. The separate two-sided Equil option from #230 is **not** integrated or validated with RAS.

The sparse-CSR Schur feasibility experiment in
[fork PR #64](https://github.com/3rav/lis/pull/64) is separate and
not integrated; this implementation retains its dense 256-DOF
interface limit.

This change can be reviewed as a narrowly scoped, **experimental,
opt-in serial** solver feature; it is not approved for production FEA
or for merge without review. The exact witness proves algebraic
correctness only in the delimited block structure above.

For a credible FEA performance comparison, use several real Newton
tangent snapshots, matched factor fill, residuals, setup/solve
times and peak memory. The existing
`./test/testhalfilut --fea-matrix path/to/contact.mtx 8 8 0.005`
runner compares only baseline with RCM-scaled ILUT, **not** RAS/Schur.

## Limits / safeguards

- **Experimental and serial-only**: ordinary real double-precision CSR, with row indices local to the single process. MPI, USER/OPERATOR and quad-precision input are explicitly rejected when the option is used; the default ILUT remains unchanged.
- Maximum interface size for the dense Schur option is 256 DOF. A larger interface or singular approximate Schur leads to an explicit error rather than an undocumented silent fallback.
- The diagonal approximation explicitly rejects interior rows with missing or weak diagonals (relative to their largest coefficient, threshold `1e-10`) and nonfinite reciprocals. The dense Schur factorization rejects pivots smaller than `1e-12` of the initial Schur maximum coefficient. These guards prevent silent invalid arithmetic; they **do not guarantee a well-conditioned or effective contact preconditioner**.
- The regression also checks that a deliberately weak interior diagonal is rejected. These tests are algebraic safeguards, not a production FEA validation.
- The nested ILUT setup may be costly in memory and time; there is currently no automatic partitioner or communication overlap.
- The adjoint/transpose RAS application is explicitly **not implemented**. Requiring an adjoint raises an error rather than returning an incorrect result.
- Domain parts own their matrices, factors and work vectors; the outer preconditioner destroys them deterministically.
- The subdomain and Schur apply workspaces are reused per preconditioner instance. Do not call **the same instance concurrently** from multiple threads without external synchronization; independent instances are separate.
- Setup OOM paths initialize the allocated subdomain records before cleanup, so a failure in the second allocation cannot traverse uninitialized pointers.
- Use flexible GMRES or GMRES for testing and inspect true residuals, iterations, setup time, factor fill and maximum RSS. For contact KKT saddle-point matrices, success is **not guaranteed**; validate across several active contact sets and scaling regimes.

The regression `testhalfilut` executes all overlaps 0/1/2 with and without Schur, checking a substantial residual decrease and explicit rejection of unsupported adjoint applies. The current code is a foundation for later **MPI-distributed** RAS plus iterative Schur and ddPQ pivot-matching work, not a claim that these have already been delivered.

Reference: https://hypre.readthedocs.io/en/latest/solvers-ilu.html

## Local ILUTP integration (stacked on upstream #229)

This RAS/Schur branch now includes the exact native ILUTP commits from
upstream PR #229 and adds genuine local ILUTP support for RAS. With `-ilutp 1`, each overlapping
subdomain builds its own pivoted ILUTP factorization instead of silently
falling back to ILUT. The local options `-ilutp_threshold`,
`-ilutp_full`, `-ilutp_full_limit`, RCM and row scaling are
forwarded from the outer solver.

For example:

```sh
-i fgmres -p ilut -scale none -ilut_ras_domains 2 \\
-ilut_ras_overlap 1 -ilut_schur 1 -ilutp 1 \\
-ilutp_threshold 0.5 -ilut_drop 0.005 -ilut_rate 8
```

The integration regression `ILUT_RAS_LOCAL_ILUTP` uses independent
zero-diagonal 2x2 sparse blocks to force dynamic column pivoting inside
RAS subdomains, with `-ilutp_full 1`. It checks both the original-system
residual and a manufactured exact solution, and rejects full-fill when
`-ilutp 0`. This does **not** mean other algorithm combinations or
production FEA problems are automatically safe or efficient.

ILUTP requires **one OpenMP solver thread**, serial real-double CSR,
and half-precision factor storage **disabled**. RAS retains its own
serial-only and Schur interface-size limits.

**Merge order:** merge upstream #229 into `develop` **before** #225.
Because #225 already contains the same #229 commits as a stacked
branch, its comparison should then narrow to the RAS/Schur-specific
changes. Confirm the new comparison, full CI and merge status again
before merging #225; neither PR should be merged out of order.
