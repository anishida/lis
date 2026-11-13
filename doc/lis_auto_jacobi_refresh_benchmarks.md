# AUTO cached-CSR Jacobi numeric refresh: regression and performance data

**PR:** [anishida/lis#167](https://github.com/anishida/lis/pull/167)

**Validation:** [8/8 PASS, 2026-11-01](https://github.com/3rav/lis/actions/runs/41286231861)

## What changed

On a value-only change with an explicitly opted-in, compatible Jacobi
preconditioner, the AUTO session used to verify the complete CSR row pattern
**twice**: once to decide whether Jacobi could survive, and again inside
generic cached CSR materialization.

The updated implementation verifies the actual row pattern once and passes
a stack-scoped proof to the *immediately following* cache refresh. It checks
the source and cache identity, CSR pointer, cache epoch and effective generation
stamps before consuming the proof. The public cached materialization API
continues to perform its own pattern verification. A failed materialization
invalidates the partially updated CSR and drops persistent solver state.

A benchmark invariant now requires equal numbers of USER row callbacks under
SAFE and JACOBI_NUMERIC policies for these controlled scenarios.

## Setup

- Three-diagonal USER row-callback matrix, dimension 3,000 / 30,000 / 100,000.
- GMRES, Jacobi, `-scale none`, zero initial vector at each update.
- Recompute RHS so the exact solution remains the all-ones vector.
- Value-only matrix updates: respectively 300 / 80 / 30 repetitions.
- `-O2 -g` builds; three runs per size and build, reversing policy order
  for the middle run.
- CPU seconds measured with C `clock()`, encompassing the source update,
  cached materialization, solver work and correctness verification.
- Linux GitHub-hosted runners; results are not deterministic dedicated-host
  microbenchmarks and do not isolate preparation time from solve time.

## Measurements

Median CPU seconds from three runs; positive change means Jacobi is slower.
A negative change means Jacobi is faster.

| Configuration / rows / updates | SAFE reset (s) | JACOBI_NUMERIC (s) | Change |
|---|---:|---:|---:|
| Real i4 / 3,000 / 300 | 0.241757 | 0.242247 | +0.2% |
| Real i4 / 30,000 / 80 | 0.616500 | 0.620120 | +0.6% |
| Real i4 / 100,000 / 30 | 0.776218 | 0.780283 | +0.5% |
| Complex i4 / 3,000 / 300 | 0.552973 | 0.565542 | +2.3% |
| Complex i4 / 30,000 / 80 | 1.416802 | 1.410896 | -0.4% |
| Complex i4 / 100,000 / 30 | 1.802347 | 1.806771 | +0.2% |
| Real OpenMP / 3,000 / 300 | 0.379167 | 0.296136 | -21.9% |
| Real OpenMP / 30,000 / 80 | 0.638296 | 0.480172 | -24.8% |
| Real OpenMP / 100,000 / 30 | 0.571398 | 0.562175 | -1.6% |

### USER row callback invocations

Both policies now perform exactly the same number of row callback invocations
in the measured cases:

| Rows / updates | SAFE | Jacobi numeric |
|---|---:|---:|
| 3,000 / 300 | 3,600,000 | 3,600,000 |
| 30,000 / 80 | 9,600,000 | 9,600,000 |
| 100,000 / 30 | 12,000,000 | 12,000,000 |

In the previous implementation, the Jacobi path performed 50% more callbacks
for these value-only USER updates.

## Regressions

The validation matrix passed 8/8 configurations:

- Real 32-bit and 64-bit LIS integer types
- Complex 32-bit and 64-bit LIS integer types
- Real long double
- Real OpenMP
- Complex ASan/UBSan
- MPI-enabled real, one rank only

Regressions include normal numeric refresh, preconditioner reuse, CSR structural
mismatch, **underreported** structural changes (incorrectly marked values-only),
callback failure after partial value writes, atomic outputs, and subsequent
recovery through rebuild.

The validation workflow and configurable benchmark variant live on a separate
fork branch and are intentionally not part of the PR diff. The production
benchmark in this PR is a fixed 3,000-row smoke benchmark.

## Interpretation and remaining limitations

The eliminated duplicate scan removes a known, measurable regression. Serial
workloads here show essentially neutral performance; some OpenMP workloads
show improvement, but the benefit is **not universal**. This is an opt-in
optimization with a limited workload-specific rationale, not a claim of general
speedup.

No distributed multi-rank MPI tests, full CalculiX integration performance
results, or independent phase-by-phase timings were obtained. Further external
review should focus on the scoped proof's invalidation guarantees and whether
the additional AUTO Jacobi-specific complexity is worth maintaining.
