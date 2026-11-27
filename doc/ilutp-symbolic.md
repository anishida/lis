# Opt-in symbolic complete-fill LU

This proposal extracts the recovered RC32 scalar symbolic factorization,
exact minimum-degree ordering, direct permuted CSR assembly and numerical
refresh. It does not import the historical benchmark archive or GLI patches.

Configuration (all new options default to 0, accept only 0/1):

    -i fgmres -p ilut -ilutp 1 -ilutp_full 1 -ilutp_threshold 0
    -ilutp_symbolic 1 -ilutp_md 1 -ilutp_reuse 1

Symbolic factors require a sorted, duplicate-free, structurally symmetric
CSR matrix in a non-MPI real-double int32 build, no scaling, half storage,
RAS, Schur, or additive Schwarz. Numerical factorization uses the original
scalar row kernel and does not dynamically pivot: symbolic mode requires
zero column-pivot threshold. Dynamic column-pivot permutations are not
cached or refreshed by this path. A zero or unusable numeric pivot fails
rather than altering the symbolic graph. Explicit symbolic
selection rejects unsupported inputs instead of silently selecting a
different factorization. OpenMP outer Krylov operations may use threads;
triangular substitutions remain serial. Minimum degree is optional;
reuse requires minimum degree and frozen-pattern sequence refresh policy.

The cache belongs to the preconditioner, not application user_data.
A factor validity flag is set only on successful construction/refresh;
failed refresh invalidates the factor and the sequence caller destroys it
before rebuilding. Solver sequence configuration snapshots include all new
options. No environment variable controls the symbolic path.

The new fields are appended to the public solver/preconditioner structs:
existing field offsets are preserved but sizeof changes. This requires
ABI review and rebuilding clients, not a claim of binary compatibility.

The grouped 2/3/4/8 kernels, DOT4, RUNS and parallel numerical variants
remain in the historical RC32 branch and are not included here.
No performance claim for the private CCX fixture is made.
The focused CI compares known solutions, fresh factors and sequence
solutions, checks state and user-data isolation, and exercises failure
recovery. Broader structural patterns and independent performance validation are still
needed before production acceptance.

## Synthetic grid scaling (informational)

The focused test generates sorted, structurally symmetric five-point CSR
grids with nonsymmetric, strictly diagonally dominant values. Six systems
per mode have the same pattern, slightly changing values and known
solutions. All four modes use FGMRES, threshold zero, the same fill cap
per grid size, a zero initial guess and tolerance 1e-12. The test checks
physical relative residual < 1e-11, relative forward error < 1e-10 and
actual factor creation/refresh counts.

Example optimized serial runner results (commit `5fce708`,
[CI log](https://github.com/3rav/lis/actions/runs/46005530008/job/140527601101)):

| DOF | Mode | First solve (ms) | Mean next 5 (ms) | L+U off-diagonal entries |
| ---: | --- | ---: | ---: | ---: |
| 1600 | Legacy | 7.881 | 6.810 | 124878 |
| 1600 | Symbolic | 2.843 | 3.080 | 124878 |
| 1600 | Symbolic MD | 2.791 | 2.861 | 39808 |
| 1600 | Symbolic MD reuse | 2.769 | 0.482 | 39808 |
| 6400 | Legacy | 104.017 | 97.504 | 1011358 |
| 6400 | Symbolic | 31.871 | 32.308 | 1011358 |
| 6400 | Symbolic MD | 15.192 | 15.592 | 242400 |
| 6400 | Symbolic MD reuse | 15.309 | 4.674 | 242400 |
| 25600 | Legacy | 1671.699 | 1642.695 | 8141118 |
| 25600 | Symbolic | 600.311 | 578.190 | 8141118 |
| 25600 | Symbolic MD | 140.821 | 138.145 | 1429828 |
| 25600 | Symbolic MD reuse | 139.902 | 81.290 | 1429828 |

At 25600 DOF the reuse mode took about 0.546 s for six solves vs
9.885 s for legacy (18.09x); factor entry count fell by 82.44%.
All tested cases passed; maximum relative residual and forward error
were below 1.46e-13. Reuse performed one creation and five numerical
refreshes; the other modes performed six creations.

These are single-run illustrative timings, not stable performance
estimates. For example, two optimized CI runs at 6400 DOF yielded
20.86x and 39.28x for mean subsequent solves. No performance threshold
is used in CI. Entry counts are not peak process memory; do not
extrapolate these synthetic results to CCX or general indefinite systems.

After building the focused test as in the CI workflow, measure each
case in a fresh process (mode 0 = legacy, 1 = symbolic, 2 = symbolic
MD, 3 = symbolic MD reuse; width 40, 80 or 160):

```sh
for rep in 0 1 2; do
  for offset in 0 1 2 3; do
    mode=$(( (rep + offset) % 4 ))
    OMP_NUM_THREADS=1 /usr/bin/time -f 'elapsed_s=%e peak_rss_kib=%M' \
      ./sequence-audit --benchmark 160 "$mode"
  done
done
```

### Repeated process measurements (25600 DOF)

Three independent processes per mode, rotating the mode order, completed
on the same optimized serial CI runner at commit `becf65a`
([CI log](https://github.com/3rav/lis/actions/runs/46006096251/job/140529509005)).
The table reports medians; brackets give observed min/max values.

| Mode | Sum of six solve times (s) | Mean next five (ms) | Process wall time (s) | Peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: |
| Legacy | 7.762 | 1284.143 | 7.78 [7.75, 7.79] | 115548 [115408, 115756] |
| Symbolic | 2.019 | 331.354 | 2.03 [2.02, 2.05] | 115272 [115256, 115272] |
| Symbolic MD | 0.520 | 85.861 | 0.53 [0.53, 0.56] | 36148 [36076, 36176] |
| Symbolic MD reuse | 0.296 | 40.752 | 0.30 [0.30, 0.32] | 35460 [35392, 35480] |

Ratios of medians (legacy / MD reuse): 26.25x for the sum of solve
times, 31.51x for the mean of the five subsequent solves, and 25.93x
for process wall time (rounded to 0.01 s by GNU time). Median peak
RSS decreased by 69.31%, from 112.84 to 34.63 MiB. The maximum
physical relative residual and forward error in all twelve processes
were below 1.46e-13. Each reuse process performed one creation and
five refreshes; each other process performed six creations.

These measurements share a single CI runner; they do not quantify
between-machine variability or establish performance on CCX. The
runner CPU model was not explicitly recorded in this historical log.
Use an optimized non-sanitized serial build, record the CPU and compiler,
and report medians/ranges across repeated processes. Peak RSS includes
all process allocations, not just LU factors. The private CCX fixture
and ABI acceptance remain outstanding.

## Release distribution

The internal ILUTP, minimum-degree and refactor includes are listed in
`src/precon/Makefile.am` and its tracked `Makefile.in`. The focused test
and this document are also in the top-level distribution manifest. The
serial CI checks archive membership and runs `make distcheck` with the
non-MPI real-double test configuration; no new installed header or
Fortran option index is needed for these internal includes.

## ABI migration gate

The serial CI compiles a matching client and a client using the pinned
upstream baseline headers, both linked to the proposed library. Only the
fixed-width ABI query is called by the baseline client; no old-layout
solver/preconditioner object is passed to the new library. The gate requires
rejection of the old layout, acceptance of the matching layout and unchanged
selected existing member offsets. This is a migration screening test, not
proof of binary compatibility or a substitute for maintainer ABI policy.
Clients must rebuild. The option remains opt-in; changing public structure
sizes still requires an explicit release/ABI decision.
