# Reproducible ILUT Matrix Market review benchmark

This change extends the existing optional **`test/testhalfilut --fea-matrix`**
harness. It does not add a new solver or change default ILUT behavior.

## Build and run

The existing test executable is built with the LIS half-enabled test profile:

```sh
./configure --enable-test --enable-half --enable-omp --disable-mpi
make -j2
OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE \
  LIS_FEA_BASELINE_THREADS=8 LIS_FEA_RCM_THREADS=8 \
  ./test/testhalfilut --fea-matrix test/testmat0.mtx 8 8 0.005
```

The original CLI `--fea-matrix MATRIX BASE_RATE RCM_RATE [DROP]` and output
variant names `baseline` and `rcm_scaled` are preserved. The optional
environment variables allow per-variant **OpenMP factorization and
application** thread counts; if unset, each variant uses the process's
existing OpenMP thread count. In a non-OpenMP build, only a thread count of
1 is accepted.

The harness constructs a deterministic known vector `x_true` and sets
`b = A x_true`. Both variants operate on the original imported matrix,
not a reordered or row-scaled surrogate. Reported **true relative
residual** is `||Ax-b||_2 / ||b||_2` on that original system.
`rcm_scaled` enables the opt-in ILUT RCM ordering and row equilibration
already provided by upstream. The test is limited to single-process
real-valued CSR matrices.

## Thread-count confounding (important)

**Use equal thread counts in both variants if the question is whether RCM +
row scaling helps the existing ILUT at a fixed parallel configuration.**
LIS's OpenMP ILUT is partitioned by thread: its create path uses
`omp_get_max_threads()` and thread-specific row ranges; the apply path
uses the corresponding ranges. Consequently, changing a factor's thread
count can change *its L/U entries, factor fill and Krylov iteration count*.
The tool deliberately supports unequal thread budgets to study
thread-dependent algorithm behavior, but such a run cannot isolate the
effect of RCM or row scaling alone. The `threads` CSV column is required
to interpret each result; `fair_fill=yes` does not eliminate this
algorithmic confounder.

For a separate thread-sensitivity experiment, hold the matrix, ordering,
drop and rate fixed and run the RCM variant twice with
`LIS_FEA_RCM_THREADS=1` and `=2`, while holding the baseline at one
thread. Compare `nnz_LU`, `iterations`, true residual, and wall times.
Never compare their timings without also considering thread budget.

This is not a demonstration of a speedup, nor evidence that a different
factor-count proxy implies equal solver quality.

### Controlled same-run demonstration

We tested `test/testmat0.mtx` (100 unknowns, 460 imported entries)
twice **in the same OpenMP GitHub Actions job**, holding the original
matrix, RHS, drop rate, ordering, row scaling and baseline thread budget
fixed. Only `LIS_FEA_RCM_THREADS` changed:

| RCM+scaled ILUT threads | Stored L+U entries | FGMRES iterations | Original-system true residual |
|---:|---:|---:|---:|
| 1 | 968 | 5 | 1.15692696e-9 |
| 2 | 776 | 17 | 9.33675937e-9 |

Source: https://github.com/3rav/lis/actions/runs/43419495035

This demonstrates that the **thread partition changes the numerical
preconditioner**. A 1-vs-2 thread comparison must *not* be called a pure
RCM/row-equilibration comparison. Both solves met `1e-8` here, but
the small bundled test is only a regression fixture, not a structural
FEA speedup result. The comparison intentionally does not assert that
all future implementations must exhibit the same thread dependence.


## Interpreting the results

CSV contains the legacy columns and additionally:

- `setup_wall_s`, `solve_wall_s`, `total_wall_s`: **wall-clock**
  measurements, unlike the existing CPU-time counters. These timings are
  approximate and can vary with load and OpenMP scheduling.
- `threads`: actual thread budget selected for each variant.
- `target_1e-8_met`, `solver_error`, `status`: distinguish a completed
  solver run from satisfying the configured `1e-8` tolerance.

**Both variant rows are printed even if either fails**. Nonconvergence
returns a nonzero process status when its true residual exceeds the legacy
`1e-7` smoke acceptance band, but a residual between `1e-8` and
`1e-7` is still explicitly reported as `not_converged` rather than
misrepresented as meeting the target. File and matrix-import failures
are reported separately.

`fair_fill=yes` means the stored L/U **entry-count proxy** differs by at
most 5%; it does **not** mean equal peak RSS. If factor statistics are
missing, fill matching is reported as `unknown`, never guessed. For fair
comparisons, tune `BASE_RATE` and `RCM_RATE` independently and measure
multiple runs. A green CI smoke check is not a performance claim.

No private contact or KKT data are bundled, and a synthetic or standard
structural matrix is not evidence of contact-KKT performance.
