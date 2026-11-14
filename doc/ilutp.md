# Native ILUTP and bounded complete-fill LU (experimental)

## What this PR delivers

An **optional, native LIS preconditioner** for serial real-double CSR
systems. The existing `-p ilut` factorization remains the default
when `-ilutp 0` (the default). No SuperLU library or plugin is needed.

- `-ilutp 1`: incomplete LU with **dynamic column pivoting** and
  `-ilutp_threshold` in `[0,1]` (default `0.5`).
  **Compatibility boundary:** in limited-fill mode, threshold `0` runs
  the exact legacy ILUT factorization with an identity column permutation.
  Any *positive* threshold uses the native ILUTP kernel, even if no column
  swap occurs. Its magnitude-based fill retention differs from legacy
  ILUT, so truncated L/U factors are not guaranteed to converge to the
  legacy factors as the threshold approaches zero. This is a deliberate
  mode boundary, not a guarantee of coefficient continuity.
- `-ilutp_full 1`: **bounded complete-fill LU**, without dropping
  **or diagonal pivot regularization**, using the same column-pivot
  selection logic. This changes both numerical policies at once.
- `-ilutp_full_limit N`: cap on stored off-diagonal L+U entries
  (not resident memory). A value of 0 uses
  `max(20*nnz(A), 4*n)`; exceeding the cap fails setup.
- `-ilut_rcm 1` and `-ilut_row_scale 1`: optional ordering and
  max-entry row scaling inside the preconditioner; GMRES still
  solves the **original** `Ax=b`.

## Optional two-sided ILUT equilibration

`-ilut_equil 1` enables two-sided maximum-absolute-value scaling in
serial real-double CSR ILUT, including the dynamic-column-pivoting and
bounded complete-fill variants. It is **disabled by default**.
The option applies one row pass and then one column pass to the
row-scaled matrix. `-ilut_row_scale 1` remains available as an
independent row-only baseline; when Equil is enabled, the row pass
is performed only once. `-ilut_rcm 1` may be combined with either.

With symmetric ordering `P`, diagonal row weights `R`, diagonal
column weights `C`, and dynamic **column** permutation `Q`:

```text
B = P A P^T
S = R B C
S Q ~= L U
forward: P^T C Q (LU)^-1 R P
adjoint: P^T R^H (LU)^-H Q^T C^H P
```

Importantly, `C` is indexed in the **pre-pivot** column coordinates,
so the outgoing forward and incoming adjoint factors use `C[Q[i]]`.
The original Krylov operator `A` remains unchanged. Row/column
maximum passes take `O(nnz)` and the column weights cost `O(n)`
additional memory. Empty rows/columns retain unit scaling; nonfinite
input or transformed coefficients are rejected. Equil is not
a substitute for numerical pivoting, and it can increase incomplete
factor fill.

The nested solver supports `-p hybrid -hybrid_p ilut -ilut_equil 1`
with and without `-ilutp 1`. Equil with a non-ILUT nested
preconditioner is rejected rather than silently ignored. The
combined path requires serial CSR input and one active OpenMP thread
when ILUTP is enabled.

To compare Equil with its controls on the **same original** sparse
Matrix Market matrix and RHS, use:

```sh
./test/benchilutp A.mtx b.mtx --equil 0 --rcm 1 --row-scale 1 --csv row.csv
./test/benchilutp A.mtx b.mtx --equil 1 --rcm 1 --row-scale 1 --csv equil.csv
```

The `equil` CSV column records the requested mode. The
`physical_relres` field is measured against the original `A`
and `b`. `tools/run_native_ilutp_contact.py` also accepts
`--equil 0|1` and forwards it to this benchmark.
For a small manufactured-system control with five profiles, run
`./test/benchilutpequil 30 5 0.005`.

Regression coverage is provided by `testilutequil` and
`testilutp`, including forward/transpose solves, sequence
invalidations, Hybrid forwarding and bounded complete-fill.
Neither lower GMRES iteration counts nor a reduction in pivot
swaps alone establishes total-time or memory superiority.

### Reproducible Stokes-type indefinite saddle point

The included `test/gen_stokes_saddle.py` manufactures a nonsingular,
**symmetric-indefinite** pressure/velocity saddle matrix with *exactly
zero* pressure diagonal. This is a Stokes-inspired finite-difference
model, **not** a real IFISS, PETSc or CalculiX export.

```text
H = [ 0    B  ]
    [ B^T  K  ]
A = D H D
```

`K` is a positive-definite five-point velocity stiffness operator,
`B` a full-row-rank discrete divergence-like constraint matrix,
and the positive diagonal `D` introduces a manufactured unit contrast.
On an `m × m` grid, `A` has `3m²` unknowns, `m²` zero
pressure pivots, `2m²` positive and `m²` negative eigenvalues.

```sh
python3 test/gen_stokes_saddle.py --side 8 --contrast 30 \
  --units mixed_checker --ordering interleaved --gamma 0.3 \
  --delta 0 --matrix stokes-A.mtx --rhs stokes-b.mtx

./test/benchilutp stokes-A.mtx stokes-b.mtx \
  --maxiter 180 --restart 70 --rtol 1e-6 --drop 0.005 \
  --rate 3 --pivot-tol 0.001 --rcm 0 --row-scale 0 \
  --equil 0 --csv no-scale.csv

./test/benchilutp stokes-A.mtx stokes-b.mtx \
  --maxiter 180 --restart 70 --rtol 1e-6 --drop 0.005 \
  --rate 3 --pivot-tol 0.001 --rcm 0 --row-scale 1 \
  --equil 0 --csv row-only.csv

./test/benchilutp stokes-A.mtx stokes-b.mtx \
  --maxiter 180 --restart 70 --rtol 1e-6 --drop 0.005 \
  --rate 3 --pivot-tol 0.001 --rcm 0 --row-scale 0 \
  --equil 1 --csv two-sided.csv
```

For the native LIS n=192 example, the following five controls all
exhausted 180 GMRES iterations: ordinary ILUT, row-scaled ILUT,
unscaled ILUTP, row-scaled ILUTP and Equil-only ILUT. By contrast,
**ILUTP+Equil converged in 9 iterations**, with original-system
relative residual **4.78e-7** and 2889 stored L+U off-diagonal
coefficients (versus at most 2734 for the controls). These controls
use the same `-ilut_rate 3` and `-ilut_drop 0.005` settings.

At n=300 (`--side 10 --gamma 1`, other parameters unchanged)
the same five controls failed at 180, while **ILUTP+Equil
converged in 11 iterations** and used 4575 L+U nonzeros (the
largest control: 4318). Original-system relative residual was
**3.14e-7**. Independent manufactured exact-solution comparisons
also confirmed relative solution errors **5.38e-4** and
**1.36e-4**, respectively, against an acceptance threshold of
`1e-3`. That optional accuracy audit resides on the experimental
branch and is not part of the generic Matrix Market benchmark.

These results demonstrate actual **numerical complementarity**:
neither pivoting nor two-sided Equil (nor row-only scaling) suffices
alone for these *selected* saddle systems. They are not claims
about general wall-clock speed, memory optimality or a real
Stokes/contact FEM dataset. The sweep, raw controls and accuracy
check are documented in the
[reproducible experiment](https://github.com/3rav/lis/blob/experiment/stokes-saddle-ilutp-equil-20261114/doc/stokes-saddle-equil-benchmark.md).

### Example: complete-fill mode

```sh
-i fgmres -p ilut -ilutp 1 -ilutp_full 1 \
-ilutp_threshold 0.1 -ilut_rcm 1 -ilut_row_scale 1 \
-ilutp_full_limit 60000000 -omp_num_threads 1 -scale none
```

On an OpenMP-enabled build, **set `OMP_NUM_THREADS=1` in the environment before creating the preconditioner**, including for direct `lis_precon_create` or `benchilutp` use. The solver option `-omp_num_threads 1` alone may not change the runtime's thread count before preconditioner setup.

This is **column** pivoting, not row-pivoted or supernodal LU.
Supported scope is serial real-double CSR with one active
OpenMP thread. MPI, complex scalars and half-precision *factor
storage* are unsupported for opt-in ILUTP; legacy ILUT's half mode
is still available. A `LIS_MATRIX_USER` or operator shell may still
use legacy ILUT through its row callback, but an explicit ILUTP
request on a non-CSR matrix **fails**, never silently falls back.

Hybrid can select the new ILUTP **on explicit CSR** using
`-p hybrid -hybrid_p ilut -ilutp 1` (with the usual inner
`-hybrid_i`, tolerances and optional `-ilutp_full 1`).
The outer opt-in ILUTP numerical parameters are passed to the
inner ILUT preconditioner. Hybrid + USER and legacy half storage
continue to work unchanged without ILUTP; unsupported combinations
with opt-in ILUTP return an error instead of using another
algorithm. Frozen-pattern factor updates and factor reuse across
different Newton matrices are not validated.

## What has been demonstrated

| Result | Observed |
| --- | --- |
| Real CCX Newton linear systems | **15/15** right-hand sides solved |
| Independent matrices in that set | **12** (four steps share one A) |
| FGMRES iterations | **1 per solve** |
| True original-system relative residual | **4.63e-15 to 2.84e-12** |
| Largest stored L+U | **51,108,311** off-diagonal entries |
| Peak observed process RSS | **799,024 KiB (~780 MiB)** |
| Longest factor setup | **65.04 CPU seconds** |

Configuration: RCM, row-max scaling, full fill, threshold 0.1,
60M factor-entry cap, FGMRES(40), serial double precision.
Every reported right-hand side was the **original** exported RHS,
and factorization was recomputed for every case. The private input
series is one FEA model with **unknown spring count**.

See [results, reproducibility and remaining evidence gap](ilutp_contact_ccx_v2_benchmarks.md).

## What is not yet demonstrated

The data do **not** establish success on the separately identified
**256-spring snapshot 03**, full end-to-end CalculiX convergence,
or an advantage over the production CCX solver. In particular,
the successful complete-fill configuration disables both dropping
and diagonal regularization; their contributions have not been
separated. A controlled test with `--pivot-tol 0` on the original
256-spring A/b pair remains outstanding.

A stored-factor limit is **not** a total process memory guarantee.

## Tests and compatibility

```sh
autoreconf -fiv
./configure --enable-test --disable-fortran --disable-f90 --disable-saamg --disable-mpi
make -j2
./test/testilutp
make -j2 check
```

`testilutp` covers column permutations, forward/transpose solves,
RCM/scale, fill limits, option replacement, Hybrid inner ILUTP,
half-mode rejection, invalid configurations and sequence behavior.
`testuserprecon` checks legacy USER/Hybrid support and fail-closed
rejection of opt-in ILUTP on USER. Dedicated native CI runs these
tests with half precision enabled and also uses ASan/LSan/UBSan
to exercise **RCM+scale after a partially allocated full-fill
factor exceeds its cap**.

**Legacy compatibility note:** while legacy ILUT's default
factorization is unchanged, this PR makes
`lis_solver_set_option_replace` and
`lis_solver_set_option_replace_if_changed` also propagate/compare
the existing `ilut_rcm` and `ilut_row_scale` settings. That can
affect legacy ILUT when options are replaced on an existing solver.

**Upstream review status:** the PR is suitable for **Ready for
review**, with current-head Linux/MSYS2 CI and dedicated native ILUTP
sanitizer/USER/half/Hybrid regressions green as of 2026-11-14.
Linux jobs are pinned to Ubuntu 24.04 after a controlled baseline/PR
comparison reproduced three unrelated eigensolver test failures
on Ubuntu 26.04 in both trees. This is a CI reproducibility
workaround, **not** a proven GCC 15 numerical fix.

Reviewers should still inspect the option-replacement and partial
allocation cleanup paths. The original 256-spring experiment and
end-to-end nonlinear CCX validation remain outstanding; neither
success should be claimed at review time.
