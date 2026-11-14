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
