# MFAMG optional FP16 preconditioner storage

MFAMG honors the existing LIS solver `-precon_half` setting. It reuses
the IEEE 754 binary16 conversion routines used by ILU.

* `off`: original full-precision storage and computation.
* `on`: compress suitable **prolongation coefficients**, **block-Jacobi
  LU factors** and **Jacobi smoother diagonals** after hierarchy
  construction and factorization, retaining full precision per array
  when any value is unsafe to encode.
* `noprobe`: skip the per-coefficient 1% relative tolerance guard, but
  still reject infinities, NaNs and nonzero-to-zero underflow.

Conversion failures and allocation failures are nonfatal: the corresponding
array stays full-precision. Compression is disabled when the build reports
FP16 unavailable, including complex builds.

Resident prolongation coefficients, block-smoother LU factors and
scalar-/block-Jacobi smoother diagonals can be compressed. Coarse CSR
matrices, pivots, terminal dense factorization and work vectors remain in
LIS scalar precision.
Restriction, prolongation and triangular solves decode coefficients and
accumulate in the solver precision. This is **not yet full mixed-precision
storage for the entire multigrid hierarchy**. Like ILU's half mode, half
storage is opt-in, separate from solver arithmetic precision.

The half-precision update path rebuilds the hierarchy transactionally
instead of attempting numeric in-place refresh of released arrays.

## Numerical range and fallback behavior

IEEE binary16 has a largest finite value of **65504**, a smallest positive
normal value of **2^-14**, and a smallest positive subnormal of **2^-24**.
MFAMG does not apply automatic per-level scaling to its stored coefficients.
For example, unscaled stiffness matrices can produce diagonals or block
factors larger than 65504. Those arrays **remain in the original precision**
instead of being silently saturated or stored as infinity. Similarly, any
nonzero coefficient that would round to zero is kept in full precision.

The `on` mode applies an **element-wise** error check (1% relative tolerance
plus a small absolute allowance for transfer/block factors, and a reciprocal
error check for Jacobi diagonals). It is **not** an end-to-end convergence
probe. The `noprobe` mode removes the rounding-error tolerance entirely;
it still rejects nonfinite values and nonzero-to-zero underflow, but may
admit more numerical error. Prefer `on` for production uses and validate
convergence against `off` on your own matrices. These settings reduce
storage of eligible arrays only when their values are representable.

## Testing

`testmfamg` now tests `off`, `on`, and `noprobe` for weighted-Jacobi
and block-Jacobi smoothers. The regressions inspect resident packed arrays
and compare V-cycle correction to the full-precision reference. The GitHub
Actions Linux `half` job configures with `--enable-half`; other jobs
exercise compatibility with default builds.

Further work: full-precision residual certification under harder matrix
condition numbers, robust fallback probing, and total resident-memory
benchmarks for large problems.

## Optional 3D Poisson benchmark (not part of default make check)

The `testmfamg` executable has an opt-in seven-point 3D Dirichlet
Laplacian benchmark. Build with `--enable-half`, then run each mode in
a **fresh process**:

```sh
./test/testmfamg --fp16-3d-bench 48 off
./test/testmfamg --fp16-3d-bench 48 on
```

The command reports iteration count, true relative residual,
eligible numeric coefficient-array MiB, and (Linux only) process-wide
`VmRSS` before/after preconditioner construction. It uses GMRES with
a deterministic manufactured right-hand side and a constant near-nullspace
vector. The `half` CI job runs 12^3 and 48^3 cases separately.

### Reproducible GitHub Actions result (48^3 = 110,592 DOF)

| Measurement | FP64 storage (off) | FP16 storage (on) |
| --- | ---: | ---: |
| Eligible coefficient arrays | 3.374802 MiB | 0.843700 MiB |
| GMRES iterations | 21 | 21 |
| True relative residual | 8.25549e-9 | 8.37126e-9 |
| Whole-process VmRSS, after preconditioner setup | 45,860 KiB | 45,792 KiB |

This demonstrates **2.531102 MiB saved (75%) in eligible resident
coefficient arrays**, without an increased Krylov iteration count on
this particular problem. The measured whole-process VmRSS difference
was only **68 KiB**, *not* 2.53 MiB. It reflects the matrix, input/output
vectors, library/workspace allocations, allocator page retention, and
process/address-space effects. These results do **not** establish
whole-preconditioner or peak-RSS savings, and are not benchmarks of
real-world FEA stiffness matrices. More problem types and independent
memory profiling are required for general performance claims.

## Diagonal safety and memory accounting

After the Jacobi diagonal has been computed and validated in solver
precision, it may be packed in FP16. The `on` policy guards the relative
error of its **reciprocal** (as used in the smoother); `noprobe` skips
that tolerance but still forbids zero/NaN/Inf. If a level's diagonal
cannot be packed, it remains at the original precision. This applies
to both weighted Jacobi and block-Jacobi's scalar fallback.

A read-only internal diagnostic, 
`lis_precon_mfamg_get_half_coefficient_bytes`, reports the full-
precision equivalent and current resident bytes of **only three
eligible coefficient categories** (prolongation, block LU, diagonals).
The test suite compares both counters, confirms no compression in
`off`/unsupported configurations, and verifies the PSD rebuild path.
These counters exclude matrix storage, work vectors, indices,
allocator metadata and overhead: they must **not** be interpreted as
whole-preconditioner resident memory or runtime measurements.
