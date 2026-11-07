# Generic FEA operator snapshots

This is a **read-only, solver-independent** audit of assembled sparse linear
systems. `tools/fea_matrix_audit.py` works with output from *any* FEA program
(or other numerical application) that follows the format below. Contact
mechanics is one use case, not a hard-coded assumption.

## Format

- `A.mtx`: the **square assembled linear operator actually passed to the
  solver**, in Matrix Market coordinate real/integer/pattern general or
  symmetric format. Real/integer general array format is also supported.
  Symmetric coordinate input is expanded; duplicate coordinates are summed.
  Complex, Hermitian and skew-symmetric files are not supported.
- Optional `b.mtx`, `x0.mtx`, `x.mtx`: RHS, initial iterate, and computed
  solution as n-by-1 Matrix Market matrices.
- Optional `B.mtx`: n-by-k candidate coarse/near-nullspace vectors, one
  vector per column, for elasticity, diffusion, mixed formulations, etc.
  Matrix Market arrays are **column-major**.
- Optional `active_dofs.csv`: active equation mapping. Required columns are
  `eq,node,dof`; optional `body` identifies a part, domain or substructure.
  `eq` is zero-based, contiguous, and occurs **exactly once** across `0..n-1`.
  `node` and `dof` are arbitrary nonempty labels; for example numerical
  IDs, `UX`, `ROT_Z`, `TEMP` or `PRESSURE`. No fixed number or ordering
  of DOFs per node is imposed. `body` can be blank if the exporter does not
  supply such groups; if present, IDs should remain stable across snapshots.

Example CSV (works for partially constrained, multiphysics or shell models):

```csv
eq,node,dof,body
0,42,UX,part-A
1,42,ROT_Z,part-A
2,101,TEMP,
```

An exporter must map its native equation numbering (possibly one-based)
to this *zero-based active equation order* and apply all relevant boundary
conditions/elimination. Do not infer a part identifier solely from a matrix
graph after introducing contact couplings.

## Usage

```sh
python3 tools/fea_matrix_audit.py \
  --matrix A.mtx --rhs b.mtx --initial x0.mtx --solution x.mtx \
  --modes B.mtx --dof-map active_dofs.csv > audit.json
```

`--expected-dofs-per-node N` is an **optional, explicitly provided physical
assumption** (such as `3` for displacement-only 3D solids, `6` for some
shell formulations, `1` for thermal problems). Without it the tool reports
a histogram of active DOFs per node and does **not** infer any missing DOFs.
With it the tool also reports how many nodes have fewer than N active DOFs.
This does not distinguish constraints from unused fields.

## Output

The JSON includes SHA-256 of A and b, operator dimension, exact numeric NNZ
(after duplicate cancellation), expanded stored coordinates, relative
Frobenius asymmetry, infinity norm, diagonal statistics, optional true
residuals `||b-A*x0||_2/||b||_2` and `||b-A*x||_2/||b||_2`,
per-mode Rayleigh quotients and normalized operator actions, and
optional equation/group/DOF-map histograms. If b is zero, residual
normalization uses denominator 1 (absolute residual).

Invalid input returns exit code **2** and a concise error, including malformed
CSV rows and nonfinite Matrix Market data. The tool is written using only
the Python standard library and stores an expanded sparse matrix in Python
dictionaries; large-scale performance measurements should instead use an
optimized sparse matrix package.

## Reproducible comparisons

Freeze the **same** A, b, x0, B and equation map for every backend. Exporter
metadata should record the originating FEA program/version, mesh, analysis
type, solver tolerances, precision, load/step/Newton iteration and any contact
active-set statistics where applicable. True residuals and convergence
histories must be checked before comparing speed. No condition number,
singular value, or nonlinear contact convergence claim is made here.

A legacy entrypoint `tools/contact_matrix_audit.py` delegates to this generic
tool for existing scripts. The library API and the CLI are FEA-independent.

Run regression tests with:

```sh
python3 -m unittest discover -s tools -p 'test_fea_matrix_audit.py' -v
```
