# Export generic LIS linear-system snapshots (opt-in)

`lis_snapshot_export_mm()` is a public helper in
`#include "lis_snapshot.h"` that writes **the matrix actually passed to LIS**
and the application's current vectors as standard Matrix Market files.
This generic interface can be reused by any application that uses LIS,
independently of its physical model, mesh, element or field representation.

```c
#include "lis_snapshot.h"
/* Call after assembly, before the solver overwrites an input vector. */
int rc = lis_snapshot_export_mm(A, b, x0, NULL,
                                near_nullspace_dim, near_nullspace,
                                "snapshots/step03_iter07");
if (rc != LIS_SUCCESS) { /* discard that prefix, report the error */ }
```

The path prefix must refer to an already existing directory. Every
successful call writes `<prefix>_A.mtx`; optional non-null `b`, `x0`, `x`
and `basis_count > 0` cause `_b.mtx`, `_x0.mtx`, `_x.mtx` and
`_modes.mtx` to be written. The modes form a column-major n-by-k array.
The distinct suffix prevents `_b`/`_B` collisions on case-insensitive filesystems.

## Scope and safety

- Explicit opt-in, no environment-triggered writes and **no changes** to
  `lis_solve` or the iterative numerical path.
- Real-valued, **single-process** matrices only in the first version.
  Compatible with native CSR, convertible explicit formats, and row-backed
  `LIS_MATRIX_USER` / `LIS_MATRIX_OPERATOR`. Operator-only matrices with
  no accessible rows, complex builds, distributed MPI and split CSR
  are unsupported and must fail instead of exporting a misleading operator.
- The input matrix and vectors are borrowed and never modified. User row
  callbacks and external CSR storage must not be mutated during capture.
- Numeric values must be finite; data are written at 17 significant digits
  for double-precision interchange. Applications requiring exact extended
  precision should use a suitable lossless format instead.
- Existing destination files (including staging filenames ending in
  `.tmp`) are rejected before capture. The exporter does not deliberately
  overwrite existing snapshots. **Use a unique prefix per call**.
- Each file is written to a temporary path, closed and renamed. If a later
  file fails, previously published files from the same call are removed
  (best effort). Multiple files are **not committed atomically**, and the
  preflight check is **not a cross-process lock**; concurrent writers must
  never share a prefix.
- Physical metadata is **optional**. When requested,
  `lis_snapshot_export_dof_map_csv()` writes
  `<prefix>_dofs.csv` using zero-based `eq,node,dof,body` columns.
  No fixed DOF count or ordering is assumed.

This captures **already assembled and eliminated** equations. A generic LIS
exporter cannot reconstruct full FE coordinates, inactive constraints or
body labels from algebraic rows. Native application adapters remain
responsible for such metadata.

## Auditing and existing LIS writers

LIS already provides `lis_output_matrix()`, `lis_output_vector()` and the
internal `lis_output_mm_*` functions for standalone Matrix Market output.
This advanced snapshot API builds an independent capture contract on top
of LIS matrices: validated serial CSR conversion (including row-backed
USER/OPERATOR matrices), optional input/solution/mode arrays, consistent
filename conventions, 17-significant-digit scalar output, and staged writes
with best-effort rollback. The existing generic writer has different
interfaces and lifecycle/error semantics; unifying the low-level formatting
paths can be considered separately.

Offline numerical auditing is **not bundled in this pull request**.
A separate [generic matrix audit utility](https://github.com/anishida/lis/pull/205)
can consume the emitted standard Matrix Market files once installed.
For example, when the auditor is present:

```sh
python3 tools/fea_matrix_audit.py \
    --matrix snapshots/step03_iter07_A.mtx \
    --rhs snapshots/step03_iter07_b.mtx \
    --initial snapshots/step03_iter07_x0.mtx \
    --modes snapshots/step03_iter07_modes.mtx \
    --dof-map snapshots/step03_iter07_dofs.csv
```

The optional auditor can calculate true residuals and operator/nullspace
metrics from identical algebraic data; this snapshot exporter makes no
convergence or performance claims.

## Optional, application-defined equation labels

The native API `lis_snapshot_export_dof_map_csv()` writes the optional
`<prefix>_dofs.csv` file accepted by the generic FEA auditor. The application
provides labels for each **active equation** with a callback; LIS does not
guess node ordering, constrain equations, or read a mesh:

```c
LIS_INT label_dof(void *ctx, LIS_INT eq, const char **node,
                  const char **dof, const char **body)
{
    /* Application uses its own eq -> mesh/entity/field map. */
    /* Assign *node, *dof and *body; body may be NULL. */
    return LIS_SUCCESS;
}

LIS_INT rc = lis_snapshot_export_dof_map_csv(
    A, label_dof, my_application_context, "snapshots/step03_iter07");
```

The callback must return nonempty node/entity and DOF/field labels; labels
may contain commas, quotes and whitespace (written as correctly quoted CSV).
`body` is optional, not a hard-coded structural partition. The helper calls
the callback once for every active serial equation in the same zero-based
ordering used by the exported operator, and writes
`eq,node,dof,body`. It is independent of Matrix Market export and works
also in complex LIS builds because labels are not numerical values.

Each DOF CSV is staged as `<prefix>_dofs.csv.tmp` and published after
successful close. An existing destination is rejected, and incomplete staged
CSV is deleted on failure. This API is separately invoked from the numeric
exporter, so the two calls are not one atomic transaction.
The caller owns its label memory and callback state throughout the call.
Use a unique prefix for each snapshot. No fixed field list, physical
dimension, number of DOFs per node or application package is required.
