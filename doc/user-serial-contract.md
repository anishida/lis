# Serial USER operator contract

`LIS_MATRIX_USER` represents an application-owned operator. It does not
own a sparse storage representation or the application callback context.
This guide describes the single-process path; it does not establish support
for distributed row-backed operators.

## Initialization and ownership

Create a matrix, set its size, then initialize USER once with
`lis_matrix_set_user(A,context,matvec,matvech)`. The descriptor must be
unassembled at initialization. Either install the needed explicit operator
callbacks or install a row callback with `lis_matrix_set_user_row`.
A row-only matrix can supply serial forward, Hermitian-adjoint and diagonal
operations. Explicit callbacks have precedence.

Keep the application context and all referenced data alive until every
matrix operation, solver, preconditioner or session using it has completed.
The USER descriptor borrows the context; it is not responsible for freeing
application storage. Use a new matrix descriptor to replace the context
or the primary callbacks. Do not reinitialize an assembled descriptor
with `lis_matrix_set_user`.

Single-row queries use local row indices and return zero-based global column
IDs. In serial these coincide with local column IDs. For a count query,
both output arrays are NULL; the fill query must reproduce the number of
entries and respect capacity. Duplicate and unsorted columns are allowed.
The adjoint means conjugate transpose for complex scalars, not ordinary
transpose. See [batched rows](user-row-batch.md) for the optional batch
interface and its fallback rules.

## Mutation and reuse

Do not mutate callback state, replace callbacks or destroy a matrix during
an active operation. After a value-only change, call
`lis_matrix_mark_values_changed(A)`. After any pattern change, call
`lis_matrix_mark_structure_changed(A)` (which also changes the value
generation). Mark each changed operand of a composed operator. Recreate
long-lived row sources after callback or data changes.

Row callbacks are not assumed thread-safe. The explicit thread-safe
opt-in applies to the currently installed single-row and batch callbacks
and their data. Replacing either requires reaffirming the opt-in.
Do not share mutable scratch across concurrent calls. The row-backed forward
kernel can use OpenMP; the row-backed Hermitian-adjoint kernel currently
accumulates contributions sequentially.

The optional `lis_matrix_auto_session_*` interface in
`lis_matrix_materialize.h` owns reusable solver and materialization state.
Configure its borrowed solver, keep the source alive, mark generations
and reset/destroy the session before destroying the source. Never destroy
a borrowed cache-owned matrix or borrowed session solver yourself.

## Supported configuration and explicit boundaries

| Operation | Direct serial USER path |
| --- | --- |
| Standard-precision Krylov solve | Requires forward operator |
| BiCG/BiCR and other adjoint-dependent methods | Requires Hermitian-adjoint capability |
| Jacobi preconditioner | Requires diagonal capability |
| Direct ILU-family row consumers | Require row access and compatible solver/preconditioner |
| Symmetric diagonal scaling | Requires diagonal capability and preconditioner none |
| QUAD / SWITCH solve | Rejected by the matrix-free solver capability check |
| Jacobi / Gauss-Seidel / SOR as standalone solvers | Require explicit matrix storage |
| Explicit format conversion through solver storage option or -use_at | Not supported for direct USER; use the documented row-backed materialization API where applicable |
| In-place diagonal or matrix shift | Not supported for USER |
| Generalized eigensolver paths needing a shift | Only the explicitly supported composed-operator paths |

Existing regressions include `testuserrow`, `testjacobiuser`,
`testshellscale`, `testoperatorrowprecon` and `testsolvercapabilities`.

This table is a capability guide, not proof that every numerical solver
or error path is covered. Test the selected combination, including the
application callbacks and preconditioner. The default-precision scalar
callback ABI cannot represent the low component of QUAD values, so
supporting true QUAD/SWITCH requires a separately specified contract,
not silent conversion to double.

## Errors and preconditioner names

Return a LIS error code when a callback fails. Its output is not valid on
failure and must not be consumed. Check the public API return and the
solver/eigensolver status separately; a successful API execution is not
necessarily a converged numerical solve. On failure do not treat the
returned solution or residual as certified. Explicitly verify the physical
residual for a successful solve.

Registered user preconditioner names must fit `LIS_PRECONNAME_MAX` bytes,
excluding the terminating NUL. Longer names are truncated by registration
and may not match the name given to the solver option parser. Use short,
distinct names, keep their callbacks alive, and do not unregister them
while a solver or preconditioner still uses them.
