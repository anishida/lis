# Optional batched row access for `LIS_MATRIX_USER`

This experimental, opt-in extension reduces per-row function-pointer dispatch in
USER matrix-vector operations and row-source consumers such as matrix
materialization and sparse preconditioner setup. It does **not** alter
ordinary CSR or the existing `LIS_MATRIX_USER_GET_ROW` contract.

## API

First initialize the USER matrix with `lis_matrix_set_user()` and install
the existing single-row callback:

```c
lis_matrix_set_user(A, user_data, NULL, NULL);
lis_matrix_set_user_row(A, get_one_row);
lis_matrix_set_user_rows_batch(A, get_multiple_rows);
lis_matrix_set_user_row_thread_safe(A, LIS_TRUE);
```

The batch callback signature is:

```c
LIS_INT get_multiple_rows(
    void *user_data, LIS_INT row_begin, LIS_INT row_end,
    LIS_INT *offsets, LIS_INT capacity,
    LIS_INT *cols, LIS_SCALAR *vals);
```

The row interval is half-open `[row_begin,row_end)`, using locally owned
row indices, while returned columns are zero-based **global** column IDs.
`offsets` contains `row_end-row_begin+1` entries, starts at zero, and
must be nondecreasing. The final offset is the number of nonzeros.

For the first, *count* call, `cols` and `vals` are `NULL` and
`capacity=0`. The callback must fill `offsets`. For the subsequent
*fill* call, `capacity` is at least the required number of nonzeros and the
callback must reproduce the same offsets while filling the arrays.

The matrix must remain unchanged between the two calls. LIS checks for
inconsistent offsets and out-of-range column indices. The application must
not retain any of the scratch pointers beyond the call.

**Mathematical equivalence is mandatory:** both installed callbacks must
describe the same rows, column/value pairs, duplicate-entry semantics and
matrix dimensions (the order of entries inside a row may differ). Depending
on locality and available memory, LIS can use scalar and batch callbacks
within the *same* matrix-vector operation. Different numerical meanings
would produce an invalid, state-dependent linear operator. A callback
replacement is a change of retrieval strategy, not a method for silently
updating matrix coefficients. If the physical operator changes between
operations, mark the matrix's value or structure generation via
`lis_matrix_mark_values_changed(A)` or
`lis_matrix_mark_structure_changed(A)` as appropriate, and recreate
long-lived row sources / materialization state before resuming.

## Performance and safety

- Each row-source instance caches up to **16 adjacent rows**. Sixteen
  is a conservative initial upper bound: the 17-entry offset array stays
  inline per worker and a sequential traversal amortizes two callback calls
  over up to 16 rows. It is **not universally optimal**.
- **Adaptive locality policy:** after two nonconsecutive misses of the
  cached window, that worker temporarily uses the ordinary scalar callback.
  After eight consecutive requests it resumes batch mode. The adaptation is
  per row-source, deterministic given the row order, and requires no public
  configuration or additional matrix storage. It avoids repeatedly fetching
  16-row windows for random-looking access while preserving batching in
  ordinary SpMV. A single random jump does not permanently disable batching.
- The locality heuristic is independent of allocation-failure fallback and
  does not suppress failures returned by the application callback in a
  batch operation that is actually attempted. On the scalar route,
  errors from the scalar callback are propagated as usual.
- If allocating the batch index/value scratch fails, the worker disables
  batching for the lifetime of that row-source object and retries through
  the original scalar callback. Insufficient memory for even the scalar
  fallback is still reported as an error.
- Errors from the application batch callback, changing offsets between
  count/fill, and out-of-range column indices are **not** silently downgraded
  to scalar callbacks. They invalidate the current batch and return an error to
  the direct row-source caller. The *solver* must also check error codes:
  some older Krylov paths on baseline `develop` ignored USER matvec errors.
  The independent correctness PR #98 fixes selected GCRO-DR paths; deploying
  #97 without those changes does **not** guarantee end-to-end solver error
  propagation for every Krylov algorithm. A new request may retry
  a failed batch after the application corrects its state.
- Nonnegative column IDs below the local dimension are required in serial
  mode. Duplicate and **unsorted** columns remain permitted, as in the
  existing scalar USER row interface; no sorting is performed.
- USER matvec and matrix materialization create worker-private row-source
  state; cache buffers and offsets are never shared across workers.
- Only serial/local USER column numbering is supported by the row-source
  implementation (same as the preexisting USER direct-row path).
- A new batch callback resets the thread-safe flag. Call
  `lis_matrix_set_user_row_thread_safe(A,LIS_TRUE)` only if **both** installed
  callbacks and their application state are safe for concurrent calls.
- Explicit application-provided matvec still takes precedence over the
  generic callback-backed matvec.
- The single-row callback remains mandatory in this first version, for
  existing direct-row API compatibility and automatic fallback.
- Do not mutate matrix coefficients or callbacks while another thread is
  constructing or using a batch cache.

## Compatibility

**No fields have been added to `LIS_MATRIX_STRUCT`.** The public
structure remains byte-for-byte unchanged relative to the upstream develop
baseline. The new exported setter is additive and stores the optional
callback in a private, globally synchronized registry keyed by the matrix
pointer. Destruction of the matrix unregisters the callback; replacing
`lis_matrix_set_user()` or the one-row callback also unregisters it, and
calling the new setter with `NULL` disables it.

There is one registry lookup per worker-local row-source **initialization**,
not per row. Worker-local row sources snapshot the registered function.
For that reason, no callback re-registration, matrix destruction, or matrix
coefficient mutation may overlap an active matvec, preconditioner operation,
or row-source object's lifetime. Destroy/recreate a row source after
replacing callbacks. Creating/destroying separate matrices can use the
registry concurrently in OpenMP builds.

Keeping the original public layout avoids changing `sizeof(LIS_MATRIX_STRUCT)`
and its field offsets. It does not eliminate the need to rebuild LIS itself
and use an ABI-compatible library release/version policy when deploying
the new exported setter symbol.

## Validation

`test/testuserbatch.c` compares scalar and batched matvec/matvech, checks
out-of-order row requests, validates malformed offset rejection, and counts
callback dispatches. The separate `user-batch-regression.yml` workflow
runs one and eight OpenMP threads with normal and ASan/UBSan builds.

The benchmark `test/testuserbatchbench.c` prints CPU time and callback
counts for four row-source cases: lightweight versus deliberately expensive
callback dispatch, each in sequential versus strided row-access order.
It additionally benchmarks **the actual generic `lis_matvec` USER path**
for both callback types. Scalar and batched paths must return identical
numerical checksums.
Timing results are not CI pass/fail thresholds: small hosted runners have
no stable benchmark isolation. A large callback-count reduction is not
evidence for equal gains in SpMV or the full FEA solver.

For release acceptance, independently measure whole SpMV and setup+solve
performance on real application matrices with two representative callback
implementations, plus a random-row preconditioner workload. Record total
time, true residuals, iterations, allocations and peak RSS against upstream
develop.
