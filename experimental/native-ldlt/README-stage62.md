# LIS BUCKLE — stage 62, candidate patch 51 (standalone prototype)

Stage 62 adds automatic symbolic analysis of a symmetric lower CSC matrix to
stage 61's multifrontal numeric factorization. C11, standard C library and libm;
no BLAS, LAPACK, OpenBLAS, METIS or SuiteSparse dependencies. This is NOT a
patch applied to the original LIS sources, nor a full CCX BUCKLE benchmark.

## API and algorithm

`symbolic_ldlt.h` exposes `p62_analyze`, `p62_free` and `p62_factorize`.
The input CSC must have n+1 column pointers, col[0]=0, col[n]=nnz, and strictly
increasing lower-triangular row indices in each column. `nnz` is the actual
row/value array length. A missing diagonal is permitted; explicit zeros are
structural entries. The numeric solver rejects non-finite values.

1. Validate the CSC, transpose the strict lower pattern to lists of previous
   neighbors per row. No global dense allocation.
2. Compute the elimination forest using ancestor path compression. Roots
   are UINT32_MAX; non-root parents have a larger index.
3. For each row, traverse the elimination-tree reach, marking visited columns
   once per row. Compute exact symbolic column counts (including the unit
   diagonal), the total strict lower structural slots and the maximum symbolic
   column front for natural-order scalar elimination. The counts ignore numerical
   cancellation and do not model delayed or block pivots.
4. Detect maximal contiguous strict supernodes: adjacent j,j+1 qualify when
   parent[j]==j+1 and count[j]==count[j+1]+1. The filled graph implies that
   the tail of column j after j+1 is a subset of column j+1's strict lower
   pattern. Equality of counts gives equality of the tails. Boundaries are
   in super[0..nsuper], with super[0]=0, super[nsuper]=n.
5. `p62_factorize` invokes stage 61 with the derived forest, keeping its
   Bunch–Kaufman strategy, delayed pivots, packing and permutation handling.

Symbolic requested memory is O(n + nnz_offdiag), independent of symbolic fill
size. The reach traversal costs O(n + nnz + symbolic slots) after computing
the tree, aside from the path-compression cost. No claim of a linear-time
column-count algorithm independent of fill is made.

## Example

```c
p62_limits sl = {128u * 1024 * 1024, 100000000};
p61_limits nl = {256, 128u * 1024 * 1024, 10000000, 1000000, 0};
p56_context *ctx = NULL;
p61_stats stats = {0};
int rc = p62_factorize(&ctx, n, nnz, col, row, values, &sl, &nl, &stats);
/* On success, call p56_apply with a separate scratch vector. */
p56_destroy(ctx);
```

`p62_analyze` permits a NULL initialized output or a previously owned result.
On failure it leaves the previous result untouched. The wrapper leaves the
previous numerical context and stats unchanged on failure. The caller must
serialize replacement and destruction against concurrent applies.

`max_work_bytes` limits requested symbolic allocation sizes (excludes allocator
metadata, caller arrays, existing result and numeric allocations). `max_slots`
limits symbolic strict lower positions; the separate stage-61 numerical limits
are required to bound delayed pivots and actual factor slots. Neither budget is
a global RSS cap. A symbolic result can be retained for repeated stage-61
factorizations of exactly the same structure; this wrapper recomputes it.

## Reproduce tests

```
make stage62-check
make clean
```

GCC 14.2 optimized and Clang 17 ASan+UBSan: 4305/4305 tests passed each.
897 new cases: 700 random structures (240 also solved numerically), 160
nonadjacent indefinite paired cases (with/without explicit zero diagonals),
34 error/limit/transaction checks, and 3 large cases. All 860 small generated
structures compared with an independent explicit graph clique-elimination
oracle for parents, counts, roots, slots, max front and strict supernode
boundaries. Numerical solves verified against independent symmetric CSC
matvec and known solution. Scales 1e-100 through 1e100 tested.
Max absolute solution error 1.554e-15; max componentwise backward error
7.619e-16 in new tests. See test-results.txt for complete actual logs.

Optimized final run (CPU times, one run):

| case | lower CSC entries | symbolic L slots | strict supernodes | max front | symbolic bytes | symbolic CPU s | auto factor CPU s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| tridiagonal 152592 | 305183 | 152591 | 152591 | 2 | 5493392 | 0.003303 | 0.114075 |
| grid 40x40 | 4720 | 62439 | 1560 | 41 | 63764 | 0.000158 | 0.559590 |
| grid 128x128 | 48896 | 2080895 | 16256 | 129 | 654420 | 0.003829 | not run |

The last test is symbolic only. For the first two, numeric factor slots and
peak front matched predictions (diagonally dominant synthetic input). Numeric
requested peak bytes: 22737252 and 1707248 respectively. These separate peaks
are not process RSS, nor a measured combined peak. Times do not represent
CCX, BUCKLE, MPI or repeated-run performance statistics.

## Remaining work

- Numerical amalgamation of strict supernodes: detection only in stage 62.
- Fill-reducing ordering (e.g. own minimum-degree or nested dissection),
  elimination-tree postordering and memory-aware scheduling. Natural ordering
  remains in effect, so difficult graphs can create large fill/fronts.
- Global allocation budgeting, symbolic cache and structural fingerprint.
- Integration with actual LIS source and CCX problem; full MPI and ABI tests.
  No original stage 52/54 source or actual CCX matrix was supplied. Stage 52
  remains the previously recommended production option.
