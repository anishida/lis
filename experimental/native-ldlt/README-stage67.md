# LIS BUCKLE stage 67 — OpenMP level-scheduled multifrontal LDLᵀ

Standalone C11 prototype, candidate patch 56. Extends stage 66 without modifying its serial numeric path. No BLAS, OpenBLAS, LAPACK, METIS or SuiteSparse. The optional parallel entry points link OpenMP (GCC libgomp in this environment), libm and libc.

## API and algorithm

- `p67_factorize`: natural order with automatic symbolic analysis and strict supernode grouping.
- `p67_ordered_factorize`: Natural/MD/ND/AUTO, combines symbolic and numeric permutations.
- `p67_numeric`: lower CSC and caller-provided valid elimination tree and optional contiguous supernode boundaries; requests 1–64 threads. One thread uses the original numeric implementation. A contracted tree with only one ready front per level falls back to the serial implementation after scheduling analysis. `levels=0` means scheduling skipped in the explicit one-thread path, not a tree with no levels.
- Stable level scheduling: independent fronts at the same level are distributed dynamically, with a barrier between levels. Each thread has its own ID map and candidate array. Each front owns its contribution; parents consume their children's contributions only after the prior levels finish. Pivot records are keyed by original ID and gathered in original topological front order before packing. This preserves the serial front ordering even when execution order differs.
- Global workspace allocation, factor-slot and delay limits protected by OpenMP critical sections. A failed front prevents later levels from running, leaves the old context and output statistics unchanged and releases partial results. Failure code may depend on which worker encounters a resource limit first. SIMD kernel and floating-point environment are inherited from the caller for the parallel region, and the worker environment is restored.
- No nested parallelism inside Schur kernels, no parallel triangular solve. Immutable factors still support concurrent independent solves with private scratch/output.

## Reproduction

`make parallel` (GCC, `-fopenmp`); `make parallel-sanitize` (GCC ASan+UBSan); `make -j4 stage67-check` (prior serial GCC tests and Clang sanitizers, new GCC OpenMP tests, separate portable SIMD test). `OMP_NUM_THREADS=4 OMP_PROC_BIND=true make parallel-bench` gives wall-clock medians of three independent runs, including ordering, symbolic analysis, numeric factorization and packing but excluding solution verification. The benchmark verifies all solutions and bitwise agreement with one-thread results. No CCX data available.

Observed: 535/535 new tests with and without ASan+UBSan; 8794/8794 cumulative normal/sanitized tests (plus 162 portable SIMD tests, counted separately). 360 randomly generated primary cases really ran with more than one worker; 1581 2x2 pivots and 1581 delays in the dedicated paired cases. Max solution error 1.332e-15, max componentwise backward error 6.393e-16. Older test suites exercise the unchanged serial API; they are not all rerun through the parallel scheduler. See `test-results.txt` and `benchmark-parallel.txt`.

## Limitations and open validation

The scheduler allocates O(n * threads) index maps and O(n) front metadata, and multiple simultaneous updates may increase memory. The numeric budget counts requested payloads of this module (including schedule and per-thread maps), excluding allocator metadata, input CSC, symbolic/reordering phases, old factor, final adapter allocations, and OpenMP runtime/stacks. On a single-chain tree the parallel planning overhead is currently paid before falling back, so peak workspace may increase even though only one thread is used. Barriers limit scaling on deep trees, and the global allocation/limit lock may be a bottleneck. No guarantee of reproducibility across different CPUs or floating-point modes.

Additional GCC ThreadSanitizer + libgomp attempt emitted race warnings, including reads/writes within the same OpenMP critical section and across OpenMP synchronization boundaries, and was terminated after 12 s (exit 124). These warnings may arise from insufficient recognition of OpenMP synchronization by this tool/runtime combination; they are NOT treated as a successful race check or dismissed as proven false positives. Raw log preserved in `thread-sanitizer.txt`. No compatible libomp/Archer detected in this environment. Further race-aware validation and review are required before production use. ASan/UBSan do not certify absence of races.

No integration with original LIS/CCX, no real BUCKLE benchmark, no MPI or ABI validation, no GitHub publication. Stage 68: reduce memory overhead, improve cache and allocations, and finish concurrency validation. Source files and logs from older stages included as historical reference.
