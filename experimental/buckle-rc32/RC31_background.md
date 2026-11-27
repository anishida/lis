# LIS BUCKLE — consolidated experimental release candidate (stages 1–31)

## Status
Local prototype against the exact uploaded `lis-2.1.25.zip`, not a port to the current GitHub develop branch. No GitHub modifications. The archive contains code/patches, test drivers and recorded results; no private CCX matrices or mode/RHS vectors. The independent reference eigenvalue table is included for certification. The full consolidated patch changes 12 files including two new internal includes and distribution lists. Apply it to pristine source, not to a tree already containing incremental patches. Rebuild library and callers: the included stage11 inner near-nullspace API changes the public LIS_ESOLVER_STRUCT layout (ABI). The functional patches 01–03 are supplied separately for potential review, not to be applied in addition to the consolidated patch. None is claimed merged or production-ready.

## Clean build (Linux, GCC/OpenMP, real double, 32-bit indices, static library)

```sh
unzip lis-2.1.25.zip
cd lis-2.1.25
patch -p1 < /path/to/LIS_2.1.25_buckle_consolidated_experimental.patch
# TARGET is an unrelated environment variable in the benchmark runtime.
env -u TARGET ./configure --enable-omp --disable-shared
make -j4
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 make -j4 check
```
If using a ZIP extraction method that loses executable modes, restore them from the archive or `chmod +x configure config/{compile,depcomp,install-sh,missing,mkinstalldirs,test-driver} test/test*.sh`. Initial configure failure due to TARGET and initial script permission failures were resolved. Final clean standard test suite rerun from scratch: 126 listed tests (125 PASS, one SKIP: OOM test requires --enable-oom-tests), zero FAIL, exit 0. The test runner lists some entries twice; counts are its reported counts, not necessarily unique test programs. Experimental environment switches were unset during standard regression. Logs: `check_uniform.log`, `test-suite.log`.

## What is included
- 01–03: Block GLI NEV certification, Block GLI linear options, selective CSR materialization.
- 04: owned inner near-nullspace prototype from stage11 (ABI change; separate coarse-basis contract on newer develop).
- 09–12: symbolic LU, paired/grouped up to eight row updates, serial triangular apply with OpenMP for other work.
- 14–16: optional level-scheduled numeric LU (disabled in recommended config), DOT4 triangular solve, integrated lazy exact minimum degree and vector permutation.
- 17–19: exact-pattern reuse, direct scatter from original CSR, integration into sequence refresh and eigensolver preconditioner rebuild. Incremental patches for these three are in `stages29_31/` for inspection.
- Not included: stage12 continuation prototype 05, singleton absorption, unsuccessful incomplete ILUTP global-drop experiment, group16. Other investigated alternatives (spectral MFAMG, smoothed aggregation, etc.) are not implemented in this package. No supernodal BLAS-3, LDL^T or true scalable parallel numeric factorization.

## Recommended experimental benchmark switches
Flags are environment-presence based: unset to disable (setting `=0` does not disable). Use them consistently throughout a session. Do not globally enable them for unrelated preconditioners or change them mid-solve.

```sh
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
export LIS_EXPERIMENT_ILUTP_SYMBOLIC=1 LIS_EXPERIMENT_ILUTP_PAIRS=1
export LIS_EXPERIMENT_ILUTP_GROUPS=1 LIS_EXPERIMENT_ILUTP_GROUP8=1
export LIS_EXPERIMENT_ILUTP_SERIAL_APPLY=1 LIS_EXPERIMENT_ILUTP_DOT4=1
export LIS_EXPERIMENT_ILUTP_MD=1 LIS_EXPERIMENT_ILUTP_REUSE=1
unset LIS_EXPERIMENT_ILUTP_PARALLEL_NUMERIC LIS_EXPERIMENT_ILUTP_RUNS
```
Linear options: `-i fgmres -p ilut -ilutp 1 -ilutp_threshold 0 -ilutp_full 1 -ilutp_full_limit 45000000 -tol 1e-10 -maxiter 40 -restart 35 -print none`.
GLI options: `-e gli -etransform buckling -shift .00107 -ss 52 -nev 6 -etol 1e-8 -emaxiter 100 -omp_num_threads 4`.
Full ILUTP is exact LU without dropping or pivoting in this path, not an incomplete preconditioner; avoid unsupported/unstable matrices. The input must be sorted, unique CSR with symmetric structure, serial real-double 32-bit. The full-fill limit counts factor off-diagonals, not peak bytes. Pivoting, distributed matrices, complex, arbitrary scaling and half-storage are not supported for reuse. For matrices outside the contract, use existing paths without the experimental flags.

## Cache ownership and limitations
`lis_solve_sequence()` uses `lis_solver_set_sequence_precon_refresh_policy(solver,LIS_SEQUENCE_PRECON_REFRESH_FROZEN_PATTERN)` and the caller must mark value/structure changes with the existing generation API. Reuse compares exact CSR pattern, checks configuration and fill limit, and refactors on stored L/U pattern; old factors are destroyed on failure in the integrated callers. A pattern mismatch falls back to a fresh build. Cache is owned by ILUT preconditioner via reserved user_data + cleanup callback; public get/set user_data rejects access while cache is attached. Explicit experimental session is also supplied for testing; borrowed factors must not outlive session, and caller must detach before destruction. Cache only reuses within the same preconditioner lifetime, not across independent nested GII solver objects. GLI's rebuild helper was tested on small matrices including a changed pattern; the full GLI run required only one factorization, so its time is not evidence of a speedup from reuse. The 3-shift linear sequence is evidence of the gain. No global cache, no new public ABI fields for cache itself. No whole-library sanitizer run; ASan/UBSan on three modified translation units and small test drivers passed 8 cases. MPI, complex, 64-bit-index, broad matrices, OOM and concurrent sessions remain untested for new features.

## Reproduction on private model
`reproduce.sh /path/to/built/lis /path/to/private/CSR /path/to/private/lanczos_rhs.bin /path/to/output`.
Input files: A.csr = D K D, B.csr = 1e9 D KG D, d.npy stores D diagonal; S.csr=A-.00107 B, S0.0009.csr and S0.00105.csr similarly. CSR binary: int32 n, nnz, rowptr[n+1], col[nnz], float64 values[nnz]. RHS file contains Lanczos vectors in float64; second RHS used. Shift .00107 maps to physical 1070000. Model 152592 DOF. The independent certificate uses SciPy/NumPy solely for validation, not solving. Portable certification script verified byte-identical CSV to the original local script on the clean run. `reproduce.sh` itself assembled from tested commands; not rerun as a single script. The driver prints an out-of-range mode diagnostic after writing all six requested vectors; solve status is success, this is a driver artifact.

## Measured results
- Fresh fully built LIS GLI 12.35 s process time, 899808 KiB peak RSS, 6/6 modes, max independent original residual 4.57075069493189e-11, K orthogonality 2.0725906031821768e-14, max relative eigenvalue error 1.872836637689307e-13. One factorization, MD integrated. Input preparation/scaling and certification excluded.
- Fresh-built three-shift public LIS sequence: no cache 18.33 s, cache 15.05 s. Max residual 2.4366982413265507e-11 both. Peak RSS 763100/789712 KiB. Single sequential trials, no median. Includes solves; different benchmark than eigenproblem.
- Historical earlier native full ILUTP 40.22 s with external MMD excluded; new 12.35 s includes MD. Historical MUMPS reference 5.19 s, different run conditions. Do not claim a controlled same-hardware ranking or superiority. Final native configuration is still slower than the MUMPS reference.
- Stage29 explicit cache and stage30 direct scatter logs included. The direct scatter avoids rematerializing permuted CSR on refresh; overall peak RSS is still governed by initial factorization. No claim of lower total peak RSS.

## Before upstream PR
Rebase to actual develop and reassess changed interfaces. Split functional fixes from performance features, replace environment flags with supported scoped options, decide ABI/version policy, complete sanitizer/OOM/MPI/complex/64-bit and broad sparse-matrix coverage. Add true supernodal numerical factorization only after a separate design and benchmark. Do not publish private matrices or vectors. No claim that all possible optimizations are exhausted.
