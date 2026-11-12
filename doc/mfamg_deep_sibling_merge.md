# Opt-in MFAMG deep-level sibling grouping

MFAMG can form a long coarse hierarchy when many short branches attach
to a shared hub. This optional, per-solver grouping step reduces the
number of aggregates that choose the same strong neighboring anchor
before the multivector QR prolongator is constructed.

## Configuration

```text
-i gmres -p mfamg -mfamg_deep_sibling_merge 1
```

The option is `0` (default; unchanged MFAMG hierarchy) or `1`.
Invalid values are rejected. Configuration is stored per solver and
participates in transactional option-profile replacement and change
detection, with no environment-variable global state.

## Algorithm and safeguards

On serial multivector hierarchy levels with `level_id >= 4`, after
the existing connected aggregation:

- An anchor edge must satisfy the **diagonal-normalized strength**
  `|a_ij| / sqrt(|a_ii| |a_jj|) >= 0.05`. Non-finite, zero, or weak
  connections are excluded.
- Among qualifying edges, the strongest **normalized** score wins.
  Exact ties use the smallest aggregate ID. This is invariant to
  positive diagonal congruence scaling in exact arithmetic, but
  not guaranteed bitwise invariant under finite-precision rounding.
- Require at least **three eligible siblings** selecting one hub.
  The hub itself joins the group, together with as many of those
  siblings as fit the group-size cap.
- If siblings exceed the cap, leave the excess aggregates separate;
  **do not** pack disconnected sibling-only groups.
- The cap is `max(existing_target, min(local_n, 8 * mode_count))`
  using overflow-safe arithmetic. The depth threshold (4), factor
  (8), and minimum fan-in (3) are heuristic choices, not proven
  optimal parameters.
- Aggregation preserves the existing local multivector QR/rank
  construction and does not explicitly discard near-nullspace modes.

Because every selected sibling has a strong edge to the included
hub, each resulting group is connected through the hub **provided
that the original constituent aggregates are themselves connected**.
This addresses the disconnected sibling-group defect of the earlier
algorithm. It does **not** prove approximation or Galerkin energy
stability, nor does it guarantee fewer Krylov iterations.

Multiple MPI ranks always take the original distributed aggregation
path. A two-rank nine-level regression and a relinked negative-control
mutation test verify this serial-only guard.

The option is off by default. Production library code remains silent.

## Numerical regression coverage

`test/testmfamgdeepmerge.c` covers:

1. A 121x121 shifted star Laplacian, verifying a connected
   hub-and-leaves group (96 unknowns) and 25 remaining singleton
   leaves: 26 groups in total. A regular nine-node line confirms
   that a two-neighbor fan-in does not trigger grouping; a separate
   six-node SPD graph verifies that diagonal-normalized anchor
   ranking selects the correct hub instead of the raw-weight maximum.
2. A **1025x1025 SPD star-of-chains** operator (32 arms, 32 unknowns per
   arm, diagonal shift `1e-3 I`) with two smooth near-nullspace modes.
   Preconditioner construction verifies that depth `>=4` is reached,
   grouping actually fires, and the next-level dimension decreases.
3. Full `lis_solve()` with GMRES/MFAMG for option states
   `0 -> 1 -> 0` on **the same LIS_SOLVER**. Each solve checks success,
   finite positive iteration count and independently evaluates the
   physical relative residual `||b-Ax||_2/||b||_2 < 1e-8`.
4. The profile transition back to zero must restore the original
   hierarchy dimensions, including level count; internal per-instance
   counters verify that no merge runs when disabled.
5. An MPI **two-rank numerical solve** on a larger SPD distributed
   line exercises a hierarchy deeper than level 4 and verifies the
   serial-only guard: zero merge attempts, identical hierarchy
   dimensions and physical residuals for `0 -> 1 -> 0`.
   A separate validation-only negative control removes the outer
   `nprocs==1` check: the test must then fail upon reaching the
   helper's independent rank guard. The negative control is never
   included in the submitted production code.

The internal `lis_precon_mfamg_get_deep_merge_info` test accessor
reports attempts, successful events, reduced aggregate count and
deepest affected level. It is **not part of the public LIS API**.

### Single-hub numerical regression (1025 unknowns)

| Setting | Hierarchy levels | Unknowns at level 5 | Groups removed | GMRES iter. | Physical relative residual |
| --- | ---: | ---: | ---: | ---: | ---: |
| Disabled | 13 | 30 | 0 | 37 | 7.63e-10 |
| Rooted grouping | 7 | 18 | 24, deepest level 5 | 38 | 7.00e-10 |
| Disabled again | 13 | 30 | 0 | 37 | 7.63e-10 |

The connected grouping does **not** guarantee that iteration counts
remain unchanged: this example needs one additional GMRES iteration.

### Single-hub benchmark (16385 unknowns, 32 arms of length 512)

Seven consecutive complete `lis_solve()` calls per option on a
single Linux CI runner (CPU time measured via `clock()`):

| Setting | Levels | Iterations | Build CPU mean +/- SD (ms) | Solve CPU mean +/- SD (ms) |
| --- | ---: | ---: | ---: | ---: |
| Disabled | 15 | 51 | 2.214 +/- 1.035 | 29.626 +/- 0.889 |
| Rooted grouping | 8 | 51 | 1.761 +/- 0.045 | 29.181 +/- 0.476 |
| Disabled again | 15 | 51 | 1.810 +/- 0.046 | 29.592 +/- 0.729 |

The difference of approximately 0.45 ms in mean solve time is
within the scale of observed variability; it is **not a demonstrated
general speedup**. Benchmark and source:
https://github.com/3rav/lis/actions/runs/43174125966

### Multi-hub benchmark (12352 unknowns)

A separate graph has **64 hubs on an 8x8 lattice**, each connected
to **24 arms of length 8**, totaling 12352 unknowns. Seven repeated
runs in one Linux CI job:

| Setting | Levels | Unknowns at level 5 | Iterations | Solve CPU mean +/- SD (ms) |
| --- | ---: | ---: | ---: | ---: |
| Disabled | 13 | 1024 | 18 | 17.321 +/- 1.663 |
| Rooted grouping | 9 | 256 | 19 | 15.668 +/- 0.034 |
| Disabled again | 13 | 1024 | 18 | 16.479 +/- 0.187 |

Deep merging removes 771 aggregate groups in this case, but requires
**one additional GMRES iteration**. CPU time decreases on this
synthetic example; differences between first/last disabled timings
and runner variability prevent a universal performance conclusion.

Benchmark and source:
https://github.com/3rav/lis/actions/runs/43174520733

The benchmarks do not prove performance improvement on arbitrary
unstructured finite-element meshes, contact problems or distributed
MPI operators. They are diagnostic, not production tuning claims.

## Memory and platform validation

The integration test destroys all preconditioners, vectors, solver and
matrix after repeated option transitions. It also requires exact
restoration of the baseline iteration count and physical residual
in the deterministic `0 -> 1 -> 0` sequence. Run it with sanitizers
(ASan/LSan/UBSan) as well as the regular Linux serial/OpenMP/quad,
MPI two-rank and Windows UCRT64 `make check` configurations to check
memory safety on the exercised paths.

This change does not modify CalculiX, the default MFAMG configuration
or the already upstreamed terminal block-Schur solver.

## ABI note

The new `mfamg_deep_sibling_merge` field is appended to the end of
the visible `LIS_SOLVER_STRUCT`, preserving offsets of existing
members. The size of this public structure nevertheless increases,
so binary compatibility with programs compiled against an older
layout is **not guaranteed**. Rebuild dependent applications against
matching LIS headers and library. No public function is added.
