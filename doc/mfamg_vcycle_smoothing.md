# Optional MFAMG V-cycle smoothing passes

The default MFAMG V-cycle performs one pre-smoothing pass from zero
and one post-smoothing correction. This PR adds a per-solver option to
allow two or three passes on **each** side:

```text
-i gmres -p mfamg -mfamg_vcycle_smooth_steps 2
```

Accepted values: `1` (default, existing behavior), `2`, `3`.
Invalid values are rejected. The option is parsed per solver, not from
an environment variable; option-profile replacement and change
detection preserve and compare it.

When more than one pre-smoothing pass is requested, only the **first**
starts from the zero iterate; subsequent passes use the current
approximation. Every post-smoothing pass corrects the prolongated
iterate. The recursive coarse solve and transfer remain unchanged.

Increasing smoothing may increase work per V-cycle, and there is no
claim of better total wall time, contact convergence, or better
spectral behavior without further benchmarks. This is strictly opt-in
and preserves the default operation of MFAMG.

The regression uses a deterministic 24x24 positive definite system to
exercise passes 1, 2 and 3 and checks the physical preconditioner
defect and valid configuration. It is run via `make check`.

## Regression evidence

The deterministic 96-unknown SPD regression deliberately builds a
hierarchy with at least four levels, so V-cycle recursion is exercised
beyond the first fine level. It checks:

- The computed correction and relative physical defect of the default
  configuration are **bitwise identical** to explicit `steps=1` on
  the tested build (not a claim of bitwise identity with historical
  releases or other compilers).
- Configurations `1`, `2`, and `3` yield finite relative physical
  defects strictly below `0.90` on this matrix.
- Additional smoothing changes the computed correction measurably.

The `0.90` bound is a regression safety check on this test matrix,
not a general convergence bound for arbitrary matrices.

The option is exercised here with GMRES. This PR does not establish
that MFAMG configured this way is a symmetric positive-definite
preconditioner suitable for conjugate gradients; users needing a
symmetry-dependent Krylov method must validate that property separately.
