# FGMRES projected-residual convergence units

FGMRES's Arnoldi/Givens least-squares estimate,
`abs(s[i+1])`, is an **absolute** residual 2-norm. The LIS
convergence condition `nrm2_r` (the only currently supported FGMRES condition) uses a **relative** tolerance.

The projected estimate is only a candidate for convergence. LIS already
recomputes the explicit physical residual `b - A*x` before accepting a
successful solve; that certification is unchanged.

## Correct normalization

For `-conv_cond nrm2_r` (default), LIS defines

```
solver->bnrm = 1 / ||r0||_2,     r0 = b - A*x0
```

Therefore multiply the absolute projected estimate by `solver->bnrm`
before monitoring, convergence-candidate selection, or the no-improvement
heuristic.

**Note:** FGMRES currently rejects `-conv_cond nrm2_b` and `nrm1_b`
at solver setup (`Option conv_cond is not implemented for solver FGMRES`).
This PR intentionally does not add support for those modes.

The least-squares coefficients in `s[]` **must not be rescaled**; doing
so would change back-substitution and the computed solution.

The unsupported `nrm1_b` mode uses a different mixed absolute/relative
condition and is deliberately out of scope.

## Scope

This is a correction of the default convergence-candidate units in both
regular and quad-precision FGMRES, on single rank and with MPI.
It does **not** introduce a new environment variable, global state,
public API, or MFAMG/CalculiX dependency.

The explicit physical true-residual check remains the final authority.
There is no claim of general performance or nonlinear-contact improvement.

## Regression

`test/testfgmresprojected.c` checks a two-equation diagonal system with:

- RHS scales `1e4` and `1e-2`;
- the supported `nrm2_r` condition;
- zero and nonzero initial guesses (`-initx_zeros false`);
- an independent calculation of the *first* Arnoldi/Givens relative
  projected residual, verified using the public solver monitor callback;
- an independently computed physical relative residual for the returned
  solution;
- both default and quad FGMRES when LIS is configured with `--enable-quad`.

The regression does not assume a fixed iteration count. It catches the
original units mismatch even when both old and new implementations
eventually produce a correct physical solution.

Run the test using `make -C test testfgmresprojected` and
`./test/testfgmresprojected`, or as part of `make check`.

## Compatibility

- No new solver options or application-visible global state.
- No changes to MFAMG, FGCRO-DR, or contact-specific logic.
- The projected-residual *candidate* may be accepted earlier or later,
  but never bypasses explicit true-residual certification.
