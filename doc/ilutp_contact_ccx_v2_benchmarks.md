# LIS native ILUTP: current results on contact matrices

This document reports the **current outcome**, not the sequence of
benchmark attempts. The implementation, options and supported modes
are described in [ILUTP](ilutp.md).

## Confirmed complete-fill result

On a private CCX Newton series with **60,265 DOF** and roughly
3.62–3.65 million stored entries per matrix, native LIS's opt-in
bounded complete-fill ILUTP achieved:

| Metric | Final observed result |
| --- | --- |
| Original A/b right-hand sides accepted | **15 / 15** |
| Distinct matrices A | **12** |
| FGMRES iterations | **1 per solve** |
| Original-system relative residual `||Ax-b||₂/||b||₂` | **4.631e-15 to 2.837e-12** |
| Largest stored L+U off-diagonal fill | **51,108,311 entries** |
| Largest peak process RSS | **799,024 KiB (about 780 MiB)** |
| Longest factor setup | **65.04 CPU seconds** |

Parameters: `-ilutp 1 -ilutp_full 1 -ilutp_threshold 0.1
-ilut_rcm 1 -ilut_row_scale 1 -ilutp_full_limit 60000000`,
serial double precision, FGMRES(40), original RHS, `rtol=1e-6`.
The independent physical-residual acceptance threshold was 1e-5.
Inputs were checked by SHA-256 and refactorized for every case,
including duplicates. The matrix in steps 1, 4, 6 and 12 was
identical, so **15 solved systems means 12 distinct A matrices**.

The original private inputs are not committed; public CI instead
tests synthetic and small numerical systems. The data are from a
single model whose **spring count is not independently verified**.

To reproduce one case with the original Matrix Market A/b pair:

```sh
python3 tools/run_native_ilutp_contact.py \
  --bench ./test/benchilutp \
  --case original /path/to/newton_A.mtx /path/to/newton_b.mtx \
  --full-limit 60000000 --full-tau 0.1 \
  --rcm 1 --row-scale 1 --restart 40 --maxiter 80 \
  --rtol 1e-6 --accept-residual 1e-5 --timeout 150 \
  --report result.json --csv result.csv
```

The factor-entry cap is **not an RSS bound**. The reported worst
fill was around 85.2% of its configured 60M cap.

## Unresolved: separately labelled 256-spring system

The original separately identified **snapshot 03** contains 60,265
unknowns and **3,639,839 stored A entries**. With the previous
limited-fill ILUTP settings, pivot regularization `1e-9` and
`drop=0.005, rate=12`, it experienced numerical breakdown and
an original-system relative residual near 1. This does **not**
prove that full fill was needed: the full mode also sets the
pivot regularization to zero. The reported `max(|D|)=1e9`
is consistent with inverse pivots hitting the regularization floor,
but is not proof of causality.

**No complete-fill success on this particular 256-spring input
has been demonstrated.** The tested 15-step series is different
and cannot be substituted for it.

Minimum discriminating comparison on the original 256-spring A/b
files, using the existing CLI without adding solver functionality:

```sh
A=/path/to/256_snapshot03_A.mtx
B=/path/to/256_snapshot03_b.mtx
OMP_NUM_THREADS=1 ./test/benchilutp "$A" "$B" \
  --maxiter 240 --restart 120 --rtol 1e-5 \
  --drop 0.005 --rate 12 --pivot-tol 0 \
  --rcm 1 --row-scale 1 --csv limited.csv

OMP_NUM_THREADS=1 ./test/benchilutp "$A" "$B" \
  --maxiter 240 --restart 120 --rtol 1e-5 \
  --drop 0.005 --rate 12 --pivot-tol 0 \
  --rcm 1 --row-scale 1 --full 1 --only-full 1 \
  --full-tau 0.1 --full-limit 60000000 --csv full.csv
```

The limited comparison emits thresholds 0.5 and 1.0, whereas the
full recipe above uses threshold 0.1. A matching-threshold control
(`--full-tau 0.5`) is needed before assigning a specific effect
to dropping/fill. Compare true residual, setup/solve time,
stored factor entries and peak RSS.

## Interpretation for upstream

**What we get now:** a working, default-off serial column-pivoted
preconditioner with optional bounded complete fill, validated on
15 real CCX Newton A/b solves (12 distinct matrices) and dedicated
regression tests.

**What we do not get yet:** verified full CCX nonlinear convergence,
a solution to the separately labelled 256-spring case, a production
solver comparison, or a guarantee of acceptable resources for all
contact models. The causal effect of fill versus pivot regularization
is still open. Current-head CI was green on 2026-11-14 after
pinning Linux runners to Ubuntu 24.04; a controlled baseline/PR
crosscheck confirmed separate eigensolver failures on Ubuntu 26.04
in both trees. The PR may proceed to **review** with those numerical
gaps explicitly stated; review is not proof of the original
256-spring case or end-to-end CCX convergence.
