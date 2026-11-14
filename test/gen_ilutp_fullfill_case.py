#!/usr/bin/env python3
"""Synthetic nonsymmetric ill-scaled stress matrix; NOT CCX contact data.

Writes original A and matching manufactured b=A*x. Python stdlib only.
"""
import argparse
import math
from pathlib import Path


def generate(n, matrix, rhs):
    if n < 40:
        raise ValueError("n must be >= 40")
    matrix.parent.mkdir(parents=True, exist_ok=True)
    rhs.parent.mkdir(parents=True, exist_ok=True)
    nnz = sum(1 + (i > 0) + (i + 1 < n) + (i >= 11) + (i + 17 < n)
              for i in range(n))

    def x(j):
        return 1.0 + 0.1 * math.sin((j + 1) * 0.017)

    with matrix.open("w") as fm, rhs.open("w") as fb:
        fm.write(f"%%MatrixMarket matrix coordinate real general\n"
                 f"% synthetic stress, NOT CCX/contact\n{n} {n} {nnz}\n")
        fb.write(f"%%MatrixMarket matrix array real general\n"
                 f"% manufactured b=A*x for original A\n{n} 1\n")
        for i in range(n):
            scale = 10.0 ** (((i * 7) % 11) - 5)
            entries = [(i, 1e-11 if i % 127 == 0
                        else 4.0 + 0.1 * math.cos(i * 0.041))]
            if i > 0:
                entries.append((i - 1, -0.65))
            if i + 1 < n:
                entries.append((i + 1, 0.73))
            if i >= 11:
                entries.append((i - 11, -0.11))
            if i + 17 < n:
                entries.append((i + 17, 0.08))
            b = 0.0
            for j, value in entries:
                av = scale * value
                fm.write(f"{i+1} {j+1} {av:.17e}\n")
                b += av * x(j)
            fb.write(f"{b:.17e}\n")
    print(f"SYNTHETIC_FULLFILL n={n} nnz={nnz} (not CCX)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--rhs", type=Path, required=True)
    args = ap.parse_args()
    generate(args.n, args.matrix, args.rhs)
