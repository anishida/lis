#!/usr/bin/env python3
"""Manufacture a symmetric-indefinite Stokes-type block saddle-point problem.

Original unscaled operator: H = [[-delta*I, gamma*B], [gamma*B.T, K]],
where K is a positive definite 5-point velocity Laplacian and B is
a one-sided difference / discrete-divergence matrix of full row rank.
The published Stokes saddle structure is standard, but this particular
matrix is *manufactured*, not a genuine IFISS or FE export.

The exported A = D H D with positive diagonal units D is congruent to
H, and hence remains symmetric indefinite with the same inertia.
The pressure-first ordering exposes exactly zero pressure pivots at
delta=0, without introducing a nullspace; the x-velocity term in B
has an identity diagonal so B has full row rank.

Outputs a Matrix Market coordinate real general matrix and the exact
manufactured RHS b=A*x, preserving the *original* matrix and RHS.
"""
import argparse
import math
from pathlib import Path


def make(side, contrast, units, ordering, delta, gamma, diag, matrix, rhs):
    if side < 2 or contrast < 1 or delta < 0 or gamma <= 0 or diag <= 4:
        raise ValueError("invalid grid or scaling parameters")
    cells = side * side
    n = 3 * cells
    pressure = lambda i: i
    ux = lambda i: cells + i
    uy = lambda i: 2 * cells + i
    entries = {}

    def add(i, j, value):
        if value:
            entries[i, j] = entries.get((i, j), 0.0) + value

    for cell in range(cells):
        x, y = cell % side, cell // side
        add(pressure(cell), pressure(cell), -delta)
        for velocity, idx, coeff in ((ux, cell, 1.0), (uy, cell, 1.0),
                                     (ux, cell - 1, -1.0 if x > 0 else 0.0),
                                     (uy, cell - side, -1.0 if y > 0 else 0.0)):
            if coeff:
                i, j = pressure(cell), velocity(idx)
                add(i, j, gamma * coeff)
                add(j, i, gamma * coeff)
        for velocity in (ux, uy):
            i = velocity(cell)
            add(i, i, diag)
            if x + 1 < side: add(i, velocity(cell + 1), -1.0)
            if x > 0: add(i, velocity(cell - 1), -1.0)
            if y + 1 < side: add(i, velocity(cell + side), -1.0)
            if y > 0: add(i, velocity(cell - side), -1.0)

    weights = [1.0] * n
    root = math.sqrt(contrast)
    for cell in range(cells):
        checker = ((cell % side) + (cell // side)) & 1
        if units == "pressure_checker":
            weights[pressure(cell)] = contrast if checker else 1.0 / contrast
        elif units == "pressure_uniform":
            weights[pressure(cell)] = contrast
        elif units == "pressure_weak":
            weights[pressure(cell)] = 1.0 / contrast
        elif units == "mixed_checker":
            weights[pressure(cell)] = contrast if checker else 1.0 / contrast
            weights[ux(cell)] = root if not checker else 1.0 / root
            weights[uy(cell)] = root if checker else 1.0 / root
        elif units == "velocity_checker":
            weights[ux(cell)] = root if checker else 1.0 / root
            weights[uy(cell)] = root if not checker else 1.0 / root
        elif units != "unity":
            raise ValueError("unknown unit mode")

    if ordering == "pressure_first":
        order = list(range(n))
    elif ordering == "velocity_first":
        order = list(range(cells, n)) + list(range(cells))
    elif ordering == "interleaved":
        order = [k for cell in range(cells)
                 for k in (pressure(cell), ux(cell), uy(cell))]
    else:
        raise ValueError("unknown ordering")

    inverse = [0] * n
    for i, old in enumerate(order):
        inverse[old] = i

    transformed = {}
    for (i, j), value in entries.items():
        new_value = value * weights[i] * weights[j]
        if not math.isfinite(new_value):
            raise ValueError("non-finite matrix coefficient")
        if new_value:
            transformed[inverse[i], inverse[j]] = new_value

    true_x = [math.cos(k * 0.17) + 0.25 * math.sin(k * 0.11) for k in range(n)]
    b = [0.0] * n
    for (i, j), value in transformed.items():
        b[i] += value * true_x[j]
    if not all(math.isfinite(value) for value in b):
        raise ValueError("non-finite original-system RHS")

    matrix, rhs = Path(matrix), Path(rhs)
    matrix.parent.mkdir(parents=True, exist_ok=True)
    rhs.parent.mkdir(parents=True, exist_ok=True)
    with matrix.open("w", encoding="ascii") as file:
        file.write("%%MatrixMarket matrix coordinate real general\n")
        file.write("% Stokes-type symmetric indefinite full-rank saddle point\n")
        file.write(f"{n} {n} {len(transformed)}\n")
        for (i, j), value in sorted(transformed.items()):
            file.write(f"{i + 1} {j + 1} {value:.17g}\n")
    with rhs.open("w", encoding="ascii") as file:
        file.write("%%MatrixMarket matrix array real general\n")
        file.write(f"{n} 1\n")
        for value in b:
            file.write(f"{value:.17g}\n")
    return {"n": n, "nnz": len(transformed), "side": side,
            "unit_contrast": contrast, "units": units,
            "ordering": ordering, "gamma": gamma, "delta": delta,
            "diag": diag}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side", type=int, default=6)
    parser.add_argument("--contrast", type=float, default=100)
    parser.add_argument("--units", choices=["unity", "pressure_checker",
                        "pressure_uniform", "pressure_weak", "mixed_checker",
                        "velocity_checker"], default="pressure_checker")
    parser.add_argument("--ordering", choices=["pressure_first",
                        "velocity_first", "interleaved"],
                        default="pressure_first")
    parser.add_argument("--delta", type=float, default=0.0)
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--diag", type=float, default=4.5)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--rhs", required=True)
    options = parser.parse_args()
    print(make(**vars(options)))
