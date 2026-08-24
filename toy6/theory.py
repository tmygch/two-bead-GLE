#!/usr/bin/env python3
"""
theory_toy6.py

Exact theoretical memory kernel for the six-bead elastic network.

Outputs
-------
Toy6_N6_mu(t)_iXXX_jYYY_genth.dat
Toy6_N6_mu(s)_iXXX_jYYY_genth.dat
"""

import argparse
import numpy as np
from scipy.linalg import expm

def remove_row_col(A, k):
    return np.delete(np.delete(A, k, axis=0), k, axis=1)
def insert_zero_rowcol(A, k):
    A = np.insert(A, k, 0.0, axis=0)
    return np.insert(A, k, 0.0, axis=1)
def laplace_transform(mu, s, t):
    """Numerical Laplace transform using the trapezoidal rule."""
    dt = np.diff(t)
    E0 = np.exp(-np.outer(s, t[:-1]))
    E1 = np.exp(-np.outer(s, t[1:]))
    return np.sum(
        0.5 * dt * (mu[:-1] * E0 + mu[1:] * E1),
        axis=1,
    )
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--i", type=int, default=2,
                        help="tagged bead i (1-based, default: 2)")
    parser.add_argument("--j", type=int, default=6,
                        help="tagged bead j (1-based, default: 6)")
    args = parser.parse_args()

    # Parameters
    N = 6
    i2, j2 = args.i, args.j       # manuscript notation: 1,...,N
    i, j = i2 - 1, j2 - 1         # Python notation: 0,...,N-1
    gamma = 1.0
    k = 1.0
    rate = k / gamma

    # Kirchhoff matrix
    L = np.array([
        [ 1, -1,  0,  0,  0,  0],
        [-1,  3, -1,  0, -1,  0],
        [ 0, -1,  4, -1, -1, -1],
        [ 0,  0, -1,  1,  0,  0],
        [ 0, -1, -1,  0,  2,  0],
        [ 0,  0, -1,  0,  0,  1],
    ], dtype=float)

    # Coordinate transformation:
    #   tilde r_i = r_i - r_j
    #   tilde r_j = (r_i + r_j)/2
    P = np.eye(N)
    P[i, j] = -1.0
    P[j, i] = 0.5
    P[j, j] = 0.5

    TL = P @ L @ np.linalg.inv(P)
    TL_d = remove_row_col(TL, i)              # tilde L'

    # Generalized inverse (tilde L')^-
    keep = [n for n in range(N) if n not in (i, j)]
    L_dd_inv = np.linalg.inv(L[np.ix_(keep, keep)])   # (L'')^-1
    j_reduced = j - 1 if j > i else j
    TL_d_ginv = insert_zero_rowcol(L_dd_inv, j_reduced)

    # Coupling row/column vectors
    Tlib_d = np.delete(TL[i, :], i)
    Tlbi_d = np.delete(TL[:, i], i)

    # Check generalized-inverse relations
    err1 = np.max(np.abs(TL_d @ TL_d_ginv @ TL_d - TL_d))
    err2 = np.max(np.abs(TL_d_ginv @ TL_d @ TL_d_ginv - TL_d_ginv))
    print(f"max|L' G L' - L'| = {err1:.3e}")
    print(f"max|G L' G - G|   = {err2:.3e}")

    # Memory kernel mu(t)
    t_values = np.logspace(-3.5, 3.9, 300)
    mu_values = np.array([
        rate * Tlib_d @ expm(-rate * t * TL_d) @ TL_d_ginv @ Tlbi_d
        for t in t_values
    ])

    np.savetxt(
        f"Toy6_N6_mu(t)_i{i2:03d}_j{j2:03d}_genth.dat",
        np.column_stack((t_values, mu_values)),
        fmt="%.17g",
    )
    # Laplace-transformed kernel mu_hat(s)
    s_values = np.logspace(-4.6, 2.6, 30)
    hat_mu = laplace_transform(mu_values, s_values, t_values)

    np.savetxt(
        f"Toy6_N6_mu(s)_i{i2:03d}_j{j2:03d}_genth.dat",
        np.column_stack((s_values, hat_mu)),
        fmt="%.17g",
    )

if __name__ == "__main__":
    main()
