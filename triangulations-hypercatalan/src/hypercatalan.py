#!/usr/bin/env python3
"""
Compute order‑2 hypercatalan numbers H(n) for dissections of the midpoint polygon.

Definition:
    H(n) = sum_{T} prod_{i=1..n} Catalan(deg_T(i)),
where T ranges over triangulations of a convex n‑gon and deg_T(i) is the number of diagonals
incident to vertex i in T.

This script uses a dynamic program over polygon intervals. Each subproblem carries only
the degrees at the two boundary vertices; internal pivot vertices are "closed" immediately
by multiplying with a Catalan factor. See paper/main.tex for the full proof of correctness.

Usage:
    python3 hypercatalan.py --n 10
    python3 hypercatalan.py --table 3 12
"""

from functools import lru_cache
from collections import defaultdict
from math import comb
import argparse

def catalan_numbers_up_to(m: int):
    """Return list C[0..m] of Catalan numbers as exact integers."""
    C = [0]*(m+1)
    for n in range(m+1):
        C[n] = comb(2*n, n) // (n+1) if n > 0 else 1
    return C

def hypercatalan_order2(n: int) -> int:
    """
    Compute H(n) = sum_T prod_{i=1..n} Catalan(deg_T(i)).
    Exact integer result. Runs fast for n up to ~25 on a laptop.
    """
    if n < 3:
        return 0

    # Degrees never exceed n-3; compute some slack.
    C = catalan_numbers_up_to(2*n)

    @lru_cache(maxsize=None)
    def dp(i: int, j: int):
        """
        For the subpolygon with consecutive vertices i..j (i < j),
        return a dict mapping (d_i, d_j) -> weight, where:
          - d_i: # diagonals from i to vertices inside [i..j]
          - d_j: same for j.
        'weight' already includes Catalan factors for all internal vertices of this subpolygon.
        """
        if j - i <= 1:
            # Edge: no diagonals, no internal vertices
            return {(0, 0): 1}

        res = defaultdict(int)

        # Choose pivot k forming triangle (i, k, j)
        for k in range(i+1, j):
            left = dp(i, k)   # maps (deg_i_left, deg_k_left) -> weight
            right = dp(k, j)  # maps (deg_k_right, deg_j_right) -> weight

            # In the full polygon, (i,k) is a diagonal unless k == i+1, and (k,j) unless k == j-1.
            add_i = 1 if k > i+1 else 0
            add_j = 1 if j > k+1 else 0
            add_k = add_i + add_j  # k gets both (unless boundary-adjacent)

            for (d_i_L, d_k_L), wL in left.items():
                for (d_k_R, d_j_R), wR in right.items():
                    d_i_tot = d_i_L + add_i
                    d_j_tot = d_j_R + add_j
                    d_k_tot = d_k_L + d_k_R + add_k

                    # Close the internal pivot k now
                    res[(d_i_tot, d_j_tot)] += wL * wR * C[d_k_tot]

        return dict(res)

    # Close the very outer endpoints 0 and n-1
    endpoint_map = dp(0, n-1)  # (deg_0, deg_{n-1}) -> weight
    total = sum(weight * C[d0] * C[dn] for (d0, dn), weight in endpoint_map.items())
    return total

def main():
    parser = argparse.ArgumentParser(description="Order‑2 hypercatalan numbers for midpoint dissections")
    parser.add_argument("--n", type=int, help="Compute H(n) for this n (n >= 3).")
    parser.add_argument("--table", nargs=2, type=int, metavar=("START", "END"),
                        help="Print H(n) for START..END (inclusive).")
    args = parser.parse_args()

    if args.n is not None:
        n = args.n
        if n < 3:
            print("n must be at least 3")
            return
        print(hypercatalan_order2(n))
        return

    if args.table is not None:
        a, b = args.table
        if a > b:
            a, b = b, a
        for n in range(a, b+1):
            print(n, hypercatalan_order2(n))
        return

    parser.print_help()

if __name__ == "__main__":
    main()
