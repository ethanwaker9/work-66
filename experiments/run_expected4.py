import csv
import os
import random
from fractions import Fraction
from math import comb, gamma, pi, sqrt
from _path import RESULTS
from latsieve.common import next_prime
from latsieve.expected import general_node_counts, pluecker_l1, lll_exact


def constants():
    rows = []
    alpha = 1 / (0.99 - 0.25)
    for t in range(2, 7):
        cs = [sqrt(comb(t, j)) * alpha ** (j * (t - j) / 4) for j in range(1, t)]
        ball = (pi * t / 4) ** (t / 2) / gamma(t / 2 + 1)
        rows.append((t, ";".join(f"{c:.3f}" for c in cs), f"{ball:.3f}"))
    with open(os.path.join(RESULTS, "constants.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["t", "c_tj", "ball_over_box"])
        wr.writerows(rows)
    print(rows)


def main():
    rng = random.Random(21)
    lo = (-8, -8, -8, 0)
    hi = (7, 7, 7, 7)
    L = [hi[i] - lo[i] for i in range(4)]
    rows = []
    for ex in (6, 9, 12, 15):
        for rep in range(2):
            r = next_prime(rng.randrange(2 ** ex, 2 ** (ex + 1)))
            lam = [rng.randrange(1, r) for _ in range(3)]
            hnf = [(r, 0, 0, 0), (lam[0], 1, 0, 0), (lam[1], 0, 1, 0), (lam[2], 0, 0, 1)]
            red = lll_exact(hnf, [Fraction(1, L[i] ** 2) for i in range(4)])
            E = Fraction(L[0] * L[1] * L[2] * L[3], r)
            for name, B in (("LLL", red), ("LS", hnf)):
                pred = [float(E * pluecker_l1(B, L, k)) for k in range(4)]
                acc = [0, 0, 0, 0]
                ntr = 24
                for _ in range(ntr):
                    s = [Fraction(rng.randrange(1 << 16), 1 << 16) for _ in range(4)]
                    tau = tuple(sum(s[j] * B[j][i] for j in range(4)) for i in range(4))
                    c = general_node_counts(B, lo, hi, tau)
                    for k in range(4):
                        acc[k] += c[k]
                meas = [a / ntr for a in acc]
                for k in range(4):
                    rows.append((ex, r, name, k + 1, pred[k], meas[k]))
                print(ex, r, name, [round(p, 2) for p in pred], [round(m, 2) for m in meas], flush=True)
    with open(os.path.join(RESULTS, "expected_nodes_t4.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["log2r", "r", "frame", "level", "predicted", "measured"])
        wr.writerows(rows)


if __name__ == "__main__":
    constants()
    main()
