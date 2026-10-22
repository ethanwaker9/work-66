import csv
import os
import random
import sys
from math import pi, sqrt
from _path import RESULTS
from latsieve import nfs_box, next_prime
from latsieve.reduction import lll3w, gauss_ints
from latsieve.zfs import zfs_frame


def flag_nodes(a, c, L, r):
    E = L[0] * L[1] * L[2] / r
    n = (a[1] * c[2] - a[2] * c[1], a[2] * c[0] - a[0] * c[2], a[0] * c[1] - a[1] * c[0])
    n2 = E * sum(abs(a[i]) / L[i] for i in range(3))
    n3 = sum(L[i] * abs(n[i]) for i in range(3)) / r
    return E, n2, n3


def ellipsoid_nodes(r, l0, l1, box):
    h = [(box.hi[i] - box.lo[i]) / 2.0 + 0.5 for i in range(3)]
    b = lll3w(((r, 0, 0), (l0, 1, 0), (l1, 0, 1)), 1 / h[0] ** 2, 1 / h[1] ** 2, 1 / h[2] ** 2)
    adj = [(b[1][1] * b[2][2] - b[1][2] * b[2][1], b[1][2] * b[2][0] - b[1][0] * b[2][2], b[1][0] * b[2][1] - b[1][1] * b[2][0]),
           (b[2][1] * b[0][2] - b[2][2] * b[0][1], b[2][2] * b[0][0] - b[2][0] * b[0][2], b[2][0] * b[0][1] - b[2][1] * b[0][0]),
           (b[0][1] * b[1][2] - b[0][2] * b[1][1], b[0][2] * b[1][0] - b[0][0] * b[1][2], b[0][0] * b[1][1] - b[0][1] * b[1][0])]
    det = abs(sum(adj[0][i] * b[0][i] for i in range(3)))
    d1 = [adj[1][i] / det for i in range(3)]
    d2 = [adj[2][i] / det for i in range(3)]
    rho2 = 3.0
    n1 = 4 * pi / 3 * rho2 ** 1.5 * h[0] * h[1] * h[2] / det
    g11 = sum(d1[i] * d1[i] * h[i] * h[i] for i in range(3))
    g22 = sum(d2[i] * d2[i] * h[i] * h[i] for i in range(3))
    g12 = sum(d1[i] * d2[i] * h[i] * h[i] for i in range(3))
    n2 = pi * rho2 * sqrt(max(g11 * g22 - g12 * g12, 0.0))
    n3 = 2 * sqrt(rho2 * g22)
    return n1, n2, n3


def main(H, per, seed):
    box = nfs_box(H, H, H)
    L = [box.hi[i] - box.lo[i] for i in range(3)]
    rng = random.Random(seed)
    rows = []
    for ex in range(4, 27, 2):
        lats = []
        for _ in range(per):
            r = next_prime(rng.randrange(2 ** ex, 2 ** (ex + 1)))
            lats.append((r, rng.randrange(1, r), rng.randrange(1, r)))
        acc = {k: [0.0, 0.0, 0.0] for k in ("LS", "PS", "SE", "ZFS")}
        Esum = 0.0
        for r, l0, l1 in lats:
            E, n2, n3 = flag_nodes((r, 0, 0), (l0, 1, 0), L, r)
            Esum += E
            for i, v in enumerate((E, n2, n3)):
                acc["LS"][i] += v
            g0, g1, k0, k1 = gauss_ints(r, 0, l0 % r, 1, 1.0 / (L[0] * L[0]), 1.0 / (L[1] * L[1]))
            E, n2, n3 = flag_nodes((g0, g1, 0), (k0, k1, 0), L, r)
            for i, v in enumerate((E, n2, n3)):
                acc["PS"][i] += v
            for i, v in enumerate(ellipsoid_nodes(r, l0, l1, box)):
                acc["SE"][i] += v
            a, c, e = zfs_frame(r, l0, l1, box)
            E, n2, n3 = flag_nodes(a, c, L, r)
            for i, v in enumerate((E, n2, n3)):
                acc["ZFS"][i] += v
        for k, v in acc.items():
            rows.append((H, ex, k, Esum / per, v[0] / per, v[1] / per, v[2] / per))
            print(H, ex, k, round(Esum / per, 3), [round(x / per, 2) for x in v], flush=True)
    with open(os.path.join(RESULTS, f"nodes_H{H}.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["H", "log2r", "method", "E", "points", "lines", "planes"])
        wr.writerows(rows)


if __name__ == "__main__":
    H = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    main(H, per, 2026 + H)
