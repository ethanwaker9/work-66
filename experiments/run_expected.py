import csv
import os
import random
from fractions import Fraction
from _path import RESULTS
from latsieve.common import next_prime, hnf_basis, cross
from latsieve.reduction import lll3w, fk_reduce, gauss2
from latsieve.zfs import select_frame, zfs_frame
from latsieve.common import Box
from latsieve.expected import predicted_nodes, count_nodes


def frames(r, l0, l1, L):
    hb = hnf_basis(r, l0, l1)
    u, v, x = lll3w(hb, 1.0 / L[0] ** 2, 1.0 / L[1] ** 2, 1.0 / L[2] ** 2)
    out = {"ZFS": zfs_frame(r, l0, l1, Box((-32, -32, 0), (31, 31, 31))), "LLL": (u, v, x), "LS": hb}
    e0, e1 = gauss2((r, 0), (l0 % r, 1), 1.0 / L[0] ** 2, 1.0 / L[1] ** 2)
    out["PS"] = ((e0[0], e0[1], 0), (e1[0], e1[1], 0), hb[2])
    return out


def main():
    rng = random.Random(11)
    lo = (-32, -32, 0)
    hi = (31, 31, 31)
    L = [hi[i] - lo[i] for i in range(3)]
    rows = []
    Q = 1 << 20
    for ex in (6, 9, 12, 15, 18, 21):
        for rep in range(3):
            r = next_prime(rng.randrange(2 ** ex, 2 ** (ex + 1)))
            l0 = rng.randrange(1, r)
            l1 = rng.randrange(1, r)
            for name, (a, c, e) in frames(r, l0, l1, L).items():
                if name in ("LS", "PS") and ex < 9:
                    pass
                pred = predicted_nodes(a, c, e, L)
                acc = [0, 0, 0]
                ntr = 120 if ex >= 12 else 40
                for _ in range(ntr):
                    t = [Fraction(rng.randrange(Q), Q) for _ in range(3)]
                    tau = tuple(t[0] * a[i] + t[1] * c[i] + t[2] * e[i] for i in range(3))
                    pl, li, po = count_nodes(a, c, e, lo, hi, tau)
                    acc[0] += pl
                    acc[1] += li
                    acc[2] += po
                meas = [v / ntr for v in acc]
                rows.append((ex, r, name, pred[0], meas[0], pred[1], meas[1], pred[2], meas[2]))
                print(ex, r, name, [round(p, 3) for p in pred], [round(m, 3) for m in meas], flush=True)
    with open(os.path.join(RESULTS, "expected_nodes.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["log2r", "r", "frame", "pred_planes", "meas_planes", "pred_lines", "meas_lines", "pred_points", "meas_points"])
        wr.writerows(rows)


if __name__ == "__main__":
    main()
