import csv
import os
import random
import sys
from _path import RESULTS
from latsieve import nfs_box, next_prime
from latsieve.reduction import lll3_ints
from latsieve.zfs import zfs_frame


def classify(r, l0, l1, box, a, c, e):
    L = [box.hi[i] - box.lo[i] for i in range(3)]
    if a == (r, 0, 0) and c == (l0, 1, 0) and e == (l1, 0, 1):
        return "LS"
    if a[2] == 0 and c[2] == 0 and e == (l1, 0, 1):
        return "coordinate planes"
    red = lll3_ints(r, 0, 0, l0, 1, 0, l1, 0, 1, 1.0 / (L[0] * L[0]), 1.0 / (L[1] * L[1]), 1.0 / (L[2] * L[2]))
    vecs = [tuple(red[0:3]), tuple(red[3:6]), tuple(red[6:9])]
    if a == vecs[0]:
        return "reduced, line u"
    return "reduced, line v or x"


def main(H, per, seed):
    box = nfs_box(H, H, H)
    L = [box.hi[i] - box.lo[i] for i in range(3)]
    rng = random.Random(seed)
    rows = []
    kinds = ["reduced, line u", "reduced, line v or x", "coordinate planes", "LS"]
    for ex in range(4, 27, 2):
        cnt = {k: 0 for k in kinds}
        nz = 0
        row = 0
        Es = 0.0
        for _ in range(per):
            r = next_prime(rng.randrange(2 ** ex, 2 ** (ex + 1)))
            l0, l1 = rng.randrange(1, r), rng.randrange(1, r)
            a, c, e = zfs_frame(r, l0, l1, box)
            cnt[classify(r, l0, l1, box, a, c, e)] += 1
            nz += sum(1 for v in a if v)
            row += 1 if a[1] == 0 and a[2] == 0 else 0
            Es += L[0] * L[1] * L[2] / r
        rows.append((H, ex, Es / per) + tuple(cnt[k] / per for k in kinds) + (nz / per, row / per))
        print(rows[-1], flush=True)
    with open(os.path.join(RESULTS, f"frames_H{H}.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["H", "log2r", "E", "reduced_u", "reduced_vx", "coordinate_planes", "LS", "nonzero_a", "line_in_row"])
        wr.writerows(rows)


if __name__ == "__main__":
    H = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    main(H, per, 2026 + H)
