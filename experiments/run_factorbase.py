import csv
import gc
import math
import os
import random
import sys
import time
import tracemalloc
from _path import RESULTS
from latsieve import METHODS, nfs_box, primes_between


def main(H, B, seed):
    box = nfs_box(H, H, H)
    rng = random.Random(seed)
    primes = primes_between(2, B)
    lats = [(r, rng.randrange(r), rng.randrange(r)) for r in primes]
    weights = [max(1, round(math.log2(r))) for r, _, _ in lats]
    rows = []
    ref = None
    for name, f in METHODS:
        tt = None
        for _ in range(2):
            S = [0] * box.size
            t0 = time.perf_counter()
            for (r, a, b), w in zip(lats, weights):
                f(r, a, b, box, S, w)
            el = time.perf_counter() - t0
            tt = el if tt is None else min(tt, el)
        gc.collect()
        peak = 0
        sample = lats[:: max(1, len(lats) // 400)]
        for r, a, b in sample:
            T = [0] * box.size
            gc.collect()
            tracemalloc.start()
            f(r, a, b, box, T)
            peak = max(peak, tracemalloc.get_traced_memory()[1])
            tracemalloc.stop()
        if ref is None and name == "LS":
            ref = S
        rows.append((name, tt, S, peak))
        print(name, round(tt, 3), peak, flush=True)
    total = sum(ref)
    thr = sorted(ref, reverse=True)[max(0, box.size // 1000)]
    cand_ref = {i for i, v in enumerate(ref) if v >= thr}
    out = []
    for name, tt, S, peak in rows:
        mass = sum(S) / total
        cand = {i for i, v in enumerate(S) if v >= thr}
        recall = len(cand & cand_ref) / max(1, len(cand_ref))
        out.append((H, B, len(lats), name, tt, mass, recall, peak))
    with open(os.path.join(RESULTS, f"factorbase_H{H}_B{B}.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["H", "B", "primes", "method", "total_time_s", "hit_mass_fraction", "candidate_recall", "peak_bytes"])
        wr.writerows(out)
    for o in out:
        print(o)


if __name__ == "__main__":
    H = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 2 ** 18
    main(H, B, 99)
