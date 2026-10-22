import csv
import gc
import os
import random
import sys
import time
import tracemalloc
from _path import RESULTS
from latsieve import METHODS, nfs_box, next_prime, reference_points


def run(H, exps, per, seed):
    box = nfs_box(H, H, H)
    rng = random.Random(seed)
    rows = []
    for ex in exps:
        lats = []
        for _ in range(per):
            r = next_prime(rng.randrange(2 ** ex, 2 ** (ex + 1)))
            lats.append((r, rng.randrange(1, r), rng.randrange(1, r)))
        refs = [set(reference_points(r, a, b, box)) for r, a, b in lats]
        for name, f in METHODS:
            tt = 0.0
            found = 0
            total = 0
            peak = 0
            for (r, a, b), ref in zip(lats, refs):
                reps = []
                for _ in range(3):
                    S = [0] * box.size
                    t0 = time.perf_counter()
                    f(r, a, b, box, S)
                    reps.append(time.perf_counter() - t0)
                tt += sorted(reps)[1]
                found += sum(1 for i in ref if S[i])
                total += len(ref)
                S = [0] * box.size
                gc.collect()
                tracemalloc.start()
                f(r, a, b, box, S)
                peak = max(peak, tracemalloc.get_traced_memory()[1])
                tracemalloc.stop()
            rows.append((H, ex, name, tt / per * 1e3, found / max(total, 1), peak, total / per))
            print(H, ex, name, round(tt / per * 1e3, 4), round(found / max(total, 1), 4), peak, flush=True)
    return rows


if __name__ == "__main__":
    H = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    exps = list(range(4, 27, 2))
    rows = run(H, exps, per, 2026 + H)
    with open(os.path.join(RESULTS, f"bench_H{H}.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["H", "log2r", "method", "time_ms", "completeness", "peak_bytes", "points"])
        wr.writerows(rows)
