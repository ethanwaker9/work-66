import random
import sys
sys.path.insert(0, __file__.rsplit("/experiments", 1)[0])
from latsieve import METHODS, nfs_box, reference_points

random.seed(7)
stats = {name: [0, 0, 0] for name, _ in METHODS}
for trial in range(300):
    box = nfs_box(random.choice([4, 8, 16]), random.choice([4, 8, 16]), random.choice([4, 8, 16]))
    r = random.choice([2, 3, 7, 31, 101, 257, 1009, 4099, 65537, 1000003])
    l0 = random.randrange(r)
    l1 = random.randrange(r)
    ref = sorted(reference_points(r, l0, l1, box))
    refset = set(ref)
    for name, f in METHODS:
        S = [0] * box.size
        f(r, l0, l1, box, S)
        got = [i for i, v in enumerate(S) for _ in range(v)]
        stats[name][0] += sum(1 for i in set(got) if i in refset)
        stats[name][1] += len(ref)
        stats[name][2] += len(got) - len(set(got)) + sum(1 for i in set(got) if i not in refset)
for name, (f, t, bad) in stats.items():
    print(name, "found", f, "of", t, "fraction", round(f / t, 4), "invalid_or_duplicate", bad)
