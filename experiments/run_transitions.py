import csv
import json
import os
import random
from _path import RESULTS
from latsieve.common import next_prime
from latsieve.transitions import transition_set, hypothesis_h, diagonal_family, terminal_vectors, positions_order, transition_set_by_queries, alpha_beta_from_hnf


def exhaustive_table():
    rows = []
    for I in range(1, 7):
        for J in range(I, 7):
            for mult in (3, 12):
                r = next_prime(mult * I * J + 1)
                best = 0
                cnt = 0
                for alpha in range(1, r):
                    for beta in range(r):
                        if not hypothesis_h(r, alpha, beta, I, J):
                            continue
                        cnt += 1
                        k = len(transition_set(r, alpha, beta, I, J))
                        best = max(best, k)
                rows.append((I, J, r, cnt, best))
                print("table", I, J, r, cnt, best, flush=True)
    with open(os.path.join(RESULTS, "transitions_max.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["I", "J", "r", "lattices", "max_T"])
        wr.writerows(rows)


def distribution():
    rng = random.Random(5)
    rows = []
    for I in (16, 32, 64):
        for ratio in (2, 16, 128):
            r = next_prime(ratio * I * I + 1)
            hist = {}
            n = 0
            while n < 600:
                alpha = rng.randrange(1, r)
                beta = rng.randrange(r)
                if not hypothesis_h(r, alpha, beta, I, I):
                    continue
                n += 1
                k = len(transition_set(r, alpha, beta, I, I))
                hist[k] = hist.get(k, 0) + 1
            for k, v in sorted(hist.items()):
                rows.append((I, ratio, r, k, v))
            print("dist", I, ratio, sorted(hist.items()), flush=True)
    with open(os.path.join(RESULTS, "transitions_hist.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["I", "ratio", "r", "T", "count"])
        wr.writerows(rows)


def family():
    rows = []
    for I in range(2, 41, 2):
        for J in (I, I + 3):
            M = 2 * (I + J)
            r = next_prime(4 * M * (I + J) + 1)
            a, b = diagonal_family(I, J, r)
            k = len(transition_set(r, a, b, I, J))
            rows.append((I, J, r, k, min(I, J) + (1 if I == J else 2), 3 * min(I, J), min(I, J) + 3))
    with open(os.path.join(RESULTS, "transitions_family.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["I", "J", "r", "T", "predicted", "upper3min", "gs_bound"])
        wr.writerows(rows)
    print("family done", flush=True)


def terminal_check():
    rng = random.Random(8)
    ok = 0
    tot = 0
    for _ in range(1500):
        I = rng.randint(1, 12)
        J = rng.randint(1, 12)
        r = next_prime(rng.randint(2 * I * J + 2, 80 * I * J + 50))
        alpha = rng.randrange(1, r)
        beta = rng.randrange(r)
        if not hypothesis_h(r, alpha, beta, I, J):
            continue
        tot += 1
        T = transition_set(r, alpha, beta, I, J)
        if T == terminal_vectors(r, alpha, beta, I, J) and len(T) <= 3 * min(I, J):
            ok += 1
    with open(os.path.join(RESULTS, "transitions_terminal.json"), "w") as fh:
        json.dump({"instances": tot, "agree": ok}, fh)
    print("terminal", tot, ok, flush=True)


def partitions():
    out = {}
    I = J = 12
    M = 2 * (I + J)
    r = next_prime(4 * M * (I + J) + 1)
    for tag, (alpha, beta) in (("family", diagonal_family(I, J, r)), ("random", (7919 % r, 3571 % r))):
        vals = positions_order(r, alpha, beta, I, J)
        n = len(vals)
        vecs = {}
        grid = [[0] * J for _ in range(I)]
        for i in range(n):
            z, a, b = vals[i]
            z2, a2, b2 = vals[(i + 1) % n]
            v = (a2 - a, b2 - b, (z2 - z) % r)
            if v not in vecs:
                vecs[v] = len(vecs)
            grid[a][b] = v
        order = sorted(vecs, key=lambda v: v[2])
        idx = {v: i for i, v in enumerate(order)}
        out[tag] = {"r": r, "alpha": alpha, "beta": beta, "grid": [[idx[grid[a][b]] for b in range(J)] for a in range(I)], "vectors": order}
    with open(os.path.join(RESULTS, "transitions_partition.json"), "w") as fh:
        json.dump(out, fh)
    print("partition", {k: len(v["vectors"]) for k, v in out.items()}, flush=True)


def query_check():
    rng = random.Random(4)
    agree = 0
    total = 0
    calls = 0
    for _ in range(300):
        I = rng.randint(1, 12)
        J = rng.randint(1, 12)
        r = next_prime(rng.randint(2 * I * J + 2, 60 * I * J + 50))
        l0 = rng.randrange(1, r)
        l1 = rng.randrange(1, r)
        alpha, beta = alpha_beta_from_hnf(r, l0, l1)
        if not hypothesis_h(r, alpha, beta, I, J):
            continue
        total += 1
        T2, q = transition_set_by_queries(r, l0, l1, I, J)
        calls += q
        if T2 == transition_set(r, alpha, beta, I, J):
            agree += 1
    with open(os.path.join(RESULTS, "transitions_queries.json"), "w") as fh:
        json.dump({"instances": total, "agree": agree, "mean_flag_enumerations": calls / max(1, total)}, fh)
    print("queries", total, agree, calls / max(1, total), flush=True)


if __name__ == "__main__":
    query_check()
    family()
    terminal_check()
    partitions()
    distribution()
    exhaustive_table()
