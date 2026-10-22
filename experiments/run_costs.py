import json
import os
import timeit
from _path import RESULTS

S = [0] * 1048576


def points(n):
    w = 1
    for q in range(5, 5 + 7 * n, 7):
        S[q] += w


def lines(n, a0, a1, a2):
    lo0, lo1, lo2 = -64, -64, 0
    hi0, hi1, hi2 = 63, 63, 63
    px, py, pz = 0, 0, 10
    sl = sh = 0
    for _ in range(n):
        if a0 > 0:
            sl = -((px - lo0) // a0)
            sh = (hi0 - px) // a0
        elif a0 < 0:
            sl = -((hi0 - px) // (-a0))
            sh = (px - lo0) // (-a0)
        else:
            sl = -(1 << 62)
            sh = 1 << 62
        if a1 > 0:
            t = -((py - lo1) // a1)
            if t > sl:
                sl = t
            t = (hi1 - py) // a1
            if t < sh:
                sh = t
        elif a1 < 0:
            t = -((hi1 - py) // (-a1))
            if t > sl:
                sl = t
            t = (py - lo1) // (-a1)
            if t < sh:
                sh = t
        if a2 > 0:
            t = -((pz - lo2) // a2)
            if t > sl:
                sl = t
            t = (hi2 - pz) // a2
            if t < sh:
                sh = t
        elif a2 < 0:
            t = -((hi2 - pz) // (-a2))
            if t > sl:
                sl = t
            t = (pz - lo2) // (-a2)
            if t < sh:
                sh = t
        if sl <= sh:
            for q in range(5, 5, 1):
                S[q] += 1
        px += 1
        py += 0
    return sl, sh


def planes(n):
    f0, t0, g0, h0 = 7919, 104729, -10 ** 9, 10 ** 9
    f1, t1, g1, h1 = 7907, -104723, -10 ** 9, 10 ** 9
    f2, t2, g2, h2 = 7901, 104711, -10 ** 9, 10 ** 9
    jl = jh = 0
    for k in range(n):
        jl = -(1 << 62)
        jh = 1 << 62
        if f0:
            t = t0 * k
            l = -((t - g0) // f0)
            h = (h0 - t) // f0
            if l > jl:
                jl = l
            if h < jh:
                jh = h
        if f1:
            t = t1 * k
            l = -((t - g1) // f1)
            h = (h1 - t) // f1
            if l > jl:
                jl = l
            if h < jh:
                jh = h
        if f2:
            t = t2 * k
            l = -((t - g2) // f2)
            h = (h2 - t) // f2
            if l > jl:
                jl = l
            if h < jh:
                jh = h
        for _ in range(0):
            pass
    return jl, jh


if __name__ == "__main__":
    out = {}
    out["point_ns"] = timeit.timeit(lambda: points(10000), number=200) / 200 / 10000 * 1e9
    for name, a in (("line_1_ns", (200, 0, 0)), ("line_2_ns", (200, 3, 0)), ("line_3_ns", (200, 3, 5))):
        out[name] = timeit.timeit(lambda: lines(10000, *a), number=50) / 50 / 10000 * 1e9
    out["line_start_ns"] = timeit.timeit(lambda: [range(5, 6, 1) for _ in range(10000)], number=50) / 50 / 10000 * 1e9
    out["plane_ns"] = timeit.timeit(lambda: planes(10000), number=50) / 50 / 10000 * 1e9
    with open(os.path.join(RESULTS, "unit_costs.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(out)
