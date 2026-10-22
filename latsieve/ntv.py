import sys
from itertools import product
from math import floor
from .reduction import lll3w, gauss2
from .line_sieve import line_sieve


def ntv(r, l0, l1, box, S, w=1, ncomb=2):
    lo, hi, I = box.lo, box.hi, box.I
    if r <= I[0] or l0 % r == 0:
        return line_sieve(r, l0, l1, box, S, w)
    level = 1 if r < I[0] * I[1] else 2
    if level == 2:
        sz = max(1.0, r / (I[0] * I[1]))
        red = lll3w(((r, 0, 0), (l0, 1, 0), (l1, 0, 1)), 1.0 / (I[0] * I[0]), 1.0 / (I[1] * I[1]), 1.0 / (sz * sz))
    else:
        e0, e1 = gauss2((r, 0), (l0 % r, 1), 1.0 / (I[0] * I[0]), 1.0 / (max(1.0, r / I[0]) ** 2))
        tx = l1 % r
        det = e0[0] * e1[1] - e0[1] * e1[0]
        ca = round(tx * e1[1] / det)
        cb = round(-e0[1] * tx / det)
        red = ((e0[0], e0[1], 0), (e1[0], e1[1], 0), (tx - ca * e0[0] - cb * e1[0], -ca * e0[1] - cb * e1[1], 1))
    rng = range(-ncomb, ncomb + 1)
    T = [[], [], []]
    SS = [[], [], []]
    seen = set()
    for c in product(rng, rng, rng):
        if c == (0, 0, 0):
            continue
        v = tuple(c[0] * red[0][i] + c[1] * red[1][i] + c[2] * red[2][i] for i in range(3))
        if v in seen:
            continue
        seen.add(v)
        if v[2] > 0:
            k = 2
        elif v[2] == 0 and v[1] > 0:
            k = 1
        elif v[2] == 0 and v[1] == 0 and v[0] > 0:
            k = 0
        else:
            continue
        SS[k].append(v)
        if all(abs(v[i]) < I[i] for i in range(k)):
            T[k].append(v)
    for k in range(3):
        T[k].sort(key=lambda v: v[k])
        SS[k].sort(key=lambda v: v[k])
    q0, q1 = gauss2((r, 0), (l0 % r, 1), 1.0 / (I[0] * I[0]), 1.0 / (I[1] * I[1]))
    det2 = q0[0] * q1[1] - q0[1] * q1[0]
    A0, A1, B0, B1, Z0, Z1 = lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]
    base = -(A0 + I[0] * (B0 + I[1] * Z0))

    def cva(t):
        fa = floor((t[0] * q1[1] - t[1] * q1[0]) / det2)
        fb = floor((q0[0] * t[1] - q0[1] * t[0]) / det2)
        return [((fa + da) * q0[0] + (fb + db) * q1[0], (fa + da) * q0[1] + (fb + db) * q1[1], 0) for da in (0, 1, -1, 2) for db in (0, 1, -1, 2)]

    def inside(p, k):
        if k >= 0 and not (A0 <= p[0] <= A1):
            return False
        if k >= 1 and not (B0 <= p[1] <= B1):
            return False
        if k >= 2 and not (Z0 <= p[2] <= Z1):
            return False
        return True

    def fallback(c, k):
        cur = c
        for _ in range(64):
            if cur[k] > hi[k]:
                break
            L = []
            for v in SS[k]:
                t = (cur[0] + v[0], cur[1] + v[1], cur[2] + v[2])
                for e in cva(t):
                    q = (t[0] - e[0], t[1] - e[1], t[2] - e[2])
                    if inside(q, 2) and q[k] > c[k]:
                        return q
                    if q[k] > cur[k]:
                        L.append(q)
            if not L:
                return None
            cur = min(L, key=lambda q: q[k])
        return None

    def walk(c, k, sg):
        while True:
            nxt = None
            for v in T[k]:
                q = (c[0] + sg * v[0], c[1] + sg * v[1], c[2] + sg * v[2])
                if inside(q, k - 1):
                    nxt = q
                    break
            if nxt is None and (k > level or k == 2) and sg == 1:
                e = fallback(c, k)
                if e is not None:
                    d = (e[0] - c[0], e[1] - c[1], e[2] - c[2])
                    if d not in T[k]:
                        T[k].append(d)
                        T[k].sort(key=lambda v: v[k])
                    nxt = e
            if nxt is None:
                return
            if (sg == 1 and nxt[k] > hi[k]) or (sg == -1 and nxt[k] < lo[k]):
                return
            c = nxt
            enum(c, k - 1)

    def enum(c, k):
        if k == 0:
            if inside(c, 2):
                S[base + c[0] + I[0] * (c[1] + I[1] * c[2])] += w
            walk(c, 0, 1)
            walk(c, 0, -1)
            return
        enum(c, k - 1)
        walk(c, k, 1)
        walk(c, k, -1)

    limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(limit, 100000))
    enum((0, 0, 0), 2)
    sys.setrecursionlimit(limit)
    return 0
