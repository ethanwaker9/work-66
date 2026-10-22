from itertools import product
from .reduction import lll3w
from .plane_sieve import PlaneWalker
from .line_sieve import line_sieve
from .plane_sieve import plane_sieve


def space_sieve(r, l0, l1, box, S, w=1, ncomb=2, keep=6):
    lo0, lo1, lo2 = box.lo
    hi0, hi1, hi2 = box.hi
    I0, I1 = box.I[0], box.I[1]
    if r <= I0 or l0 % r == 0:
        return line_sieve(r, l0, l1, box, S, w)
    if r <= I0 * I1:
        return plane_sieve(r, l0, l1, box, S, w)
    wz = I0 * I1 / r
    red = lll3w(((r, 0, 0), (l0, 1, 0), (l1, 0, 1)), 1.0 / (I0 * I0), 1.0 / (I1 * I1), wz * wz)
    red = [v if v[2] >= 0 else (-v[0], -v[1], -v[2]) for v in red]
    cand = set()
    rng = range(-ncomb, ncomb + 1)
    for c in product(rng, rng, rng):
        if c == (0, 0, 0):
            continue
        v = tuple(c[0] * red[0][i] + c[1] * red[1][i] + c[2] * red[2][i] for i in range(3))
        if v[2] < 0:
            v = (-v[0], -v[1], -v[2])
        if abs(v[0]) < I0 and abs(v[1]) < I1:
            cand.add(v)
    flat = sorted((v for v in cand if v[2] == 0), key=lambda v: abs(v[0]) / I0 + abs(v[1]) / I1)
    lz0 = []
    for v in flat:
        if all(v[0] * u[1] - v[1] * u[0] != 0 for u in lz0):
            lz0.append(v)
        if len(lz0) == 2:
            break
    lnz = sorted((v for v in cand if v[2] > 0), key=lambda v: v[2])[:keep]
    pw = None
    base = -(lo0 + I0 * (lo1 + I1 * lo2))
    s = (0, 0, 0)
    while s[2] <= hi2:
        lc = [s]
        for v in lz0:
            for sg in (1, -1):
                px, py = s[0] + sg * v[0], s[1] + sg * v[1]
                while lo0 <= px <= hi0 and lo1 <= py <= hi1:
                    if (px, py, s[2]) not in lc:
                        lc.append((px, py, s[2]))
                    px += sg * v[0]
                    py += sg * v[1]
        for p in lc:
            S[base + p[0] + I0 * (p[1] + I1 * p[2])] += w
        st = None
        for sn in lc:
            for v in lnz:
                qx, qy, qz = sn[0] + v[0], sn[1] + v[1], sn[2] + v[2]
                if lo0 <= qx <= hi0 and lo1 <= qy <= hi1 and qz <= hi2:
                    if st is None or qz < st[2]:
                        st = (qx, qy, qz)
                    break
        if st is None:
            if pw is None:
                pw = PlaneWalker(r, l0, l1, box)
            sx, sy = lc[0][0], lc[0][1]
            z = s[2]
            while z < hi2:
                z += 1
                sx, sy = pw.next_start(sx, sy)
                hits = pw.plane_hits(sx, sy)
                if hits:
                    hb = min(hits, key=lambda h: h[1])
                    st = (hb[0], hb[1], z)
                    break
            if st is None:
                break
            lnz.append((st[0] - s[0], st[1] - s[1], st[2] - s[2]))
            lnz.sort(key=lambda u: u[2])
        s = st
    return 0
