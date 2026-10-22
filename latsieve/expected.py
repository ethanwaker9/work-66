from fractions import Fraction
from math import floor, ceil
from .common import cross, dot


def predicted_nodes(a, c, e, L):
    vol = L[0] * L[1] * L[2]
    det = abs(dot(cross(a, c), e))
    n = cross(a, c)
    p1 = Fraction(vol, det)
    p2 = p1 * sum(Fraction(abs(a[i]), L[i]) for i in range(3))
    p3 = Fraction(sum(L[i] * abs(n[i]) for i in range(3)), det)
    return float(p3), float(p2), float(p1)


def count_nodes(a, c, e, lo, hi, tau):
    n = cross(a, c)
    D = dot(n, e)
    if D < 0:
        e = (-e[0], -e[1], -e[2])
        D = -D
    m = cross(e, a)
    lo = [Fraction(v) - tau[i] for i, v in enumerate(lo)]
    hi = [Fraction(v) - tau[i] for i, v in enumerate(hi)]

    def ext(v):
        s = sum(v[i] * lo[i] if v[i] > 0 else v[i] * hi[i] for i in range(3))
        t = sum(v[i] * hi[i] if v[i] > 0 else v[i] * lo[i] for i in range(3))
        return s, t

    nm, nM = ext(n)
    k0 = ceil(nm / D)
    k1 = floor(nM / D)
    slabs = []
    for i in range(3):
        if n[i] == 0:
            continue
        s = tuple(m[i] * n[j] - n[i] * m[j] for j in range(3))
        g, h = ext(s)
        slabs.append((n[i], m[i], g, h))
    planes = 0
    lines = 0
    points = 0
    for k in range(k0, k1 + 1):
        planes += 1
        jl = None
        jh = None
        for ni, mi, g, h in slabs:
            t = mi * D * k
            if ni > 0:
                l = ceil((t - h) / (ni * D))
                u = floor((t - g) / (ni * D))
            else:
                l = ceil((g - t) / (-ni * D))
                u = floor((h - t) / (-ni * D))
            jl = l if jl is None or l > jl else jl
            jh = u if jh is None or u < jh else jh
        if jl is None or jl > jh:
            continue
        for j in range(jl, jh + 1):
            lines += 1
            p = [j * c[i] + k * e[i] for i in range(3)]
            sl = None
            sh = None
            ok = True
            for i in range(3):
                if a[i] > 0:
                    l = ceil((lo[i] - p[i]) / a[i])
                    u = floor((hi[i] - p[i]) / a[i])
                elif a[i] < 0:
                    l = ceil((hi[i] - p[i]) / a[i])
                    u = floor((lo[i] - p[i]) / a[i])
                else:
                    if not (lo[i] <= p[i] <= hi[i]):
                        ok = False
                    continue
                sl = l if sl is None or l > sl else sl
                sh = u if sh is None or u < sh else sh
            if ok and sl is not None and sl <= sh:
                points += sh - sl + 1
    return planes, lines, points


def _det(M):
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    det = Fraction(1)
    for i in range(n):
        piv = next((k for k in range(i, n) if A[k][i] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != i:
            A[i], A[piv] = A[piv], A[i]
            det = -det
        det *= A[i][i]
        for k in range(i + 1, n):
            f = A[k][i] / A[i][i]
            if f:
                A[k] = [A[k][j] - f * A[i][j] for j in range(n)]
    return det


def inverse_rows(B):
    t = len(B)
    M = [[Fraction(B[j][i]) for j in range(t)] for i in range(t)]
    A = [M[i] + [Fraction(int(i == j)) for j in range(t)] for i in range(t)]
    for i in range(t):
        piv = next(k for k in range(i, t) if A[k][i] != 0)
        A[i], A[piv] = A[piv], A[i]
        p = A[i][i]
        A[i] = [x / p for x in A[i]]
        for k in range(t):
            if k != i and A[k][i] != 0:
                f = A[k][i]
                A[k] = [A[k][j] - f * A[i][j] for j in range(2 * t)]
    return [row[t:] for row in A]


def cross_general(vectors, m):
    out = []
    for i in range(m):
        minor = [[v[j] for j in range(m) if j != i] for v in vectors]
        out.append((-1) ** i * _det(minor) if minor else Fraction((-1) ** i))
    return out


def zonotope_hrep(gens, z0):
    from itertools import combinations
    m = len(z0)
    slabs = []
    seen = set()
    if m == 1:
        lo = z0[0] + sum(min(0, g[0]) for g in gens)
        hi = z0[0] + sum(max(0, g[0]) for g in gens)
        return [((Fraction(1),), lo, hi)]
    for R in combinations(range(len(gens)), m - 1):
        nu = cross_general([gens[i] for i in R], m)
        if all(x == 0 for x in nu):
            continue
        key = tuple(nu)
        neg = tuple(-x for x in nu)
        if key in seen or neg in seen:
            continue
        seen.add(key)
        base = sum(nu[j] * z0[j] for j in range(m))
        vals = [sum(nu[j] * g[j] for j in range(m)) for g in gens]
        slabs.append((key, base + sum(min(0, v) for v in vals), base + sum(max(0, v) for v in vals)))
    return slabs


def general_node_counts(B, lo, hi, tau):
    from math import floor, ceil
    t = len(B)
    Dm = inverse_rows(B)
    hreps = []
    for k in range(t):
        rows = Dm[k:]
        m = t - k
        gens = [[rows[j][i] * (hi[i] - lo[i]) for j in range(m)] for i in range(t)]
        z0 = [sum(rows[j][i] * (Fraction(lo[i]) - tau[i]) for i in range(t)) for j in range(m)]
        hreps.append(zonotope_hrep(gens, z0))
    counts = [0] * t

    def rec(k, tail):
        sl = None
        sh = None
        for nu, a, b in hreps[k]:
            rest = sum(nu[j + 1] * tail[j] for j in range(len(tail)))
            if nu[0] == 0:
                if not (a <= rest <= b):
                    return
                continue
            u = (a - rest) / nu[0]
            v = (b - rest) / nu[0]
            if u > v:
                u, v = v, u
            sl = u if sl is None or u > sl else sl
            sh = v if sh is None or v < sh else sh
        if sl is None or sl > sh:
            return
        for n in range(ceil(sl), floor(sh) + 1):
            counts[k] += 1
            if k > 0:
                rec(k - 1, [n] + tail)

    rec(t - 1, [])
    return counts


def pluecker_l1(B, L, j):
    from itertools import combinations
    t = len(B)
    if j == 0:
        return Fraction(1)
    total = Fraction(0)
    for U in combinations(range(t), j):
        M = [[Fraction(B[c][i], L[i]) for c in range(j)] for i in U]
        total += abs(_det(M))
    return total


def lll_exact(B, w, delta=Fraction(99, 100)):
    b = [list(v) for v in B]
    n = len(b)

    def ip(u, v):
        return sum(Fraction(w[i]) * u[i] * v[i] for i in range(len(u)))

    def gso():
        bs = []
        mu = [[Fraction(0)] * n for _ in range(n)]
        Bn = []
        for i in range(n):
            v = [Fraction(x) for x in b[i]]
            for j in range(i):
                mu[i][j] = ip(b[i], bs[j]) / Bn[j]
                v = [v[k] - mu[i][j] * bs[j][k] for k in range(len(v))]
            bs.append(v)
            Bn.append(ip(v, v))
        return mu, Bn

    k = 1
    mu, Bn = gso()
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                b[k] = [b[k][i] - q * b[j][i] for i in range(len(b[k]))]
                mu, Bn = gso()
        if Bn[k] >= (delta - mu[k][k - 1] ** 2) * Bn[k - 1]:
            k += 1
        else:
            b[k], b[k - 1] = b[k - 1], b[k]
            mu, Bn = gso()
            k = max(k - 1, 1)
    return [tuple(v) for v in b]
