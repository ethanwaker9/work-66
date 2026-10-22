from .reduction import fk_reduce
from .line_sieve import line_sieve


def _sgn(a):
    return (a > 0) - (a < 0)


def _add(a, b, k=1):
    return (a[0] + k * b[0], a[1] + k * b[1], a[2] + k * b[2])


def _neg(a):
    return (-a[0], -a[1], -a[2])


def _is_opposite(u, v, w):
    if u[0] == 0:
        return _sgn(v[0]) != _sgn(w[0])
    g = u[1] / u[0]
    return _sgn(g * v[0] - v[1]) != _sgn(g * w[0] - w[1])


def _radiate(U):
    u, v, w = U
    if _is_opposite(u, v, w):
        if not _is_opposite(v, w, u):
            u = _neg(u)
    elif _is_opposite(v, w, u):
        v = _neg(v)
    else:
        w = _neg(w)
    return [u, v, w]


def _reduce1(U, I0, I1, cap):
    x, y = U[0], U[1]
    if x[2] > y[2] or (x[2] == y[2] and x[1] > y[1]):
        x, y = y, x
    w = U[2]
    n = 0
    while n < cap:
        while abs(w[0] + x[0]) < I0 and n < cap:
            if abs(w[1]) < abs(x[1]) or abs(w[1]) < I1:
                U[2] = w
                return
            w = _add(w, x)
            n += 1
        if abs(w[1]) < abs(y[1]) or abs(w[1]) < I1:
            U[2] = w
            return
        w = _add(w, y)
        n += 1
    U[2] = w


def _reduce2(U, I0, I1, cap):
    if _sgn(U[0][1]) != _sgn(U[2][1]):
        x, y = U[0], U[1]
    else:
        x, y = U[1], U[0]
    w = U[2]
    n = 0
    while abs(w[1]) >= I1 and n < cap:
        w = _add(w, x)
        n += 1
        while abs(w[0]) >= I0 and n < cap:
            w = _add(w, y)
            n += 1
        if abs(w[1]) < abs(U[0][1]) or abs(w[1]) < abs(U[1][1]):
            break
    U[2] = w


def _adjust(U, I0, I1):
    for _ in range(64):
        changed = False
        for a in range(3):
            for b in range(3):
                if a == b:
                    continue
                z = _add(U[a], U[b], -1)
                if abs(z[0]) < I0 and abs(z[1]) < I1 and z != (0, 0, 0):
                    if abs(U[a][2]) >= abs(U[b][2]):
                        U[a] = z
                    else:
                        U[b] = z
                    U[:] = _radiate(U)
                    changed = True
                    break
            if changed:
                break
        if not changed:
            break


def hakt_basis(r, l0, l1, I0, I1):
    (al, be), (ga, de) = fk_reduce(r, l0, I0)
    a = (-l1 * pow(l0, -1, r)) % r
    if a > 0:
        a -= r
    U = [(al, be, 0), (ga, de, 0), (0, a, 1)]
    cap = 8 * (I0 + I1) + 64
    it = 0
    while abs(U[2][1]) >= I1:
        it += 1
        if it > 64:
            return None
        U[:] = _radiate(U)
        if _sgn(U[0][1]) == _sgn(U[1][1]):
            _reduce1(U, I0, I1, cap)
        else:
            _reduce2(U, I0, I1, cap)
        if abs(U[0][1]) > abs(U[1][1]):
            U[2], U[0] = U[0], U[2]
        else:
            U[2], U[1] = U[1], U[2]
    U = [v if v[2] >= 0 else _neg(v) for v in U]
    _adjust(U, I0, I1)
    return [v if v[2] >= 0 else _neg(v) for v in U]


def hakt(r, l0, l1, box, S, w=1, cap=10 ** 6):
    lo0, lo1, lo2 = box.lo
    hi0, hi1, hi2 = box.hi
    I0, I1 = box.I[0], box.I[1]
    if r <= I0 or l0 % r == 0:
        return line_sieve(r, l0, l1, box, S, w)
    U = hakt_basis(r, l0, l1, I0, I1)
    if U is None or U[0][2] + U[1][2] + U[2][2] == 0:
        return 0
    a, b, c = U

    def stop(s, d):
        if s[0] > hi0 or s[1] > hi1 or s[2] > hi2:
            return True
        if s[0] < lo0 and d[0] <= 0 and d[2] == 0:
            return True
        if s[1] < lo1 and d[1] <= 0 and d[2] == 0:
            return True
        return False

    steps = 0
    p = (0, 0, 0)
    base = -(lo0 + I0 * (lo1 + I1 * lo2))
    while steps < cap:
        S[base + p[0] + I0 * (p[1] + I1 * p[2])] += w
        found = None
        for d in (a, b, c):
            t = _add(p, d)
            steps += 1
            if lo0 <= t[0] <= hi0 and lo1 <= t[1] <= hi1 and lo2 <= t[2] <= hi2:
                found = t
                break
        q = p
        while found is None and steps < cap:
            rr = q
            while found is None and steps < cap:
                s2 = rr
                while steps < cap:
                    s2 = _add(s2, a)
                    steps += 1
                    if lo0 <= s2[0] <= hi0 and lo1 <= s2[1] <= hi1 and lo2 <= s2[2] <= hi2:
                        found = s2
                        break
                    if stop(s2, a):
                        break
                if found is not None:
                    break
                rr = _add(rr, b)
                if stop(rr, b):
                    break
            if found is not None:
                break
            q = _add(q, c)
            if q[2] > hi2 or (c[2] == 0 and stop(q, c)):
                break
        if found is None:
            break
        p = found
    return 0
