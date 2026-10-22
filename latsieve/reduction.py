def lll3w(B, w0, w1, w2, delta=0.99):
    a0, a1, a2, b0, b1, b2, c0, c1, c2 = lll3_ints(B[0][0], B[0][1], B[0][2], B[1][0], B[1][1], B[1][2], B[2][0], B[2][1], B[2][2], w0, w1, w2, delta)
    return (a0, a1, a2), (b0, b1, b2), (c0, c1, c2)


def lll3_ints(a0, a1, a2, b0, b1, b2, c0, c1, c2, w0, w1, w2, delta=0.99):
    k = 1
    while True:
        A = w0 * a0 * a0 + w1 * a1 * a1 + w2 * a2 * a2
        ab = w0 * a0 * b0 + w1 * a1 * b1 + w2 * a2 * b2
        q = round(ab / A)
        if q:
            b0 -= q * a0
            b1 -= q * a1
            b2 -= q * a2
            ab -= q * A
        m = ab / A
        Bb = w0 * b0 * b0 + w1 * b1 * b1 + w2 * b2 * b2
        Bs = Bb - m * ab
        if k == 1:
            if Bs < (delta - m * m) * A:
                a0, a1, a2, b0, b1, b2 = b0, b1, b2, a0, a1, a2
                continue
            k = 2
        ac = w0 * a0 * c0 + w1 * a1 * c1 + w2 * a2 * c2
        bc = w0 * b0 * c0 + w1 * b1 * c1 + w2 * b2 * c2
        m20 = ac / A
        m21 = (bc - m20 * ab) / Bs
        q = round(m21)
        if q:
            c0 -= q * b0
            c1 -= q * b1
            c2 -= q * b2
            ac -= q * ab
            bc -= q * Bb
            m20 = ac / A
            m21 = (bc - m20 * ab) / Bs
        q = round(m20)
        if q:
            c0 -= q * a0
            c1 -= q * a1
            c2 -= q * a2
            ac -= q * A
            bc -= q * ab
            m20 = ac / A
            m21 = (bc - m20 * ab) / Bs
        Cc = w0 * c0 * c0 + w1 * c1 * c1 + w2 * c2 * c2
        Cs = Cc - m20 * ac - m21 * (bc - m20 * ab)
        if Cs < (delta - m21 * m21) * Bs:
            b0, b1, b2, c0, c1, c2 = c0, c1, c2, b0, b1, b2
            k = 1
            continue
        return a0, a1, a2, b0, b1, b2, c0, c1, c2


def lll3(B, w, delta=0.99):
    return list(lll3w(B, w[0], w[1], w[2], delta))


def fk_reduce(r, lam, I):
    a0, b0 = -r, 0
    a1, b1 = lam % r, 1
    while a1 >= I:
        k = (-a0) // a1
        a0 += k * a1
        b0 += k * b1
        if -a0 < I:
            break
        k = a1 // (-a0)
        a1 += k * a0
        b1 += k * b0
    if a1 >= I:
        k = (a1 - I) // (-a0) + 1
        a1 += k * a0
        b1 += k * b0
    elif -a0 >= I:
        k = (-a0 - I) // a1 + 1
        a0 += k * a1
        b0 += k * b1
    return (a0, b0), (a1, b1)


def gauss2(u, v, wx, wy):
    def n(p):
        return wx * p[0] * p[0] + wy * p[1] * p[1]

    def ip(p, q):
        return wx * p[0] * q[0] + wy * p[1] * q[1]

    if n(u) > n(v):
        u, v = v, u
    while True:
        q = round(ip(u, v) / n(u))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if n(v) >= n(u):
            return u, v
        u, v = v, u


def gauss_ints(u0, u1, v0, v1, wx, wy):
    nu = wx * u0 * u0 + wy * u1 * u1
    nv = wx * v0 * v0 + wy * v1 * v1
    if nu > nv:
        u0, u1, v0, v1, nu, nv = v0, v1, u0, u1, nv, nu
    while True:
        q = round((wx * u0 * v0 + wy * u1 * v1) / nu)
        v0 -= q * u0
        v1 -= q * u1
        nv = wx * v0 * v0 + wy * v1 * v1
        if nv >= nu:
            return u0, u1, v0, v1
        u0, u1, v0, v1, nu, nv = v0, v1, u0, u1, nv, nu
