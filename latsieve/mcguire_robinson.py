from .reduction import lll3w


def mcguire_robinson(r, l0, l1, box, S, w=1):
    lo0, lo1, lo2 = box.lo
    hi0, hi1, hi2 = box.hi
    lo, hi = box.lo, box.hi
    I0, I1 = box.I[0], box.I[1]
    red = sorted(lll3w(((r, 0, 0), (l0, 1, 0), (l1, 0, 1)), 1.0, 1.0, 1.0), key=lambda v: v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    v1, v2, v3 = red
    n = (v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0])
    D = n[0] * v3[0] + n[1] * v3[1] + n[2] * v3[2]
    if D < 0:
        v3 = (-v3[0], -v3[1], -v3[2])
        D = -D
    nmin = sum(n[i] * lo[i] if n[i] > 0 else n[i] * hi[i] for i in range(3))
    nmax = sum(n[i] * hi[i] if n[i] > 0 else n[i] * lo[i] for i in range(3))
    kmin = -((-nmin) // D)
    kmax = nmax // D
    m = (v1[1] * n[2] - v1[2] * n[1], v1[2] * n[0] - v1[0] * n[2], v1[0] * n[1] - v1[1] * n[0])
    dm = m[0] * v2[0] + m[1] * v2[1] + m[2] * v2[2]
    off1 = v1[0] + I0 * (v1[1] + I1 * v1[2])
    base = -(lo0 + I0 * (lo1 + I1 * lo2))
    corners = [(x, y, z) for x in (lo0, hi0) for y in (lo1, hi1) for z in (lo2, hi2)]

    def inside(px, py, pz):
        return lo0 <= px <= hi0 and lo1 <= py <= hi1 and lo2 <= pz <= hi2

    def arange(px, py, pz):
        sl = -(1 << 62)
        sh = 1 << 62
        for a, p, l, h in ((v1[0], px, lo0, hi0), (v1[1], py, lo1, hi1), (v1[2], pz, lo2, hi2)):
            if a > 0:
                sl = max(sl, -((p - l) // a))
                sh = min(sh, (h - p) // a)
            elif a < 0:
                sl = max(sl, -((h - p) // (-a)))
                sh = min(sh, (p - l) // (-a))
            elif not (l <= p <= h):
                return None
        return (sl, sh) if sl <= sh else None

    def ilp(R):
        vals = [m[0] * (c[0] - R[0]) + m[1] * (c[1] - R[1]) + m[2] * (c[2] - R[2]) for c in corners]
        if dm > 0:
            bl, bh = -((-min(vals)) // dm), max(vals) // dm
        else:
            bl, bh = -((max(vals)) // (-dm)), (-min(vals)) // (-dm)
        for b in range(bl, bh + 1):
            px, py, pz = R[0] + b * v2[0], R[1] + b * v2[1], R[2] + b * v2[2]
            ar = arange(px, py, pz)
            if ar is not None:
                a = ar[0]
                return (px + a * v1[0], py + a * v1[1], pz + a * v1[2])
        return None

    def walkline(px, py, pz):
        idx = base + px + I0 * (py + I1 * pz)
        x, y, z, q = px, py, pz, idx
        while lo0 <= x <= hi0 and lo1 <= y <= hi1 and lo2 <= z <= hi2:
            S[q] += w
            x += v1[0]
            y += v1[1]
            z += v1[2]
            q += off1
        x, y, z, q = px - v1[0], py - v1[1], pz - v1[2], idx - off1
        while lo0 <= x <= hi0 and lo1 <= y <= hi1 and lo2 <= z <= hi2:
            S[q] += w
            x -= v1[0]
            y -= v1[1]
            z -= v1[2]
            q -= off1

    for k in range(kmin, kmax + 1):
        p0 = ilp((k * v3[0], k * v3[1], k * v3[2]))
        if p0 is None:
            continue
        walkline(*p0)
        for sg in (1, -1):
            qx, qy, qz = p0
            while True:
                qx += sg * v2[0]
                qy += sg * v2[1]
                qz += sg * v2[2]
                if not inside(qx, qy, qz):
                    if inside(qx + v1[0], qy + v1[1], qz + v1[2]):
                        qx, qy, qz = qx + v1[0], qy + v1[1], qz + v1[2]
                    elif inside(qx - v1[0], qy - v1[1], qz - v1[2]):
                        qx, qy, qz = qx - v1[0], qy - v1[1], qz - v1[2]
                    else:
                        break
                walkline(qx, qy, qz)
    return 0
