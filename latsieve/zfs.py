from .reduction import lll3w, lll3_ints, gauss_ints


def select_frame(u, v, x, hnf, L0, L1, L2, c_line, c_plane):
    q = c_line * L0 * L1 * L2
    lu = q * (abs(u[0]) / L0 + abs(u[1]) / L1 + abs(u[2]) / L2)
    lv = q * (abs(v[0]) / L0 + abs(v[1]) / L1 + abs(v[2]) / L2)
    lx = q * (abs(x[0]) / L0 + abs(x[1]) / L1 + abs(x[2]) / L2)
    n = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
    puv = c_plane * (L0 * abs(n[0]) + L1 * abs(n[1]) + L2 * abs(n[2]))
    n = (u[1] * x[2] - u[2] * x[1], u[2] * x[0] - u[0] * x[2], u[0] * x[1] - u[1] * x[0])
    pux = c_plane * (L0 * abs(n[0]) + L1 * abs(n[1]) + L2 * abs(n[2]))
    n = (v[1] * x[2] - v[2] * x[1], v[2] * x[0] - v[0] * x[2], v[0] * x[1] - v[1] * x[0])
    pvx = c_plane * (L0 * abs(n[0]) + L1 * abs(n[1]) + L2 * abs(n[2]))
    cands = ((lu + puv, u, v, x), (lu + pux, u, x, v), (lv + puv, v, u, x), (lv + pvx, v, x, u), (lx + pux, x, u, v), (lx + pvx, x, v, u))
    best = min(cands, key=lambda t: t[0])
    if hnf is not None:
        a, c, e = hnf
        nh = (a[1] * c[2] - a[2] * c[1], a[2] * c[0] - a[0] * c[2], a[0] * c[1] - a[1] * c[0])
        ch = q * (abs(a[0]) / L0 + abs(a[1]) / L1 + abs(a[2]) / L2) + c_plane * (L0 * abs(nh[0]) + L1 * abs(nh[1]) + L2 * abs(nh[2]))
        if ch < best[0]:
            return a, c, e
    return best[1], best[2], best[3]


def zfs_frame(r, l0, l1, box, c_line=4.0, c_coord=1.5, c_plane=25.0, c_spread=1.0, c_row=0.4):
    lo0, lo1, lo2 = box.lo
    hi0, hi1, hi2 = box.hi
    L0 = hi0 - lo0 or 1
    L1 = hi1 - lo1 or 1
    L2 = hi2 - lo2 or 1
    u0, u1, u2, v0, v1, v2, x0, x1, x2 = lll3_ints(r, 0, 0, l0, 1, 0, l1, 0, 1, 1.0 / (L0 * L0), 1.0 / (L1 * L1), 1.0 / (L2 * L2))
    V = L0 * L1 * L2
    if V < r:
        a0, a1, a2, c0, c1, c2, e0, e1, e2 = u0, u1, u2, v0, v1, v2, x0, x1, x2
    else:
        q = c_line + c_coord * ((u0 != 0) + (u1 != 0) + (u2 != 0))
        lu = V * (abs(u0) / L0 + abs(u1) / L1 + abs(u2) / L2) * q + (c_spread + c_row * (u1 != 0 or u2 != 0)) * V
        q = c_line + c_coord * ((v0 != 0) + (v1 != 0) + (v2 != 0))
        lv = V * (abs(v0) / L0 + abs(v1) / L1 + abs(v2) / L2) * q + (c_spread + c_row * (v1 != 0 or v2 != 0)) * V
        q = c_line + c_coord * ((x0 != 0) + (x1 != 0) + (x2 != 0))
        lx = V * (abs(x0) / L0 + abs(x1) / L1 + abs(x2) / L2) * q + (c_spread + c_row * (x1 != 0 or x2 != 0)) * V
        puv = c_plane * (L0 * abs(u1 * v2 - u2 * v1) + L1 * abs(u2 * v0 - u0 * v2) + L2 * abs(u0 * v1 - u1 * v0))
        pux = c_plane * (L0 * abs(u1 * x2 - u2 * x1) + L1 * abs(u2 * x0 - u0 * x2) + L2 * abs(u0 * x1 - u1 * x0))
        pvx = c_plane * (L0 * abs(v1 * x2 - v2 * x1) + L1 * abs(v2 * x0 - v0 * x2) + L2 * abs(v0 * x1 - v1 * x0))
        sel = 0
        best = lu + puv
        if lu + pux < best:
            best = lu + pux
            sel = 1
        if lv + puv < best:
            best = lv + puv
            sel = 2
        if lv + pvx < best:
            best = lv + pvx
            sel = 3
        if lx + pux < best:
            best = lx + pux
            sel = 4
        if lx + pvx < best:
            best = lx + pvx
            sel = 5
        q = (L1 * L2 * (c_line + c_coord) + c_plane * L2) * r
        if q < best:
            best = q
            sel = 6
        if (c_spread + c_row) * V >= c_plane * L2 * r:
            g0, g1, h0, h1 = gauss_ints(r, 0, l0 % r, 1, 1.0 / (L0 * L0), 1.0 / (L1 * L1))
            if V * (abs(g0) / L0 + abs(g1) / L1) * (c_line + c_coord * ((g0 != 0) + (g1 != 0))) + c_row * (g1 != 0) * V + c_plane * L2 * r < best:
                sel = 7
        del lu, lv, lx, puv, pux, pvx, best, q
        if sel == 0:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = u0, u1, u2, v0, v1, v2, x0, x1, x2
        elif sel == 1:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = u0, u1, u2, x0, x1, x2, v0, v1, v2
        elif sel == 2:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = v0, v1, v2, u0, u1, u2, x0, x1, x2
        elif sel == 3:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = v0, v1, v2, x0, x1, x2, u0, u1, u2
        elif sel == 4:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = x0, x1, x2, u0, u1, u2, v0, v1, v2
        elif sel == 5:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = x0, x1, x2, v0, v1, v2, u0, u1, u2
        elif sel == 6:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = r, 0, 0, l0, 1, 0, l1, 0, 1
        else:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = g0, g1, 0, h0, h1, 0, l1, 0, 1
    return (a0, a1, a2), (c0, c1, c2), (e0, e1, e2)


def zfs(r, l0, l1, box, S, w=1, c_line=4.0, c_coord=1.5, c_plane=25.0, c_spread=1.0, c_row=0.4):
    lo0, lo1, lo2 = box.lo
    hi0, hi1, hi2 = box.hi
    L0 = hi0 - lo0 or 1
    L1 = hi1 - lo1 or 1
    L2 = hi2 - lo2 or 1
    u0, u1, u2, v0, v1, v2, x0, x1, x2 = lll3_ints(r, 0, 0, l0, 1, 0, l1, 0, 1, 1.0 / (L0 * L0), 1.0 / (L1 * L1), 1.0 / (L2 * L2))
    V = L0 * L1 * L2
    if V < r:
        a0, a1, a2, c0, c1, c2, e0, e1, e2 = u0, u1, u2, v0, v1, v2, x0, x1, x2
    else:
        q = c_line + c_coord * ((u0 != 0) + (u1 != 0) + (u2 != 0))
        lu = V * (abs(u0) / L0 + abs(u1) / L1 + abs(u2) / L2) * q + (c_spread + c_row * (u1 != 0 or u2 != 0)) * V
        q = c_line + c_coord * ((v0 != 0) + (v1 != 0) + (v2 != 0))
        lv = V * (abs(v0) / L0 + abs(v1) / L1 + abs(v2) / L2) * q + (c_spread + c_row * (v1 != 0 or v2 != 0)) * V
        q = c_line + c_coord * ((x0 != 0) + (x1 != 0) + (x2 != 0))
        lx = V * (abs(x0) / L0 + abs(x1) / L1 + abs(x2) / L2) * q + (c_spread + c_row * (x1 != 0 or x2 != 0)) * V
        puv = c_plane * (L0 * abs(u1 * v2 - u2 * v1) + L1 * abs(u2 * v0 - u0 * v2) + L2 * abs(u0 * v1 - u1 * v0))
        pux = c_plane * (L0 * abs(u1 * x2 - u2 * x1) + L1 * abs(u2 * x0 - u0 * x2) + L2 * abs(u0 * x1 - u1 * x0))
        pvx = c_plane * (L0 * abs(v1 * x2 - v2 * x1) + L1 * abs(v2 * x0 - v0 * x2) + L2 * abs(v0 * x1 - v1 * x0))
        sel = 0
        best = lu + puv
        if lu + pux < best:
            best = lu + pux
            sel = 1
        if lv + puv < best:
            best = lv + puv
            sel = 2
        if lv + pvx < best:
            best = lv + pvx
            sel = 3
        if lx + pux < best:
            best = lx + pux
            sel = 4
        if lx + pvx < best:
            best = lx + pvx
            sel = 5
        q = (L1 * L2 * (c_line + c_coord) + c_plane * L2) * r
        if q < best:
            best = q
            sel = 6
        if (c_spread + c_row) * V >= c_plane * L2 * r:
            g0, g1, h0, h1 = gauss_ints(r, 0, l0 % r, 1, 1.0 / (L0 * L0), 1.0 / (L1 * L1))
            if V * (abs(g0) / L0 + abs(g1) / L1) * (c_line + c_coord * ((g0 != 0) + (g1 != 0))) + c_row * (g1 != 0) * V + c_plane * L2 * r < best:
                sel = 7
        del lu, lv, lx, puv, pux, pvx, best, q
        if sel == 0:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = u0, u1, u2, v0, v1, v2, x0, x1, x2
        elif sel == 1:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = u0, u1, u2, x0, x1, x2, v0, v1, v2
        elif sel == 2:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = v0, v1, v2, u0, u1, u2, x0, x1, x2
        elif sel == 3:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = v0, v1, v2, x0, x1, x2, u0, u1, u2
        elif sel == 4:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = x0, x1, x2, u0, u1, u2, v0, v1, v2
        elif sel == 5:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = x0, x1, x2, v0, v1, v2, u0, u1, u2
        elif sel == 6:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = r, 0, 0, l0, 1, 0, l1, 0, 1
        else:
            a0, a1, a2, c0, c1, c2, e0, e1, e2 = g0, g1, 0, h0, h1, 0, l1, 0, 1
    del u0, u1, u2, v0, v1, v2, x0, x1, x2, V, L0, L1, L2
    n0 = a1 * c2 - a2 * c1
    n1 = a2 * c0 - a0 * c2
    n2 = a0 * c1 - a1 * c0
    D = n0 * e0 + n1 * e1 + n2 * e2
    if D < 0:
        e0, e1, e2, D = -e0, -e1, -e2, -D
    klo = -((-((n0 * lo0 if n0 > 0 else n0 * hi0) + (n1 * lo1 if n1 > 0 else n1 * hi1) + (n2 * lo2 if n2 > 0 else n2 * hi2))) // D)
    khi = ((n0 * hi0 if n0 > 0 else n0 * lo0) + (n1 * hi1 if n1 > 0 else n1 * lo1) + (n2 * hi2 if n2 > 0 else n2 * lo2)) // D
    if klo > khi:
        return 0
    m0 = e1 * a2 - e2 * a1
    m1 = e2 * a0 - e0 * a2
    m2 = e0 * a1 - e1 * a0
    f0 = f1 = f2 = t0 = t1 = t2 = g0 = g1 = g2 = h0 = h1 = h2 = 0
    if n0:
        s1 = m0 * n1 - n0 * m1
        s2 = m0 * n2 - n0 * m2
        g0 = (s1 * lo1 if s1 > 0 else s1 * hi1) + (s2 * lo2 if s2 > 0 else s2 * hi2)
        h0 = (s1 * hi1 if s1 > 0 else s1 * lo1) + (s2 * hi2 if s2 > 0 else s2 * lo2)
        if n0 > 0:
            f0, t0, g0, h0 = n0 * D, -m0 * D, -h0, -g0
        else:
            f0, t0 = -n0 * D, m0 * D
    if n1:
        s1 = m1 * n0 - n1 * m0
        s2 = m1 * n2 - n1 * m2
        g1 = (s1 * lo0 if s1 > 0 else s1 * hi0) + (s2 * lo2 if s2 > 0 else s2 * hi2)
        h1 = (s1 * hi0 if s1 > 0 else s1 * lo0) + (s2 * hi2 if s2 > 0 else s2 * lo2)
        if n1 > 0:
            f1, t1, g1, h1 = n1 * D, -m1 * D, -h1, -g1
        else:
            f1, t1 = -n1 * D, m1 * D
    if n2:
        s1 = m2 * n0 - n2 * m0
        s2 = m2 * n1 - n2 * m1
        g2 = (s1 * lo0 if s1 > 0 else s1 * hi0) + (s2 * lo1 if s2 > 0 else s2 * hi1)
        h2 = (s1 * hi0 if s1 > 0 else s1 * lo0) + (s2 * hi1 if s2 > 0 else s2 * lo1)
        if n2 > 0:
            f2, t2, g2, h2 = n2 * D, -m2 * D, -h2, -g2
        else:
            f2, t2 = -n2 * D, m2 * D
    I0, I1 = box.I[0], box.I[1]
    offa = a0 + I0 * (a1 + I1 * a2)
    offc = c0 + I0 * (c1 + I1 * c2)
    base0 = -(lo0 + I0 * (lo1 + I1 * lo2))
    del n0, n1, n2, m0, m1, m2, D, s1, s2
    for k in range(klo, khi + 1):
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
        if jl > jh:
            continue
        px = jl * c0 + k * e0
        py = jl * c1 + k * e1
        pz = jl * c2 + k * e2
        idx = base0 + px + I0 * (py + I1 * pz)
        for _ in range(jl, jh + 1):
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
                t = idx + sl * offa
                for q in range(t, t + (sh - sl + 1) * offa, offa):
                    S[q] += w
            px += c0
            py += c1
            pz += c2
            idx += offc
    return 0


def zonotope_fibers(a, c, e, lo, hi):
    n = (a[1] * c[2] - a[2] * c[1], a[2] * c[0] - a[0] * c[2], a[0] * c[1] - a[1] * c[0])
    D = n[0] * e[0] + n[1] * e[1] + n[2] * e[2]
    if D < 0:
        e = (-e[0], -e[1], -e[2])
        D = -D
    m = (e[1] * a[2] - e[2] * a[1], e[2] * a[0] - e[0] * a[2], e[0] * a[1] - e[1] * a[0])

    def lo_hi(v):
        s = sum(v[i] * lo[i] if v[i] > 0 else v[i] * hi[i] for i in range(3))
        t = sum(v[i] * hi[i] if v[i] > 0 else v[i] * lo[i] for i in range(3))
        return s, t

    nm, nM = lo_hi(n)
    klo = -((-nm) // D)
    khi = nM // D
    slabs = []
    for i in range(3):
        if n[i] == 0:
            continue
        s = tuple(m[i] * n[j] - n[i] * m[j] for j in range(3))
        g, h = lo_hi(s)
        if n[i] > 0:
            slabs.append((n[i] * D, -m[i] * D, -h, -g))
        else:
            slabs.append((-n[i] * D, m[i] * D, g, h))
    return e, klo, khi, slabs


def flag_points(basis, box, offset=(0, 0, 0), frame=None):
    lo = tuple(box.lo[i] - offset[i] for i in range(3))
    hi = tuple(box.hi[i] - offset[i] for i in range(3))
    if frame is None:
        L = [max(hi[i] - lo[i], 1) for i in range(3)]
        u, v, x = lll3w(basis, 1.0 / L[0] ** 2, 1.0 / L[1] ** 2, 1.0 / L[2] ** 2)
        a, c, e = select_frame(u, v, x, None, L[0], L[1], L[2], 13.0, 25.0)
    else:
        a, c, e = frame
    e, klo, khi, slabs = zonotope_fibers(a, c, e, lo, hi)
    out = []
    planes = 0
    lines = 0
    for k in range(klo, khi + 1):
        planes += 1
        jl = -(1 << 62)
        jh = 1 << 62
        for f, tc, g, h in slabs:
            t = tc * k
            jl = max(jl, -((t - g) // f))
            jh = min(jh, (h - t) // f)
        for j in range(jl, jh + 1):
            lines += 1
            p = [j * c[i] + k * e[i] for i in range(3)]
            sl = -(1 << 62)
            sh = 1 << 62
            for i in range(3):
                if a[i] > 0:
                    sl = max(sl, -((p[i] - lo[i]) // a[i]))
                    sh = min(sh, (hi[i] - p[i]) // a[i])
                elif a[i] < 0:
                    sl = max(sl, -((hi[i] - p[i]) // (-a[i])))
                    sh = min(sh, (p[i] - lo[i]) // (-a[i]))
            for s in range(sl, sh + 1):
                out.append(tuple(p[i] + s * a[i] + offset[i] for i in range(3)))
    return out, planes, lines
