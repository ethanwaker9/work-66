from math import sqrt, floor, ceil
from .reduction import lll3w


def sphere_enum(r, l0, l1, box, S, w=1):
    lo0, lo1, lo2 = box.lo
    hi0, hi1, hi2 = box.hi
    I0, I1 = box.I[0], box.I[1]
    h = [(box.hi[i] - box.lo[i]) / 2.0 + 0.5 for i in range(3)]
    cen = [(box.hi[i] + box.lo[i]) / 2.0 for i in range(3)]
    wt = [1.0 / (h[i] * h[i]) for i in range(3)]
    b = lll3w(((r, 0, 0), (l0, 1, 0), (l1, 0, 1)), wt[0], wt[1], wt[2])

    def ip(u, v):
        return wt[0] * u[0] * v[0] + wt[1] * u[1] * v[1] + wt[2] * u[2] * v[2]

    bs = []
    mu = [[0.0] * 3 for _ in range(3)]
    Bn = [0.0] * 3
    for i in range(3):
        v = [float(x) for x in b[i]]
        for j in range(i):
            mu[i][j] = ip(b[i], bs[j]) / Bn[j]
            v = [v[k] - mu[i][j] * bs[j][k] for k in range(3)]
        bs.append(v)
        Bn[i] = ip(v, v)
    adj = [(b[1][1] * b[2][2] - b[1][2] * b[2][1], b[1][2] * b[2][0] - b[1][0] * b[2][2], b[1][0] * b[2][1] - b[1][1] * b[2][0]),
           (b[2][1] * b[0][2] - b[2][2] * b[0][1], b[2][2] * b[0][0] - b[2][0] * b[0][2], b[2][0] * b[0][1] - b[2][1] * b[0][0]),
           (b[0][1] * b[1][2] - b[0][2] * b[1][1], b[0][2] * b[1][0] - b[0][0] * b[1][2], b[0][0] * b[1][1] - b[0][1] * b[1][0])]
    det = sum(adj[0][i] * b[0][i] for i in range(3))
    t = [sum(adj[j][i] * cen[i] for i in range(3)) / det for j in range(3)]
    R2 = 3.0
    rad2 = sqrt(R2 / Bn[2])
    base = -(lo0 + I0 * (lo1 + I1 * lo2))
    for y2 in range(ceil(t[2] - rad2), floor(t[2] + rad2) + 1):
        d2 = (y2 - t[2]) ** 2 * Bn[2]
        rem = R2 - d2
        if rem < 0:
            continue
        c1 = t[1] - mu[2][1] * (y2 - t[2])
        rad1 = sqrt(rem / Bn[1])
        for y1 in range(ceil(c1 - rad1), floor(c1 + rad1) + 1):
            rem1 = rem - (y1 - c1) ** 2 * Bn[1]
            if rem1 < 0:
                continue
            c0 = t[0] - mu[1][0] * (y1 - t[1]) - mu[2][0] * (y2 - t[2])
            rad0 = sqrt(rem1 / Bn[0])
            px = y1 * b[1][0] + y2 * b[2][0]
            py = y1 * b[1][1] + y2 * b[2][1]
            pz = y1 * b[1][2] + y2 * b[2][2]
            for y0 in range(ceil(c0 - rad0), floor(c0 + rad0) + 1):
                x = px + y0 * b[0][0]
                yy = py + y0 * b[0][1]
                z = pz + y0 * b[0][2]
                if lo0 <= x <= hi0 and lo1 <= yy <= hi1 and lo2 <= z <= hi2:
                    S[base + x + I0 * (yy + I1 * z)] += w
    return 0
