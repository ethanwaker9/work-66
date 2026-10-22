from .reduction import fk_reduce, gauss2
from .line_sieve import line_sieve


class PlaneWalker:
    def __init__(self, r, l0, l1, box):
        I0, I1 = box.I[0], box.I[1]
        self.box = box
        (al, be), (ga, de) = fk_reduce(r, l0, I0)
        self.u = (al, be)
        self.v = (ga, de)
        wx = 1.0 / (I0 * I0)
        wy = 1.0 / (I1 * I1)
        e0, e1 = gauss2((r, 0), (l0 % r, 1), wx, wy)
        tx = l1 % r
        det = e0[0] * e1[1] - e0[1] * e1[0]
        a = (tx * e1[1]) / det
        b = (-e0[1] * tx) / det
        fa = int(a // 1)
        fb = int(b // 1)
        if (a - fa) + (b - fb) <= 1:
            tri = ((fa, fb), (fa + 1, fb), (fa, fb + 1))
        else:
            tri = ((fa + 1, fb + 1), (fa + 1, fb), (fa, fb + 1))
        self.zvecs = tuple((tx - ca * e0[0] - cb * e1[0], -ca * e0[1] - cb * e1[1]) for ca, cb in tri)
        self.big = (al, be) if -al >= ga else (ga, de)

    def next_start(self, sx, sy):
        lo, hi = self.box.lo, self.box.hi
        A0, A1, B0, B1 = lo[0], hi[0], lo[1], hi[1]
        for dx, dy in self.zvecs:
            nx = sx + dx
            ny = sy + dy
            if A0 <= nx <= A1 and B0 <= ny <= B1:
                return nx, ny
        sx += self.zvecs[0][0]
        sy += self.zvecs[0][1]
        bx, by = self.big
        if bx < 0:
            bx, by = -bx, -by
        while sx > A1:
            sx -= bx
            sy -= by
        while sx < A0:
            sx += bx
            sy += by
        return sx, sy

    def plane_hits(self, sx, sy):
        lo, hi = self.box.lo, self.box.hi
        A0, A1, B0, B1 = lo[0], hi[0], lo[1], hi[1]
        al, be = self.u
        ga, de = self.v
        out = []
        x, y = sx, sy
        while y <= B1:
            if y >= B0:
                out.append((x, y))
            if x >= A0 - al:
                x += al
                y += be
            elif x + ga <= A1:
                x += ga
                y += de
            else:
                x += al + ga
                y += be + de
        x, y = sx, sy
        while True:
            if x <= A1 + al:
                x -= al
                y -= be
            elif x - ga >= A0:
                x -= ga
                y -= de
            else:
                x -= al + ga
                y -= be + de
            if y < B0:
                break
            if y <= B1:
                out.append((x, y))
        return out


def plane_sieve(r, l0, l1, box, S, w=1):
    lo0, lo1, lo2 = box.lo
    hi0, hi1, hi2 = box.hi
    I0, I1 = box.I[0], box.I[1]
    if r <= I0 or l0 % r == 0:
        return line_sieve(r, l0, l1, box, S, w)
    pw = PlaneWalker(r, l0, l1, box)
    al, be = pw.u
    ga, de = pw.v
    offu = al + I0 * be
    offv = ga + I0 * de
    offuv = offu + offv
    x0 = (l1 * lo2) % r
    sx = x0 - r * ((x0 - lo0) // r)
    sy = 0
    while sx > hi0:
        if al < 0:
            sx += al
            sy += be
        else:
            sx -= ga
            sy -= de
    for z in range(lo2, hi2 + 1):
        if z > lo2:
            sx, sy = pw.next_start(sx, sy)
        x, y = sx, sy
        idx = (x - lo0) + I0 * ((y - lo1) + I1 * (z - lo2))
        while y <= hi1:
            if y >= lo1:
                S[idx] += w
            if x >= lo0 - al:
                x += al
                y += be
                idx += offu
            elif x + ga <= hi0:
                x += ga
                y += de
                idx += offv
            else:
                x += al + ga
                y += be + de
                idx += offuv
        x, y = sx, sy
        idx = (x - lo0) + I0 * ((y - lo1) + I1 * (z - lo2))
        while True:
            if x <= hi0 + al:
                x -= al
                y -= be
                idx -= offu
            elif x - ga >= lo0:
                x -= ga
                y -= de
                idx -= offv
            else:
                x -= al + ga
                y -= be + de
                idx -= offuv
            if y < lo1:
                break
            if y <= hi1:
                S[idx] += w
    return 0
