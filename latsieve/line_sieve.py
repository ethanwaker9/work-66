def line_sieve(r, l0, l1, box, S, w=1):
    lo0, lo1, lo2 = box.lo
    hi0, hi1, hi2 = box.hi
    I0, I1 = box.I[0], box.I[1]
    for z in range(lo2, hi2 + 1):
        x0 = (l0 * lo1 + l1 * z) % r
        base = I0 * I1 * (z - lo2) - lo0
        for y in range(lo1, hi1 + 1):
            x = x0 + r * ((lo0 - x0 + r - 1) // r)
            if x <= hi0:
                t = base + x
                for q in range(t, t + ((hi0 - x) // r + 1) * r, r):
                    S[q] += w
            x0 += l0
            if x0 >= r:
                x0 -= r
            base += I0
    return 0
