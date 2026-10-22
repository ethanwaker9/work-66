from .common import is_prime


def positions_order(r, alpha, beta, I, J):
    vals = []
    for a in range(I):
        for b in range(J):
            vals.append(((alpha * a + beta * b) % r, a, b))
    vals.sort()
    return vals


def hypothesis_h(r, alpha, beta, I, J, factor=1):
    for dy in range(-factor * J + 1, factor * J):
        for dx in range(-factor * I + 1, factor * I):
            if (dx, dy) != (0, 0) and (alpha * dx + beta * dy) % r == 0:
                return False
    return True


def transition_set(r, alpha, beta, I, J):
    vals = positions_order(r, alpha, beta, I, J)
    n = len(vals)
    T = set()
    for i in range(n):
        z, a, b = vals[i]
        z2, a2, b2 = vals[(i + 1) % n]
        g = (z2 - z) % r
        if g == 0:
            g = r
        T.add((a2 - a, b2 - b, g))
    return T


def transition_links(r, alpha, beta, I, J):
    vals = positions_order(r, alpha, beta, I, J)
    n = len(vals)
    links = []
    for i in range(n):
        z, a, b = vals[i]
        z2, a2, b2 = vals[(i + 1) % n]
        g = (z2 - z) % r or r
        links.append(((a, b), (a2, b2), g, z))
    return links


def terminal_vectors(r, alpha, beta, I, J):
    links = transition_links(r, alpha, beta, I, J)
    last_col = sorted(((alpha * (I - 1) + beta * b) % r) for b in range(J))
    out = set()
    for p, q, g, z in links:
        if p[0] == 0 or q[0] == 0:
            out.add((q[0] - p[0], q[1] - p[1], g))
            continue
        start = (z - alpha) % r
        hit = False
        for v in last_col:
            d = (v - start) % r
            if 0 < d < g:
                hit = True
                break
        if hit:
            out.add((q[0] - p[0], q[1] - p[1], g))
    return out


def diagonal_alpha_beta(I, J):
    M = 2 * (I + J)
    if I <= J:
        return M, 1 - M
    return 1 - M, M


def diagonal_family(I, J, r):
    alpha, beta = diagonal_alpha_beta(I, J)
    return alpha % r, beta % r


def lattice_from_alpha_beta(r, alpha, beta):
    C = (-pow(alpha, -1, r)) % r
    B = (beta * pow(alpha, -1, r)) % r
    l0 = (-B) % r
    l1 = (-C) % r
    return l0, l1


def alpha_beta_from_hnf(r, l0, l1):
    inv = pow(l1, -1, r)
    alpha = inv % r
    beta = (-l0 * inv) % r
    return alpha, beta


def max_transition_count(I, J, r):
    best = 0
    arg = None
    for alpha in range(1, r):
        for beta in range(0, r):
            if not hypothesis_h(r, alpha, beta, I, J):
                continue
            k = len(transition_set(r, alpha, beta, I, J))
            if k > best:
                best = k
                arg = (alpha, beta)
    return best, arg


def _strip(basis, frame, a0, I, b0, J, zlo, zhi):
    from .common import Box
    from .zfs import flag_points
    pts, _, _ = flag_points(basis, Box((a0, b0, zlo), (a0 + I - 1, b0 + J - 1, zhi)), frame=frame)
    return pts


def _neighbor(basis, frame, a0, I, b0, J, z, up, queries):
    G = 1
    while True:
        queries[0] += 1
        if up:
            pts = _strip(basis, frame, a0, I, b0, J, z + 1, z + G)
            if pts:
                return min(pts, key=lambda p: p[2])
        else:
            pts = _strip(basis, frame, a0, I, b0, J, z - G, z - 1)
            if pts:
                return max(pts, key=lambda p: p[2])
        G *= 2


def transition_set_by_queries(r, l0, l1, I, J, a0=0, b0=0):
    from .reduction import lll3w
    from .zfs import select_frame
    basis = ((r, 0, 0), (l0, 1, 0), (l1, 0, 1))
    swap = J > I
    h = max(1, r // (I * J))
    u, v, x = lll3w(basis, 1.0 / I ** 2, 1.0 / J ** 2, 1.0 / h ** 2)
    frame = select_frame(u, v, x, None, I, J, h, 13.0, 25.0)
    inv = pow(l1, -1, r)

    def zof(px, py):
        return ((px - l0 * py) * inv) % r

    queries = [0]
    links = set()
    if not swap:
        first = [(a0, b0 + k) for k in range(J)]
        last = [(a0 + I - 1, b0 + k) for k in range(J)]
        shift = (1, 0)
    else:
        first = [(a0 + k, b0) for k in range(I)]
        last = [(a0 + k, b0 + J - 1) for k in range(I)]
        shift = (0, 1)
    for p in first:
        z = zof(*p)
        q = _neighbor(basis, frame, a0, I, b0, J, z, True, queries)
        links.add((q[0] - p[0], q[1] - p[1], q[2] - z))
        u0 = _neighbor(basis, frame, a0, I, b0, J, z, False, queries)
        links.add((p[0] - u0[0], p[1] - u0[1], z - u0[2]))
    for s in last:
        zs = zof(s[0] + shift[0], s[1] + shift[1])
        q = _neighbor(basis, frame, a0, I, b0, J, zs - 1, True, queries)
        if q[2] == zs:
            continue
        u0 = _neighbor(basis, frame, a0, I, b0, J, zs, False, queries)
        links.add((q[0] - u0[0], q[1] - u0[1], q[2] - u0[2]))
    return links, queries[0]
