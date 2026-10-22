class Box:
    def __init__(self, lo, hi):
        self.lo = tuple(lo)
        self.hi = tuple(hi)
        self.I = tuple(h - l + 1 for l, h in zip(self.lo, self.hi))
        self.size = self.I[0] * self.I[1] * self.I[2]

    def index(self, x):
        return (x[0] - self.lo[0]) + self.I[0] * ((x[1] - self.lo[1]) + self.I[1] * (x[2] - self.lo[2]))

    def point(self, i):
        x = i % self.I[0]
        i //= self.I[0]
        y = i % self.I[1]
        z = i // self.I[1]
        return (x + self.lo[0], y + self.lo[1], z + self.lo[2])

    def contains(self, x):
        return all(l <= v <= h for l, v, h in zip(self.lo, x, self.hi))


def nfs_box(h0, h1, h2):
    return Box((-h0, -h1, 0), (h0 - 1, h1 - 1, h2 - 1))


def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def next_prime(n):
    while not is_prime(n):
        n += 1
    return n


def primes_between(a, b):
    sieve = bytearray([1]) * (b + 1)
    sieve[0:2] = b"\x00\x00"
    i = 2
    while i * i <= b:
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
        i += 1
    return [p for p in range(max(a, 2), b + 1) if sieve[p]]


def hnf_basis(r, l0, l1):
    return ((r, 0, 0), (l0, 1, 0), (l1, 0, 1))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def det3(a, b, c):
    return dot(cross(a, b), c)


def reference_points(r, l0, l1, box):
    out = []
    lo, hi = box.lo, box.hi
    for z in range(lo[2], hi[2] + 1):
        for y in range(lo[1], hi[1] + 1):
            x0 = (l0 * y + l1 * z) % r
            x = x0 + r * ((lo[0] - x0 + r - 1) // r)
            while x <= hi[0]:
                out.append(box.index((x, y, z)))
                x += r
    return out
