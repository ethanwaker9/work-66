from .common import Box, nfs_box, is_prime, next_prime, primes_between, hnf_basis, reference_points
from .zfs import zfs, flag_points
from .line_sieve import line_sieve
from .plane_sieve import plane_sieve
from .space_sieve import space_sieve
from .hakt import hakt
from .ntv import ntv
from .mcguire_robinson import mcguire_robinson
from .sphere_enum import sphere_enum

METHODS = (
    ("LS", line_sieve),
    ("PS", plane_sieve),
    ("GGV", space_sieve),
    ("HAKT", hakt),
    ("NTV", ntv),
    ("MR", mcguire_robinson),
    ("SE", sphere_enum),
    ("ZFS", zfs),
)
