# Exhaustive lattice sieving in three dimensions

This repository contains our implementation and the experimental pipeline of our research
*Transition Vectors, Zonotope Fibers and Exhaustive Lattice Sieving in Three Dimensions*. The sieving step of the number field sieve enumerates, for every prime ideal of norm `r`, the points of the lattice
```
Λ = { x in Z^3 : x0 ≡ λ0·x1 + λ1·x2 (mod r) }
```
inside a box `H = [m0,M0) × [m1,M1) × [m2,M2)` and adds a weight to the corresponding cells of a sieve array.
The package implements the zonotope flag sieve (ZFS) of our research together with seven existing methods,
all on the same primitives, and the tools used to study the transition vectors of the three dimensional space sieve.

## Contents

```
latsieve/            library
  common.py          boxes, primes, sieve lattices, exhaustive reference enumeration
  reduction.py       LLL reduction in dimension 3, Franke–Kleinjung reduction, Gauss reduction
  zfs.py             ZFS, generic flag enumeration, zonotope fibers
  line_sieve.py      LS, line sieve
  plane_sieve.py     PS, plane sieve of Gaudry, Grémy and Videau
  space_sieve.py     GGV, line, plane and space sieve of Gaudry, Grémy and Videau
  hakt.py            HAKT, three dimensional lattice sieve of Hayasaka, Aoki, Kobayashi and Takagi
  ntv.py             NTV, nearly transition vector sieve of Grémy (localntv framework)
  mcguire_robinson.py MR, planes spanned by successive minima of McGuire and Robinson
  sphere_enum.py     SE, enumeration of the ellipsoid containing the box (De Micheli, Gaudry and Pierrot)
  transitions.py     transition sets, terminal links
  expected.py        expected node counts and exact node counting for rational translates, in all dimensions
experiments/         scripts that produce every table and figure
results/             CSV and JSON outputs of the experiments
```

Every sieving function has the signature
```python
method(r, l0, l1, box, S, w=1)
```
where `box` is a `latsieve.Box`, `S` is a Python list of length `box.size` indexed by
`(x0-m0) + I0*((x1-m1) + I1*(x2-m2))`, and `w` is the weight added to every hit.

## Running the experiments
It requires Python 3.10 or later and matplotlib 3.7 or later, for the figures.
```python
from latsieve import nfs_box, zfs, reference_points

box = nfs_box(64, 64, 64)
S = [0] * box.size
r, l0, l1 = 1000003, 123457, 765431
zfs(r, l0, l1, box, S)
hits = sorted(i for i, v in enumerate(S) if v)
assert hits == sorted(reference_points(r, l0, l1, box))
```
All commands are run from the `experiments` directory.

```bash
cd experiments
python3 check_exactness.py
```
checks every method against the exhaustive reference on 300 random lattices and boxes and prints the fraction of points found and the number of invalid or repeated hits.
```bash
./run_all.sh
```
runs the complete pipeline. On a single core of an Apple M2 the whole pipeline takes about forty minutes; the factor base sweep over all 23000 prime ideals of norm at most 2^18 is the longest step. Timings are the median of three runs per lattice, and the peak working memory is measured with `tracemalloc` during one call, excluding the sieve array, which is allocated before the call.

