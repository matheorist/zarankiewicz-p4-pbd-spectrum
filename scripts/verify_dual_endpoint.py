"""Verify the modules and arithmetic used at the 5-space endpoint."""

from collections import Counter
from itertools import combinations


BASE = (0, 1, 4, 6)


def sigma(m):
    return (m * (m - 1) // 2) % 5


def cyclic_blocks(n):
    return [tuple(sorted((x + a) % n for a in BASE)) for x in range(n)]


def verify_packing(blocks, point_count):
    covered = set()
    replications = Counter()
    for block in blocks:
        assert len(block) == 4 and len(set(block)) == 4
        replications.update(block)
        for pair in combinations(block, 2):
            pair = tuple(sorted(pair))
            assert pair not in covered
            covered.add(pair)
    assert all(0 <= x < point_count for block in blocks for x in block)
    return replications, covered


def dual_floor_blocks(m, n4):
    """Return the primal packing dual to disjoint 5-edges and singleton edges."""
    assert m % 4 == 0
    h = n4 - m // 4
    assert 0 <= 5 * h <= n4

    dual_edges = []
    for j in range(h):
        dual_edges.append(tuple(range(5 * j, 5 * j + 5)))
    dual_degrees = Counter(x for edge in dual_edges for x in edge)
    for x in range(n4):
        dual_edges.extend([(x,)] * (4 - dual_degrees[x]))
    assert len(dual_edges) == m

    primal = []
    for x in range(n4):
        primal.append(tuple(i for i, edge in enumerate(dual_edges) if x in edge))
    return primal


# The cyclic module is linear and 4-regular for every order at least 13.
for n4 in range(13, 301):
    blocks = cyclic_blocks(n4)
    assert len(set(blocks)) == n4
    reps, pairs = verify_packing(blocks, n4)
    assert set(reps.values()) == {4}
    assert len(pairs) == 6 * n4


# The exact constant floors in the m = 1 (mod 4) classes.
nu_offsets = {1: 15, 5: 15, 9: 16, 13: 13, 17: 16}
for residue, expected in nu_offsets.items():
    for m in range(residue or 20, 1000, 20):
        nu = next(t for t in range(13, 18) if t % 5 == sigma(m))
        assert nu == expected


# The dual module realizes every raw integer in its interval; after the
# pair-count filter, its leave is K5-divisible.
mu_offsets = {0: 0, 4: 0, 8: 1, 12: 3, 16: 1}
for m in range(40, 401, 4):
    mu = next(t for t in range(m // 4, m // 4 + 5) if t % 5 == sigma(m))
    assert mu - m // 4 == mu_offsets[m % 20]
    for n4 in range(m // 4, 5 * m // 16 + 1):
        blocks = dual_floor_blocks(m, n4)
        reps, pairs = verify_packing(blocks, m)
        assert set(reps.values()).issubset({1, 5})
        assert len(pairs) == 6 * n4
        if n4 % 5 == sigma(m):
            leave_degrees = {m - 1 - 3 * reps[x] for x in range(m)}
            assert all(d % 4 == 0 for d in leave_degrees)
            assert (m * (m - 1) // 2 - 6 * n4) % 10 == 0


# Check the two Steiner-core replication classes and block increments.
for v in range(4, 2000):
    if v % 12 in (1, 4):
        replication = (v - 1) // 3
        assert replication % 4 == (0 if v % 12 == 1 else 1)
        blocks = v * (v - 1) // 12
        next_blocks = (v + 12) * (v + 11) // 12
        assert next_blocks - blocks == 2 * v + 11


def residue_interval(lower, upper, residue):
    lower += (residue - lower) % 5
    upper -= (upper - residue) % 5
    return None if lower > upper else (lower, upper)


def verify_amplified_overlap(m):
    """Check representative residue-filtered core bands at delta = 49/50."""
    residue = sigma(m)
    bands = []
    if m % 4 == 1:
        nu = next(t for t in range(13, 18) if t % 5 == residue)
        bands.append(residue_interval(nu, m, residue))
        core_orders = range(13, m // 50 + 1, 12)
        for v in core_orders:
            b = v * (v - 1) // 12
            bands.append(residue_interval(b + 13, b + m - v, residue))
    else:
        mu = next(t for t in range(m // 4, m // 4 + 5) if t % 5 == residue)
        bands.append(residue_interval(mu, 5 * m // 16, residue))
        core_orders = range(4, m // 50 + 1, 12)
        for v in core_orders:
            b = v * (v - 1) // 12
            remaining = m - v
            bands.append(
                residue_interval(
                    b + remaining // 4, b + 5 * remaining // 16, residue
                )
            )
    bands = sorted(band for band in bands if band is not None)
    reach = bands[0][1]
    for lower, upper in bands[1:]:
        assert lower <= reach + 5
        reach = max(reach, upper)


for order_residue in (0, 1, 4, 5, 8, 9, 12, 13, 16, 17):
    verify_amplified_overlap(100_000 + order_residue)


print(
    "Dual endpoint modules, floors, congruences, core increments, "
    "and amplified overlaps verified."
)
