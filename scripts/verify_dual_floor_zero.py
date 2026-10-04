"""Verify the cyclic 4-regular quadruple modules for the dual endpoint."""

from collections import Counter
from itertools import combinations


BASE = (0, 1, 4, 6)


def verify(n):
    blocks = [tuple(sorted((x + a) % n for a in BASE)) for x in range(n)]
    assert len(set(blocks)) == n
    degrees = Counter(x for block in blocks for x in block)
    assert set(degrees.values()) == {4}
    covered = set()
    for block in blocks:
        assert len(set(block)) == 4
        for pair in combinations(block, 2):
            pair = tuple(sorted(pair))
            assert pair not in covered
            covered.add(pair)
    assert len(covered) == 6 * n


for n in range(15, 301, 5):
    verify(n)

for v in range(25, 1000):
    if v % 60 in (1, 25):
        assert v % 12 == 1
        assert ((v - 1) // 3) % 4 == 0
        assert (v * (v - 1) // 12) % 5 == 0

for m in range(21, 1000):
    if m % 20 in (1, 5):
        pair_count = m * (m - 1) // 2
        assert pair_count % 10 == 0
        for n4 in range(15, m + 1, 5):
            assert (pair_count - 6 * n4) % 10 == 0
            assert (m - 13) % 4 == 0

print("Dual floor-zero cyclic modules and core congruences verified through 1000.")
