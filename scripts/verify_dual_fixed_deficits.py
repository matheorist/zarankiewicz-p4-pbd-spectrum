"""Verify excess-eight certificates at the m == 1 (mod 4) endpoint."""

from collections import Counter
from itertools import combinations


MODULES = {
    1: (
        [(0, 2, 4, 6), (1, 3, 5, 7)],
        [(0, 1), (2, 3), (4, 5), (6, 7)],
    ),
    3: (
        [
            (0, 2, 18, 19), (1, 3, 18, 20), (0, 4, 8, 13),
            (0, 1, 9, 14), (1, 5, 10, 15), (0, 6, 11, 16),
            (1, 7, 12, 17), (2, 3, 4, 5, 6, 7),
            (8, 9, 10, 11, 12, 18), (13, 14, 15, 16, 17, 18),
        ],
        [(19, 20)],
    ),
    5: (
        [
            (2, 4, 9, 24), (1, 2, 5, 10), (0, 1, 6, 24),
            (0, 3, 14, 19), (1, 11, 18, 20), (0, 2, 15, 21),
            (0, 7, 16, 23), (3, 12, 17, 22), (1, 3, 8, 23),
            (2, 3, 13, 18), (4, 5, 6, 7, 8, 18),
            (9, 10, 11, 12, 13, 23), (14, 15, 16, 17, 18, 24),
            (19, 20, 21, 22, 23, 24),
        ],
        [],
    ),
    6: (
        [
            (0, 2, 10, 16), (1, 3, 5, 11), (2, 3, 6, 17),
            (0, 3, 7, 12), (1, 4, 18, 22), (0, 1, 13, 21),
            (2, 4, 8, 14), (0, 4, 15, 19), (1, 2, 9, 20),
            (3, 4, 21, 23), (5, 6, 7, 8, 9, 21),
            (10, 11, 12, 13, 14, 15), (16, 17, 18, 19, 20, 21),
        ],
        [(22, 23)],
    ),
    8: (
        [
            (1, 3, 17, 23), (0, 2, 11, 18), (1, 2, 5, 19),
            (0, 3, 16, 24), (0, 4, 5, 7), (2, 6, 8, 12),
            (3, 5, 9, 13), (1, 4, 10, 14), (4, 6, 16, 20),
            (0, 1, 6, 21), (2, 3, 4, 22), (5, 6, 15, 22),
            (7, 8, 9, 10, 16, 22), (11, 12, 13, 14, 15, 16),
            (17, 18, 19, 20, 21, 22),
        ],
        [(23, 24)],
    ),
    10: (
        [
            (0, 2, 7, 16), (0, 4, 6, 9), (1, 2, 8, 14),
            (0, 3, 8, 15), (0, 1, 5, 10), (1, 3, 6, 19),
            (1, 4, 7, 17), (2, 3, 4, 11), (2, 5, 6, 12),
            (3, 5, 7, 18), (4, 5, 8, 19), (6, 7, 8, 13),
            (9, 10, 11, 12, 13, 19),
        ],
        [(14, 15), (16, 17), (18, 19)],
    ),
    11: (
        [
            (0, 6, 8, 10), (0, 3, 9, 13), (1, 6, 7, 12),
            (1, 2, 14, 17), (0, 2, 16, 18), (2, 4, 7, 15),
            (0, 5, 7, 11), (1, 4, 9, 10), (1, 3, 8, 11),
            (2, 3, 5, 12), (4, 5, 8, 13), (3, 4, 6, 14),
            (5, 6, 9, 15), (7, 8, 9, 16),
            (10, 11, 12, 13, 14, 15, 16),
        ],
        [(17, 18)],
    ),
}


def verify(t, blocks, uncovered):
    covered = set()
    weighted = Counter()
    for block in blocks:
        weighted.update({point: len(block) - 1 for point in block})
        for pair in combinations(sorted(block), 2):
            assert pair not in covered
            covered.add(pair)
    for edge in uncovered:
        pair = tuple(sorted(edge))
        assert pair not in covered
        covered.add(pair)
        weighted.update(edge)
    assert all(value % 4 == 0 for value in weighted.values())
    correction = sum(len(block) * (len(block) - 1) // 2 - 10 for block in blocks)
    correction += len(uncovered)
    excess = sum((len(block) - 4) * (len(block) - 5) for block in blocks)
    excess += 2 * len(uncovered)
    assert correction == -4 * t
    assert excess == 8
    return max(weighted) + 1, Counter(map(len, blocks)), len(uncovered)


for target, (packing, leave) in MODULES.items():
    points, signature, uncovered_count = verify(target, packing, leave)
    print(
        f"t={target}: E=8 certificate on {points} points; "
        f"blocks={dict(signature)}, uncovered={uncovered_count}"
    )
