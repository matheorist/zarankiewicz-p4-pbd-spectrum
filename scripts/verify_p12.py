"""Verify the 3-regular (20,5,1)-packing P_12 used in Theorem 3."""

from collections import Counter
from itertools import combinations


BLOCKS = [
    (11, 4, 17, 12, 14),
    (5, 3, 0, 2, 8),
    (3, 15, 4, 7, 13),
    (12, 13, 8, 10, 6),
    (3, 14, 18, 6, 16),
    (4, 19, 9, 6, 2),
    (2, 13, 17, 18, 1),
    (14, 5, 7, 10, 19),
    (9, 12, 1, 15, 5),
    (17, 9, 10, 16, 0),
    (11, 0, 19, 15, 18),
    (16, 8, 1, 11, 7),
]


def main():
    assert len(BLOCKS) == 12
    assert all(len(block) == 5 and len(set(block)) == 5 for block in BLOCKS)
    assert set().union(*map(set, BLOCKS)) == set(range(20))

    degrees = Counter(point for block in BLOCKS for point in block)
    assert {degrees[point] for point in range(20)} == {3}

    covered_pairs = set()
    for block in BLOCKS:
        for pair in combinations(sorted(block), 2):
            assert pair not in covered_pairs
            covered_pairs.add(pair)

    assert len(covered_pairs) == 120
    print(
        "P_12 verified: 12 quintuples on 20 points, "
        "all degrees 3, and all 120 covered pairs are distinct."
    )


if __name__ == "__main__":
    main()
