"""Verify the endpoint-deficit modules and the t=3 E=8 obstruction."""

from itertools import combinations

from search_defect8_incidence import solve
from search_defect8_modules import count_options, verify


MODULES = {
    # target: (expected excess, blocks, uncovered pairs)
    3: (
        16,
        [
            (0, 3, 18), (5, 8, 13), (1, 9, 13),
            (0, 2, 6), (4, 8, 14), (0, 9, 14),
            (1, 2, 7, 10, 15),
            (3, 4, 7, 11, 16),
            (5, 6, 7, 12, 17),
            (8, 9, 10, 11, 12, 18),
            (13, 14, 15, 16, 17, 18),
        ],
        [],
    ),
    6: (
        8,
        [
            (2, 9, 13), (7, 10, 17), (3, 4, 10), (0, 1, 14),
            (0, 5, 8, 10, 18),
            (1, 5, 11, 15, 19),
            (2, 6, 8, 15, 20),
            (3, 6, 11, 16, 18),
            (4, 5, 6, 12, 21),
            (7, 8, 12, 16, 19),
            (9, 10, 11, 12, 20),
            (13, 14, 15, 16, 21),
            (17, 18, 19, 20, 21),
        ],
        [],
    ),
    9: (
        8,
        [
            (1, 2, 3), (0, 5, 7),
            (0, 4, 9, 10, 16),
            (1, 4, 8, 11, 17),
            (2, 6, 8, 12, 18),
            (3, 4, 6, 13, 19),
            (5, 6, 9, 14, 20),
            (7, 8, 9, 15, 21),
            (10, 11, 12, 13, 14, 15),
            (16, 17, 18, 19, 20, 21),
        ],
        [],
    ),
}


DUAL_MODULES = {
    5: [
        (0, 1, 2),
        (2, 4, 10, 16), (0, 5, 11, 17), (1, 6, 12, 18),
        (0, 3, 7, 13), (1, 3, 14, 19), (3, 8, 15, 20),
        (2, 3, 9, 21),
        (4, 5, 6, 7, 8, 9),
        (10, 11, 12, 13, 14, 15),
        (16, 17, 18, 19, 20, 21),
    ],
    10: [
        (0, 4, 8), (1, 5, 6), (2, 3, 7),
        (0, 3, 6, 9), (1, 4, 7, 10), (2, 5, 8, 11),
        (0, 1, 2, 12), (3, 4, 5, 13), (6, 7, 8, 14),
        (9, 10, 11, 12, 13, 14),
    ],
}


def verify_modules():
    for target, (excess, blocks, uncovered) in MODULES.items():
        order = 1 + max(point for block in blocks for point in block)
        verify(order, target, blocks, uncovered, excess)
        print(f"t={target}: E={excess} module verified on {order} points")


def verify_dual_modules():
    for target, blocks in DUAL_MODULES.items():
        seen = set()
        local = {}
        excess = 0
        correction = 0
        for block in blocks:
            size = len(block)
            excess += (size - 4) * (size - 5)
            correction += size * (size - 1) // 2 - 10
            for point in block:
                local[point] = local.get(point, 0) + size - 1
            for pair in combinations(block, 2):
                assert pair not in seen
                seen.add(pair)
        assert excess == 8
        assert correction == -4 * target
        assert all(value % 4 == 0 for value in local.values())
        print(
            f"dual t={target}: E=8 module verified on {len(local)} points"
        )


def verify_t3_obstruction():
    options = list(count_options(3))
    assert len(options) == 14
    for index, (counts, uncovered_count) in enumerate(options):
        module, status, _ = solve(3, counts, uncovered_count, time_limit=60)
        assert module is None and status == "INFEASIBLE", (index, counts, status)
    print("t=3: all 14 E=8 incidence signatures are infeasible")


if __name__ == "__main__":
    verify_modules()
    verify_dual_modules()
    verify_t3_obstruction()
