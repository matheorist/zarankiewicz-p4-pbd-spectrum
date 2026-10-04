"""Verify the module families for the 5-block transition deficits."""

from math import ceil


def represent(total, maximum):
    count = ceil(total / maximum)
    parts = [1] * count
    remaining = total - count
    for index in range(count):
        add = min(maximum - 1, remaining)
        parts[index] += add
        remaining -= add
    assert remaining == 0
    return parts


def alpha_one_module(index):
    # Four triples with an intersection graph having y=index-1 edges.
    y = index - 1
    private = 12 - 2 * y
    shared = y
    excess = 4 * 2
    correction = 4 * (-3)
    local_defect = -private + shared
    delta = local_defect - 5 * (correction // 4)
    return excess, delta


def alpha_two_module(index):
    # j triples, 4-j sextuples, and y pairwise block intersections.
    j = ceil((index - 1) / 3)
    y = 1 + 3 * j - index
    assert 0 <= j <= 4 and 0 <= y <= 2
    private = 24 - 3 * j - 2 * y
    excess = 4 * 2
    correction = -3 * j + 9 * (4 - j)
    local_defect = 2 * private + y
    delta = local_defect - 5 * (correction // 4)
    return excess, delta


def verify():
    for index in range(1, 8):
        assert alpha_one_module(index) == (8, 3 * index)
    for index in range(1, 14):
        assert alpha_two_module(index) == (8, 3 * index)

    for alpha, maximum, denominator, module in (
        (1, 7, 21, alpha_one_module),
        (2, 13, 39, alpha_two_module),
    ):
        for delta in range(3, 180, 3):
            parts = represent(delta // 3, maximum)
            excess = sum(module(part)[0] for part in parts)
            realised_delta = sum(module(part)[1] for part in parts)
            assert realised_delta == delta
            assert excess == 8 * ceil(delta / denominator)
        print(
            f"alpha={alpha}: every delta in {{3,6,...,177}} verified "
            f"with denominator {denominator}"
        )


if __name__ == "__main__":
    verify()
