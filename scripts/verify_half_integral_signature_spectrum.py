"""Verify the finite arithmetic signature spectrum at the 5-space endpoint.

This script checks the support bound and enumerates the non-4 block-count
signatures for fixed defect d.  It does not decide the separate incidence
realisation problem.
"""

from math import floor


def support(d):
    budget = 4 + 8 * d
    limit = floor((9 + (17 + 32 * d) ** 0.5) / 2)
    return [
        size
        for size in range(1, limit + 1)
        if size not in (4, 5) and (size - 4) * (size - 5) <= budget
    ]


def signatures(d):
    budget = 4 + 8 * d
    sizes = support(d)
    result = []

    def visit(index, remaining, counts):
        if index == len(sizes):
            if remaining < 0 or remaining % 2:
                return
            uncovered = remaining // 2
            if uncovered < 0:
                return
            correction = uncovered
            for size, count in counts.items():
                correction += count * (size * (size - 1) // 2 - 10)
            result.append((uncovered, tuple(sorted(counts.items())), correction))
            return

        size = sizes[index]
        cost = (size - 4) * (size - 5)
        maximum = remaining // cost
        for count in range(maximum + 1):
            if count:
                counts[size] = count
            visit(index + 1, remaining - count * cost, counts)
            counts.pop(size, None)

    visit(0, budget, {})
    return result


def verify_zero_layer():
    rows = []
    for uncovered, counts, correction in signatures(0):
        counts = dict(counts)
        rows.append((counts.get(3, 0), counts.get(6, 0), uncovered, correction))
    expected_rows = {
        (0, 0, 2, 2),
        (0, 1, 1, 6),
        (0, 2, 0, 10),
        (1, 0, 1, -6),
        (1, 1, 0, -2),
        (2, 0, 0, -14),
    }
    assert set(rows) == expected_rows
    assert all(set(dict(counts)) <= {3, 6} for _, counts, _ in signatures(0))


def main():
    verify_zero_layer()
    print("d=0: six and only six numerical signatures verified.")
    for d in range(0, 9):
        options = signatures(d)
        print(
            f"d={d}: support={support(d)}, "
            f"non-4 signatures={len(options)}, "
            f"uncovered-edge range=[0,{2 + 4 * d}]"
        )


if __name__ == "__main__":
    main()
