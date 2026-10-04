"""Search for cyclic 5-regular linear triple systems on 18, 21, and 33 points."""

from collections import Counter
from itertools import combinations


def candidate(order, orbit_offsets, class_offsets):
    g = order // 3
    a, b = orbit_offsets
    c, d = class_offsets
    blocks = [
        tuple(sorted((x, (x + a) % order, (x + b) % order)))
        for x in range(order)
    ]
    blocks += [(3 * i, 3 * i + 1, 3 * i + 2) for i in range(g)]
    blocks += [
        tuple(sorted((3 * i, 3 * ((i + c) % g) + 1,
                      3 * ((i + d) % g) + 2)))
        for i in range(g)
    ]
    return blocks


def verify(order, blocks):
    if len(blocks) != 5 * order // 3 or len(set(blocks)) != len(blocks):
        return False
    if any(len(set(block)) != 3 for block in blocks):
        return False
    degrees = Counter(point for block in blocks for point in block)
    if {degrees[point] for point in range(order)} != {5}:
        return False
    pairs = set()
    for block in blocks:
        for pair in combinations(block, 2):
            if pair in pairs:
                return False
            pairs.add(pair)
    return True


def search(order):
    g = order // 3
    for a in range(1, order):
        for b in range(a + 1, order):
            for c in range(1, g):
                for d in range(1, g):
                    if c == d:
                        continue
                    blocks = candidate(order, (a, b), (c, d))
                    if verify(order, blocks):
                        return (a, b), (c, d), blocks
    return None


if __name__ == "__main__":
    for order in (18, 21, 33):
        result = search(order)
        if result is None:
            print(f"order {order}: no module in the searched cyclic family")
        else:
            orbit, parallel, _ = result
            print(f"order {order}: orbit={orbit}, second class={parallel}")
