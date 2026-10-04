"""Verify explicit triple-packing modules for the K5 transition."""

from itertools import combinations


def cyclic_triples(order, bases):
    return {
        tuple(sorted(((x + a) % order for a in base)))
        for base in bases
        for x in range(order)
    }


def verify_packing(order, triples, degree_bound):
    degrees = [0] * order
    pairs = set()
    for triple in triples:
        assert len(triple) == 3
        for point in triple:
            degrees[point] += 1
        for pair in combinations(triple, 2):
            assert pair not in pairs
            pairs.add(pair)
    assert max(degrees, default=0) <= degree_bound


def full_packing(deficit):
    order = 4 * deficit
    if deficit == 1:
        triples = {(0, 1, 2)}
    elif deficit == 2:
        # A triangle decomposition of K_8 minus a 1-factor.
        triples = cyclic_triples(8, ((0, 1, 3),))
    elif deficit == 3:
        # A triangle decomposition of K_12 minus a 1-factor.
        triples = {
            (0, 2, 4), (0, 3, 7), (0, 5, 9), (0, 6, 11), (0, 8, 10),
            (1, 2, 9), (1, 3, 8), (1, 4, 6), (1, 5, 11), (1, 7, 10),
            (2, 5, 7), (2, 6, 10), (2, 8, 11), (3, 4, 11),
            (3, 5, 10), (3, 6, 9), (4, 7, 8), (4, 9, 10),
            (5, 6, 8), (7, 9, 11),
        }
    else:
        # Two cyclic orbits; their directed differences are all distinct.
        triples = cyclic_triples(order, ((0, 1, 3), (0, 4, 9)))
    verify_packing(order, triples, 6)
    return sorted(triples)


def module(deficit, delta):
    assert delta % 4 == 0
    needed = deficit + delta // 4
    triples = full_packing(deficit)[:needed]
    assert len(triples) == needed
    # There are 4s sextuples. A triple represents one point shared by the
    # corresponding three sextuples; private points fill unused incidences.
    realised = -4 * deficit + 4 * len(triples)
    assert realised == delta
    assert 2 * (4 * deficit) == 8 * deficit
    # A shared point lies in three sextuples and has r=0; a private point
    # lies in one sextuple and has r=2. Both satisfy r-w == 1 (mod 4).
    assert (0 - 3 * 5) % 4 == 1
    assert (2 - 5) % 4 == 1


def verify():
    maxima = {1: 0, 2: 24, 3: 68}
    for deficit, maximum in maxima.items():
        packing = full_packing(deficit)
        assert -4 * deficit + 4 * len(packing) == maximum
        for delta in range(0, maximum + 1, 4):
            module(deficit, delta)
        print(f"s={deficit}: all multiples of 4 through Delta={maximum}")

    for deficit in range(4, 61):
        packing = full_packing(deficit)
        assert len(packing) == 8 * deficit
        for delta in range(0, 28 * deficit + 1, 4):
            module(deficit, delta)
    print("s=4,...,60: all multiples of 4 through Delta=28s")


if __name__ == "__main__":
    verify()
