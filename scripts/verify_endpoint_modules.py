"""Verify the dual endpoint modules used near the 4-space threshold."""

from collections import Counter
from itertools import combinations


def triple_module(order):
    """A 5-regular linear triple system on order vertices, for 3 | order."""
    assert order % 3 == 0 and order >= 15
    g = order // 3
    blocks = [
        tuple(sorted((x, (x + 3) % order, (x + 7) % order)))
        for x in range(order)
    ]
    blocks += [(3 * i, 3 * i + 1, 3 * i + 2) for i in range(g)]
    blocks += [
        tuple(sorted((3 * i, 3 * ((i + 3) % g) + 1,
                      3 * ((i + 1) % g) + 2)))
        for i in range(g)
    ]
    return blocks


def one_four_module(order, vertices):
    """A degree-5 dual with order edges of sizes 1 and 4."""
    high_edges = (5 * vertices - order) // 3
    assert 5 * vertices - order == 3 * high_edges
    assert 0 <= 4 * high_edges <= vertices
    edges = [tuple(range(4 * i, 4 * i + 4)) for i in range(high_edges)]
    high_degree = Counter(point for edge in edges for point in edge)
    for point in range(vertices):
        edges += [(point,)] * (5 - high_degree[point])
    assert len(edges) == order
    return edges


def verify_linear_regular(vertices, edges, degree, sizes):
    assert all(len(edge) in sizes and len(set(edge)) == len(edge) for edge in edges)
    incidence = Counter(point for edge in edges for point in edge)
    assert {incidence[point] for point in range(vertices)} == {degree}
    pairs = set()
    for edge in edges:
        for pair in combinations(sorted(edge), 2):
            assert pair not in pairs
            pairs.add(pair)


def mu(order):
    return min(value for value in range((order + 4) // 5, order + 1)
               if value % 3 == 1)


def main():
    for vertices in range(15, 301, 3):
        edges = triple_module(vertices)
        assert len(edges) == 5 * vertices // 3
        verify_linear_regular(vertices, edges, 5, {3})

    tested = 0
    for order in range(101, 601):
        if order % 12 not in {5, 8}:
            continue
        for vertices in range(mu(order), 4 * order // 17 + 1, 3):
            edges = one_four_module(order, vertices)
            verify_linear_regular(vertices, edges, 5, {1, 4})
            tested += 1

    print("Triple modules verified for every multiple of 3 from 15 through 300.")
    print(f"One-four modules verified for {tested} endpoint parameter pairs.")


if __name__ == "__main__":
    main()
