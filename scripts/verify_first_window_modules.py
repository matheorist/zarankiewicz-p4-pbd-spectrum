"""Finite checks for the nonintegral first-transition module table.

The script checks the restricted triple-packing numbers used for beta=1,
the small beta=2 and beta=3 module lists, and the arithmetic exceptional
table printed in the manuscript.  It is an ancillary verification only;
the asymptotic completion is supplied by the dense K_4-decomposition theorem.
"""

from itertools import combinations
from math import ceil

from ortools.sat.python import cp_model


TAU = (0, 0, 1, 1, 2, 4, 7, 8, 12, 13, 17, 20)

# Explicit 6-regular certificates for the six finite orders not covered by
# the two cyclic orbits used for k >= 19.  Each tuple is a triple of labels.
FINITE_TAU_CERTIFICATES = {
    13: [(0, 1, 9), (0, 2, 11), (0, 3, 10), (0, 4, 7), (0, 5, 8),
         (0, 6, 12), (1, 2, 4), (1, 3, 12), (1, 5, 7), (1, 6, 11),
         (1, 8, 10), (2, 3, 9), (2, 5, 10), (2, 6, 8), (2, 7, 12),
         (3, 4, 6), (3, 5, 11), (3, 7, 8), (4, 5, 12), (4, 8, 11),
         (4, 9, 10), (5, 6, 9), (6, 7, 10), (7, 9, 11), (8, 9, 12),
         (10, 11, 12)],
    14: [(0, 1, 12), (0, 2, 11), (0, 3, 7), (0, 4, 5), (0, 6, 8),
         (0, 10, 13), (1, 2, 9), (1, 4, 10), (1, 5, 7), (1, 6, 13),
         (1, 8, 11), (2, 3, 8), (2, 4, 12), (2, 6, 10), (2, 7, 13),
         (3, 4, 6), (3, 5, 13), (3, 9, 11), (3, 10, 12), (4, 7, 8),
         (4, 9, 13), (5, 6, 9), (5, 8, 12), (5, 10, 11), (6, 7, 11),
         (7, 9, 12), (8, 9, 10), (11, 12, 13)],
    15: [(0, 1, 14), (0, 2, 4), (0, 3, 8), (0, 5, 12), (0, 6, 11),
         (0, 7, 10), (1, 2, 6), (1, 3, 11), (1, 4, 5), (1, 7, 12),
         (1, 8, 9), (2, 7, 9), (2, 8, 14), (2, 10, 12), (2, 11, 13),
         (3, 4, 14), (3, 5, 9), (3, 6, 12), (3, 10, 13), (4, 7, 13),
         (4, 8, 11), (4, 9, 12), (5, 7, 11), (5, 8, 13), (5, 10, 14),
         (6, 7, 14), (6, 8, 10), (6, 9, 13), (9, 10, 11), (12, 13, 14)],
    16: [(0, 1, 9), (0, 2, 5), (0, 3, 13), (0, 4, 6), (0, 7, 8),
         (0, 10, 14), (1, 2, 11), (1, 3, 15), (1, 5, 10), (1, 6, 8),
         (1, 12, 13), (2, 3, 8), (2, 4, 14), (2, 6, 13), (2, 7, 15),
         (3, 4, 5), (3, 7, 10), (3, 11, 14), (4, 7, 12), (4, 9, 11),
         (4, 10, 15), (5, 6, 11), (5, 7, 14), (5, 9, 12), (6, 7, 9),
         (6, 12, 15), (8, 9, 15), (8, 11, 13), (8, 12, 14), (9, 10, 13),
         (10, 11, 12), (13, 14, 15)],
    17: [(0, 2, 3), (0, 4, 11), (0, 5, 14), (0, 7, 16), (0, 8, 12),
         (0, 10, 13), (1, 3, 14), (1, 4, 5), (1, 6, 7), (1, 9, 13),
         (1, 11, 16), (1, 12, 15), (2, 5, 16), (2, 6, 13), (2, 7, 8),
         (2, 9, 15), (2, 10, 12), (3, 4, 7), (3, 5, 11), (3, 6, 16),
         (3, 10, 15), (4, 6, 15), (4, 9, 10), (4, 13, 14), (5, 6, 12),
         (5, 8, 9), (6, 8, 11), (7, 9, 11), (7, 10, 14), (8, 10, 16),
         (8, 13, 15), (9, 12, 14), (11, 12, 13), (14, 15, 16)],
    18: [(0, 1, 15), (0, 2, 17), (0, 3, 5), (0, 6, 7), (0, 10, 14),
         (0, 11, 13), (1, 2, 16), (1, 3, 14), (1, 4, 9), (1, 5, 10),
         (1, 6, 12), (2, 4, 13), (2, 5, 15), (2, 6, 8), (2, 9, 14),
         (3, 4, 12), (3, 6, 10), (3, 9, 15), (3, 13, 17), (4, 7, 16),
         (4, 10, 17), (4, 11, 14), (5, 6, 16), (5, 7, 13), (5, 8, 17),
         (6, 11, 15), (7, 8, 14), (7, 9, 12), (7, 11, 17), (8, 9, 13),
         (8, 10, 11), (8, 12, 15), (9, 11, 16), (10, 12, 16),
         (12, 13, 14), (15, 16, 17)],
}


def verify_linear_packing(k, triples):
    assert len(triples) == 2 * k
    assert all(len(set(triple)) == 3 for triple in triples)
    assert all(0 <= vertex < k for triple in triples for vertex in triple)
    assert len({tuple(sorted(triple)) for triple in triples}) == len(triples)
    pair_sets = [set(triple) for triple in triples]
    for first, second in combinations(range(len(pair_sets)), 2):
        assert len(pair_sets[first] & pair_sets[second]) <= 1
    degrees = [sum(vertex in triple for triple in triples) for vertex in range(k)]
    assert degrees == [6] * k


def cyclic_certificate(k):
    return [
        tuple(sorted((shift + vertex) % k for vertex in base))
        for base in ((0, 1, 3), (0, 4, 9))
        for shift in range(k)
    ]


def triple_capacity(k):
    triples = list(combinations(range(k), 3))
    model = cp_model.CpModel()
    chosen = [model.new_bool_var(f"t_{i}") for i in range(len(triples))]
    for pair in combinations(range(k), 2):
        model.add(
            sum(chosen[i] for i, triple in enumerate(triples) if set(pair) <= set(triple))
            <= 1
        )
    for vertex in range(k):
        model.add(
            sum(chosen[i] for i, triple in enumerate(triples) if vertex in triple)
            <= 6
        )
    model.maximize(sum(chosen))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 20
    solver.parameters.num_search_workers = 8
    status = solver.solve(model)
    if status != cp_model.OPTIMAL:
        raise RuntimeError(f"triple packing optimization failed for k={k}")
    return round(solver.objective_value)


def check_tau():
    for k, expected in enumerate(TAU, start=1):
        actual = triple_capacity(k)
        assert actual == expected, (k, actual, expected)
    for k, certificate in FINITE_TAU_CERTIFICATES.items():
        verify_linear_packing(k, certificate)
    for k in (19, 20, 31, 64):
        verify_linear_packing(k, cyclic_certificate(k))
    print("beta=1 finite capacities and k >= 19 cyclic certificates: OK")


def module_feasible(k, beta, target):
    """Check a module made from uncovered pairs, triples, and sextuples."""
    sizes = {"D": 2, "T": 3, "S": 6}
    weights = {"D": 1, "T": 2, "S": 5}
    corrections = {"D": 1, "T": -3, "S": 5}
    for d in range(k + 1):
        for t in range(k - d + 1):
            s = k - d - t
            kinds = ["D"] * d + ["T"] * t + ["S"] * s
            objects = [(sizes[kind], weights[kind]) for kind in kinds]
            point_types = []
            for mask in range(1, 1 << k):
                members = tuple(i for i in range(k) if mask & (1 << i))
                weight = sum(objects[i][1] for i in members)
                value = beta - (beta + weight) % 4
                point_types.append((members, value))

            model = cp_model.CpModel()
            variables = []
            for index, (members, _) in enumerate(point_types):
                upper = max(objects[i][0] for i in members) if len(members) == 1 else 1
                variables.append(model.new_int_var(0, upper, f"x_{index}"))
            for object_index, (size, _) in enumerate(objects):
                model.add(
                    sum(
                        variables[i]
                        for i, (members, _) in enumerate(point_types)
                        if object_index in members
                    )
                    == size
                )
            for first, second in combinations(range(k), 2):
                model.add(
                    sum(
                        variables[i]
                        for i, (members, _) in enumerate(point_types)
                        if first in members and second in members
                    )
                    <= 1
                )
            correction = d - 3 * t + 5 * s
            model.add(
                correction
                + sum(value * variables[i] for i, (_, value) in enumerate(point_types))
                == target
            )
            solver = cp_model.CpSolver()
            solver.parameters.max_time_in_seconds = 2
            solver.parameters.num_search_workers = 8
            status = solver.solve(model)
            if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                return True
    return False


def check_small_modules():
    beta2 = {1: [3], 2: [2, 6], 3: [1, 5, 9]}
    beta3 = {
        1: [3, 7, 23],
        2: [2, 6, 10, 14, 22, 26, 30, 42, 46],
        3: list(range(1, 70, 4)),
    }
    for beta, data in ((2, beta2), (3, beta3)):
        for k, targets in data.items():
            for target in targets:
                assert module_feasible(k, beta, target), (beta, k, target)
    print("beta=2,3 finite module lists: OK")


def arithmetic_exceptions():
    denominators = {1: 28, 2: 22, 3: 92}
    exceptions = {
        1: {
            1: {3: 1, 7: 2, 11: 2, 15: 2, 19: 2, 23: 2,
                27: 2, 31: 2, 35: 2, 43: 3, 47: 3, 51: 3,
                55: 3, 59: 3, 63: 3},
            2: {2: 1, 6: 1, 10: 1, 14: 2, 18: 2, 22: 2,
                30: 2, 34: 2, 38: 2, 42: 2, 50: 3, 54: 3,
                58: 3, 62: 3, 66: 3, 70: 3},
            3: {5: 1, 9: 1, 13: 1, 17: 1, 21: 1, 29: 2,
                33: 2, 37: 2, 41: 2, 45: 2, 49: 2, 61: 3,
                65: 3, 69: 3, 73: 3, 77: 3},
        },
        2: {2: {10: 1}, 3: {13: 1}},
        3: {1: {11: 1, 15: 1, 19: 1}, 2: {18: 1, 34: 1, 38: 1}},
    }
    for beta, denominator in denominators.items():
        count = 0
        for delta in range(1, 1680):
            if delta % 4 == 0:
                continue
            theta = (-delta) % 4
            j = ceil(delta / denominator - theta / 4 - 1e-12)
            rho = exceptions.get(beta, {}).get(theta, {}).get(delta, j)
            assert rho >= j
            count += 1
        assert count == 1260
    print("first-window arithmetic classes: 1260 per beta")


def check_generic_layer_arithmetic():
    """Check the interval reductions used in the generic module proof."""
    for theta in (1, 2, 3):
        for s in range(1, 20):
            k = theta + 4 * s

            # beta=1: Delta=4T-k and the generic interval requires
            # T in {2k-6,...,2k}.
            beta1 = [
                delta
                for delta in range(1, 2000)
                if 28 * (s - 1) + 7 * theta < delta <= 28 * s + 7 * theta
                and (-delta) % 4 == theta % 4
            ]
            assert {((delta + k) // 4) for delta in beta1} <= set(
                range(2 * k - 6, 2 * k + 1)
            )

            # beta=2: Delta=-k+4e and a 6-regular graph supplies all
            # edge counts up to 3k once k >= 7.
            beta2 = [
                delta
                for delta in range(1, 2000)
                if 22 * (s - 1) + 11 * theta / 2 < delta <= 22 * s + 11 * theta / 2
                and (-delta) % 4 == theta % 4
            ]
            if k >= 7:
                edges = {(delta + k) // 4 for delta in beta2}
                assert edges <= set(range(3 * k + 1))

            # beta=3: Delta=23k-4(5d+e), and the needed reductions are
            # exactly 0,...,23.
            reductions = {
                (23 * k - delta) // 4
                for delta in range(1, 2000)
                if 92 * (s - 1) + 23 * theta < delta <= 92 * s + 23 * theta
                and (-delta) % 4 == theta % 4
            }
            assert reductions <= set(range(24))
            assert {
                5 * d + e
                for d in range(5)
                for e in range(5 if d < 4 else 4)
            } == set(range(24))
    print("generic first-window layer arithmetic: OK")


if __name__ == "__main__":
    check_tau()
    check_small_modules()
    arithmetic_exceptions()
    check_generic_layer_arithmetic()
