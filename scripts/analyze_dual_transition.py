"""Classify the small-excess modules at the K5 side of the p=4 range.

For m == 0 (mod 4), a finite collection P of non-5-blocks and a graph D
leaves prescribed degrees for the background 4-block packing.  This script
enumerates every signature with excess E <= 24 and maximises the transition
parameter Delta.  Points are encoded by the labelled special blocks and
uncovered edges that contain them, so pairwise linearity is imposed exactly.
"""

from itertools import combinations

from ortools.sat.python import cp_model


BLOCK_EXCESS = {1: 12, 2: 6, 3: 2, 6: 2, 7: 6, 8: 12, 9: 20}
BLOCK_CORRECTION = {k: k * (k - 1) // 2 - 10 for k in BLOCK_EXCESS}


def signatures(excess):
    items = list(BLOCK_EXCESS.items()) + [("D", 2)]

    def generate(index, remaining, counts):
        if index == len(items):
            if remaining == 0:
                yield dict(counts)
            return
        kind, contribution = items[index]
        for number in range(remaining // contribution + 1):
            if number:
                counts[kind] = number
            yield from generate(
                index + 1, remaining - number * contribution, counts
            )
            counts.pop(kind, None)

    yield from generate(0, excess, {})


def labelled_objects(signature):
    objects = []
    for kind in sorted(k for k in signature if k != "D"):
        for index in range(signature[kind]):
            objects.append((f"B{kind}_{index}", kind, kind - 1))
    for index in range(signature.get("D", 0)):
        objects.append((f"D_{index}", 2, 1))
    return objects


def point_value(objects, members):
    weight = sum(objects[index][2] for index in members)
    residual_degree = (1 + weight) % 4
    return 1 - residual_degree


def build_model(signature):
    objects = labelled_objects(signature)
    point_types = []
    for mask in range(1, 1 << len(objects)):
        members = tuple(i for i in range(len(objects)) if mask & (1 << i))
        point_types.append((members, point_value(objects, members)))

    model = cp_model.CpModel()
    multiplicities = []
    for index, (members, _) in enumerate(point_types):
        upper = max(objects[i][1] for i in members) if len(members) == 1 else 1
        multiplicities.append(model.new_int_var(0, upper, f"x_{index}"))

    for object_index, (_, size, _) in enumerate(objects):
        model.add(
            sum(
                multiplicities[i]
                for i, (members, _) in enumerate(point_types)
                if object_index in members
            )
            == size
        )

    for first, second in combinations(range(len(objects)), 2):
        model.add(
            sum(
                multiplicities[i]
                for i, (members, _) in enumerate(point_types)
                if first in members and second in members
            )
            <= 1
        )

    point_contribution = sum(
        value * multiplicities[i]
        for i, (_, value) in enumerate(point_types)
    )
    correction = sum(
        signature.get(k, 0) * value for k, value in BLOCK_CORRECTION.items()
    ) + signature.get("D", 0)
    return model, multiplicities, point_types, point_contribution, correction


def maximise_delta(signature, time_limit=10):
    model, multiplicities, point_types, point_contribution, correction = (
        build_model(signature)
    )
    model.maximize(point_contribution)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    if status != cp_model.OPTIMAL:
        return None, solver.status_name(status), None

    delta = correction + round(solver.objective_value)
    selected = [
        members
        for i, (members, _) in enumerate(point_types)
        for _ in range(solver.value(multiplicities[i]))
    ]
    return delta, solver.status_name(status), selected


def exceeds_delta(signature, threshold, time_limit=30):
    model, _, _, point_contribution, correction = build_model(signature)
    model.add(point_contribution + correction >= threshold)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    return solver.status_name(status)


def crude_upper_bound(signature):
    # Every point contributes at most 1 to sum_x(1-r(x)).
    correction = sum(
        signature.get(k, 0) * value for k, value in BLOCK_CORRECTION.items()
    ) + signature.get("D", 0)
    incidences = sum(
        signature.get(k, 0) * k for k in BLOCK_CORRECTION
    ) + 2 * signature.get("D", 0)
    return correction + incidences


def classify(max_deficit=3, time_limit=60):
    expected = {1: 0, 2: 24, 3: 68}
    for deficit in range(1, max_deficit + 1):
        signature_list = list(signatures(8 * deficit))
        candidates = [
            signature
            for signature in signature_list
            if crude_upper_bound(signature) > expected[deficit]
        ]
        for signature in candidates:
            status = exceeds_delta(
                signature, expected[deficit] + 1, time_limit=time_limit
            )
            if status != "INFEASIBLE":
                raise RuntimeError(
                    f"unresolved signature for s={deficit}: {signature} ({status})"
                )
        print(
            f"s={deficit}: all {len(signature_list)} signatures certified; "
            f"maximum Delta={expected[deficit]} "
            f"({len(candidates)} nontrivial CP-SAT checks)"
        )


if __name__ == "__main__":
    classify()
