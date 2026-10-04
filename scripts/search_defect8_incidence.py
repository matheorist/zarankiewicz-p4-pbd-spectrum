"""Exact incidence-type search for the finite E=8 reduction.

Each point is represented by the subset of labelled special blocks containing
it.  Pairwise linearity says that every pair of block labels occurs in at most
one point type.  Uncovered edges are represented as labelled size-2
pseudo-blocks; this gives exactly the same local congruence and intersection
conditions as the graph D in Proposition 9.
"""

import argparse
from itertools import combinations

from ortools.sat.python import cp_model

from search_defect8_modules import count_options, verify


def labelled_blocks(counts, uncovered_count):
    labels = []
    real = []
    for size in sorted(counts):
        for index in range(counts[size]):
            labels.append((f"B{size}_{index}", size, False))
            real.append(len(labels) - 1)
    uncovered = []
    for index in range(uncovered_count):
        labels.append((f"D_{index}", 2, True))
        uncovered.append(len(labels) - 1)
    return labels, real, uncovered


def admissible_types(labels):
    block_count = len(labels)
    weights = [size - 1 for _, size, _ in labels]
    types = []
    for mask in range(1, 1 << block_count):
        members = tuple(i for i in range(block_count) if mask & (1 << i))
        if sum(weights[i] for i in members) % 3 == 0:
            types.append(members)
    return types


def solve(target, counts, uncovered_count, time_limit=60, expected_excess=8):
    labels, real_indices, uncovered_indices = labelled_blocks(
        counts, uncovered_count
    )
    types = admissible_types(labels)
    model = cp_model.CpModel()
    multiplicity = []
    for index, point_type in enumerate(types):
        upper = max(labels[i][1] for i in point_type) if len(point_type) == 1 else 1
        multiplicity.append(model.new_int_var(0, upper, f"x_{index}"))

    for block_index, (_, size, _) in enumerate(labels):
        model.add(
            sum(
                multiplicity[i]
                for i, point_type in enumerate(types)
                if block_index in point_type
            )
            == size
        )

    for first, second in combinations(range(len(labels)), 2):
        model.add(
            sum(
                multiplicity[i]
                for i, point_type in enumerate(types)
                if first in point_type and second in point_type
            )
            <= 1
        )

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    solver.parameters.symmetry_level = 2
    status = solver.solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None, solver.status_name(status), len(types)

    block_points = [[] for _ in labels]
    point = 0
    for index, point_type in enumerate(types):
        for _ in range(solver.value(multiplicity[index])):
            for block_index in point_type:
                block_points[block_index].append(point)
            point += 1
    blocks = [tuple(block_points[i]) for i in real_indices]
    uncovered = [tuple(block_points[i]) for i in uncovered_indices]
    verify(point, target, blocks, uncovered, expected_excess)
    return (blocks, uncovered), solver.status_name(status), len(types)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=int, required=True, choices=(3, 6, 9))
    parser.add_argument("--signature", type=int, required=True)
    parser.add_argument("--time-limit", type=float, default=60)
    args = parser.parse_args()

    options = list(count_options(args.target))
    if not 0 <= args.signature < len(options):
        parser.error(f"signature must lie in [0, {len(options) - 1}]")
    counts, uncovered_count = options[args.signature]
    module, status, type_count = solve(
        args.target, counts, uncovered_count, args.time_limit
    )
    print(
        f"t={args.target}, signature={args.signature}, counts={counts}, "
        f"uncovered={uncovered_count}, types={type_count}: {status}"
    )
    if module is not None:
        blocks, uncovered = module
        for block in blocks:
            print(block)
        print("uncovered:", uncovered)


if __name__ == "__main__":
    main()
