"""Search for finite E(H)=8 modules for the n5=3,6,9 exceptions.

A selected packing is completed by K4s in its dense leave.  The local weighted
degrees are constrained to be divisible by 3, making that leave K4-divisible.
"""

import argparse
from itertools import combinations


def pairs(block):
    return combinations(block, 2)


def verify(order, target, blocks, uncovered, expected_excess=8):
    seen = set()
    local = [0] * order
    excess = 2 * len(uncovered)
    correction = len(uncovered)
    for block in blocks:
        size = len(block)
        excess += (size - 4) * (size - 5)
        correction += size * (size - 1) // 2 - 6
        for point in block:
            local[point] += size - 1
        for pair in pairs(block):
            assert pair not in seen
            seen.add(pair)
    for pair in uncovered:
        assert pair not in seen
        for point in pair:
            local[point] += 1
    assert excess == expected_excess
    assert correction == 4 * target
    assert all(value % 3 == 0 for value in local)


def search(order, target, counts, uncovered_count=0, time_limit=20):
    from ortools.sat.python import cp_model

    model = cp_model.CpModel()
    candidates = []
    allowed_sizes = (2, 3, 5, 6, 7)
    by_size = {size: [] for size in allowed_sizes}
    by_pair = {pair: [] for pair in combinations(range(order), 2)}
    by_point = [[] for _ in range(order)]

    for size in allowed_sizes:
        if not counts.get(size, 0):
            continue
        for block in combinations(range(order), size):
            index = len(candidates)
            candidates.append(block)
            by_size[size].append(index)
            for pair in pairs(block):
                by_pair[pair].append(index)
            for point in block:
                by_point[point].append(index)

    use = [model.new_bool_var(f"b_{i}") for i in range(len(candidates))]
    miss = {
        pair: model.new_bool_var(f"d_{pair[0]}_{pair[1]}")
        for pair in by_pair
    } if uncovered_count else {}
    for size, count in counts.items():
        model.add(sum(use[i] for i in by_size[size]) == count)
    for pair, indices in by_pair.items():
        model.add(sum(use[i] for i in indices) + (miss[pair] if miss else 0) <= 1)
    if miss:
        model.add(sum(miss.values()) == uncovered_count)

    weights = {size: size - 1 for size in allowed_sizes}
    sizes = [len(block) for block in candidates]
    max_local = sum(weights[size] * count for size, count in counts.items()) + uncovered_count
    for point, indices in enumerate(by_point):
        quotient = model.new_int_var(0, max_local // 3, f"q_{point}")
        missing_degree = (
            sum(var for pair, var in miss.items() if point in pair) if miss else 0
        )
        model.add(
            sum(weights[sizes[i]] * use[i] for i in indices) + missing_degree
            == 3 * quotient
        )

    # Vertex transitivity lets us prescribe one block of the largest used size.
    anchor_size = max(size for size, count in counts.items() if count)
    anchor = tuple(range(anchor_size))
    anchor_index = candidates.index(anchor)
    model.add(use[anchor_index] == 1)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    status = solver.solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None, solver.status_name(status)
    blocks = [candidates[i] for i in range(len(candidates)) if solver.value(use[i])]
    uncovered = [pair for pair, var in miss.items() if solver.value(var)]
    verify(order, target, blocks, uncovered)
    return (blocks, uncovered), solver.status_name(status)


def count_options(target, max_uncovered=4):
    """Enumerate every E=8 signature with block sizes at least two.

    Sizes outside 2,3,4,5,6,7 already contribute more than eight to E;
    size-4 blocks are irrelevant to a minimal fixed module.
    """
    for uncovered in range(max_uncovered + 1):
        for pairs_count in range(2):
            for septuples in range(2):
                remaining = 4 - uncovered - 3 * (pairs_count + septuples)
                if remaining < 0:
                    continue
                for triples in range(remaining + 1):
                    sextuples = remaining - triples
                    numerator = (
                        4 * target
                        - uncovered
                        + 5 * pairs_count
                        + 3 * triples
                        - 9 * sextuples
                        - 15 * septuples
                    )
                    if numerator >= 0 and numerator % 4 == 0:
                        counts = {
                            2: pairs_count,
                            3: triples,
                            5: numerator // 4,
                            6: sextuples,
                            7: septuples,
                        }
                        yield {k: v for k, v in counts.items() if v}, uncovered


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=int, choices=(3, 6, 9))
    parser.add_argument("--order", type=int)
    parser.add_argument("--signature", type=int, help="zero-based signature index")
    parser.add_argument("--time-limit", type=float, default=20)
    parser.add_argument("--list", action="store_true", help="list finite signatures")
    args = parser.parse_args()

    targets = (args.target,) if args.target else (3, 6, 9)
    signatures = {
        target: list(count_options(target))
        for target in targets
    }
    if args.list or args.order is None:
        for target in targets:
            for index, (counts, uncovered_count) in enumerate(signatures[target]):
                print(
                    f"t={target}, signature={index}, counts={counts}, "
                    f"uncovered={uncovered_count}"
                )
        return

    if args.target is None or args.signature is None:
        parser.error("--order requires both --target and --signature")
    options = signatures[args.target]
    if not 0 <= args.signature < len(options):
        parser.error(f"signature must lie in [0, {len(options) - 1}]")
    counts, uncovered_count = options[args.signature]
    module, status = search(
        args.order,
        args.target,
        counts,
        uncovered_count,
        args.time_limit,
    )
    print(
        f"t={args.target}, order={args.order}, counts={counts}, "
        f"uncovered={uncovered_count}: {status}"
    )
    if module is not None:
        blocks, uncovered = module
        for block in blocks:
            print(block)
        print("uncovered:", uncovered)


if __name__ == "__main__":
    main()
