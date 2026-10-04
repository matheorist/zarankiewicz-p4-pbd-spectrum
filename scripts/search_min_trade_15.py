"""Search for the minimum pure {4,5}-trade on 15 points.

The 15 translates of {0,1,4,6} modulo 15 decompose K_15 minus the
difference-7 Hamilton cycle into K_4s.  This script asks whether the same
90-edge graph has a K_5-decomposition.
"""

from itertools import combinations

from ortools.sat.python import cp_model


V = range(15)
BASE4 = (0, 1, 4, 6)


def pairs(block):
    return {tuple(sorted(pair)) for pair in combinations(block, 2)}


quadruples = [tuple(sorted((x + a) % 15 for a in BASE4)) for x in V]
target = set().union(*(pairs(block) for block in quadruples))
assert len(target) == 90

candidates = [block for block in combinations(V, 5) if pairs(block) <= target]
by_pair = {pair: [] for pair in target}
for i, block in enumerate(candidates):
    for pair in pairs(block):
        by_pair[pair].append(i)

model = cp_model.CpModel()
use = [model.new_bool_var(f"b_{i}") for i in range(len(candidates))]
for ids in by_pair.values():
    model.add_exactly_one(use[i] for i in ids)

solver = cp_model.CpSolver()
solver.parameters.num_search_workers = 8
solver.parameters.max_time_in_seconds = 300
status = solver.solve(model)

print(f"candidate K5 blocks: {len(candidates)}")
print(f"status: {solver.status_name(status)}")
if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    quintuples = [candidates[i] for i in range(len(candidates)) if solver.value(use[i])]
    assert len(quintuples) == 9
    assert set().union(*(pairs(block) for block in quintuples)) == target
    print("quadruples:")
    for block in quadruples:
        print(block)
    print("quintuples:")
    for block in quintuples:
        print(block)
