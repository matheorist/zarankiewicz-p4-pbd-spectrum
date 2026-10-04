"""Joint CP-SAT search for a 15-point pure {4,5}-trade of volume (15,9)."""

from itertools import combinations

from ortools.sat.python import cp_model


V = tuple(range(15))
PAIRS = tuple(combinations(V, 2))
BLOCKS4 = tuple(combinations(V, 4))
BLOCKS5 = tuple(combinations(V, 5))


def incidence(blocks, objects):
    return {
        obj: [i for i, block in enumerate(blocks) if set(obj) <= set(block)]
        for obj in objects
    }


pair4 = incidence(BLOCKS4, PAIRS)
pair5 = incidence(BLOCKS5, PAIRS)
point4 = incidence(BLOCKS4, ((x,) for x in V))
point5 = incidence(BLOCKS5, ((x,) for x in V))

model = cp_model.CpModel()
x4 = [model.new_bool_var(f"q_{i}") for i in range(len(BLOCKS4))]
x5 = [model.new_bool_var(f"p_{i}") for i in range(len(BLOCKS5))]

for pair in PAIRS:
    lhs = sum(x4[i] for i in pair4[pair])
    rhs = sum(x5[i] for i in pair5[pair])
    model.add(lhs == rhs)
    model.add(lhs <= 1)

for point in V:
    model.add(sum(x4[i] for i in point4[(point,)]) == 4)
    model.add(sum(x5[i] for i in point5[(point,)]) == 3)

# Vertex transitivity lets us prescribe one quintuple.
model.add(x5[BLOCKS5.index((0, 1, 2, 3, 4))] == 1)

solver = cp_model.CpSolver()
solver.parameters.num_search_workers = 8
solver.parameters.max_time_in_seconds = 600
solver.parameters.log_search_progress = True
status = solver.solve(model)

print(f"status: {solver.status_name(status)}")
if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    quads = [BLOCKS4[i] for i in range(len(BLOCKS4)) if solver.value(x4[i])]
    quints = [BLOCKS5[i] for i in range(len(BLOCKS5)) if solver.value(x5[i])]
    assert len(quads) == 15 and len(quints) == 9
    print("quadruples:")
    for block in quads:
        print(block)
    print("quintuples:")
    for block in quints:
        print(block)
