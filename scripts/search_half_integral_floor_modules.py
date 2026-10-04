"""Finite incidence search for the six d=0 half-integral signatures.

This exploratory CP-SAT model searches for finite linear packings. It is
evidence for or against small-module realizability, not by itself a proof of
eventual realizability.
"""

from itertools import combinations
from ortools.sat.python import cp_model

SIGNATURES = ((0, 0, 2), (0, 1, 1), (0, 2, 0),
              (1, 0, 1), (1, 1, 0), (2, 0, 0))


def solve(h, signature, time_limit=20.0):
    q3, q6, e = signature
    q4 = h - 2 * q3 + q6
    if q4 < 0:
        return None

    labels = [f"Q4_{i}" for i in range(q4)]
    labels += [f"Q3_{i}" for i in range(q3)]
    labels += [f"Q6_{i}" for i in range(q6)]
    labels += [f"D_{i}" for i in range(e)]
    weights = [4] * q4 + [3] * q3 + [6] * q6 + [2] * e
    n = len(labels)

    types = []
    for mask in range(1, 1 << n):
        support = tuple(i for i in range(n) if mask & (1 << i))
        if sum(weights[i] for i in support) % 3 == 0:
            types.append(support)

    model = cp_model.CpModel()
    xs = []
    for j, support in enumerate(types):
        upper = max(weights[i] for i in support) if len(support) == 1 else 1
        xs.append(model.NewIntVar(0, upper, f"x_{j}"))

    for i, weight in enumerate(weights):
        model.Add(sum(xs[j] for j, support in enumerate(types) if i in support) == weight)

    for i, k in combinations(range(n), 2):
        model.Add(sum(xs[j] for j, support in enumerate(types)
                      if i in support and k in support) <= 1)

    model.Maximize(sum(xs))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    status = solver.Solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return {"status": solver.StatusName(status), "q4": q4}

    active = []
    for j, support in enumerate(types):
        value = solver.Value(xs[j])
        if value:
            active.append((value, tuple(labels[i] for i in support)))
    return {"status": solver.StatusName(status), "q4": q4,
            "points": sum(value for value, _ in active), "active": active}


def main():
    for h in range(1, 13):
        print(f"h={h}")
        for signature in SIGNATURES:
            result = solve(h, signature)
            if result is None:
                print(f"  {signature}: skipped (q4<0)")
            elif result["status"] in ("OPTIMAL", "FEASIBLE"):
                print(f"  {signature}: {result['status']}, q4={result['q4']}, points={result['points']}")
                for value, support in result["active"]:
                    print(f"    {value}x {', '.join(support)}")
            else:
                print(f"  {signature}: {result['status']}, q4={result['q4']}")


if __name__ == "__main__":
    main()
