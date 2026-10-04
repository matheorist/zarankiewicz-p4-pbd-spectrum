"""
Verify the fixed-degree bad-residue modules.

For fixed R == 2 mod 3 and q == 0,9 mod 12, the manuscript constructs
v=5q points and N=Rq quintuples. Each point lies in exactly R quintuples.
"""

from collections import Counter
from itertools import combinations


def build(q, R):
    dual_vertices = [(x, y) for x in range(q) for y in range(R)]
    index = {xy: i for i, xy in enumerate(dual_vertices)}
    dual_edges = []
    for i in range(5):
        for b in range(q):
            dual_edges.append(tuple(index[((b + i * y) % q, y)] for y in range(R)))

    # Primal quintuples are dual vertices; primal points are dual edges.
    blocks = [set() for _ in dual_vertices]
    for point, edge in enumerate(dual_edges):
        for dual_vertex in edge:
            blocks[dual_vertex].add(point)
    return [tuple(sorted(block)) for block in blocks]


def verify(q, R, blocks):
    v = 5 * q
    n5 = R * q
    if len(blocks) != n5:
        return False, "wrong number of quintuples"
    if any(len(block) != 5 for block in blocks):
        return False, "some block does not have size 5"

    pts = sorted({x for block in blocks for x in block})
    if pts != list(range(v)):
        return False, f"point set mismatch: got {len(pts)} points"

    pair_seen = {}
    for block in blocks:
        for pair in combinations(block, 2):
            if pair in pair_seen:
                return False, f"pair {pair} repeated"
            pair_seen[pair] = block

    deg = Counter()
    for block in blocks:
        for x in block:
            deg[x] += 1
    deg_values = sorted({deg[x] for x in range(v)})
    if deg_values != [R]:
        return False, f"bad point degrees: {deg_values}"

    leave_degrees = [v - 1 - 4 * deg[x] for x in range(v)]
    if any(d % 3 for d in leave_degrees):
        return False, "leave degree not divisible by 3"
    leave_edges = v * (v - 1) // 2 - 10 * n5
    if leave_edges % 6:
        return False, "leave edge count not divisible by 6"

    return True, f"R={R}, q={q}, v={v}, n5={n5}, degree-values={deg_values}"


if __name__ == "__main__":
    for R in (5, 8, 11, 14, 17):
        for q in range(max(36, 4 * (R - 1) + 1), 181):
            if q % 12 not in (0, 9):
                continue
            blocks = build(q, R)
            ok, msg = verify(q, R, blocks)
            print(("OK   " if ok else "FAIL ") + msg)
            if not ok:
                raise SystemExit(1)
