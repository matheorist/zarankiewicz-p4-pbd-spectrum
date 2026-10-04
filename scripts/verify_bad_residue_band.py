"""
Verify the bad-residue lower band construction from the manuscript.

For v == 0,9 mod 12 and N == 0 mod 3 with either mu_2(v) <= N <= v-15 or
N = v, the proof constructs N quintuple blocks whose point-degrees are in
{2,5}.  The leave is then K4-divisible and dense.
"""

from collections import Counter
from itertools import combinations


D5 = (0, 1, 3, 7, 12)


def mu2(v):
    t = (2 * v + 4) // 5
    while t % 3 != 0:
        t += 1
    return t


def cyclic_edges(n, h):
    return [tuple((s + d) % n for d in D5) for s in range(h)]


def havel_hakimi_bipartite(left_degrees, right_size, right_degree=5):
    rem = list(left_degrees)
    edges = []
    for b in range(right_size):
        choices = sorted(((-d, i) for i, d in enumerate(rem) if d > 0))[:right_degree]
        if len(choices) < right_degree:
            raise ValueError("not enough left degree")
        for _, i in choices:
            rem[i] -= 1
            edges.append((i, b))
    if any(rem):
        raise ValueError(f"left degrees remain: {rem}")
    return edges


def build(v, N):
    if N == v:
        blocks = [set() for _ in range(N)]
        for point in range(v):
            for d in D5:
                blocks[(point + d) % v].add(point)
        return [tuple(sorted(block)) for block in blocks]

    h = (5 * N - 2 * v) // 3
    a_size = (N + h) // 2
    b_size = (N - h) // 2
    f_edges = cyclic_edges(a_size, h)

    deg_f = [0] * a_size
    for edge in f_edges:
        for x in edge:
            deg_f[x] += 1
    left_degrees = [5 - d for d in deg_f]
    graph_edges = havel_hakimi_bipartite(left_degrees, b_size)

    # Dual vertices are future quintuples: A vertices 0..a_size-1 and
    # B vertices a_size..a_size+b_size-1.
    blocks = [set() for _ in range(N)]
    point = 0
    for edge in f_edges:
        for a in edge:
            blocks[a].add(point)
        point += 1
    for a, b in graph_edges:
        blocks[a].add(point)
        blocks[a_size + b].add(point)
        point += 1

    return [tuple(sorted(block)) for block in blocks]


def verify(v, N, blocks):
    if len(blocks) != N:
        return False, "wrong number of blocks"
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
    if any(d not in (2, 5) for d in deg_values):
        return False, f"bad point degrees: {deg_values}"
    if sum(deg.values()) != 5 * N:
        return False, "degree sum mismatch"

    leave_degrees = [v - 1 - 4 * deg[x] for x in range(v)]
    if any(d % 3 for d in leave_degrees):
        return False, "leave degree not divisible by 3"
    leave_edges = v * (v - 1) // 2 - 10 * N
    if leave_edges % 6:
        return False, "leave edge count not divisible by 6"

    return True, (
        f"v={v}, N={N}, point-degree-values={deg_values}, "
        f"min-leave-degree={min(leave_degrees)}"
    )


if __name__ == "__main__":
    for v in range(120, 241):
        if v % 12 not in (0, 9):
            continue
        candidates = list(range(mu2(v), v - 14)) + [v]
        for N in candidates:
            if N % 3:
                continue
            blocks = build(v, N)
            ok, msg = verify(v, N, blocks)
            print(("OK   " if ok else "FAIL ") + msg)
            if not ok:
                raise SystemExit(1)
