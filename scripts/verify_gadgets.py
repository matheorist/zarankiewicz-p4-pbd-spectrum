"""
Computer verification of the four small {4,5}-PBD gadgets used in the note
"Zarankiewicz numbers in the p=4 region".

A linear {4,5}-space on v points = a PBD(v,{4,5}): blocks of size 4 and 5 such
that every pair of points lies in exactly one block.

Gadgets required:
  16 points, block-type (n4,n5) = (15,3)   [Lemma B]
  16 points, block-type (n4,n5) = (10,6)   [Lemma B]
  25 points, block-type (n4,n5) = (5,27)   [Lemma B']
  25 points, block-type (n4,n5) = (10,24)  [Lemma B']

Each is constructed by randomized exact-cover backtracking (branch on the first
uncovered pair, randomized block order, restart on node-budget overflow) and then
independently verified.
"""

import sys
import random
from itertools import combinations

sys.setrecursionlimit(100000)


class Budget(Exception):
    pass


def attempt(v, n4, n5, seed, node_cap):
    rng = random.Random(seed)
    covered = [[False] * v for _ in range(v)]
    blocks = []
    budget = [n4, n5]
    nodes = [0]

    def first_uncovered():
        for i in range(v):
            row = covered[i]
            for j in range(i + 1, v):
                if not row[j]:
                    return (i, j)
        return None

    def mark(block, val):
        for x, y in combinations(block, 2):
            covered[x][y] = val

    def cliques(cands, need):
        res = []
        for combo in combinations(cands, need):
            if all(not covered[x][y] for x, y in combinations(combo, 2)):
                res.append(combo)
        return res

    def recurse():
        nodes[0] += 1
        if nodes[0] > node_cap:
            raise Budget()
        pair = first_uncovered()
        if pair is None:
            return budget[0] == 0 and budget[1] == 0
        a, b = pair
        cands = [c for c in range(v)
                 if c != a and c != b
                 and not covered[a][c] and not covered[b][c]]
        # build candidate blocks covering {a,b}
        cand_blocks = []
        for size, bi in ((4, 0), (5, 1)):
            if budget[bi] == 0:
                continue
            for extra in cliques(cands, size - 2):
                cand_blocks.append((bi, tuple(sorted((a, b) + extra))))
        rng.shuffle(cand_blocks)
        for bi, block in cand_blocks:
            mark(block, True)
            blocks.append(block)
            budget[bi] -= 1
            if recurse():
                return True
            budget[bi] += 1
            blocks.pop()
            mark(block, False)
        return False

    try:
        if recurse():
            return list(blocks), nodes[0]
    except Budget:
        return None, nodes[0]
    return None, nodes[0]


def solve(v, n4, n5, seeds=range(4000), node_cap=4000):
    assert 6 * n4 + 10 * n5 == v * (v - 1) // 2, "pair count mismatch"
    for seed in seeds:
        design, nodes = attempt(v, n4, n5, seed, node_cap)
        if design is not None:
            return design, seed, nodes
    return None, None, None


def verify(v, n4, n5, blocks):
    if blocks is None:
        return False, "no design found"
    s4 = sum(1 for b in blocks if len(b) == 4)
    s5 = sum(1 for b in blocks if len(b) == 5)
    if (s4, s5) != (n4, n5):
        return False, f"block-type {(s4, s5)} != {(n4, n5)}"
    if any(len(b) not in (4, 5) for b in blocks):
        return False, "block of illegal size"
    seen = {}
    for b in blocks:
        if len(set(b)) != len(b):
            return False, f"repeated point in block {b}"
        for x, y in combinations(sorted(b), 2):
            if (x, y) in seen:
                return False, f"pair {(x, y)} covered twice by {seen[(x,y)]} and {b}"
            seen[(x, y)] = b
    total = v * (v - 1) // 2
    if len(seen) != total:
        return False, f"only {len(seen)}/{total} pairs covered"
    return True, f"OK: {n4} quadruples + {n5} quintuples cover all {total} pairs exactly once"


def main():
    cases = [(16, 15, 3), (16, 10, 6), (25, 5, 27), (25, 10, 24)]
    all_ok = True
    for v, n4, n5 in cases:
        print(f"\n=== PBD({v},{{4,5}}) block-type (n4,n5)=({n4},{n5}) ===", flush=True)
        design, seed, nodes = solve(v, n4, n5)
        ok, msg = verify(v, n4, n5, design)
        print(f"  [seed={seed}, nodes={nodes}] {msg}", flush=True)
        if ok:
            for b in design:
                print("   ", b, flush=True)
        else:
            all_ok = False
    print("\n" + ("ALL FOUR GADGETS VERIFIED" if all_ok else "FAILURE"), flush=True)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
