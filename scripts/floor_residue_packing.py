"""
Construct and verify low-degree 5-uniform linear packings for v == 2 (mod 3).

For a PBD(v,{4,5}) in the residue classes v == 5,8 (mod 12), the leave after
the quintuple packing is K4-divisible exactly when every vertex has
    r5(x) == 1 (mod 3).

This script builds a packing with degrees in {1,4}: choose h in {0,...,4} such
that v + 3h is divisible by 5.  Then h "high" vertices have degree 4 and all
other vertices have degree 1.  The construction keeps high-high pairs repeated
at most once, and fills the remaining incidences by fresh low vertices.
"""

from itertools import combinations


def high_count(v):
    for h in range(5):
        if (v + 3 * h) % 5 == 0:
            return h
    raise AssertionError("unreachable")


def high_patterns(h):
    highs = list(range(h))
    if h == 0:
        return []
    if h == 1:
        return [(0,), (0,), (0,), (0,)]
    if h == 2:
        return [(0, 1), (0,), (0,), (0,), (1,), (1,), (1,)]
    if h == 3:
        return [(0, 1, 2)] + [(0,), (0,), (0,), (1,), (1,), (1,), (2,), (2,), (2,)]
    if h == 4:
        return [(0, 1, 2, 3)] + [
            (0,), (0,), (0,), (1,), (1,), (1,),
            (2,), (2,), (2,), (3,), (3,), (3,)
        ]
    raise ValueError(h)


def build(v):
    h = high_count(v)
    patterns = high_patterns(h)
    n5 = (v + 3 * h) // 5
    if n5 < len(patterns):
        raise ValueError(f"v={v} is too small for h={h}: n5={n5}, patterns={len(patterns)}")

    lows = iter(range(h, v))
    blocks = []
    for pat in patterns:
        need = 5 - len(pat)
        block = list(pat) + [next(lows) for _ in range(need)]
        blocks.append(tuple(block))

    for _ in range(n5 - len(patterns)):
        blocks.append(tuple(next(lows) for _ in range(5)))

    return h, blocks


def verify(v, blocks):
    deg = [0] * v
    seen = {}
    for b in blocks:
        if len(b) != 5 or len(set(b)) != 5:
            return False, f"bad block {b}"
        for x in b:
            deg[x] += 1
        for p in combinations(sorted(b), 2):
            if p in seen:
                return False, f"pair {p} repeated in {seen[p]} and {b}"
            seen[p] = b

    if any(d % 3 != 1 for d in deg):
        return False, f"bad residues: {sorted(set(deg))}"
    if sum(deg) != 5 * len(blocks):
        return False, "degree sum mismatch"
    if len(blocks) % 3 != 1:
        return False, f"n5={len(blocks)} is not 1 mod 3"

    leave_degrees = [v - 1 - 4 * d for d in deg]
    if any(x % 3 != 0 for x in leave_degrees):
        return False, f"leave degrees not 0 mod 3: {sorted(set(leave_degrees))}"
    leave_edges = v * (v - 1) // 2 - 10 * len(blocks)
    if leave_edges % 6 != 0:
        return False, f"leave edge count {leave_edges} not divisible by 6"

    return True, (
        f"v={v}, h={sum(1 for d in deg if d == 4)}, n5={len(blocks)}, "
        f"r5-values={sorted(set(deg))}, min leave degree={min(leave_degrees)}"
    )


if __name__ == "__main__":
    for v in range(17, 121):
        if v % 12 not in (5, 8):
            continue
        h, blocks = build(v)
        ok, msg = verify(v, blocks)
        print(("OK   " if ok else "FAIL ") + msg)
        if not ok:
            raise SystemExit(1)
