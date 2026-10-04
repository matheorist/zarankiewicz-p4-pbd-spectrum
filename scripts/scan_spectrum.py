"""
Fast bitset exact-cover search for PBD(v,{4,5}) with prescribed block-type.
Scans which (v, n5) are realizable for the small orders relevant to the note,
to determine the true small-design spectrum.
"""
import random, time

def search(v, n4, n5, time_budget=4.0, node_cap=200000):
    assert 6*n4 + 10*n5 == v*(v-1)//2
    FULL = (1 << v) - 1
    t0 = time.time()

    def attempt(rng):
        adj = [0]*v          # adj[i] bit j set if pair covered
        need = [n4, n5]
        blocks = []
        nodes = 0
        def rec():
            nonlocal nodes
            nodes += 1
            if nodes > node_cap:
                return "cap"
            # first uncovered pair
            a = -1; b = -1
            for i in range(v):
                un = (~adj[i]) & FULL & ~(1 << i)
                un &= ~((1 << (i+1)) - 1)   # only j > i
                if un:
                    a = i
                    b = (un & -un).bit_length() - 1
                    break
            if a == -1:
                return need[0] == 0 and need[1] == 0
            avail = (~adj[a]) & (~adj[b]) & FULL & ~(1 << a) & ~(1 << b)
            cands = [c for c in range(v) if (avail >> c) & 1]
            cblocks = []
            if need[0] > 0:
                for ci in range(len(cands)):
                    c = cands[ci]
                    for d in cands[ci+1:]:
                        if not (adj[c] >> d) & 1:
                            cblocks.append((0, (a, b, c, d)))
            if need[1] > 0:
                for ci in range(len(cands)):
                    c = cands[ci]
                    for di in range(ci+1, len(cands)):
                        d = cands[di]
                        if (adj[c] >> d) & 1:
                            continue
                        for e in cands[di+1:]:
                            if not (adj[c] >> e) & 1 and not (adj[d] >> e) & 1:
                                cblocks.append((1, (a, b, c, d, e)))
            rng.shuffle(cblocks)
            for bi, blk in cblocks:
                for x in range(len(blk)):
                    for y in range(x+1, len(blk)):
                        adj[blk[x]] |= (1 << blk[y]); adj[blk[y]] |= (1 << blk[x])
                need[bi] -= 1; blocks.append(blk)
                r = rec()
                if r is True:
                    return True
                if r == "cap":
                    return "cap"
                need[bi] += 1; blocks.pop()
                for x in range(len(blk)):
                    for y in range(x+1, len(blk)):
                        adj[blk[x]] &= ~(1 << blk[y]); adj[blk[y]] &= ~(1 << blk[x])
            return False
        r = rec()
        return (list(blocks) if r is True else None)

    seed = 0
    while time.time() - t0 < time_budget:
        rng = random.Random(seed); seed += 1
        res = attempt(rng)
        if res is not None:
            return res, seed, time.time()-t0
    return None, seed, time.time()-t0

def repl_feasible(v, n5):
    """Is there a nonneg (r4,r5) with 3 r4 + 4 r5 = v-1 and r5<=n5 for some vertex
    distribution? Returns the set of allowed r5 values per vertex."""
    allowed = []
    for r5 in range(0, (v-1)//4 + 1):
        if (v-1-4*r5) % 3 == 0 and (v-1-4*r5) >= 0:
            allowed.append(r5)
    return allowed

def verify(v, n4, n5, blocks):
    if blocks is None: return False
    from itertools import combinations
    seen=set()
    for b in blocks:
        for x,y in combinations(sorted(b),2):
            if (x,y) in seen: return False
            seen.add((x,y))
    return len(seen)==v*(v-1)//2 and sum(len(b)==4 for b in blocks)==n4 and sum(len(b)==5 for b in blocks)==n5

print("vertex r5-values allowed by replication 3r4+4r5=v-1:")
for v in (16, 25, 28, 37, 40):
    print(f"  v={v}: r5 in {repl_feasible(v, None)}")

print("\n--- v=16 spectrum (n5 = 0,3,6,9,12) ---")
for n5 in (0,3,6,9,12):
    n4 = (16*15//2 - 10*n5)//6
    if (16*15//2 - 10*n5) % 6 != 0:
        print(f"  n5={n5}: n4 non-integer, skip"); continue
    blk, tries, dt = search(16, n4, n5, time_budget=6.0)
    ok = verify(16, n4, n5, blk)
    print(f"  n5={n5:2d} (n4={n4:2d}): {'FOUND' if ok else 'not found'}  [{tries} tries, {dt:.1f}s]")

print("\n--- v=25 spectrum (n4 = 0,5,10,15,20 ; n5 = 30,27,24,21,18) ---")
for n4 in (0,5,10,15,20):
    n5 = (25*24//2 - 6*n4)//10
    if (25*24//2 - 6*n4) % 10 != 0:
        print(f"  n4={n4}: n5 non-integer, skip"); continue
    blk, tries, dt = search(25, n4, n5, time_budget=6.0)
    ok = verify(25, n4, n5, blk)
    print(f"  n4={n4:2d} (n5={n5:2d}): {'FOUND' if ok else 'not found'}  [{tries} tries, {dt:.1f}s]")
