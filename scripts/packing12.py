"""
Does the spectrum gap have minimum exactly 12? The value n5=12 requires a
3-regular (20,5,1)-packing: 12 quintuples on 20 points, every point in exactly 3,
any two quintuples meeting in <=1 point. We search for it directly (bitset).
If found, n5=12 is at least supported at the quintuple-substructure level.
We then also try to extend to a full PBD(v,{4,5}) with n5=12 for small admissible v.
"""
import random
from itertools import combinations

def regular_packing(v=20, k=5, r=3, time_budget_seeds=200000):
    """Find b = v*r/k blocks of size k on v points, each point in exactly r blocks,
    pairwise intersection <= 1."""
    b = v*r//k
    for seed in range(time_budget_seeds):
        rng = random.Random(seed)
        adj = [0]*v          # covered pairs
        deg = [0]*v
        blocks = []
        ok = True
        # greedy with restart: repeatedly pick a vertex of min remaining capacity
        for _ in range(b):
            # candidate block: 5 points, all pairwise uncovered, all with deg<r
            avail = [x for x in range(v) if deg[x] < r]
            rng.shuffle(avail)
            found = None
            # try to extend greedily
            def ext(block, pool):
                if len(block) == k:
                    return list(block)
                for idx, c in enumerate(pool):
                    if all(not (adj[c] >> y) & 1 for y in block):
                        res = ext(block+[c], pool[idx+1:])
                        if res: return res
                return None
            found = ext([], avail)
            if not found:
                ok = False; break
            for x, y in combinations(found, 2):
                adj[x] |= 1 << y; adj[y] |= 1 << x
            for x in found:
                deg[x] += 1
            blocks.append(tuple(found))
        if ok and all(deg[x] == r for x in range(v)) and len(blocks) == b:
            # verify pairwise <=1
            good = True
            for i in range(len(blocks)):
                for j in range(i+1, len(blocks)):
                    if len(set(blocks[i]) & set(blocks[j])) > 1:
                        good = False
            if good:
                return blocks, seed
    return None, None

blocks, seed = regular_packing()
if blocks:
    print(f"FOUND 3-regular (20,5,1)-packing (seed={seed}): 12 quintuples on 20 points")
    for b in blocks:
        print("  ", b)
    # check degrees
    from collections import Counter
    cnt = Counter()
    for b in blocks:
        for x in b: cnt[x]+=1
    print("  degrees:", dict(sorted(cnt.items())))
    print("  => the size-5 substructure for n5=12 EXISTS.")
else:
    print("no 3-regular (20,5,1)-packing found in budget")
