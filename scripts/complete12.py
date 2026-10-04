"""
Try to complete the 3-regular (20,5,1)-packing of 12 quintuples to a full
PBD(v,{4,5}) with exactly n5=12, by covering all remaining pairs with quadruples.
We try v = 25, 28, 37, 40 (admissible: v = 1,4 mod 12).
"""
import random
from itertools import combinations

QUINTS = [(11,4,17,12,14),(5,3,0,2,8),(3,15,4,7,13),(12,13,8,10,6),
          (3,14,18,6,16),(4,19,9,6,2),(2,13,17,18,1),(14,5,7,10,19),
          (9,12,1,15,5),(17,9,10,16,0),(11,0,19,15,18),(16,8,1,11,7)]

def complete(v, seeds=300000, node_cap=40000):
    FULL = (1<<v)-1
    base = [0]*v
    for q in QUINTS:
        for x,y in combinations(q,2):
            base[x] |= 1<<y; base[y] |= 1<<x
    # required quadruples count
    P = v*(v-1)//2
    n4 = (P - 12*10)//6
    if (P - 120) % 6 != 0:
        return None, f"v={v}: pair count not divisible"
    for seed in range(seeds):
        rng = random.Random(seed)
        adj = list(base); need=[n4]; blocks=[]; nodes=[0]
        def rec():
            nodes[0]+=1
            if nodes[0]>node_cap: raise StopIteration
            a=-1;b=-1
            for i in range(v):
                un=(~adj[i])&FULL&~(1<<i)&~((1<<(i+1))-1)
                if un:
                    a=i;b=(un&-un).bit_length()-1;break
            if a==-1:
                return need[0]==0
            if need[0]==0:
                return False  # uncovered pair but no quads left
            avail=(~adj[a])&(~adj[b])&FULL&~(1<<a)&~(1<<b)
            cands=[c for c in range(v) if (avail>>c)&1]
            cb=[]
            for ci in range(len(cands)):
                c=cands[ci]
                for d in cands[ci+1:]:
                    if not (adj[c]>>d)&1:
                        cb.append((a,b,c,d))
            rng.shuffle(cb)
            for blk in cb:
                for x in range(4):
                    for y in range(x+1,4):
                        adj[blk[x]]|=1<<blk[y]; adj[blk[y]]|=1<<blk[x]
                need[0]-=1; blocks.append(blk)
                if rec(): return True
                need[0]+=1; blocks.pop()
                for x in range(4):
                    for y in range(x+1,4):
                        adj[blk[x]]&=~(1<<blk[y]); adj[blk[y]]&=~(1<<blk[x])
            return False
        try:
            if rec():
                return blocks, seed
        except StopIteration:
            pass
    return None, None

for v in (25, 28, 37, 40):
    res, seed = complete(v, seeds=20000, node_cap=30000)
    if isinstance(seed, str):
        print(seed); continue
    if res:
        # verify full design
        allb = [tuple(q) for q in QUINTS] + [tuple(b) for b in res]
        seen=set(); okv=True
        for bb in allb:
            for x,y in combinations(sorted(bb),2):
                if (x,y) in seen: okv=False
                seen.add((x,y))
        full = (len(seen)==v*(v-1)//2)
        print(f"v={v}: FOUND full PBD({v},{{4,5}}) with n5=12, n4={len(res)} "
              f"[verify pairs={'OK' if full and okv else 'FAIL'}, seed={seed}]")
        break
    else:
        print(f"v={v}: not completed in budget")
