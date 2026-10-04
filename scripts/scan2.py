"""Focused search: smallest n5>0 realizable for v=25,28 (v=1,4 mod 12).
Lower node_cap => more restarts (restarts beat deep backtracking here)."""
import random, time
from itertools import combinations

def search(v, n4, n5, time_budget, node_cap):
    FULL = (1 << v) - 1
    t0 = time.time(); seed = 0
    while time.time() - t0 < time_budget:
        rng = random.Random(seed); seed += 1
        adj = [0]*v; need=[n4,n5]; blocks=[]; nodes=[0]
        def rec():
            nodes[0]+=1
            if nodes[0] > node_cap: raise StopIteration
            a=-1;b=-1
            for i in range(v):
                un=(~adj[i])&FULL&~(1<<i)&~((1<<(i+1))-1)
                if un:
                    a=i; b=(un&-un).bit_length()-1; break
            if a==-1: return need[0]==0 and need[1]==0
            avail=(~adj[a])&(~adj[b])&FULL&~(1<<a)&~(1<<b)
            cands=[c for c in range(v) if (avail>>c)&1]
            cb=[]
            if need[0]>0:
                for ci in range(len(cands)):
                    c=cands[ci]
                    for d in cands[ci+1:]:
                        if not (adj[c]>>d)&1: cb.append((0,(a,b,c,d)))
            if need[1]>0:
                for ci in range(len(cands)):
                    c=cands[ci]
                    for di in range(ci+1,len(cands)):
                        d=cands[di]
                        if (adj[c]>>d)&1: continue
                        for e in cands[di+1:]:
                            if not (adj[c]>>e)&1 and not (adj[d]>>e)&1:
                                cb.append((1,(a,b,c,d,e)))
            rng.shuffle(cb)
            for bi,blk in cb:
                for x in range(len(blk)):
                    for y in range(x+1,len(blk)):
                        adj[blk[x]]|=1<<blk[y]; adj[blk[y]]|=1<<blk[x]
                need[bi]-=1; blocks.append(blk)
                if rec(): return True
                need[bi]+=1; blocks.pop()
                for x in range(len(blk)):
                    for y in range(x+1,len(blk)):
                        adj[blk[x]]&=~(1<<blk[y]); adj[blk[y]]&=~(1<<blk[x])
            return False
        try:
            if rec():
                # verify
                seen=set()
                for bb in blocks:
                    for x,y in combinations(sorted(bb),2):
                        assert (x,y) not in seen; seen.add((x,y))
                assert len(seen)==v*(v-1)//2
                return blocks, seed, time.time()-t0
        except StopIteration:
            pass
    return None, seed, time.time()-t0

for v in (25, 28):
    P = v*(v-1)//2
    print(f"\n=== v={v} (P={P} pairs) : smallest n5>0 ? ===", flush=True)
    found_any=False
    for n5 in range(12, P//10+1, 3):
        if (P-10*n5)%6: continue
        n4=(P-10*n5)//6
        if n4<0: break
        blk,tries,dt=search(v,n4,n5,time_budget=12.0,node_cap=20000)
        st = "FOUND" if blk else "not found"
        print(f"  n5={n5:2d} (n4={n4:2d}): {st}  [{tries} restarts, {dt:.1f}s]", flush=True)
        if blk:
            found_any=True
            print(f"      blocks: {[tuple(b) for b in blk]}", flush=True)
            break   # smallest found; stop
    if not found_any:
        print("  (no n5>0 found in budget)", flush=True)
