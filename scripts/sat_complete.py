"""
SAT-decide: can the 12-quintuple 3-regular (20,5,1)-packing be completed to a
full PBD(v,{4,5}) with exactly n5=12 (rest quadruples)? I.e. does the leave graph
L = K_v minus the twelve K_5's have a K_4-decomposition?
If SAT for some admissible v, then n5=12 is realizable -> spectrum gap is tight
at the lower end (n5 = 0 or >=12, and 12 is attained).
"""
from itertools import combinations
from pysat.solvers import Cadical153 as Solver
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

QUINTS = [(11,4,17,12,14),(5,3,0,2,8),(3,15,4,7,13),(12,13,8,10,6),
          (3,14,18,6,16),(4,19,9,6,2),(2,13,17,18,1),(14,5,7,10,19),
          (9,12,1,15,5),(17,9,10,16,0),(11,0,19,15,18),(16,8,1,11,7)]

def try_v(v):
    covered = set()                      # pairs covered by quintuples
    for q in QUINTS:
        for p in combinations(sorted(q), 2):
            covered.add(p)
    # leave pairs = all pairs not covered by a quintuple
    leave = [p for p in combinations(range(v), 2) if p not in covered]
    leaveset = set(leave)
    # candidate 4-blocks: all 4-subsets whose 6 pairs are all in the leave
    cand = []
    for B in combinations(range(v), 4):
        if all((B[i], B[j]) in leaveset for i, j in combinations(range(4), 2)):
            cand.append(B)
    # SAT: var per candidate block; each leave pair covered by exactly one block
    pool = IDPool()
    var = {B: pool.id(('b', B)) for B in cand}
    pair_blocks = {p: [] for p in leave}
    for B in cand:
        for p in combinations(B, 2):
            pair_blocks[p].append(var[B])
    s = Solver()
    feasible = True
    for p in leave:
        lits = pair_blocks[p]
        if not lits:
            feasible = False
            break
        # exactly-one(lits)
        cnf = CardEnc.equals(lits=lits, bound=1, vpool=pool, encoding=EncType.seqcounter)
        for cl in cnf.clauses:
            s.add_clause(cl)
    if not feasible:
        return None
    ok = s.solve()
    if not ok:
        s.delete()
        return False
    model = set(l for l in s.get_model() if l > 0)
    blocks = [B for B in cand if var[B] in model]
    s.delete()
    return blocks

if __name__ == "__main__":
    for v in (25, 28, 37, 40):
        P = v*(v-1)//2
        if (P - 120) % 6 != 0:
            print(f"v={v}: pair count not divisible, skip", flush=True); continue
        print(f"v={v}: encoding & solving (n4 needed = {(P-120)//6}) ...", flush=True)
        res = try_v(v)
        if res is None:
            print(f"  v={v}: infeasible (a leave pair has no candidate block)", flush=True)
        elif res is False:
            print(f"  v={v}: UNSAT (no K4-decomposition of the leave)", flush=True)
        else:
            # verify full design
            allb = [tuple(sorted(q)) for q in QUINTS] + [tuple(sorted(b)) for b in res]
            seen = set(); ok = True
            for b in allb:
                for p in combinations(sorted(b), 2):
                    if p in seen: ok = False
                    seen.add(p)
            full = (len(seen) == P)
            print(f"  v={v}: SAT! full PBD({v},{{4,5}}) with n5=12, n4={len(res)}, "
                  f"pairs {'all covered once' if (ok and full) else 'ERROR'}", flush=True)
            if ok and full:
                print(f"     => n5=12 IS realizable; spectrum gap tight at lower end.", flush=True)
                break
