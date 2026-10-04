"""
Joint SAT search for a PBD(v,{4,5}) with EXACTLY n5=12 quintuples.
Quintuples (n5=12, sum of r5 = 60, each vertex r5 in {0,3,6,...}) touch <= 20
vertices, so WLOG they lie within {0,...,19}; quadruples range over all v points.
Variables: 5-subsets of {0..19}, and 4-subsets of {0..v-1}.
Constraints: every pair covered exactly once; exactly 12 quintuples.
UNSAT => n5=12 not realizable on v points (definitive).
SAT   => explicit design; gap is tight (n5=12 attained).
"""
import sys, time
from itertools import combinations
from pysat.solvers import Cadical153 as Solver
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

def try_v(v, n5=12, qpts=20):
    t0 = time.time()
    pool = IDPool()
    quint = list(combinations(range(qpts), 5))
    quad  = list(combinations(range(v), 4))
    qv = {B: pool.id(('q5', B)) for B in quint}
    dv = {B: pool.id(('q4', B)) for B in quad}
    pair_blocks = {p: [] for p in combinations(range(v), 2)}
    for B in quint:
        for p in combinations(B, 2):
            pair_blocks[p].append(qv[B])
    for B in quad:
        for p in combinations(B, 2):
            pair_blocks[p].append(dv[B])
    s = Solver()
    ncl = 0
    for p, lits in pair_blocks.items():
        if not lits:
            s.delete(); return ('infeasible', p)
        cnf = CardEnc.equals(lits=lits, bound=1, vpool=pool, encoding=EncType.seqcounter)
        for cl in cnf.clauses:
            s.add_clause(cl); ncl += 1
    cnf = CardEnc.equals(lits=list(qv.values()), bound=n5, vpool=pool,
                         encoding=EncType.seqcounter)
    for cl in cnf.clauses:
        s.add_clause(cl); ncl += 1
    print(f"  v={v}: vars~{pool.top}, clauses~{ncl}, built in {time.time()-t0:.1f}s; solving...",
          flush=True)
    ok = s.solve()
    dt = time.time() - t0
    if not ok:
        s.delete(); return ('UNSAT', dt)
    model = set(l for l in s.get_model() if l > 0)
    q5 = [B for B in quint if qv[B] in model]
    q4 = [B for B in quad if dv[B] in model]
    s.delete()
    return ('SAT', q5, q4, dt)

def verify(v, q5, q4):
    seen = set()
    for B in list(q5) + list(q4):
        for p in combinations(sorted(B), 2):
            if p in seen: return False
            seen.add(p)
    return len(seen) == v*(v-1)//2 and len(q5) == 12

if __name__ == "__main__":
    vs = [int(x) for x in sys.argv[1:]] or [25]
    for v in vs:
        P = v*(v-1)//2
        if (P - 120) % 6 != 0:
            print(f"v={v}: not admissible", flush=True); continue
        res = try_v(v)
        tag = res[0]
        if tag == 'SAT':
            _, q5, q4, dt = res
            good = verify(v, q5, q4)
            print(f"v={v}: SAT in {dt:.1f}s  [verify {'OK' if good else 'FAIL'}]", flush=True)
            print(f"   n5=12 quintuples: {q5}", flush=True)
            print(f"   n4={len(q4)} quadruples: {q4}", flush=True)
            print(f"   ==> n5=12 REALIZABLE on v={v}; spectrum gap tight at lower end.", flush=True)
            break
        elif tag == 'UNSAT':
            print(f"v={v}: UNSAT in {res[1]:.1f}s  (no PBD(v,{{4,5}}) with n5=12)", flush=True)
        else:
            print(f"v={v}: {tag} {res[1]}", flush=True)
