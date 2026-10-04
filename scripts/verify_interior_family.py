"""
Verify the interior exact-value construction:
  TD(5,q) truncated to a 5th group of size w; four groups of size q filled with
  S(2,4,q); the truncated group of size w filled with S(2,4,w).
Result should be a linear {4,5}-space on m=4q+w points with n5=q*w blocks of size 5
and total degree U_4 = (1/4)C(m,2)+(5/2)n.
"""
from itertools import combinations

def mols(q):
    # q-1 MOLS of order q for q prime: L_a[i][j] = (a*i+j) mod q, a=1..q-1
    return [[[(a*i + j) % q for j in range(q)] for i in range(q)] for a in range(1, q)]

def TD5(q):
    Ls = mols(q)
    blocks = []
    for i in range(q):
        for j in range(q):
            blocks.append([('G0', i), ('G1', j),
                           ('G2', Ls[0][i][j]), ('G3', Ls[1][i][j]), ('G4', Ls[2][i][j])])
    groups = {g: [(g, s) for s in range(q)] for g in ('G0','G1','G2','G3','G4')}
    return blocks, groups

def S24(points):
    pts = list(points); v = len(pts)
    cov = [[False]*v for _ in range(v)]
    out = []
    def first():
        for a in range(v):
            for b in range(a+1, v):
                if not cov[a][b]: return a, b
        return None
    def mark(t, val):
        for x, y in combinations(t, 2): cov[x][y] = val
    def rec():
        p = first()
        if p is None: return True
        a, b = p
        cand = [c for c in range(v) if c not in (a, b)
                and not cov[min(a,c)][max(a,c)] and not cov[min(b,c)][max(b,c)]]
        for ci in range(len(cand)):
            c = cand[ci]
            for d in cand[ci+1:]:
                if not cov[min(c,d)][max(c,d)]:
                    t = (a, b, c, d); mark(t, True); out.append(t)
                    if rec(): return True
                    out.pop(); mark(t, False)
        return False
    assert rec(), f"no S(2,4,{v})"
    return [tuple(pts[k] for k in t) for t in out]

def build(q, w):
    td, groups = TD5(q)
    deleted = set(('G4', s) for s in range(w, q))
    blocks = [tuple(p for p in blk if p not in deleted) for blk in td]
    for g in ('G0','G1','G2','G3'):
        blocks += S24(groups[g])
    blocks += S24([('G4', s) for s in range(w)])
    return blocks, 4*q + w

def verify(blocks, m, q, w):
    seen = {}; allpts = set()
    for b in blocks:
        allpts.update(b)
        if len(set(b)) != len(b): return False, f"repeat in {b}"
        for x, y in combinations(b, 2):
            key = tuple(sorted((str(x), str(y))))
            if key in seen: return False, f"pair {key} covered twice"
            seen[key] = b
    if len(allpts) != m: return False, f"{len(allpts)} pts != {m}"
    if len(seen) != m*(m-1)//2: return False, f"{len(seen)}/{m*(m-1)//2} pairs"
    if any(len(b) not in (4, 5) for b in blocks): return False, "bad block size"
    n5 = sum(len(b) == 5 for b in blocks); n4 = sum(len(b) == 4 for b in blocks)
    n = n4 + n5; z = 4*n4 + 5*n5
    U4num = m*(m-1)//2 + 10*n
    if U4num % 4: return False, "U4 not integer"
    U4 = U4num // 4
    if z != U4: return False, f"z={z} != U4={U4}"
    if n5 != q*w: return False, f"n5={n5} != qw={q*w}"
    return True, f"m={m}, n={n}, (n4,n5)=({n4},{n5}=qw), z=U4={z}=m(m-1)/3-5qw/3"

if __name__ == "__main__":
    for q, w in [(13, 13)]:
        blocks, m = build(q, w)
        ok, msg = verify(blocks, m, q, w)
        print(f"q={q} w={w} m={m} (m mod12={m%12}):", "PASS" if ok else "FAIL", msg, flush=True)
