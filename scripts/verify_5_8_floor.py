"""
For v=5,8 mod 12 (v=2 mod 3) every vertex of a {4,5}-space has r_5 = 1 mod 3 >= 1,
so n_5 >= v/5: the spectrum floor mu(v) = Theta(v), NOT O(1).
We verify the minimal construction at v=17 (=5 mod 12): mu(17)=ceil(17/5)=4.
Construction: center 0, 16 petals; 4 quintuples Q_i = {0} U (group i of 4 petals),
sharing only the center (a sunflower). Leave on the 16 petals = K_16 minus 4
disjoint K_4 = complete 4-partite K_{4,4,4,4}, which has a K_4-decomposition into
transversals via 2 MOLS(4) (= TD(4,4)). The union is a {4,5}-space with n_5=4.
"""
from itertools import combinations

def mols4():
    # 2 MOLS of order 4 over GF(4)? Use Z_2^2 / known pair. Simple known pair:
    L1 = [[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]]
    L2 = [[0,1,2,3],[2,3,0,1],[3,2,1,0],[1,0,3,2]]
    # check orthogonality
    pairs=set()
    for i in range(4):
        for j in range(4):
            pairs.add((L1[i][j],L2[i][j]))
    assert len(pairs)==16, "not orthogonal"
    return L1,L2

def build_17():
    center=0
    # petals 1..16 in 4 groups of 4: group g = {1+4g,...,4+4g}
    groups=[[1+4*g+r for r in range(4)] for g in range(4)]
    quints=[tuple([center]+groups[g]) for g in range(4)]
    # leave on petals = K_{4,4,4,4} with parts = groups; transversal K_4 via 2 MOLS(4)
    L1,L2=mols4()
    quads=[]
    # A transversal K_4: one petal from each group, indexed by (row r) with column via MOLS
    # Standard TD(4,4): blocks {(0,a),(1,b),(2,c),(3,d)} with b,c,d determined by a,row.
    # Use: for each cell (i,j) in 4x4, block = { group0[i], group1[j], group2[L1[i][j]], group3[L2[i][j]] }
    for i in range(4):
        for j in range(4):
            b=(groups[0][i],groups[1][j],groups[2][L1[i][j]],groups[3][L2[i][j]])
            quads.append(b)
    return quints, quads

def verify(v, quints, quads):
    seen={}
    for B in list(quints)+list(quads):
        for p in combinations(sorted(B),2):
            if p in seen: return False, f"pair {p} twice ({seen[p]} & {B})"
            seen[p]=B
    total=v*(v-1)//2
    if len(seen)!=total: return False, f"{len(seen)}/{total} pairs"
    n5=sum(len(B)==5 for B in quints+quads); n4=sum(len(B)==4 for B in quints+quads)
    if any(len(B) not in (4,5) for B in quints+quads): return False,"bad size"
    # replication check
    from collections import Counter
    r5=Counter(); r4=Counter()
    for B in quints:
        for x in B: r5[x]+=1
    for B in quads:
        for x in B: r4[x]+=1
    bad=[x for x in range(v) if 3*r4[x]+4*r5[x]!=v-1]
    return True, (f"OK v={v}: n5={n5}, n4={n4}, all pairs once; "
                  f"r5 values={sorted(set(r5[x] for x in range(v)))}, "
                  f"replication-ok={not bad}")

if __name__=="__main__":
    q5,q4=build_17()
    ok,msg=verify(17,q5,q4)
    print(msg if ok else "FAIL: "+msg)
    if ok:
        print(f"  => mu(17)=4=ceil(17/5); sub-v/5 values n5 in {{1}} are NOT realizable.")
        print(f"  quintuples: {q5}")
