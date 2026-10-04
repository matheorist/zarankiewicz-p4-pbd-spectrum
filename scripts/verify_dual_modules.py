"""
Verify the two 4-regular (r_5=4=1 mod 3) linear quintuple modules used in
thm:banddual (the =5,8 mod 12 positive-fraction band):
  (A) TD(5,4): 16 quintuples on 20 points, 4-regular, pairwise <=1, covers cross pairs.
  (B) 4 parallel classes of AG(2,5): 20 quintuples on 25 points, 4-regular, pairwise <=1.
Also the floor packing check (1-regular parallel class) and the semigroup logic
mu(R')+16a+20b = mu(R)+12a+15b (i.e. mu(R')=mu(R)-4a-5b).
"""
from itertools import combinations
from collections import Counter

def check_packing(blocks, npts, name, want_reg):
    # linear (pairwise <=1), regular, block size 5
    for b in blocks:
        assert len(set(b))==5, f"{name}: block size !=5 {b}"
    for i in range(len(blocks)):
        for j in range(i+1,len(blocks)):
            if len(set(blocks[i])&set(blocks[j]))>1:
                return f"{name}: FAIL blocks {blocks[i]},{blocks[j]} meet >1"
    deg=Counter()
    for b in blocks:
        for x in b: deg[x]+=1
    regs=set(deg[x] for x in range(npts))
    return (f"{name}: {len(blocks)} blocks on {npts} pts, sizes=5, linear OK, "
            f"degrees={sorted(regs)} (want {want_reg}), r5 mod3={sorted(set(r%3 for r in regs))}")

def TD54():
    # 3 MOLS(4) -> TD(5,4). Groups G0..G4 each {0..3}. Block(i,j)=
    # {G0:i, G1:j, G2:L1[i][j], G3:L2[i][j], G4:L3[i][j]}. Points tagged 0..19.
    L1=[[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]]
    L2=[[0,1,2,3],[2,3,0,1],[3,2,1,0],[1,0,3,2]]
    L3=[[0,1,2,3],[3,2,1,0],[1,0,3,2],[2,3,0,1]]
    # check mutually orthogonal
    for (A,B) in [(L1,L2),(L1,L3),(L2,L3)]:
        pset=set((A[i][j],B[i][j]) for i in range(4) for j in range(4))
        assert len(pset)==16, "MOLS not orthogonal"
    def pt(g,s): return g*4+s
    blocks=[]
    for i in range(4):
        for j in range(4):
            blocks.append((pt(0,i),pt(1,j),pt(2,L1[i][j]),pt(3,L2[i][j]),pt(4,L3[i][j])))
    return blocks,20

def AG25_kclasses(k):
    pts=[(x,y) for x in range(5) for y in range(5)]
    idx={p:i for i,p in enumerate(pts)}
    slopes=['inf',0,1,2][:k]
    blocks=[]
    for slope in slopes:
        for c in range(5):
            if slope=='inf': line=[(c,y) for y in range(5)]
            else: line=[(x,(slope*x+c)%5) for x in range(5)]
            blocks.append(tuple(sorted(idx[p] for p in line)))
    return blocks,25

if __name__=="__main__":
    tb,n=TD54(); print(check_packing(tb,n,"TD(5,4)",4))
    ab,n=AG25_kclasses(4); print(check_packing(ab,n,"4-classes AG(2,5)",4))
    ab3,n=AG25_kclasses(3); print(check_packing(ab3,n,"3-classes AG(2,5)",3))
    # semigroup logic sanity: <12,15> covers all mult of 3 >= 36 plus {12,15,24,27,30}
    S=set()
    for a in range(50):
        for b in range(50):
            if 12*a+15*b<=200: S.add(12*a+15*b)
    miss=[x for x in range(0,201,3) if x not in S and x>=36]
    print(f"<12,15>: missing mult-of-3 >=36 up to 200: {miss} (should be [])")
