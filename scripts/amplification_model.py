"""
Asymptotic bookkeeping model for the truncated-TD amplification.

This is not a proof.  It records the residue arithmetic behind the heuristic
in the manuscript:

    beta -> 4/5 + beta/5

The model treats the base band for a group of size s as an interval of admissible
quintuple counts from its lower floor window (O(1) for s=1,4 mod 12 and O(s) for
s=5,8 mod 12) up to beta*M(s), where M(s)=binom(s,2)/10.  It then enumerates
decompositions m=4q+w with 1<=w<=q and q,w in the residue classes where the
current base bands are available.
"""

from collections import defaultdict


BASE_CLASSES = {1, 4, 5, 8}  # mod 12, covered by band2 and banddual near threshold
TARGET_CLASSES = {0, 1, 4, 5, 8, 9}  # m == 0 or 1 mod 4


def M(v):
    return v * (v - 1) / 20.0


def interval_ratio(m, q, w, beta):
    # Asymptotic upper endpoint:
    # cross TD contributes q*w quintuples; five group fillers contribute beta*M.
    high = q * w + beta * (4 * M(q) + M(w))
    return high / M(m)


def scan(beta=0.01, lo=500, hi=2000):
    best = defaultdict(lambda: (0.0, None))
    for m in range(lo, hi + 1):
        if m % 12 not in TARGET_CLASSES:
            continue
        for q in range(1, m):
            w = m - 4 * q
            if not (1 <= w <= q):
                continue
            if q % 12 not in BASE_CLASSES or w % 12 not in BASE_CLASSES:
                continue
            r = interval_ratio(m, q, w, beta)
            key = m % 12
            if r > best[key][0]:
                best[key] = (r, (m, q, w, q % 12, w % 12))
    return best


def iterate(beta, rounds):
    vals = [beta]
    for _ in range(rounds):
        beta = 4 / 5 + beta / 5
        vals.append(beta)
    return vals


if __name__ == "__main__":
    print("beta iteration:")
    for i, b in enumerate(iterate(0.001, 8)):
        print(f"  round {i}: beta={b:.6f}")

    print("\nResidue scan using beta=0.001:")
    best = scan(beta=0.001)
    for cls in sorted(TARGET_CLASSES):
        r, data = best[cls]
        print(f"  m mod 12 = {cls:2d}: max ratio {r:.6f} via {data}")

    print("\nResidue scan using beta after one amplification:")
    beta1 = 4 / 5 + 0.001 / 5
    best = scan(beta=beta1)
    for cls in sorted(TARGET_CLASSES):
        r, data = best[cls]
        print(f"  m mod 12 = {cls:2d}: max ratio {r:.6f} via {data}")
