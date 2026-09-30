# leaky_sage.sage -- Sage usage of the leaky forcing checker.
# Run from this folder:   sage leaky_sage.sage
# or inside a Sage session: load("leaky_sage.sage")
# (leaky.py accepts Sage Graph objects directly; lf.c is compiled on first use.)
import sys, os
sys.path.insert(0, os.getcwd())
from leaky import Z_leaky, is_leaky_forcing_set

def sage_GP_labels(n, k):
    """Sage's graphs.GeneralizedPetersenGraph(n,k) has outer vertices 0..n-1 (i ~ i+1),
    spokes i ~ n+i, inner n+i ~ n+(i+k mod n). Map the manuscript's y_i, x_i (1-indexed)."""
    y = lambda i: (i - 1) % n
    x = lambda i: n + (i - 1) % n
    return x, y

def thm414_set_sage(n, k):
    x, y = sage_GP_labels(n, k)
    B = 4*k + 4
    return [y(j) for i in range(ceil(n/B)) for j in range(B*i + 1, B*i + 2*k + 3) if j <= n]

# ---- checks from the referee report ----
P = graphs.PetersenGraph()
print("Petersen graph Z^(0), Z^(1), Z^(2):", [Z_leaky(P, l) for l in (0, 1, 2)])   # expect [5, 5, 5]

# Table 1 / Concern 3: P(7,2) and P(7,3) are isomorphic, so Z must agree
G72, G73 = graphs.GeneralizedPetersenGraph(7, 2), graphs.GeneralizedPetersenGraph(7, 3)
print("P(7,2) ~= P(7,3):", G72.is_isomorphic(G73),
      "  Z(P(7,2)) =", Z_leaky(G72, 0), "  Z(P(7,3)) =", Z_leaky(G73, 0))   # expect True, 6, 6

# Theorem 4.14 construction, e.g. k = 7, n = 87 (Figure 3)
n, k = 87, 7
G = graphs.GeneralizedPetersenGraph(n, k)
print("Thm 4.14 set 2-leaky for P(87,7):", is_leaky_forcing_set(G, 2, thm414_set_sage(n, k)))

# Concern 2: Theorem 5.2 first example, two stars with centres joined
for l in (1, 2, 3):
    for a in (l + 2, l + 1):                     # S_{l+2} as printed vs S_{l+1}
        T = graphs.StarGraph(a - 1).disjoint_union(graphs.StarGraph(a - 1))
        c1, c2 = (0, 0), (1, 0)                  # centres of the two stars
        G = T.copy(); G.add_edge(c1, c2)
        print(f"l={l}, S_{a}: Z(G) - Z(G-e) =", Z_leaky(G, l) - Z_leaky(T, l))
