"""
Reproduces every computational claim in the Q402 follow-up referee report.
Run from the leaky_forcing folder:   python3 checks/run_all.py        (needs networkx)
Total runtime ~1-2 minutes; the Z^(2)(P(n,2)) sequence for n=10..13 takes the longest.
"""
import sys, os, random, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import networkx as nx
from networkx.generators.atlas import graph_atlas_g
from leaky import Z_leaky as Z, is_leaky_forcing_set, closure, GP, thm414_set

# ---------------------------------------------------------------- Section 3 formulas
def cycle_order(G):
    C = nx.cycle_basis(G)[0]; H = G.subgraph(C); order = [C[0]]
    while len(order) < len(C):
        for w in H[order[-1]]:
            if w not in order and (len(order) < 2 or w != order[-2]):
                order.append(w); break
    return order

def t36_as_stated(G):
    d = [G.degree(v) for v in cycle_order(G)]
    R, F, D = sum(x >= 4 for x in d), sum(x == 3 for x in d), sum(x == 2 for x in d)
    if R >= 2 or all(x >= 3 for x in d): return 0
    if F <= 1 and D >= 2: return 2
    return 1

def t36_as_proved(G):                      # t = 2 iff R empty and |F| <= 1
    d = [G.degree(v) for v in cycle_order(G)]
    R, F = sum(x >= 4 for x in d), sum(x == 3 for x in d)
    if R >= 2 or all(x >= 3 for x in d): return 0
    if R == 0 and F <= 1: return 2
    return 1

def t37(G):
    c = cycle_order(G); g = len(c); d = [G.degree(v) for v in c]
    hi = [x >= 3 for x in d]
    E = [i for i in range(g) if hi[i] and hi[(i + 1) % g]]
    disjoint = any((b - a) % g >= 2 and (a - b) % g >= 2 for a in E for b in E)
    cond2 = any(d[i] >= 4 and (d[i-1] >= 4 or d[(i+1) % g] >= 4 or (d[i-1] == 3 and d[(i+1) % g] == 3))
                for i in range(g))
    if disjoint or cond2: return 0
    F = [i for i in range(g) if d[i] == 3]
    if len(F) == 2 and (F[1] - F[0]) % g in (1, g - 1) and all(d[i] == 2 for i in range(g) if i not in F):
        return 2
    dist = lambda a, b: min((a - b) % g, (b - a) % g)
    if all(x <= 3 for x in d) and all(dist(a, b) >= 3 for a, b in itertools.combinations(F, 2)):
        return 2
    return 1

def t38(G):
    c = cycle_order(G); g = len(c); d = [G.degree(v) for v in c]
    F = [i for i in range(g) if d[i] == 3]
    if g in (3, 4) and all(x == 3 for x in d): return 2
    if g == 3 and len(F) == 2: return 1
    if g == 4 and ((len(F) == 2 and any((a - b) % 4 == 2 for a in F for b in F)) or len(F) == 3): return 1
    return 0

def check_unicyclic(G, t36):
    q = sum(1 for v in G if G.degree(v) == 1)
    r = sum(1 for v in G if G.degree(v) <= 2)
    g = len(cycle_order(G))
    res = {1: Z(G, 1) == q + (t36(G) if g == 3 else t37(G)),
           2: Z(G, 2) == r + t38(G)}
    for l in (3, 4):
        res[l] = Z(G, l) == sum(1 for v in G if G.degree(v) <= l)
    return res, g

def random_unicyclic(nv):
    g = random.randint(3, min(7, nv)); G = nx.cycle_graph(g)
    for v in range(g, nv): G.add_edge(v, random.randrange(v))
    return G

print("== Section 3: 4000 random unicyclic graphs on 4..13 vertices ==")
random.seed(1)
fails = {"3.6 as stated": 0, "3.6 as proved": 0, "3.7": 0, "3.8": 0, "3.9": 0}
for _ in range(4000):
    G = random_unicyclic(random.randint(4, 13))
    r_stated, g = check_unicyclic(G, t36_as_stated)
    r_proved, _ = check_unicyclic(G, t36_as_proved)
    if g == 3:
        fails["3.6 as stated"] += not r_stated[1]; fails["3.6 as proved"] += not r_proved[1]
    else:
        fails["3.7"] += not r_stated[1]
    fails["3.8"] += not r_stated[2]; fails["3.9"] += (not r_stated[3]) + (not r_stated[4])
print("  mismatches:", fails)
G = nx.Graph([(0, 1), (1, 2), (0, 2), (0, 3), (0, 4)])  # triangle + two leaves on c1
print("  counterexample to Thm 3.6 statement: Z1 =", Z(G, 1), " statement predicts", 2 + t36_as_stated(G))

print("== All connected graphs on <= 7 vertices: Thms 6.1, 6.4 and unicyclic formulas ==")
b61 = b64 = bu = 0
for G in graph_atlas_g()[1:]:
    n = G.number_of_nodes()
    if n < 2 or not nx.is_connected(G): continue
    z1 = Z(G, 1)
    path_or_cycle = max(d for _, d in G.degree()) <= 2
    b61 += (z1 == 2) != path_or_cycle
    if n >= 3:
        deg = sorted(d for _, d in G.degree()); m = G.number_of_edges()
        in_list = (m == n*(n-1)//2 or (m == n-1 and deg[-1] == n-1) or m == n*(n-1)//2 - 1 or
                   (deg[0] == 1 and m == (n-1)*(n-2)//2 + 1 and
                    nx.is_isomorphic(G.subgraph([v for v in G if G.degree(v) > 1]), nx.complete_graph(n-1))))
        b64 += (z1 == n - 1) != in_list
    if G.number_of_edges() == n:
        res, _ = check_unicyclic(G, t36_as_proved); bu += not all(res.values())
print(f"  failures: Thm 6.1 {b61}, Thm 6.4 {b64}, unicyclic (3.6 as proved, 3.7-3.9) {bu}")

# ---------------------------------------------------------------- Section 4
print("== Table 1 (rows: [Z, Z1, Z2]) ==")
for k, ns in ((2, range(5, 10)), (3, range(7, 12))):
    for n in ns:
        print(f"  P({n},{k}):", [Z(GP(n, k), l) for l in (0, 1, 2)])
print("  P(7,2) isomorphic to P(7,3):", nx.is_isomorphic(GP(7, 2), GP(7, 3)))
print("  Z2(P(n,1)), n=3..10:", [Z(GP(n, 1), 2) for n in range(3, 11)])

bad = [(n, k) for k in range(1, 9) for n in range(2*k + 1, 2*k + 40)
       if not is_leaky_forcing_set(GP(n, k), 1, [("x", i) for i in range(1, 2*k + 1)] + [("y", k), ("y", k + 1)])[0]]
print("== Thm 4.7 set (1-leaky), k=1..8: failures", bad)
bad = [(n, k) for k in range(1, 10) for n in range(6*k + 6, 10*k + 14)
       if not is_leaky_forcing_set(GP(n, k), 1, thm414_set(n, k))[0]]
print("== Lemma 4.12 (union of A_i is 1-leaky, n >= 6k+6), k=1..9: failures", bad)

print("== Thm 4.14 set (2-leaky) ==")
for k in range(1, 11):
    lo = 6*k + 6; hi = 10*k + 10 + 2*(4*k + 4) if k <= 8 else 10*k + 18
    fails = [(n, is_leaky_forcing_set(GP(n, k), 2, thm414_set(n, k))[1]) for n in range(lo, hi + 1)]
    fails = [f for f in fails if f[1]]
    print(f"  k={k:2d}, n={lo}..{hi}: {len(fails)} failing n"
          + (f"; e.g. n={fails[0][0]}, leaks {fails[0][1]}" if fails else ""))

print("== Z2(P(n,2)) growth, n=5..13 (slowest step) ==")
print("  ", [Z(GP(n, 2), 2) for n in range(5, 14)])

# ---------------------------------------------------------------- Section 5
print("== Thm 5.2 / 5.3 examples ==")
def double_star(a):                  # two copies of S_a (a vertices each), centres joined by e
    G = nx.Graph([("c1", "c2")])
    for i in range(a - 1): G.add_edge("c1", ("l1", i)); G.add_edge("c2", ("l2", i))
    return G
for l in (1, 2, 3):
    for a, name in ((l + 2, "S_{l+2} (as printed)"), (l + 1, "S_{l+1}")):
        G = double_star(a); H = G.copy(); H.remove_edge("c1", "c2")
        print(f"  l={l} {name}: Z(G)-Z(G-e) = {Z(G, l) - Z(H, l)}")
G = nx.path_graph(6); G.add_edge(1, 4)
print("  P6 + e, l=1: Z1(G) =", Z(G, 1), " Z1(P6) =", Z(nx.path_graph(6), 1))
def fig5(l):   # with the typo fixed: ceil(l/2) leaves on v12 and v22
    G = nx.complete_graph(["v11", "v12", "v21", "v22"]); fl, cl = l // 2, (l + 1) // 2
    for v, c in (("v11", fl), ("v21", fl), ("v12", cl), ("v22", cl)):
        for i in range(c): G.add_edge(v, (v, "leaf", i))
    for p in ("v1", "v2"):
        for i in range(fl): G.add_edge((p, "w", i), p + "1"); G.add_edge((p, "w", i), p + "2")
    return G
for l in (2, 3, 4, 5):
    G = fig5(l); H = G.copy(); H.remove_edge("v11", "v21")
    print(f"  Figure 5, l={l}: Z(G)-Z(G-e) = {Z(G, l) - Z(H, l)}")
for l in (1, 2, 3):
    H1 = nx.Graph()
    for c in (0, 1):
        K = [(c, i) for i in range(l)]
        for i in range(l):
            H1.add_edge("v", K[i]); H1.add_edges_from((K[i], K[j]) for j in range(i + 1, l))
    H2 = H1.copy(); H2.add_edge("v", "w")
    T = nx.Graph()
    for i in range(l + 1):
        T.add_edge("r", ("c", i)); T.add_edges_from((("c", i), ("c", i, j)) for j in range(l))
    T2 = T.copy(); T2.remove_node("r")
    print(f"  Thm 5.3, l={l}: upper example diff {Z(H2, l) - Z(H1, l)} (want 1);"
          f" lower example diff {Z(T, l) - Z(T2, l)} (want -{l + 1})")

# ---------------------------------------------------------------- cross-check C against pure Python
random.seed(7); disagree = 0
for _ in range(300):
    G = nx.gnp_random_graph(9, 0.35); V = list(G)
    S = random.sample(V, random.randint(1, 5)); L = random.sample(V, random.randint(0, 2))
    ok_c = is_leaky_forcing_set(G, len(L), S)[0]
    ok_py = all(len(closure(G, S, LL)) == len(V) for j in range(len(L) + 1) for LL in itertools.combinations(V, j))
    disagree += ok_c != ok_py
print("== C vs pure-Python cross-check on 300 random instances: disagreements", disagree)
