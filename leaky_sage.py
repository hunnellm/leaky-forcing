"""
leaky_sage.py -- pure-Sage implementation of the leaky forcing tools in leaky.py
(no C compiler or lf.c needed).

Load directly from GitHub in Sage:

    sage: load("https://raw.githubusercontent.com/hunnellm/leaky-forcing/main/leaky_sage.py")
    sage: G = graphs.PetersenGraph()
    sage: Z_leaky(G, 1)
    sage: is_leaky_forcing_set(G, 1, [0, 1, 2, 3])

    Z_leaky(G, ell)                 -> Z^(ell)(G), exact, by exhaustive search
    Z_leaky(G, ell, return_sets=True) -> all minimum-size ell-leaky forcing sets
    is_leaky_forcing_set(G, ell, S) -> (True, []) or (False, first_failing_leak_set)
                                       ([] means S is not even a zero forcing set)

Conventions: a leak may be any vertex; a leaky vertex can be forced but never forces.
S is an ell-leaky forcing set iff it forces the whole graph for EVERY placement of
ell leaks (more leaks only make forcing harder, so exactly min(ell, n) leaks suffice).
Z_leaky only searches supersets of the vertices of degree <= ell (Lemma 2.2).
"""
from itertools import combinations
from sage.all import Graph


def _setup(G):
    """Return (vertex list, neighbour bitmasks) for a Sage graph."""
    V = list(G.vertices(sort=False))
    idx = {v: i for i, v in enumerate(V)}
    nbr = [0] * len(V)
    for u, v in G.edges(labels=False):
        if u != v:
            nbr[idx[u]] |= 1 << idx[v]
            nbr[idx[v]] |= 1 << idx[u]
    return V, idx, nbr


def _closure_mask(nbr, blue, leaks):
    """Blue bitmask after exhaustively applying the leaky colour change rule."""
    blue = int(blue); leaks = int(leaks)
    n = len(nbr)
    changed = True
    while changed:
        changed = False
        for v in range(n):
            if (blue >> v) & 1 and not (leaks >> v) & 1:
                w = nbr[v] & ~blue
                if w and not (w & (w - 1)):      # exactly one white neighbour
                    blue |= w
                    changed = True
    return blue


def _first_failing_leaks(nbr, S, ell):
    """None if S is ell-leaky forcing, else the first failing tuple of leak indices."""
    n = len(nbr)
    full = (1 << n) - 1
    for leaks in combinations(range(n), min(ell, n)):
        lm = 0
        for x in leaks:
            lm |= 1 << x
        if _closure_mask(nbr, S, lm) != full:
            return leaks
    return None


def is_leaky_forcing_set(G, ell, S):
    """Test whether S is an ell-leaky forcing set.
    Returns (True, []) or (False, leaks); leaks == [] means S is not a zero forcing set."""
    V, idx, nbr = _setup(G)
    Sm = 0
    for v in S:
        Sm |= 1 << idx[v]
    full = (1 << len(V)) - 1
    if _closure_mask(nbr, Sm, 0) != full:
        return False, []
    bad = _first_failing_leaks(nbr, Sm, int(ell))
    if bad is None:
        return True, []
    return False, [V[i] for i in bad]


def Z_leaky(G, ell, return_sets=False):
    """Return Z^(ell)(G), or all minimum-size ell-leaky forcing sets when requested.

    With return_sets=True, returns a list of sets of Sage graph vertices.
    """
    ell = int(ell)
    V, idx, nbr = _setup(G)
    n = len(V)
    forced = 0
    for i in range(n):
        if bin(nbr[i]).count("1") <= ell:
            forced |= 1 << i
    base = bin(forced).count("1")
    rest = [i for i in range(n) if not (forced >> i) & 1]
    for extra in range(len(rest) + 1):
        forcing_sets = []
        for add in combinations(rest, extra):
            Sm = forced
            for i in add:
                Sm |= 1 << i
            if _first_failing_leaks(nbr, Sm, ell) is None:
                if not return_sets:
                    return base + extra
                forcing_sets.append({V[i] for i in range(n) if (Sm >> i) & 1})
        if forcing_sets:
            return forcing_sets
    return [] if return_sets else n


def closure(G, S, leaks=()):
    """Blue set after exhaustively applying the leaky colour change rule."""
    V, idx, nbr = _setup(G)
    Sm = sum(1 << idx[v] for v in set(S))
    Lm = sum(1 << idx[v] for v in set(leaks))
    b = _closure_mask(nbr, Sm, Lm)
    return {V[i] for i in range(len(V)) if (b >> i) & 1}


# ---------- generalized Petersen graphs in the manuscript's labelling ----------
def GP(n, k):
    """P(n,k) as a Sage Graph: outer ('y',i)~('y',i+1), spokes ('x',i)~('y',i),
    inner ('x',i)~('x',i+k); 1-indexed mod n."""
    G = Graph()
    for i in range(1, n + 1):
        G.add_edge(("y", i), ("y", i % n + 1))
        G.add_edge(("x", i), ("y", i))
        G.add_edge(("x", i), ("x", (i + k - 1) % n + 1))
    return G


def thm414_set(n, k):
    """Union of the A_i from Section 4.2: outer vertices y_j with j in the first 2k+2
    positions of each block of 4k+4 consecutive indices (last block truncated)."""
    B, S = 4 * k + 4, []
    for i in range(-(-n // B)):
        S += [("y", j) for j in range(B * i + 1, B * i + 2 * k + 3) if j <= n]
    return S
