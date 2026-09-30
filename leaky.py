"""
leaky.py -- Python/Sage front end for the brute-force leaky forcing checker lf.c.

Works with either a networkx Graph or a Sage Graph (anything with .vertices()/.edges()
in Sage's style, or .nodes()/.edges() in networkx's style).

    Z_leaky(G, ell)                -> Z^(ell)(G)           (exact; |V(G)| <= 64)
    is_leaky_forcing_set(G, ell, S) -> (True, [])  or (False, first_failing_leak_set)
                                                           (|V(G)| <= 1000, ell <= 2)

Conventions (as in the manuscript): a leak may be any vertex; a leaky vertex can be
forced but never forces. S is an ell-leaky forcing set iff it forces the whole graph for
EVERY placement of ell leaks. Z_leaky only searches supersets of the vertices of
degree <= ell (Lemma 2.2 of the manuscript / [3, Lemma 2.4]), which is a pure speed-up.

The C source is compiled automatically on first use (needs gcc or cc on the PATH).
"""
import os, subprocess, shutil

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "lf.c")
_BIN = os.path.join(_HERE, "lf")


def _binary():
    if not os.path.exists(_BIN) or os.path.getmtime(_BIN) < os.path.getmtime(_SRC):
        cc = shutil.which("gcc") or shutil.which("cc") or shutil.which("clang")
        if cc is None:
            raise RuntimeError("No C compiler found; compile lf.c by hand: gcc -O2 -o lf lf.c")
        subprocess.run([cc, "-O2", "-o", _BIN, _SRC], check=True)
    return _BIN


def _vertices_edges(G):
    """Return (vertex list, edge list of vertex pairs) for networkx or Sage graphs."""
    if hasattr(G, "nodes"):                      # networkx
        V = list(G.nodes())
        E = [(u, v) for u, v in G.edges()]
    else:                                        # Sage
        try:
            V = G.vertices(sort=False)
        except TypeError:                        # older Sage
            V = G.vertices()
        E = [(u, v) for u, v in G.edges(labels=False)]
    E = [(u, v) for u, v in E if u != v]         # simple graphs only
    return V, E


def _run(G, ell, mode, S=None):
    V, E = _vertices_edges(G)
    idx = {v: i for i, v in enumerate(V)}
    lines = [f"{len(V)} {len(E)} {ell} {mode}"]
    lines += [f"{idx[u]} {idx[v]}" for u, v in E]
    if mode == 1:
        lines.append(str(len(S)))
        lines.append(" ".join(str(idx[v]) for v in S))
    out = subprocess.run([_binary()], input="\n".join(lines) + "\n",
                         capture_output=True, text=True, check=True).stdout.split()
    return V, out


def Z_leaky(G, ell):
    """The ell-leaky forcing number Z^(ell)(G) by exhaustive search (|V| <= 64)."""
    V, _ = _vertices_edges(G)
    if len(V) > 64:
        raise ValueError("Z_leaky needs |V(G)| <= 64")
    _, out = _run(G, int(ell), 0)
    return int(out[0])


def is_leaky_forcing_set(G, ell, S):
    """Test whether S is an ell-leaky forcing set (ell <= 2, |V| <= 1000).
    Returns (True, []) or (False, leaks) where leaks is the first failing leak set
    ([] means S is not even a zero forcing set)."""
    if ell > 2:
        raise ValueError("is_leaky_forcing_set supports ell <= 2")
    V, out = _run(G, int(ell), 1, list(S))
    if out[0] == "OK":
        return True, []
    return False, [V[int(x)] for x in out[1:] if x != "none"]


# ---------- pure-Python reference (slow; for cross-checking the C code) ----------
def closure(G, S, leaks=()):
    """Blue set after exhaustively applying the leaky colour change rule."""
    V, E = _vertices_edges(G)
    nbr = {v: set() for v in V}
    for u, v in E:
        nbr[u].add(v); nbr[v].add(u)
    blue, leaks = set(S), set(leaks)
    changed = True
    while changed:
        changed = False
        for v in list(blue):
            if v in leaks:
                continue
            W = nbr[v] - blue
            if len(W) == 1:
                blue |= W; changed = True
    return blue


# ---------- generalized Petersen graphs in the manuscript's labelling ----------
def GP(n, k):
    """networkx P(n,k): outer ('y',i)~('y',i+1), spokes ('x',i)~('y',i), inner ('x',i)~('x',i+k); 1-indexed mod n."""
    import networkx as nx
    G = nx.Graph()
    for i in range(1, n + 1):
        G.add_edge(("y", i), ("y", i % n + 1))
        G.add_edge(("x", i), ("y", i))
        G.add_edge(("x", i), ("x", (i + k - 1) % n + 1))
    return G


def thm414_set(n, k):
    """S = union of the A_i from Section 4.2: outer vertices y_j with j in the first 2k+2
    positions of each block of 4k+4 consecutive indices (last block truncated)."""
    B, S = 4 * k + 4, []
    for i in range(-(-n // B)):
        S += [("y", j) for j in range(B * i + 1, B * i + 2 * k + 3) if j <= n]
    return S
