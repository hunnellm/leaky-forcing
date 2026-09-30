# Leaky forcing checker (Q402 referee checks)

Brute-force computation of ℓ-leaky forcing numbers, used for the computational claims in the Q402 follow-up referee report.

## Files

- `lf.c`: C core. Exact Z^(ℓ)(G) for graphs on up to 64 vertices. Also tests whether a given set is an ℓ-leaky forcing set (ℓ ≤ 2, up to 1000 vertices).
- `leaky.py`: Python front end. Accepts a **networkx or Sage** `Graph`. Compiles `lf.c` automatically on first use (needs `gcc`/`cc`).
  - `Z_leaky(G, ell)`
  - `is_leaky_forcing_set(G, ell, S)`: returns `(True, [])`, or `(False, first failing leak set)`
  - `closure(G, S, leaks)`: pure-Python reference implementation, used for cross-checking
  - `GP(n, k)`, `thm414_set(n, k)`: P(n,k) in the manuscript's labelling, and the Section 4.2 set ∪A_i
- `leaky_sage.sage`: Sage examples using `graphs.GeneralizedPetersenGraph` and `graphs.StarGraph`.
- `checks/run_all.py`: reproduces every number in the report. Takes about 1.5 min.

## Usage

```
python3 checks/run_all.py            # needs networkx
sage leaky_sage.sage                 # from this folder
```

Inside a Sage session:

```
sage: import sys; sys.path.insert(0, "/path/to/leaky_forcing")
sage: from leaky import Z_leaky, is_leaky_forcing_set
sage: Z_leaky(graphs.PetersenGraph(), 1)
5
```

## Conventions

- **Colour change rule:** the standard zero forcing rule, except that leaky vertices never force. Any vertex can be a leak, including an initially blue one.
- **Forcing sets:** a set is ℓ-leaky forcing only if it forces the whole graph for every leak set of size at most ℓ.
- **Search:** `Z_leaky` searches only supersets of the vertices of degree ≤ ℓ, using Lemma 2.2 of the manuscript. This speeds up the search and cannot change the answer.
- **Sage labels:** in `graphs.GeneralizedPetersenGraph(n, k)`, outer vertex y_i is `i-1` and inner vertex x_i is `n+i-1`.

The Python interface has been tested with networkx. It has not been run under a real Sage install, but a mock of Sage's `Graph` API returns the right values.
