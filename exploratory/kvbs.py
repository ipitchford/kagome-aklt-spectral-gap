"""Virtual-space (valence-bond) machinery for the spin-2 AKLT model on the kagome lattice.

Conventions
-----------
* Kagome Bravais vectors a1=(2,0), a2=(1,sqrt3); sites A=R, B=R+(1,0), C=R+(1/2,sqrt3/2).
* Up triangle U[i,j] = {A_ij, B_ij, C_ij}; down triangle D[i,j] = {B_ij, A_{i+1,j}, C_{i+1,j-1}}.
* Every site carries four virtual spin-1/2 ("half-edges"), two in its up triangle and two in its
  down triangle.  A spin-2 site is the symmetric subspace of its four qubits.
* Spin basis index k = number of down spins, so |k> has S^z = s-k.
* A region X (set of sites, with all kagome edges joining two sites of X) has the valence-bond map
  Psi_X : V_X -> (C^5)^{X}, where V_X is the tensor product over boundary sites a of the spin-(z'_a/2)
  space of the z'_a dangling half-edges of a (fused, i.e. symmetrised).  Psi_X is built by placing a
  singlet on every internal edge, attaching the fused dangling legs, and projecting each site on
  spin 2.  The range of Psi_X is the kernel of the AKLT Hamiltonian of X (standard; checked
  numerically for small regions in tests).

All arithmetic here is floating point (numpy, float64).  Exact rational checks live elsewhere.
"""
from __future__ import annotations

import itertools
from functools import lru_cache
from math import comb, sqrt

import numpy as np
import opt_einsum as oe

# ----------------------------------------------------------------------------- geometry

def up(i, j):
    return [("A", i, j), ("B", i, j), ("C", i, j)]


def down(i, j):
    return [("B", i, j), ("A", i + 1, j), ("C", i + 1, j - 1)]


def site_triangles(s):
    """(up index, down index) of a site."""
    t, i, j = s
    if t == "A":
        return (i, j), (i - 1, j)
    if t == "B":
        return (i, j), (i, j)
    if t == "C":
        return (i, j), (i - 1, j + 1)
    raise ValueError(s)


def site_edges(s):
    """The four kagome edges at s, each as a frozenset, in a fixed order:
    two up-triangle edges then two down-triangle edges."""
    (ui, uj), (di, dj) = site_triangles(s)
    U = up(ui, uj)
    D = down(di, dj)
    out = []
    for tri in (U, D):
        for o in tri:
            if o != s:
                out.append(frozenset((s, o)))
    return out


def edges_of_triangle(tri):
    return [frozenset(p) for p in itertools.combinations(tri, 2)]


def induced_edges(sites):
    S = set(sites)
    E = set()
    for s in S:
        for e in site_edges(s):
            if e <= S:
                E.add(e)
    return E


def position(s):
    t, i, j = s
    x = 2 * i + j
    y = sqrt(3) * j
    if t == "B":
        x += 1
    elif t == "C":
        x += 0.5
        y += sqrt(3) / 2
    return x, y

# ----------------------------------------------------------------------------- local tensors

@lru_cache(None)
def sym_embed(n):
    """Isometry (2,)*n x (n+1): column k is the normalised symmetric n-qubit state with k ones."""
    T = np.zeros((2 ** n, n + 1))
    for bits in itertools.product((0, 1), repeat=n):
        k = sum(bits)
        idx = int("".join(map(str, bits)), 2) if n else 0
        T[idx, k] = 1.0 / sqrt(comb(n, k))
    return T.reshape((2,) * n + (n + 1,))


SINGLET = np.array([[0.0, 1.0], [-1.0, 0.0]]) / sqrt(2.0)


def spin_ops(dim):
    """S^z, S^+ in the basis k = number of down spins (S^z = s-k)."""
    s = (dim - 1) / 2.0
    Sz = np.diag([s - k for k in range(dim)])
    Sp = np.zeros((dim, dim))
    for k in range(1, dim):
        m = s - k
        Sp[k - 1, k] = sqrt(s * (s + 1) - m * (m + 1))
    return Sz, Sp

# ----------------------------------------------------------------------------- regions

class Region:
    """A set of kagome sites together with all kagome edges internal to it."""

    def __init__(self, sites, name=""):
        self.sites = sorted(set(sites))
        self.name = name
        self.S = set(self.sites)
        self.edges = sorted(induced_edges(self.sites), key=lambda e: sorted(e))
        # dangling half-edges per site: list of (site, edge)
        self.dangling = {}
        for s in self.sites:
            d = [e for e in site_edges(s) if not e <= self.S]
            if d:
                self.dangling[s] = d
        self.leg_sites = [s for s in self.sites if s in self.dangling]
        self.leg_dims = [len(self.dangling[s]) + 1 for s in self.leg_sites]
        self.vdim = int(np.prod(self.leg_dims)) if self.leg_dims else 1

    def __repr__(self):
        return f"Region({self.name}, n={len(self.sites)}, edges={len(self.edges)}, legs={self.leg_dims}, vdim={self.vdim})"

    # ------------------------------------------------------------------ Psi
    def psi(self):
        """Dense Psi_X with shape (5,)*n + tuple(leg_dims)."""
        n = len(self.sites)
        idx = {}
        counter = itertools.count()

        def q(s, e):
            key = (s, e)
            if key not in idx:
                idx[key] = next(counter)
            return idx[key]

        P = sym_embed(4).reshape(16, 5).T.reshape(5, 2, 2, 2, 2)  # phys <- 4 qubits
        operands = []
        phys_ind = []
        for a, s in enumerate(self.sites):
            p = 10_000 + a
            phys_ind.append(p)
            operands += [P, [p] + [q(s, e) for e in site_edges(s)]]
        for e in self.edges:
            s, t = sorted(e)
            operands += [SINGLET, [q(s, e), q(t, e)]]
        leg_ind = []
        for b, s in enumerate(self.leg_sites):
            d = self.dangling[s]
            v = 20_000 + b
            leg_ind.append(v)
            operands += [sym_embed(len(d)), [q(s, e) for e in d] + [v]]
        out = oe.contract(*operands, phys_ind + leg_ind, optimize="greedy")
        return out

    def gram(self):
        Psi = self.psi()
        n = len(self.sites)
        M = Psi.reshape(5 ** n, self.vdim)
        return M.T @ M

    def leg_sz(self):
        """Total S^z of each virtual basis vector (flattened)."""
        sz = np.zeros(1)
        for d in self.leg_dims:
            s = (d - 1) / 2.0
            sz = (sz[:, None] + np.array([s - k for k in range(d)])[None, :]).ravel()
        return sz


def contraction_map(big_parts, ext_region_sites=None, contract_edges=None):
    """Map C : V_ext -> (tensor over parts of V_part), for the union of `big_parts` (list of Regions,
    disjoint site sets).  Edges joining two different parts are contracted with singlets; the union's
    dangling legs are the external legs.  Returns (C as array with shape part_legs... + ext_legs,
    ext Region)."""
    union = Region([s for R in big_parts for s in R.sites], name="union")
    if contract_edges is not None:
        # legs of the union must include half-edges of non-contracted cross edges
        union.edges = [e for e in union.edges if any(e <= R.S for R in big_parts) or e in contract_edges]
        E = set(union.edges)
        union.dangling = {}
        for s in union.sites:
            d = [e for e in site_edges(s) if e not in E]
            if d:
                union.dangling[s] = d
        union.leg_sites = [s for s in union.sites if s in union.dangling]
        union.leg_dims = [len(union.dangling[s]) + 1 for s in union.leg_sites]
        union.vdim = int(np.prod(union.leg_dims))
    counter = itertools.count()
    idx = {}

    def q(s, e):
        key = (s, e)
        if key not in idx:
            idx[key] = next(counter)
        return idx[key]

    operands = []
    part_ind = []
    for pi, R in enumerate(big_parts):
        for b, s in enumerate(R.leg_sites):
            d = R.dangling[s]
            v = 10_000 + 100 * pi + b
            part_ind.append(v)
            # bra: Emb^dagger on all dangling qubits of s in R
            operands += [sym_embed(len(d)), [q(s, e) for e in d] + [v]]
    # singlets on edges joining two parts
    for e in union.edges:
        s, t = sorted(e)
        ps = [i for i, R in enumerate(big_parts) if s in R.S][0]
        pt = [i for i, R in enumerate(big_parts) if t in R.S][0]
        if ps != pt and (contract_edges is None or e in contract_edges):
            operands += [SINGLET, [q(s, e), q(t, e)]]
    ext_ind = []
    for b, s in enumerate(union.leg_sites):
        d = union.dangling[s]
        v = 20_000 + b
        ext_ind.append(v)
        operands += [sym_embed(len(d)), [q(s, e) for e in d] + [v]]
    C = oe.contract(*operands, part_ind + ext_ind, optimize="greedy")
    return C, union


def sqrtm_psd(G):
    w, V = np.linalg.eigh((G + G.T) / 2)
    if w.min() <= 0:
        raise ValueError(f"Gram not positive definite: min eig {w.min()}")
    return (V * np.sqrt(w)) @ V.T, (V / np.sqrt(w)) @ V.T, w


def orth_basis(D):
    """Orthonormal basis of range(D) via symmetric orthonormalisation D (D^T D)^{-1/2}."""
    G = D.T @ D
    _, Gmh, w = sqrtm_psd(G)
    return D @ Gmh, w
