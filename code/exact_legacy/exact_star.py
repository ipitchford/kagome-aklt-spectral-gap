"""Reconstruct and certify a local spin-2 star obstruction using rational arithmetic.

Usage: python derive_star.py [--output report.json]
Python optimization flags do not disable checks. SymPy is the sole dependency.
This proves an eigenvalue of the five-spin star, not a lattice spectral gap.
"""
import argparse
import json
from itertools import product
from math import comb
from pathlib import Path
import sympy as s


def require(condition, message):
    if not condition:
        raise ValueError(message)


def orbit_keys(total):
    return sorted({(k[0],) + tuple(sorted(k[1:]))
                   for k in product(range(5), repeat=5) if sum(k) == total})


def members(key):
    return [k for k in product(range(5), repeat=5)
            if (k[0],) + tuple(sorted(k[1:])) == key]


def principal_minors(M):
    return [M[:i, :i].det() for i in range(1, M.rows + 1)]


def reconstruct():
    # e_k=|m=2-k>/sqrt(binomial(4,k)); S+ e_k=k e_(k-1).
    a, b = orbit_keys(3), orbit_keys(2)  # magnetizations 7 and 8
    H, U, G = s.zeros(len(a)), s.zeros(len(b), len(a)), s.zeros(len(a))
    for j, key in enumerate(a):
        ms = members(key)
        G[j, j] = s.Rational(len(ms), s.prod(comb(4, k) for k in key))
        outH, outU = {}, {}
        for k in ms:
            for leaf in range(1, 5):
                t = k[0] + k[leaf]
                for p in range(max(0, t - 4), min(4, t) + 1):
                    q = list(k)
                    q[0], q[leaf] = p, t - p
                    q = tuple(q)
                    outH[q] = outH.get(q, 0) + s.Rational(
                        comb(4, p) * comb(4, t - p), comb(8, t))
            for site in range(5):
                if k[site]:
                    q = list(k)
                    q[site] -= 1
                    q = tuple(q)
                    outU[q] = outU.get(q, 0) + k[site]
        for i, key2 in enumerate(a):
            vals = {outH.get(k, 0) for k in members(key2)}
            require(len(vals) == 1, 'star action must preserve leaf symmetry')
            H[i, j] = vals.pop()
        for i, key2 in enumerate(b):
            vals = {outU.get(k, 0) for k in members(key2)}
            require(len(vals) == 1, 'raising action must preserve leaf symmetry')
            U[i, j] = vals.pop()
    r = s.Rational
    T = s.Matrix([[r(1,2),r(3,4),-r(1,4)],
                  [-r(1,2),-r(1,4),0], [1,0,0], [0,-r(3,2),r(3,4)],
                  [0,1,0], [0,0,-r(3,4)], [0,0,1]])
    A = s.Matrix([[r(5,2),r(3,2),0], [r(3,7),r(13,7),-r(6,7)],
                  [r(1,7),-r(3,14),r(3,14)]])
    B = s.Matrix([[r(7,16),r(7,16),-r(1,8)],
                  [r(7,16),r(17,16),-r(3,8)],
                  [-r(1,8),-r(3,8),r(1,2)]])
    return a, b, H, U, G, T, A, B


def validate(H, U, G, T, A, B):
    require(U.rank() == 4 and T.rank() == 3, 'highest-weight dimension mismatch')
    require(U*T == s.zeros(4, 3), 'T is not highest weight')
    require(H*T == T*A, 'invariant-block action mismatch')
    require(B == T.T*G*T, 'Gram reconstruction mismatch')
    require(B == B.T, 'Gram matrix must be symmetric')
    require(all(x > 0 for x in principal_minors(B)), 'Gram matrix must be positive')
    require(B*A == A.T*B, 'block must be Gram-self-adjoint')
    require(all(x > 0 for x in principal_minors(B*A)), 'block must be positive')
    x = s.Symbol('x')
    expected = (2*x - 3)*(14*x*x - 43*x + 2)
    require(s.expand(28*A.charpoly(x).as_expr() - expected) == 0,
            'characteristic polynomial mismatch')
    # The smaller root alpha of 14x^2-43x+2 is positive and below 1/20.
    alpha = (43 - 3*s.sqrt(193))/28
    require(alpha > 0 and alpha < s.Rational(1,20), 'root bound failed')
    eigvec = s.Matrix([3, 2*alpha - 5, 3*alpha - 9])
    require((A*eigvec-alpha*eigvec).applyfunc(s.simplify) == s.zeros(3,1),
            'exact eigenvector failed')
    # An entirely rational upper bound, independent of evaluating radicals.
    v = s.Matrix([3,-5,-9])
    quotient = (v.T*B*A*v)[0]/(v.T*B*v)[0]
    require(quotient == s.Rational(82,1729) and quotient < s.Rational(1,20),
            'rational Rayleigh witness failed')
    return alpha, quotient


def run():
    a,b,H,U,G,T,A,B = reconstruct()
    alpha, quotient = validate(H,U,G,T,A,B)
    rejected = []
    for label in ['block_entry', 'gram_entry']:
        Ac, Bc = A.copy(), B.copy()
        if label == 'block_entry':
            Ac[0,0] += 1
        else:
            Bc[0,0] += 1
        try:
            validate(H,U,G,T,Ac,Bc)
        except ValueError:
            rejected.append(label)
        else:
            raise ValueError('negative control was not rejected: ' + label)
    return {
        'status':'PASS', 'scope':'local star obstruction only; no lattice gap proved',
        'arithmetic':'exact SymPy rational and algebraic arithmetic',
        'sympy_version':s.__version__,
        'basis_orbits':a,'raising_target_orbits':b,
        'H_orbit':str(H),'raising_U':str(U),'orbit_Gram_G':str(G),'T':str(T),
        'Gram_B':str(B),'block_A':str(A),
        'Gram_principal_minors':[str(x) for x in principal_minors(B)],
        'positive_block_principal_minors':[str(x) for x in principal_minors(B*A)],
        'primitive_charpoly':'(2*x - 3)*(14*x**2 - 43*x + 2)',
        'exact_positive_star_eigenvalue':str(alpha),
        'exact_interval':'0 < alpha < 1/20',
        'rayleigh_coefficients':[3,-5,-9],'rayleigh_value':str(quotient),
        'negative_controls_rejected':rejected
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end='')
