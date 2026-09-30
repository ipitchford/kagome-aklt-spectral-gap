"""Exact triangle spectrum and strict triangle-cover overlap obstruction.
Requires Python 3 and SymPy. Prints deterministic JSON; writes no files.
At one site use b_k=sqrt(binomial(4,k)) |m=2-k>, k=0,...,4.
Thus matrix(h)=D^-1 h D, physical coefficients=D*x, and Gram W=D^2.
This verifies a failed sufficient criterion, not a thermodynamic gap.
"""
import itertools
import json
import sympy as s


def require(condition, message):
    if not condition:
        raise ValueError(message)


c = [s.binomial(4, i) for i in range(5)]


def weight(a):
    return s.prod(c[i] for i in a)


def projector_entry(a, b, i, j):
    if any(a[t] != b[t] for t in range(len(a)) if t not in (i, j)):
        return s.S.Zero
    if a[i] + a[j] != b[i] + b[j]:
        return s.S.Zero
    return c[b[i]] * c[b[j]] / s.binomial(8, b[i] + b[j])


def hamiltonian(states, edges):
    return s.Matrix(len(states), len(states), lambda i, j: sum(
        projector_entry(states[i], states[j], a, b) for a, b in edges))


def validate_projector(P, W):
    require(P * P == P, 'projector is not idempotent')
    require(W * P == P.T * W, 'projector is not Gram-self-adjoint')
    require(s.trace(P) == 9, 'two-site projector rank is not nine')


pair_states = list(itertools.product(range(5), repeat=2))
P = hamiltonian(pair_states, [(0, 1)])
Wp = s.diag(*[weight(a) for a in pair_states])
validate_projector(P, Wp)
corrupted = P.copy()
corrupted[0, 0] += s.Rational(1, 10)
try:
    validate_projector(corrupted, Wp)
except ValueError:
    corrupted_projector_rejected = True
else:
    raise ValueError('negative control accepted a corrupted projector')

blocks = {}
fullspec = {}
for k in range(13):
    states = [a for a in itertools.product(range(5), repeat=3) if sum(a) == k]
    H = hamiltonian(states, [(0, 1), (1, 2), (0, 2)])
    W = s.diag(*[weight(a) for a in states])
    require(W * H == H.T * W, 'triangle Hamiltonian not Gram-self-adjoint')
    ev = H.eigenvals()
    require(sum(ev.values()) == len(states), 'incomplete exact triangle spectrum')
    for e, m in ev.items():
        fullspec[str(e)] = fullspec.get(str(e), 0) + int(m)
    nn = H.nullspace()
    if nn:
        N = s.Matrix.hstack(*nn)
        Q = N * (N.T * W * N).inv() * N.T * W
    else:
        Q = s.zeros(len(states))
    require(Q * Q == Q and H * Q == s.zeros(len(states)) and W * Q == Q.T * W,
            'invalid triangle ground projector')
    blocks[k] = (states, {a: i for i, a in enumerate(states)}, Q)
expected_spectrum = {'0': 27, '11/14': 18, '6/7': 7, '69/70': 10,
                     '36/35': 5, '15/14': 14, '10/7': 9, '3/2': 22, '3': 13}
require(fullspec == expected_spectrum, 'triangle spectrum differs from claimed spectrum')

# The fifteen-dimensional degree-four triangle block has a three-dimensional
# kernel. Verify the explicit cyclic basis used in the hand derivation.
st4, idx4, Q4 = blocks[4]
v0 = {(0, 1, 3): -18, (0, 2, 2): 8, (1, 0, 3): 18,
      (1, 2, 1): -3, (2, 0, 2): -8, (2, 1, 1): 3}
N4 = s.Matrix(len(st4), 3, lambda i, j: v0.get(st4[i][j:] + st4[i][:j], 0))
W4 = s.diag(*[weight(a) for a in st4])
G4 = N4.T * W4 * N4
require(G4 == 19872 * s.eye(3) - 3168 * s.ones(3), 'tiny Gram matrix mismatch')
require(Q4 == N4 * G4.inv() * N4.T * W4, 'cyclic basis does not span ground block')

states = [a for a in itertools.product(range(5), repeat=5) if sum(a) == 5]
n = len(states)
W = s.diag(*[weight(a) for a in states])


def ground_entry(a, b, sites):
    if any(a[t] != b[t] for t in range(5) if t not in sites):
        return s.S.Zero
    aa = tuple(a[t] for t in sites)
    bb = tuple(b[t] for t in sites)
    if sum(aa) != sum(bb):
        return s.S.Zero
    _, indices, Q = blocks[sum(aa)]
    return Q[indices[aa], indices[bb]]


E = s.Matrix(n, n, lambda i, j: ground_entry(states[i], states[j], (0, 1, 2)))
F = s.Matrix(n, n, lambda i, j: ground_entry(states[i], states[j], (0, 3, 4)))
require(E * E == E and F * F == F and W * E == E.T * W and W * F == F.T * W,
        'invalid embedded ground projectors')


def perm_sign(t):
    return (-1) ** sum(t[i] > t[j] for i in range(3) for j in range(i + 1, 3))


x = s.Matrix([perm_sign(a[:3]) if sorted(a[:3]) == [0, 1, 2] and a[3:] == (0, 2)
              else 0 for a in states])
require(E * x == x, 'witness is not in first triangle ground space')
intersection = (2 * s.eye(n) - E - F).nullspace()
require(not intersection, 'witness sector contains common ground states')
num = (x.T * W * F * x)[0]
den = (x.T * W * x)[0]
ratio = num / den


def validate_obstruction(proposed_ratio):
    require(proposed_ratio == ratio, 'claimed ratio differs from exact contraction')
    require(proposed_ratio == s.Rational(55, 414), 'unexpected contraction value')
    require(proposed_ratio - s.Rational(1, 9) == s.Rational(1, 46),
            'strict overlap obstruction not established')


validate_obstruction(ratio)
try:
    validate_obstruction(s.Rational(1, 10))
except ValueError:
    corrupted_ratio_rejected = True
else:
    raise ValueError('negative control accepted corrupted overlap ratio')
diagonals = []
for cc in range(3):
    _, idx, Q = blocks[cc + 2]
    z = (cc, 0, 2)
    diagonals.append(Q[idx[z], idx[z]])
require(diagonals == [0, s.Rational(1, 6), s.Rational(16, 69)], 'diagonal values mismatch')

out = {
    'status': 'exact_rational_verification_passed',
    'scope': 'triangle-cover obstruction only; no thermodynamic gap proved',
    'basis': 'b_k=sqrt(binomial(4,k))|m=2-k>; operator D^-1 h D; physical vector D*x; Gram W=D^2',
    'triangle_spectrum': expected_spectrum,
    'triangle_gap': '11/14',
    'sector_total_down': 5, 'sector_total_magnetization': 5, 'sector_dimension': n,
    'ranks': [int(s.trace(E)), int(s.trace(F))], 'common_range_dimension': len(intersection),
    'F_diagonal_entries': [str(v) for v in diagonals],
    'witness_nonzero': [{'state': a, 'coefficient': str(v)} for a, v in zip(states, x) if v],
    'numerator': str(num), 'denominator': str(den), 'rayleigh_ratio': str(ratio),
    'eta_squared_lower_bound': str(ratio), 'threshold_squared': '1/9',
    'strict_obstruction_difference': str(ratio - s.Rational(1, 9)),
    'tiny_Gram_matrix': [[str(z) for z in row] for row in G4.tolist()],
    'negative_controls': {'corrupted_projector_rejected': corrupted_projector_rejected,
                          'corrupted_ratio_rejected': corrupted_ratio_rejected}}
print(json.dumps(out, indent=2, sort_keys=True))
