"""Exact rational evidence for the five-site cross LTI obstruction.

Standard-library-only checker.  The representation-theoretic kernel proof is
in ../local_sdp_review.md; this script verifies its finite support geometry and
constructs the negative vector directly in the 5-spin tensor product.
"""

from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb
import json


def require(condition, message):
    if not condition:
        raise ValueError(message)


C = ((0, 0), (1, 0), (0, 1), (-1, 0), (0, -1))
CS = frozenset(C)


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def translate(support, shift):
    return frozenset(add(x, shift) for x in support)


def canonical(support):
    # Lexicographic minimum is covariant under translations.
    x = min(support)
    return tuple(sorted((y[0] - x[0], y[1] - x[1]) for y in support))


differences = Counter((x[0] - y[0], x[1] - y[1]) for x in C for y in C)
require(differences[(0, 0)] == 5, "cross cardinality")
require(max(n for shift, n in differences.items() if shift != (0, 0)) == 2,
        "nonzero cross intersection bound")
for shift, count in differences.items():
    require(len(CS & translate(CS, shift)) == count, "intersection multiplicity")

support_classes = {}
for size in range(1, 6):
    classes = defaultdict(list)
    for support in combinations(C, size):
        classes[canonical(support)].append(frozenset(support))
    support_classes[size] = classes
    if size >= 3:
        require(all(len(members) == 1 for members in classes.values()),
                "unique higher-support classes")

pair_class_sizes = sorted(len(members) for members in support_classes[2].values())
require(pair_class_sizes == [1, 1, 2, 2, 2, 2], "pair class counts")
for members in support_classes[2].values():
    if len(members) == 2:
        support, other = members
        require(other == frozenset((-x[0], -x[1]) for x in support),
                "repeated pair class is central inversion")
        a, b = tuple(support)
        shift = add(a, b)
        # Translating the inverted pair back reverses its two endpoints.
        require(add((-a[0], -a[1]), shift) == b, "first endpoint reversal")
        require(add((-b[0], -b[1]), shift) == a, "second endpoint reversal")

# e_k = |m=2-k>/sqrt(binomial(4,k)); its norm squared is 1/binomial(4,k).
# The first site is the center.  Remaining sites are four leaves.
orbit_coefficients = {
    (0, 0, 0, 0, 3): 0,
    (0, 0, 0, 1, 2): 1,
    (0, 0, 1, 1, 1): -8,
    (1, 0, 0, 0, 2): -3,
    (1, 0, 0, 1, 1): 12,
    (2, 0, 0, 0, 1): -15,
    (3, 0, 0, 0, 0): 20,
}
basis = [k for k in product(range(5), repeat=5) if sum(k) == 3]
psi = {k: Q(orbit_coefficients[(k[0],) + tuple(sorted(k[1:]))]) for k in basis}


def apply_star(vector):
    out = defaultdict(Q)
    for k, value in vector.items():
        for leaf in range(1, 5):
            total = k[0] + k[leaf]
            for p in range(max(0, total - 4), min(4, total) + 1):
                target = list(k)
                target[0], target[leaf] = p, total - p
                coefficient = Q(comb(4, p) * comb(4, total - p), comb(8, total))
                out[tuple(target)] += coefficient * value
    return dict(out)


def inner(left, right):
    out = Q(0)
    for k, value in left.items():
        denominator = 1
        for number in k:
            denominator *= comb(4, number)
        out += value * right.get(k, Q(0)) / denominator
    return out


raised = defaultdict(Q)
for k, value in psi.items():
    for site in range(5):
        if k[site]:
            target = list(k)
            target[site] -= 1
            raised[tuple(target)] += k[site] * value
require(all(value == 0 for value in raised.values()), "highest-weight condition")

kpsi = apply_star(psi)
k2psi = apply_star(kpsi)
norm = inner(psi, psi)
first = inner(psi, kpsi)
second = inner(kpsi, kpsi)
require(second == inner(psi, k2psi), "K-squared consistency")
require(norm == 157, "witness norm")
require(first == Q(115, 14), "witness first moment")
require(second == Q(311, 196), "witness second moment")
g0 = second - first / 2
require(g0 == -Q(247, 98), "negative constant")
require(-first / 2 == -Q(115, 28), "negative delta coefficient")
require(g0 < 0 and first > 0, "strict negativity for nonnegative delta")

# Controls check that corrupted data fail ordinary validation even under -O.
rejected_controls = []
try:
    require(pair_class_sizes == [1, 1, 1, 2, 2, 2], "corrupted pair class")
except ValueError:
    rejected_controls.append("corrupted_pair_class")
corrupted_psi = dict(psi)
corrupted_psi[(3, 0, 0, 0, 0)] += 1
try:
    require(inner(corrupted_psi, corrupted_psi) == 157, "corrupted witness")
except ValueError:
    rejected_controls.append("corrupted_witness")
require(len(rejected_controls) == 2, "all negative controls rejected")

result = {
    "status": "PASS",
    "arithmetic": "Python fractions.Fraction; no numerical tolerances",
    "maximum_nonzero_cross_intersection": 2,
    "pair_translation_class_sizes": pair_class_sizes,
    "pair_translation_classes": [
        {"canonical": [list(site) for site in key],
         "representatives": [[list(site) for site in sorted(support)]
                             for support in members]}
        for key, members in sorted(support_classes[2].items())
    ],
    "support_counts_by_size": {
        str(size): {"subsets": sum(len(members) for members in classes.values()),
                    "translation_classes": len(classes)}
        for size, classes in support_classes.items()
    },
    "higher_support_translation_classes_are_singletons": True,
    "witness_orbit_coefficients": [
        {"orbit": list(orbit), "coefficient": coefficient}
        for orbit, coefficient in orbit_coefficients.items()
    ],
    "witness_is_highest_weight_m7": True,
    "witness_norm_squared": str(norm),
    "witness_K_expectation": str(first),
    "witness_K_squared_expectation": str(second),
    "witness_g_delta_expectation": "-247/98 - (115/28)*delta",
    "normalized_g0_expectation": str(g0 / norm),
    "rational_parameter_upper_bound": str(2 * second / first - 1),
    "negative_controls_rejected": rejected_controls,
    "scope": "cross-supported generator after discarding all disjoint-edge anticommutators; not a lattice-gap obstruction",
}
print(json.dumps(result, indent=2))
