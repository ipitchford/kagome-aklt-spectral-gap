"""Exact replay of the sharp two-bond spin-2 anticommutator constant.

Prints a JSON receipt to stdout; never writes package files. Requires SymPy.
All checks use explicit exceptions and remain active under python -O.
"""
import json
import sympy as s
from sympy.physics.wigner import wigner_6j


def require(value, message):
    if not value:
        raise ValueError(message)


def validate(entries):
    r = s.Rational
    expected = {2:r(1,70),3:-r(1,14),4:r(3,14),5:-r(1,2),6:s.Integer(1)}
    require(entries == expected, 'exact recoupling table mismatch')
    eta = r(1,2)
    checks = {}
    for J in range(2,6):
        c = abs(entries[J])
        require(0 < c < 1, 'invalid principal cosine')
        # On the two-projector block, S=P+Q has eigenvalues 1 +/- c.
        # {P,Q}+eta(P+Q) = S^2+(eta-1)S.
        vals = [s.factor((1-c)*(eta-c)),s.factor((1+c)*(eta+c))]
        require(all(x >= 0 for x in vals), 'anticommutator bound fails')
        checks[J] = [str(x) for x in vals]
    require(entries[6] == 1, 'common maximal-spin range mismatch')
    # J=6 has P=Q=1: anticommutator+eta(P+Q)=2+2eta.
    require(2+2*eta > 0, 'common range bound failed')
    # The J=5 sum eigenvalue 1-c=1/2 gives ratio anti/sum = -c.
    c = abs(entries[5])
    sum_eigenvalue = 1-c
    anticommutator_eigenvalue = sum_eigenvalue**2-sum_eigenvalue
    require(sum_eigenvalue > 0 and
            -anticommutator_eigenvalue/sum_eigenvalue == eta,
            'sharpness witness failed')
    return checks


def main():
    entries = {J:s.simplify(9*wigner_6j(2,2,4,2,J,4)) for J in range(2,7)}
    checks = validate(entries)
    corrupted = dict(entries)
    corrupted[5] += s.Rational(1,100)
    try:
        validate(corrupted)
    except ValueError:
        negative_control = 'rejected'
    else:
        raise ValueError('corrupted recoupling entry was accepted')
    print(json.dumps({
        'status':'PASS',
        'scope':'sharp eta=1/2 for two adjacent spin-2 maximal-spin projectors',
        'sympy_version':s.__version__,
        'signed_9_times_6j':{str(J):str(v) for J,v in entries.items()},
        'optimal_eta':'1/2',
        'bound_operator_eigenvalues_on_nontrivial_blocks':checks,
        'J5_sum_eigenvalue':'1/2','J5_anticommutator_eigenvalue':'-1/4',
        'corrupted_entry_negative_control':negative_control,
        'lattice_implication':'edgewise scalar estimate on degree-four graph yields only H^2 >= -2 H'
    },indent=2))


if __name__ == '__main__':
    main()
