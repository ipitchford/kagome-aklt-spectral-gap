#!/usr/bin/env python3
"""Exact local identities and finite combinatorial checks; NOT an AKLT gap proof."""
from fractions import Fraction as F
from math import comb
from itertools import combinations
import json
import sympy as s

def require(value, message):
    if not value:
        raise ValueError(message)

def projector():
    # e_i=sqrt(binomial(4,i)) |m=2-i>; metric D_i=binomial(4,i).
    D=s.diag(*[comb(4,i)*comb(4,j) for i in range(5) for j in range(5)])
    P=s.zeros(25)
    for i in range(5):
        for j in range(5):
            for k in range(5):
                l=i+j-k
                if 0 <= l < 5:
                    P[5*i+j,5*k+l]=s.Rational(comb(4,k)*comb(4,l),comb(8,i+j))
    E=s.zeros(5); G=s.zeros(5)
    for i in range(1,5): E[i-1,i]=5-i
    for i in range(4): G[i+1,i]=i+1
    Z=s.diag(2,1,0,-1,-2)
    X=s.kronecker_product(Z,Z)+(s.kronecker_product(E,G)+s.kronecker_product(G,E))/2
    I=s.eye(25)
    poly=X*(X+3*I)*(X+5*I)*(X+6*I)/2520
    require(P==poly, 'polynomial projector does not match Clebsch-Gordan construction')
    require(P*P==P, 'projector idempotency')
    require(D*P==P.T*D,'weighted self-adjointness')
    require(s.trace(P)==9,'spin-four rank')
    return {'dimension':25,'rank':9,'polynomial_denominator':2520,
            'idempotent':True,'metric_self_adjoint':True,'independent_constructions_agree':True}

def patch(N):
    edges={}
    for x in range(N):
        for y in range(N):
            if x+1<N: edges[((x,y),(x+1,y))]=int(y not in (0,N-1))
            if y+1<N: edges[((x,y),(x,y+1))]=int(x not in (0,N-1))
    for t in range(N):
        for e in [((-1,t),(0,t)),((N-1,t),(N,t)),((t,-1),(t,0)),((t,N-1),(t,N))]:
            edges[e]=0
    return edges

def signature(pair):
    x0=min(p[0] for e in pair for p in e)
    y0=min(p[1] for e in pair for p in e)
    return tuple(sorted(tuple((x-x0,y-y0) for x,y in e) for e in pair))

def patch_check(N,a=F(7,5),corrupt=False):
    edges=patch(N)
    vertices={v for e in edges for v in e}
    f=2*(2*N-1)+(N-1)*(N-2)*a
    g=2*N+2*(N-2)*a+(N-2)**2*a*a
    totals=[F(0),F(0)]
    for e,power in edges.items():
        totals[e[0][0]==e[1][0]]+=a**power
    require(totals==[f,f], 'one-edge translation coefficients')
    require(len(vertices)==N*N+4*N,'patch vertex count')
    require(len(edges)==2*N*(N+1),'patch edge count')
    accumulated={}
    for (e,p),(h,q) in combinations(edges.items(),2):
        key=signature((e,h))
        accumulated[key]=accumulated.get(key,F(0))+a**(p+q)
    if corrupt:
        key=next(k for k in accumulated if set(k[0])&set(k[1]))
        accumulated[key]+=F(1,1000)
    adjacent=[]; disjoint=[]
    for key,weight in accumulated.items():
        if set(key[0])&set(key[1]): adjacent.append(weight)
        else: disjoint.append(weight)
    require(len(adjacent)==6, 'six unoriented adjacent edge-pair translation classes')
    require(all(v==g for v in adjacent), 'adjacent anticommutator coefficient must equal g')
    require(all(v<=g for v in disjoint), 'disjoint product coefficients cannot exceed g')
    threshold=(2*(2*N-1)+(N-1)*(N-2)*a*a-g)/f
    return {'N':N,'a':str(a),'sites':len(vertices),'edges':len(edges),
            'f':str(f),'g':str(g),'threshold':str(threshold),
            'adjacent_classes':len(adjacent),'disjoint_classes':len(disjoint),
            'maximum_disjoint_weight':str(max(disjoint))}

def negative_controls():
    caught=[]
    try: patch_check(4,corrupt=True)
    except ValueError: caught.append('corrupt_adjacent_coefficient')
    require(caught==['corrupt_adjacent_coefficient'],'negative control did not fail')
    # A wrong normalization changes a rank-nine projector into a nonprojector.
    P=s.diag(*([s.Rational(9,10)]*9+[0]*16))
    try: require(P*P==P,'corrupted projector')
    except ValueError: caught.append('wrong_projector_normalization')
    require(len(caught)==2,'normalization negative control did not fail')
    return caught

if __name__=='__main__':
    result={'claim_ceiling':'exact local identities and finite patch counting only; no uniform gap proved',
            'projector':projector(), 'patches':[patch_check(n) for n in range(4,11)],
            'negative_controls':negative_controls()}
    print(json.dumps(result,indent=2))
