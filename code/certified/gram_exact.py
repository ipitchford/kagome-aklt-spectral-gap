"""Independent integer-coefficient Schwinger polynomial Gram construction."""
from fractions import Fraction as F
from collections import defaultdict
from itertools import product
from math import comb,prod
import numpy as np

def exact_gram(n,edges):
 if not (n>0 and all(0<=a<b<n for a,b in edges) and len(set(edges))==len(edges)):
  raise ArithmeticError('simple canonical graph required')
 # A coefficient of the singlet product has l1 norm at most 2**m.
 # Every signed partial Gram sum is bounded by 12**n * 4**m.
 if 12**n*4**len(edges)>=2**62:raise ArithmeticError('int64 exact accumulation bound exceeded')
 deg=[0]*n
 for a,b in edges:deg[a]+=1;deg[b]+=1
 if any(d>4 for d in deg):raise ArithmeticError("degree exceeds four")
 dangling=[4-d for d in deg]; dims=[d+1 for d in dangling]
 poly={(0,)*n:1}
 for a,b in edges:
  new=defaultdict(int)
  for exp,c in poly.items():
   ea=list(exp);ea[b]+=1;new[tuple(ea)]+=c
   eb=list(exp);eb[a]+=1;new[tuple(eb)]-=c
  poly={e:c for e,c in new.items() if c}
 rows=defaultdict(dict);bs=list(product(*(range(d) for d in dims)))
 for j,src in enumerate(bs):
  for exp,c in poly.items():
   dst=tuple(x+y for x,y in zip(src,exp))
   rows[dst][j]=c
 N=len(bs);D=12**n;G=np.zeros((N,N),dtype=np.int64)
 for s,row in rows.items():
  den=prod(comb(4,k) for k in s)
  if D%den:raise ArithmeticError("nonintegral common denominator")
  f=D//den;items=list(row.items())
  for a,ca in items:
   for b,cb in items:G[a,b]+=f*ca*cb
 return G,D,dims,bs

