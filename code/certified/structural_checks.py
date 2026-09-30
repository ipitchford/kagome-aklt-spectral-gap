"""Exact coverage, normalization, geometry and scalar comparison checks."""
from itertools import product,permutations,combinations
from collections import Counter
from math import comb,prod
from pathlib import Path
import json
import numpy as np
import sympy as sp
from gram_exact import exact_gram
from rational_hub import bases,exbasis,mags
from block_gap import triangle_map
from validated import require
from paths import EVIDENCE

# Corrected Schwinger-polynomial identification.
sts=list(product(range(5),repeat=2));P=sp.Matrix(25,25,lambda i,j:sp.Rational(comb(4,sts[i][0])*comb(4,sts[i][1]),comb(8,sum(sts[i]))) if sum(sts[i])==sum(sts[j]) else 0)
v=sp.zeros(25,1);v[sts.index((0,4))]=1;v[sts.index((1,3))]=-1
bad=v.copy();bad[sts.index((1,3))]=-16
require(P*v==sp.zeros(25,1),'correct polynomial singlet is annihilated')
require(P*bad!=sp.zeros(25,1),'incorrect normalization is rejected')
psi,W,G,omega=triangle_map();require(psi.rank()==27,'injective triangle map')
for a,b in [(0,1),(0,2),(1,2)]:
 H=sp.zeros(125)
 for j,src in enumerate(product(range(5),repeat=3)):
  s=src[a]+src[b]
  for k in range(max(0,s-4),min(4,s)+1):
   dst=list(src);dst[a]=k;dst[b]=s-k;i=int(np.ravel_multi_index(dst,(5,5,5)))
   H[i,j]=sp.Rational(comb(4,k)*comb(4,s-k),comb(8,s))
 require(H*psi==sp.zeros(125,27),'triangle map in every exact edge kernel')

edges=[(0,1),(0,2),(1,2),(3,4),(3,5),(4,5),(2,5)]
GL,DL,ds,bs=exact_gram(6,edges);UL=bases(True);UE=bases(False)
pair={s:np.kron(exbasis(3,s),np.eye(2,dtype=np.int64)) for s in (1,-1)}
UL['y']=np.kron(pair[-1],pair[1]);UE['y']=np.kron(exbasis(3,-1),exbasis(3,1))
for U,n in [(UL,324),(UE,81)]:
 allU=np.concatenate(list(U.values()),axis=1);B=allU.T@allU
 require(allU.shape==(n,n) and np.all(np.diag(B)>0) and np.array_equal(B,np.diag(np.diag(B))),'complete orthogonal symmetry basis')
idx={x:i for i,x in enumerate(bs)}
for perm in [(1,0,2,3,4,5),(0,1,2,4,3,5),(3,4,5,0,1,2)]:
 pi=[idx[tuple(s[p] for p in perm)] for s in bs]
 require(np.array_equal(GL[np.ix_(pi,pi)],GL),'exact D8 Gram commutation')
for t,u in combinations(UL,2):require(np.count_nonzero(UL[t].T@GL@UL[u])==0,'exact orthogonality of all private symmetry sectors')
C=np.zeros((324,3,81),dtype=np.int64)
for ext in product(range(3),repeat=4):
 j=np.ravel_multi_index(ext,(3,)*4)
 for a,b,h,c in [(1,1,0,1),(0,0,2,1),(0,1,1,-1),(1,0,1,-1)]:
  i=np.ravel_multi_index((ext[0],ext[1],a,ext[2],ext[3],b),ds);C[i,h,j]=c
for t in UL:
 u,e=UL[t],UE[t];norm=np.sum(u*u,axis=0)
 for h in range(3):
  raw=C[:,h,:]@e;coords=sp.diag(*[sp.Rational(1,int(q)) for q in norm])*sp.Matrix(u.T@raw)
  require(sp.Matrix(u)*coords==sp.Matrix(raw),'all six exact gluing sectors, including mixed partner')
# Exchange of the two private triangles gives a signed monomial equivalence x <-> y.
pi=[idx[(s[3],s[4],s[5],s[0],s[1],s[2])] for s in bs]
require(np.array_equal(UL['x'][pi].T@GL@UL['x'][pi],UL['x'].T@GL@UL['x']),'mixed partner Gram equivalence')
MH,_,_,mhstates=exact_gram(3,[(0,1),(0,2),(1,2)]);mi={s:i for i,s in enumerate(mhstates)}
for p in permutations(range(3)):
 pi=[mi[tuple(s[k] for k in p)] for s in mhstates]
 require(np.array_equal(MH[np.ix_(pi,pi)],MH),'S3 hub metric invariance')

def counts(ds):
 c=Counter({0:1})
 for d in ds:
  new=Counter()
  for a,n in c.items():
   for j in range(d):new[a+d-1-2*j]+=n
  c=new
 return c

def combine(*cs):
 c=Counter({0:1})
 for d in cs:
  new=Counter()
  for a,n in c.items():
   for b,m in d.items():new[a+b]+=n*m
  c=new
 return c
L=counts(ds);E=counts([3]*4);T=counts([3]*3)
full=combine(L,L,L,T)[0];require(full==106935090,'full localized Mz0 count')
pairdim=combine(E,counts([3,3]),L)[0];paircommon=combine(E,E,counts([3]))[0]
require((pairdim,paircommon)==(34606,3139),'unreduced pair coordinates and common count')
require(counts([4]*4)[0]==44,'star kernel weight-zero dimension')
require(counts([3]*6)[0]==141,'nine-site block kernel weight-zero dimension')

# Exact combinatorics, rectangular tori and arbitrary block-union edge patches.
def up(i,j):return {('A',i,j),('B',i,j),('C',i,j)}
def down(i,j):return {('B',i,j),('A',i+1,j),('C',i+1,j-1)}
def block(i,j):
 hubs={(i,j),(i+1,j),(i+1,j-1)};S=set().union(*(up(*q) for q in hubs))
 edges=set(frozenset(p) for tr in [down(i,j)]+[up(*q) for q in hubs] for p in combinations(tr,2))
 return S,edges,hubs

def geometry(m,n):
 red=lambda s:(s[0],s[1]%m,s[2]%n)
 blocks={};allE={frozenset(red(s) for s in e) for i,j in product(range(m),range(n)) for tr in [up(i,j),down(i,j)] for e in combinations(tr,2)}
 ok=len(allE)==6*m*n
 for i,j in product(range(m),range(n)):
  S,E,h=block(i,j);S={red(s) for s in S};E={frozenset(red(s) for s in e) for e in E};h={(a%m,b%n) for a,b in h}
  ok=ok and len(S)==9 and len(E)==12 and E=={e for e in allE if e<=S} and len(h)==3
  blocks[i,j]=(S,E,h)
 for (S,E,h),(T,F,k) in combinations(blocks.values(),2):
  if S&T:ok=ok and len(h&k)==1 and S&T=={red(s) for s in up(*next(iter(h&k)))}
 for i,j in product(range(m),range(n)):
  p=set().union(*(block(a,b)[0] for a,b in [(i,j),(i-1,j),(i-1,j+1)]))
  ok=ok and len(p)==21 and len({red(s) for s in p})==21
 return {'m':m,'n':n,'ok':bool(ok),'sites':3*m*n,'edges':len(allE)}
rects=[geometry(m,n) for m,n in product(range(3,9),repeat=2)]
require(all(x['ok'] for x in rects),'rectangular torus geometry')
require(not geometry(2,2)['ok'],'small-torus negative control')
base=list(product(range(3),repeat=2));patches=0
for mask in range(1,1<<len(base)):
 chosen=[block(*q) for j,q in enumerate(base) if (mask>>j)&1];hubs=Counter(h for _,_,hs in chosen for h in hs)
 require(all(1<=v<=3 for v in hubs.values()),'partial hub size')
 require(sum(hubs.values())==3*len(chosen),'partial cover incidence identity')
 for S,E,h in chosen:
  require(len(S)==9 and len(E)==12,'patch blocks')
 for (S,E,h),(T,F,k) in combinations(chosen,2):
  require((not S&T) or (len(h&k)==1 and S&T==up(*next(iter(h&k)))),'patch unique-overlap identity')
 patches+=1
# The small rational lower constant is certified using only square comparisons.
b1=sp.Rational(2)-sp.Rational(42,11)*sp.Rational(7,100)
b2=sp.Rational(73)-80*sp.Rational(691,1000)
require(b1>0 and b1*b1>3,'gamma_0 > 7/100')
require(b2>0 and b2*b2>313,'mu_3 > 691/1000')
require(sp.Rational(7,100)*(3*sp.Rational(691,1000)-2)>sp.Rational(51,10000),'uniform gap > 0.0051')
require((sp.Rational(43)-28*sp.Rational(6,125))**2<9*193,'star alpha < 6/125')
require(sp.Rational(33,2246)<sp.Rational(15,1000),'square threshold N=34,a=2')
res={'status':'PASS','normalization_counterexample_rejected':True,'triangle_rank':27,'D8_private_dimensions':{k:v.shape[1] for k,v in UL.items()},'D8_external_dimensions':{k:v.shape[1] for k,v in UE.items()},'full_localized_Mz0_dimension':full,'pair_compressed_Mz0_dimension':pairdim,'pair_common_Mz0_dimension':paircommon,'triple_compressed_Mz0_dimension':combine(counts([3]*4),counts([3]*4),counts([3]),L)[0],'triple_common_Mz0_dimension':counts([3]*12)[0],'rectangular_tori':rects,'patches_checked':patches,'gap_exact':'11*(2-sqrt(3))*(59-3*sqrt(313))/3360','gap_strict_rational_lower':'51/10000','limits':'finite geometry checks corroborate the general counting proof; they do not replace it'}
(EVIDENCE/'structural_checks.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
