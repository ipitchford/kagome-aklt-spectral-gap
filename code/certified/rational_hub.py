"""Exact rational boundary frames; no square roots or numerical rank choices."""
from fractions import Fraction as F
from math import comb,prod
from itertools import product
import numpy as np, sympy as sp, time,pickle,os
from gram_exact import exact_gram

def require(test,message):
 if not test:raise ArithmeticError(message)

def exbasis(d,sgn):
 cols=[]
 for i in range(d):
  for j in range(i,d):
   if i==j and sgn<0:continue
   v=np.zeros((d,d),dtype=np.int64);v[i,j]=1
   if i!=j:v[j,i]=sgn
   cols.append(v.ravel())
 return np.array(cols,dtype=np.int64).T

def bases(half):
 pair={s:exbasis(3,s) for s in (1,-1)}
 if half:pair={s:np.kron(U,np.eye(2,dtype=np.int64)) for s,U in pair.items()}
 out={}
 for s,tag in [(1,'p'),(-1,'m')]:
  for t,sg in [(1,'p'),(-1,'m')]:out[tag+sg]=np.kron(pair[s],pair[s])@exbasis(pair[s].shape[1],t)
 out['x']=np.kron(pair[1],pair[-1])
 return out

def mags(ds):
 return np.array([sum(d-1-2*k for d,k in zip(ds,s))//2 for s in product(*(range(d) for d in ds))],dtype=int)

def rational_frame(G,m):
 """Exact Sz-preserving congruence to a diagonal metric in [1,4)."""
 n=G.rows;W=sp.zeros(n);Winv=sp.zeros(n);diag=[None]*n
 for z in sorted(set(m)):
  ix=list(map(int,np.where(m==z)[0]));B=G.extract(ix,ix)
  L,D=B.LDLdecomposition(hermitian=False)
  require(L*D*L.T==B,'exact LDL reconstruction')
  scales=[]
  for i in range(len(ix)):
   d=D[i,i];require(d>0,'positive exact pivot')
   p=sp.Integer(1)
   while d*p*p<1:p*=2
   while d*p*p>=4:p/=2
   scales.append(p);diag[ix[i]]=d*p*p
  Q=L.T.inv()*sp.diag(*scales)
  Qinv=sp.diag(*[1/p for p in scales])*L.T
  require(Q.T*B*Q==sp.diag(*[diag[i] for i in ix]),'diagonal frame')
  for a,i in enumerate(ix):
   for b,j in enumerate(ix):W[i,j]=Q[a,b];Winv[i,j]=Qinv[a,b]
 require(W*Winv==sp.eye(n),'frame inverse')
 return W,Winv,diag

def construct():
 t=time.time()
 edges=[(0,1),(0,2),(1,2),(3,4),(3,5),(4,5),(2,5)]
 GL,DL,ds,bs=exact_gram(6,edges);MH,DH,_,_=exact_gram(3,[(0,1),(0,2),(1,2)])
 # Positive scalar normalisation is irrelevant to orthogonal projectors.
 M=sp.Matrix(MH)*sp.Rational(27,int(np.trace(MH)))
 ULS=bases(True);UES=bases(False);ml=mags(ds);me=mags([3]*4)
 C=np.zeros((324,3,81),dtype=np.int64)
 for ext in product(range(3),repeat=4):
  e=np.ravel_multi_index(ext,(3,)*4)
  for a,b,h,coef in [(1,1,0,1),(0,0,2,1),(0,1,1,-1),(1,0,1,-1)]:
   l=np.ravel_multi_index((ext[0],ext[1],a,ext[2],ext[3],b),tuple(ds));C[l,h,e]=coef
 out={}
 for tag in ULS:
  U,V=ULS[tag],UES[tag];d,e=U.shape[1],V.shape[1]
  norms=np.sum(U*U,axis=0)
  mL=np.array([ml[np.flatnonzero(U[:,i])[0]] for i in range(d)])
  mE=np.array([me[np.flatnonzero(V[:,i])[0]] for i in range(e)])
  require(np.all((ml[:,None]-mL[None,:])*U==0),'private Sz basis')
  require(np.all((me[:,None]-mE[None,:])*V==0),'external Sz basis')
  require(12**6*4**7*int(np.max(np.sum(np.abs(U),axis=0)))**2 < 2**62,"projected int64 Gram bound")
  gs=sp.Matrix(U.T@GL@U)/DL
  # Ground tensor must remain wholly in the matching D8 symmetry sector.
  Cs=[]
  for h in range(3):
   raw=C[:,h,:]@V
   ch=sp.diag(*[sp.Rational(1,int(v)) for v in norms])*sp.Matrix(U.T@raw)
   require(sp.Matrix(U)*ch==sp.Matrix(raw),'gluing equivariance')
   Cs.append(ch)
  WL,WLinv,diag=rational_frame(gs,mL)
  Cs=[WLinv*c for c in Cs];G=sp.diag(*diag)
  E=Cs[0].T*G*Cs[0]+Cs[1].T*G*Cs[1]/2+Cs[2].T*G*Cs[2]
  WE,WEinv,ediag=rational_frame(E,mE)
  Cs=[c*WE for c in Cs]
  D=np.empty((e,e,3,3),dtype=object)
  for h,hp in product(range(3),repeat=2):
   chan=Cs[h].T*G*Cs[hp]
   D[:,:,h,hp]=np.array(chan,dtype=object)
  for h,c in enumerate(Cs):
   for i,j in product(range(d),range(e)):
    if mL[i]+1-h!=mE[j]:require(c[i,j]==0,'exact Sz selection')
  out[tag]={'d':d,'e':e,'mL':mL,'mE':mE,'diag':diag,'C':Cs,'D':D,'WL':WL,'WE':WE}
  print('rational',tag,d,e,'conditioned metrics',min(diag),max(diag),'t',time.time()-t,flush=True)
 return out,M

def make_local(data,M,tag):
 """Exact V=W C and exact right inverse factor K=V(C^T W C)^-1."""
 t=time.time();s=data[tag];d,e=s['d'],s['e'];n=d*27;r=e*9
 C=sp.MutableSparseMatrix(n,r,{})
 for i,h,k in product(range(d),range(3),range(e)):
  v=s['C'][h][i,k]
  if not v:continue
  for x,y in product(range(3),repeat=2):C[i*27+x*9+y*3+h,k*9+x*3+y]=v
 # The d-dimensional private metric is diagonal, so multiplication is cheap.
 W=sp.kronecker_product(sp.SparseMatrix(sp.diag(*s['diag'])),sp.SparseMatrix(M))
 V=W*C;B=C.T*V
 print(" small Gram ready",tag,time.time()-t,flush=True)
 mrow=(s['mL'][:,None]+mags([3]*3)[None,:]).ravel()
 mcol=(s['mE'][:,None]+mags([3]*2)[None,:]).ravel()
 require(all(V[i,j]==0 for i in range(n) for j in range(r) if mrow[i]!=mcol[j]),'global projector factor Sz selection')
 require(all(B[i,j]==0 for i in range(r) for j in range(r) if mcol[i]!=mcol[j]),'global projector Gram Sz selection')
 K=sp.zeros(n,r)
 for z in sorted(set(mcol)):
  rows=list(map(int,np.where(mrow==z)[0]));cols=list(map(int,np.where(mcol==z)[0]))
  b=B.extract(cols,cols);v=V.extract(rows,cols)
  require(b==b.T,'small projector Gram symmetric')
  # Bareiss/rational inversion; no floating inverses enter the target operator.
  inv=b.inv(method='DM');kk=v*inv
  require(kk*b==v,'exact projector right-factor solve')
  for i,ri in enumerate(rows):
   for j,cj in enumerate(cols):K[ri,cj]=kk[i,j]
 print('local exact',tag,n,r,'t',time.time()-t,flush=True)
 return {'V':V,'K':K,'B':B,'Cblock':C}

def prepare(work):
 """Build all exact inputs from polynomial definitions, overwriting generated caches."""
 from pathlib import Path
 work=Path(work);work.mkdir(parents=True,exist_ok=True)
 data,M=construct()
 with open(work/'rational_data.pkl','wb') as f:pickle.dump((data,M),f)
 for tag in data:
  with open(work/('local_'+tag+'.pkl'),'wb') as f:pickle.dump(make_local(data,M,tag),f)
 return data,M

if __name__=='__main__':
 from paths import WORK
 prepare(WORK)
 print('EXACT DATA REBUILT',flush=True)
