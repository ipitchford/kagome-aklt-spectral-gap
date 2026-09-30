"""Exact-input certification of the weighted nine-site block via two ground spaces."""
import itertools,time,json,os,math
from fractions import Fraction as F
from collections import defaultdict
import sympy as sp
import numpy as np
from scipy.linalg import cholesky,solve_triangular,eigh
from rational_hub import exbasis,mags,rational_frame
from validated import *
from paths import EVIDENCE

def triangle_map():
 poly={(0,0,0):1}
 for a,b in [(0,1),(0,2),(1,2)]:
  nxt=defaultdict(int)
  for ex,c in poly.items():
   t=list(ex);t[b]+=1;nxt[tuple(t)]+=c
   t=list(ex);t[a]+=1;nxt[tuple(t)]-=c
  poly={x:c for x,c in nxt.items() if c}
 psi=sp.zeros(125,27)
 for ex in itertools.product(range(3),repeat=3):
  j=np.ravel_multi_index(ex,(3,3,3))
  for p,c in poly.items():
   dst=tuple(a+b for a,b in zip(ex,p));i=np.ravel_multi_index(dst,(5,5,5));psi[i,j]+=c
 w=[sp.Rational(1,math.prod(math.comb(4,k) for k in x)) for x in itertools.product(range(5),repeat=3)]
 W=sp.diag(*w);G=psi.T*W*psi
 return psi,W,G,poly

def setup():
 t=time.time();psi,W,G,omega=triangle_map();exactR=W*psi*G.inv(method='DM')*psi.T*W
 R,er=exact_float(exactR)
 # Order each row/column spin pair adjacently, ready for three contractions.
 R=R.reshape([5]*6).transpose(0,3,1,4,2,5).reshape(25,625)
 out={}
 for tag,sgn in [('p',1),('m',-1)]:
  UE=sp.Matrix(exbasis(3,sgn));U=sp.kronecker_product(sp.eye(3),UE);d,e=U.cols,UE.cols
  me=np.array([mags([3,3])[next(i for i in range(9) if UE[i,j])] for j in range(e)])
  ml=(np.array([1,0,-1])[:,None]+me[None,:]).ravel()
  Q,Qi,g=rational_frame(U.T*G*U,ml);V=psi*U*Q
  D=np.empty((d,d,5,5),dtype=object)
  for a,b in itertools.product(range(5),repeat=2):
   VA=V[a*25:(a+1)*25,:];VB=V[b*25:(b+1)*25,:]
   weights=sp.diag(*[sp.Rational(1,math.comb(4,x)*math.comb(4,y)) for x,y in itertools.product(range(5),repeat=2)])
   D[:,:,a,b]=np.array(VA.T*weights*VB,dtype=object)
  df,ed=exact_float(D);gf,eg=exact_float(g)
  out[tag]={'d':d,'e':e,'m':ml,'me':me,'D':df.reshape(d*d,25),'ed':ed,'g':gf,'eg':eg,'Qinv':np.array(Qi,dtype=float)}
 print('block exact tensors',time.time()-t,'Rradius',er,flush=True)
 return out,R,er,omega

def contraction(A,B,C,R,er):
 d1,d2,d3=A['d'],B['d'],C['d']
 et=matmul_error(A['D'],A['ed'],R,er);T=A['D']@R
 T=T.reshape(d1*d1,25,25).transpose(0,2,1).reshape(d1*d1*25,25)
 eu=matmul_error(T,et,B['D'].T,B['ed']);V=T@B['D'].T
 V=V.reshape(d1*d1,25,d2*d2).transpose(0,2,1).reshape(d1*d1*d2*d2,25)
 ev=matmul_error(V,eu,C['D'].T,C['ed']);X=V@C['D'].T
 X=X.reshape(d1,d1,d2,d2,d3,d3).transpose(0,2,4,1,3,5).reshape(d1*d2*d3,d1*d2*d3)
 return X,ev

def check(data,R,er,omega,tags,q=F(2,3)):
 start=time.time();A,B,C=[data[t] for t in tags];d1,d2,d3=A['d'],B['d'],C['d'];e1,e2,e3=A['e'],B['e'],C['e']
 X,ex=contraction(A,B,C,R,er)
 mm=(A['m'][:,None,None]+B['m'][None,:,None]+C['m'][None,None,:]).ravel();ix=np.where(mm==0)[0]
 H=X[np.ix_(ix,ix)].copy();del X
 n=len(ix);cerr=add(mul(n,ex),mul(gam(1),mul(.5,add(norminf(H),normone(H)))));symmetrise(H)
 states=np.array(list(itertools.product(range(d1),range(d2),range(d3))))[ix]
 g=np.kron(np.kron(A['g'],B['g']),C['g'])[ix]
 # Diagonal product conversion and two multiplication errors.
 maxg=max(g);maxprod=mul(max(A['g']),mul(max(B['g']),max(C['g'])));eb=mul(gam(2),maxprod)
 for I,J,K in [(A,B,C),(B,A,C),(C,A,B)]:eb=add(eb,mul(I['eg'],mul(add(max(J['g']),J['eg']),add(max(K['g']),K['eg']))))
 gs=[s for s in itertools.product(range(e1),range(e2),range(e3)) if A['me'][s[0]]+B['me'][s[1]]+C['me'][s[2]]==0]
 gsidx=np.array([np.ravel_multi_index(s,(e1,e2,e3)) for s in gs]);k=len(gs)
 Z=np.zeros((n,k))
 for (h1,h2,h3),coef in omega.items():
  if not coef:continue
  ZA=A['Qinv'][:,h1*e1:(h1+1)*e1];ZB=B['Qinv'][:,h2*e2:(h2+1)*e2];ZC=C['Qinv'][:,h3*e3:(h3+1)*e3]
  Z+=coef*np.kron(np.kron(ZA,ZB),ZC)[np.ix_(ix,gsidx)]
 BZ=g[:,None]*Z;gram=Z.T@BZ;gram=(gram+gram.T)*.5
 L=cholesky(gram,lower=True);Y=solve_triangular(L,BZ.T,lower=True).T*np.sqrt(2.)
 qf=float(q);eq=up(float(abs(F.from_float(qf)-q)))
 err=add(cerr,add(mul(add(abs(qf),eq),eb),mul(eq,maxg)))
 err=add(err,mul(gam(3),add(norminf(H),mul(abs(qf),maxg))))
 # Optional diagnostic, not part of acceptance.
 Q0=Z@np.linalg.inv(cholesky(gram,lower=True).T)
 V=H/np.sqrt(g[:,None]*g[None,:]);QQ=np.sqrt(g[:,None])*Q0
 estimate=eigh((V-QQ@QQ.T+V.T-QQ@QQ.T)*.5,subset_by_index=[n-1,n-1],eigvals_only=True)[0]
 H*=-1;H[np.diag_indices(n)]+=qf*g
 require(np.all(np.isfinite(Y)),"finite dyadic rank lift")
 yn=div(add(float(np.max(np.sum(Y*Y,axis=1))),TINY),math.nextafter(1-gam(2*k),-INF))
 err=add(err,add(mul(n,add(mul(gam(k),yn),TINY)),mul(gam(1),add(norminf(H),mul(n,yn)))))
 H+=Y@Y.T
 err=add(err,mul(gam(1),mul(.5,add(norminf(H),normone(H)))));symmetrise(H)
 result=shifted_cholesky(H,err)
 result.update({'mode':'nine_site_block_angle','tags':tags,'threshold_squared':str(q),'known_common_dimension':k,'rank_lift_columns':k,'overlap_squared_estimate_not_used_in_acceptance':float(estimate),'elapsed_seconds':time.time()-start})
 print('BLOCK',json.dumps(result),flush=True)
 return result

if __name__=='__main__':
 data,R,er,omega=setup();results=[]
 for tags in itertools.combinations_with_replacement(['p','m'],3):results.append(check(data,R,er,omega,tags))
 path=EVIDENCE/'block_gap_validated.json'
 with open(path,'w') as f:json.dump({'status':'PASS' if all(x['status']=='PASS' for x in results) else 'REJECT','results':results},f,indent=2)

 require(all(x["status"]=="PASS" for x in results),"block-angle certificate rejected")
