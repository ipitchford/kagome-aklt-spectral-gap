"""Exploratory symmetry-reduced intersection-angle computation. Float64, not a certificate."""
import sys,itertools,time
import numpy as np
from scipy.linalg import eigh
# kvbs.py is the original v0.2 implementation copied beside this file.
from kvbs import Region,up,down,sqrtm_psd

def exchange_basis(d,sign):
 cols=[]
 for i in range(d):
  for j in range(i,d):
   if i==j:
    if sign<0:continue
    v=np.zeros((d,d));v[i,j]=1
   else:
    v=np.zeros((d,d));v[i,j]=1/np.sqrt(2);v[j,i]=sign/np.sqrt(2)
   cols.append(v.ravel())
 return np.array(cols).T

def symmetry_bases(half):
 pair={s:exchange_basis(3,s) for s in (1,-1)}
 if half:pair={s:np.kron(U,np.eye(2)) for s,U in pair.items()}
 out={}
 for s,tag in [(1,'p'),(-1,'m')]:
  d=pair[s].shape[1]
  for t,sg in [(1,'p'),(-1,'m')]:out[tag+sg]=np.kron(pair[s],pair[s])@exchange_basis(d,t)
 out['x']=np.kron(pair[1],pair[-1])
 return out

def mags(ds):
 v=np.zeros(1)
 for d in ds:v=(v[:,None]+(d-1)/2-np.arange(d)[None,:]).ravel()
 return v

def setup():
 L=Region(up(1,0)+up(1,-1));hub=Region(up(0,0))
 conn=[s for s in L.leg_sites if len(L.dangling[s])==1]
 groups=[]
 for v in conn:
  tri=next(t for t in [up(1,0),up(1,-1)] if v in t)
  groups.append(sorted(s for s in tri if s!=v)+[v])
 order=sum(groups,[]); perm=[L.leg_sites.index(s) for s in order]
 G=L.gram().reshape(L.leg_dims*2).transpose(perm+[6+i for i in perm]).reshape(324,324)
 M=hub.gram();M/=np.trace(M)/27
 T,_,_=sqrtm_psd(M)
 UL=symmetry_bases(True);UE=symmetry_bases(False)
 # A: external four spin-one legs -> private plus one hub spin-one leg.
 S=np.zeros((2,2,3));S[1,1,0]=S[0,0,2]=1/np.sqrt(3)
 S[0,1,1]=S[1,0,1]=-1/np.sqrt(6)
 C=np.zeros((3,3,2,3,3,2,3,3,3,3,3))
 for i,j,k,l in itertools.product(range(3),repeat=4):C[i,j,:,k,l,:,:,i,j,k,l]=S
 C=C.reshape(324,3,81)
 ml=mags([3,3,2,3,3,2]);me=mags([3]*4)
 out={}
 for tag in UL:
  U=UL[tag];V=UE[tag]
  gl=U.T@G@U; X,_,ev=sqrtm_psd(gl)
  A=np.einsum('li,lhk,ke->ihe',U,C,V,optimize=True)
  A=np.einsum('ij,jhk->ihk',X,A,optimize=True)
  Aflat=A.reshape(-1,V.shape[1]);XX,inv,_=sqrtm_psd(Aflat.T@Aflat)
  A=(Aflat@inv).reshape(U.shape[1],3,V.shape[1])
  mL=np.rint(np.diag(U.T@(ml[:,None]*U))).astype(int)
  mE=np.rint(np.diag(V.T@(me[:,None]*V))).astype(int)
  # Force exact magnetisation sparsity; record removed roundoff.
  mask=mL[:,None,None]+np.array([1,0,-1])[None,:,None]!=mE[None,None,:]
  defect=np.max(abs(A[mask]));A[mask]=0
  chan=np.einsum('iha,ikb->abhk',A,A,optimize=True)
  out[tag]={'A':A,'D':chan,'mL':mL,'mE':mE,'d':U.shape[1],'e':V.shape[1],'gram_min':ev[0], 'm_defect':float(defect)}
  # Check raw symmetry-sector invariance and exact dimension counting independently later.
  print(tag,U.shape[1],V.shape[1],ev[0],ev[-1],defect,flush=True)
 return out,M,T

def remaining_projector(R,T,leg):
 # V on private L x hub, column external e x the two other hub legs.
 A=R['A'];d,_,e=A.shape
 V=np.zeros((d,3,3,3,e,3,3))
 other=[i for i in range(3) if i!=leg]
 for u,v in itertools.product(range(3),repeat=2):
  sl=[slice(None)]*7;sl[1+other[0]]=u;sl[1+other[1]]=v;sl[5]=u;sl[6]=v
  V[tuple(sl)]=A
 V=V.reshape(d,27,-1);V=np.einsum('uv,ivk->iuk',T,V,optimize=True).reshape(d*27,-1)
 Q=(V@sqrtm_psd(V.T@V)[1]).reshape(d,27,-1)
 Z=np.einsum('uv,ivk->iuk',T,Q,optimize=True).reshape(d*27,-1)
 # Z Z^T = T G T; arrangement h1,h2,h3,L on each side
 Rmat=(Z@Z.T).reshape(d,3,3,3,d,3,3,3).transpose(1,2,3,0,5,6,7,4)
 return Rmat

def two_map_gram(A,B,M):
 # Gram of (A1 x A2 x I_h3) under hub metric M
 D1,D2=A['D'],B['D']; e1,e2=A['e'],B['e']
 G=np.einsum('abik,cdjl,ijxkly->acxbdy',D1,D2,M.reshape([3]*6),optimize=True)
 G=G.reshape(e1*e2*3,e1*e2*3)
 return (G+G.T)/2

def triple_angle(data,M,T,tags=('mm','mm','mm'),do_eigen=True):
 t=time.time();A,B,C=[data[x] for x in tags]; e1,e2,e3=A['e'],B['e'],C['e'];d=C['d']
 G=two_map_gram(A,B,M);X,Xinv,_=sqrtm_psd(G)
 # Preserve exact Sz support in whitening.
 mg=(A['mE'][:,None,None]+B['mE'][None,:,None]+np.array([1,0,-1])[None,None,:]).ravel()
 X[mg[:,None]!=mg[None,:]]=0;Xinv[mg[:,None]!=mg[None,:]]=0
 R=remaining_projector(C,T,2)
 # Reduced Mz=0 layout: first indices (a,b), then z,private c.
 mz=np.array([1,0,-1])[:,None]+C['mL'][None,:]
 n=0;groups={};states=[]
 for a,b in itertools.product(range(e1),range(e2)):
  zs,cs=np.where(mz==-A['mE'][a]-B['mE'][b]); ids=np.arange(n,n+len(zs));n+=len(zs)
  groups[a,b]=(ids,zs,cs)
  states.extend((a,b,int(z),int(c)) for z,c in zip(zs,cs))
 states=np.array(states);print('triple',tags,'N0',n,'t',time.time()-t,flush=True)
 V=np.zeros((n,n));D1=A['D'];D2=B['D']
 # R axes (h1,h2,z,c,h1prime,h2prime,zprime,cprime)
 for a in range(e1):
  for ap in range(e1):
   coeff=D1[a,ap];ii,jj=np.where(abs(coeff)>1e-14)
   W=sum((coeff[i,j]*R[i,:,:,:,j,:,:,:] for i,j in zip(ii,jj)),np.zeros((3,3,d,3,3,d)))
   # W h2,z,c,h2prime,zprime,cprime
   for b in range(e2):
    ix,z,c=groups[a,b]
    if not len(ix):continue
    for bp in range(e2):
     iy,zp,cp=groups[ap,bp]
     if not len(iy):continue
     coeff2=D2[b,bp]; kk,ll=np.where(abs(coeff2)>1e-14)
     val=np.zeros((len(ix),len(iy)))
     for k,l in zip(kk,ll):val+=coeff2[k,l]*W[k,z[:,None],c[:,None],l,zp[None,:],cp[None,:]]
     V[np.ix_(ix,iy)]=val
  if a%5==0:print(' build row',a,'/',e1,'t',time.time()-t,flush=True)
 print('raw matrix ready',time.time()-t,'sym',np.max(abs(V-V.T)),flush=True)
 # Whitening groups for each private basis vector c.
 for c in range(d):
  rows=np.where(states[:,3]==c)[0]
  gidx=(states[rows,0]*e2+states[rows,1])*3+states[rows,2]
  W=Xinv[np.ix_(gidx,gidx)]
  V[rows,:]=W@V[rows,:]
 for c in range(d):
  rows=np.where(states[:,3]==c)[0]
  gidx=(states[rows,0]*e2+states[rows,1])*3+states[rows,2]
  W=Xinv[np.ix_(gidx,gidx)]
  V[:,rows]=V[:,rows]@W
 # Common space: X acting on a,b,z applied to I_e1 I_e2 A3(c,z,k)
 gstates=[(a,b,k) for a,b,k in itertools.product(range(e1),range(e2),range(e3)) if A['mE'][a]+B['mE'][b]+C['mE'][k]==0]
 Bc=np.zeros((n,len(gstates)))
 for j,(a,b,k) in enumerate(gstates):
  rows,z,c=groups[a,b];Bc[rows,j]=C['A'][c,z,k]
 for c in range(d):
  rows=np.where(states[:,3]==c)[0];gidx=(states[rows,0]*e2+states[rows,1])*3+states[rows,2]
  Bc[rows,:]=X[np.ix_(gidx,gidx)]@Bc[rows,:]
 if len(gstates):
  Bc=Bc@sqrtm_psd(Bc.T@Bc)[1]
  containment=np.linalg.norm(V@Bc-Bc,ord='fro');print('common',len(gstates),'resid',containment,flush=True)
  V-=Bc@Bc.T
 else:containment=0
 V=(V+V.T)/2
 from scipy.sparse.linalg import eigsh
 ev=eigsh(V,k=1,which='LA',tol=1e-10,return_eigenvectors=False)[0]
 print('theta',np.sqrt(ev),'theta2',ev,'total',time.time()-t,flush=True)
 cp=None
 return {'tags':tags,'N0':n,'common':len(gstates),'theta_estimate':float(np.sqrt(ev)),'theta2_estimate':float(ev),'float_cholesky_1_16':cp,'common_residual_F':containment,'elapsed_seconds':time.time()-t},V

if __name__=='__main__':
 import json
 data,M,T=setup()
 tags=tuple(sys.argv[1:4]) if len(sys.argv)>1 else ('mm','mm','mm')
 out,V=triple_angle(data,M,T,tags)
 print(json.dumps(out))
