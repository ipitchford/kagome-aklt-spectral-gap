"""Deterministic candidate certificate generator and checker for hub angles."""
import os,sys,pickle,time,json,itertools,gc
from fractions import Fraction as F
import numpy as np
from scipy.linalg import cholesky,solve_triangular
from validated import *
from rational_hub import mags
from paths import WORK,EVIDENCE
ROOT=str(WORK)
_R_CACHE={}

def load():
 exact,M=pickle.load(open(ROOT+'/rational_data.pkl','rb'));M,em=exact_float(M)
 data={}
 for tag,s in exact.items():
  D,ed=exact_float(s['D']);G,eg=exact_float(s['diag']);cs=[]
  for c in s['C']:cs.append(exact_float(c)[0])
  data[tag]={**s,'D':D,'ed':ed,'g':G,'eg':eg,'Cs':cs}
 return data,M,em

def local_R(s,tag):
 if tag in _R_CACHE:return _R_CACHE[tag]
 t=time.time();v=pickle.load(open(ROOT+'/local_'+tag+'.pkl','rb'))
 V,ev=exact_float(v['V']);K,ek=exact_float(v['K']);del v
 er=matmul_error(K,ek,V.T,ev);R=K@V.T
 er=add(er,add(mul(gam(1),maxabs(R)),TINY));symmetrise(R)
 d=s['d'];R=R.reshape(d,3,3,3,d,3,3,3).transpose(1,2,3,0,5,6,7,4)
 print(' R ball',tag,'radius',er,'max',maxabs(R.reshape(27*d,27*d)),'t',time.time()-t,flush=True)
 _R_CACHE[tag]=(R,er)
 return R,er

def one_gram(A,M,em):
 e=A['e'];D=A['D'];G=np.zeros((e,3,3,e,3,3));err=0.;mm=maxabs(M);MM=M.reshape([3]*6)
 for a,ap in itertools.product(range(e),repeat=2):
  ii,jj=np.nonzero(D[a,ap]);cs=D[a,ap,ii,jj]
  for i,j,c in zip(ii,jj,cs):G[a,:,:,ap,:,:]+=c*MM[i,:,:,j,:,:]
  err=max(err,linear_error(cs,A['ed'],mm,em))
 G=G.reshape(e*9,e*9);err=add(err,add(mul(gam(1),maxabs(G)),TINY));symmetrise(G)
 return G,err

def two_gram(A,B,M,em):
 ea,eb=A['e'],B['e'];D1,D2=A['D'],B['D'];G=np.zeros((ea,eb,3,ea,eb,3));err=0.;mm=maxabs(M);MM=M.reshape([3]*6)
 for a,ap in itertools.product(range(ea),repeat=2):
  ii,jj=np.nonzero(D1[a,ap]);cs=D1[a,ap,ii,jj];W=np.zeros([3]*4)
  for i,j,c in zip(ii,jj,cs):W+=c*MM[i,:,:,j,:,:]
  ew=linear_error(cs,A['ed'],mm,em);mw=maxabs(W)
  for b,bp in itertools.product(range(eb),repeat=2):
   kk,ll=np.nonzero(D2[b,bp]);co=D2[b,bp,kk,ll]
   for k,l,c in zip(kk,ll,co):G[a,b,:,ap,bp,:]+=c*W[k,:,l,:]
   err=max(err,linear_error(co,B['ed'],mw,ew))
 G=G.reshape(ea*eb*3,ea*eb*3);err=add(err,add(mul(gam(1),maxabs(G)),TINY));symmetrise(G)
 return G,err

def finish(C,errrow,errcol,states,gs,G,eg,B,q,mode,tags,start):
 """Build q metric-C, add an arbitrary rank-k dyadic lift, and validate PD."""
 n=len(C);m=G.shape[0];d=B['d']
 divisions=math.nextafter(1-gam(max(1,n)),-INF)
 # Each row/column radius was a positive sum of at most n block contributions.
 er=div(float(np.max(errrow)),divisions);ec=div(float(np.max(errcol)),divisions)
 cnorm=norminf(C);cone=normone(C)
 cerr=add(mul(.5,add(er,ec)),mul(gam(1),mul(.5,add(cnorm,cone))))
 symmetrise(C)
 cnorm=norminf(C);gnorm=norminf(G);gmax=max(map(abs,B['g']));qf=float(q);eq=float(abs(F.from_float(qf)-q));eq=up(eq)
 # Full metric is a direct sum of g_c times principal submatrices of G.
 bnorm=mul(add(gmax,B['eg']),add(gnorm,mul(m,eg)))
 berr=add(mul(add(gmax,B['eg']),mul(m,eg)),mul(B['eg'],gnorm))
 targeterr=add(cerr,add(mul(add(abs(qf),eq),berr),mul(eq,bnorm)))
 targeterr=add(targeterr,mul(gam(4),add(cnorm,mul(abs(qf),mul(gmax,gnorm)))))
 # Common-coordinate matrix: gs is the list of boundary monomials in Mz=0.
 k=len(gs);Z=np.zeros((n,k));lookup={tuple(s):i for i,s in enumerate(states)}
 for col,source in enumerate(gs):
  if mode=='triple':
   a,b,ext=source
   for z,c in itertools.product(range(3),range(d)):
    row=lookup.get((a,b,z,c));v=B['Cs'][z][c,ext]
    if row is not None and v:Z[row,col]=v
  else:
   a,ext,y=source
   for z,c in itertools.product(range(3),range(d)):
    row=lookup.get((a,y,z,c));v=B['Cs'][z][c,ext]
    if row is not None and v:Z[row,col]=v
 BZ=np.zeros_like(Z);C*=-1
 for c in range(d):
  rows=np.where(states[:,3]==c)[0]
  if mode=='triple':idx=(states[rows,0]*B['_e2']+states[rows,1])*3+states[rows,2]
  else:idx=states[rows,0]*9+states[rows,1]*3+states[rows,2]
  block=G[np.ix_(idx,idx)]*B['g'][c]
  BZ[rows]=block@Z[rows]
  C[np.ix_(rows,rows)]+=qf*block
 gram=Z.T@BZ;gram=(gram+gram.T)/2
 chol=cholesky(gram,lower=True,check_finite=False)
 # Y is deliberately arbitrary dyadic data, not a claimed exact ground-space basis.
 Y=solve_triangular(chol,BZ.T,lower=True,check_finite=False).T*np.sqrt(2.)
 del Z,BZ,gram,chol,lookup
 ynorm=0.
 for i in range(0,n,192):ynorm=max(ynorm,float(np.max(np.sum(Y[i:i+192]**2,axis=1))))
 require(np.all(np.isfinite(Y)),'finite dyadic rank lift')
 ynorm=div(add(ynorm,TINY),math.nextafter(1-gam(2*k),-INF))
 anorm=norminf(C)
 targeterr=add(targeterr,add(mul(n,add(mul(gam(k),ynorm),TINY)),mul(gam(1),add(anorm,mul(n,ynorm)))))
 for i in range(0,n,192):C[i:i+192]+=Y[i:i+192]@Y.T
 del Y
 targeterr=add(targeterr,mul(gam(1),mul(.5,add(norminf(C),normone(C)))))
 symmetrise(C)
 print(' PSD target',tags,n,'common',k,'input norm radius',targeterr,'maxentry',maxabs(C),'t',time.time()-start,flush=True)
 result=shifted_cholesky(C,targeterr)
 result.update({'tags':tags,'mode':mode,'threshold_squared':str(q),'known_common_dimension':k,'rank_lift_columns':k,
 'elapsed_seconds':time.time()-start,'rational_tensor_input':True,'random_start':False})
 print(' CERTIFICATE',json.dumps(result),flush=True)
 return result

def triple(data,M,em,tags,q=F(1,20)):
 t=time.time();A,B,C=[data[x] for x in tags];ea,eb=A['e'],B['e'];d=C['d']
 G,eg=two_gram(A,B,M,em);R,er=local_R(C,tags[2]);mr=maxabs(R.reshape(27*d,27*d))
 groups={};states=[];n=0;mz=np.array([1,0,-1])[:,None]+C['mL'][None,:]
 for a,b in itertools.product(range(ea),range(eb)):
  z,c=np.where(mz==-A['mE'][a]-B['mE'][b]);ix=np.arange(n,n+len(z));n+=len(z);groups[a,b]=(ix,z,c)
  states.extend((a,b,int(zz),int(cc)) for zz,cc in zip(z,c))
 states=np.array(states);V=np.zeros((n,n));rr=np.zeros(n);cc=np.zeros(n);D1,D2=A['D'],B['D']
 for a,ap in itertools.product(range(ea),repeat=2):
  ii,jj=np.nonzero(D1[a,ap]);co=D1[a,ap,ii,jj]
  if not len(co):continue
  W=np.zeros((3,3,d,3,3,d))
  for i,j,v in zip(ii,jj,co):W+=v*R[i,:,:,:,j,:,:,:]
  ew=linear_error(co,A['ed'],mr,er);mw=maxabs(W)
  for b,bp in itertools.product(range(eb),repeat=2):
   ix,z,c=groups[a,b];iy,zp,cp=groups[ap,bp]
   if not len(ix) or not len(iy):continue
   kk,ll=np.nonzero(D2[b,bp]);co2=D2[b,bp,kk,ll]
   if not len(co2):continue
   val=np.zeros((len(ix),len(iy)))
   for k,l,v in zip(kk,ll,co2):val+=v*W[k,z[:,None],c[:,None],l,zp[None,:],cp[None,:]]
   V[np.ix_(ix,iy)]=val
   ev=linear_error(co2,B['ed'],mw,ew)
   rr[ix]+=mul(len(iy),ev);cc[iy]+=mul(len(ix),ev)
 gs=[(a,b,k) for a,b,k in itertools.product(range(ea),range(eb),range(C['e'])) if A['mE'][a]+B['mE'][b]+C['mE'][k]==0]
 del R,W
 return finish(V,rr,cc,states,gs,G,eg,{**C,'_e2':eb},q,'triple',tags,t)

def pair(data,M,em,tags,q=F(49,1600)):
 t=time.time();A,B=[data[x] for x in tags];e=A['e'];d=B['d']
 G,eg=one_gram(A,M,em);R,er=local_R(B,tags[1]);mr=maxabs(R.reshape(27*d,27*d))
 groups={};states=[];n=0;mz=mags([3,3])[:,None]+B['mL'][None,:]
 for a in range(e):
  yz,c=np.where(mz==-A['mE'][a]);y,z=yz//3,yz%3;ix=np.arange(n,n+len(c));n+=len(c);groups[a]=(ix,y,z,c)
  states.extend((a,int(yy),int(zz),int(v)) for yy,zz,v in zip(y,z,c))
 states=np.array(states);V=np.zeros((n,n));rr=np.zeros(n);cc=np.zeros(n);D=A['D']
 for a,ap in itertools.product(range(e),repeat=2):
  ix,y,z,c=groups[a];iy,yp,zp,cp=groups[ap]
  ii,jj=np.nonzero(D[a,ap]);co=D[a,ap,ii,jj]
  if not len(co):continue
  val=np.zeros((len(ix),len(iy)))
  for i,j,v in zip(ii,jj,co):val+=v*R[i,y[:,None],z[:,None],c[:,None],j,yp[None,:],zp[None,:],cp[None,:]]
  V[np.ix_(ix,iy)]=val;ev=linear_error(co,A['ed'],mr,er)
  rr[ix]+=mul(len(iy),ev);cc[iy]+=mul(len(ix),ev)
 gs=[(a,b,y) for a,b,y in itertools.product(range(e),range(B['e']),range(3)) if A['mE'][a]+B['mE'][b]+1-y==0]
 del R
 return finish(V,rr,cc,states,gs,G,eg,B,q,'pair',tags,t)

if __name__=='__main__':
 data,M,em=load();mode=sys.argv[1] if len(sys.argv)>1 else 'pair';tags=tuple(sys.argv[2:]) or ('mm','mm')
 result=pair(data,M,em,tags) if mode=='pair' else triple(data,M,em,tags)
 dest=EVIDENCE/('validated_'+mode+'_'+'_'.join(tags)+'.json')
 with open(dest,'w') as f:json.dump(result,f,indent=2)
