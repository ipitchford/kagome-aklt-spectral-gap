"""Complete exact bowtie spectrum by highest-weight and leaf-exchange blocks.
No FLINT, floating roots, tolerance exclusions or stored characteristic factors.
"""
from itertools import product
from math import comb,prod
from collections import defaultdict
from pathlib import Path
import time,json
import sympy as s
from block_gap import triangle_map
from paths import EVIDENCE

def require(v,msg):
 if not v:raise ArithmeticError(msg)

def orbit_basis(total,signs):
 reps=[];orbits=[]
 for v in product(range(5),repeat=5):
  if sum(v)!=total or v[1]>v[2] or v[3]>v[4]:continue
  if (signs[0]<0 and v[1]==v[2]) or (signs[1]<0 and v[3]==v[4]):continue
  out={}
  for a in range(2 if v[1]!=v[2] else 1):
   for b in range(2 if v[3]!=v[4] else 1):
    w=list(v)
    if a:w[1],w[2]=w[2],w[1]
    if b:w[3],w[4]=w[4],w[3]
    out[tuple(w)]=(signs[0] if a else 1)*(signs[1] if b else 1)
  reps.append(v);orbits.append(out)
 return reps,orbits

def run():
 t=time.time();psi,W,G,_=triangle_map();Q=psi*G.inv(method='DM')*psi.T*W
 require(Q*Q==Q and W*Q==Q.T*W,'triangle ground projector')
 rankE=sum(Q[a*25+b*5+c,a*25+b*5+c] for a,b,c,d,e in product(range(5),repeat=5) if a+b+c+d+e==10)
 require(rankE==105,'exact Mz0 trace of embedded triangle projector')
 st3=list(product(range(5),repeat=3));colmap={st3[j]:[(st3[i],Q[i,j]) for i in range(125) if Q[i,j]] for j in range(125)}
 x=s.Symbol('x');factors={};blocks=[];dim0=0;dimfull=0
 for signs in product([1,-1],repeat=2):
  for J in range(11):
   reps,orbits=orbit_basis(10-J,signs);targets,targetorbits=orbit_basis(9-J,signs)
   n=len(reps)
   if not n:continue
   H=s.zeros(n);R=s.zeros(len(targets),n);metric=s.diag(*[s.Rational(len(o),prod(comb(4,k) for k in v)) for o,v in zip(orbits,reps)])
   for j,orb in enumerate(orbits):
    h=defaultdict(lambda:s.S.Zero);raiseout=defaultdict(int)
    for src,coef in orb.items():
     for sites in [(0,1,2),(0,3,4)]:
      inp=tuple(src[k] for k in sites)
      for dest,a in colmap[inp]:
       dst=list(src)
       for k,b in zip(sites,dest):dst[k]=b
       h[tuple(dst)]+=coef*a
     for site in range(5):
      if src[site]:
       dst=list(src);dst[site]-=1;raiseout[tuple(dst)]+=coef*src[site]
    for i,(rep,orb2) in enumerate(zip(reps,orbits)):
     H[i,j]=h[rep]
     require(all(h[v]==sign*H[i,j] for v,sign in orb2.items()),'leaf symmetry of ground operator')
    for i,(rep,orb2) in enumerate(zip(targets,targetorbits)):
     R[i,j]=raiseout[rep]
     require(all(raiseout[v]==sign*R[i,j] for v,sign in orb2.items()),'leaf symmetry of raising operator')
   null=R.nullspace();d=len(null)
   if not d:continue
   T=s.Matrix.hstack(*null);piv=list(T.T.rref()[1]);A=T.extract(piv,range(d)).inv()*((H*T).extract(piv,range(d)))
   gram=T.T*metric*T
   require(H*T==T*A,'exact highest-weight invariance')
   require(R*T==s.zeros(R.rows,d),'highest-weight condition')
   require(gram*A==A.T*gram,'exact self-adjointness')
   cp=A.charpoly(x).as_expr();fac=s.factor_list(cp,x)[1]
   rec=[]
   for p,e in fac:
    pp=s.Poly(p,x).clear_denoms()[1].primitive()[1]
    if pp.LC()<0:pp=-pp
    key=str(pp.as_expr());factors.setdefault(key,{'poly':pp,'mz_multiplicity':0,'physical_multiplicity':0})
    factors[key]['mz_multiplicity']+=int(e);factors[key]['physical_multiplicity']+=(2*J+1)*int(e)
    rec.append({'coefficients':[str(c) for c in pp.all_coeffs()],'exponent':int(e)})
   blocks.append({'J':J,'leaf_exchange_signs':signs,'dimension':d,'factors':rec,'matrix':[[str(z) for z in row] for row in A.tolist()]})
   dim0+=d;dimfull+=(2*J+1)*d
  print('bowtie signs complete',signs,'t',time.time()-t,flush=True)
 require((dim0,dimfull)==(381,3125),'complete SU(2)/exchange coverage')
 top=s.Poly(69*x*x-185*x+108,x);noninteger=[];count=0
 for key,record in factors.items():
  p=record['poly']
  if p.degree()==1 and any(p.eval(k)==0 for k in (0,1,2)):continue
  require(all(p.eval(k)!=0 for k in (0,1,2)),'exact integer-root exclusion')
  require(p.eval(s.Rational(91,50))!=0,'strict rational ordering endpoint')
  nr=int(p.count_roots(s.Rational(91,50),s.Rational(2)))
  require(nr==0 or p==top,'no competing upper root')
  count+=nr;noninteger.append({'coefficients':[str(c) for c in p.all_coeffs()],'roots_above_91_over_50':nr,'mz_multiplicity':record['mz_multiplicity']})
 require(count==1,'unique distinct upper noninteger root')
 require(top.count_roots(s.Rational(91,50),s.Rational(183,100))==1,'rational top root isolator')
 def mult(k,which):
  return sum(r[which] for r in factors.values() if r['poly'].degree()==1 and r['poly'].eval(k)==0)
 require(mult(2,'mz_multiplicity')==19 and mult(0,'mz_multiplicity')==190,'common range/kernel counts')
 require(mult(2,'physical_multiplicity')==81,'full common ground dimension')
 return {'status':'PASS','arithmetic':'exact rational highest-weight blocks and rational Sturm counts','dimension_Mz0':dim0,'dimension_physical':dimfull,'largest_block':max(b['dimension'] for b in blocks),'common_range_Mz0':19,'common_kernel_Mz0':190,'rank_E_Mz0':105,'rank_F_Mz0':105,'eta':'(47+sqrt(4417))/138','top_sum_projector_polynomial':'69*x^2-185*x+108','top_sum_projector_isolator':['91/50','183/100'],'noninteger_factors':noninteger,'blocks':blocks,'elapsed_seconds':time.time()-t}

if __name__=='__main__':
 ans=run();(EVIDENCE/'bowtie_certified.json').write_text(json.dumps(ans,indent=2));print('BOWTIE EXACT PASS',ans['largest_block'],ans['elapsed_seconds'],flush=True)
