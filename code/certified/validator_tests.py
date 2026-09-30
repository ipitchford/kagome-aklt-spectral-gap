"""Adversarial arithmetic tests and a small all-rational kagome certificate."""
import json,pickle,subprocess,itertools,math
from fractions import Fraction as F
import numpy as np
import sympy as sp
import validated as v
from validated import require
from paths import WORK,EVIDENCE

# Every entry of every test product is compared to a rational exact product.
rng=np.random.default_rng(20260930);count=0
for k in [1,2,5,11]:
 for trial in range(5):
  A=np.array([[F(int(rng.integers(-100,101)),int(rng.integers(1,40))) for j in range(k)] for i in range(3)],dtype=object)
  B=np.array([[F(int(rng.integers(-100,101)),int(rng.integers(1,40))) for j in range(4)] for i in range(k)],dtype=object)
  af,ea=v.exact_float(A);bf,eb=v.exact_float(B);err=v.matmul_error(af,ea,bf,eb);c=af@bf;truth=A@B
  require(all(abs(F.from_float(c[ix])-truth[ix])<=F.from_float(err) for ix in np.ndindex(c.shape)),'every exact product entry enclosed');count+=c.size
for A,B in [(np.array([[F(2)**-600]],dtype=object),np.array([[F(2)**-600]],dtype=object)),(np.array([[F(2)**40,F(1),-F(2)**40]],dtype=object),np.array([[F(1)],[F(2)**-40],[F(1)]],dtype=object))]:
 af,ea=v.exact_float(A);bf,eb=v.exact_float(B);err=v.matmul_error(af,ea,bf,eb);c=af@bf
 require(abs(F.from_float(c[0,0])-(A@B)[0,0])<=F.from_float(err),'underflow/cancellation enclosure');count+=1
require(v.shifted_cholesky(np.diag([1.,2.,3.]),0.)['status']=='PASS','positive control')
require(v.shifted_cholesky(np.diag([1.,-1.,3.]),0.)['status']!='PASS','indefinite rejected')
require(v.shifted_cholesky(np.diag([-.5,-.5,1.])+np.array([[2.,0,0],[0,0,0],[0,0,0]]),0.)['status']!='PASS','insufficient-rank lift rejected')
old=v.cholesky
try:
 v.cholesky=lambda A,**kw:np.eye(len(A))
 fake=v.shifted_cholesky(np.diag([1.,-1.,3.]),0.)
 require(fake['status']=='REJECT_ERROR_BUDGET','bogus Cholesky factor caught by full residual')
finally:v.cholesky=old
try:v.shifted_cholesky(np.array([[1.,.25],[0.,1.]]),0.)
except ArithmeticError:pass
else:raise ArithmeticError('asymmetric input accepted')

# The end-to-end numerical checker must reject a deliberately too-strong angle bound.
import certify_hub as h
data,M,em=h.load();neg=h.triple(data,M,em,('mm','mm','mm'),q=F(3,100))
require(neg['status']!='PASS','overstrong triple threshold rejected')
pos=h.triple(data,M,em,('mm','mm','mm'))
require(pos['status']=='PASS','small triple positive control')
# Rebuild the same small target independently, ENTRY BY ENTRY over rationals;
# verify inertia over integers, using no floating-point acceptance at all.
exact,hub=pickle.load(open(WORK/'rational_data.pkl','rb'));s=exact['mm'];local=pickle.load(open(WORK/'local_mm.pkl','rb'))
R=local['K']*local['V'].T;d,e=s['d'],s['e'];D=s['D'];g=s['diag']
states=[(a,b,z,c) for a,b,z,c in itertools.product(range(e),range(e),range(3),range(d)) if s['mE'][a]+s['mE'][b]+1-z+s['mL'][c]==0]
n=len(states);require(n==83,'small exact dimension')
chan={(a,ap):[(i,j,D[a,ap,i,j]) for i,j in itertools.product(range(3),repeat=2) if D[a,ap,i,j]] for a,ap in itertools.product(range(e),repeat=2)}
B=sp.zeros(n);C=sp.zeros(n)
for row,(a,b,z,c) in enumerate(states):
 for col,(ap,bp,zp,cp) in enumerate(states):
  terms=[(i,ip,j,jp,aa*bb) for i,ip,aa in chan[a,ap] for j,jp,bb in chan[b,bp]]
  C[row,col]=sum(co*R[c*27+i*9+j*3+z,cp*27+ip*9+jp*3+zp] for i,ip,j,jp,co in terms)
  if c==cp:B[row,col]=g[c]*sum(co*hub[i*9+j*3+z,ip*9+jp*3+zp] for i,ip,j,jp,co in terms)
require(B==B.T and C==C.T,'independent exact compressed forms symmetric')
gs=[(a,b,k) for a,b,k in itertools.product(range(e),repeat=3) if s['mE'][a]+s['mE'][b]+s['mE'][k]==0]
Z=sp.Matrix(n,len(gs),lambda i,j:s['C'][states[i][2]][states[i][3],gs[j][2]] if states[i][:2]==gs[j][:2] else 0)
require(len(gs)==7 and Z.rank()==7,'seven independent exact common vectors')
require(C*Z==B*Z,'exact compressed projector equals identity on common range')
BZ=B*Z;lift=2*BZ*(Z.T*BZ).inv()*BZ.T
checks=[]
for q in [sp.Rational(1,20),sp.Rational(3,100)]:
 target=q*B-C+lift;den=sp.ilcm(*[x.q for x in target]);integers=target*den
 text=str(n)+'\n'+'\n'.join(' '.join(str(int(x)) for x in row) for row in integers.tolist())+'\n'
 run=subprocess.run([str(WORK/'inertia')],input=text,text=True,capture_output=True,check=True);r=json.loads(run.stdout)
 if q==sp.Rational(1,20):require((r['positive'],r['negative'],r['zero'])==(83,0,0),'all-rational triple certificate')
 else:require(r['negative']>0,'all-rational countercheck to stronger threshold')
 r.update(threshold_squared=str(q));checks.append(r)
res={'status':'PASS','arithmetic_product_entries_checked':count,'fake_factor_rejected':True,'asymmetric_input_rejected':True,'insufficient_rank_lift_rejected':True,'overstrong_threshold_rejected':True,'small_all_rational_certificates':checks,'python_optimization_does_not_disable_checks':True}
(EVIDENCE/'validator_tests.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
