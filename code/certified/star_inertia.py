from fractions import Fraction as F
from itertools import product
from math import comb, prod
from pathlib import Path
import subprocess,json,hashlib,time
from paths import WORK,EVIDENCE
from validated import require
C=[comb(4,k) for k in range(5)]
sts=[s for s in product(range(5),repeat=5) if sum(s)==10]; ix={s:i for i,s in enumerate(sts)};n=len(sts)
W=[prod(C[x] for x in s) for s in sts]; D=7000; shift=F(6,125)
A=[[0]*n for _ in sts]
for j,s in enumerate(sts):
 for l in range(1,5):
  tot=s[0]+s[l]
  for k in range(max(0,tot-4),min(4,tot)+1):
   t=list(s);t[0]=k;t[l]=tot-k;i=ix[tuple(t)]
   r=D*W[i]*F(C[s[0]]*C[s[l]],comb(8,tot))
   require(r.denominator==1,"integer matrix construction")
   A[i][j]+=r.numerator
 A[j][j]-=int(D*W[j]*shift)
text=str(n)+'\n'+'\n'.join(' '.join(map(str,row)) for row in A)+'\n'
(WORK/'star_matrix.txt').write_text(text)
t=time.time();run=subprocess.run([str(WORK/'inertia')],input=text,text=True,capture_output=True,check=True);data=json.loads(run.stdout)
require((data['negative'],data['positive'],data['zero'])==(45,336,0),'star inertia count')
data.update({'source':'independently reconstructed rational physical projector, not stored factorisation','matrix_sha256':hashlib.sha256(text.encode()).hexdigest(),'threshold':'6/125','scale':D,'elapsed_seconds':time.time()-t,'ground_dimension_exact':44,'positive_eigenvalues_below_threshold':1})
(EVIDENCE/'star_inertia.json').write_text(json.dumps(data,indent=2));print(json.dumps(data))
