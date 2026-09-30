"""Outward scalar error bounds and an all-entry shifted-Cholesky PSD check.

The Cholesky routine is an untrusted candidate generator: acceptance checks its
full residual and does not assume that successful factorisation proves positivity.
IEEE binary64 arithmetic, conventional (non-Strassen) BLAS dot products.
"""
import math,ctypes
from fractions import Fraction as F
import numpy as np
from scipy.linalg import cholesky
U=2.0**-53;TINY=2.0**-900;INF=float('inf')

def require(x,msg):
 if not x:raise ArithmeticError(msg)

def up(x):
 require(math.isfinite(x),'nonfinite scalar bound')
 return math.nextafter(float(x),INF)

def add(a,b):return up(a+b)
def mul(a,b):return up(a*b)
def div(a,b):
 require(b>0,'positive divisor');return up(a/b)
def gam(k):
 require(0<=k<2**48,'dot product length supported')
 return 4*k*U

def sum_bound(values):
 values=np.asarray(values);require(np.all(np.isfinite(values)),'finite summands')
 s=float(np.sum(np.abs(values)))
 den=math.nextafter(1-gam(values.size),-INF)
 return div(s,den)

def maxabs(A):
 ans=0.
 for i in range(0,len(A),256):
  v=np.asarray(A[i:i+256]);require(np.all(np.isfinite(v)),'finite matrix')
  if v.size:ans=max(ans,float(np.max(np.abs(v))))
 return ans

def norminf(A):
 ans=0.;n=A.shape[1]
 for i in range(0,len(A),256):
  v=np.asarray(A[i:i+256]);require(np.all(np.isfinite(v)),'finite matrix')
  if v.size:ans=max(ans,float(np.max(np.sum(np.abs(v),axis=1))))
 return div(ans,math.nextafter(1-gam(n),-INF))

def normone(A):return norminf(A.T)

def exact_float(A):
 """Conversion verified against each exact rational, not assumed correctly rounded."""
 A=np.asarray(A,dtype=object);B=np.empty(A.shape,dtype=float);err=F(0)
 for ix in np.ndindex(A.shape):
  a=A[ix];x=float(a);require(math.isfinite(x),'finite rational conversion')
  B[ix]=x
  if a:
   require(x!=0,'nonzero rational must not underflow in conversion')
   p,q=a.as_numer_denom() if hasattr(a,'as_numer_denom') else (a.numerator,a.denominator) if isinstance(a,F) else (int(a),1)
   er=abs(F.from_float(x)-F(int(p),int(q)));err=max(err,er)
  else:require(x==0,'zero conversion')
 eps=float(err)
 if F.from_float(eps)<err:eps=math.nextafter(eps,INF)
 return B,eps

def matmul_error(A,ea,B,eb):
 k=A.shape[1];ra=norminf(A);cb=normone(B);ma=maxabs(A);mb=maxabs(B)
 round_term=mul(gam(k),min(mul(ra,mb),mul(ma,cb)))
 return add(add(add(round_term,mul(ea,cb)),mul(eb,ra)),add(mul(k,mul(ea,eb)),TINY))

def linear_error(coeff,ec,mx,ex):
 """Uniform radius for sum_i coeff_i X_i, |X_i|<=mx, radius<=ex."""
 coeff=np.asarray(coeff);s=coeff.size
 if not s:return 0.
 c1=sum_bound(coeff)
 return add(add(mul(gam(s),mul(c1,mx)),mul(c1,ex)),add(mul(s,mul(ec,add(mx,ex))),TINY))

def symmetrise(A):
 n=len(A);step=192
 for i in range(0,n,step):
  for j in range(i,n,step):
   v=(A[i:i+step,j:j+step]+A[j:j+step,i:i+step].T)*0.5
   A[i:i+step,j:j+step]=v;A[j:j+step,i:i+step]=v.T
 return A

def check_environment():
 require(np.finfo(float).eps==2**-52,'binary64 epsilon')
 require(ctypes.CDLL(None).fegetround()==0,'round-to-nearest required')
 # Backend identity is included in the release's environment receipt.

def shifted_cholesky(A,target_error,shift=None):
 """Certify the exact symmetric target K>0 from ||K-A||_2<=target_error.

Returns PASS only if shift > input error + diagonal subtraction error + a
validated bound on ||(A-shift I)-L L^T||_infinity. Every residual entry is checked.
The source A is overwritten by its diagonal shift; L is not exported.
"""
 check_environment();n=len(A);require(A.shape==(n,n),'square matrix')
 ma=maxabs(A);require(target_error>=0,'nonnegative input radius')
 # Pick a dyadic shift larger than the budget; its adequacy is checked below.
 if shift is None:
  preliminary=add(target_error,mul(n*gam(n),max(1.,ma)))
  shift=2.0**math.ceil(math.log2(max(2**-40,16*preliminary)))
 require(shift>0 and math.isfinite(shift),'finite positive shift')
 require(np.array_equal(A,A.T),'exactly symmetric floating candidate required')
 diag=np.diag_indices(n);A[diag]-=shift
 subtraction_error=add(mul(gam(1),add(ma,shift)),TINY)
 try:L=cholesky(A,lower=True,check_finite=False)
 except np.linalg.LinAlgError:
  return {'status':'REJECT_CHOLESKY','dimension':n,'shift':shift,'target_error_bound':target_error}
 require(np.all(np.isfinite(L)),'finite Cholesky factor')
 require(np.all(np.diag(L)>0),'positive Cholesky pivots')
 # Exact row energy <= computed sum of squares /(1-gamma).
 maxenergy=0.
 for i in range(0,n,192):
  b=L[i:i+192];v=float(np.max(np.sum(b*b,axis=1)));maxenergy=max(maxenergy,v)
 maxenergy=div(add(maxenergy,TINY),math.nextafter(1-gam(2*n),-INF))
 residual_rows=0.;mashift=add(ma,shift)
 for i in range(0,n,192):
  v=A[i:i+192]-L[i:i+192]@L.T
  require(np.all(np.isfinite(v)),'finite full residual')
  residual_rows=max(residual_rows,float(np.max(np.sum(np.abs(v),axis=1))))
 residual_rows=div(residual_rows,math.nextafter(1-gam(n),-INF))
 dot_error=mul(n,add(mul(gam(n),maxenergy),TINY))
 subtract_error=mul(n,add(mul(gam(1),add(mashift,mul(1+gam(n),maxenergy))),TINY))
 residual_bound=add(residual_rows,add(dot_error,subtract_error))
 budget=add(target_error,add(subtraction_error,residual_bound))
 answer={'status':'PASS' if shift>budget else 'REJECT_ERROR_BUDGET','dimension':n,'shift':shift,
  'target_error_bound':target_error,'shift_roundoff_bound':subtraction_error,
  'residual_center_infinity_upper':residual_rows,'residual_roundoff_upper':add(dot_error,subtract_error),
  'residual_total_upper':residual_bound,'total_error_upper':budget,
  'positive_margin_lower':math.nextafter(shift-budget,-INF),'max_factor_row_energy_upper':maxenergy,
  'method':'all-entry residual; rank-lift; exact-rational tensor input'}
 del L
 return answer
