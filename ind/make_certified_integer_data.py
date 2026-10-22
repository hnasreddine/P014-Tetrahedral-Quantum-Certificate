"""Creates purely integer approximants to exact SOS linear system and rational inverse.
The algebraic entries were independently rebuilt in fast_exact_rows.json.gz.
All error bounds in the companion report follow from isqrt enclosures.
"""
import gzip,json,math,os
import numpy as np
from fractions import Fraction as Fr
P=os.path.dirname(__file__)
z=json.load(gzip.open(os.path.join(P,'fast_exact_rows.json.gz'),'rt'))
D=z['W_DEN'];QQ=6*D*D;SA=10**12;SR=10**12;SB=10**12;SY=10**11
BIG=10**40
r2=math.isqrt(2*BIG*BIG);r3=math.isqrt(3*BIG*BIG);r6=math.isqrt(6*BIG*BIG)
roots=(r2,r3,r6)
def approx_tuple(v,den,scale):
 # rigorous interval for f = a+b sqrt2+c sqrt3+d sqrt6 over denominator den
 a,b,c,d=v;lo=a*BIG;hi=a*BIG
 for coef,root in zip((b,c,d),roots):
  if coef>=0:lo+=coef*root;hi+=coef*(root+1)
  else:lo+=coef*(root+1);hi+=coef*root
 assert hi>=lo
 # nearest integer to interval midpoint; strict error <= 1/(2scale) + width/den/BIG
 value=(lo+hi)*scale//(2*den*BIG)
 # exact upper bound error in unit of scale
 err=max(abs(lo*scale-value*den*BIG),abs(hi*scale-value*den*BIG))
 return int(value),Fr(err,den*BIG*scale)
keep=z['used'];piv=z['pivots'];n=len(keep)
A=np.zeros((n,n),dtype=np.float64);maxerrA=Fr(0);rowerr=[Fr(0)]*n
A_sparse=[]
for j,i in enumerate(keep):
 nz=[]
 for k,field in z['allrows'][i][0]:
  val,err=approx_tuple(field,QQ,SA)
  if val:nz.append((k,val));A[j,k]=val/SA
  rowerr[j]+=err
 A_sparse.append(nz)
maxerrA=max(rowerr)
print('A max certified infinity error:',float(maxerrA),'A norm',np.max(np.sum(abs(A),axis=1)))
# b exact components have common denominator
BB=math.lcm(*(Fr(y).denominator for i in keep for y in z['allrows'][i][1]))
bint=[];maxerrB=Fr(0)
for i in keep:
 comp=[int(Fr(y)*BB) for y in z['allrows'][i][1]]
 val,err=approx_tuple(comp,QQ*BB,SB)
 bint.append(val);maxerrB=max(maxerrB,err)
print('b certified infinity error:',float(maxerrB),'max value',max(map(abs,bint)))
# get exact rational rounded-grams xhat (not numpy-recreated)
raw=open(os.path.join(P,'..','prior','rounded_grams.txt')).read().split();it=iter(raw)
nblock=int(next(it));Xhat=[];margin=int(5e-5*SY)
for name in z['blocks']:
 label=next(it);nr=int(next(it));r=nr//2;assert label==name
 H=[[int(next(it)) for j in range(nr)] for k in range(nr)]
 for j in range(r):Xhat.append(H[j][j]+margin)
 for j in range(r):
  for k in range(j+1,r):Xhat.append(H[j][k]);Xhat.append(H[r+j][k])
assert len(Xhat)==3175
xhat=[Xhat[k] for k in piv]
err=0
for i in range(n):err=max(err,abs(bint[i]/SB-sum(A[i,j]*xhat[j]/SY for j in range(n))))
print('nearby pivot residual estimate',err,'max |xhat|',max(map(abs,xhat))/SY)
R=np.linalg.inv(A)
Rint=np.rint(R*SR).astype(np.int64)
print('float inverse rho',np.max(np.sum(abs(np.eye(n)-Rint.astype(float)@A/SR),axis=1)))
# ensure no INT64 overflow
assert max(abs(int(v)) for row in Rint for v in row)<2**61
assert max(abs(v) for row in A_sparse for _,v in row)<2**62
fname=os.path.join(P,'integer_certificate_data.txt')
with open(fname,'w') as f:
 f.write(f'{n} {SA} {SR} {SB} {SY}\n')
 f.write(str(maxerrA.numerator)+' '+str(maxerrA.denominator)+' '+str(maxerrB.numerator)+' '+str(maxerrB.denominator)+'\n')
 for row in A_sparse:
  f.write(str(len(row))+' '+' '.join(f'{j} {v}' for j,v in row)+'\n')
 for row in Rint:
  f.write(' '.join(str(int(v)) for v in row)+'\n')
 f.write(' '.join(map(str,bint))+'\n')
 f.write(' '.join(map(str,xhat))+'\n')
print('saved',fname,'size',os.stat(fname).st_size)
